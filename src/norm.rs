mod atomize;
mod reference;
mod split;

use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use noodles_vcf::{Header, variant::RecordBuf};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::format::{HeaderMode, OutputFormat, Reader, RecordScratch, VariantWriter, Writer};
pub(crate) use reference::MismatchPolicy;
use reference::{Outcome, ReferenceNormalizer};

#[derive(Clone, Debug)]
pub(crate) struct Options {
    pub(crate) reference: Option<PathBuf>,
    pub(crate) split_multiallelic: bool,
    pub(crate) mismatch_policy: MismatchPolicy,
    pub(crate) atomize: bool,
    pub(crate) keep_sum_ad: bool,
    pub(crate) output_format: OutputFormat,
    pub(crate) site_window: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct Summary {
    pub(crate) read: u64,
    pub(crate) written: u64,
    pub(crate) changed: u64,
    pub(crate) unchanged: u64,
    pub(crate) unsupported: u64,
    pub(crate) split: u64,
    pub(crate) skipped: u64,
    pub(crate) atomized: u64,
    pub(crate) output_format: OutputFormat,
}

struct Pending {
    reference: usize,
    position: usize,
    serial: u64,
    input: u64,
    record: RecordBuf,
}

impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for Pending {}

impl PartialOrd for Pending {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Pending {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key().cmp(&other.key())
    }
}

impl Pending {
    fn key(&self) -> (usize, usize, u64) {
        (self.reference, self.position, self.serial)
    }
}

pub(crate) fn write(input: &Path, options: &Options, output: impl Write) -> Result<Summary> {
    if options.site_window == 0 {
        return Err(RsomicsError::ConfigError(
            "--site-window must be at least 1".to_owned(),
        ));
    }

    let mut reader = Reader::open(input)?;
    let (header, _, _) = reader.read_header()?;
    let mut writer = Writer::new(output, options.output_format);
    writer.write_header(&header, HeaderMode::Full)?;
    let summary = normalize_stream(&mut reader, &header, options, &mut writer)?;
    writer.finish()?;
    Ok(summary)
}

fn normalize_stream(
    reader: &mut Reader,
    header: &Header,
    options: &Options,
    writer: &mut impl VariantWriter,
) -> Result<Summary> {
    let reference_order: HashMap<_, _> = header
        .contigs()
        .keys()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    let mut normalizer = options
        .reference
        .as_deref()
        .map(|path| ReferenceNormalizer::open(path, options.mismatch_policy))
        .transpose()?;
    split::validate(header, options.keep_sum_ad)?;
    let mut scratch = RecordScratch::default();
    let mut pending = BinaryHeap::new();
    let mut seen = HashSet::new();
    let mut input_position = None;
    let mut output_position = None;
    let mut summary = Summary {
        read: 0,
        written: 0,
        changed: 0,
        unchanged: 0,
        unsupported: 0,
        split: 0,
        skipped: 0,
        atomized: 0,
        output_format: options.output_format,
    };

    loop {
        let number = summary.read + 1;
        let Some(record) = reader.read_record(header, &mut scratch, number)? else {
            break;
        };
        summary.read += 1;

        let reference = reference_order
            .get(record.reference_sequence_name())
            .copied()
            .ok_or_else(|| invalid(number, "reference sequence is absent from the header"))?;
        let position = record
            .variant_start()
            .map(usize::from)
            .ok_or_else(|| invalid(number, "variant position is missing"))?;
        validate_input_order(number, reference, position, &mut input_position, &mut seen)?;

        let records = if options.split_multiallelic {
            split::split(header, &record, options.keep_sum_ad)?
        } else {
            vec![record]
        };
        if records.len() > 1 {
            summary.split += 1;
        }
        for mut record in records {
            if let Some(normalizer) = &mut normalizer {
                match normalizer.normalize(&mut record)? {
                    Outcome::Changed => summary.changed += 1,
                    Outcome::Unchanged => summary.unchanged += 1,
                    Outcome::Unsupported => summary.unsupported += 1,
                    Outcome::Skipped => {
                        summary.skipped += 1;
                        continue;
                    }
                }
            }
            let (records, atomized) = if options.atomize {
                atomize::atomize(record)?
            } else {
                (vec![record], false)
            };
            summary.atomized += u64::from(atomized);
            for record in records {
                let normalized_position = record
                    .variant_start()
                    .map(usize::from)
                    .ok_or_else(|| invalid(number, "normalized variant position is missing"))?;
                let serial = summary
                    .written
                    .checked_add(u64::try_from(pending.len()).map_err(|_| {
                        RsomicsError::InvalidInput("pending record count exceeds u64".to_owned())
                    })?)
                    .ok_or_else(|| {
                        RsomicsError::InvalidInput("output record count exceeds u64".to_owned())
                    })?;
                pending.push(Reverse(Pending {
                    reference,
                    position: normalized_position,
                    serial,
                    input: number,
                    record,
                }));
            }
        }

        let threshold = position.saturating_sub(options.site_window);
        flush_ready(
            &mut pending,
            header,
            writer,
            reference,
            threshold,
            &mut output_position,
            &mut summary,
        )?;
    }

    flush_all(
        &mut pending,
        header,
        writer,
        &mut output_position,
        &mut summary,
    )?;
    Ok(summary)
}

fn validate_input_order(
    number: u64,
    reference: usize,
    position: usize,
    previous: &mut Option<(usize, usize)>,
    seen: &mut HashSet<usize>,
) -> Result<()> {
    if let Some((previous_reference, previous_position)) = *previous {
        if reference < previous_reference
            || (reference == previous_reference && position < previous_position)
        {
            return Err(invalid(number, "input records are not coordinate sorted"));
        }
        if reference != previous_reference && !seen.insert(previous_reference) {
            return Err(invalid(
                number,
                "reference sequence blocks are not contiguous",
            ));
        }
    }
    *previous = Some((reference, position));
    Ok(())
}

fn flush_ready(
    pending: &mut BinaryHeap<Reverse<Pending>>,
    header: &Header,
    writer: &mut impl VariantWriter,
    reference: usize,
    threshold: usize,
    output_position: &mut Option<(usize, usize)>,
    summary: &mut Summary,
) -> Result<()> {
    while pending.peek().is_some_and(|Reverse(record)| {
        record.reference < reference
            || (record.reference == reference && record.position <= threshold)
    }) {
        let Reverse(record) = pending.pop().unwrap();
        write_pending(record, header, writer, output_position, summary)?;
    }
    Ok(())
}

fn flush_all(
    pending: &mut BinaryHeap<Reverse<Pending>>,
    header: &Header,
    writer: &mut impl VariantWriter,
    output_position: &mut Option<(usize, usize)>,
    summary: &mut Summary,
) -> Result<()> {
    while let Some(Reverse(record)) = pending.pop() {
        write_pending(record, header, writer, output_position, summary)?;
    }
    Ok(())
}

fn write_pending(
    pending: Pending,
    header: &Header,
    writer: &mut impl VariantWriter,
    output_position: &mut Option<(usize, usize)>,
    summary: &mut Summary,
) -> Result<()> {
    let position = (pending.reference, pending.position);
    if output_position.is_some_and(|previous| position < previous) {
        return Err(invalid(
            pending.input,
            "normalized position exceeds --site-window; increase the window",
        ));
    }
    writer.write_record(header, &pending.record, pending.serial + 1)?;
    *output_position = Some(position);
    summary.written += 1;
    Ok(())
}

fn invalid(number: u64, message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("normalizing variant record {number}: {message}"))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    fn fixture() -> (tempfile::TempDir, PathBuf, PathBuf) {
        let directory = tempfile::tempdir().unwrap();
        let reference = directory.path().join("reference.fa");
        let input = directory.path().join("input.vcf");
        fs::write(&reference, b">chr1\nAAAAAACGTACGT\n").unwrap();
        fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
        fs::write(
            &input,
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tA\tC\t.\tPASS\t.\n\
chr1\t4\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.\n",
        )
        .unwrap();
        (directory, reference, input)
    }

    #[test]
    fn normalizes_and_reorders_records_within_the_site_window() {
        let (_directory, reference, input) = fixture();
        let options = Options {
            reference: Some(reference),
            split_multiallelic: false,
            mismatch_policy: MismatchPolicy::Exit,
            atomize: false,
            keep_sum_ad: false,
            output_format: OutputFormat::Vcf,
            site_window: 1000,
        };
        let mut output = Vec::new();
        let summary = write(&input, &options, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert_eq!(summary.read, 3);
        assert_eq!(summary.changed, 2);
        assert_eq!(summary.unchanged, 1);
        let records: Vec<_> = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect();
        assert_eq!(records[0].split('\t').nth(1), Some("1"));
        assert_eq!(records[1].split('\t').nth(1), Some("2"));
        assert_eq!(records[2].split('\t').nth(1), Some("11"));
    }

    #[test]
    fn rejects_unsorted_input_without_partial_success() {
        let (_directory, reference, input) = fixture();
        fs::write(
            &input,
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.\n\
chr1\t4\t.\tA\tAA\t.\tPASS\t.\n",
        )
        .unwrap();
        let options = Options {
            reference: Some(reference),
            split_multiallelic: false,
            mismatch_policy: MismatchPolicy::Exit,
            atomize: false,
            keep_sum_ad: false,
            output_format: OutputFormat::Vcf,
            site_window: 1000,
        };
        let error = write(&input, &options, Vec::new()).unwrap_err().to_string();
        assert!(error.contains("not coordinate sorted"), "{error}");
    }
}
