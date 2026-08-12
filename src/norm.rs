mod atomize;
mod cardinality;
mod duplicate;
mod merge;
mod reference;
mod split;

use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use noodles_vcf::{Header, variant::RecordBuf, variant::io::Write as _};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::expression::Compiled;
use crate::filter::Logic;
use crate::format::{
    HeaderMode, OutputFormat, Reader, RecordScratch, VariantWriter, Writer, trim_line_ending,
};
use crate::regions::RegionSelection;
pub(crate) use duplicate::Policy as DuplicatePolicy;
pub(crate) use merge::Policy as JoinPolicy;
pub(crate) use reference::MismatchPolicy;
use reference::{Outcome, ReferenceNormalizer};

#[derive(Clone, Debug)]
pub(crate) struct Options {
    pub(crate) reference: Option<PathBuf>,
    pub(crate) expression: Option<String>,
    pub(crate) expression_logic: Logic,
    pub(crate) targets: Option<RegionSelection>,
    pub(crate) split_multiallelic: bool,
    pub(crate) join_multiallelic: Option<JoinPolicy>,
    pub(crate) strict_filter: bool,
    pub(crate) split_overlaps_missing: bool,
    pub(crate) mismatch_policy: MismatchPolicy,
    pub(crate) atomize: bool,
    pub(crate) atom_overlaps_star: bool,
    pub(crate) old_record_tag: Option<String>,
    pub(crate) duplicate_policy: Option<DuplicatePolicy>,
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
    pub(crate) joined: u64,
    pub(crate) not_selected: u64,
    pub(crate) target_filtered: u64,
    pub(crate) skipped: u64,
    pub(crate) atomized: u64,
    pub(crate) duplicates: u64,
    pub(crate) output_format: OutputFormat,
}

struct Pending {
    reference: usize,
    position: usize,
    serial: u64,
    input: u64,
    selected: bool,
    record: RecordBuf,
}

struct OutputState {
    position: Option<(usize, usize)>,
    duplicates: duplicate::State,
}

struct OutputOptions<'a> {
    old_record_tag: Option<&'a str>,
    join_multiallelic: Option<JoinPolicy>,
    strict_filter: bool,
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
    let (mut header, _, _) = reader.read_header()?;
    if let Some(tag) = options.old_record_tag.as_deref() {
        atomize::prepare_header(&mut header, tag)?;
    }
    let expression = options
        .expression
        .as_deref()
        .map(|source| {
            Compiled::bind(source, &header).map_err(|error| {
                RsomicsError::ConfigError(format!("invalid norm expression: {error}"))
            })
        })
        .transpose()?;
    let mut writer = Writer::new(output, options.output_format);
    writer.write_header(&header, HeaderMode::Full)?;
    let summary = normalize_stream(
        &mut reader,
        &header,
        options,
        expression.as_ref(),
        &mut writer,
    )?;
    writer.finish()?;
    Ok(summary)
}

fn normalize_stream(
    reader: &mut Reader,
    header: &Header,
    options: &Options,
    expression: Option<&Compiled>,
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
    let mut output_state = OutputState {
        position: None,
        duplicates: duplicate::State::new(options.duplicate_policy),
    };
    let output_options = OutputOptions {
        old_record_tag: options.old_record_tag.as_deref(),
        join_multiallelic: options.join_multiallelic,
        strict_filter: options.strict_filter,
    };
    let mut serial = 0;
    let mut summary = Summary {
        read: 0,
        written: 0,
        changed: 0,
        unchanged: 0,
        unsupported: 0,
        split: 0,
        joined: 0,
        not_selected: 0,
        target_filtered: 0,
        skipped: 0,
        atomized: 0,
        duplicates: 0,
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

        if options
            .targets
            .as_ref()
            .is_some_and(|targets| !targets.keeps(&record))
        {
            summary.target_filtered += 1;
            let threshold = position.saturating_sub(options.site_window);
            flush_ready(
                &mut pending,
                header,
                writer,
                &output_options,
                (reference, threshold),
                &mut output_state,
                &mut summary,
            )?;
            continue;
        }

        let selected = expression
            .map(|expression| {
                expression
                    .evaluate(header, &record)
                    .map(|truth| options.expression_logic.accepts(truth.site_passes()))
                    .map_err(|error| invalid(number, &format!("evaluating expression: {error}")))
            })
            .transpose()?
            .unwrap_or(true);
        summary.not_selected += u64::from(!selected);
        let origin = (selected && options.split_multiallelic && options.old_record_tag.is_some())
            .then(|| record.clone());
        let records = if selected && options.split_multiallelic {
            split::split(
                header,
                &record,
                options.keep_sum_ad,
                options.split_overlaps_missing,
            )?
        } else {
            vec![record]
        };
        if selected && records.len() > 1 {
            summary.split += 1;
        }
        for (split_index, mut record) in records.into_iter().enumerate() {
            if selected && let Some(normalizer) = &mut normalizer {
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
            let (records, atomized) = if selected && options.atomize {
                atomize::atomize(
                    header,
                    record,
                    options.atom_overlaps_star,
                    options.old_record_tag.as_deref(),
                    origin.as_ref().map(|record| atomize::Origin {
                        record,
                        alternate: split_index + 1,
                    }),
                )?
            } else {
                (vec![record], false)
            };
            summary.atomized += u64::from(atomized);
            for record in records {
                let normalized_position = record
                    .variant_start()
                    .map(usize::from)
                    .ok_or_else(|| invalid(number, "normalized variant position is missing"))?;
                pending.push(Reverse(Pending {
                    reference,
                    position: normalized_position,
                    serial,
                    input: number,
                    selected,
                    record,
                }));
                serial = serial.checked_add(1).ok_or_else(|| {
                    RsomicsError::InvalidInput("output record count exceeds u64".to_owned())
                })?;
            }
        }

        let threshold = position.saturating_sub(options.site_window);
        flush_ready(
            &mut pending,
            header,
            writer,
            &output_options,
            (reference, threshold),
            &mut output_state,
            &mut summary,
        )?;
    }

    flush_all(
        &mut pending,
        header,
        writer,
        &output_options,
        &mut output_state,
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
    options: &OutputOptions<'_>,
    through: (usize, usize),
    state: &mut OutputState,
    summary: &mut Summary,
) -> Result<()> {
    while pending.peek().is_some_and(|Reverse(record)| {
        record.reference < through.0
            || (record.reference == through.0 && record.position <= through.1)
    }) {
        let records = pop_coordinate(pending);
        write_coordinate(records, header, writer, options, state, summary)?;
    }
    Ok(())
}

fn flush_all(
    pending: &mut BinaryHeap<Reverse<Pending>>,
    header: &Header,
    writer: &mut impl VariantWriter,
    options: &OutputOptions<'_>,
    state: &mut OutputState,
    summary: &mut Summary,
) -> Result<()> {
    while !pending.is_empty() {
        let records = pop_coordinate(pending);
        write_coordinate(records, header, writer, options, state, summary)?;
    }
    Ok(())
}

fn pop_coordinate(pending: &mut BinaryHeap<Reverse<Pending>>) -> Vec<Pending> {
    let Reverse(first) = pending.pop().unwrap();
    let coordinate = (first.reference, first.position);
    let mut records = vec![first];
    while pending
        .peek()
        .is_some_and(|Reverse(record)| (record.reference, record.position) == coordinate)
    {
        records.push(pending.pop().unwrap().0);
    }
    records
}

fn write_coordinate(
    mut records: Vec<Pending>,
    header: &Header,
    writer: &mut impl VariantWriter,
    options: &OutputOptions<'_>,
    state: &mut OutputState,
    summary: &mut Summary,
) -> Result<()> {
    let selected_indices = records
        .iter()
        .enumerate()
        .filter_map(|(index, record)| record.selected.then_some(index))
        .collect::<Vec<_>>();
    if let Some(policy) = options
        .join_multiallelic
        .filter(|_| selected_indices.len() > 1)
    {
        let (joined, count) = merge::join(
            policy,
            options.strict_filter,
            header,
            selected_indices.iter().map(|index| &records[*index].record),
        )?;
        let mut sources = records.into_iter().map(Some).collect::<Vec<_>>();
        let merged = joined
            .into_iter()
            .map(|(index, record)| {
                let mut pending = sources[selected_indices[index]].take().unwrap();
                pending.record = record;
                pending
            })
            .collect::<Vec<_>>();
        for index in selected_indices {
            sources[index] = None;
        }
        records = sources.into_iter().flatten().chain(merged).collect();
        summary.joined += count;
    }
    for record in records {
        write_pending(
            record,
            header,
            writer,
            options.old_record_tag,
            state,
            summary,
        )?;
    }
    Ok(())
}

fn write_pending(
    pending: Pending,
    header: &Header,
    writer: &mut impl VariantWriter,
    old_record_tag: Option<&str>,
    state: &mut OutputState,
    summary: &mut Summary,
) -> Result<()> {
    let position = (pending.reference, pending.position);
    if state.position.is_some_and(|previous| position < previous) {
        return Err(invalid(
            pending.input,
            "normalized position exceeds --site-window; increase the window",
        ));
    }
    if pending.selected && state.duplicates.remove(position, &pending.record) {
        summary.duplicates += 1;
        return Ok(());
    }
    let number = pending.serial + 1;
    if let Some(tag) = old_record_tag
        .filter(|tag| writer.supports_vcf_records() && pending.record.info().get(*tag).is_some())
    {
        let mut record = Vec::new();
        noodles_vcf::io::Writer::new(&mut record)
            .write_variant_record(header, &pending.record)
            .map_err(|error| {
                RsomicsError::Io(std::io::Error::new(
                    error.kind(),
                    format!("rendering variant record {number}: {error}"),
                ))
            })?;
        trim_line_ending(&mut record);
        restore_origin_commas(&mut record, tag)?;
        writer.write_vcf_record(&record, number)?;
    } else {
        writer.write_record(header, &pending.record, number)?;
    }
    state.position = Some(position);
    summary.written += 1;
    Ok(())
}

fn restore_origin_commas(record: &mut Vec<u8>, tag: &str) -> Result<()> {
    let info_start = record
        .iter()
        .enumerate()
        .filter(|(_, byte)| **byte == b'\t')
        .nth(6)
        .map(|(index, _)| index + 1)
        .ok_or_else(|| {
            RsomicsError::InvalidInput("rendered VCF record has no INFO field".to_owned())
        })?;
    let info_end = record[info_start..]
        .iter()
        .position(|byte| *byte == b'\t')
        .map_or(record.len(), |offset| info_start + offset);
    let prefix = format!("{tag}=");
    let mut field_start = info_start;
    let (value_start, value_end) = loop {
        let field_end = record[field_start..info_end]
            .iter()
            .position(|byte| *byte == b';')
            .map_or(info_end, |offset| field_start + offset);
        if record[field_start..field_end].starts_with(prefix.as_bytes()) {
            break (field_start + prefix.len(), field_end);
        }
        if field_end == info_end {
            return Err(RsomicsError::InvalidInput(format!(
                "rendered VCF record is missing INFO/{tag}"
            )));
        }
        field_start = field_end + 1;
    };
    let mut restored = Vec::with_capacity(record.len());
    restored.extend_from_slice(&record[..value_start]);
    let mut remaining = &record[value_start..value_end];
    while let Some(offset) = remaining.windows(3).position(|window| window == b"%2C") {
        restored.extend_from_slice(&remaining[..offset]);
        restored.push(b',');
        remaining = &remaining[offset + 3..];
    }
    restored.extend_from_slice(remaining);
    restored.extend_from_slice(&record[value_end..]);
    *record = restored;
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
            expression: None,
            expression_logic: Logic::Include,
            targets: None,
            split_multiallelic: false,
            join_multiallelic: None,
            strict_filter: false,
            split_overlaps_missing: false,
            mismatch_policy: MismatchPolicy::Exit,
            atomize: false,
            atom_overlaps_star: true,
            old_record_tag: None,
            duplicate_policy: None,
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
            expression: None,
            expression_logic: Logic::Include,
            targets: None,
            split_multiallelic: false,
            join_multiallelic: None,
            strict_filter: false,
            split_overlaps_missing: false,
            mismatch_policy: MismatchPolicy::Exit,
            atomize: false,
            atom_overlaps_star: true,
            old_record_tag: None,
            duplicate_policy: None,
            keep_sum_ad: false,
            output_format: OutputFormat::Vcf,
            site_window: 1000,
        };
        let error = write(&input, &options, Vec::new()).unwrap_err().to_string();
        assert!(error.contains("not coordinate sorted"), "{error}");
    }
}
