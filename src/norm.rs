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
use crate::regions::{IndexedRecords, RegionSelection, RegionSet};
pub(crate) use duplicate::Policy as DuplicatePolicy;
pub(crate) use merge::Policy as JoinPolicy;
pub(crate) use reference::MismatchPolicy;
use reference::{Outcome, ReferenceNormalizer};

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum SortOrder {
    Position,
    Lexicographic,
}

#[derive(Clone, Debug)]
pub(crate) struct Options {
    pub(crate) reference: Option<PathBuf>,
    pub(crate) expression: Option<String>,
    pub(crate) expression_logic: Logic,
    pub(crate) regions: Option<RegionSet>,
    pub(crate) targets: Option<RegionSelection>,
    pub(crate) sort: SortOrder,
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
    sort: SortOrder,
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

struct Normalizer<'a, W> {
    header: &'a Header,
    options: &'a Options,
    expression: Option<&'a Compiled>,
    writer: &'a mut W,
    reference_order: HashMap<&'a str, usize>,
    reference_normalizer: Option<ReferenceNormalizer>,
    pending: BinaryHeap<Reverse<Pending>>,
    seen: HashSet<usize>,
    input_position: Option<(usize, usize)>,
    output_state: OutputState,
    output_options: OutputOptions<'a>,
    serial: u64,
    summary: Summary,
}

impl PartialEq for Pending {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
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
        self.coordinate()
            .cmp(&other.coordinate())
            .then_with(|| self.sort.cmp(&other.sort))
            .then_with(|| {
                if self.sort == SortOrder::Lexicographic {
                    compare_alleles(&self.record, &other.record)
                } else {
                    Ordering::Equal
                }
            })
            .then_with(|| self.serial.cmp(&other.serial))
    }
}

fn compare_alleles(left: &RecordBuf, right: &RecordBuf) -> Ordering {
    compare_ascii_case_insensitive(
        left.reference_bases().as_bytes(),
        right.reference_bases().as_bytes(),
    )
    .then_with(|| {
        let left = left.alternate_bases();
        let right = right.alternate_bases();
        left.as_ref()
            .iter()
            .zip(right.as_ref())
            .find_map(|(left, right)| {
                let order = compare_ascii_case_insensitive(left.as_bytes(), right.as_bytes());
                (order != Ordering::Equal).then_some(order)
            })
            .unwrap_or_else(|| left.as_ref().len().cmp(&right.as_ref().len()))
    })
}

fn compare_ascii_case_insensitive(left: &[u8], right: &[u8]) -> Ordering {
    left.iter()
        .map(u8::to_ascii_lowercase)
        .cmp(right.iter().map(u8::to_ascii_lowercase))
}

impl Pending {
    fn coordinate(&self) -> (usize, usize) {
        (self.reference, self.position)
    }
}

pub(crate) fn write(input: &Path, options: &Options, output: impl Write) -> Result<Summary> {
    if options.site_window == 0 {
        return Err(RsomicsError::ConfigError(
            "--site-window must be at least 1".to_owned(),
        ));
    }
    let mut writer = Writer::new(output, options.output_format);
    let summary = if let Some(regions) = &options.regions {
        if input == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "indexed regions require a named input".to_owned(),
            ));
        }
        normalize_indexed(input, regions, options, &mut writer)?
    } else {
        normalize_stream(input, options, &mut writer)?
    };
    writer.finish()?;
    Ok(summary)
}

fn normalize_stream(
    input: &Path,
    options: &Options,
    writer: &mut impl VariantWriter,
) -> Result<Summary> {
    let mut reader = Reader::open(input)?;
    let (mut header, _, _) = reader.read_header()?;
    prepare_header(&mut header, options)?;
    let expression = bind_expression(&header, options)?;
    let mut normalizer = Normalizer::new(&header, options, expression.as_ref(), writer)?;
    normalizer.write_header()?;
    let mut scratch = RecordScratch::default();
    loop {
        let number = normalizer.summary.read + 1;
        let Some(record) = reader.read_record(&header, &mut scratch, number)? else {
            break;
        };
        normalizer.push(record, number)?;
    }
    normalizer.finish()?;
    Ok(normalizer.summary)
}

fn normalize_indexed(
    input: &Path,
    regions: &RegionSet,
    options: &Options,
    writer: &mut impl VariantWriter,
) -> Result<Summary> {
    let mut reader = IndexedRecords::open(input, regions)?;
    let mut header = reader.header().clone();
    prepare_header(&mut header, options)?;
    let expression = bind_expression(&header, options)?;
    let mut normalizer = Normalizer::new(&header, options, expression.as_ref(), writer)?;
    normalizer.write_header()?;
    let read = reader.visit(|_, record, number| normalizer.push(record, number))?;
    normalizer.summary.read = read;
    normalizer.finish()?;
    Ok(normalizer.summary)
}

fn prepare_header(header: &mut Header, options: &Options) -> Result<()> {
    if let Some(tag) = options.old_record_tag.as_deref() {
        atomize::prepare_header(header, tag)?;
    }
    Ok(())
}

fn bind_expression(header: &Header, options: &Options) -> Result<Option<Compiled>> {
    options
        .expression
        .as_deref()
        .map(|source| {
            Compiled::bind(source, header).map_err(|error| {
                RsomicsError::ConfigError(format!("invalid norm expression: {error}"))
            })
        })
        .transpose()
}

impl<'a, W: VariantWriter> Normalizer<'a, W> {
    fn new(
        header: &'a Header,
        options: &'a Options,
        expression: Option<&'a Compiled>,
        writer: &'a mut W,
    ) -> Result<Self> {
        split::validate(header, options.keep_sum_ad)?;
        let reference_normalizer = options
            .reference
            .as_deref()
            .map(|path| ReferenceNormalizer::open(path, options.mismatch_policy))
            .transpose()?;
        Ok(Self {
            header,
            options,
            expression,
            writer,
            reference_order: header
                .contigs()
                .keys()
                .enumerate()
                .map(|(index, name)| (name.as_str(), index))
                .collect(),
            reference_normalizer,
            pending: BinaryHeap::new(),
            seen: HashSet::new(),
            input_position: None,
            output_state: OutputState {
                position: None,
                duplicates: duplicate::State::new(options.duplicate_policy),
            },
            output_options: OutputOptions {
                old_record_tag: options.old_record_tag.as_deref(),
                join_multiallelic: options.join_multiallelic,
                strict_filter: options.strict_filter,
            },
            serial: 0,
            summary: Summary {
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
            },
        })
    }

    fn write_header(&mut self) -> Result<()> {
        self.writer.write_header(self.header, HeaderMode::Full)
    }

    fn push(&mut self, record: RecordBuf, number: u64) -> Result<()> {
        self.summary.read = self.summary.read.max(number);
        let reference = self
            .reference_order
            .get(record.reference_sequence_name())
            .copied()
            .ok_or_else(|| invalid(number, "reference sequence is absent from the header"))?;
        let position = record
            .variant_start()
            .map(usize::from)
            .ok_or_else(|| invalid(number, "variant position is missing"))?;
        validate_input_order(
            number,
            reference,
            position,
            &mut self.input_position,
            &mut self.seen,
        )?;

        if self
            .options
            .targets
            .as_ref()
            .is_some_and(|targets| !targets.keeps(&record))
        {
            self.summary.target_filtered += 1;
            let threshold = position.saturating_sub(self.options.site_window);
            flush_ready(
                &mut self.pending,
                self.header,
                self.writer,
                &self.output_options,
                (reference, threshold),
                &mut self.output_state,
                &mut self.summary,
            )?;
            return Ok(());
        }

        let selected = self
            .expression
            .map(|expression| {
                expression
                    .evaluate(self.header, &record)
                    .map(|truth| self.options.expression_logic.accepts(truth.site_passes()))
                    .map_err(|error| invalid(number, &format!("evaluating expression: {error}")))
            })
            .transpose()?
            .unwrap_or(true);
        self.summary.not_selected += u64::from(!selected);
        let origin =
            (selected && self.options.split_multiallelic && self.options.old_record_tag.is_some())
                .then(|| record.clone());
        let records = if selected && self.options.split_multiallelic {
            split::split(
                self.header,
                &record,
                self.options.keep_sum_ad,
                self.options.split_overlaps_missing,
            )?
        } else {
            vec![record]
        };
        if selected && records.len() > 1 {
            self.summary.split += 1;
        }
        for (split_index, mut record) in records.into_iter().enumerate() {
            if selected && let Some(normalizer) = &mut self.reference_normalizer {
                match normalizer.normalize(&mut record)? {
                    Outcome::Changed => self.summary.changed += 1,
                    Outcome::Unchanged => self.summary.unchanged += 1,
                    Outcome::Unsupported => self.summary.unsupported += 1,
                    Outcome::Skipped => {
                        self.summary.skipped += 1;
                        continue;
                    }
                }
            }
            let (records, atomized) = if selected && self.options.atomize {
                atomize::atomize(
                    self.header,
                    record,
                    self.options.atom_overlaps_star,
                    self.options.old_record_tag.as_deref(),
                    origin.as_ref().map(|record| atomize::Origin {
                        record,
                        alternate: split_index + 1,
                    }),
                )?
            } else {
                (vec![record], false)
            };
            self.summary.atomized += u64::from(atomized);
            for record in records {
                let normalized_position = record
                    .variant_start()
                    .map(usize::from)
                    .ok_or_else(|| invalid(number, "normalized variant position is missing"))?;
                self.pending.push(Reverse(Pending {
                    reference,
                    position: normalized_position,
                    serial: self.serial,
                    input: number,
                    selected,
                    sort: self.options.sort,
                    record,
                }));
                self.serial = self.serial.checked_add(1).ok_or_else(|| {
                    RsomicsError::InvalidInput("output record count exceeds u64".to_owned())
                })?;
            }
        }

        let threshold = position.saturating_sub(self.options.site_window);
        flush_ready(
            &mut self.pending,
            self.header,
            self.writer,
            &self.output_options,
            (reference, threshold),
            &mut self.output_state,
            &mut self.summary,
        )
    }

    fn finish(&mut self) -> Result<()> {
        flush_all(
            &mut self.pending,
            self.header,
            self.writer,
            &self.output_options,
            &mut self.output_state,
            &mut self.summary,
        )
    }
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
            regions: None,
            targets: None,
            sort: SortOrder::Position,
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
            regions: None,
            targets: None,
            sort: SortOrder::Position,
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
