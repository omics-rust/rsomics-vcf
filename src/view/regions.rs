use std::io::Write;
use std::path::Path;

use noodles_core::{Position, Region, region::Interval};
use noodles_util::variant::io::indexed_reader;
use noodles_vcf::variant::{RecordBuf, record::info::field::key, record_buf::info::field::Value};
use rsomics_common::{Context, Result, RsomicsError};

use crate::format::Writer;

use super::{HeaderMode, Options, Summary, samples, selection};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum OverlapMode {
    Position,
    #[default]
    Record,
    Variant,
}

#[derive(Clone, Debug)]
pub struct RegionSet {
    regions: Vec<Region>,
    overlap: OverlapMode,
}

impl RegionSet {
    pub fn new(regions: Vec<Region>, overlap: OverlapMode) -> Result<Self> {
        if regions.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "region list is empty".to_owned(),
            ));
        }
        Ok(Self { regions, overlap })
    }

    pub fn parse(values: impl IntoIterator<Item = String>, overlap: OverlapMode) -> Result<Self> {
        let regions = values
            .into_iter()
            .map(|value| {
                value.parse().map_err(|error| {
                    RsomicsError::InvalidInput(format!("invalid region {value:?}: {error}"))
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Self::new(regions, overlap)
    }

    pub(super) fn matches(&self, record: &RecordBuf) -> bool {
        self.regions.iter().any(|region| {
            let name: &[u8] = region.name().as_ref();
            name == record.reference_sequence_name().as_bytes()
                && overlaps(record, region.interval(), self.overlap)
        })
    }
}

pub(super) fn write_indexed(
    input: &Path,
    options: &Options,
    regions: &RegionSet,
    output: impl Write,
) -> Result<Summary> {
    let mut reader = indexed_reader::Builder::default()
        .build_from_path(input)
        .rs_with_context(|| format!("opening indexed variant input {}", input.display()))?;
    let header = reader
        .read_header()
        .rs_with_context(|| format!("reading variant header {}", input.display()))?;
    let query_regions = regions.merged(&header)?;
    let projection = samples::Projection::new(&header, options)?;
    let output_header = projection.header(&header, options);
    let mut writer = Writer::new(output, options.output_format);
    writer.write_header(&output_header, options.header)?;

    if options.header == HeaderMode::HeaderOnly {
        writer.finish()?;
        return Ok(Summary {
            read: 0,
            written: 0,
            output_format: options.output_format,
        });
    }

    let mut read = 0;
    let mut written = 0;
    for (region_index, region) in query_regions.iter().enumerate() {
        let records = reader
            .query(&header, region)
            .rs_with_context(|| format!("querying region {region}"))?;
        for result in records {
            let number = read + 1;
            let record =
                result.rs_with_context(|| format!("reading indexed variant record {number}"))?;
            let record = RecordBuf::try_from_variant_record(&header, record.as_ref())
                .map_err(|error| super::invalid(input, number, "decoding indexed input", error))?;
            read += 1;

            if !overlaps(&record, region.interval(), regions.overlap)
                || query_regions[..region_index].iter().any(|previous| {
                    previous.name() == region.name()
                        && overlaps(&record, previous.interval(), regions.overlap)
                })
                || !selection::keep(&record, options)
                || options
                    .targets
                    .as_ref()
                    .is_some_and(|targets| !targets.matches(&record))
            {
                continue;
            }

            let record = projection.apply(record, options)?;
            writer.write_record(&output_header, &record, number)?;
            written += 1;
        }
    }
    writer.finish()?;

    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}

impl RegionSet {
    fn merged(&self, header: &noodles_vcf::Header) -> Result<Vec<Region>> {
        let mut values = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            let name = std::str::from_utf8(region.name().as_ref()).map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "region reference name is not UTF-8: {}",
                    region.name()
                ))
            })?;
            let reference = header.contigs().get_index_of(name).ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "region reference is absent from the VCF header: {name}"
                ))
            })?;
            values.push((
                reference,
                name.to_owned(),
                region.interval().start().unwrap_or(Position::MIN),
                region.interval().end().unwrap_or(Position::MAX),
            ));
        }
        values.sort_by_key(|(reference, _, start, end)| (*reference, *start, *end));

        let mut merged: Vec<(usize, String, Position, Position)> = Vec::new();
        for (reference, name, start, end) in values {
            if let Some((last_reference, _, _, last_end)) = merged.last_mut()
                && *last_reference == reference
                && start <= last_end.checked_add(1).unwrap_or(Position::MAX)
            {
                *last_end = (*last_end).max(end);
                continue;
            }
            merged.push((reference, name, start, end));
        }
        Ok(merged
            .into_iter()
            .map(|(_, name, start, end)| Region::new(name, start..=end))
            .collect())
    }
}

fn overlaps(record: &RecordBuf, interval: Interval, mode: OverlapMode) -> bool {
    let Some(start) = record.variant_start() else {
        return false;
    };
    match mode {
        OverlapMode::Position => interval.contains(start),
        OverlapMode::Record => interval.intersects(record_interval(record, start)),
        OverlapMode::Variant => variant_intervals(record, start)
            .into_iter()
            .any(|variant| interval.intersects(variant)),
    }
}

fn record_interval(record: &RecordBuf, start: Position) -> Interval {
    let end = record
        .info()
        .get(key::END_POSITION)
        .flatten()
        .and_then(|value| match value {
            Value::Integer(value) => usize::try_from(*value).ok(),
            _ => None,
        })
        .and_then(|value| Position::try_from(value).ok())
        .filter(|end| *end >= start)
        .unwrap_or_else(|| {
            start
                .checked_add(record.reference_bases().len().max(1) - 1)
                .unwrap_or(Position::MAX)
        });
    Interval::from(start..=end)
}

fn variant_intervals(record: &RecordBuf, start: Position) -> Vec<Interval> {
    let reference = record.reference_bases().as_bytes();
    let mut intervals = Vec::new();
    for alternate in record.alternate_bases().as_ref() {
        let alternate = alternate.as_bytes();
        if alternate.starts_with(b"<") {
            intervals.push(record_interval(record, start));
            continue;
        }
        if alternate == b"*" || alternate.contains(&b'[') || alternate.contains(&b']') {
            intervals.push(Interval::from(start..=start));
            continue;
        }

        let mut prefix = 0;
        while prefix < reference.len()
            && prefix < alternate.len()
            && reference[prefix].eq_ignore_ascii_case(&alternate[prefix])
        {
            prefix += 1;
        }
        let mut reference_end = reference.len();
        let mut alternate_end = alternate.len();
        while reference_end > prefix
            && alternate_end > prefix
            && reference[reference_end - 1].eq_ignore_ascii_case(&alternate[alternate_end - 1])
        {
            reference_end -= 1;
            alternate_end -= 1;
        }
        if reference_end == prefix && alternate_end == prefix {
            intervals.push(Interval::from(start..=start));
            continue;
        }
        let offset = prefix;
        let variant_start = start.checked_add(offset).unwrap_or(Position::MAX);
        let changed_reference = reference_end.saturating_sub(prefix);
        let variant_end = variant_start
            .checked_add(changed_reference.saturating_sub(1))
            .unwrap_or(variant_start);
        intervals.push(Interval::from(variant_start..=variant_end));
    }
    if intervals.is_empty() {
        intervals.push(Interval::from(start..=start));
    }
    intervals
}
