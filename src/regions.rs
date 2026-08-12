use std::fs::File;
use std::path::{Path, PathBuf};

use noodles_bgzf as bgzf;
use noodles_core::{Position, Region, region::Interval};
use noodles_util::variant::io::indexed_reader;
use noodles_vcf::variant::{RecordBuf, record::info::field::key, record_buf::info::field::Value};
use rsomics_common::{Context, Result, RsomicsError};

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

#[derive(Clone, Debug)]
pub struct RegionSelection {
    regions: RegionSet,
    exclude: bool,
}

pub(crate) struct IndexedRecords {
    reader: indexed_reader::IndexedReader<bgzf::io::Reader<File>>,
    header: noodles_vcf::Header,
    regions: Vec<Region>,
    overlap: OverlapMode,
    input: PathBuf,
}

impl IndexedRecords {
    pub(crate) fn open(input: &Path, regions: &RegionSet) -> Result<Self> {
        let mut reader = indexed_reader::Builder::default()
            .build_from_path(input)
            .rs_with_context(|| format!("opening indexed variant input {}", input.display()))?;
        let header = reader
            .read_header()
            .rs_with_context(|| format!("reading variant header {}", input.display()))?;
        let query_regions = regions.merged(&header)?;
        Ok(Self {
            reader,
            header,
            regions: query_regions,
            overlap: regions.overlap(),
            input: input.to_path_buf(),
        })
    }

    pub(crate) fn header(&self) -> &noodles_vcf::Header {
        &self.header
    }

    pub(crate) fn visit(
        &mut self,
        mut visit: impl FnMut(&noodles_vcf::Header, RecordBuf, u64) -> Result<()>,
    ) -> Result<u64> {
        let Self {
            reader,
            header,
            regions,
            overlap,
            input,
        } = self;
        let mut read = 0;
        for (region_index, region) in regions.iter().enumerate() {
            let records = reader
                .query(header, region)
                .rs_with_context(|| format!("querying region {region}"))?;
            for result in records {
                let number = read + 1;
                let record = result
                    .rs_with_context(|| format!("reading indexed variant record {number}"))?;
                let record = RecordBuf::try_from_variant_record(header, record.as_ref()).map_err(
                    |error| {
                        RsomicsError::InvalidInput(format!(
                            "{}: decoding indexed variant record {number}: {error}",
                            input.display()
                        ))
                    },
                )?;
                read += 1;
                if !overlaps(&record, region.interval(), *overlap)
                    || regions[..region_index].iter().any(|previous| {
                        previous.name() == region.name()
                            && overlaps(&record, previous.interval(), *overlap)
                    })
                {
                    continue;
                }
                visit(header, record, number)?;
            }
        }
        Ok(read)
    }
}

impl RegionSelection {
    pub fn new(regions: RegionSet, exclude: bool) -> Self {
        Self { regions, exclude }
    }

    pub fn parse(
        values: impl IntoIterator<Item = String>,
        overlap: OverlapMode,
        exclude: bool,
    ) -> Result<Self> {
        RegionSet::parse(values, overlap).map(|regions| Self::new(regions, exclude))
    }

    pub(crate) fn keeps(&self, record: &RecordBuf) -> bool {
        self.regions.matches(record) != self.exclude
    }
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
            .map(|value| parse_region(&value))
            .collect::<Result<Vec<_>>>()?;
        Self::new(regions, overlap)
    }

    pub(crate) fn matches(&self, record: &RecordBuf) -> bool {
        self.regions.iter().any(|region| {
            let name: &[u8] = region.name().as_ref();
            name == record.reference_sequence_name().as_bytes()
                && overlaps(record, region.interval(), self.overlap)
        })
    }

    pub(crate) fn merged(&self, header: &noodles_vcf::Header) -> Result<Vec<Region>> {
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

    pub(crate) fn overlap(&self) -> OverlapMode {
        self.overlap
    }
}

fn parse_region(value: &str) -> Result<Region> {
    let interval = value.rsplit_once(':').map(|(_, interval)| interval);
    let open_end = interval
        .and_then(|interval| interval.strip_suffix('-'))
        .is_some_and(is_coordinate);
    let normalized = if open_end {
        value.strip_suffix('-').unwrap()
    } else {
        value
    };
    let region: Region = normalized.parse().map_err(|error| {
        RsomicsError::InvalidInput(format!("invalid region {value:?}: {error}"))
    })?;
    if interval.is_some_and(is_coordinate) {
        let start = region.interval().start().ok_or_else(|| {
            RsomicsError::InvalidInput(format!("invalid region {value:?}: missing position"))
        })?;
        Ok(Region::new(region.name().to_owned(), start..=start))
    } else {
        Ok(region)
    }
}

fn is_coordinate(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || byte == b',')
}

pub(crate) fn overlaps(record: &RecordBuf, interval: Interval, mode: OverlapMode) -> bool {
    let Some(start) = record.variant_start() else {
        return false;
    };
    match mode {
        OverlapMode::Position => interval.contains(start),
        OverlapMode::Record => interval.intersects(record_interval(record, start)),
        OverlapMode::Variant => {
            variant_interval(record, start).is_some_and(|variant| interval.intersects(variant))
        }
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

fn variant_interval(record: &RecordBuf, start: Position) -> Option<Interval> {
    let alternates = record.alternate_bases();
    if alternates.as_ref().is_empty() {
        return None;
    }
    let reference = record.reference_bases().as_bytes();
    let end = record_interval(record, start).end().unwrap();
    let mut offset = usize::from(end) - usize::from(start) + 1;
    for alternate in alternates.as_ref() {
        let prefix = reference
            .iter()
            .zip(alternate.as_bytes())
            .take_while(|(reference, alternate)| reference == alternate)
            .count();
        offset = offset.min(prefix);
        if offset == 0 {
            break;
        }
    }
    let variant_start = start.checked_add(offset)?;
    (variant_start <= end).then(|| Interval::from(variant_start..=end))
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;

    fn record(value: &[u8]) -> RecordBuf {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let record = vcf::Record::try_from(value).unwrap();
        RecordBuf::try_from_variant_record(&header, &record).unwrap()
    }

    #[test]
    fn literal_variant_overlap_keeps_the_reference_suffix() {
        let record = record(b"chr1\t10\t.\tAACGT\tAATGT\t10\tPASS\t.");
        assert!(!overlaps(
            &record,
            "12-14".parse().unwrap(),
            OverlapMode::Position
        ));
        assert!(overlaps(
            &record,
            "14".parse().unwrap(),
            OverlapMode::Variant
        ));
    }

    #[test]
    fn literal_insertions_and_absent_alternates_have_no_reference_span() {
        for record in [
            record(b"chr1\t25\t.\tA\tAT\t10\tPASS\t."),
            record(b"chr1\t25\t.\tA\t.\t10\tPASS\t."),
        ] {
            assert!(!overlaps(
                &record,
                "25-26".parse().unwrap(),
                OverlapMode::Variant
            ));
        }
    }

    #[test]
    fn region_selection_can_include_or_exclude_matches() {
        let record = record(b"chr1\t10\t.\tA\tC\t10\tPASS\t.");
        let regions = RegionSet::parse(["chr1:10".to_owned()], OverlapMode::Position).unwrap();

        assert!(RegionSelection::new(regions.clone(), false).keeps(&record));
        assert!(!RegionSelection::new(regions, true).keeps(&record));
    }

    #[test]
    fn single_coordinate_regions_are_not_open_ended() {
        let exact = RegionSet::parse(["chr1:10".to_owned()], OverlapMode::Position).unwrap();
        let open = RegionSet::parse(["chr1:10-".to_owned()], OverlapMode::Position).unwrap();
        let later = record(b"chr1\t12\t.\tA\tC\t10\tPASS\t.");

        assert!(!exact.matches(&later));
        assert!(open.matches(&later));
    }
}
