use noodles_vcf::{
    self as vcf,
    variant::{RecordBuf, record_buf::info::field::Value},
};
use rsomics_common::{Result, RsomicsError};

use super::columns::{BoundColumns, Column, MatchField};
use super::source::{AnnotationRecord, AnnotationSource, Payload};
use crate::variant_type;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum PairLogic {
    Snps,
    Indels,
    Both,
    All,
    #[default]
    Some,
    Exact,
    Id,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct OverlapFractions {
    pub(crate) annotation: f64,
    pub(crate) target: f64,
}

#[derive(Debug, PartialEq)]
pub(crate) struct Matched<'a> {
    pub(crate) source: &'a AnnotationRecord,
    pub(crate) allele_map: Vec<Option<usize>>,
}

#[derive(Clone, Copy)]
struct Constraints {
    alleles: bool,
    id: bool,
    end: bool,
}

impl AnnotationSource {
    pub(crate) fn first_match(
        &mut self,
        target: &RecordBuf,
        logic: PairLogic,
        fractions: OverlapFractions,
    ) -> Result<Option<Matched<'_>>> {
        fractions.validate()?;
        let target_contig = *self
            .contigs
            .get(target.reference_sequence_name())
            .ok_or_else(|| {
                invalid(format!(
                    "target uses unknown contig {:?}",
                    target.reference_sequence_name()
                ))
            })?;
        let (target_start, target_end, target_info_end) = target_span(target)?;
        let coordinate = (target_contig, target_start);
        if self
            .last_target_coordinate
            .is_some_and(|previous| coordinate < previous)
        {
            return Err(invalid("target is out of coordinate order"));
        }
        self.last_target_coordinate = Some(coordinate);

        self.active.retain(|record| {
            record.contig == target_contig && record.inclusive_end() >= target_start
        });
        while self.next.as_ref().is_some_and(|record| {
            record.contig < target_contig
                || (record.contig == target_contig && record.inclusive_start() <= target_end)
        }) {
            let record = self
                .read_record()?
                .expect("lookahead record exists while advancing");
            if record.contig == target_contig && record.inclusive_end() >= target_start {
                self.active.push_back(record);
            }
        }

        let constraints = Constraints::from_columns(self.columns());
        for source in &self.active {
            if !overlaps(source, target_contig, target_start, target_end)
                || !fractions.keeps(source, target_start, target_end)
            {
                continue;
            }
            let allele_map = match &source.payload {
                Payload::Variant(_) => pair_match(source, target, logic)?,
                Payload::Tabular(_) => {
                    tabular_match(source, target, target_start, target_info_end, constraints)?
                }
            };
            if let Some(allele_map) = allele_map {
                return Ok(Some(Matched { source, allele_map }));
            }
        }
        Ok(None)
    }

    pub(crate) fn active_len(&self) -> usize {
        self.active.len()
    }
}

impl OverlapFractions {
    pub(crate) fn validate(self) -> Result<()> {
        if [self.annotation, self.target]
            .into_iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(&value))
        {
            return Err(RsomicsError::ConfigError(
                "annotation overlap fractions must be finite values between 0 and 1".to_owned(),
            ));
        }
        Ok(())
    }

    fn keeps(self, source: &AnnotationRecord, target_start: usize, target_end: usize) -> bool {
        let start = source.inclusive_start().max(target_start);
        let end = source.inclusive_end().min(target_end);
        if start > end {
            return false;
        }
        let intersection = (end - start + 1) as f64;
        let annotation = (source.inclusive_end() - source.inclusive_start() + 1) as f64;
        let target = (target_end - target_start + 1) as f64;
        intersection / annotation >= self.annotation && intersection / target >= self.target
    }
}

impl Constraints {
    fn from_columns(columns: &BoundColumns) -> Self {
        let mut reference = false;
        let mut alternate = false;
        for column in columns.spec().fields() {
            match column {
                Column::Match(MatchField::Ref) => reference = true,
                Column::Match(MatchField::Alt) => alternate = true,
                Column::Match(_) | Column::Transfer(_) => {}
            }
        }
        let tabular = columns.source_kind() == super::columns::SourceKind::Tabular;
        debug_assert!(!tabular || reference == alternate);
        Self {
            alleles: tabular && reference && alternate,
            id: columns.spec().matches_id(),
            end: columns.spec().matches_end(),
        }
    }
}

fn overlaps(
    source: &AnnotationRecord,
    target_contig: usize,
    target_start: usize,
    target_end: usize,
) -> bool {
    source.contig == target_contig
        && source.inclusive_start() <= target_end
        && source.inclusive_end() >= target_start
}

fn tabular_match(
    source: &AnnotationRecord,
    target: &RecordBuf,
    target_start: usize,
    target_info_end: Option<usize>,
    constraints: Constraints,
) -> Result<Option<Vec<Option<usize>>>> {
    if constraints.id && !ids_overlap(source.id.as_deref(), target) {
        return Ok(None);
    }
    if constraints.end && target_info_end.is_some() && source.info_end != target_info_end {
        return Ok(None);
    }
    if !constraints.alleles {
        return Ok(Some(Vec::new()));
    }
    if source.inclusive_start() != target_start {
        return Ok(None);
    }
    let map = allele_map(source, target);
    let matches = if source.alternates.is_empty() {
        target.alternate_bases().as_ref().is_empty() && map.first() == Some(&Some(0))
    } else {
        map.iter().skip(1).any(Option::is_some)
    };
    Ok(matches.then_some(map))
}

fn pair_match(
    source: &AnnotationRecord,
    target: &RecordBuf,
    logic: PairLogic,
) -> Result<Option<Vec<Option<usize>>>> {
    let target_start = target.variant_start().map_or(0, usize::from);
    if source.start != target_start {
        return Ok(None);
    }
    let target_info_end = info_end(target, "target")?;
    let shared = shared_allele(source, target, target_info_end);
    let exact = exact_alleles(source, target, target_info_end);
    let source_record = match &source.payload {
        Payload::Variant(record) => record.as_ref(),
        Payload::Tabular(_) => return Ok(None),
    };
    let source_mask = variant_type::record_mask(source_record);
    let target_mask = variant_type::record_mask(target);
    let compatible = match logic {
        PairLogic::Snps => shared || ref_type_pair(source_mask, target_mask, TypeClass::Snp),
        PairLogic::Indels => shared || ref_type_pair(source_mask, target_mask, TypeClass::Indel),
        PairLogic::Both => shared || both_type_pair(source_mask, target_mask),
        PairLogic::All => true,
        PairLogic::Some => shared,
        PairLogic::Exact => exact,
        PairLogic::Id => exact && identical_ids(source.id.as_deref(), target),
    };
    Ok(compatible.then(|| allele_map(source, target)))
}

fn allele_map(source: &AnnotationRecord, target: &RecordBuf) -> Vec<Option<usize>> {
    let Some(reference) = source.reference.as_deref() else {
        return Vec::new();
    };
    let Some(delta) = RefDelta::between(target.reference_bases(), reference) else {
        return vec![None; source.alternates.len() + 1];
    };
    let mut map = Vec::with_capacity(source.alternates.len() + 1);
    map.push(Some(0));
    for alternate in &source.alternates {
        let index = target
            .alternate_bases()
            .as_ref()
            .iter()
            .position(|target_alternate| compatible_alternate(target_alternate, alternate, delta))
            .map(|index| index + 1);
        map.push(index);
    }
    map
}

fn target_span(record: &RecordBuf) -> Result<(usize, usize, Option<usize>)> {
    let start = record.variant_start().map_or(0, usize::from);
    let info_end = info_end(record, "target")?;
    let end = match info_end {
        Some(end) => end,
        None => start
            .checked_add(record.reference_bases().len().max(1) - 1)
            .ok_or_else(|| invalid("target span exceeds usize"))?,
    };
    if end < start {
        return Err(invalid(format!("target end {end} precedes start {start}")));
    }
    Ok((start, end, info_end))
}

fn info_end(record: &RecordBuf, label: &str) -> Result<Option<usize>> {
    match record
        .info()
        .get(vcf::variant::record::info::field::key::END_POSITION)
    {
        None => Ok(None),
        Some(Some(Value::Integer(value))) => usize::try_from(*value)
            .map(Some)
            .map_err(|_| invalid(format!("{label} has invalid END value {value}"))),
        Some(_) => Err(invalid(format!("{label} has invalid END value"))),
    }
}

fn shared_allele(source: &AnnotationRecord, target: &RecordBuf, target_end: Option<usize>) -> bool {
    if !source
        .reference
        .as_deref()
        .is_some_and(|reference| reference.eq_ignore_ascii_case(target.reference_bases()))
    {
        return false;
    }
    if source.alternates.is_empty() && target.alternate_bases().as_ref().is_empty() {
        return true;
    }
    source.alternates.iter().any(|source_alternate| {
        target
            .alternate_bases()
            .as_ref()
            .iter()
            .any(|target_alternate| {
                source_alternate.eq_ignore_ascii_case(target_alternate)
                    && symbolic_end_matches(source_alternate, source.info_end, target_end)
            })
    })
}

fn exact_alleles(source: &AnnotationRecord, target: &RecordBuf, target_end: Option<usize>) -> bool {
    if source.alternates.len() != target.alternate_bases().as_ref().len()
        || !source
            .reference
            .as_deref()
            .is_some_and(|reference| reference.eq_ignore_ascii_case(target.reference_bases()))
    {
        return false;
    }
    let mut used = vec![false; target.alternate_bases().as_ref().len()];
    source.alternates.iter().all(|source_alternate| {
        let index = target
            .alternate_bases()
            .as_ref()
            .iter()
            .enumerate()
            .position(|(index, target_alternate)| {
                !used[index]
                    && source_alternate.eq_ignore_ascii_case(target_alternate)
                    && symbolic_end_matches(source_alternate, source.info_end, target_end)
            });
        if let Some(index) = index {
            used[index] = true;
            true
        } else {
            false
        }
    })
}

fn symbolic_end_matches(allele: &str, source: Option<usize>, target: Option<usize>) -> bool {
    !allele.starts_with('<') || source == target
}

fn ids_overlap(source: Option<&str>, target: &RecordBuf) -> bool {
    let Some(source) = source else {
        return false;
    };
    source
        .split(';')
        .any(|source| target.ids().as_ref().iter().any(|target| source == target))
}

fn identical_ids(source: Option<&str>, target: &RecordBuf) -> bool {
    match source {
        Some(source) => source
            .split(';')
            .eq(target.ids().as_ref().iter().map(String::as_str)),
        None => target.ids().as_ref().is_empty(),
    }
}

#[derive(Clone, Copy)]
enum TypeClass {
    Snp,
    Indel,
}

fn ref_type_pair(source: u32, target: u32, class: TypeClass) -> bool {
    let bit = match class {
        TypeClass::Snp => variant_type::SNP | variant_type::MNP,
        TypeClass::Indel => variant_type::INDEL,
    };
    (source == variant_type::REF && target & bit != 0)
        || (target == variant_type::REF && source & bit != 0)
}

fn both_type_pair(source: u32, target: u32) -> bool {
    let snp = variant_type::SNP | variant_type::MNP;
    (source & snp != 0 && target & snp != 0)
        || (source & variant_type::INDEL != 0 && target & variant_type::INDEL != 0)
        || ref_type_pair(source, target, TypeClass::Snp)
        || ref_type_pair(source, target, TypeClass::Indel)
}

#[derive(Clone, Copy)]
enum RefDelta<'a> {
    Equal,
    TargetLonger(&'a [u8]),
    SourceLonger(&'a [u8]),
}

impl<'a> RefDelta<'a> {
    fn between(target: &'a str, source: &'a str) -> Option<Self> {
        let target = target.as_bytes();
        let source = source.as_bytes();
        let shared = target
            .iter()
            .zip(source)
            .take_while(|(target, source)| target.eq_ignore_ascii_case(source))
            .count();
        match (shared == target.len(), shared == source.len()) {
            (true, true) => Some(Self::Equal),
            (false, true) => Some(Self::TargetLonger(&target[shared..])),
            (true, false) => Some(Self::SourceLonger(&source[shared..])),
            (false, false) => None,
        }
    }
}

fn compatible_alternate(target: &str, source: &str, delta: RefDelta<'_>) -> bool {
    let target = target.as_bytes();
    let source = source.as_bytes();
    let shared = target
        .iter()
        .zip(source)
        .take_while(|(target, source)| target.eq_ignore_ascii_case(source))
        .count();
    let target = &target[shared..];
    let source = &source[shared..];
    match delta {
        RefDelta::Equal => target.is_empty() && source.is_empty(),
        RefDelta::TargetLonger(expected) => source.is_empty() && equal_bytes(target, expected),
        RefDelta::SourceLonger(expected) => target.is_empty() && equal_bytes(source, expected),
    }
}

fn equal_bytes(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(left, right)| left.eq_ignore_ascii_case(right))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    use std::path::Path;

    use noodles_core::Position;
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;
    use crate::annotate::columns::ColumnSpec;
    use crate::annotate::source::{AnnotationRecord, AnnotationSource, variant_record};

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=2000000>\n\
##contig=<ID=chr2,length=2000000>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    struct TestSource {
        _directory: tempfile::TempDir,
        source: AnnotationSource,
    }

    fn header() -> vcf::Header {
        HEADER.parse().unwrap()
    }

    fn record(line: &str) -> RecordBuf {
        let raw = vcf::Record::try_from(line.as_bytes()).unwrap();
        RecordBuf::try_from_variant_record(&header(), &raw).unwrap()
    }

    fn annotation(line: &str) -> AnnotationRecord {
        variant_record(record(line), 1).unwrap()
    }

    fn source(name: &str, content: &str, columns: &str) -> TestSource {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join(name);
        std::fs::write(&path, content).unwrap();
        let spec = ColumnSpec::parse(columns).unwrap();
        let source = AnnotationSource::open(&path, &header(), spec).unwrap();
        TestSource {
            _directory: directory,
            source,
        }
    }

    fn matches_pair(source: &str, target: &str, logic: PairLogic) -> bool {
        pair_match(&annotation(source), &record(target), logic)
            .unwrap()
            .is_some()
    }

    #[test]
    fn matches_bed_and_tabular_intervals_with_distinct_coordinates() {
        let target = record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let fractions = OverlapFractions::default();

        let mut bed = source("regions.bed", "chr1\t9\t20\tBED\n", "CHROM,FROM,TO,INFO/X");
        let matched = bed
            .source
            .first_match(&target, PairLogic::Some, fractions)
            .unwrap()
            .unwrap();
        assert_eq!(matched.source.serial, 1);
        assert_eq!(matched.source.inclusive_start(), 10);

        let mut tab = source("regions.tsv", "chr1\t10\t20\tTAB\n", "CHROM,FROM,TO,INFO/X");
        assert!(
            tab.source
                .first_match(&target, PairLogic::Some, fractions)
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn target_ref_and_end_spans_match_source_points() {
        let mut input = source(
            "points.tsv",
            "chr1\t13\tREF\nchr1\t20\tEND\n",
            "CHROM,POS,INFO/X",
        );

        let by_ref = record("chr1\t10\t.\tAAAA\tA\t.\tPASS\t.");
        let first = input
            .source
            .first_match(&by_ref, PairLogic::Some, OverlapFractions::default())
            .unwrap()
            .unwrap();
        assert_eq!(first.source.serial, 1);

        let by_end = record("chr1\t14\t.\tA\t<DEL>\t.\tPASS\tEND=20");
        let second = input
            .source
            .first_match(&by_end, PairLogic::Some, OverlapFractions::default())
            .unwrap()
            .unwrap();
        assert_eq!(second.source.serial, 2);
    }

    #[test]
    fn reciprocal_overlap_thresholds_include_exact_boundaries() {
        let target = record("chr1\t15\t.\tAAAAA\tA\t.\tPASS\t.");
        let fractions = OverlapFractions {
            annotation: 0.5,
            target: 1.0,
        };
        let mut exact = source("exact.tsv", "chr1\t10\t19\tX\n", "CHROM,FROM,TO,INFO/X");
        assert!(
            exact
                .source
                .first_match(&target, PairLogic::Some, fractions)
                .unwrap()
                .is_some()
        );

        let mut above = source("above.tsv", "chr1\t10\t19\tX\n", "CHROM,FROM,TO,INFO/X");
        assert!(
            above
                .source
                .first_match(
                    &target,
                    PairLogic::Some,
                    OverlapFractions {
                        annotation: 0.500_001,
                        target: 1.0,
                    },
                )
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn enforces_tabular_ref_alt_id_and_end_constraints() {
        let columns = "CHROM,POS,REF,ALT,~ID,~INFO/END,INFO/X";
        let line = "chr1\t10\tA\tC\trs1\t20\tX\n";

        for (target, expected) in [
            ("chr1\t10\trs1\tA\tC\t.\tPASS\tEND=20", true),
            ("chr1\t10\trs1\tA\tC\t.\tPASS\t.", true),
            ("chr1\t10\trs2\tA\tC\t.\tPASS\tEND=20", false),
            ("chr1\t10\trs1\tA\tG\t.\tPASS\tEND=20", false),
            ("chr1\t10\trs1\tA\tC\t.\tPASS\tEND=21", false),
        ] {
            let mut input = source("constraints.tsv", line, columns);
            let target = record(target);
            assert_eq!(
                input
                    .source
                    .first_match(&target, PairLogic::Some, OverlapFractions::default())
                    .unwrap()
                    .is_some(),
                expected,
                "{target:?}"
            );
        }
    }

    #[test]
    fn tabular_symbolic_alleles_need_end_only_when_requested() {
        let mut input = source(
            "symbolic.tsv",
            "chr1\t10\tN\t<DEL>\tX\n",
            "CHROM,POS,REF,ALT,INFO/X",
        );
        let target = record("chr1\t10\t.\tN\t<DEL>\t.\tPASS\tEND=20");
        assert!(
            input
                .source
                .first_match(&target, PairLogic::Some, OverlapFractions::default())
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn pair_logic_matches_bcftools_1_24_categories() {
        let snp_source = "chr1\t10\tz\tA\tG\t.\tPASS\t.";
        let snp_target = "chr1\t10\ta\tA\tC\t.\tPASS\t.";
        assert!(!matches_pair(snp_source, snp_target, PairLogic::Snps));
        assert!(!matches_pair(snp_source, snp_target, PairLogic::Indels));
        assert!(matches_pair(snp_source, snp_target, PairLogic::Both));
        assert!(matches_pair(snp_source, snp_target, PairLogic::All));
        assert!(!matches_pair(snp_source, snp_target, PairLogic::Some));
        assert!(!matches_pair(snp_source, snp_target, PairLogic::Exact));

        let shared_source = "chr1\t10\tn\tA\tG,T\t.\tPASS\t.";
        let shared_target = "chr1\t10\tm\tA\tC,G\t.\tPASS\t.";
        for logic in [
            PairLogic::Snps,
            PairLogic::Indels,
            PairLogic::Both,
            PairLogic::All,
            PairLogic::Some,
        ] {
            assert!(
                matches_pair(shared_source, shared_target, logic),
                "{logic:?}"
            );
        }
        assert!(!matches_pair(
            shared_source,
            shared_target,
            PairLogic::Exact
        ));

        let reordered = "chr1\t10\tn\tA\tT,G\t.\tPASS\t.";
        assert!(matches_pair(shared_source, reordered, PairLogic::Exact));
        assert!(matches_pair(shared_source, reordered, PairLogic::Id));
        let different_id = "chr1\t10\tx\tA\tT,G\t.\tPASS\t.";
        assert!(!matches_pair(shared_source, different_id, PairLogic::Id));
    }

    #[test]
    fn variant_stream_uses_pair_logic_and_returns_the_allele_map() {
        let content = format!("{HEADER}chr1\t10\tz\tA\tG\t.\tPASS\t.\n");
        let mut input = source("annotation.vcf", &content, "ID");
        let target = record("chr1\t10\ta\tA\tC\t.\tPASS\t.");
        let matched = input
            .source
            .first_match(&target, PairLogic::Both, OverlapFractions::default())
            .unwrap()
            .unwrap();

        assert_eq!(matched.source.serial, 1);
        assert_eq!(matched.allele_map, [Some(0), None]);
    }

    #[test]
    fn pair_logic_handles_symbolic_breakend_overlap_ref_and_mixed_records() {
        let symbolic = "chr1\t10\ts\tN\t<DEL>\t.\tPASS\tEND=20";
        assert!(matches_pair(symbolic, symbolic, PairLogic::Some));
        assert!(!matches_pair(
            symbolic,
            "chr1\t10\ts\tN\t<DEL>\t.\tPASS\tEND=21",
            PairLogic::Some,
        ));
        assert!(matches_pair(
            symbolic,
            "chr1\t10\ts\tN\t<DEL>\t.\tPASS\tEND=21",
            PairLogic::All,
        ));

        let breakend = "chr1\t10\tb\tN\tN]chr2:1]\t.\tPASS\t.";
        assert!(matches_pair(breakend, breakend, PairLogic::Some));
        assert!(!matches_pair(
            breakend,
            "chr1\t10\tb\tN\tN]chr2:2]\t.\tPASS\t.",
            PairLogic::Some,
        ));

        let overlap = "chr1\t10\to\tA\t*\t.\tPASS\t.";
        assert!(matches_pair(overlap, overlap, PairLogic::Some));

        let reference = "chr1\t10\tr\tA\t.\t.\tPASS\t.";
        let snp = "chr1\t10\ts\tA\tC\t.\tPASS\t.";
        let indel = "chr1\t10\ti\tA\tAT\t.\tPASS\t.";
        assert!(matches_pair(reference, snp, PairLogic::Snps));
        assert!(matches_pair(reference, indel, PairLogic::Indels));
        assert!(!matches_pair(reference, snp, PairLogic::Some));

        let mixed_source = "chr1\t10\tm\tA\tG,AT\t.\tPASS\t.";
        let mixed_target = "chr1\t10\tn\tA\tC,AG\t.\tPASS\t.";
        assert!(matches_pair(mixed_source, mixed_target, PairLogic::Both));
    }

    #[test]
    fn allele_maps_follow_source_order_and_reference_extension() {
        let source = annotation("chr1\t10\t.\tA\tG,C\t.\tPASS\t.");
        let target = record("chr1\t10\t.\tA\tC,G\t.\tPASS\t.");
        let map = pair_match(&source, &target, PairLogic::Some)
            .unwrap()
            .unwrap();
        assert_eq!(map, [Some(0), Some(2), Some(1)]);

        let source = annotation("chr1\t10\t.\tA\tG\t.\tPASS\t.");
        let target = record("chr1\t10\t.\tAT\tGT\t.\tPASS\t.");
        let map = allele_map(&source, &target);
        assert_eq!(map, [Some(0), Some(1)]);
    }

    #[test]
    fn repeated_coordinates_preserve_source_order_across_contigs() {
        let mut input = source(
            "ordered.tsv",
            "chr1\t10\tA\tG\tfirst\n\
chr1\t10\tA\tC\tsecond\n\
chr2\t1\tA\tT\tthird\n",
            "CHROM,POS,REF,ALT,INFO/X",
        );

        let target = record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let matched = input
            .source
            .first_match(&target, PairLogic::Some, OverlapFractions::default())
            .unwrap()
            .unwrap();
        assert_eq!(matched.source.serial, 2);

        let target = record("chr2\t1\t.\tA\tT\t.\tPASS\t.");
        let matched = input
            .source
            .first_match(&target, PairLogic::Some, OverlapFractions::default())
            .unwrap()
            .unwrap();
        assert_eq!(matched.source.serial, 3);
    }

    #[test]
    fn rejects_target_coordinate_regressions_and_invalid_fractions() {
        let mut input = source("ordered.tsv", "chr1\t1\t10\tX\n", "CHROM,FROM,TO,INFO/X");
        let second = record("chr1\t2\t.\tA\tC\t.\tPASS\t.");
        input
            .source
            .first_match(&second, PairLogic::Some, OverlapFractions::default())
            .unwrap();
        let first = record("chr1\t1\t.\tA\tC\t.\tPASS\t.");
        let error = input
            .source
            .first_match(&first, PairLogic::Some, OverlapFractions::default())
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("target is out of coordinate order")
        );

        let mut input = source("fractions.tsv", "chr1\t1\t10\tX\n", "CHROM,FROM,TO,INFO/X");
        let error = input
            .source
            .first_match(
                &first,
                PairLogic::Some,
                OverlapFractions {
                    annotation: f64::NAN,
                    target: 0.0,
                },
            )
            .unwrap_err();
        assert!(error.to_string().contains("overlap fraction"));
    }

    #[test]
    fn forward_join_discards_expired_annotations() {
        const RECORDS: usize = 1_000_000;
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("million.tsv");
        write_generated(&path, RECORDS);
        let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/X").unwrap();
        let mut source = AnnotationSource::open(&path, &header(), spec).unwrap();
        let mut target = RecordBuf::builder()
            .set_reference_sequence_name("chr1")
            .set_variant_start(Position::MIN)
            .set_reference_bases("A")
            .build();

        for position in 1..=RECORDS {
            *target.variant_start_mut() = Position::try_from(position).ok();
            assert!(
                source
                    .first_match(&target, PairLogic::Some, OverlapFractions::default())
                    .unwrap()
                    .is_some()
            );
            assert!(
                source.active_len() <= 5,
                "{position}: {}",
                source.active_len()
            );
        }
    }

    fn write_generated(path: &Path, records: usize) {
        let mut writer = BufWriter::with_capacity(1024 * 1024, File::create(path).unwrap());
        for position in 1..=records {
            writeln!(writer, "chr1\t{position}\t{}\tX", position + 4).unwrap();
        }
        writer.flush().unwrap();
    }
}
