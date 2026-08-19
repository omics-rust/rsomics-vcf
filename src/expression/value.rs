use std::fmt;

use noodles_vcf::{
    self as vcf,
    variant::{
        RecordBuf,
        record::samples::series::value::genotype::Phasing,
        record_buf::{
            info::field::{Value as InfoValue, value::Array as InfoArray},
            samples::sample::{Value as SampleValue, value::Array as SampleArray},
        },
    },
};

use super::bind::{BoundField, CalculatedField, FieldKind, FixedField};

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Values<'a> {
    Site(Vec<Atom<'a>>),
    Samples(SampleValues<'a>),
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SampleValues<'a> {
    pub values: Vec<Vec<Atom<'a>>>,
    pub selected: Box<[bool]>,
}

impl<'a> SampleValues<'a> {
    pub(crate) fn all(values: Vec<Vec<Atom<'a>>>) -> Self {
        Self {
            selected: vec![true; values.len()].into_boxed_slice(),
            values,
        }
    }
}

pub(crate) fn genotype_indices(record: &RecordBuf) -> Result<Vec<Vec<usize>>, ValueError> {
    let sample_count = record.samples().values().count();
    let Some(genotypes) = record.samples().select("GT") else {
        return Ok(vec![Vec::new(); sample_count]);
    };
    (0..sample_count)
        .map(|index| match genotypes.get(index).flatten() {
            Some(SampleValue::Genotype(genotype)) => {
                let mut indices = Vec::new();
                for allele in genotype.as_ref() {
                    if let Some(position) = allele.position()
                        && !indices.contains(&position)
                    {
                        indices.push(position);
                    }
                }
                indices.sort_unstable();
                Ok(indices)
            }
            Some(_) => Err(ValueError::new("FORMAT/GT is not a genotype")),
            None => Ok(Vec::new()),
        })
        .collect()
}

pub(crate) fn diploid_genotype_indices(
    record: &RecordBuf,
) -> Result<Vec<Option<[usize; 2]>>, ValueError> {
    let sample_count = record.samples().values().count();
    let Some(genotypes) = record.samples().select("GT") else {
        return Ok(vec![None; sample_count]);
    };
    (0..sample_count)
        .map(|index| match genotypes.get(index).flatten() {
            Some(SampleValue::Genotype(genotype)) => {
                let mut alleles = genotype.as_ref().iter();
                let first = alleles.next().and_then(|allele| allele.position());
                let second = alleles.next().and_then(|allele| allele.position());
                Ok(first.zip(second).map(|(first, second)| [first, second]))
            }
            Some(_) => Err(ValueError::new("FORMAT/GT is not a genotype")),
            None => Ok(None),
        })
        .collect()
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Atom<'a> {
    Absent,
    Missing,
    Number(f64),
    Flag,
    Text(&'a str),
    OwnedText(String),
    Genotype(Genotype),
    Filter(Filter<'a>),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Genotype(Vec<Allele>);

impl Genotype {
    pub(super) fn spelling(&self) -> String {
        let mut output = String::new();
        for (index, allele) in self.0.iter().enumerate() {
            if index > 0 {
                output.push(match allele.phasing {
                    Phasing::Phased => '|',
                    Phasing::Unphased => '/',
                });
            }
            match allele.position {
                Some(position) => output.push_str(&position.to_string()),
                None => output.push('.'),
            }
        }
        output
    }

    pub(super) fn matches_class(&self, pattern: &str) -> Option<bool> {
        let complete = self.0.iter().all(|allele| allele.position.is_some());
        let positions: Vec<_> = self.0.iter().filter_map(|allele| allele.position).collect();
        let ploidy = self.0.len();
        let reference = complete && positions.iter().all(|position| *position == 0);
        let alternate = complete && positions.iter().any(|position| *position > 0);
        let homogeneous = complete
            && ploidy > 1
            && positions
                .first()
                .is_some_and(|first| positions.iter().all(|position| position == first));
        let heterogeneous = complete && ploidy > 1 && !homogeneous;

        if pattern.eq_ignore_ascii_case("mis") {
            return Some(!complete);
        }
        if pattern.eq_ignore_ascii_case("ref") {
            return Some(reference);
        }
        if pattern.eq_ignore_ascii_case("alt") {
            return Some(alternate);
        }
        if pattern.eq_ignore_ascii_case("hom") {
            return Some(homogeneous);
        }
        if pattern.eq_ignore_ascii_case("het") {
            return Some(heterogeneous);
        }
        if pattern.eq_ignore_ascii_case("hap") {
            return Some(complete && ploidy == 1);
        }

        match pattern {
            "R" | "r" => Some(complete && ploidy == 1 && reference),
            "A" | "a" => Some(complete && ploidy == 1 && alternate),
            _ if pattern.len() == 2 && pattern.bytes().all(|byte| matches!(byte, b'R' | b'r')) => {
                Some(reference && ploidy > 1)
            }
            _ if pattern.len() == 2
                && pattern
                    .bytes()
                    .all(|byte| matches!(byte, b'R' | b'r' | b'A' | b'a'))
                && pattern.bytes().any(|byte| matches!(byte, b'R' | b'r'))
                && pattern.bytes().any(|byte| matches!(byte, b'A' | b'a')) =>
            {
                Some(complete && ploidy > 1 && reference_alleles(&positions) && alternate)
            }
            _ if pattern.len() == 2 && pattern.bytes().all(|byte| matches!(byte, b'A' | b'a')) => {
                let bytes = pattern.as_bytes();
                let same_case = bytes[0].is_ascii_uppercase() == bytes[1].is_ascii_uppercase();
                Some(
                    complete
                        && ploidy > 1
                        && positions.iter().all(|position| *position > 0)
                        && if same_case {
                            homogeneous
                        } else {
                            heterogeneous
                        },
                )
            }
            _ => None,
        }
    }
}

fn reference_alleles(positions: &[usize]) -> bool {
    positions.contains(&0)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Allele {
    pub position: Option<usize>,
    pub phasing: Phasing,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Filter<'a> {
    Pass,
    Missing,
    Failed(Vec<&'a str>),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ValueError(String);

impl ValueError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for ValueError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for ValueError {}

pub(crate) fn read<'a>(
    field: &BoundField,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Values<'a>, ValueError> {
    match &field.kind {
        FieldKind::Fixed(fixed) => read_fixed(*fixed, record),
        FieldKind::Info(name) => read_info(name, record),
        FieldKind::Format(name) => read_format(name, record),
        FieldKind::Calculated(calculated) => read_calculated(*calculated, header, record),
    }
}

fn read_fixed<'a>(field: FixedField, record: &'a RecordBuf) -> Result<Values<'a>, ValueError> {
    let values = match field {
        FixedField::Chrom => vec![Atom::Text(record.reference_sequence_name())],
        FixedField::Position => vec![
            record
                .variant_start()
                .map(|position| Atom::Number(usize::from(position) as f64))
                .unwrap_or(Atom::Missing),
        ],
        FixedField::Id => {
            let values: Vec<_> = record
                .ids()
                .as_ref()
                .iter()
                .map(|id| Atom::Text(id))
                .collect();
            if values.is_empty() {
                vec![Atom::Missing]
            } else {
                values
            }
        }
        FixedField::Reference => vec![Atom::Text(record.reference_bases())],
        FixedField::Alternate => {
            let values: Vec<_> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|alternate| Atom::Text(alternate))
                .collect();
            if values.is_empty() {
                vec![Atom::Missing]
            } else {
                values
            }
        }
        FixedField::Quality => vec![
            record
                .quality_score()
                .map(|quality| Atom::Number(f64::from(quality)))
                .unwrap_or(Atom::Missing),
        ],
        FixedField::Filter => {
            let filter = if record.filters().is_pass() {
                Filter::Pass
            } else {
                let values: Vec<_> = record
                    .filters()
                    .as_ref()
                    .iter()
                    .map(String::as_str)
                    .collect();
                if values.is_empty() {
                    Filter::Missing
                } else {
                    Filter::Failed(values)
                }
            };
            vec![Atom::Filter(filter)]
        }
        FixedField::Type => {
            let mut output = Vec::new();
            let alternates = record.alternate_bases().as_ref().join(",");
            crate::variant_type::write(
                &mut output,
                record.reference_bases().as_bytes(),
                alternates.as_bytes(),
            );
            vec![Atom::OwnedText(
                String::from_utf8(output)
                    .map_err(|_| ValueError::new("variant type is not UTF-8"))?,
            )]
        }
    };
    Ok(Values::Site(values))
}

fn read_info<'a>(name: &str, record: &'a RecordBuf) -> Result<Values<'a>, ValueError> {
    let values = match record.info().as_ref().get(name) {
        Some(Some(value)) => info_atoms(value),
        Some(None) | None => vec![Atom::Absent],
    };
    Ok(Values::Site(values))
}

fn info_atoms(value: &InfoValue) -> Vec<Atom<'_>> {
    match value {
        InfoValue::Integer(value) => vec![Atom::Number(f64::from(*value))],
        InfoValue::Float(value) => vec![Atom::Number(f64::from(*value))],
        InfoValue::Flag => vec![Atom::Flag],
        InfoValue::Character(value) => vec![Atom::OwnedText(value.to_string())],
        InfoValue::String(value) => vec![Atom::Text(value)],
        InfoValue::Array(values) => match values {
            InfoArray::Integer(values) => optional_numbers(values, |value| f64::from(*value)),
            InfoArray::Float(values) => optional_numbers(values, |value| f64::from(*value)),
            InfoArray::Character(values) => values
                .iter()
                .map(|value| {
                    value
                        .map(|value| Atom::OwnedText(value.to_string()))
                        .unwrap_or(Atom::Missing)
                })
                .collect(),
            InfoArray::String(values) => values
                .iter()
                .map(|value| value.as_deref().map(Atom::Text).unwrap_or(Atom::Missing))
                .collect(),
        },
    }
}

fn read_format<'a>(name: &str, record: &'a RecordBuf) -> Result<Values<'a>, ValueError> {
    if record
        .samples()
        .keys()
        .as_ref()
        .get_index_of(name)
        .is_none()
    {
        return Ok(Values::Samples(SampleValues::all(
            record
                .samples()
                .values()
                .map(|_| vec![Atom::Missing])
                .collect(),
        )));
    }
    let samples = record
        .samples()
        .values()
        .map(|sample| match sample.get(name).flatten() {
            Some(value) => sample_atoms(value),
            None if name == "GT" => Ok(vec![Atom::Genotype(Genotype(vec![Allele {
                position: None,
                phasing: Phasing::Phased,
            }]))]),
            None => Ok(vec![Atom::Missing]),
        })
        .collect::<Result<Vec<_>, ValueError>>()?;
    Ok(Values::Samples(SampleValues::all(samples)))
}

fn sample_atoms(value: &SampleValue) -> Result<Vec<Atom<'_>>, ValueError> {
    match value {
        SampleValue::Integer(value) => Ok(vec![Atom::Number(f64::from(*value))]),
        SampleValue::Float(value) => Ok(vec![Atom::Number(f64::from(*value))]),
        SampleValue::Character(value) => Ok(vec![Atom::OwnedText(value.to_string())]),
        SampleValue::String(value) => Ok(vec![Atom::Text(value)]),
        SampleValue::Genotype(value) => Ok(vec![Atom::Genotype(Genotype(
            value
                .as_ref()
                .iter()
                .map(|allele| Allele {
                    position: allele.position(),
                    phasing: allele.phasing(),
                })
                .collect(),
        ))]),
        SampleValue::Array(values) => Ok(match values {
            SampleArray::Integer(values) => optional_numbers(values, |value| f64::from(*value)),
            SampleArray::Float(values) => optional_numbers(values, |value| f64::from(*value)),
            SampleArray::Character(values) => values
                .iter()
                .map(|value| {
                    value
                        .map(|value| Atom::OwnedText(value.to_string()))
                        .unwrap_or(Atom::Missing)
                })
                .collect(),
            SampleArray::String(values) => values
                .iter()
                .map(|value| value.as_deref().map(Atom::Text).unwrap_or(Atom::Missing))
                .collect(),
        }),
    }
}

fn optional_numbers<'a, T>(
    values: &'a [Option<T>],
    convert: impl Fn(&'a T) -> f64,
) -> Vec<Atom<'a>> {
    values
        .iter()
        .map(|value| {
            value
                .as_ref()
                .map(|value| Atom::Number(convert(value)))
                .unwrap_or(Atom::Missing)
        })
        .collect()
}

fn read_calculated<'a>(
    field: CalculatedField,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Values<'a>, ValueError> {
    let values = match field {
        CalculatedField::AlternateCount => {
            vec![Atom::Number(record.alternate_bases().as_ref().len() as f64)]
        }
        CalculatedField::SampleCount => vec![Atom::Number(header.sample_names().len() as f64)],
        CalculatedField::AlleleCount => genotype_metrics(record)?.alternate_counts(),
        CalculatedField::MinorAlleleCount => {
            vec![
                genotype_metrics(record)?
                    .minor_count()
                    .map(|count| Atom::Number(count as f64))
                    .unwrap_or(Atom::Missing),
            ]
        }
        CalculatedField::AlleleFrequency => genotype_metrics(record)?.alternate_frequencies(),
        CalculatedField::MinorAlleleFrequency => {
            let metrics = genotype_metrics(record)?;
            vec![
                metrics
                    .minor_count()
                    .and_then(|count| metrics.frequency(count))
                    .map(Atom::Number)
                    .unwrap_or(Atom::Missing),
            ]
        }
        CalculatedField::TotalAlleleCount => {
            vec![Atom::Number(genotype_metrics(record)?.total as f64)]
        }
        CalculatedField::MissingSampleCount => {
            vec![
                missing_sample_metrics(record)?
                    .map(|(missing, _)| Atom::Number(missing as f64))
                    .unwrap_or(Atom::Missing),
            ]
        }
        CalculatedField::MissingSampleFraction => {
            let value = missing_sample_metrics(record)?.map(|(missing, samples)| {
                Atom::Number(if samples == 0 {
                    0.0
                } else {
                    missing as f64 / samples as f64
                })
            });
            vec![value.unwrap_or(Atom::Missing)]
        }
        CalculatedField::IndelLength => {
            let reference = record.reference_bases().len();
            let values: Vec<_> = record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(|alternate| {
                    if alternate.starts_with('<') {
                        Atom::Missing
                    } else {
                        Atom::Number(alternate.len() as f64 - reference as f64)
                    }
                })
                .collect();
            if values.is_empty() {
                vec![Atom::Missing]
            } else {
                values
            }
        }
    };
    Ok(Values::Site(values))
}

struct GenotypeMetrics {
    counts: Vec<u64>,
    total: u64,
}

impl GenotypeMetrics {
    fn alternate_counts(&self) -> Vec<Atom<'static>> {
        if self.total == 0 {
            return vec![Atom::Missing];
        }
        if self.counts.len() == 1 {
            return vec![Atom::Number(0.0)];
        }
        self.counts
            .iter()
            .skip(1)
            .map(|count| Atom::Number(*count as f64))
            .collect()
    }

    fn alternate_frequencies(&self) -> Vec<Atom<'static>> {
        if self.total == 0 {
            return vec![Atom::Missing];
        }
        if self.counts.len() == 1 {
            return vec![Atom::Number(0.0)];
        }
        self.counts
            .iter()
            .skip(1)
            .map(|count| {
                self.frequency(*count)
                    .map(Atom::Number)
                    .unwrap_or(Atom::Missing)
            })
            .collect()
    }

    fn minor_count(&self) -> Option<u64> {
        let count = self
            .counts
            .iter()
            .copied()
            .filter(|count| *count > 0)
            .min()?;
        Some(if count == self.total { 0 } else { count })
    }

    fn frequency(&self, count: u64) -> Option<f64> {
        (self.total > 0).then(|| count as f64 / self.total as f64)
    }
}

fn genotype_metrics(record: &RecordBuf) -> Result<GenotypeMetrics, ValueError> {
    let alleles = record.alternate_bases().as_ref().len() + 1;
    let mut metrics = GenotypeMetrics {
        counts: vec![0; alleles],
        total: 0,
    };
    let Some(genotypes) = record.samples().select("GT") else {
        return Ok(metrics);
    };
    for index in 0..record.samples().values().count() {
        let Some(Some(SampleValue::Genotype(genotype))) = genotypes.get(index) else {
            continue;
        };
        for allele in genotype.as_ref() {
            let Some(position) = allele.position() else {
                continue;
            };
            let count = metrics.counts.get_mut(position).ok_or_else(|| {
                ValueError::new(format!(
                    "genotype allele index {position} exceeds {} ALT alleles",
                    alleles - 1
                ))
            })?;
            *count = count
                .checked_add(1)
                .ok_or_else(|| ValueError::new("allele count overflow"))?;
            metrics.total = metrics
                .total
                .checked_add(1)
                .ok_or_else(|| ValueError::new("total allele count overflow"))?;
        }
    }
    Ok(metrics)
}

fn missing_sample_metrics(record: &RecordBuf) -> Result<Option<(u64, u64)>, ValueError> {
    let samples = record.samples().values().count() as u64;
    if samples == 0 {
        return Ok(Some((0, 0)));
    }
    let Some(genotypes) = record.samples().select("GT") else {
        return Ok(None);
    };
    let mut missing = 0u64;
    for index in 0..samples as usize {
        let sample_missing = match genotypes.get(index).flatten() {
            Some(SampleValue::Genotype(genotype)) => genotype
                .as_ref()
                .iter()
                .any(|allele| allele.position().is_none()),
            Some(_) => return Err(ValueError::new("FORMAT/GT is not a genotype")),
            None => true,
        };
        if sample_missing {
            missing = missing
                .checked_add(1)
                .ok_or_else(|| ValueError::new("missing sample count overflow"))?;
        }
    }
    Ok(Some((missing, samples)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expression::{bind::bind, syntax::parse};

    fn parse_record(header: &vcf::Header, line: &[u8]) -> RecordBuf {
        let raw = vcf::Record::try_from(line).unwrap();
        RecordBuf::try_from_variant_record(header, &raw).unwrap()
    }

    fn record() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1>\n\
##FILTER=<ID=LowQual,Description=\"low\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"frequency\">\n\
##INFO=<ID=FLAG,Number=0,Type=Flag,Description=\"flag\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            .parse()
            .unwrap();
        let line = b"chr1\t7\trs1;rs2\tA\tC,G\t50\tLowQual\tDP=12;AF=0.5,.;FLAG\tGT:DP:AD\t0/1:8:4,4,0\t1/2:.:0,3,5\t0/.:2:2,0,.";
        let record = parse_record(&header, line);
        (header, record)
    }

    fn values<'a>(source: &str, header: &vcf::Header, record: &'a RecordBuf) -> Values<'a> {
        let expression = bind(parse(source).unwrap(), header).unwrap();
        let super::super::bind::BoundExpression::Value(super::super::bind::BoundValue::Field(
            field,
        )) = expression
        else {
            panic!("expected field");
        };
        read(&field, header, record).unwrap()
    }

    #[test]
    fn fixed_and_info_fields_preserve_typed_vectors_and_missing_values() {
        let (header, record) = record();
        assert_eq!(
            values("POS", &header, &record),
            Values::Site(vec![Atom::Number(7.0)])
        );
        assert_eq!(
            values("ID", &header, &record),
            Values::Site(vec![Atom::Text("rs1"), Atom::Text("rs2")])
        );
        assert_eq!(
            values("INFO/AF", &header, &record),
            Values::Site(vec![Atom::Number(0.5), Atom::Missing])
        );
        assert_eq!(
            values("INFO/FLAG", &header, &record),
            Values::Site(vec![Atom::Flag])
        );
        assert_eq!(
            values("FILTER", &header, &record),
            Values::Site(vec![Atom::Filter(Filter::Failed(vec!["LowQual"]))])
        );
    }

    #[test]
    fn format_fields_retain_sample_and_value_dimensions() {
        let (header, record) = record();
        assert_eq!(
            values("FMT/DP", &header, &record),
            Values::Samples(SampleValues::all(vec![
                vec![Atom::Number(8.0)],
                vec![Atom::Missing],
                vec![Atom::Number(2.0)],
            ]))
        );
        assert_eq!(
            values("FMT/AD", &header, &record),
            Values::Samples(SampleValues::all(vec![
                vec![Atom::Number(4.0), Atom::Number(4.0), Atom::Number(0.0)],
                vec![Atom::Number(0.0), Atom::Number(3.0), Atom::Number(5.0)],
                vec![Atom::Number(2.0), Atom::Number(0.0), Atom::Missing],
            ]))
        );
    }

    #[test]
    fn genotype_calculations_match_bcftools_1_24_semantics() {
        let (header, record) = record();
        assert_eq!(
            values("AC", &header, &record),
            Values::Site(vec![Atom::Number(2.0), Atom::Number(1.0)])
        );
        assert_eq!(
            values("AN", &header, &record),
            Values::Site(vec![Atom::Number(5.0)])
        );
        assert_eq!(
            values("MAC", &header, &record),
            Values::Site(vec![Atom::Number(1.0)])
        );
        assert_eq!(
            values("N_MISSING", &header, &record),
            Values::Site(vec![Atom::Number(1.0)])
        );
        assert_eq!(
            values("F_MISSING", &header, &record),
            Values::Site(vec![Atom::Number(1.0 / 3.0)])
        );
    }

    #[test]
    fn calculated_fields_preserve_empty_monomorphic_and_symbolic_semantics() {
        let (header, _) = record();
        let monomorphic = parse_record(&header, b"chr1\t8\t.\tA\tC\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1");
        assert_eq!(
            values("MAC", &header, &monomorphic),
            Values::Site(vec![Atom::Number(0.0)])
        );
        assert_eq!(
            values("MAF", &header, &monomorphic),
            Values::Site(vec![Atom::Number(0.0)])
        );

        let missing = parse_record(&header, b"chr1\t9\t.\tA\tC\t.\tPASS\t.\tGT\t./.\t./.\t./.");
        assert_eq!(
            values("AC", &header, &missing),
            Values::Site(vec![Atom::Missing])
        );
        assert_eq!(
            values("MAC", &header, &missing),
            Values::Site(vec![Atom::Missing])
        );

        let symbolic = parse_record(
            &header,
            b"chr1\t10\t.\tAC\t<DEL>,ACT\t.\tPASS\t.\tGT\t0/1\t0/2\t0/0",
        );
        assert_eq!(
            values("ILEN", &header, &symbolic),
            Values::Site(vec![Atom::Missing, Atom::Number(1.0)])
        );

        let reference_only =
            parse_record(&header, b"chr1\t12\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0");
        assert_eq!(
            values("AC", &header, &reference_only),
            Values::Site(vec![Atom::Number(0.0)])
        );
        assert_eq!(
            genotype_metrics(&reference_only)
                .unwrap()
                .alternate_frequencies(),
            vec![Atom::Number(0.0)]
        );
        assert_eq!(
            values("ILEN", &header, &reference_only),
            Values::Site(vec![Atom::Missing])
        );
    }

    #[test]
    fn missing_sample_metrics_require_a_genotype_column() {
        let (header, _) = record();
        let without_gt = parse_record(&header, b"chr1\t11\t.\tA\tC\t.\tPASS\t.\tDP\t1\t.\t3");
        assert_eq!(
            values("N_MISSING", &header, &without_gt),
            Values::Site(vec![Atom::Missing])
        );
        assert_eq!(
            values("F_MISSING", &header, &without_gt),
            Values::Site(vec![Atom::Missing])
        );
    }
}
