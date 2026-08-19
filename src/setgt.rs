mod random;
mod replacement;
mod target;

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format,
    variant::{
        RecordBuf,
        record_buf::samples::sample::{
            Value as SampleValue,
            value::{Array as SampleArray, Genotype},
        },
    },
};
use rsomics_common::{Context, Result, RsomicsError};

use crate::{
    expression::{Compiled, binomial_two_sided},
    genotype::{Change, InfoPolicy, MissingPolicy, edit_selected, reconcile_ac_an},
};
use random::Random48;
use replacement::Replacement;
use target::{Comparison, Principal, Target};

pub(crate) struct Program {
    target: Target,
    replacement: Replacement,
    expression: Option<Compiled>,
    exclude: bool,
    random: Option<Random48>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Query {
    Include(String),
    Exclude(String),
}

impl Program {
    pub(crate) fn bind(
        header: &vcf::Header,
        target: Target,
        replacement: Replacement,
        query: Option<Query>,
        seed: i64,
    ) -> Result<Self> {
        if !header.sample_names().is_empty() {
            let schema = header.formats().get("GT").ok_or_else(|| {
                config("setGT requires FORMAT/GT in a header that declares samples")
            })?;
            if schema.number() != format::Number::Count(1) || schema.ty() != format::Type::String {
                return Err(config("setGT requires FORMAT/GT Number=1,Type=String"));
            }
        }
        let is_query = matches!(target.principal, Principal::Query);
        if is_query != query.is_some() {
            return Err(config(if is_query {
                "setGT target q requires exactly one include or exclude expression"
            } else {
                "setGT expressions require target q"
            }));
        }
        if let Principal::Binomial(binomial) = &target.principal {
            let schema = header.formats().get(&binomial.tag).ok_or_else(|| {
                config(format!(
                    "setGT binomial target requires FORMAT/{} in the header",
                    binomial.tag
                ))
            })?;
            if schema.number() != format::Number::ReferenceAlternateBases
                || schema.ty() != format::Type::Integer
            {
                return Err(config(format!(
                    "setGT binomial target requires FORMAT/{} Number=R,Type=Integer",
                    binomial.tag
                )));
            }
        }
        replacement.validate(header)?;

        let (expression, exclude) = match query {
            Some(Query::Include(source)) => (Some(compile(&source, header)?), false),
            Some(Query::Exclude(source)) => (Some(compile(&source, header)?), true),
            None => (None, false),
        };
        let random = target.random_fraction.map(|_| Random48::new(seed));
        Ok(Self {
            target,
            replacement,
            expression,
            exclude,
            random,
        })
    }

    fn select(
        &mut self,
        header: &vcf::Header,
        record: &RecordBuf,
        number: u64,
    ) -> Result<Vec<bool>> {
        let genotypes = read_genotypes(record, number)?;
        let mut selected = match &self.target.principal {
            Principal::AnyMissing => genotypes
                .iter()
                .map(|genotype| {
                    genotype
                        .as_ref()
                        .iter()
                        .any(|allele| allele.position().is_none())
                })
                .collect(),
            Principal::PartialMissing => genotypes
                .iter()
                .map(|genotype| {
                    let mut missing = false;
                    let mut called = false;
                    for allele in genotype.as_ref() {
                        missing |= allele.position().is_none();
                        called |= allele.position().is_some();
                    }
                    missing && called
                })
                .collect(),
            Principal::CompleteMissing => genotypes
                .iter()
                .map(|genotype| {
                    !genotype.as_ref().is_empty()
                        && genotype
                            .as_ref()
                            .iter()
                            .all(|allele| allele.position().is_none())
                })
                .collect(),
            Principal::All => vec![true; genotypes.len()],
            Principal::Query => self.select_query(header, record, genotypes.len(), number)?,
            Principal::Binomial(binomial) => select_binomial(record, &genotypes, binomial, number)?,
        };

        if let Some(fraction) = self.target.random_fraction {
            let random = self
                .random
                .as_mut()
                .ok_or_else(|| invalid(number, "random state is not initialized"))?;
            for value in &mut selected {
                if *value {
                    *value = random.next() <= fraction;
                }
            }
        }
        Ok(selected)
    }

    pub(crate) fn apply(
        &mut self,
        header: &vcf::Header,
        record: &mut RecordBuf,
        number: u64,
    ) -> Result<Change> {
        let selected = self.select(header, record, number)?;
        if !selected.iter().any(|selected| *selected) {
            return Ok(Change::default());
        }
        let prepared = self
            .replacement
            .prepare(record, &selected)
            .rs_with_context(|| format!("record {number}"))?;
        let change = edit_selected(
            record,
            &selected,
            MissingPolicy::Error,
            |sample, genotype| prepared.resolve(sample, genotype),
        )
        .rs_with_context(|| format!("record {number}"))?;
        if change.genotypes > 0 {
            reconcile_ac_an(header, record, InfoPolicy::Strict)
                .rs_with_context(|| format!("record {number}"))?;
        }
        Ok(change)
    }

    fn select_query(
        &self,
        header: &vcf::Header,
        record: &RecordBuf,
        sample_count: usize,
        number: u64,
    ) -> Result<Vec<bool>> {
        let expression = self
            .expression
            .as_ref()
            .ok_or_else(|| invalid(number, "query expression is not initialized"))?;
        let truth = expression
            .evaluate(header, record)
            .map_err(|error| invalid(number, format!("evaluating setGT expression: {error}")))?;
        match (truth.sample_passes(), truth.sample_selection()) {
            (Some(passes), Some(selection)) if passes.len() == sample_count => Ok(passes
                .iter()
                .zip(selection)
                .map(|(passes, selected)| {
                    if self.exclude {
                        *selected && !*passes
                    } else {
                        *passes
                    }
                })
                .collect()),
            (Some(_), Some(_)) => Err(invalid(
                number,
                "expression sample count differs from record",
            )),
            (None, None) => Ok(vec![
                if self.exclude {
                    !truth.site_passes()
                } else {
                    truth.site_passes()
                };
                sample_count
            ]),
            _ => Err(invalid(number, "expression sample selection is incomplete")),
        }
    }
}

fn read_genotypes(record: &RecordBuf, number: u64) -> Result<Vec<Genotype>> {
    let sample_count = record.samples().values().count();
    if sample_count == 0 {
        return Ok(Vec::new());
    }
    let genotypes = record
        .samples()
        .select("GT")
        .ok_or_else(|| invalid(number, "record has no FORMAT/GT field"))?;
    (0..sample_count)
        .map(|sample| match genotypes.get(sample) {
            Some(Some(SampleValue::Genotype(genotype))) => {
                validate_alleles(record, genotype, number, sample)?;
                Ok(genotype.clone())
            }
            Some(Some(_)) => Err(invalid(
                number,
                format!(
                    "sample {} FORMAT/GT is not encoded as a genotype",
                    sample + 1
                ),
            )),
            _ => Err(invalid(
                number,
                format!("sample {} has no FORMAT/GT value", sample + 1),
            )),
        })
        .collect()
}

fn validate_alleles(
    record: &RecordBuf,
    genotype: &Genotype,
    number: u64,
    sample: usize,
) -> Result<()> {
    let alternate_count = record.alternate_bases().as_ref().len();
    if let Some(position) = genotype
        .as_ref()
        .iter()
        .filter_map(|allele| allele.position())
        .find(|position| *position > alternate_count)
    {
        return Err(invalid(
            number,
            format!(
                "sample {} genotype allele index {position} exceeds {alternate_count} ALT alleles",
                sample + 1
            ),
        ));
    }
    Ok(())
}

fn select_binomial(
    record: &RecordBuf,
    genotypes: &[Genotype],
    binomial: &target::Binomial,
    number: u64,
) -> Result<Vec<bool>> {
    let values = record.samples().select(&binomial.tag).ok_or_else(|| {
        invalid(
            number,
            format!("record has no FORMAT/{} field", binomial.tag),
        )
    })?;
    let expected = record.alternate_bases().as_ref().len() + 1;
    genotypes
        .iter()
        .enumerate()
        .map(|(sample, genotype)| {
            let positions = genotype
                .as_ref()
                .iter()
                .filter_map(|allele| allele.position())
                .collect::<Vec<_>>();
            if genotype.as_ref().len() != 2 || positions.len() != 2 || positions[0] == positions[1]
            {
                return Ok(false);
            }
            let Some(Some(SampleValue::Array(SampleArray::Integer(counts)))) = values.get(sample)
            else {
                return Err(invalid(
                    number,
                    format!(
                        "sample {} FORMAT/{} is missing or not an integer array",
                        sample + 1,
                        binomial.tag
                    ),
                ));
            };
            if counts.len() != expected {
                return Err(invalid(
                    number,
                    format!(
                        "sample {} FORMAT/{} has {} values, expected {expected}",
                        sample + 1,
                        binomial.tag,
                        counts.len()
                    ),
                ));
            }
            let count = |position: usize| {
                counts[position].filter(|value| *value >= 0).ok_or_else(|| {
                    invalid(
                        number,
                        format!(
                            "sample {} FORMAT/{} allele {position} is missing or negative",
                            sample + 1,
                            binomial.tag
                        ),
                    )
                })
            };
            let probability = binomial_two_sided(count(positions[0])?, count(positions[1])?);
            Ok(probability
                .is_some_and(|value| compare(value, binomial.threshold, binomial.comparison)))
        })
        .collect()
}

fn compare(left: f64, right: f64, comparison: Comparison) -> bool {
    match comparison {
        Comparison::Less => left < right,
        Comparison::LessEqual => left <= right,
        Comparison::Equal => left == right,
        Comparison::GreaterEqual => left >= right,
        Comparison::Greater => left > right,
    }
}

fn compile(source: &str, header: &vcf::Header) -> Result<Compiled> {
    Compiled::bind(source, header)
        .map_err(|error| config(format!("invalid setGT expression: {error}")))
}

fn config(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

fn invalid(number: u64, message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(format!("record {number}: {}", message.into()))
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        self as vcf,
        variant::{
            RecordBuf,
            record_buf::{
                Samples,
                info::field::{Value as InfoValue, value::Array as InfoArray},
                samples::sample::Value as SampleValue,
            },
        },
    };

    use super::*;

    fn fixture() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"allele number\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depths\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(
            b"chr1\t7\t.\tA\tC,G\t50\tPASS\tAC=3,2;AN=10\tGT:AD\t0/1:4,4,0\t./.:.,.,.\t1/.:0,3,.\t1/2:0,3,5\t0|0:8,0,0\t2:0,0,7".as_slice(),
        )
        .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn sample_fixture(samples: &[&str], info: &str) -> (vcf::Header, RecordBuf) {
        let names = (1..=samples.len())
            .map(|index| format!("S{index}"))
            .collect::<Vec<_>>()
            .join("\t");
        let header: vcf::Header = format!(
            "##fileformat=VCFv4.3\n\
             ##contig=<ID=chr1>\n\
             ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count\">\n\
             ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"allele number\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
             ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depths\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{names}\n"
        )
        .parse()
        .unwrap();
        let line = format!(
            "chr1\t7\t.\tA\tC,G\t50\tPASS\t{info}\tGT:AD\t{}",
            samples.join("\t")
        );
        let raw = vcf::Record::try_from(line.as_bytes()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn program(
        header: &vcf::Header,
        targets: &[&str],
        replacement: &str,
        query: Option<Query>,
        seed: i64,
    ) -> Program {
        Program::bind(
            header,
            target::Target::parse(
                &targets
                    .iter()
                    .map(|value| (*value).to_owned())
                    .collect::<Vec<_>>(),
            )
            .unwrap(),
            replacement::Replacement::parse(replacement).unwrap(),
            query,
            seed,
        )
        .unwrap()
    }

    #[test]
    fn selects_every_missing_state_and_all_genotypes() {
        let (header, record) = fixture();
        for (target, expected) in [
            (".", &[false, true, true, false, false, false][..]),
            ("./x", &[false, false, true, false, false, false][..]),
            ("./.", &[false, true, false, false, false, false][..]),
            ("a", &[true, true, true, true, true, true][..]),
        ] {
            let mut program = program(&header, &[target], "p", None, 0);
            assert_eq!(program.select(&header, &record, 7).unwrap(), expected);
        }
    }

    #[test]
    fn selects_site_sample_include_exclude_and_selected_sample_queries() {
        let (header, record) = fixture();
        for (query, expected) in [
            (
                Query::Include("QUAL > 10".to_owned()),
                &[true, true, true, true, true, true][..],
            ),
            (
                Query::Exclude("QUAL > 10".to_owned()),
                &[false, false, false, false, false, false][..],
            ),
            (
                Query::Include("GT = 'het'".to_owned()),
                &[true, false, false, true, false, false][..],
            ),
            (
                Query::Exclude("GT = 'het'".to_owned()),
                &[false, true, true, false, true, true][..],
            ),
            (
                Query::Exclude("GT[0,3] = 'het'".to_owned()),
                &[false, false, false, false, false, false][..],
            ),
        ] {
            let mut program = program(&header, &["q"], "p", Some(query), 0);
            assert_eq!(program.select(&header, &record, 7).unwrap(), expected);
        }
    }

    #[test]
    fn selects_binomial_probabilities_and_composes_random_after_principal() {
        let (header, record) = fixture();
        let mut binomial = program(&header, &["b:AD<0.8"], "p", None, 0);
        assert_eq!(
            binomial.select(&header, &record, 7).unwrap(),
            [false, false, false, true, false, false]
        );

        let mut random = program(&header, &["a", "r:0.3"], "p", None, 7);
        assert_eq!(
            random.select(&header, &record, 7).unwrap(),
            [true, false, true, true, false, true]
        );

        let mut composed = program(&header, &[".", "r:0.3"], "p", None, 7);
        assert_eq!(
            composed.select(&header, &record, 7).unwrap(),
            [false, true, false, false, false, false]
        );

        let boundary = format!("r:{}", f64::from_bits(0x3fd10d6bf5d44040));
        let mut inclusive = program(&header, &["a", &boundary], "p", None, 7);
        assert!(inclusive.select(&header, &record, 7).unwrap()[0]);
    }

    #[test]
    fn query_contract_is_bound_before_records_are_read() {
        let (header, _) = fixture();
        let target = |value: &str| target::Target::parse(&[value.to_owned()]).unwrap();
        let replacement = || replacement::Replacement::parse("p").unwrap();

        assert!(Program::bind(&header, target("q"), replacement(), None, 0).is_err());
        assert!(
            Program::bind(
                &header,
                target("a"),
                replacement(),
                Some(Query::Include("QUAL > 0".to_owned())),
                0,
            )
            .is_err()
        );
        assert!(
            Program::bind(
                &header,
                target("q"),
                replacement(),
                Some(Query::Include("UNKNOWN > 0".to_owned())),
                0,
            )
            .is_err()
        );
    }

    fn spelling(record: &RecordBuf, sample: usize) -> String {
        let Some(Some(SampleValue::Genotype(genotype))) =
            record.samples().get_index(sample).unwrap().get("GT")
        else {
            panic!("missing genotype")
        };
        let mut output = String::new();
        for (index, allele) in genotype.as_ref().iter().enumerate() {
            if index > 0 {
                output.push(match allele.phasing() {
                    noodles_vcf::variant::record::samples::series::value::genotype::Phasing::Phased => '|',
                    noodles_vcf::variant::record::samples::series::value::genotype::Phasing::Unphased => '/',
                });
            }
            match allele.position() {
                Some(position) => output.push_str(&position.to_string()),
                None => output.push('.'),
            }
        }
        output
    }

    #[test]
    fn applies_missing_and_reconciles_existing_counts() {
        let (header, mut record) = fixture();
        let mut program = program(&header, &["a"], ".", None, 0);
        let change = program.apply(&header, &mut record, 7).unwrap();

        assert_eq!(change.genotypes, 5);
        assert_eq!(change.alleles, 8);
        assert_eq!(
            (0..6)
                .map(|sample| spelling(&record, sample))
                .collect::<Vec<_>>(),
            ["./.", "./.", "./.", "./.", "./.", "."]
        );
        assert!(matches!(
            record.info().get("AC"),
            Some(Some(InfoValue::Array(InfoArray::Integer(values))))
                if values == &[Some(0), Some(0)]
        ));
        assert!(matches!(
            record.info().get("AN"),
            Some(Some(InfoValue::Integer(0)))
        ));
    }

    #[test]
    fn phases_unphases_sorts_and_inverts_without_fixed_ploidy() {
        let samples = ["1/0:3,4,0", "2|0:2,0,5", "0/1/2:2,3,4"];
        let (header, record) = sample_fixture(&samples, "AC=2,2;AN=7");

        let mut phased = record.clone();
        program(&header, &["a"], "p", None, 0)
            .apply(&header, &mut phased, 7)
            .unwrap();
        assert_eq!(spelling(&phased, 0), "1|0");
        assert_eq!(spelling(&phased, 2), "0|1|2");

        let mut unphased = record.clone();
        program(&header, &["a"], "u", None, 0)
            .apply(&header, &mut unphased, 7)
            .unwrap();
        assert_eq!(
            (0..3)
                .map(|sample| spelling(&unphased, sample))
                .collect::<Vec<_>>(),
            ["0/1", "0/2", "0/1/2"]
        );

        let mut inverted = record;
        let change = program(&header, &["a"], "i", None, 0)
            .apply(&header, &mut inverted, 7)
            .unwrap();
        assert_eq!(spelling(&inverted, 0), "0/1");
        assert_eq!(spelling(&inverted, 1), "0|2");
        assert_eq!(spelling(&inverted, 2), "0/1/2");
        assert_eq!(change.genotypes, 2);
    }

    #[test]
    fn resolves_major_minor_depth_and_custom_templates() {
        let (header, record) = fixture();

        let mut major = record.clone();
        program(&header, &["a"], "M", None, 0)
            .apply(&header, &mut major, 7)
            .unwrap();
        assert_eq!(spelling(&major, 0), "0/0");

        let mut minor = record.clone();
        program(&header, &["a"], "m", None, 0)
            .apply(&header, &mut minor, 7)
            .unwrap();
        assert_eq!(spelling(&minor, 0), "1/1");

        let mut depth = record.clone();
        program(&header, &["a"], "X", None, 0)
            .apply(&header, &mut depth, 7)
            .unwrap();
        assert_eq!(
            (0..6)
                .map(|sample| spelling(&depth, sample))
                .collect::<Vec<_>>(),
            ["0/0", "./.", "1/1", "2/2", "0/0", "2"]
        );

        let mut custom = record;
        program(
            &header,
            &["q"],
            "c:0/m|M/X",
            Some(Query::Include("GT[0] = 'het'".to_owned())),
            0,
        )
        .apply(&header, &mut custom, 7)
        .unwrap();
        assert_eq!(spelling(&custom, 0), "0/1|0/0");
        assert_eq!(spelling(&custom, 1), "./.");
    }

    #[test]
    fn custom_out_of_range_positions_become_missing_and_replace_ploidy() {
        let (header, mut record) = fixture();
        let change = program(
            &header,
            &["q"],
            "c:7|0",
            Some(Query::Include("GT[5] = 'A'".to_owned())),
            0,
        )
        .apply(&header, &mut record, 7)
        .unwrap();

        assert_eq!(spelling(&record, 5), ".|0");
        assert_eq!(change.genotypes, 1);
        assert_eq!(change.alleles, 2);
    }

    #[test]
    fn selected_malformed_depth_and_unresolvable_counts_fail_with_context() {
        for sample in ["0/1:4,2", "0/1:4,-1,0"] {
            let (header, mut record) = sample_fixture(&[sample], "AC=1,0;AN=2");
            let error = program(&header, &["a"], "X", None, 0)
                .apply(&header, &mut record, 9)
                .unwrap_err();
            let message = error.to_string();
            assert!(message.contains("record 9"), "{message}");
            assert!(message.contains("sample 1"), "{message}");
        }

        let (header, mut record) = sample_fixture(&["./.:.,.,."], "AC=0,0;AN=0");
        let error = program(&header, &["a"], "M", None, 0)
            .apply(&header, &mut record, 9)
            .unwrap_err();
        assert!(error.to_string().contains("record 9"));
    }

    #[test]
    fn genotype_and_binomial_structure_fail_loud_with_record_context() {
        let (header, mut wrong_type) = fixture();
        let keys = wrong_type.samples().keys().clone();
        let genotype = keys.as_ref().get_index_of("GT").unwrap();
        let mut rows = wrong_type
            .samples()
            .values()
            .map(|sample| sample.values().to_vec())
            .collect::<Vec<_>>();
        rows[0][genotype] = Some(SampleValue::String("0/1".to_owned()));
        *wrong_type.samples_mut() = Samples::new(keys, rows);
        let error = program(&header, &["a"], "p", None, 0)
            .apply(&header, &mut wrong_type, 11)
            .unwrap_err();
        assert!(error.to_string().contains("record 11"));
        assert!(error.to_string().contains("sample 1"));

        let (header, mut out_of_range) = sample_fixture(&["3/1:0,2,3"], "AC=1,0;AN=2");
        let error = program(&header, &["a"], "p", None, 0)
            .apply(&header, &mut out_of_range, 12)
            .unwrap_err();
        assert!(error.to_string().contains("allele index 3"));

        for sample in ["0/1:4,2", "0/1:4,-1,0"] {
            let (header, mut record) = sample_fixture(&[sample], "AC=1,0;AN=2");
            let error = program(&header, &["b:AD<1"], "p", None, 0)
                .apply(&header, &mut record, 13)
                .unwrap_err();
            assert!(error.to_string().contains("record 13"));
            assert!(error.to_string().contains("sample 1"));
        }
    }

    #[test]
    fn sites_only_records_are_unchanged_without_a_gt_schema() {
        let sample_header: vcf::Header = "##fileformat=VCFv4.3\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        assert!(
            Program::bind(
                &sample_header,
                target::Target::parse(&["a".to_owned()]).unwrap(),
                replacement::Replacement::parse("p").unwrap(),
                None,
                0,
            )
            .is_err()
        );

        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(b"chr1\t1\t.\tA\tC\t10\tPASS\t.".as_slice()).unwrap();
        let mut record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let mut program = Program::bind(
            &header,
            target::Target::parse(&["a".to_owned()]).unwrap(),
            replacement::Replacement::parse("p").unwrap(),
            None,
            0,
        )
        .unwrap();

        assert_eq!(
            program.apply(&header, &mut record, 1).unwrap(),
            Change::default()
        );
    }
}
