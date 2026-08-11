use noodles_vcf::{
    self as vcf,
    header::record::value::{Map, map::Filter as HeaderFilter},
    variant::{
        RecordBuf,
        record::samples::series::value::genotype::Phasing,
        record_buf::{
            Filters, Samples,
            info::field::{Value as InfoValue, value::Array},
            samples::sample::Value as SampleValue,
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use crate::expression::Compiled;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum Logic {
    #[default]
    Include,
    Exclude,
}

impl Logic {
    fn accepts(self, value: bool) -> bool {
        match self {
            Self::Include => value,
            Self::Exclude => !value,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum AnnotationMode {
    #[default]
    Replace,
    Add,
    ResetPassed,
    AddAndResetPassed,
}

impl AnnotationMode {
    fn adds(self) -> bool {
        matches!(self, Self::Add | Self::AddAndResetPassed)
    }

    fn resets_passed(self) -> bool {
        matches!(self, Self::ResetPassed | Self::AddAndResetPassed)
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct Options {
    pub logic: Logic,
    pub soft_filter: Option<String>,
    pub mode: AnnotationMode,
    pub set_genotypes: Option<GenotypeReplacement>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum GenotypeReplacement {
    Missing,
    Reference,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RetainFailed {
    No,
    Yes,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Disposition {
    Keep,
    Drop,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Decision {
    disposition: Disposition,
    site_passes: bool,
    sample_passes: Option<Vec<bool>>,
}

impl Decision {
    pub(crate) fn disposition(&self) -> Disposition {
        self.disposition
    }

    pub(crate) fn site_passes(&self) -> bool {
        self.site_passes
    }

    pub(crate) fn sample_passes(&self) -> Option<&[bool]> {
        self.sample_passes.as_deref()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Program {
    expression: Compiled,
    logic: Logic,
    failure_filter: Option<String>,
    mode: AnnotationMode,
    set_genotypes: Option<GenotypeReplacement>,
}

impl Program {
    pub(crate) fn bind(header: &mut vcf::Header, source: &str, options: Options) -> Result<Self> {
        let expression = Compiled::bind(source, header).map_err(|error| {
            RsomicsError::ConfigError(format!("invalid filter expression: {error}"))
        })?;
        let failure_filter = options
            .soft_filter
            .map(|name| register_filter(header, name, source, options.logic));
        Ok(Self {
            expression,
            logic: options.logic,
            failure_filter,
            mode: options.mode,
            set_genotypes: options.set_genotypes,
        })
    }

    pub(crate) fn failure_filter(&self) -> Option<&str> {
        self.failure_filter.as_deref()
    }

    pub(crate) fn apply(
        &self,
        header: &vcf::Header,
        record: &mut RecordBuf,
        retain_failed: RetainFailed,
    ) -> Result<Decision> {
        let truth = self.expression.evaluate(header, record).map_err(|error| {
            RsomicsError::InvalidInput(format!("evaluating filter expression: {error}"))
        })?;
        let site_passes = self.logic.accepts(truth.site_passes());
        let sample_passes = truth.sample_passes().map(|passes| {
            passes
                .iter()
                .map(|passes| self.logic.accepts(*passes))
                .collect()
        });

        if site_passes {
            if self.mode.resets_passed() || record.filters().as_ref().is_empty() {
                *record.filters_mut() = Filters::pass();
            }
        } else if let Some(filter) = &self.failure_filter {
            if self.mode.adds() {
                let filters = record.filters_mut().as_mut();
                filters.shift_remove("PASS");
                filters.insert(filter.clone());
            } else {
                *record.filters_mut() = [filter.clone()].into_iter().collect();
            }
        }
        if let Some(replacement) = self.set_genotypes {
            replace_failed_genotypes(record, site_passes, sample_passes.as_deref(), replacement)?;
        }

        let disposition = if site_passes
            || self.failure_filter.is_some()
            || self.set_genotypes.is_some()
            || retain_failed == RetainFailed::Yes
        {
            Disposition::Keep
        } else {
            Disposition::Drop
        };
        Ok(Decision {
            disposition,
            site_passes,
            sample_passes,
        })
    }
}

fn replace_failed_genotypes(
    record: &mut RecordBuf,
    site_passes: bool,
    sample_passes: Option<&[bool]>,
    replacement: GenotypeReplacement,
) -> Result<()> {
    let keys = record.samples().keys().clone();
    let Some(genotype_index) = keys.as_ref().get_index_of("GT") else {
        return Ok(());
    };
    let mut samples: Vec<_> = record
        .samples()
        .values()
        .map(|sample| sample.values().to_vec())
        .collect();
    if sample_passes.is_some_and(|passes| passes.len() != samples.len()) {
        return Err(RsomicsError::InvalidInput(
            "filter sample count does not match the record".to_owned(),
        ));
    }

    let alternate_count = record.alternate_bases().as_ref().len();
    let mut ac = valid_ac(record, alternate_count);
    let mut an = valid_an(record);
    for (sample_index, sample) in samples.iter_mut().enumerate() {
        let passes = sample_passes.map_or(site_passes, |passes| passes[sample_index]);
        if passes {
            continue;
        }
        let Some(Some(value)) = sample.get_mut(genotype_index) else {
            continue;
        };
        let SampleValue::Genotype(genotype) = value else {
            return Err(RsomicsError::InvalidInput(
                "FORMAT/GT is not encoded as a genotype".to_owned(),
            ));
        };
        for allele in genotype.as_mut() {
            let position = allele.position();
            if let Some(position) = position {
                if position > alternate_count {
                    return Err(RsomicsError::InvalidInput(format!(
                        "genotype allele index {position} exceeds {alternate_count} ALT alleles"
                    )));
                }
                if position > 0 {
                    decrement_ac(&mut ac, position - 1)?;
                }
            }
            match replacement {
                GenotypeReplacement::Missing => {
                    if position.is_some() {
                        adjust_an(&mut an, -1)?;
                    }
                    *allele.position_mut() = None;
                }
                GenotypeReplacement::Reference => {
                    if position.is_none() {
                        adjust_an(&mut an, 1)?;
                    }
                    *allele.position_mut() = Some(0);
                }
            }
            *allele.phasing_mut() = Phasing::Unphased;
        }
    }
    *record.samples_mut() = Samples::new(keys, samples);
    if let Some(ac) = ac {
        record.info_mut().insert(
            "AC".to_owned(),
            Some(InfoValue::Array(Array::Integer(
                ac.into_iter().map(Some).collect(),
            ))),
        );
    }
    if let Some(an) = an {
        record
            .info_mut()
            .insert("AN".to_owned(), Some(InfoValue::Integer(an)));
    }
    Ok(())
}

fn valid_ac(record: &RecordBuf, alternate_count: usize) -> Option<Vec<i32>> {
    let Some(Some(InfoValue::Array(Array::Integer(values)))) = record.info().get("AC") else {
        return None;
    };
    (values.len() == alternate_count)
        .then(|| values.iter().copied().collect::<Option<Vec<_>>>())
        .flatten()
}

fn valid_an(record: &RecordBuf) -> Option<i32> {
    match record.info().get("AN") {
        Some(Some(InfoValue::Integer(value))) => Some(*value),
        _ => None,
    }
}

fn decrement_ac(ac: &mut Option<Vec<i32>>, index: usize) -> Result<()> {
    let Some(ac) = ac else {
        return Ok(());
    };
    ac[index] = ac[index]
        .checked_sub(1)
        .ok_or_else(|| RsomicsError::InvalidInput("INFO/AC count underflow".to_owned()))?;
    Ok(())
}

fn adjust_an(an: &mut Option<i32>, difference: i32) -> Result<()> {
    let Some(an) = an else {
        return Ok(());
    };
    *an = an
        .checked_add(difference)
        .ok_or_else(|| RsomicsError::InvalidInput("INFO/AN count overflow".to_owned()))?;
    Ok(())
}

fn register_filter(
    header: &mut vcf::Header,
    requested: String,
    source: &str,
    logic: Logic,
) -> String {
    let name = if requested == "+" {
        (1..)
            .map(|index| format!("Filter{index}"))
            .find(|name| !header.filters().contains_key(name))
            .unwrap()
    } else {
        requested
    };
    let description = match logic {
        Logic::Include => format!("Set if not true: {source}"),
        Logic::Exclude => format!("Set if true: {source}"),
    };
    header
        .filters_mut()
        .entry(name.clone())
        .or_insert_with(|| Map::<HeaderFilter>::new(description));
    name
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        self as vcf,
        variant::{
            RecordBuf,
            record_buf::{
                info::field::{Value as InfoValue, value::Array},
                samples::sample::Value as SampleValue,
            },
        },
    };

    use super::*;

    fn fixture() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=LowQual,Description=\"low quality\">\n\
##FILTER=<ID=Filter1,Description=\"existing\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            .parse()
            .unwrap();
        let record =
            vcf::Record::try_from(b"chr1\t1\t.\tA\tC\t10\tLowQual\t.\tDP\t8\t20".as_slice())
                .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &record).unwrap();
        (header, record)
    }

    fn genotype_fixture(record: &[u8]) -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"allele number\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            .parse()
            .unwrap();
        let record = vcf::Record::try_from(record).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &record).unwrap();
        (header, record)
    }

    fn genotype(record: &RecordBuf, sample: usize) -> Vec<Option<usize>> {
        let value = record
            .samples()
            .get_index(sample)
            .unwrap()
            .get("GT")
            .unwrap()
            .unwrap();
        let SampleValue::Genotype(value) = value else {
            panic!("GT is not a genotype")
        };
        value
            .as_ref()
            .iter()
            .map(|allele| allele.position())
            .collect()
    }

    fn integer_info(record: &RecordBuf, key: &str) -> Option<i32> {
        match record.info().get(key) {
            Some(Some(InfoValue::Integer(value))) => Some(*value),
            _ => None,
        }
    }

    fn integer_array_info(record: &RecordBuf, key: &str) -> Vec<Option<i32>> {
        match record.info().get(key) {
            Some(Some(InfoValue::Array(Array::Integer(values)))) => values.clone(),
            _ => Vec::new(),
        }
    }

    #[test]
    fn hard_include_and_exclude_choose_record_disposition() {
        let (mut header, record) = fixture();
        let include = Program::bind(
            &mut header,
            "QUAL >= 20",
            Options {
                logic: Logic::Include,
                ..Options::default()
            },
        )
        .unwrap();
        let mut included = record.clone();
        let decision = include
            .apply(&header, &mut included, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.disposition(), Disposition::Drop);
        assert!(!decision.site_passes());

        let exclude = Program::bind(
            &mut header,
            "QUAL < 20",
            Options {
                logic: Logic::Exclude,
                ..Options::default()
            },
        )
        .unwrap();
        let mut excluded = record;
        let decision = exclude
            .apply(&header, &mut excluded, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.disposition(), Disposition::Drop);

        let decision = include
            .apply(&header, &mut excluded, RetainFailed::Yes)
            .unwrap();
        assert_eq!(decision.disposition(), Disposition::Keep);
    }

    #[test]
    fn sample_decisions_follow_include_and_exclude_logic() {
        let (mut header, record) = fixture();
        let include = Program::bind(
            &mut header,
            "FMT/DP >= 10",
            Options {
                logic: Logic::Include,
                ..Options::default()
            },
        )
        .unwrap();
        let mut included = record.clone();
        let decision = include
            .apply(&header, &mut included, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.sample_passes(), Some(&[false, true][..]));

        let exclude = Program::bind(
            &mut header,
            "FMT/DP >= 10",
            Options {
                logic: Logic::Exclude,
                ..Options::default()
            },
        )
        .unwrap();
        let mut excluded = record;
        let decision = exclude
            .apply(&header, &mut excluded, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.sample_passes(), Some(&[true, false][..]));
    }

    #[test]
    fn soft_filter_replaces_adds_and_resets_filters() {
        let (mut header, record) = fixture();
        let replace = Program::bind(
            &mut header,
            "QUAL >= 20",
            Options {
                soft_filter: Some("ExprFail".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        let mut replaced = record.clone();
        let decision = replace
            .apply(&header, &mut replaced, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.disposition(), Disposition::Keep);
        assert_eq!(
            replaced.filters().as_ref().iter().collect::<Vec<_>>(),
            ["ExprFail"]
        );

        let add = Program::bind(
            &mut header,
            "QUAL >= 20",
            Options {
                soft_filter: Some("Added".to_owned()),
                mode: AnnotationMode::Add,
                ..Options::default()
            },
        )
        .unwrap();
        let mut added = record.clone();
        add.apply(&header, &mut added, RetainFailed::No).unwrap();
        assert_eq!(
            added.filters().as_ref().iter().collect::<Vec<_>>(),
            ["LowQual", "Added"]
        );

        let reset = Program::bind(
            &mut header,
            "QUAL >= 5",
            Options {
                soft_filter: Some("Unused".to_owned()),
                mode: AnnotationMode::ResetPassed,
                ..Options::default()
            },
        )
        .unwrap();
        let mut passed = record;
        reset.apply(&header, &mut passed, RetainFailed::No).unwrap();
        assert!(passed.filters().is_pass());

        let add_and_reset = Program::bind(
            &mut header,
            "QUAL >= 5",
            Options {
                soft_filter: Some("Unused2".to_owned()),
                mode: AnnotationMode::AddAndResetPassed,
                ..Options::default()
            },
        )
        .unwrap();
        let mut passed = replaced;
        add_and_reset
            .apply(&header, &mut passed, RetainFailed::No)
            .unwrap();
        assert!(passed.filters().is_pass());
    }

    #[test]
    fn generated_soft_filter_name_skips_existing_header_ids() {
        let (mut header, _) = fixture();
        let program = Program::bind(
            &mut header,
            "QUAL >= 20",
            Options {
                soft_filter: Some("+".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(program.failure_filter(), Some("Filter2"));
        assert!(header.filters().contains_key("Filter2"));
    }

    #[test]
    fn failed_samples_can_be_set_missing_and_update_valid_ac_an() {
        let (mut header, mut record) =
            genotype_fixture(b"chr1\t1\t.\tA\tC,G\t10\tPASS\tAC=4,3;AN=10\tGT:DP\t0|1:8\t2/2:20");
        let program = Program::bind(
            &mut header,
            "FMT/DP >= 10",
            Options {
                set_genotypes: Some(GenotypeReplacement::Missing),
                ..Options::default()
            },
        )
        .unwrap();
        let decision = program
            .apply(&header, &mut record, RetainFailed::No)
            .unwrap();
        assert_eq!(decision.disposition(), Disposition::Keep);
        assert_eq!(genotype(&record, 0), [None, None]);
        assert_eq!(genotype(&record, 1), [Some(2), Some(2)]);
        assert_eq!(integer_array_info(&record, "AC"), [Some(3), Some(3)]);
        assert_eq!(integer_info(&record, "AN"), Some(8));
    }

    #[test]
    fn failed_samples_can_be_set_reference_and_count_missing_alleles() {
        let (mut header, mut record) =
            genotype_fixture(b"chr1\t1\t.\tA\tC,G\t10\tPASS\tAC=4,3;AN=10\tGT:DP\t./2:8\t1/1:20");
        let program = Program::bind(
            &mut header,
            "FMT/DP >= 10",
            Options {
                set_genotypes: Some(GenotypeReplacement::Reference),
                ..Options::default()
            },
        )
        .unwrap();
        program
            .apply(&header, &mut record, RetainFailed::No)
            .unwrap();
        assert_eq!(genotype(&record, 0), [Some(0), Some(0)]);
        assert_eq!(genotype(&record, 1), [Some(1), Some(1)]);
        assert_eq!(integer_array_info(&record, "AC"), [Some(4), Some(2)]);
        assert_eq!(integer_info(&record, "AN"), Some(11));
    }

    #[test]
    fn malformed_ac_an_are_preserved_while_genotypes_change() {
        let (mut header, mut record) =
            genotype_fixture(b"chr1\t1\t.\tA\tC,G\t10\tPASS\tAC=4;AN=.\tGT:DP\t0/1:8\t2/2:20");
        let original_info = record.info().clone();
        let program = Program::bind(
            &mut header,
            "FMT/DP >= 10",
            Options {
                set_genotypes: Some(GenotypeReplacement::Missing),
                ..Options::default()
            },
        )
        .unwrap();
        program
            .apply(&header, &mut record, RetainFailed::No)
            .unwrap();
        assert_eq!(genotype(&record, 0), [None, None]);
        assert_eq!(record.info(), &original_info);
    }

    #[test]
    fn site_failure_rewrites_every_sample_and_single_alt_ac() {
        let (mut header, mut record) =
            genotype_fixture(b"chr1\t1\t.\tA\tC\t10\tPASS\tAC=3;AN=4\tGT:DP\t0/1:8\t1/1:20");
        let program = Program::bind(
            &mut header,
            "QUAL >= 20",
            Options {
                set_genotypes: Some(GenotypeReplacement::Missing),
                ..Options::default()
            },
        )
        .unwrap();
        program
            .apply(&header, &mut record, RetainFailed::No)
            .unwrap();
        assert_eq!(genotype(&record, 0), [None, None]);
        assert_eq!(genotype(&record, 1), [None, None]);
        assert_eq!(integer_array_info(&record, "AC"), [Some(0)]);
        assert_eq!(integer_info(&record, "AN"), Some(0));
    }
}
