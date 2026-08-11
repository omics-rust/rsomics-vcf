mod gaps;

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

use crate::{expression::Compiled, regions::RegionSet};

#[derive(Clone, Debug)]
pub(crate) struct Mask {
    regions: RegionSet,
    negate: bool,
}

impl Mask {
    pub(crate) fn new(regions: RegionSet, negate: bool) -> Self {
        Self { regions, negate }
    }

    fn passes(&self, record: &RecordBuf) -> bool {
        self.regions.matches(record) == self.negate
    }
}

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

#[derive(Clone, Debug, Default)]
pub(crate) struct Options {
    pub logic: Logic,
    pub mask: Option<Mask>,
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

#[derive(Clone, Debug)]
pub(crate) struct Program {
    expression: Option<Compiled>,
    mask: Option<Mask>,
    logic: Logic,
    failure_filter: Option<String>,
    mode: AnnotationMode,
    set_genotypes: Option<GenotypeReplacement>,
}

pub(crate) struct Processor {
    program: Option<Program>,
    gaps: Option<gaps::GapBuffer>,
}

impl Processor {
    pub(crate) fn bind(
        header: &mut vcf::Header,
        source: Option<&str>,
        options: Options,
        mut gap_options: gaps::Options,
    ) -> Result<Self> {
        let has_program = source.is_some() || options.mask.is_some();
        let has_gaps = gap_options.snp_gap.is_some() || gap_options.indel_gap.is_some();
        if !has_program && !has_gaps {
            return Err(RsomicsError::ConfigError(
                "filter requires an expression, mask, SnpGap, or IndelGap".to_owned(),
            ));
        }
        if !has_program && options.set_genotypes.is_some() {
            return Err(RsomicsError::ConfigError(
                "setting failed genotypes requires an expression or mask".to_owned(),
            ));
        }

        gap_options.soft = options.soft_filter.is_some();
        let program = has_program
            .then(|| Program::bind(header, source, options))
            .transpose()?;
        let gaps = has_gaps
            .then(|| gaps::GapBuffer::new(header, gap_options))
            .transpose()?;
        Ok(Self { program, gaps })
    }

    pub(crate) fn push<F>(
        &mut self,
        header: &vcf::Header,
        mut record: RecordBuf,
        emit: &mut F,
    ) -> Result<()>
    where
        F: FnMut(RecordBuf) -> Result<()>,
    {
        if let Some(program) = &self.program
            && program
                .apply(header, &mut record, RetainFailed::No)?
                .disposition()
                == Disposition::Drop
        {
            return Ok(());
        }
        let Some(gaps) = &mut self.gaps else {
            return emit(record);
        };
        for output in gaps.push(record)? {
            if output.disposition() == Disposition::Keep {
                emit(output.into_record())?;
            }
        }
        Ok(())
    }

    pub(crate) fn finish<F>(&mut self, emit: &mut F) -> Result<()>
    where
        F: FnMut(RecordBuf) -> Result<()>,
    {
        let Some(gaps) = &mut self.gaps else {
            return Ok(());
        };
        for output in gaps.finish()? {
            if output.disposition() == Disposition::Keep {
                emit(output.into_record())?;
            }
        }
        Ok(())
    }
}

impl Program {
    pub(crate) fn bind(
        header: &mut vcf::Header,
        source: Option<&str>,
        options: Options,
    ) -> Result<Self> {
        if source.is_none() && options.mask.is_none() {
            return Err(RsomicsError::ConfigError(
                "filter requires an expression or mask".to_owned(),
            ));
        }
        if options.mask.is_some() && options.soft_filter.is_none() {
            return Err(RsomicsError::ConfigError(
                "mask requires a soft filter".to_owned(),
            ));
        }
        let expression = source
            .map(|source| {
                Compiled::bind(source, header).map_err(|error| {
                    RsomicsError::ConfigError(format!("invalid filter expression: {error}"))
                })
            })
            .transpose()?;
        let description = source.map_or_else(
            || "Record masked by region".to_owned(),
            |source| match options.logic {
                Logic::Include => format!("Set if not true: {source}"),
                Logic::Exclude => format!("Set if true: {source}"),
            },
        );
        let failure_filter = options
            .soft_filter
            .map(|name| register_filter(header, name, description));
        Ok(Self {
            expression,
            mask: options.mask,
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
        let (expression_passes, sample_passes) = if let Some(expression) = &self.expression {
            let truth = expression.evaluate(header, record).map_err(|error| {
                RsomicsError::InvalidInput(format!("evaluating filter expression: {error}"))
            })?;
            (
                self.logic.accepts(truth.site_passes()),
                truth.sample_passes().map(|passes| {
                    passes
                        .iter()
                        .map(|passes| self.logic.accepts(*passes))
                        .collect::<Vec<_>>()
                }),
            )
        } else {
            (true, None)
        };
        let mask_passes = self.mask.as_ref().is_none_or(|mask| mask.passes(record));
        let site_passes = expression_passes && mask_passes;

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

fn register_filter(header: &mut vcf::Header, requested: String, description: String) -> String {
    let name = if requested == "+" {
        (1..)
            .map(|index| format!("Filter{index}"))
            .find(|name| !header.filters().contains_key(name))
            .unwrap()
    } else {
        requested
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
    use crate::regions::OverlapMode;

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

    fn variant_fixture() -> (vcf::Header, RecordBuf) {
        let (header, _) = fixture();
        let record =
            vcf::Record::try_from(b"chr1\t10\t.\tAACGT\tAATGT\t10\tPASS\t.\tDP\t8\t20".as_slice())
                .unwrap();
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

    fn process(
        mut processor: Processor,
        header: &vcf::Header,
        records: &[&[u8]],
    ) -> Vec<RecordBuf> {
        let mut output = Vec::new();
        for record in records {
            let record = vcf::Record::try_from(*record).unwrap();
            let record = RecordBuf::try_from_variant_record(header, &record).unwrap();
            processor
                .push(header, record, &mut |record| {
                    output.push(record);
                    Ok(())
                })
                .unwrap();
        }
        processor
            .finish(&mut |record| {
                output.push(record);
                Ok(())
            })
            .unwrap();
        output
    }

    #[test]
    fn hard_include_and_exclude_choose_record_disposition() {
        let (mut header, record) = fixture();
        let include = Program::bind(
            &mut header,
            Some("QUAL >= 20"),
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
            Some("QUAL < 20"),
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
            Some("FMT/DP >= 10"),
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
            Some("FMT/DP >= 10"),
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
            Some("QUAL >= 20"),
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
            Some("QUAL >= 20"),
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
            Some("QUAL >= 5"),
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
            Some("QUAL >= 5"),
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
            Some("QUAL >= 20"),
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
            Some("FMT/DP >= 10"),
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
            Some("FMT/DP >= 10"),
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
            Some("FMT/DP >= 10"),
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
            Some("QUAL >= 20"),
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

    #[test]
    fn mask_only_failure_rewrites_every_sample() {
        let (mut header, mut record) =
            genotype_fixture(b"chr1\t1\t.\tA\tC\t10\tPASS\tAC=3;AN=4\tGT:DP\t0/1:8\t1/1:20");
        let regions = RegionSet::parse(["chr1:1-1".to_owned()], OverlapMode::Position).unwrap();
        let program = Program::bind(
            &mut header,
            None,
            Options {
                mask: Some(Mask::new(regions, false)),
                soft_filter: Some("Masked".to_owned()),
                set_genotypes: Some(GenotypeReplacement::Missing),
                ..Options::default()
            },
        )
        .unwrap();
        let decision = program
            .apply(&header, &mut record, RetainFailed::No)
            .unwrap();
        assert!(!decision.site_passes());
        assert!(decision.sample_passes().is_none());
        assert_eq!(genotype(&record, 0), [None, None]);
        assert_eq!(genotype(&record, 1), [None, None]);
        assert_eq!(integer_array_info(&record, "AC"), [Some(0)]);
        assert_eq!(integer_info(&record, "AN"), Some(0));
    }

    #[test]
    fn masks_share_position_record_and_variant_overlap_semantics() {
        let cases = [
            (OverlapMode::Position, "chr1:10-10", false),
            (OverlapMode::Position, "chr1:14-14", true),
            (OverlapMode::Record, "chr1:14-14", false),
            (OverlapMode::Record, "chr1:15-15", true),
            (OverlapMode::Variant, "chr1:10-10", true),
            (OverlapMode::Variant, "chr1:14-14", false),
            (OverlapMode::Variant, "chr1:12-12", false),
        ];
        for (overlap, region, expected) in cases {
            let (mut header, mut record) = variant_fixture();
            let regions = RegionSet::parse([region.to_owned()], overlap).unwrap();
            let program = Program::bind(
                &mut header,
                None,
                Options {
                    mask: Some(Mask::new(regions, false)),
                    soft_filter: Some("Masked".to_owned()),
                    ..Options::default()
                },
            )
            .unwrap();
            let decision = program
                .apply(&header, &mut record, RetainFailed::No)
                .unwrap();
            assert_eq!(decision.site_passes(), expected, "{overlap:?} {region}");
        }
    }

    #[test]
    fn negated_masks_and_expressions_both_must_pass() {
        let (mut header, record) = fixture();
        let regions = RegionSet::parse(["chr1:1-1".to_owned()], OverlapMode::Position).unwrap();
        let expression_passes = Program::bind(
            &mut header,
            Some("FMT/DP >= 10"),
            Options {
                mask: Some(Mask::new(regions.clone(), false)),
                soft_filter: Some("Masked".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        let mut masked = record.clone();
        let decision = expression_passes
            .apply(&header, &mut masked, RetainFailed::No)
            .unwrap();
        assert!(!decision.site_passes());
        assert_eq!(decision.sample_passes(), Some(&[false, true][..]));

        let expression_fails = Program::bind(
            &mut header,
            Some("QUAL >= 20"),
            Options {
                mask: Some(Mask::new(regions.clone(), true)),
                soft_filter: Some("ExprFail".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        let mut failed = record.clone();
        assert!(
            !expression_fails
                .apply(&header, &mut failed, RetainFailed::No)
                .unwrap()
                .site_passes()
        );

        let both_pass = Program::bind(
            &mut header,
            Some("QUAL >= 5"),
            Options {
                mask: Some(Mask::new(regions, true)),
                soft_filter: Some("Unused".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        let mut passed = record;
        assert!(
            both_pass
                .apply(&header, &mut passed, RetainFailed::No)
                .unwrap()
                .site_passes()
        );
    }

    #[test]
    fn masks_require_soft_filters_and_mask_only_headers_are_explicit() {
        let (mut header, _) = fixture();
        let error = Program::bind(&mut header, None, Options::default()).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("filter requires an expression or mask")
        );

        let regions = RegionSet::parse(["chr1:1-1".to_owned()], OverlapMode::Position).unwrap();
        let error = Program::bind(
            &mut header,
            None,
            Options {
                mask: Some(Mask::new(regions.clone(), false)),
                ..Options::default()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("mask requires a soft filter"));

        Program::bind(
            &mut header,
            None,
            Options {
                mask: Some(Mask::new(regions, false)),
                soft_filter: Some("Masked".to_owned()),
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(
            header.filters().get("Masked").unwrap().description(),
            "Record masked by region"
        );
    }

    #[test]
    fn processor_filters_records_before_gap_context() {
        let (mut header, _) = fixture();
        let processor = Processor::bind(
            &mut header,
            Some("QUAL >= 10"),
            Options::default(),
            gaps::Options {
                snp_gap: Some(gaps::SnpGap {
                    distance: 3,
                    types: crate::variant_type::INDEL,
                }),
                ..gaps::Options::default()
            },
        )
        .unwrap();
        let output = process(
            processor,
            &header,
            &[
                b"chr1\t10\t.\tA\tC\t10\tPASS\t.\tDP\t8\t20",
                b"chr1\t12\t.\tA\tAT\t1\tPASS\t.\tDP\t8\t20",
            ],
        );

        assert_eq!(
            output
                .iter()
                .map(|record| usize::from(record.variant_start().unwrap()))
                .collect::<Vec<_>>(),
            [10]
        );
        assert!(output[0].filters().is_pass());
    }

    #[test]
    fn processor_keeps_soft_failures_in_gap_context() {
        let (mut header, _) = fixture();
        let processor = Processor::bind(
            &mut header,
            Some("QUAL >= 10"),
            Options {
                soft_filter: Some("ExprFail".to_owned()),
                ..Options::default()
            },
            gaps::Options {
                snp_gap: Some(gaps::SnpGap {
                    distance: 3,
                    types: crate::variant_type::INDEL,
                }),
                ..gaps::Options::default()
            },
        )
        .unwrap();
        let output = process(
            processor,
            &header,
            &[
                b"chr1\t10\t.\tA\tC\t10\tPASS\t.\tDP\t8\t20",
                b"chr1\t12\t.\tA\tAT\t1\tPASS\t.\tDP\t8\t20",
            ],
        );

        assert_eq!(output.len(), 2);
        assert_eq!(
            output[0].filters().as_ref().iter().collect::<Vec<_>>(),
            ["SnpGap"]
        );
        assert_eq!(
            output[1].filters().as_ref().iter().collect::<Vec<_>>(),
            ["ExprFail"]
        );
    }

    #[test]
    fn processor_accepts_gap_only_filters() {
        let (mut header, _) = fixture();
        Processor::bind(
            &mut header,
            None,
            Options::default(),
            gaps::Options {
                indel_gap: Some(2),
                ..gaps::Options::default()
            },
        )
        .unwrap();

        assert!(header.filters().contains_key("IndelGap"));
    }
}
