use noodles_vcf::{
    self as vcf,
    header::record::value::{Map, map::Filter as HeaderFilter},
    variant::{RecordBuf, record_buf::Filters},
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

        let disposition =
            if site_passes || self.failure_filter.is_some() || retain_failed == RetainFailed::Yes {
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
    use noodles_vcf::{self as vcf, variant::RecordBuf};

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
}
