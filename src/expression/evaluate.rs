use std::{collections::HashSet, fmt};

use super::{
    bind::{BoundExpression, BoundValue},
    syntax::{BinaryOperator, UnaryOperator},
    value::{self, Atom, Values},
};
use noodles_vcf::{self as vcf, variant::RecordBuf};

mod arithmetic;
mod comparison;
mod function;
mod logical;
mod select;

use arithmetic::{arithmetic, negate};
use comparison::compare;
use logical::logical;

pub(crate) use function::binomial_two_sided;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Truth {
    site: bool,
    samples: Option<SampleTruth>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SampleTruth {
    passes: Vec<bool>,
    selected: Box<[bool]>,
}

impl Truth {
    pub(crate) fn site_passes(&self) -> bool {
        self.site
    }

    pub(crate) fn sample_passes(&self) -> Option<&[bool]> {
        self.samples
            .as_ref()
            .map(|samples| samples.passes.as_slice())
    }

    pub(crate) fn sample_selection(&self) -> Option<&[bool]> {
        self.samples
            .as_ref()
            .map(|samples| samples.selected.as_ref())
    }

    fn site(value: bool) -> Self {
        Self {
            site: value,
            samples: None,
        }
    }

    #[cfg(test)]
    fn samples(passes: Vec<bool>) -> Self {
        let selected = vec![true; passes.len()].into_boxed_slice();
        Self::selected_samples(passes, selected)
    }

    fn selected_samples(passes: Vec<bool>, selected: Box<[bool]>) -> Self {
        let site = passes
            .iter()
            .zip(&selected)
            .any(|(passes, selected)| *passes && *selected);
        Self::with_samples(site, passes, selected)
    }

    fn with_samples(site: bool, mut passes: Vec<bool>, selected: Box<[bool]>) -> Self {
        debug_assert_eq!(passes.len(), selected.len());
        for (passes, selected) in passes.iter_mut().zip(&selected) {
            *passes &= *selected;
        }
        Self {
            site,
            samples: Some(SampleTruth { passes, selected }),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct EvaluateError(String);

impl EvaluateError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for EvaluateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for EvaluateError {}

enum Evaluated<'a> {
    Values(Values<'a>),
    Truth(Truth),
    Set(&'a HashSet<String>),
}

pub(crate) fn evaluate<'a>(
    expression: &'a BoundExpression,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Truth, EvaluateError> {
    match evaluate_node(expression, header, record)? {
        Evaluated::Truth(truth) => Ok(truth),
        Evaluated::Values(_) => Err(EvaluateError::new(
            "filter expression does not produce a truth value",
        )),
        Evaluated::Set(_) => Err(EvaluateError::new(
            "value file must be used in an equality comparison",
        )),
    }
}

fn evaluate_node<'a>(
    expression: &'a BoundExpression,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Evaluated<'a>, EvaluateError> {
    match expression {
        BoundExpression::Value(value) => evaluate_value(value, header, record),
        BoundExpression::Unary {
            operator,
            expression,
        } => {
            let values = require_values(evaluate_node(expression, header, record)?)?;
            match operator {
                UnaryOperator::Negate => negate(values).map(Evaluated::Values),
            }
        }
        BoundExpression::Binary {
            operator,
            left,
            right,
        } => {
            let left = evaluate_node(left, header, record)?;
            let right = evaluate_node(right, header, record)?;
            if operator.is_arithmetic() {
                arithmetic(require_values(left)?, require_values(right)?, *operator)
                    .map(Evaluated::Values)
            } else if operator.is_comparison() {
                compare_evaluated(left, right, *operator).map(Evaluated::Truth)
            } else if operator.is_logical() {
                logical(require_truth(left)?, require_truth(right)?, *operator)
                    .map(Evaluated::Truth)
            } else {
                Err(EvaluateError::new("operator is not implemented"))
            }
        }
        BoundExpression::Function { kind, arguments } => {
            let evaluated_arguments = arguments
                .iter()
                .map(|argument| evaluate_node(argument, header, record))
                .collect::<Result<_, _>>()?;
            function::evaluate(*kind, evaluated_arguments, arguments, record)
        }
    }
}

fn evaluate_value<'a>(
    value: &'a BoundValue,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Evaluated<'a>, EvaluateError> {
    let values = match value {
        BoundValue::Number(value) => Values::Site(vec![Atom::Number(*value)]),
        BoundValue::String(value) => Values::Site(vec![Atom::OwnedText(value.clone())]),
        BoundValue::Missing => Values::Site(vec![Atom::Missing]),
        BoundValue::File(values) => return Ok(Evaluated::Set(values)),
        BoundValue::Field(field) => {
            let values = value::read(field, header, record)
                .map_err(|error| EvaluateError::new(error.to_string()))?;
            match &field.subscript {
                Some(subscript) => select::apply(subscript, values, record)?,
                None => values,
            }
        }
    };
    Ok(Evaluated::Values(values))
}

fn compare_evaluated(
    left: Evaluated<'_>,
    right: Evaluated<'_>,
    operator: BinaryOperator,
) -> Result<Truth, EvaluateError> {
    match (left, right) {
        (Evaluated::Values(left), Evaluated::Values(right)) => compare(left, right, operator),
        (Evaluated::Values(values), Evaluated::Set(set))
        | (Evaluated::Set(set), Evaluated::Values(values)) => {
            comparison::compare_set(values, set, operator)
        }
        (Evaluated::Set(_), Evaluated::Set(_)) => {
            Err(EvaluateError::new("cannot compare two value files"))
        }
        _ => Err(EvaluateError::new(
            "truth value cannot be used as a comparison operand",
        )),
    }
}

fn require_values(value: Evaluated<'_>) -> Result<Values<'_>, EvaluateError> {
    match value {
        Evaluated::Values(values) => Ok(values),
        Evaluated::Truth(_) => Err(EvaluateError::new(
            "truth value cannot be used as an arithmetic operand",
        )),
        Evaluated::Set(_) => Err(EvaluateError::new(
            "value file cannot be used as an arithmetic operand",
        )),
    }
}

fn require_truth(value: Evaluated<'_>) -> Result<Truth, EvaluateError> {
    match value {
        Evaluated::Truth(truth) => Ok(truth),
        Evaluated::Values(_) => Err(EvaluateError::new(
            "value cannot be used as a logical operand",
        )),
        Evaluated::Set(_) => Err(EvaluateError::new(
            "value file cannot be used as a logical operand",
        )),
    }
}

fn number(atom: &Atom<'_>) -> Result<Option<f64>, EvaluateError> {
    match atom {
        Atom::Absent | Atom::Missing => Ok(None),
        Atom::Number(value) => Ok(Some(*value)),
        Atom::Flag => Ok(Some(1.0)),
        _ => Err(EvaluateError::new("expected a numeric value")),
    }
}

fn sample_width(samples: &[Vec<Atom<'_>>]) -> usize {
    samples.iter().map(Vec::len).max().unwrap_or(0)
}

trait OperatorKind {
    fn is_arithmetic(&self) -> bool;
    fn is_comparison(&self) -> bool;
    fn is_logical(&self) -> bool;
}

impl OperatorKind for BinaryOperator {
    fn is_arithmetic(&self) -> bool {
        matches!(
            self,
            Self::Add | Self::Subtract | Self::Multiply | Self::Divide | Self::Modulo
        )
    }

    fn is_comparison(&self) -> bool {
        matches!(
            self,
            Self::Equal
                | Self::NotEqual
                | Self::Less
                | Self::LessEqual
                | Self::Greater
                | Self::GreaterEqual
                | Self::Regex
                | Self::NotRegex
        )
    }

    fn is_logical(&self) -> bool {
        matches!(
            self,
            Self::SampleAnd | Self::SiteAnd | Self::SampleOr | Self::SiteOr
        )
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;
    use crate::expression::{bind::bind, syntax::parse};

    fn fixture() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"frequency\">\n\
##INFO=<ID=B,Number=2,Type=Integer,Description=\"binomial counts\">\n\
##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"strand counts\">\n\
##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"forward allele depth\">\n\
##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"reverse allele depth\">\n\
##INFO=<ID=R,Number=R,Type=Integer,Description=\"allele values\">\n\
##INFO=<ID=M,Number=2,Type=Integer,Description=\"values with missing\">\n\
##INFO=<ID=Z,Number=2,Type=Integer,Description=\"zero counts\">\n\
##INFO=<ID=X,Number=1,Type=Integer,Description=\"optional\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">\n\
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"strand counts\">\n\
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"forward allele depth\">\n\
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"reverse allele depth\">\n\
##FORMAT=<ID=ST,Number=1,Type=String,Description=\"status\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            .parse()
            .unwrap();
        let line = b"chr1\t7\trs1;rs2\tA\tC,G\t50\tPASS\tDP=12;AF=0.1,0.9;B=3,5;DP4=8,2,1,5;ADF=8,1,0;ADR=2,5,0;R=1,2,3;M=1,.;Z=0,0\tGT:DP:AD:DP4:ADF:ADR:ST\t0/1:8:4,4,0:8,2,1,5:8,2,0:1,5,0:alpha\t1/2:.:0,3,5:3,1,1,3:0,3,1:0,1,3:BETA\t0/.:2:2,0,.:8,2,1,.:2,0,.:1,0,.:.";
        let raw = vcf::Record::try_from(line.as_slice()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn truth(source: &str, header: &vcf::Header, record: &RecordBuf) -> Truth {
        let expression = bind(parse(source).unwrap(), header).unwrap();
        evaluate(&expression, header, record).unwrap()
    }

    fn genotype_fixture() -> (vcf::Header, RecordBuf) {
        let source = include_str!("../../tests/fixtures/expression-genotypes.vcf");
        let (header, line) = source.trim_end().rsplit_once('\n').unwrap();
        let header: vcf::Header = format!("{header}\n").parse().unwrap();
        let raw = vcf::Record::try_from(line.as_bytes()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn genotype_samples(indices: &[usize]) -> Truth {
        let mut passes = vec![false; 16];
        for index in indices {
            passes[*index] = true;
        }
        Truth::samples(passes)
    }

    #[test]
    fn site_arithmetic_broadcasts_scalars_and_comparisons_match_any_pair() {
        let (header, record) = fixture();
        assert_eq!(
            truth("INFO/DP + 3 >= 15", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("-(INFO/DP - 2) / 2 = -5", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("INFO/DP % 5 = 2", &header, &record),
            Truth::site(true)
        );
        assert_eq!(truth("AF > 0.5", &header, &record), Truth::site(true));
        assert_eq!(truth("AF = 0.1", &header, &record), Truth::site(true));
        assert_eq!(truth("AF < 0", &header, &record), Truth::site(false));
    }

    #[test]
    fn format_arithmetic_retains_per_sample_truth() {
        let (header, record) = fixture();
        assert_eq!(
            truth("FMT/DP + 2 >= 10", &header, &record),
            Truth::samples(vec![true, false, false])
        );
        assert_eq!(
            truth("FMT/AD + FMT/AD > 9", &header, &record),
            Truth::samples(vec![false, true, false])
        );
        assert_eq!(
            truth("FMT/AD > 4", &header, &record),
            Truth::samples(vec![false, true, false])
        );
        assert_eq!(
            truth("FMT/AD > FMT/AD", &header, &record),
            Truth::samples(vec![false, false, false])
        );
    }

    #[test]
    fn equality_has_explicit_missing_and_string_semantics() {
        let (header, record) = fixture();
        assert_eq!(truth("X = \".\"", &header, &record), Truth::site(true));
        assert_eq!(truth("X != \".\"", &header, &record), Truth::site(false));
        assert_eq!(truth("AF != \".\"", &header, &record), Truth::site(true));
        assert_eq!(truth("CHROM = 'chr1'", &header, &record), Truth::site(true));
        assert_eq!(truth("ID = 'rs2'", &header, &record), Truth::site(true));
    }

    #[test]
    fn genotype_named_classes_match_typed_samples() {
        let (header, record) = genotype_fixture();
        assert_eq!(
            truth("GT = 'het'", &header, &record),
            genotype_samples(&[0, 2, 3, 11, 13, 14])
        );
        assert_eq!(
            truth("GT = 'hom'", &header, &record),
            genotype_samples(&[1, 8, 12, 15])
        );
        assert_eq!(
            truth("GT = 'ref'", &header, &record),
            genotype_samples(&[6, 8, 15])
        );
        assert_eq!(
            truth("GT = 'ALT'", &header, &record),
            genotype_samples(&[0, 1, 2, 3, 7, 11, 12, 13, 14])
        );
        assert_eq!(
            truth("GT = 'mis'", &header, &record),
            genotype_samples(&[4, 5, 9, 10])
        );
        assert_eq!(
            truth("GT = 'hap'", &header, &record),
            genotype_samples(&[6, 7])
        );
    }

    #[test]
    fn genotype_symbolic_classes_preserve_case_patterns() {
        let (header, record) = genotype_fixture();
        for source in ["GT = 'AA'", "GT = 'aa'"] {
            assert_eq!(truth(source, &header, &record), genotype_samples(&[1, 12]));
        }
        for source in ["GT = 'Aa'", "GT = 'aA'"] {
            assert_eq!(truth(source, &header, &record), genotype_samples(&[2, 13]));
        }
        for (source, expected) in [
            ("GT = 'RR'", &[8, 15][..]),
            ("GT = 'RA'", &[0, 3, 11, 14][..]),
            ("GT = 'AR'", &[0, 3, 11, 14][..]),
            ("GT = 'R'", &[6][..]),
            ("GT = 'A'", &[7][..]),
        ] {
            assert_eq!(truth(source, &header, &record), genotype_samples(expected));
        }
    }

    #[test]
    fn genotype_exact_spelling_regex_and_inequality_are_typed() {
        let (header, record) = genotype_fixture();
        for (source, expected) in [
            ("GT = './.'", 4),
            ("GT = '.|.'", 5),
            ("GT = '.'", 10),
            ("GT = '0/1/2'", 3),
            ("GT = '0/1|2'", 11),
            ("GT ~ '^1[|/]0$'", 0),
        ] {
            assert_eq!(
                truth(source, &header, &record),
                genotype_samples(&[expected])
            );
        }
        assert_eq!(
            truth("GT != 'mis'", &header, &record),
            genotype_samples(&[0, 1, 2, 3, 6, 7, 8, 11, 12, 13, 14, 15])
        );

        let expression = bind(parse("GT > 'het'").unwrap(), &header).unwrap();
        assert!(evaluate(&expression, &header, &record).is_err());
    }

    #[test]
    fn incompatible_arithmetic_vectors_fail() {
        let (header, record) = fixture();
        let expression = bind(parse("AF + R > 0").unwrap(), &header).unwrap();
        assert!(evaluate(&expression, &header, &record).is_err());
    }

    #[test]
    fn logical_operators_keep_sample_and_site_semantics_distinct() {
        let (header, record) = fixture();
        let first = "FMT/DP >= 8";
        let second = "FMT/AD > 4";
        assert_eq!(
            truth(&format!("{first} & {second}"), &header, &record),
            Truth::samples(vec![false, false, false])
        );
        assert_eq!(
            truth(&format!("{first} && {second}"), &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth(&format!("{first} | {second}"), &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth(&format!("{first} || {second}"), &header, &record),
            Truth::samples(vec![true, true, true])
        );
    }

    #[test]
    fn site_truth_controls_mixed_logical_sample_masks() {
        let (header, record) = fixture();
        let sample = "FMT/DP >= 8";
        for operator in ["&", "&&"] {
            assert_eq!(
                truth(&format!("QUAL > 40 {operator} {sample}"), &header, &record,),
                Truth::samples(vec![true, false, false])
            );
            assert_eq!(
                truth(&format!("QUAL < 40 {operator} {sample}"), &header, &record,),
                Truth::samples(vec![false, false, false])
            );
        }
        assert_eq!(
            truth("QUAL > 40 | FMT/DP < 0", &header, &record),
            Truth::with_samples(
                true,
                vec![false, false, false],
                vec![true, true, true].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("QUAL > 40 || FMT/DP < 0", &header, &record),
            Truth::samples(vec![true, true, true])
        );
        assert_eq!(
            truth("QUAL < 40 || FMT/DP >= 8", &header, &record),
            Truth::samples(vec![true, false, false])
        );
    }

    #[test]
    fn regex_comparisons_support_vectors_samples_and_case_suffixes() {
        let (header, record) = fixture();
        assert_eq!(
            truth("ID ~ '^RS[12]$/i'", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("ID ~ '^RS3$/i'", &header, &record),
            Truth::site(false)
        );
        assert_eq!(truth("ID !~ '^rs1$'", &header, &record), Truth::site(true));
        assert_eq!(
            truth("FMT/ST ~ '^beta$/i'", &header, &record),
            Truth::samples(vec![false, true, false])
        );
    }

    #[test]
    fn invalid_or_non_string_regexes_fail() {
        let (header, record) = fixture();
        for source in ["ID ~ '['", "INFO/DP ~ '1'"] {
            let expression = bind(parse(source).unwrap(), &header).unwrap();
            assert!(evaluate(&expression, &header, &record).is_err(), "{source}");
        }
    }

    #[test]
    fn info_subscripts_select_indices_ranges_and_missing_values() {
        let (header, record) = fixture();
        assert_eq!(truth("AF[0] = 0.1", &header, &record), Truth::site(true));
        assert_eq!(truth("AF[1] = 0.9", &header, &record), Truth::site(true));
        assert_eq!(truth("AF[0-1] > 0.5", &header, &record), Truth::site(true));
        assert_eq!(truth("AF[2] = '.'", &header, &record), Truth::site(true));
        assert_eq!(
            truth("INFO/DP[9] = 12", &header, &record),
            Truth::site(true)
        );
    }

    #[test]
    fn format_subscripts_keep_sample_and_value_axes_distinct() {
        let (header, record) = fixture();
        assert_eq!(
            truth("FMT/DP[0] >= 8", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/DP[1-] > 1", &header, &record),
            Truth::selected_samples(
                vec![false, false, true],
                vec![false, true, true].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/AD[:1] >= 4", &header, &record),
            Truth::samples(vec![true, false, false])
        );
        assert_eq!(
            truth("FMT/AD[1:2] > 4", &header, &record),
            Truth::selected_samples(
                vec![false, true, false],
                vec![false, true, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/AD[0,2:0-1] > 3", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, true].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/DP[0] = '.'", &header, &record),
            Truth::selected_samples(
                vec![false, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/DP[1] = '.'", &header, &record),
            Truth::selected_samples(
                vec![false, true, false],
                vec![false, true, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/DP[0] + 2 >= 10", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/ST[1] ~ '^beta$/i'", &header, &record),
            Truth::selected_samples(
                vec![false, true, false],
                vec![false, true, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("FMT/DP[0] >= 8 | FMT/AD[1:2] > 4", &header, &record,),
            Truth::selected_samples(
                vec![true, true, false],
                vec![true, true, false].into_boxed_slice()
            )
        );
    }

    #[test]
    fn genotype_subscripts_select_called_alleles_per_sample() {
        let (header, record) = fixture();
        assert_eq!(
            truth("FMT/AD[GT] > 4", &header, &record),
            Truth::samples(vec![false, true, false])
        );
        assert_eq!(
            truth("FMT/AD[0:GT] > 3", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
    }

    #[test]
    fn sample_files_are_resolved_during_binding() {
        use std::io::Write;

        let (header, record) = fixture();
        let mut samples = tempfile::NamedTempFile::new_in(".").unwrap();
        writeln!(samples, "S2\nS3").unwrap();
        let source = format!("FMT/DP[@{}] > 1", samples.path().display());
        let expression = bind(parse(&source).unwrap(), &header).unwrap();
        drop(samples);
        assert_eq!(
            evaluate(&expression, &header, &record).unwrap(),
            Truth::selected_samples(
                vec![false, false, true],
                vec![false, true, true].into_boxed_slice()
            )
        );
    }

    #[test]
    fn global_numeric_functions_reduce_sites_and_selected_samples() {
        let (header, record) = fixture();
        for source in [
            "MAX(AF) = 0.9",
            "MIN(AF) = 0.1",
            "SUM(AF) = 1",
            "MEAN(AF) = 0.5",
            "AVG(AF) = 0.5",
            "MEDIAN(AF) = 0.5",
            "STDEV(AF) > 0.39 & STDEV(AF) < 0.41",
            "STDEV(INFO/DP) = 0",
            "SUM(FMT/AD) = 18",
            "MAX(FMT/AD) = 5",
            "SUM(FMT/AD[0]) = 8",
            "SUM(X) = '.'",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
    }

    #[test]
    fn sample_numeric_functions_reduce_each_selected_sample() {
        let (header, record) = fixture();
        assert_eq!(
            truth("SMPL_SUM(FMT/AD) = 8", &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth("sMEAN(FMT/AD) > 2", &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth("SMPL_MEDIAN(FMT/AD) > 2", &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth("sSTDEV(FMT/AD) > 2", &header, &record),
            Truth::samples(vec![false, true, false])
        );
        assert_eq!(
            truth("SMPL_SUM(FMT/AD[0]) = 8", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("SMPL_SUM(FMT/DP) = '.'", &header, &record),
            Truth::samples(vec![false, true, false])
        );
    }

    #[test]
    fn elementwise_functions_preserve_vectors_samples_and_missing_values() {
        let (header, record) = fixture();
        assert_eq!(
            truth("ABS(-AF[0]) = 0.1", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("ABS(FMT/DP - 10) = 2", &header, &record),
            Truth::samples(vec![true, false, false])
        );
        assert_eq!(truth("ABS(X) = '.'", &header, &record), Truth::site(true));
        assert_eq!(
            truth("STRLEN(CHROM) = 4", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("STRLEN('abc') = 3", &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth("STRLEN('.') = 0", &header, &record),
            Truth::site(true)
        );
    }

    #[test]
    fn phred_maps_numeric_values_and_preserves_shape() {
        let (header, record) = fixture();
        for source in ["PHRED(0.01) = 20", "PHRED(AF[0]) = 10", "PHRED(X) = '.'"] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
        assert_eq!(
            truth("PHRED(FMT/DP[0]) < -9", &header, &record),
            Truth::selected_samples(
                vec![true, false, false],
                vec![true, false, false].into_boxed_slice()
            )
        );
        assert_eq!(
            truth("PHRED(FMT/DP[1]) = '.'", &header, &record),
            Truth::selected_samples(
                vec![false, true, false],
                vec![false, true, false].into_boxed_slice()
            )
        );
    }

    #[test]
    fn binomial_supports_site_sample_and_explicit_count_pairs() {
        let (header, record) = fixture();
        for source in [
            "BINOM(B) > 0.7 & BINOM(B) < 0.8",
            "BINOM(B[0], B[1]) > 0.7",
            "BINOM(Z) = '.'",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
        assert_eq!(
            truth("BINOM(FMT/AD) > 0.7", &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth("BINOM(FMT/AD[:0], FMT/AD[:1]) > 0.2", &header, &record),
            Truth::samples(vec![true, true, true])
        );
    }

    #[test]
    fn fisher_supports_four_counts_paired_vectors_and_genotype_selection() {
        let (header, record) = fixture();
        for source in [
            "FISHER(INFO/DP4) > 0.03 & FISHER(INFO/DP4) < 0.04",
            "FISHER(INFO/ADF[0,1], INFO/ADR[0,1]) < 0.04",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
        for source in [
            "FISHER(FMT/DP4) < 0.04",
            "FISHER(FMT/ADF[:0,1], FMT/ADR[:0,1]) < 0.04",
            "FISHER(FMT/ADF, FMT/ADR) < 0.04",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::samples(vec![true, false, false]),
                "{source}"
            );
        }
    }

    #[test]
    fn fisher_distinguishes_unindexed_number_r_from_explicit_pairs() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"forward allele depth\">\n\
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"reverse allele depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(
            b"chr1\t7\t.\tA\tC\t.\tPASS\t.\tGT:ADF:ADR\t0/0:8,2:1,5".as_slice(),
        )
        .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        assert_eq!(
            truth("FISHER(FMT/ADF, FMT/ADR) = 1", &header, &record),
            Truth::samples(vec![true])
        );
        assert_eq!(
            truth(
                "FISHER(FMT/ADF[:0,1], FMT/ADR[:0,1]) < 0.04",
                &header,
                &record
            ),
            Truth::samples(vec![true])
        );
    }

    #[test]
    fn value_files_match_id_and_format_strings_as_sets() {
        let (header, record) = fixture();
        let mut file = tempfile::NamedTempFile::new_in(".").unwrap();
        writeln!(file, "rs2 annotation\nBETA").unwrap();
        let path = file.path().display();
        assert_eq!(
            truth(&format!("ID=@{path}"), &header, &record),
            Truth::site(true)
        );
        assert_eq!(
            truth(&format!("ID!=@{path}"), &header, &record),
            Truth::site(false)
        );
        assert_eq!(
            truth(&format!("FMT/ST=@{path}"), &header, &record),
            Truth::samples(vec![false, true, false])
        );
        assert_eq!(
            truth(&format!("FMT/ST!=@{path}"), &header, &record),
            Truth::samples(vec![true, false, true])
        );
        assert!(bind(parse(&format!("ID>@{path}")).unwrap(), &header).is_err());
    }

    #[test]
    fn missing_value_files_fail_during_binding() {
        let (header, _) = fixture();
        let directory = tempfile::tempdir_in(".").unwrap();
        let source = format!("ID=@{}", directory.path().join("missing").display());
        assert!(bind(parse(&source).unwrap(), &header).is_err());
    }

    #[test]
    fn count_functions_distinguish_site_slots_sample_values_and_absence() {
        let (header, record) = fixture();
        for source in [
            "COUNT(AF) = 2",
            "COUNT(M) = 2",
            "COUNT(M[0-1]) = 2",
            "COUNT(M[1]) = 0",
            "COUNT(M[3-]) = 0",
            "COUNT(X) = 0",
            "COUNT(AF[2]) = 0",
            "COUNT(INFO/DP[9]) = 1",
            "COUNT(FMT/AD) = 8",
            "COUNT(FMT/AD[:2]) = 2",
            "COUNT(FMT/DP) = 2",
            "COUNT(FMT/AD[0]) = 3",
            "COUNT(FMT/DP >= 8) = 1",
            "COUNT(FMT/DP[0] >= 8) = 1",
            "COUNT(SUM(X)) = 0",
            "COUNT(SUM(M)) = 1",
            "COUNT(SUM(M[1])) = 0",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
        assert_eq!(
            truth("SMPL_COUNT(FMT/AD) = 3", &header, &record),
            Truth::samples(vec![true, true, false])
        );
        assert_eq!(
            truth("sCOUNT(FMT/DP) = 1", &header, &record),
            Truth::samples(vec![true, false, true])
        );
        assert_eq!(
            truth("sCOUNT(FMT/AD[:2]) = 1", &header, &record),
            Truth::samples(vec![true, true, false])
        );
    }

    #[test]
    fn passing_sample_functions_count_and_fraction_all_samples() {
        let (header, record) = fixture();
        for source in ["N_PASS(FMT/DP >= 8) = 1", "F_PASS(FMT/DP >= 8) = 1 / 3"] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
    }

    #[test]
    fn passing_sample_fraction_uses_the_expression_selection() {
        let (header, record) = fixture();
        for source in [
            "N_PASS(FMT/DP[0] >= 8) = 1",
            "F_PASS(FMT/DP[0] >= 8) = 1",
            "N_PASS(FMT/DP[0] >= 8 | FMT/DP[2] > 2) = 1",
            "F_PASS(FMT/DP[0] >= 8 | FMT/DP[2] > 2) = 0.5",
        ] {
            assert_eq!(
                truth(source, &header, &record),
                Truth::site(true),
                "{source}"
            );
        }
    }
}
