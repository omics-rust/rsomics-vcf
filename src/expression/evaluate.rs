use std::fmt;

use noodles_vcf::{self as vcf, variant::RecordBuf};
use regex::{Regex, RegexBuilder};

use super::{
    bind::{BoundExpression, BoundValue},
    syntax::{BinaryOperator, UnaryOperator},
    value::{self, Atom, Values},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Truth {
    pub site: bool,
    pub samples: Option<Vec<bool>>,
}

impl Truth {
    fn site(value: bool) -> Self {
        Self {
            site: value,
            samples: None,
        }
    }

    fn samples(values: Vec<bool>) -> Self {
        Self::with_samples(values.iter().any(|value| *value), values)
    }

    fn with_samples(site: bool, values: Vec<bool>) -> Self {
        Self {
            site,
            samples: Some(values),
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
}

pub(crate) fn evaluate(
    expression: &BoundExpression,
    header: &vcf::Header,
    record: &RecordBuf,
) -> Result<Truth, EvaluateError> {
    match evaluate_node(expression, header, record)? {
        Evaluated::Truth(truth) => Ok(truth),
        Evaluated::Values(_) => Err(EvaluateError::new(
            "filter expression does not produce a truth value",
        )),
    }
}

fn evaluate_node<'a>(
    expression: &BoundExpression,
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
                compare(require_values(left)?, require_values(right)?, *operator)
                    .map(Evaluated::Truth)
            } else if operator.is_logical() {
                logical(require_truth(left)?, require_truth(right)?, *operator)
                    .map(Evaluated::Truth)
            } else {
                Err(EvaluateError::new("operator is not implemented"))
            }
        }
        BoundExpression::Function { .. } => Err(EvaluateError::new(
            "function expression evaluation is not implemented",
        )),
    }
}

fn evaluate_value<'a>(
    value: &BoundValue,
    header: &vcf::Header,
    record: &'a RecordBuf,
) -> Result<Evaluated<'a>, EvaluateError> {
    let values = match value {
        BoundValue::Number(value) => Values::Site(vec![Atom::Number(*value)]),
        BoundValue::String(value) => Values::Site(vec![Atom::OwnedText(value.clone())]),
        BoundValue::Missing => Values::Site(vec![Atom::Missing]),
        BoundValue::File(path) => {
            return Err(EvaluateError::new(format!(
                "value file {} is not loaded",
                path.display()
            )));
        }
        BoundValue::Field(field) => {
            if field.subscript.is_some() {
                return Err(EvaluateError::new(
                    "field subscript evaluation is not implemented",
                ));
            }
            value::read(field, header, record)
                .map_err(|error| EvaluateError::new(error.to_string()))?
        }
    };
    Ok(Evaluated::Values(values))
}

fn require_values(value: Evaluated<'_>) -> Result<Values<'_>, EvaluateError> {
    match value {
        Evaluated::Values(values) => Ok(values),
        Evaluated::Truth(_) => Err(EvaluateError::new(
            "truth value cannot be used as an arithmetic operand",
        )),
    }
}

fn require_truth(value: Evaluated<'_>) -> Result<Truth, EvaluateError> {
    match value {
        Evaluated::Truth(truth) => Ok(truth),
        Evaluated::Values(_) => Err(EvaluateError::new(
            "value cannot be used as a logical operand",
        )),
    }
}

fn logical(left: Truth, right: Truth, operator: BinaryOperator) -> Result<Truth, EvaluateError> {
    let sample_count = logical_sample_count(&left, &right)?;
    match operator {
        BinaryOperator::SampleAnd | BinaryOperator::SiteAnd => {
            if !left.site || !right.site {
                return Ok(sample_count
                    .map(|count| Truth::with_samples(false, vec![false; count]))
                    .unwrap_or_else(|| Truth::site(false)));
            }
            match (left.samples, right.samples) {
                (None, None) => Ok(Truth::site(true)),
                (Some(samples), None) | (None, Some(samples)) => {
                    Ok(Truth::with_samples(true, samples))
                }
                (Some(left), Some(right)) => {
                    let samples = left
                        .iter()
                        .zip(right)
                        .map(|(left, right)| match operator {
                            BinaryOperator::SampleAnd => *left && right,
                            BinaryOperator::SiteAnd => *left || right,
                            _ => unreachable!(),
                        })
                        .collect();
                    Ok(match operator {
                        BinaryOperator::SampleAnd => Truth::samples(samples),
                        BinaryOperator::SiteAnd => Truth::with_samples(true, samples),
                        _ => unreachable!(),
                    })
                }
            }
        }
        BinaryOperator::SampleOr | BinaryOperator::SiteOr => {
            if !left.site && !right.site {
                return Ok(sample_count
                    .map(|count| Truth::with_samples(false, vec![false; count]))
                    .unwrap_or_else(|| Truth::site(false)));
            }
            match (left.samples, right.samples) {
                (None, None) => Ok(Truth::site(true)),
                (Some(samples), None) => {
                    if operator == BinaryOperator::SiteOr && right.site {
                        Ok(Truth::with_samples(true, vec![true; samples.len()]))
                    } else {
                        Ok(Truth::with_samples(true, samples))
                    }
                }
                (None, Some(samples)) => {
                    if operator == BinaryOperator::SiteOr && left.site {
                        Ok(Truth::with_samples(true, vec![true; samples.len()]))
                    } else {
                        Ok(Truth::with_samples(true, samples))
                    }
                }
                (Some(left), Some(right)) => {
                    if operator == BinaryOperator::SiteOr {
                        Ok(Truth::with_samples(true, vec![true; left.len()]))
                    } else {
                        let samples = left
                            .into_iter()
                            .zip(right)
                            .map(|(left, right)| left || right)
                            .collect();
                        Ok(Truth::with_samples(true, samples))
                    }
                }
            }
        }
        _ => Err(EvaluateError::new("operator is not logical")),
    }
}

fn logical_sample_count(left: &Truth, right: &Truth) -> Result<Option<usize>, EvaluateError> {
    match (&left.samples, &right.samples) {
        (Some(left), Some(right)) if left.len() != right.len() => Err(EvaluateError::new(format!(
            "incompatible sample counts in logical expression: {} vs {}",
            left.len(),
            right.len()
        ))),
        (Some(samples), _) | (_, Some(samples)) => Ok(Some(samples.len())),
        (None, None) => Ok(None),
    }
}

fn negate(values: Values<'_>) -> Result<Values<'_>, EvaluateError> {
    map_unary(values, |atom| match number(atom)? {
        Some(value) => Ok(Atom::Number(-value)),
        None => Ok(Atom::Missing),
    })
}

fn map_unary<'a>(
    values: Values<'a>,
    operation: impl Fn(&Atom<'a>) -> Result<Atom<'a>, EvaluateError>,
) -> Result<Values<'a>, EvaluateError> {
    match values {
        Values::Site(values) => values
            .iter()
            .map(&operation)
            .collect::<Result<_, _>>()
            .map(Values::Site),
        Values::Samples(samples) => samples
            .iter()
            .map(|values| values.iter().map(&operation).collect::<Result<_, _>>())
            .collect::<Result<_, _>>()
            .map(Values::Samples),
    }
}

fn arithmetic<'a>(
    left: Values<'a>,
    right: Values<'a>,
    operator: BinaryOperator,
) -> Result<Values<'a>, EvaluateError> {
    match (left, right) {
        (Values::Site(left), Values::Site(right)) => {
            arithmetic_vectors(&left, &right, operator).map(Values::Site)
        }
        (Values::Samples(left), Values::Samples(right)) => {
            if left.len() != right.len() {
                return Err(EvaluateError::new(format!(
                    "incompatible sample counts in arithmetic: {} vs {}",
                    left.len(),
                    right.len()
                )));
            }
            let left_width = sample_width(&left);
            let right_width = sample_width(&right);
            let width = broadcast_width(left_width, right_width)?;
            left.iter()
                .zip(&right)
                .map(|(left, right)| {
                    arithmetic_sample(left, right, left_width, right_width, width, operator)
                })
                .collect::<Result<_, _>>()
                .map(Values::Samples)
        }
        (Values::Samples(samples), Values::Site(site)) => {
            arithmetic_samples_site(samples, site, operator, false)
        }
        (Values::Site(site), Values::Samples(samples)) => {
            arithmetic_samples_site(samples, site, operator, true)
        }
    }
}

fn arithmetic_vectors<'a>(
    left: &[Atom<'a>],
    right: &[Atom<'a>],
    operator: BinaryOperator,
) -> Result<Vec<Atom<'a>>, EvaluateError> {
    let width = broadcast_width(left.len(), right.len())?;
    (0..width)
        .map(|index| {
            let left = broadcast_atom(left, index);
            let right = broadcast_atom(right, index);
            arithmetic_atoms(left, right, operator)
        })
        .collect()
}

fn arithmetic_sample<'a>(
    left: &[Atom<'a>],
    right: &[Atom<'a>],
    left_width: usize,
    right_width: usize,
    width: usize,
    operator: BinaryOperator,
) -> Result<Vec<Atom<'a>>, EvaluateError> {
    (0..width)
        .map(|index| {
            let left = sample_atom(left, left_width, index);
            let right = sample_atom(right, right_width, index);
            arithmetic_atoms(left, right, operator)
        })
        .collect()
}

fn arithmetic_samples_site<'a>(
    samples: Vec<Vec<Atom<'a>>>,
    site: Vec<Atom<'a>>,
    operator: BinaryOperator,
    site_first: bool,
) -> Result<Values<'a>, EvaluateError> {
    if site.len() != 1 {
        return Err(EvaluateError::new(format!(
            "sample arithmetic requires a scalar site value, found {} values",
            site.len()
        )));
    }
    let site = &site[0];
    samples
        .iter()
        .map(|sample| {
            sample
                .iter()
                .map(|atom| {
                    if site_first {
                        arithmetic_atoms(site, atom, operator)
                    } else {
                        arithmetic_atoms(atom, site, operator)
                    }
                })
                .collect::<Result<_, _>>()
        })
        .collect::<Result<_, _>>()
        .map(Values::Samples)
}

fn arithmetic_atoms<'a>(
    left: &Atom<'a>,
    right: &Atom<'a>,
    operator: BinaryOperator,
) -> Result<Atom<'a>, EvaluateError> {
    let (Some(left), Some(right)) = (number(left)?, number(right)?) else {
        return Ok(Atom::Missing);
    };
    let value = match operator {
        BinaryOperator::Add => left + right,
        BinaryOperator::Subtract => left - right,
        BinaryOperator::Multiply => left * right,
        BinaryOperator::Divide => left / right,
        BinaryOperator::Modulo => {
            let left = left as i64;
            let right = right as i64;
            if right == 0 {
                return Err(EvaluateError::new("modulo by zero"));
            }
            (left % right) as f64
        }
        _ => return Err(EvaluateError::new("operator is not arithmetic")),
    };
    Ok(Atom::Number(value))
}

fn compare(
    left: Values<'_>,
    right: Values<'_>,
    operator: BinaryOperator,
) -> Result<Truth, EvaluateError> {
    if matches!(operator, BinaryOperator::Regex | BinaryOperator::NotRegex) {
        return compare_regex(left, right, operator);
    }
    match (left, right) {
        (Values::Site(left), Values::Site(right)) => {
            compare_cross(&left, &right, operator).map(Truth::site)
        }
        (Values::Samples(left), Values::Samples(right)) => {
            if left.len() != right.len() {
                return Err(EvaluateError::new(format!(
                    "incompatible sample counts in comparison: {} vs {}",
                    left.len(),
                    right.len()
                )));
            }
            let left_width = sample_width(&left);
            let right_width = sample_width(&right);
            if left_width != right_width {
                return Err(EvaluateError::new(format!(
                    "incompatible per-sample value counts in comparison: {left_width} vs {right_width}"
                )));
            }
            left.iter()
                .zip(&right)
                .map(|(left, right)| compare_pairs(left, right, left_width, operator))
                .collect::<Result<_, _>>()
                .map(Truth::samples)
        }
        (Values::Samples(samples), Values::Site(site)) => samples
            .iter()
            .map(|sample| compare_cross(sample, &site, operator))
            .collect::<Result<_, _>>()
            .map(Truth::samples),
        (Values::Site(site), Values::Samples(samples)) => samples
            .iter()
            .map(|sample| compare_cross(&site, sample, operator))
            .collect::<Result<_, _>>()
            .map(Truth::samples),
    }
}

fn compare_regex(
    values: Values<'_>,
    pattern: Values<'_>,
    operator: BinaryOperator,
) -> Result<Truth, EvaluateError> {
    let Values::Site(pattern) = pattern else {
        return Err(EvaluateError::new("regex pattern must be a site value"));
    };
    if pattern.len() != 1 {
        return Err(EvaluateError::new(format!(
            "regex pattern must be scalar, found {} values",
            pattern.len()
        )));
    }
    let pattern = match &pattern[0] {
        Atom::Text(value) => *value,
        Atom::OwnedText(value) => value,
        Atom::Missing => ".",
        _ => return Err(EvaluateError::new("regex pattern must be a string")),
    };
    let (pattern, case_insensitive) = pattern
        .strip_suffix("/i")
        .map(|pattern| (pattern, true))
        .unwrap_or((pattern, false));
    let regex = RegexBuilder::new(pattern)
        .case_insensitive(case_insensitive)
        .build()
        .map_err(|error| EvaluateError::new(format!("invalid regex: {error}")))?;
    let negate = operator == BinaryOperator::NotRegex;
    match values {
        Values::Site(values) => regex_values(&values, &regex, negate).map(Truth::site),
        Values::Samples(samples) => samples
            .iter()
            .map(|values| regex_values(values, &regex, negate))
            .collect::<Result<_, _>>()
            .map(Truth::samples),
    }
}

fn regex_values(values: &[Atom<'_>], regex: &Regex, negate: bool) -> Result<bool, EvaluateError> {
    for value in values {
        let passes = match value {
            Atom::Missing => regex.is_match(".") != negate,
            Atom::Text(value) => regex.is_match(value) != negate,
            Atom::OwnedText(value) => regex.is_match(value) != negate,
            Atom::Filter(filter) => match filter {
                value::Filter::Pass => regex.is_match("PASS") != negate,
                value::Filter::Missing => regex.is_match(".") != negate,
                value::Filter::Failed(filters) => {
                    filters.iter().any(|value| regex.is_match(value) != negate)
                }
            },
            _ => return Err(EvaluateError::new("regex operand must be a string")),
        };
        if passes {
            return Ok(true);
        }
    }
    Ok(false)
}

fn compare_cross(
    left: &[Atom<'_>],
    right: &[Atom<'_>],
    operator: BinaryOperator,
) -> Result<bool, EvaluateError> {
    for left in left {
        for right in right {
            if compare_atoms(left, right, operator)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn compare_pairs(
    left: &[Atom<'_>],
    right: &[Atom<'_>],
    width: usize,
    operator: BinaryOperator,
) -> Result<bool, EvaluateError> {
    for index in 0..width {
        let left = left.get(index).unwrap_or(&Atom::Missing);
        let right = right.get(index).unwrap_or(&Atom::Missing);
        if compare_atoms(left, right, operator)? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn compare_atoms(
    left: &Atom<'_>,
    right: &Atom<'_>,
    operator: BinaryOperator,
) -> Result<bool, EvaluateError> {
    match (left, right) {
        (Atom::Missing, Atom::Missing) => Ok(operator == BinaryOperator::Equal),
        (Atom::Missing, _) | (_, Atom::Missing) => Ok(operator == BinaryOperator::NotEqual),
        _ => match (number(left), number(right)) {
            (Ok(Some(left)), Ok(Some(right))) => compare_numbers(left, right, operator),
            (Err(_), Err(_)) => {
                let left = text(left)?;
                let right = text(right)?;
                match operator {
                    BinaryOperator::Equal => Ok(left == right),
                    BinaryOperator::NotEqual => Ok(left != right),
                    _ => Err(EvaluateError::new(
                        "strings support only equality and regex comparisons",
                    )),
                }
            }
            _ => Err(EvaluateError::new("cannot compare strings and numbers")),
        },
    }
}

fn compare_numbers(left: f64, right: f64, operator: BinaryOperator) -> Result<bool, EvaluateError> {
    if left > 16_777_216.0 || right > 16_777_216.0 {
        return compare_ordered(left, right, operator);
    }
    compare_ordered(left as f32, right as f32, operator)
}

fn compare_ordered<T: PartialEq + PartialOrd>(
    left: T,
    right: T,
    operator: BinaryOperator,
) -> Result<bool, EvaluateError> {
    match operator {
        BinaryOperator::Equal => Ok(left == right),
        BinaryOperator::NotEqual => Ok(left != right),
        BinaryOperator::Less => Ok(left < right),
        BinaryOperator::LessEqual => Ok(left <= right),
        BinaryOperator::Greater => Ok(left > right),
        BinaryOperator::GreaterEqual => Ok(left >= right),
        _ => Err(EvaluateError::new("operator is not a comparison")),
    }
}

fn number(atom: &Atom<'_>) -> Result<Option<f64>, EvaluateError> {
    match atom {
        Atom::Missing => Ok(None),
        Atom::Number(value) => Ok(Some(*value)),
        Atom::Flag => Ok(Some(1.0)),
        _ => Err(EvaluateError::new("expected a numeric value")),
    }
}

fn text<'a>(atom: &'a Atom<'_>) -> Result<&'a str, EvaluateError> {
    match atom {
        Atom::Text(value) => Ok(value),
        Atom::OwnedText(value) => Ok(value),
        _ => Err(EvaluateError::new("expected a string value")),
    }
}

fn broadcast_width(left: usize, right: usize) -> Result<usize, EvaluateError> {
    if left == right || left == 1 || right == 1 {
        Ok(left.max(right))
    } else {
        Err(EvaluateError::new(format!(
            "incompatible arithmetic vector lengths: {left} vs {right}"
        )))
    }
}

fn broadcast_atom<'a, 'v>(values: &'v [Atom<'a>], index: usize) -> &'v Atom<'a> {
    if values.len() <= 1 {
        values.first().unwrap_or(&Atom::Missing)
    } else {
        &values[index]
    }
}

fn sample_atom<'a, 'v>(values: &'v [Atom<'a>], width: usize, index: usize) -> &'v Atom<'a> {
    if width == 1 {
        values.first().unwrap_or(&Atom::Missing)
    } else {
        values.get(index).unwrap_or(&Atom::Missing)
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
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;
    use crate::expression::{bind::bind, syntax::parse};

    fn fixture() -> (vcf::Header, RecordBuf) {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"frequency\">\n\
##INFO=<ID=R,Number=R,Type=Integer,Description=\"allele values\">\n\
##INFO=<ID=X,Number=1,Type=Integer,Description=\"optional\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">\n\
##FORMAT=<ID=ST,Number=1,Type=String,Description=\"status\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            .parse()
            .unwrap();
        let line = b"chr1\t7\trs1;rs2\tA\tC,G\t50\tPASS\tDP=12;AF=0.1,0.9;R=1,2,3\tGT:DP:AD:ST\t0/1:8:4,4,0:alpha\t1/2:.:0,3,5:BETA\t0/.:2:2,0,.:.";
        let raw = vcf::Record::try_from(line.as_slice()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, record)
    }

    fn truth(source: &str, header: &vcf::Header, record: &RecordBuf) -> Truth {
        let expression = bind(parse(source).unwrap(), header).unwrap();
        evaluate(&expression, header, record).unwrap()
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
            Truth {
                site: true,
                samples: Some(vec![false, false, false]),
            }
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
}
