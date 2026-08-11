use regex::{Regex, RegexBuilder};

use crate::expression::{
    syntax::BinaryOperator,
    value::{self, Atom, Values},
};

use super::{EvaluateError, Truth, number, sample_width};

pub(super) fn compare(
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
            if left.values.len() != right.values.len() {
                return Err(EvaluateError::new(format!(
                    "incompatible sample counts in comparison: {} vs {}",
                    left.values.len(),
                    right.values.len()
                )));
            }
            let left_width = sample_width(&left.values);
            let right_width = sample_width(&right.values);
            if left_width != right_width {
                return Err(EvaluateError::new(format!(
                    "incompatible per-sample value counts in comparison: {left_width} vs {right_width}"
                )));
            }
            left.values
                .iter()
                .zip(&right.values)
                .zip(left.selected.iter().zip(&right.selected))
                .map(|((left, right), (left_selected, right_selected))| {
                    if *left_selected && *right_selected {
                        compare_pairs(left, right, left_width, operator)
                    } else {
                        Ok(false)
                    }
                })
                .collect::<Result<_, _>>()
                .map(Truth::samples)
        }
        (Values::Samples(samples), Values::Site(site)) => samples
            .values
            .iter()
            .zip(&samples.selected)
            .map(|(sample, selected)| {
                if *selected {
                    compare_cross(sample, &site, operator)
                } else {
                    Ok(false)
                }
            })
            .collect::<Result<_, _>>()
            .map(Truth::samples),
        (Values::Site(site), Values::Samples(samples)) => samples
            .values
            .iter()
            .zip(&samples.selected)
            .map(|(sample, selected)| {
                if *selected {
                    compare_cross(&site, sample, operator)
                } else {
                    Ok(false)
                }
            })
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
        Atom::Absent | Atom::Missing => ".",
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
            .values
            .iter()
            .zip(&samples.selected)
            .map(|(values, selected)| {
                if *selected {
                    regex_values(values, &regex, negate)
                } else {
                    Ok(false)
                }
            })
            .collect::<Result<_, _>>()
            .map(Truth::samples),
    }
}

fn regex_values(values: &[Atom<'_>], regex: &Regex, negate: bool) -> Result<bool, EvaluateError> {
    for value in values {
        let passes = match value {
            Atom::Absent | Atom::Missing => regex.is_match(".") != negate,
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
        (left, right) if is_missing(left) && is_missing(right) => {
            Ok(operator == BinaryOperator::Equal)
        }
        (left, right) if is_missing(left) || is_missing(right) => {
            Ok(operator == BinaryOperator::NotEqual)
        }
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

fn is_missing(atom: &Atom<'_>) -> bool {
    matches!(atom, Atom::Absent | Atom::Missing)
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

fn text<'a>(atom: &'a Atom<'_>) -> Result<&'a str, EvaluateError> {
    match atom {
        Atom::Text(value) => Ok(value),
        Atom::OwnedText(value) => Ok(value),
        _ => Err(EvaluateError::new("expected a string value")),
    }
}
