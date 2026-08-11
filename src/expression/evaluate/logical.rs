use crate::expression::syntax::BinaryOperator;

use super::{EvaluateError, Truth};

pub(super) fn logical(
    left: Truth,
    right: Truth,
    operator: BinaryOperator,
) -> Result<Truth, EvaluateError> {
    let sample_count = sample_count(&left, &right)?;
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

fn sample_count(left: &Truth, right: &Truth) -> Result<Option<usize>, EvaluateError> {
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
