use crate::expression::syntax::BinaryOperator;

use super::{EvaluateError, Truth};

pub(super) fn logical(
    left: Truth,
    right: Truth,
    operator: BinaryOperator,
) -> Result<Truth, EvaluateError> {
    let selected = selected_union(&left, &right)?;
    match operator {
        BinaryOperator::SampleAnd | BinaryOperator::SiteAnd => {
            if !left.site || !right.site {
                return Ok(false_truth(selected));
            }
            match (left.samples, right.samples) {
                (None, None) => Ok(Truth::site(true)),
                (Some(samples), None) | (None, Some(samples)) => {
                    Ok(Truth::with_samples(true, samples.passes, samples.selected))
                }
                (Some(left), Some(right)) => {
                    let selected = selected.expect("sample truth has a selection");
                    let passes = left
                        .passes
                        .into_iter()
                        .zip(right.passes)
                        .map(|(left, right)| match operator {
                            BinaryOperator::SampleAnd => left && right,
                            BinaryOperator::SiteAnd => left || right,
                            _ => unreachable!(),
                        })
                        .collect();
                    Ok(match operator {
                        BinaryOperator::SampleAnd => Truth::selected_samples(passes, selected),
                        BinaryOperator::SiteAnd => Truth::with_samples(true, passes, selected),
                        _ => unreachable!(),
                    })
                }
            }
        }
        BinaryOperator::SampleOr | BinaryOperator::SiteOr => {
            if !left.site && !right.site {
                return Ok(false_truth(selected));
            }
            match (left.samples, right.samples) {
                (None, None) => Ok(Truth::site(true)),
                (Some(samples), None) => {
                    let passes = if operator == BinaryOperator::SiteOr && right.site {
                        samples.selected.to_vec()
                    } else {
                        samples.passes
                    };
                    Ok(Truth::with_samples(true, passes, samples.selected))
                }
                (None, Some(samples)) => {
                    let passes = if operator == BinaryOperator::SiteOr && left.site {
                        samples.selected.to_vec()
                    } else {
                        samples.passes
                    };
                    Ok(Truth::with_samples(true, passes, samples.selected))
                }
                (Some(left), Some(right)) => {
                    let selected = selected.expect("sample truth has a selection");
                    let passes = if operator == BinaryOperator::SiteOr {
                        selected.to_vec()
                    } else {
                        left.passes
                            .into_iter()
                            .zip(right.passes)
                            .map(|(left, right)| left || right)
                            .collect()
                    };
                    Ok(Truth::with_samples(true, passes, selected))
                }
            }
        }
        _ => Err(EvaluateError::new("operator is not logical")),
    }
}

fn selected_union(left: &Truth, right: &Truth) -> Result<Option<Box<[bool]>>, EvaluateError> {
    match (&left.samples, &right.samples) {
        (Some(left), Some(right)) if left.selected.len() != right.selected.len() => {
            Err(EvaluateError::new(format!(
                "incompatible sample counts in logical expression: {} vs {}",
                left.selected.len(),
                right.selected.len()
            )))
        }
        (Some(left), Some(right)) => Ok(Some(
            left.selected
                .iter()
                .zip(&right.selected)
                .map(|(left, right)| *left || *right)
                .collect(),
        )),
        (Some(samples), None) | (None, Some(samples)) => Ok(Some(samples.selected.clone())),
        (None, None) => Ok(None),
    }
}

fn false_truth(selected: Option<Box<[bool]>>) -> Truth {
    match selected {
        Some(selected) => Truth::with_samples(false, vec![false; selected.len()], selected),
        None => Truth::site(false),
    }
}
