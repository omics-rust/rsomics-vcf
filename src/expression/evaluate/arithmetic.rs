use crate::expression::{
    syntax::BinaryOperator,
    value::{Atom, SampleValues, Values},
};

use super::{EvaluateError, number, sample_width};

pub(super) fn negate(values: Values<'_>) -> Result<Values<'_>, EvaluateError> {
    map_unary(values, |atom| {
        if matches!(atom, Atom::Absent) {
            return Ok(Atom::Absent);
        }
        match number(atom)? {
            Some(value) => Ok(Atom::Number(-value)),
            None => Ok(Atom::Missing),
        }
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
        Values::Samples(mut samples) => {
            samples.values = samples
                .values
                .iter()
                .map(|values| values.iter().map(&operation).collect::<Result<_, _>>())
                .collect::<Result<_, _>>()?;
            Ok(Values::Samples(samples))
        }
    }
}

pub(super) fn arithmetic<'a>(
    left: Values<'a>,
    right: Values<'a>,
    operator: BinaryOperator,
) -> Result<Values<'a>, EvaluateError> {
    match (left, right) {
        (Values::Site(left), Values::Site(right)) => {
            arithmetic_vectors(&left, &right, operator).map(Values::Site)
        }
        (Values::Samples(left), Values::Samples(right)) => {
            if left.values.len() != right.values.len() {
                return Err(EvaluateError::new(format!(
                    "incompatible sample counts in arithmetic: {} vs {}",
                    left.values.len(),
                    right.values.len()
                )));
            }
            let left_width = sample_width(&left.values);
            let right_width = sample_width(&right.values);
            let width = broadcast_width(left_width, right_width)?;
            let values = left
                .values
                .iter()
                .zip(&right.values)
                .map(|(left, right)| {
                    arithmetic_sample(left, right, left_width, right_width, width, operator)
                })
                .collect::<Result<_, _>>()?;
            let selected = left
                .selected
                .iter()
                .zip(&right.selected)
                .map(|(left, right)| *left && *right)
                .collect();
            Ok(Values::Samples(SampleValues { values, selected }))
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
    mut samples: SampleValues<'a>,
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
    samples.values = samples
        .values
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
        .collect::<Result<_, _>>()?;
    Ok(Values::Samples(samples))
}

fn arithmetic_atoms<'a>(
    left: &Atom<'a>,
    right: &Atom<'a>,
    operator: BinaryOperator,
) -> Result<Atom<'a>, EvaluateError> {
    if matches!(left, Atom::Absent) || matches!(right, Atom::Absent) {
        return Ok(Atom::Absent);
    }
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
