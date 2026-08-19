use crate::expression::{
    bind::{BoundExpression, FunctionKind},
    value::{Atom, Filter, Values},
};
use noodles_vcf::variant::RecordBuf;

use super::{EvaluateError, Evaluated};

mod binomial;
mod fisher;

#[cfg(test)]
pub(crate) use binomial::binomial_two_sided;

pub(super) fn evaluate<'a>(
    kind: FunctionKind,
    mut arguments: Vec<Evaluated<'a>>,
    expressions: &[BoundExpression],
    record: &RecordBuf,
) -> Result<Evaluated<'a>, EvaluateError> {
    if matches!(kind, FunctionKind::Binomial | FunctionKind::Fisher) {
        let arguments = arguments
            .into_iter()
            .map(super::require_values)
            .collect::<Result<_, _>>()?;
        let values = match kind {
            FunctionKind::Binomial => binomial::evaluate(arguments, record),
            FunctionKind::Fisher => fisher::evaluate(arguments, expressions, record),
            _ => unreachable!(),
        }?;
        return Ok(Evaluated::Values(values));
    }
    if arguments.len() != 1 {
        return Err(EvaluateError::new(
            "function evaluation for this arity is not implemented",
        ));
    }
    let argument = arguments.pop().expect("single function argument");
    match (kind, argument) {
        (FunctionKind::Count | FunctionKind::SampleCount, Evaluated::Values(values)) => {
            count(kind, values).map(Evaluated::Values)
        }
        (FunctionKind::Count, Evaluated::Truth(truth)) => {
            let count = truth
                .samples
                .map(|samples| {
                    samples
                        .passes
                        .into_iter()
                        .zip(samples.selected)
                        .filter(|(passes, selected)| *passes && *selected)
                        .count()
                })
                .unwrap_or(usize::from(truth.site));
            Ok(Evaluated::Values(Values::Site(vec![Atom::Number(
                count as f64,
            )])))
        }
        (
            kind @ (FunctionKind::PassingSampleCount | FunctionKind::PassingSampleFraction),
            Evaluated::Truth(truth),
        ) => passing_samples(kind, truth).map(Evaluated::Values),
        (_, Evaluated::Values(values)) => apply_values(kind, values).map(Evaluated::Values),
        (_, Evaluated::Truth(_)) => Err(EvaluateError::new(
            "function requires values rather than a truth expression",
        )),
        (_, Evaluated::Set(_)) => Err(EvaluateError::new(
            "value file cannot be used as a function argument",
        )),
    }
}

fn numeric_count(value: &Atom<'_>, function: &str) -> Result<Option<i32>, EvaluateError> {
    match value {
        Atom::Absent | Atom::Missing => Ok(None),
        Atom::Number(value) if value.is_finite() && *value >= 0.0 && *value <= i32::MAX as f64 => {
            Ok(Some(*value as i32))
        }
        Atom::Flag => Ok(Some(1)),
        Atom::Number(_) => Err(EvaluateError::new(format!(
            "{function} counts must be finite nonnegative 32-bit values"
        ))),
        _ => Err(EvaluateError::new(format!(
            "{function} received a nonnumeric value"
        ))),
    }
}

fn passing_samples<'a>(
    kind: FunctionKind,
    truth: super::Truth,
) -> Result<Values<'a>, EvaluateError> {
    let Some(samples) = truth.samples else {
        return Err(EvaluateError::new(
            "passing-sample function requires a FORMAT expression",
        ));
    };
    let selected = samples
        .selected
        .iter()
        .filter(|selected| **selected)
        .count();
    let passes = samples
        .passes
        .iter()
        .zip(&samples.selected)
        .filter(|(passes, selected)| **passes && **selected)
        .count();
    let value = match kind {
        FunctionKind::PassingSampleCount => passes as f64,
        FunctionKind::PassingSampleFraction if selected == 0 => 0.0,
        FunctionKind::PassingSampleFraction => passes as f64 / selected as f64,
        _ => {
            return Err(EvaluateError::new(
                "function is not a passing-sample operation",
            ));
        }
    };
    Ok(Values::Site(vec![Atom::Number(value)]))
}

pub(super) fn apply_values<'a>(
    kind: FunctionKind,
    values: Values<'a>,
) -> Result<Values<'a>, EvaluateError> {
    match kind {
        FunctionKind::Absolute => return map_values(values, absolute),
        FunctionKind::Phred => return map_values(values, phred),
        FunctionKind::StringLength => return map_values(values, string_length),
        _ => {}
    }
    if is_sample_reduction(kind)
        && let Values::Samples(mut samples) = values
    {
        let kind = global_kind(kind);
        samples.values = samples
            .values
            .iter()
            .map(|values| reduce_atoms(kind, values.iter()).map(|value| vec![value]))
            .collect::<Result<_, _>>()?;
        return Ok(Values::Samples(samples));
    }
    let kind = global_kind(kind);
    let mut value = match values {
        Values::Site(values) => reduce_atoms(kind, values.iter())?,
        Values::Samples(samples) => reduce_atoms(
            kind,
            samples
                .values
                .iter()
                .zip(&samples.selected)
                .filter(|(_, selected)| **selected)
                .flat_map(|(values, _)| values),
        )?,
    };
    if matches!(value, Atom::Missing) {
        value = Atom::Absent;
    }
    Ok(Values::Site(vec![value]))
}

fn count<'a>(kind: FunctionKind, values: Values<'a>) -> Result<Values<'a>, EvaluateError> {
    match (kind, values) {
        (FunctionKind::Count, Values::Site(values)) => Ok(Values::Site(vec![Atom::Number(
            values
                .iter()
                .filter(|value| !matches!(value, Atom::Absent))
                .count() as f64,
        )])),
        (FunctionKind::Count, Values::Samples(samples)) => {
            let count = samples
                .values
                .iter()
                .zip(&samples.selected)
                .filter(|(_, selected)| **selected)
                .flat_map(|(values, _)| values)
                .filter(|value| !matches!(value, Atom::Absent | Atom::Missing))
                .count();
            Ok(Values::Site(vec![Atom::Number(count as f64)]))
        }
        (FunctionKind::SampleCount, Values::Samples(mut samples)) => {
            samples.values = samples
                .values
                .iter()
                .map(|values| {
                    vec![Atom::Number(
                        values
                            .iter()
                            .filter(|value| !matches!(value, Atom::Absent | Atom::Missing))
                            .count() as f64,
                    )]
                })
                .collect();
            Ok(Values::Samples(samples))
        }
        (FunctionKind::SampleCount, values @ Values::Site(_)) => count(FunctionKind::Count, values),
        _ => Err(EvaluateError::new("function is not a count operation")),
    }
}

fn map_values<'a>(
    values: Values<'a>,
    operation: fn(&Atom<'_>) -> Result<Atom<'static>, EvaluateError>,
) -> Result<Values<'a>, EvaluateError> {
    match values {
        Values::Site(values) => values
            .iter()
            .map(operation)
            .collect::<Result<_, _>>()
            .map(Values::Site),
        Values::Samples(mut samples) => {
            samples.values = samples
                .values
                .iter()
                .map(|values| values.iter().map(operation).collect::<Result<_, _>>())
                .collect::<Result<_, _>>()?;
            Ok(Values::Samples(samples))
        }
    }
}

fn absolute(atom: &Atom<'_>) -> Result<Atom<'static>, EvaluateError> {
    match atom {
        Atom::Absent => Ok(Atom::Absent),
        Atom::Missing => Ok(Atom::Missing),
        Atom::Number(value) => Ok(Atom::Number(value.abs())),
        Atom::Flag => Ok(Atom::Number(1.0)),
        _ => Err(EvaluateError::new("ABS received a nonnumeric value")),
    }
}

fn phred(atom: &Atom<'_>) -> Result<Atom<'static>, EvaluateError> {
    match atom {
        Atom::Absent => Ok(Atom::Absent),
        Atom::Missing => Ok(Atom::Missing),
        Atom::Number(value) => Ok(Atom::Number(-4.342_944_819_03 * value.ln())),
        Atom::Flag => Ok(Atom::Number(-0.0)),
        _ => Err(EvaluateError::new("PHRED received a nonnumeric value")),
    }
}

fn string_length(atom: &Atom<'_>) -> Result<Atom<'static>, EvaluateError> {
    let length = match atom {
        Atom::Absent => return Ok(Atom::Absent),
        Atom::Missing => 0,
        Atom::Text(value) => value.len(),
        Atom::OwnedText(value) => value.len(),
        Atom::Filter(Filter::Pass) => 4,
        Atom::Filter(Filter::Missing) => 0,
        Atom::Filter(Filter::Failed(filters)) => {
            filters.iter().map(|filter| filter.len()).sum::<usize>()
                + filters.len().saturating_sub(1)
        }
        _ => return Err(EvaluateError::new("STRLEN received a nonstring value")),
    };
    Ok(Atom::Number(length as f64))
}

fn reduce_atoms<'a, 'v>(
    kind: FunctionKind,
    atoms: impl Iterator<Item = &'a Atom<'v>>,
) -> Result<Atom<'static>, EvaluateError>
where
    'v: 'a,
{
    let mut values = Vec::new();
    for atom in atoms {
        match atom {
            Atom::Absent | Atom::Missing => {}
            Atom::Number(value) => values.push(*value),
            Atom::Flag => values.push(1.0),
            _ => {
                return Err(EvaluateError::new(
                    "numeric function received a nonnumeric value",
                ));
            }
        }
    }
    if values.is_empty() {
        return Ok(Atom::Missing);
    }
    let value = match kind {
        FunctionKind::Max => values.into_iter().fold(f64::NEG_INFINITY, f64::max),
        FunctionKind::Min => values.into_iter().fold(f64::INFINITY, f64::min),
        FunctionKind::Mean => values.iter().sum::<f64>() / values.len() as f64,
        FunctionKind::Median => {
            values.sort_unstable_by(f64::total_cmp);
            let middle = values.len() / 2;
            if values.len() % 2 == 0 {
                (values[middle - 1] + values[middle]) / 2.0
            } else {
                values[middle]
            }
        }
        FunctionKind::StandardDeviation => {
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            (values
                .iter()
                .map(|value| (value - mean) * (value - mean))
                .sum::<f64>()
                / values.len() as f64)
                .sqrt()
        }
        FunctionKind::Sum => values.iter().sum(),
        _ => return Err(EvaluateError::new("function is not a numeric reduction")),
    };
    Ok(Atom::Number(value))
}

fn is_sample_reduction(kind: FunctionKind) -> bool {
    matches!(
        kind,
        FunctionKind::SampleMax
            | FunctionKind::SampleMin
            | FunctionKind::SampleMean
            | FunctionKind::SampleMedian
            | FunctionKind::SampleStandardDeviation
            | FunctionKind::SampleSum
    )
}

fn global_kind(kind: FunctionKind) -> FunctionKind {
    match kind {
        FunctionKind::SampleMax => FunctionKind::Max,
        FunctionKind::SampleMin => FunctionKind::Min,
        FunctionKind::SampleMean => FunctionKind::Mean,
        FunctionKind::SampleMedian => FunctionKind::Median,
        FunctionKind::SampleStandardDeviation => FunctionKind::StandardDeviation,
        FunctionKind::SampleSum => FunctionKind::Sum,
        kind => kind,
    }
}
