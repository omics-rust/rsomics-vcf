use crate::expression::{
    bind::FunctionKind,
    value::{Atom, Values},
};

use super::EvaluateError;

pub(super) fn reduce<'a>(
    kind: FunctionKind,
    values: Values<'a>,
) -> Result<Values<'a>, EvaluateError> {
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
    let value = match values {
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
    Ok(Values::Site(vec![value]))
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
            Atom::Missing => {}
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
