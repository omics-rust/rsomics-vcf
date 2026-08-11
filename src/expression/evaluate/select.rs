use std::collections::BTreeSet;

use noodles_vcf::variant::RecordBuf;

use crate::expression::{
    bind::{BoundSubscript, SampleSelector, ValueSelector},
    syntax::IndexSelector,
    value::{self, Atom, SampleValues, Values},
};

use super::EvaluateError;

pub(super) fn apply<'a>(
    subscript: &BoundSubscript,
    values: Values<'a>,
    record: &RecordBuf,
) -> Result<Values<'a>, EvaluateError> {
    match (subscript, values) {
        (BoundSubscript::Values(selector), Values::Site(values)) => {
            select_site_values(values, selector).map(Values::Site)
        }
        (BoundSubscript::SampleValues { samples, values }, Values::Samples(sample_values)) => {
            select_samples(sample_values, samples, values, record).map(Values::Samples)
        }
        (BoundSubscript::Values(_), Values::Samples(_)) => Err(EvaluateError::new(
            "site subscript used with a FORMAT field",
        )),
        (BoundSubscript::SampleValues { .. }, Values::Site(_)) => Err(EvaluateError::new(
            "sample subscript used with a site field",
        )),
    }
}

fn select_samples<'a>(
    mut samples: SampleValues<'a>,
    sample_selector: &SampleSelector,
    value_selector: &ValueSelector,
    record: &RecordBuf,
) -> Result<SampleValues<'a>, EvaluateError> {
    if let SampleSelector::Selected(selected) = sample_selector {
        if selected.len() != samples.values.len() {
            return Err(EvaluateError::new(format!(
                "subscript has {} samples but record has {}",
                selected.len(),
                samples.values.len()
            )));
        }
        for (current, selected) in samples.selected.iter_mut().zip(selected) {
            *current &= selected;
        }
    }
    let genotypes = if matches!(value_selector, ValueSelector::Genotype) {
        Some(
            value::genotype_indices(record)
                .map_err(|error| EvaluateError::new(error.to_string()))?,
        )
    } else {
        None
    };
    for (index, values) in samples.values.iter_mut().enumerate() {
        let genotype = genotypes
            .as_ref()
            .map(|genotypes| genotypes[index].as_slice());
        *values = select_sample_values(std::mem::take(values), value_selector, genotype)?;
    }
    Ok(samples)
}

fn select_site_values<'a>(
    values: Vec<Atom<'a>>,
    selector: &ValueSelector,
) -> Result<Vec<Atom<'a>>, EvaluateError> {
    match selector {
        ValueSelector::All => Ok(values),
        ValueSelector::Genotype => Err(EvaluateError::new(
            "GT selection requires per-sample genotypes",
        )),
        ValueSelector::Indices(ranges)
            if matches!(ranges.as_slice(), [IndexSelector { end: None, .. }]) =>
        {
            if values.len() == 1 && !matches!(values[0], Atom::Absent) {
                return Ok(values);
            }
            let index = ranges[0].start;
            Ok(match values.get(index) {
                Some(Atom::Absent | Atom::Missing) | None => vec![Atom::Absent],
                Some(value) => vec![value.clone()],
            })
        }
        ValueSelector::Indices(ranges) => {
            let selected: Vec<_> = expand_indices(ranges, values.len())
                .into_iter()
                .filter_map(|index| values.get(index).cloned())
                .collect();
            Ok(if selected.is_empty() {
                vec![Atom::Absent]
            } else {
                selected
            })
        }
    }
}

fn select_sample_values<'a>(
    values: Vec<Atom<'a>>,
    selector: &ValueSelector,
    genotype: Option<&[usize]>,
) -> Result<Vec<Atom<'a>>, EvaluateError> {
    let indices = match selector {
        ValueSelector::All => return Ok(values),
        ValueSelector::Genotype => genotype
            .ok_or_else(|| EvaluateError::new("GT selection requires per-sample genotypes"))?
            .to_vec(),
        ValueSelector::Indices(ranges) => expand_indices(ranges, values.len()),
    };
    let selected: Vec<_> = indices
        .into_iter()
        .map(|index| values.get(index).cloned().unwrap_or(Atom::Missing))
        .collect();
    Ok(if selected.is_empty() {
        vec![Atom::Missing]
    } else {
        selected
    })
}

fn expand_indices(ranges: &[IndexSelector], value_count: usize) -> Vec<usize> {
    let mut indices = BTreeSet::new();
    for range in ranges {
        let end = match range.end {
            Some(usize::MAX) => value_count.saturating_sub(1),
            Some(end) => end,
            None => range.start,
        };
        if range.start <= end {
            indices.extend(range.start..=end);
        }
    }
    indices.into_iter().collect()
}
