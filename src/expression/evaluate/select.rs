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
            select_values(values, selector, None).map(Values::Site)
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
        *values = select_values(std::mem::take(values), value_selector, genotype)?;
    }
    Ok(samples)
}

fn select_values<'a>(
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
