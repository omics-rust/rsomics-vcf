use std::collections::HashSet;
use std::path::Path;

use noodles_vcf::{
    Header,
    header::record::value::map::format::{self, Number, Type},
    variant::{
        RecordBuf,
        record_buf::{
            Samples,
            samples::sample::{Value, value::Array},
        },
    },
};
use rsomics_common::Result;

use super::{Editor, State, choose, invalid, remap_array};
use crate::{
    annotate::{
        columns::{Destination, SourceField, WriteMode},
        matching::Matched,
        source::Payload,
    },
    norm::cardinality::infer_ploidy,
};

#[derive(Clone)]
pub(super) struct FormatPlan {
    pub(super) number: Number,
    pub(super) ty: Type,
    pub(super) genotype: bool,
    pub(super) context: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SampleSelection {
    source_to_target: Vec<(usize, usize)>,
}

impl SampleSelection {
    pub(crate) fn bind(
        source: &Header,
        target: &Header,
        list: Option<&str>,
        file: Option<&Path>,
    ) -> Result<Self> {
        if list.is_some() && file.is_some() {
            return Err(invalid(
                "annotation samples and samples file are mutually exclusive",
            ));
        }
        let (exclude, requested) = match (list, file) {
            (Some(list), None) => parse_sample_list(list)?,
            (None, Some(path)) => {
                let raw = path.to_string_lossy();
                let (exclude, path) = match raw.strip_prefix('^') {
                    Some(path) => (true, Path::new(path)),
                    None => (false, path),
                };
                let content = std::fs::read_to_string(path).map_err(|error| {
                    invalid(format!(
                        "reading annotation samples {}: {error}",
                        path.display()
                    ))
                })?;
                let names = content
                    .lines()
                    .map(str::trim)
                    .filter(|line| !line.is_empty() && !line.starts_with('#'))
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                if names.is_empty() {
                    return Err(invalid("annotation samples file is empty"));
                }
                (exclude, names)
            }
            (None, None) => (false, Vec::new()),
            (Some(_), Some(_)) => unreachable!(),
        };
        let explicit = list.is_some() || file.is_some();
        let requested = requested.into_iter().collect::<HashSet<_>>();
        let mut pairs = Vec::new();
        for (source_index, name) in source.sample_names().iter().enumerate() {
            let selected = if !explicit {
                true
            } else if exclude {
                !requested.contains(name)
            } else {
                requested.contains(name)
            };
            if !selected {
                continue;
            }
            if let Some(target_index) = target.sample_names().get_index_of(name) {
                pairs.push((source_index, target_index));
            } else if explicit {
                return Err(invalid(format!(
                    "annotation sample {name:?} is not present in the target"
                )));
            }
        }
        if explicit {
            for name in &requested {
                if !source.sample_names().contains(name) {
                    return Err(invalid(format!(
                        "annotation sample {name:?} is not present in the source"
                    )));
                }
            }
        }
        if pairs.is_empty() {
            return Err(invalid("annotation sample selection has no target samples"));
        }
        Ok(Self {
            source_to_target: pairs,
        })
    }

    pub(crate) fn pairs(&self) -> &[(usize, usize)] {
        &self.source_to_target
    }
}

fn parse_sample_list(raw: &str) -> Result<(bool, Vec<String>)> {
    let (exclude, raw) = raw
        .strip_prefix('^')
        .map_or((false, raw), |raw| (true, raw));
    let names = raw
        .split(',')
        .map(str::trim)
        .map(|name| {
            if name.is_empty() {
                Err(invalid("annotation samples contain an empty name"))
            } else {
                Ok(name.to_owned())
            }
        })
        .collect::<Result<Vec<_>>>()?;
    if names.is_empty() {
        return Err(invalid("annotation samples must not be empty"));
    }
    Ok((exclude, names))
}

impl Editor {
    pub(crate) fn set_samples(&mut self, selection: SampleSelection) {
        self.samples = Some(selection);
    }

    pub(crate) fn apply_samples(
        &self,
        header: &Header,
        matched: &Matched<'_>,
        target: &mut RecordBuf,
    ) -> Result<bool> {
        let format_transfers = self
            .transfers
            .iter()
            .filter(|bound| bound.format.is_some())
            .collect::<Vec<_>>();
        if format_transfers.is_empty() {
            return Ok(false);
        }
        let selection = self
            .samples
            .as_ref()
            .ok_or_else(|| invalid("FORMAT transfer has no sample selection"))?;
        let Payload::Variant(source) = &matched.source.payload else {
            return Err(invalid(
                "FORMAT transfer requires a variant annotation source",
            ));
        };
        let original = target.samples().clone();
        let (mut keys, mut values) = original.clone().into();
        if values.len() != header.sample_names().len() {
            return Err(invalid("target sample count does not match its header"));
        }
        let key_count = keys.as_ref().len();
        for row in &mut values {
            row.resize(key_count, None);
        }

        for bound in format_transfers {
            let plan = bound.format.as_ref().expect("FORMAT transfer has a schema");
            let SourceField::Format(source_key) = &bound.transfer.source else {
                unreachable!("FORMAT transfer has a FORMAT source");
            };
            let Destination::Format(target_key) = &bound.transfer.destination else {
                unreachable!("FORMAT transfer has a FORMAT destination");
            };
            let Some(source_key_index) = source.samples().keys().as_ref().get_index_of(source_key)
            else {
                continue;
            };
            if bound.transfer.mode == WriteMode::ReplaceExisting
                && !keys.as_ref().contains(target_key)
            {
                continue;
            }
            let (target_key_index, existed) = match keys.as_ref().get_index_of(target_key) {
                Some(index) => (index, true),
                None => {
                    let index = keys.as_ref().len();
                    keys.as_mut().insert(target_key.clone());
                    for row in &mut values {
                        row.resize(index + 1, None);
                    }
                    (index, false)
                }
            };

            for &(source_sample, target_sample) in selection.pairs() {
                let source_row = source.samples().get_index(source_sample).ok_or_else(|| {
                    invalid(format!(
                        "annotation record {} is missing source sample {}",
                        matched.source.serial,
                        source_sample + 1
                    ))
                })?;
                let source_value = source_row
                    .values()
                    .get(source_key_index)
                    .cloned()
                    .flatten()
                    .map_or(State::Missing, State::Value);
                let source_value = match source_value {
                    State::Value(value) => State::Value(remap_value(
                        plan,
                        value,
                        matched,
                        target.alternate_bases().as_ref().len() + 1,
                    )?),
                    state => state,
                };
                let row = values.get_mut(target_sample).ok_or_else(|| {
                    invalid(format!("target is missing sample {}", target_sample + 1))
                })?;
                row.resize(keys.as_ref().len(), None);
                let current = if existed {
                    row.get(target_key_index)
                        .cloned()
                        .flatten()
                        .map_or(State::Missing, State::Value)
                } else {
                    State::Absent
                };
                let next = if matches!(
                    bound.transfer.mode,
                    WriteMode::Replace | WriteMode::ReplaceExisting
                ) && matches!(source_value, State::Missing)
                {
                    State::Missing
                } else {
                    choose(bound.transfer.mode, source_value, current.clone())?
                };
                row[target_key_index] = match next {
                    State::Value(value) => Some(value),
                    State::Absent | State::Missing | State::Remove => None,
                };
            }
        }

        let updated = Samples::new(keys, values);
        if updated == original {
            Ok(false)
        } else {
            *target.samples_mut() = updated;
            Ok(true)
        }
    }
}

pub(super) fn prepare_definition(
    target: &mut Header,
    key: &str,
    source: &noodles_vcf::header::record::value::Map<format::Format>,
) -> Result<()> {
    match target.formats().get(key) {
        Some(existing) if existing.number() != source.number() || existing.ty() != source.ty() => {
            Err(invalid(format!("FORMAT/{key} has an incompatible schema")))
        }
        Some(_) => Ok(()),
        None => {
            target.formats_mut().insert(key.to_owned(), source.clone());
            Ok(())
        }
    }
}

fn remap_value(
    plan: &FormatPlan,
    value: Value,
    matched: &Matched<'_>,
    target_alleles: usize,
) -> Result<Value> {
    if let Value::Genotype(mut genotype) = value {
        if !plan.genotype {
            return Err(invalid(format!(
                "{} has an unexpected genotype value",
                plan.context
            )));
        }
        for allele in genotype.as_mut() {
            let Some(source) = allele.position() else {
                continue;
            };
            let target = matched
                .allele_map
                .get(source)
                .copied()
                .flatten()
                .ok_or_else(|| {
                    invalid(format!(
                        "{} allele {source} cannot be represented in the target",
                        plan.context
                    ))
                })?;
            *allele.position_mut() = Some(target);
        }
        return Ok(Value::Genotype(genotype));
    }
    validate_value(&value, plan, matched.allele_map.len())?;
    let number = match plan.number {
        Number::AlternateBases => {
            noodles_vcf::header::record::value::map::info::Number::AlternateBases
        }
        Number::ReferenceAlternateBases => {
            noodles_vcf::header::record::value::map::info::Number::ReferenceAlternateBases
        }
        Number::Samples => noodles_vcf::header::record::value::map::info::Number::Samples,
        _ => return Ok(value),
    };
    macro_rules! remap {
        ($variant:ident, $values:expr) => {
            remap_array(number, $values, matched, target_alleles, &plan.context)
                .map(Array::$variant)
                .map(Value::Array)
        };
    }
    match value {
        Value::Array(Array::Integer(values)) => remap!(Integer, values),
        Value::Array(Array::Float(values)) => remap!(Float, values),
        Value::Array(Array::Character(values)) => remap!(Character, values),
        Value::Array(Array::String(values)) => remap!(String, values),
        _ => Err(invalid(format!(
            "{} allele-indexed value must be an array",
            plan.context
        ))),
    }
}

fn validate_value(value: &Value, plan: &FormatPlan, alleles: usize) -> Result<()> {
    let (ty, count) = match value {
        Value::Integer(_) => (Type::Integer, 1),
        Value::Float(_) => (Type::Float, 1),
        Value::Character(_) => (Type::Character, 1),
        Value::String(_) => (Type::String, 1),
        Value::Array(Array::Integer(values)) => (Type::Integer, values.len()),
        Value::Array(Array::Float(values)) => (Type::Float, values.len()),
        Value::Array(Array::Character(values)) => (Type::Character, values.len()),
        Value::Array(Array::String(values)) => (Type::String, values.len()),
        Value::Genotype(_) => {
            return Err(invalid(format!(
                "{} has an unexpected genotype value",
                plan.context
            )));
        }
    };
    if ty != plan.ty {
        return Err(invalid(format!(
            "{} has type {ty}, expected {}",
            plan.context, plan.ty
        )));
    }
    let valid = match plan.number {
        Number::Count(expected) => count == expected,
        Number::AlternateBases => count == alleles - 1,
        Number::ReferenceAlternateBases => count == alleles,
        Number::Samples => infer_ploidy(alleles, count).is_some(),
        Number::Unknown => true,
        Number::LocalAlternateBases
        | Number::LocalReferenceAlternateBases
        | Number::LocalSamples
        | Number::Ploidy
        | Number::BaseModifications => {
            return Err(invalid(format!(
                "{} uses an unsupported cardinality",
                plan.context
            )));
        }
    };
    if !valid {
        return Err(invalid(format!(
            "{} has invalid cardinality {count}",
            plan.context
        )));
    }
    Ok(())
}
