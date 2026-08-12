use noodles_vcf::{
    Header,
    header::record::value::map::{format, info},
    variant::{
        RecordBuf,
        record::samples::series::value::genotype::Phasing,
        record_buf::{
            Samples,
            info::field::{Value as InfoValue, value::Array as InfoArray},
            samples::sample::{
                Value as SampleValue,
                value::{Array as SampleArray, genotype::allele::Allele},
            },
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use super::super::cardinality::{combinations, genotype_index, infer_ploidy, visit_genotypes};

pub(super) fn merge_info(
    header: &Header,
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    output: &mut RecordBuf,
) -> Result<()> {
    for (key, schema) in header.infos() {
        if records[0].info().get(key).is_none() {
            continue;
        }
        let sources = records
            .iter()
            .map(|record| record.info().get(key).flatten().cloned())
            .collect::<Vec<_>>();
        let value = merge_info_value(key, schema.number(), &sources, records, mappings)?;
        output.info_mut().insert(key.clone(), value);
    }
    Ok(())
}

fn merge_info_value(
    key: &str,
    number: info::Number,
    sources: &[Option<InfoValue>],
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
) -> Result<Option<InfoValue>> {
    let Some(first) = sources[0].as_ref() else {
        return Ok(None);
    };
    if matches!(number, info::Number::Count(_) | info::Number::Unknown) {
        return Ok(Some(first.clone()));
    }
    let target_alleles = mappings
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(1, |maximum| maximum + 1);
    match first {
        InfoValue::Array(InfoArray::Integer(_)) => merge_info_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::info(number, target_alleles, key, true),
            |value| match value {
                InfoValue::Array(InfoArray::Integer(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(InfoArray::Integer)
        .map(InfoValue::Array)
        .map(Some),
        InfoValue::Array(InfoArray::Float(_)) => merge_info_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::info(number, target_alleles, key, true),
            |value| match value {
                InfoValue::Array(InfoArray::Float(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(InfoArray::Float)
        .map(InfoValue::Array)
        .map(Some),
        InfoValue::Array(InfoArray::Character(_)) => merge_info_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::info(number, target_alleles, key, false),
            |value| match value {
                InfoValue::Array(InfoArray::Character(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(InfoArray::Character)
        .map(InfoValue::Array)
        .map(Some),
        InfoValue::Array(InfoArray::String(_)) => merge_info_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::info(number, target_alleles, key, false),
            |value| match value {
                InfoValue::Array(InfoArray::String(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(InfoArray::String)
        .map(InfoValue::Array)
        .map(Some),
        _ => Ok(Some(first.clone())),
    }
}

fn merge_info_arrays<T: Clone>(
    sources: &[Option<InfoValue>],
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    merge: ArrayMerge,
    extract: impl Fn(&InfoValue) -> Option<Vec<Option<T>>>,
) -> Result<Vec<Option<T>>> {
    merge_arrays(
        sources
            .iter()
            .map(|source| source.as_ref().and_then(&extract)),
        records,
        mappings,
        merge,
    )
}

pub(super) fn merge_samples(
    header: &Header,
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    output: &mut RecordBuf,
) -> Result<()> {
    let keys = records[0].samples().keys().clone();
    let mut samples = Vec::new();
    for sample_index in 0..records[0].samples().values().count() {
        let mut values = records[0]
            .samples()
            .values()
            .nth(sample_index)
            .unwrap()
            .values()
            .to_vec();
        for (output_index, key) in keys.as_ref().iter().enumerate() {
            let sources = records
                .iter()
                .map(|record| {
                    record
                        .samples()
                        .keys()
                        .as_ref()
                        .get_index_of(key)
                        .and_then(|index| {
                            record
                                .samples()
                                .values()
                                .nth(sample_index)?
                                .values()
                                .get(index)
                        })
                        .cloned()
                        .flatten()
                })
                .collect::<Vec<_>>();
            values[output_index] = if key == "GT" {
                merge_genotypes(&sources, mappings, sample_index)?
            } else {
                let schema = header
                    .formats()
                    .get(key)
                    .ok_or_else(|| invalid(format!("FORMAT/{key} is absent from the header")))?;
                merge_sample_value(
                    key,
                    schema.number(),
                    &sources,
                    records,
                    mappings,
                    sample_index,
                )?
            };
        }
        samples.push(values);
    }
    *output.samples_mut() = Samples::new(keys, samples);
    Ok(())
}

fn merge_sample_value(
    key: &str,
    number: format::Number,
    sources: &[Option<SampleValue>],
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    sample: usize,
) -> Result<Option<SampleValue>> {
    let Some(first) = sources[0].as_ref() else {
        return Ok(None);
    };
    if !matches!(
        number,
        format::Number::AlternateBases
            | format::Number::ReferenceAlternateBases
            | format::Number::Samples
    ) {
        return Ok(Some(first.clone()));
    }
    let target_alleles = mappings
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(1, |maximum| maximum + 1);
    let context = format!("FORMAT/{key} sample {}", sample + 1);
    match first {
        SampleValue::Array(SampleArray::Integer(_)) => merge_sample_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::sample(number, target_alleles, &context, true),
            |value| match value {
                SampleValue::Array(SampleArray::Integer(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(SampleArray::Integer)
        .map(SampleValue::Array)
        .map(Some),
        SampleValue::Array(SampleArray::Float(_)) => merge_sample_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::sample(number, target_alleles, &context, true),
            |value| match value {
                SampleValue::Array(SampleArray::Float(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(SampleArray::Float)
        .map(SampleValue::Array)
        .map(Some),
        SampleValue::Array(SampleArray::Character(_)) => merge_sample_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::sample(number, target_alleles, &context, false),
            |value| match value {
                SampleValue::Array(SampleArray::Character(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(SampleArray::Character)
        .map(SampleValue::Array)
        .map(Some),
        SampleValue::Array(SampleArray::String(_)) => merge_sample_arrays(
            sources,
            records,
            mappings,
            ArrayMerge::sample(number, target_alleles, &context, false),
            |value| match value {
                SampleValue::Array(SampleArray::String(values)) => Some(values.clone()),
                _ => None,
            },
        )
        .map(SampleArray::String)
        .map(SampleValue::Array)
        .map(Some),
        _ => Err(invalid(format!("{context} is not encoded as an array"))),
    }
}

fn merge_sample_arrays<T: Clone>(
    sources: &[Option<SampleValue>],
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    merge: ArrayMerge,
    extract: impl Fn(&SampleValue) -> Option<Vec<Option<T>>>,
) -> Result<Vec<Option<T>>> {
    merge_arrays(
        sources
            .iter()
            .map(|source| source.as_ref().and_then(&extract)),
        records,
        mappings,
        merge,
    )
}

#[derive(Clone, Copy)]
enum Number {
    Alternate,
    ReferenceAlternate,
    Genotype,
}

#[derive(Clone, Copy)]
enum GenotypePolicy {
    Info,
    NumericFormat,
    StringFormat,
}

struct ArrayMerge {
    number: Number,
    target_alleles: usize,
    context: String,
    overwrite: bool,
    genotype_policy: GenotypePolicy,
}

impl ArrayMerge {
    fn info(number: info::Number, target_alleles: usize, key: &str, overwrite: bool) -> Self {
        Self {
            number: number.into(),
            target_alleles,
            context: format!("INFO/{key}"),
            overwrite,
            genotype_policy: GenotypePolicy::Info,
        }
    }

    fn sample(
        number: format::Number,
        target_alleles: usize,
        context: &str,
        overwrite: bool,
    ) -> Self {
        Self {
            number: number.into(),
            target_alleles,
            context: context.to_owned(),
            overwrite,
            genotype_policy: if overwrite {
                GenotypePolicy::NumericFormat
            } else {
                GenotypePolicy::StringFormat
            },
        }
    }
}

impl From<info::Number> for Number {
    fn from(number: info::Number) -> Self {
        match number {
            info::Number::AlternateBases => Self::Alternate,
            info::Number::ReferenceAlternateBases => Self::ReferenceAlternate,
            info::Number::Samples => Self::Genotype,
            _ => unreachable!(),
        }
    }
}

impl From<format::Number> for Number {
    fn from(number: format::Number) -> Self {
        match number {
            format::Number::AlternateBases => Self::Alternate,
            format::Number::ReferenceAlternateBases => Self::ReferenceAlternate,
            format::Number::Samples => Self::Genotype,
            _ => unreachable!(),
        }
    }
}

fn merge_arrays<T: Clone>(
    sources: impl IntoIterator<Item = Option<Vec<Option<T>>>>,
    records: &[&RecordBuf],
    mappings: &[Vec<usize>],
    merge: ArrayMerge,
) -> Result<Vec<Option<T>>> {
    let sources: Vec<_> = sources.into_iter().collect();
    let (source_ploidies, target_ploidy) = genotype_ploidies(
        &sources,
        records,
        merge.number,
        &merge.context,
        merge.genotype_policy,
    )?;
    let target = match merge.number {
        Number::Alternate => merge.target_alleles - 1,
        Number::ReferenceAlternate => merge.target_alleles,
        Number::Genotype => {
            combinations(merge.target_alleles + target_ploidy - 1, target_ploidy)
                .ok_or_else(|| invalid(format!("{} cardinality overflows", merge.context)))?
        }
    };
    let mut output = vec![None; target];
    for (record_index, source) in sources.into_iter().enumerate() {
        let Some(source) = source else {
            continue;
        };
        let source_alleles = records[record_index].alternate_bases().as_ref().len() + 1;
        let source_ploidy = source_ploidies[record_index];
        let expected = match merge.number {
            Number::Alternate => source_alleles - 1,
            Number::ReferenceAlternate => source_alleles,
            Number::Genotype => combinations(source_alleles + source_ploidy - 1, source_ploidy)
                .ok_or_else(|| invalid(format!("{} cardinality overflows", merge.context)))?,
        };
        if source.len() != expected
            && !matches!(
                (merge.number, merge.genotype_policy),
                (Number::Genotype, GenotypePolicy::StringFormat)
            )
        {
            return Err(invalid(format!(
                "{} has {} values, expected {expected}",
                merge.context,
                source.len()
            )));
        }
        for (source_index, target_index) in
            field_mapping(merge.number, &mappings[record_index], source_ploidy)?
        {
            if merge.overwrite || output[target_index].is_none() {
                output[target_index] = source.get(source_index).cloned().flatten();
            }
        }
    }
    Ok(output)
}

fn genotype_ploidies<T>(
    sources: &[Option<Vec<Option<T>>>],
    records: &[&RecordBuf],
    number: Number,
    context: &str,
    policy: GenotypePolicy,
) -> Result<(Vec<usize>, usize)> {
    if !matches!(number, Number::Genotype) {
        return Ok((vec![0; sources.len()], 0));
    }
    if matches!(policy, GenotypePolicy::Info) {
        return Ok((vec![2; sources.len()], 2));
    }
    let inferred = sources
        .iter()
        .zip(records)
        .map(|(source, record)| {
            source.as_ref().and_then(|values| {
                infer_ploidy(record.alternate_bases().as_ref().len() + 1, values.len())
                    .filter(|ploidy| matches!(ploidy, 1 | 2))
            })
        })
        .collect::<Vec<_>>();
    let first = inferred[0].ok_or_else(|| {
        invalid(format!(
            "{context} cardinality is neither haploid nor diploid"
        ))
    })?;
    let target = match policy {
        GenotypePolicy::NumericFormat => inferred.iter().flatten().copied().max().unwrap_or(first),
        GenotypePolicy::StringFormat => first,
        GenotypePolicy::Info => unreachable!(),
    };
    let source_ploidies = match policy {
        GenotypePolicy::NumericFormat => inferred
            .into_iter()
            .map(|ploidy| ploidy.unwrap_or(target))
            .collect(),
        GenotypePolicy::StringFormat => vec![target; sources.len()],
        GenotypePolicy::Info => unreachable!(),
    };
    Ok((source_ploidies, target))
}

fn field_mapping(number: Number, alleles: &[usize], ploidy: usize) -> Result<Vec<(usize, usize)>> {
    match number {
        Number::Alternate => Ok((1..alleles.len())
            .map(|source| (source - 1, alleles[source] - 1))
            .collect()),
        Number::ReferenceAlternate => Ok(alleles.iter().copied().enumerate().collect()),
        Number::Genotype => {
            let mut mapping = Vec::new();
            visit_genotypes(alleles.len(), ploidy, |genotype| {
                let source = genotype_index(genotype).unwrap();
                let mut mapped: Vec<_> = genotype.iter().map(|&allele| alleles[allele]).collect();
                mapped.sort_unstable();
                mapping.push((source, genotype_index(&mapped).unwrap()));
            });
            Ok(mapping)
        }
    }
}

fn merge_genotypes(
    sources: &[Option<SampleValue>],
    mappings: &[Vec<usize>],
    sample: usize,
) -> Result<Option<SampleValue>> {
    let Some(SampleValue::Genotype(mut output)) = sources[0].clone() else {
        return Ok(sources[0].clone());
    };
    remap_genotype(&mut output, &mappings[0], sample)?;
    for (source, mapping) in sources.iter().skip(1).zip(&mappings[1..]) {
        let Some(SampleValue::Genotype(source)) = source else {
            continue;
        };
        while output.as_ref().len() < source.as_ref().len() {
            output.as_mut().push(Allele::new(None, Phasing::Unphased));
        }
        for (slot, allele) in source.as_ref().iter().enumerate() {
            let Some(position) = allele.position() else {
                continue;
            };
            let mapped = *mapping.get(position).ok_or_else(|| {
                invalid(format!(
                    "FORMAT/GT sample {} allele {position} exceeds {} ALT alleles",
                    sample + 1,
                    mapping.len() - 1
                ))
            })?;
            let destination = &mut output.as_mut()[slot];
            if position == 0 && destination.position().is_some() {
                continue;
            }
            if destination.position().is_none_or(|position| position == 0) {
                *destination.position_mut() = Some(mapped);
                continue;
            }
            if let Some(destination) = output
                .as_mut()
                .iter_mut()
                .find(|allele| allele.position().is_none_or(|position| position == 0))
            {
                *destination.position_mut() = Some(mapped);
                *destination.phasing_mut() = Phasing::Unphased;
            }
        }
    }
    Ok(Some(SampleValue::Genotype(output)))
}

fn remap_genotype(
    genotype: &mut noodles_vcf::variant::record_buf::samples::sample::value::Genotype,
    mapping: &[usize],
    sample: usize,
) -> Result<()> {
    for allele in genotype.as_mut() {
        let Some(position) = allele.position() else {
            continue;
        };
        *allele.position_mut() = Some(*mapping.get(position).ok_or_else(|| {
            invalid(format!(
                "FORMAT/GT sample {} allele {position} exceeds {} ALT alleles",
                sample + 1,
                mapping.len() - 1
            ))
        })?);
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
