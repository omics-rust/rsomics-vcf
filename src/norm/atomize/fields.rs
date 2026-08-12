use noodles_vcf::{
    Header,
    header::record::value::map::{format, info},
    variant::{
        RecordBuf,
        record_buf::{
            Samples,
            info::field::{Value as InfoValue, value::Array as InfoArray},
            samples::sample::{Value as SampleValue, value::Array as SampleArray},
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use crate::norm::cardinality::{combinations, infer_ploidy};

pub(super) fn remap_genotypes(
    source: &RecordBuf,
    output: &mut RecordBuf,
    mapping: &[usize],
    star_allele: bool,
) -> Result<()> {
    let keys = output.samples().keys().clone();
    let Some(genotype_index) = keys.as_ref().get_index_of("GT") else {
        return Ok(());
    };
    let mut samples = Vec::new();
    for (sample_index, (source_sample, output_sample)) in source
        .samples()
        .values()
        .zip(output.samples().values())
        .enumerate()
    {
        let mut values = output_sample.values().to_vec();
        let genotype = source_sample
            .values()
            .get(genotype_index)
            .cloned()
            .ok_or_else(|| invalid(format!("FORMAT/GT sample {} is missing", sample_index + 1)))?;
        values[genotype_index] = remap_genotype(genotype, mapping, sample_index, star_allele)?;
        samples.push(values);
    }
    *output.samples_mut() = Samples::new(keys, samples);
    Ok(())
}

pub(super) fn extend_allele_fields(
    header: &Header,
    record: &mut RecordBuf,
    star: bool,
) -> Result<()> {
    let mut remove = Vec::new();
    for (key, schema) in header.infos() {
        let Some(value) = record.info_mut().get_mut(key) else {
            continue;
        };
        match schema.number() {
            info::Number::AlternateBases => {
                if star {
                    append_info_missing(value, 1, 1, &format!("INFO/{key}"))?;
                }
            }
            info::Number::ReferenceAlternateBases => {
                if star {
                    append_info_missing(value, 2, 1, &format!("INFO/{key}"))?;
                }
            }
            info::Number::Samples => remove.push(key.clone()),
            info::Number::Count(_) | info::Number::Unknown => {}
        }
    }
    for key in remove {
        record.info_mut().as_mut().shift_remove(&key);
    }

    let keys = record.samples().keys().clone();
    let mut samples = Vec::new();
    for (sample_index, sample) in record.samples().values().enumerate() {
        let mut values = sample.values().to_vec();
        for (index, key) in keys.as_ref().iter().enumerate() {
            let Some(value) = values.get_mut(index) else {
                continue;
            };
            let schema = header
                .formats()
                .get(key)
                .ok_or_else(|| invalid(format!("FORMAT/{key} is absent from the header")))?;
            let context = format!("FORMAT/{key} sample {}", sample_index + 1);
            match schema.number() {
                format::Number::AlternateBases => {
                    if star {
                        append_sample_missing(value, 1, 1, &context)?;
                    }
                }
                format::Number::ReferenceAlternateBases => {
                    if star {
                        append_sample_missing(value, 2, 1, &context)?;
                    }
                }
                format::Number::Samples if star => extend_sample_genotypes(value, &context)?,
                format::Number::Samples => {}
                format::Number::Count(_)
                | format::Number::LocalAlternateBases
                | format::Number::LocalReferenceAlternateBases
                | format::Number::LocalSamples
                | format::Number::Ploidy
                | format::Number::BaseModifications
                | format::Number::Unknown => {}
            }
        }
        samples.push(values);
    }
    *record.samples_mut() = Samples::new(keys, samples);
    Ok(())
}

fn remap_genotype(
    value: Option<SampleValue>,
    mapping: &[usize],
    sample: usize,
    star_allele: bool,
) -> Result<Option<SampleValue>> {
    let Some(SampleValue::Genotype(mut genotype)) = value else {
        return Ok(value);
    };
    for allele in genotype.as_mut() {
        let Some(position) = allele.position() else {
            continue;
        };
        if position == 0 {
            continue;
        }
        let mapped = *mapping.get(position - 1).ok_or_else(|| {
            invalid(format!(
                "FORMAT/GT sample {} allele {position} exceeds {} ALT alleles",
                sample + 1,
                mapping.len()
            ))
        })?;
        *allele.position_mut() = (mapped != 2 || star_allele).then_some(mapped);
    }
    Ok(Some(SampleValue::Genotype(genotype)))
}

fn append_info_missing(
    value: &mut Option<InfoValue>,
    expected: usize,
    count: usize,
    context: &str,
) -> Result<()> {
    let Some(InfoValue::Array(array)) = value else {
        return Err(invalid(format!("{context} is not encoded as an array")));
    };
    match array {
        InfoArray::Integer(values) => append_missing(values, expected, count, context),
        InfoArray::Float(values) => append_missing(values, expected, count, context),
        InfoArray::Character(values) => append_missing(values, expected, count, context),
        InfoArray::String(values) => append_missing(values, expected, count, context),
    }
}

fn append_sample_missing(
    value: &mut Option<SampleValue>,
    expected: usize,
    count: usize,
    context: &str,
) -> Result<()> {
    let Some(SampleValue::Array(array)) = value else {
        return Err(invalid(format!("{context} is not encoded as an array")));
    };
    match array {
        SampleArray::Integer(values) => append_missing(values, expected, count, context),
        SampleArray::Float(values) => append_missing(values, expected, count, context),
        SampleArray::Character(values) => append_missing(values, expected, count, context),
        SampleArray::String(values) => append_missing(values, expected, count, context),
    }
}

fn extend_sample_genotypes(value: &mut Option<SampleValue>, context: &str) -> Result<()> {
    let Some(SampleValue::Array(array)) = value else {
        return Err(invalid(format!("{context} is not encoded as an array")));
    };
    match array {
        SampleArray::Integer(values) => extend_genotypes(values, context),
        SampleArray::Float(values) => extend_genotypes(values, context),
        SampleArray::Character(values) => extend_genotypes(values, context),
        SampleArray::String(values) => extend_genotypes(values, context),
    }
}

fn append_missing<T>(
    values: &mut Vec<Option<T>>,
    expected: usize,
    count: usize,
    context: &str,
) -> Result<()> {
    if values.len() != expected {
        return Err(invalid(format!(
            "{context} has {} values, expected {expected}",
            values.len()
        )));
    }
    values.extend((0..count).map(|_| None));
    Ok(())
}

fn extend_genotypes<T>(values: &mut Vec<Option<T>>, context: &str) -> Result<()> {
    if values.len() == 1 && values[0].is_none() {
        return Ok(());
    }
    let ploidy = infer_ploidy(2, values.len()).ok_or_else(|| {
        invalid(format!(
            "{context} cardinality {} does not identify a ploidy for 2 alleles",
            values.len()
        ))
    })?;
    let target = combinations(3 + ploidy - 1, ploidy)
        .ok_or_else(|| invalid(format!("{context} cardinality overflows")))?;
    values.extend((values.len()..target).map(|_| None));
    Ok(())
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
