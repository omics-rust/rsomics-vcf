use noodles_vcf::variant::{
    RecordBuf,
    record_buf::{
        AlternateBases, Samples,
        info::field::{Value as InfoValue, value::Array as InfoArray},
        samples::sample::Value as SampleValue,
    },
};
use rsomics_common::Result;

use super::{invalid, is_ordinary, iupac_first};

pub(super) fn fix_reference(
    record: &mut RecordBuf,
    alleles: &mut [Vec<u8>],
    expected: &[u8],
) -> Result<()> {
    let alternate = alleles
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, allele)| {
            is_ordinary(allele)
                && allele.len() == expected.len()
                && allele.eq_ignore_ascii_case(expected)
        })
        .map(|(index, _)| index);
    if let Some(alternate) = alternate {
        alleles.swap(0, alternate);
        set_alleles(record, alleles);
        let reference_copies = swap_genotypes(record, alternate)?;
        update_allele_count(record, alternate, reference_copies)?;
        return Ok(());
    }

    let old_reference = alleles[0].clone();
    for allele in alleles
        .iter_mut()
        .skip(1)
        .filter(|allele| is_ordinary(allele))
    {
        for ((base, old), new) in allele.iter_mut().zip(&old_reference).zip(expected) {
            if base.eq_ignore_ascii_case(old) {
                *base = *new;
            }
        }
    }
    alleles[0].clone_from_slice(expected);
    set_alleles(record, alleles);
    Ok(())
}

pub(super) fn set_alleles(record: &mut RecordBuf, alleles: &[Vec<u8>]) {
    *record.reference_bases_mut() = String::from_utf8(alleles[0].clone()).unwrap();
    *record.alternate_bases_mut() = AlternateBases::from(
        alleles[1..]
            .iter()
            .map(|allele| String::from_utf8(allele.clone()).unwrap())
            .collect::<Vec<_>>(),
    );
}

pub(super) fn remove_reference_alternates(
    record: &mut RecordBuf,
    alleles: &mut Vec<Vec<u8>>,
) -> Result<bool> {
    let mut mapping = Vec::with_capacity(alleles.len());
    mapping.push(0);
    let mut retained = Vec::with_capacity(alleles.len());
    retained.push(alleles[0].clone());
    for allele in alleles.iter().skip(1) {
        if allele.eq_ignore_ascii_case(&alleles[0]) {
            mapping.push(0);
        } else {
            mapping.push(retained.len());
            retained.push(allele.clone());
        }
    }
    if retained.len() == alleles.len() {
        return Ok(false);
    }
    *alleles = retained;
    set_alleles(record, alleles);
    remap_genotypes(record, &mapping)?;
    Ok(true)
}

pub(super) fn clean_iupac(allele: &mut [u8]) -> bool {
    if !is_ordinary(allele) {
        return false;
    }
    let mut changed = false;
    for base in allele {
        let replacement = iupac_first(*base);
        let replacement = if base.is_ascii_lowercase() {
            replacement.to_ascii_lowercase()
        } else {
            replacement
        };
        if replacement != *base {
            *base = replacement;
            changed = true;
        }
    }
    changed
}

pub(super) fn resolve_unknown_reference_bases(alleles: &mut [Vec<u8>], expected: &[u8]) -> bool {
    let mut changed = false;
    for (index, expected) in expected.iter().copied().enumerate() {
        if !alleles[0][index].eq_ignore_ascii_case(&b'N') || expected == b'N' {
            continue;
        }
        alleles[0][index] = expected;
        changed = true;
        for allele in alleles
            .iter_mut()
            .skip(1)
            .filter(|allele| is_ordinary(allele))
        {
            if allele
                .get(index)
                .is_some_and(|base| base.eq_ignore_ascii_case(&b'N'))
            {
                allele[index] = expected;
            }
        }
    }
    changed
}

fn swap_genotypes(record: &mut RecordBuf, alternate: usize) -> Result<usize> {
    let keys = record.samples().keys().clone();
    let Some(genotype_index) = keys.as_ref().get_index_of("GT") else {
        return Ok(0);
    };
    let alternate_count = record.alternate_bases().as_ref().len();
    let mut reference_copies = 0usize;
    let mut samples = Vec::new();
    for (sample_index, sample) in record.samples().values().enumerate() {
        let mut values = sample.values().to_vec();
        let Some(value) = values.get_mut(genotype_index) else {
            return Err(invalid(
                record,
                &format!("FORMAT/GT sample {} is missing", sample_index + 1),
            ));
        };
        if let Some(SampleValue::Genotype(genotype)) = value {
            for allele in genotype.as_mut() {
                let Some(position) = allele.position() else {
                    continue;
                };
                if position > alternate_count {
                    return Err(invalid(
                        record,
                        &format!(
                            "FORMAT/GT sample {} allele {position} exceeds {alternate_count} ALT alleles",
                            sample_index + 1
                        ),
                    ));
                }
                if position == 0 {
                    *allele.position_mut() = Some(alternate);
                    reference_copies += 1;
                } else if position == alternate {
                    *allele.position_mut() = Some(0);
                }
            }
        }
        samples.push(values);
    }
    *record.samples_mut() = Samples::new(keys, samples);
    Ok(reference_copies)
}

fn remap_genotypes(record: &mut RecordBuf, mapping: &[usize]) -> Result<()> {
    let keys = record.samples().keys().clone();
    let Some(genotype_index) = keys.as_ref().get_index_of("GT") else {
        return Ok(());
    };
    let alternate_count = mapping.len() - 1;
    let mut samples = Vec::new();
    for (sample_index, sample) in record.samples().values().enumerate() {
        let mut values = sample.values().to_vec();
        let Some(value) = values.get_mut(genotype_index) else {
            return Err(invalid(
                record,
                &format!("FORMAT/GT sample {} is missing", sample_index + 1),
            ));
        };
        if let Some(SampleValue::Genotype(genotype)) = value {
            for allele in genotype.as_mut() {
                let Some(position) = allele.position() else {
                    continue;
                };
                let mapped = mapping.get(position).ok_or_else(|| {
                    invalid(
                        record,
                        &format!(
                            "FORMAT/GT sample {} allele {position} exceeds {alternate_count} ALT alleles",
                            sample_index + 1
                        ),
                    )
                })?;
                *allele.position_mut() = Some(*mapped);
            }
        }
        samples.push(values);
    }
    *record.samples_mut() = Samples::new(keys, samples);
    Ok(())
}

fn update_allele_count(record: &mut RecordBuf, alternate: usize, count: usize) -> Result<()> {
    let count =
        i32::try_from(count).map_err(|_| invalid(record, "reference allele count exceeds i32"))?;
    let Some(value) = record.info_mut().get_mut("AC") else {
        return Ok(());
    };
    match value {
        Some(InfoValue::Integer(value)) if alternate == 1 => *value = count,
        Some(InfoValue::Array(InfoArray::Integer(values))) => {
            if let Some(value) = values.get_mut(alternate - 1) {
                *value = Some(count);
            }
        }
        _ => {}
    }
    Ok(())
}
