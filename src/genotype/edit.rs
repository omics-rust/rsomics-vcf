use noodles_vcf::variant::{
    RecordBuf,
    record_buf::{
        Samples,
        samples::sample::{Value as SampleValue, value::Genotype},
    },
};
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MissingPolicy {
    Ignore,
    Error,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Change {
    pub(crate) genotypes: u64,
    pub(crate) alleles: u64,
}

pub(crate) fn edit_selected<F>(
    record: &mut RecordBuf,
    selected: &[bool],
    missing: MissingPolicy,
    mut edit: F,
) -> Result<Change>
where
    F: FnMut(usize, &Genotype) -> Result<Genotype>,
{
    let sample_count = record.samples().values().count();
    if selected.len() != sample_count {
        return Err(RsomicsError::InvalidInput(format!(
            "genotype selection has {} samples but the record has {sample_count}",
            selected.len()
        )));
    }

    let keys = record.samples().keys().clone();
    let Some(genotype_index) = keys.as_ref().get_index_of("GT") else {
        return match missing {
            MissingPolicy::Ignore => Ok(Change::default()),
            MissingPolicy::Error => Err(RsomicsError::InvalidInput(
                "record has no FORMAT/GT field".to_owned(),
            )),
        };
    };
    let mut rows = record
        .samples()
        .values()
        .map(|sample| sample.values().to_vec())
        .collect::<Vec<_>>();
    let mut change = Change::default();

    for (sample, (row, selected)) in rows.iter_mut().zip(selected).enumerate() {
        if !selected {
            continue;
        }
        let Some(value) = row.get_mut(genotype_index) else {
            if missing == MissingPolicy::Error {
                return Err(missing_value(sample));
            }
            continue;
        };
        let Some(value) = value else {
            if missing == MissingPolicy::Error {
                return Err(missing_value(sample));
            }
            continue;
        };
        let SampleValue::Genotype(genotype) = value else {
            return Err(RsomicsError::InvalidInput(format!(
                "sample {} FORMAT/GT is not encoded as a genotype",
                sample + 1
            )));
        };

        let replacement = edit(sample, genotype)?;
        if replacement != *genotype {
            change.genotypes = change.genotypes.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("changed genotype count exceeds u64".to_owned())
            })?;
            let changed = changed_alleles(genotype, &replacement);
            change.alleles = change
                .alleles
                .checked_add(u64::try_from(changed).map_err(|_| {
                    RsomicsError::InvalidInput("changed allele count exceeds u64".to_owned())
                })?)
                .ok_or_else(|| {
                    RsomicsError::InvalidInput("changed allele count exceeds u64".to_owned())
                })?;
            *genotype = replacement;
        }
    }

    *record.samples_mut() = Samples::new(keys, rows);
    Ok(change)
}

fn changed_alleles(left: &Genotype, right: &Genotype) -> usize {
    let shared = left.as_ref().len().min(right.as_ref().len());
    left.as_ref()[..shared]
        .iter()
        .zip(&right.as_ref()[..shared])
        .filter(|(left, right)| left != right)
        .count()
        + left.as_ref().len().abs_diff(right.as_ref().len())
}

fn missing_value(sample: usize) -> RsomicsError {
    RsomicsError::InvalidInput(format!("sample {} has no FORMAT/GT value", sample + 1))
}
