use noodles_vcf::{
    self as vcf,
    header::record::value::map::info::{Number, Type},
    variant::{
        RecordBuf,
        record_buf::{
            info::field::{Value as InfoValue, value::Array},
            samples::sample::Value as SampleValue,
        },
    },
};
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum InfoPolicy {
    BestEffort,
    Strict,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct AlleleCounts {
    pub(crate) counts: Vec<u64>,
    pub(crate) total: u64,
}

pub(crate) fn allele_counts(record: &RecordBuf) -> Result<AlleleCounts> {
    let allele_count = record
        .alternate_bases()
        .as_ref()
        .len()
        .checked_add(1)
        .ok_or_else(|| RsomicsError::InvalidInput("allele count exceeds usize".to_owned()))?;
    let Some(genotypes) = record.samples().select("GT") else {
        return Err(RsomicsError::InvalidInput(
            "record has no FORMAT/GT field".to_owned(),
        ));
    };
    let mut counts = vec![0u64; allele_count];
    let mut total = 0u64;

    for sample in 0..record.samples().values().count() {
        let Some(Some(value)) = genotypes.get(sample) else {
            continue;
        };
        let SampleValue::Genotype(genotype) = value else {
            return Err(RsomicsError::InvalidInput(format!(
                "sample {} FORMAT/GT is not encoded as a genotype",
                sample + 1
            )));
        };
        for allele in genotype.as_ref() {
            let Some(position) = allele.position() else {
                continue;
            };
            let Some(count) = counts.get_mut(position) else {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample {} genotype allele index {position} exceeds {} ALT alleles",
                    sample + 1,
                    allele_count - 1
                )));
            };
            *count = count
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("allele count exceeds u64".to_owned()))?;
            total = total
                .checked_add(1)
                .ok_or_else(|| RsomicsError::InvalidInput("allele total exceeds u64".to_owned()))?;
        }
    }

    Ok(AlleleCounts { counts, total })
}

pub(crate) fn reconcile_ac_an(
    header: &vcf::Header,
    record: &mut RecordBuf,
    policy: InfoPolicy,
) -> Result<()> {
    let has_ac = record.info().as_ref().contains_key("AC");
    let has_an = record.info().as_ref().contains_key("AN");
    if !has_ac && !has_an {
        return Ok(());
    }

    let update_ac = has_ac && validate_ac(header, record, policy)?;
    let update_an = has_an && validate_an(header, record, policy)?;
    if !update_ac && !update_an {
        return Ok(());
    }
    let counts = match allele_counts(record) {
        Ok(counts) => counts,
        Err(_) if policy == InfoPolicy::BestEffort => return Ok(()),
        Err(error) => return Err(error),
    };

    let ac = update_ac
        .then(|| {
            counts.counts[1..]
                .iter()
                .map(|count| {
                    i32::try_from(*count).map(Some).map_err(|_| {
                        RsomicsError::InvalidInput("INFO/AC count exceeds int32".to_owned())
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .transpose();
    let an = update_an.then(|| {
        i32::try_from(counts.total)
            .map_err(|_| RsomicsError::InvalidInput("INFO/AN count exceeds int32".to_owned()))
    });
    let (ac, an) = match (ac, an.transpose()) {
        (Ok(ac), Ok(an)) => (ac, an),
        (Err(error), _) | (_, Err(error)) if policy == InfoPolicy::Strict => return Err(error),
        (ac, an) => (ac.ok().flatten(), an.ok().flatten()),
    };

    if has_ac && let Some(ac) = ac {
        record
            .info_mut()
            .insert("AC".to_owned(), Some(InfoValue::Array(Array::Integer(ac))));
    }
    if has_an && let Some(an) = an {
        record
            .info_mut()
            .insert("AN".to_owned(), Some(InfoValue::Integer(an)));
    }
    Ok(())
}

fn validate_ac(header: &vcf::Header, record: &RecordBuf, policy: InfoPolicy) -> Result<bool> {
    let definition = header.infos().get("AC");
    let definition_valid = definition.is_some_and(|definition| {
        definition.number() == Number::AlternateBases && definition.ty() == Type::Integer
    });
    let value_valid = matches!(
        record.info().get("AC"),
        Some(Some(InfoValue::Array(Array::Integer(values))))
            if values.len() == record.alternate_bases().as_ref().len()
                && values.iter().all(|value| value.is_some_and(|value| value >= 0))
    );
    validate_contract(
        definition_valid && value_valid,
        policy,
        "INFO/AC must have Number=A, Type=Integer, and one nonnegative value per ALT allele",
    )
}

fn validate_an(header: &vcf::Header, record: &RecordBuf, policy: InfoPolicy) -> Result<bool> {
    let definition = header.infos().get("AN");
    let definition_valid = definition.is_some_and(|definition| {
        definition.number() == Number::Count(1) && definition.ty() == Type::Integer
    });
    let value_valid = matches!(
        record.info().get("AN"),
        Some(Some(InfoValue::Integer(value))) if *value >= 0
    );
    validate_contract(
        definition_valid && value_valid,
        policy,
        "INFO/AN must have Number=1, Type=Integer, and a nonnegative integer value",
    )
}

fn validate_contract(valid: bool, policy: InfoPolicy, message: &str) -> Result<bool> {
    if valid {
        Ok(true)
    } else if policy == InfoPolicy::BestEffort {
        Ok(false)
    } else {
        Err(RsomicsError::InvalidInput(message.to_owned()))
    }
}
