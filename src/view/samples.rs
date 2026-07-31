use std::collections::HashSet;

use noodles_vcf::{
    self as vcf,
    header::record::value::{Map, map::Info},
    variant::{
        record::info::field::key,
        record_buf::{
            Samples,
            info::field::{Value as InfoValue, value::Array},
            samples::sample::Value as SampleValue,
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use super::Options;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SampleSelection {
    pub names: Vec<String>,
    pub exclude: bool,
    pub force: bool,
}

pub(super) struct Projection {
    indices: Option<Vec<usize>>,
}

impl Projection {
    pub(super) fn new(header: &vcf::Header, options: &Options) -> Result<Self> {
        let Some(selection) = &options.samples else {
            return Ok(Self { indices: None });
        };
        let requested: HashSet<_> = selection.names.iter().map(String::as_str).collect();
        if !selection.exclude && requested.len() != selection.names.len() {
            return Err(RsomicsError::InvalidInput(
                "sample inclusion list contains a duplicate name".to_owned(),
            ));
        }

        let mut missing = Vec::new();
        for name in &selection.names {
            if !header.sample_names().contains(name) {
                missing.push(name.as_str());
            }
        }
        if !missing.is_empty() && !selection.force {
            return Err(RsomicsError::InvalidInput(format!(
                "unknown sample{}: {}",
                if missing.len() == 1 { "" } else { "s" },
                missing.join(", ")
            )));
        }

        let indices = if selection.exclude {
            header
                .sample_names()
                .iter()
                .enumerate()
                .filter_map(|(index, name)| (!requested.contains(name.as_str())).then_some(index))
                .collect()
        } else {
            selection
                .names
                .iter()
                .filter_map(|name| header.sample_names().get_index_of(name))
                .collect()
        };
        Ok(Self {
            indices: Some(indices),
        })
    }

    pub(super) fn header(&self, header: &vcf::Header, options: &Options) -> vcf::Header {
        let mut output = header.clone();
        if options.drop_genotypes {
            output.sample_names_mut().clear();
            output.formats_mut().clear();
        } else if let Some(indices) = &self.indices {
            let names: Vec<_> = indices
                .iter()
                .filter_map(|&index| header.sample_names().get_index(index))
                .cloned()
                .collect();
            output.sample_names_mut().clear();
            output.sample_names_mut().extend(names);
            if indices.is_empty() {
                output.formats_mut().clear();
            }
        }

        if self.indices.is_some() && options.update_info {
            output
                .infos_mut()
                .entry(key::ALLELE_COUNT.into())
                .or_insert_with(|| Map::<Info>::from(key::ALLELE_COUNT));
            output
                .infos_mut()
                .entry(key::TOTAL_ALLELE_COUNT.into())
                .or_insert_with(|| Map::<Info>::from(key::TOTAL_ALLELE_COUNT));
        }
        output
    }

    pub(super) fn apply(
        &self,
        mut record: vcf::variant::RecordBuf,
        options: &Options,
    ) -> Result<vcf::variant::RecordBuf> {
        let Some(indices) = &self.indices else {
            if options.drop_genotypes {
                record.samples_mut().clear();
            }
            return Ok(record);
        };

        let keys = record.samples().keys().clone();
        let values: Vec<_> = indices
            .iter()
            .map(|&index| {
                record
                    .samples()
                    .get_index(index)
                    .map(|sample| sample.values().to_vec())
                    .unwrap_or_default()
            })
            .collect();
        *record.samples_mut() = Samples::new(keys, values);

        if options.update_info && record.samples().select("GT").is_some() {
            update_ac_an(&mut record)?;
        }
        if options.drop_genotypes {
            record.samples_mut().clear();
        }
        Ok(record)
    }
}

fn update_ac_an(record: &mut vcf::variant::RecordBuf) -> Result<()> {
    let allele_count = record.alternate_bases().as_ref().len() + 1;
    let mut counts = vec![0i32; allele_count];

    if let Some(genotypes) = record.samples().select("GT") {
        for index in 0..record.samples().values().count() {
            let Some(Some(value)) = genotypes.get(index) else {
                continue;
            };
            let SampleValue::Genotype(genotype) = value else {
                return Err(RsomicsError::InvalidInput(
                    "FORMAT/GT is not encoded as a genotype".to_owned(),
                ));
            };
            for allele in genotype.as_ref() {
                if let Some(position) = allele.position() {
                    let count = counts.get_mut(position).ok_or_else(|| {
                        RsomicsError::InvalidInput(format!(
                            "genotype allele index {position} exceeds {} ALT alleles",
                            allele_count - 1
                        ))
                    })?;
                    *count = count.checked_add(1).ok_or_else(|| {
                        RsomicsError::InvalidInput("allele count exceeds int32".to_owned())
                    })?;
                }
            }
        }
    }

    let an = counts.iter().try_fold(0i32, |total, count| {
        total.checked_add(*count).ok_or_else(|| {
            RsomicsError::InvalidInput("total allele count exceeds int32".to_owned())
        })
    })?;
    let ac = counts.into_iter().skip(1).map(Some).collect::<Vec<_>>();
    if ac.is_empty() {
        record.info_mut().as_mut().shift_remove(key::ALLELE_COUNT);
    } else {
        record.info_mut().insert(
            key::ALLELE_COUNT.into(),
            Some(InfoValue::Array(Array::Integer(ac))),
        );
    }
    record
        .info_mut()
        .insert(key::TOTAL_ALLELE_COUNT.into(), Some(InfoValue::Integer(an)));
    Ok(())
}
