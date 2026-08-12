mod fields;

use noodles_vcf::{
    Header,
    variant::{
        RecordBuf,
        record_buf::{AlternateBases, Filters},
    },
};
use rsomics_common::{Result, RsomicsError};

use crate::variant_type;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Policy {
    Snps,
    Indels,
    Both,
    Any,
}

pub(super) fn join<'a>(
    policy: Policy,
    header: &Header,
    records: impl IntoIterator<Item = &'a RecordBuf>,
) -> Result<(Vec<(usize, RecordBuf)>, u64)> {
    let records: Vec<_> = records.into_iter().collect();
    if records.is_empty() {
        return Err(invalid("cannot join an empty record group"));
    }
    if records.len() == 1 {
        return Ok((vec![(0, records[0].clone())], 0));
    }
    if matches!(policy, Policy::Any) {
        return Ok((vec![(0, join_group(header, &records)?)], 1));
    }

    let mut order = records
        .iter()
        .enumerate()
        .map(|(index, record)| (index, *record, category(record)))
        .collect::<Vec<_>>();
    order.sort_by_key(|entry| entry.2);
    let targets: &[u32] = match policy {
        Policy::Snps => &[variant_type::SNP],
        Policy::Indels => &[variant_type::INDEL],
        Policy::Both => &[
            variant_type::SNP,
            variant_type::MNP,
            variant_type::INDEL,
            variant_type::OTHER,
        ],
        Policy::Any => unreachable!(),
    };
    let mut output = Vec::new();
    let mut joined = 0u64;
    let mut cursor = 0;
    while cursor < order.len() {
        let mut end = cursor;
        for target in targets {
            end = cursor;
            while end < order.len() && (order[end].2 == 0 || order[end].2 == *target) {
                end += 1;
            }
            if end > cursor {
                break;
            }
        }
        if end == cursor {
            end += 1;
        }
        let group = &order[cursor..end];
        let record = if group.len() == 1 {
            group[0].1.clone()
        } else {
            joined += 1;
            join_group(
                header,
                &group.iter().map(|entry| entry.1).collect::<Vec<_>>(),
            )?
        };
        output.push((group[0].0, record));
        cursor = end;
    }
    Ok((output, joined))
}

fn join_group(header: &Header, records: &[&RecordBuf]) -> Result<RecordBuf> {
    let first = records[0];

    let (alleles, mappings) = merge_alleles(records)?;
    let mut output = first.clone();
    *output.reference_bases_mut() = alleles[0].clone();
    *output.alternate_bases_mut() = AlternateBases::from(alleles[1..].to_vec());
    output.ids_mut().extend(
        records
            .iter()
            .skip(1)
            .flat_map(|record| record.ids().as_ref().iter().cloned()),
    );
    *output.quality_score_mut() = records
        .iter()
        .filter_map(|record| record.quality_score())
        .reduce(f32::max);
    merge_filters(records, &mut output);
    fields::merge_info(header, records, &mappings, &mut output)?;
    fields::merge_samples(header, records, &mappings, &mut output)?;
    Ok(output)
}

fn category(record: &RecordBuf) -> u32 {
    variant_type::record_mask(record) & !variant_type::REF
}

fn merge_alleles(records: &[&RecordBuf]) -> Result<(Vec<String>, Vec<Vec<usize>>)> {
    let mut alleles = std::iter::once(records[0].reference_bases().to_owned())
        .chain(records[0].alternate_bases().as_ref().iter().cloned())
        .collect::<Vec<_>>();
    let mut mappings = vec![(0..alleles.len()).collect()];
    for record in &records[1..] {
        let mut source = std::iter::once(record.reference_bases().to_owned())
            .chain(record.alternate_bases().as_ref().iter().cloned())
            .collect::<Vec<_>>();
        let prefix_length = source[0].len().min(alleles[0].len());
        let source_prefix = &source[0][..prefix_length];
        let target_prefix = &alleles[0][..prefix_length];
        if source_prefix != target_prefix {
            if !source_prefix.eq_ignore_ascii_case(target_prefix) {
                return Err(invalid(format!(
                    "cannot join incompatible REF alleles {} and {} at {}:{}",
                    source[0],
                    alleles[0],
                    record.reference_sequence_name(),
                    record.variant_start().map_or(0, usize::from)
                )));
            }
            alleles
                .iter_mut()
                .for_each(|allele| allele.make_ascii_uppercase());
            source
                .iter_mut()
                .for_each(|allele| allele.make_ascii_uppercase());
        }
        if source[0].len() > alleles[0].len() {
            let suffix = source[0][alleles[0].len()..].to_owned();
            expand_alleles(&mut alleles, &suffix);
        } else if alleles[0].len() > source[0].len() {
            let suffix = alleles[0][source[0].len()..].to_owned();
            expand_alleles(&mut source, &suffix);
        }
        if source[0] != alleles[0] {
            return Err(invalid(format!(
                "cannot join incompatible REF alleles {} and {} at {}:{}",
                source[0],
                alleles[0],
                record.reference_sequence_name(),
                record.variant_start().map_or(0, usize::from)
            )));
        }
        let mut mapping = vec![0];
        for alternate in source.into_iter().skip(1) {
            let index = alleles
                .iter()
                .enumerate()
                .skip(1)
                .find_map(|(index, allele)| {
                    allele.eq_ignore_ascii_case(&alternate).then_some(index)
                })
                .unwrap_or_else(|| {
                    alleles.push(alternate);
                    alleles.len() - 1
                });
            mapping.push(index);
        }
        mappings.push(mapping);
    }
    Ok((alleles, mappings))
}

fn expand_alleles(alleles: &mut [String], suffix: &str) {
    for allele in alleles {
        if !allele.starts_with('<') && !allele.starts_with('*') {
            allele.push_str(suffix);
        }
    }
}

fn merge_filters(records: &[&RecordBuf], output: &mut RecordBuf) {
    let mut filters: Filters = records[0]
        .filters()
        .as_ref()
        .iter()
        .filter(|filter| filter.as_str() != "PASS")
        .cloned()
        .collect();
    filters.extend(
        records
            .iter()
            .skip(1)
            .flat_map(|record| record.filters().as_ref())
            .filter(|filter| filter.as_str() != "PASS")
            .cloned(),
    );
    *output.filters_mut() = if filters.as_ref().is_empty() {
        Filters::pass()
    } else {
        filters
    };
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
