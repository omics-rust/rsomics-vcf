mod fields;

use noodles_vcf::{
    Header,
    variant::{
        RecordBuf,
        record_buf::{AlternateBases, Filters},
    },
};
use rsomics_common::{Result, RsomicsError};

pub(super) fn join<'a>(
    header: &Header,
    records: impl IntoIterator<Item = &'a RecordBuf>,
) -> Result<RecordBuf> {
    let records: Vec<_> = records.into_iter().collect();
    let Some(first) = records.first().copied() else {
        return Err(invalid("cannot join an empty record group"));
    };
    if records.len() == 1 {
        return Ok(first.clone());
    }

    let (alleles, mappings) = merge_alleles(&records)?;
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
    merge_filters(&records, &mut output);
    fields::merge_info(header, &records, &mappings, &mut output)?;
    fields::merge_samples(header, &records, &mappings, &mut output)?;
    Ok(output)
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
