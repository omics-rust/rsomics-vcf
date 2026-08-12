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

use super::cardinality::{combinations, infer_ploidy};

#[derive(Clone, Copy)]
struct SampleProjection {
    alternate_count: usize,
    alternate: usize,
    ploidy: Option<usize>,
    sample: usize,
    keep_sum_ad: bool,
}

pub(super) fn validate(header: &Header, keep_sum_ad: bool) -> Result<()> {
    if !keep_sum_ad {
        return Ok(());
    }
    let schema = header
        .formats()
        .get("AD")
        .ok_or_else(|| invalid("--keep-sum AD requires FORMAT/AD in the header"))?;
    if schema.number() != format::Number::ReferenceAlternateBases
        || schema.ty() != format::Type::Integer
    {
        return Err(invalid(
            "--keep-sum AD requires FORMAT/AD Number=R,Type=Integer",
        ));
    }
    Ok(())
}

pub(super) fn split(
    header: &Header,
    record: &RecordBuf,
    keep_sum_ad: bool,
) -> Result<Vec<RecordBuf>> {
    let alternates = record.alternate_bases().as_ref();
    if alternates.len() < 2 {
        return Ok(vec![record.clone()]);
    }

    (0..alternates.len())
        .map(|alternate| split_one(header, record, alternate, keep_sum_ad))
        .collect()
}

pub(super) fn split_one(
    header: &Header,
    source: &RecordBuf,
    alternate: usize,
    keep_sum_ad: bool,
) -> Result<RecordBuf> {
    let alternate_count = source.alternate_bases().as_ref().len();
    let mut record = source.clone();
    *record.alternate_bases_mut() =
        vec![source.alternate_bases().as_ref()[alternate].clone()].into();

    for (key, schema) in header.infos() {
        let Some(value) = record.info_mut().get_mut(key) else {
            continue;
        };
        *value = project_info(
            value.take(),
            schema.number(),
            alternate_count,
            alternate,
            key,
        )?;
    }

    let keys = source.samples().keys().clone();
    let genotype_index = keys.as_ref().get_index_of("GT");
    let mut samples = Vec::new();
    for (sample_index, sample) in source.samples().values().enumerate() {
        let mut values = sample.values().to_vec();
        let ploidy = genotype_index
            .and_then(|index| values.get(index))
            .and_then(Option::as_ref)
            .map(|value| genotype_ploidy(value, alternate_count, sample_index))
            .transpose()?;

        for (index, key) in keys.as_ref().iter().enumerate() {
            let Some(value) = values.get_mut(index) else {
                continue;
            };
            if key == "GT" {
                remap_genotype(value, alternate_count, alternate, sample_index)?;
                continue;
            }
            let schema = header
                .formats()
                .get(key)
                .ok_or_else(|| invalid(format!("FORMAT/{key} is absent from the header")))?;
            *value = project_sample(
                value.take(),
                schema.number(),
                key,
                SampleProjection {
                    alternate_count,
                    alternate,
                    ploidy,
                    sample: sample_index,
                    keep_sum_ad,
                },
            )?;
        }
        samples.push(values);
    }
    *record.samples_mut() = Samples::new(keys, samples);
    Ok(record)
}

fn project_info(
    value: Option<InfoValue>,
    number: info::Number,
    alternate_count: usize,
    alternate: usize,
    key: &str,
) -> Result<Option<InfoValue>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let indexes = match number {
        info::Number::AlternateBases => Some(vec![alternate]),
        info::Number::ReferenceAlternateBases => Some(vec![0, alternate + 1]),
        info::Number::Samples => None,
        info::Number::Count(_) | info::Number::Unknown => return Ok(Some(value)),
    };
    let context = format!("INFO/{key}");
    match (value, indexes) {
        (InfoValue::Array(array), Some(indexes)) => project_info_array(
            array,
            &indexes,
            expected(number, alternate_count)?,
            &context,
        )
        .map(InfoValue::Array)
        .map(Some),
        (InfoValue::Array(array), None) => {
            project_info_genotypes(array, alternate_count, alternate, &context)
                .map(InfoValue::Array)
                .map(Some)
        }
        _ => Err(invalid(format!("{context} is not encoded as an array"))),
    }
}

fn project_sample(
    value: Option<SampleValue>,
    number: format::Number,
    key: &str,
    projection: SampleProjection,
) -> Result<Option<SampleValue>> {
    let Some(value) = value else {
        return Ok(None);
    };
    if projection.keep_sum_ad && key == "AD" {
        return project_ad_sum(value, projection).map(Some);
    }
    let indexes = match number {
        format::Number::AlternateBases => Some(vec![projection.alternate]),
        format::Number::ReferenceAlternateBases => Some(vec![0, projection.alternate + 1]),
        format::Number::Samples => None,
        format::Number::Count(_)
        | format::Number::LocalAlternateBases
        | format::Number::LocalReferenceAlternateBases
        | format::Number::LocalSamples
        | format::Number::Ploidy
        | format::Number::BaseModifications
        | format::Number::Unknown => return Ok(Some(value)),
    };
    let context = format!("FORMAT/{key} sample {}", projection.sample + 1);
    match (value, indexes) {
        (SampleValue::Array(array), Some(indexes)) => project_sample_array(
            array,
            &indexes,
            expected_format(number, projection.alternate_count)?,
            &context,
        )
        .map(SampleValue::Array)
        .map(Some),
        (SampleValue::Array(array), None) => project_sample_genotypes(
            array,
            projection.alternate_count,
            projection.alternate,
            projection.ploidy,
            &context,
        )
        .map(SampleValue::Array)
        .map(Some),
        _ => Err(invalid(format!("{context} is not encoded as an array"))),
    }
}

fn project_ad_sum(value: SampleValue, projection: SampleProjection) -> Result<SampleValue> {
    let SampleValue::Array(SampleArray::Integer(values)) = value else {
        return Err(invalid(format!(
            "FORMAT/AD sample {} is not encoded as an integer array",
            projection.sample + 1
        )));
    };
    let expected = projection.alternate_count + 1;
    if values.len() != expected {
        return Err(invalid(format!(
            "FORMAT/AD sample {} has {} values, expected {expected}",
            projection.sample + 1,
            values.len()
        )));
    }
    let selected = projection.alternate + 1;
    let reference = values
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != selected)
        .filter_map(|(_, value)| *value)
        .try_fold(0i32, |sum, value| {
            sum.checked_add(value).ok_or_else(|| {
                invalid(format!(
                    "FORMAT/AD sample {} sum exceeds int32",
                    projection.sample + 1
                ))
            })
        })?;
    Ok(SampleValue::Array(SampleArray::Integer(vec![
        Some(reference),
        values[selected],
    ])))
}

fn remap_genotype(
    value: &mut Option<SampleValue>,
    alternate_count: usize,
    alternate: usize,
    sample: usize,
) -> Result<()> {
    let Some(value) = value else {
        return Ok(());
    };
    let SampleValue::Genotype(genotype) = value else {
        return Err(invalid(format!(
            "FORMAT/GT sample {} is not encoded as a genotype",
            sample + 1
        )));
    };
    for allele in genotype.as_mut() {
        let Some(position) = allele.position() else {
            continue;
        };
        if position > alternate_count {
            return Err(invalid(format!(
                "FORMAT/GT sample {} allele {position} exceeds {alternate_count} ALT alleles",
                sample + 1
            )));
        }
        *allele.position_mut() = Some(usize::from(position == alternate + 1));
    }
    Ok(())
}

fn genotype_ploidy(value: &SampleValue, alternate_count: usize, sample: usize) -> Result<usize> {
    let SampleValue::Genotype(genotype) = value else {
        return Err(invalid(format!(
            "FORMAT/GT sample {} is not encoded as a genotype",
            sample + 1
        )));
    };
    for allele in genotype.as_ref() {
        if allele
            .position()
            .is_some_and(|position| position > alternate_count)
        {
            return Err(invalid(format!(
                "FORMAT/GT sample {} has an allele outside 0..={alternate_count}",
                sample + 1
            )));
        }
    }
    Ok(genotype.as_ref().len())
}

fn project_info_array(
    array: InfoArray,
    indexes: &[usize],
    expected: usize,
    context: &str,
) -> Result<InfoArray> {
    match array {
        InfoArray::Integer(values) => {
            project(values, indexes, expected, context).map(InfoArray::Integer)
        }
        InfoArray::Float(values) => {
            project(values, indexes, expected, context).map(InfoArray::Float)
        }
        InfoArray::Character(values) => {
            project(values, indexes, expected, context).map(InfoArray::Character)
        }
        InfoArray::String(values) => {
            project(values, indexes, expected, context).map(InfoArray::String)
        }
    }
}

fn project_sample_array(
    array: SampleArray,
    indexes: &[usize],
    expected: usize,
    context: &str,
) -> Result<SampleArray> {
    match array {
        SampleArray::Integer(values) => {
            project(values, indexes, expected, context).map(SampleArray::Integer)
        }
        SampleArray::Float(values) => {
            project(values, indexes, expected, context).map(SampleArray::Float)
        }
        SampleArray::Character(values) => {
            project(values, indexes, expected, context).map(SampleArray::Character)
        }
        SampleArray::String(values) => {
            project(values, indexes, expected, context).map(SampleArray::String)
        }
    }
}

fn project_info_genotypes(
    array: InfoArray,
    alternate_count: usize,
    alternate: usize,
    context: &str,
) -> Result<InfoArray> {
    match array {
        InfoArray::Integer(values) => {
            project_genotypes(values, alternate_count, alternate, None, context)
                .map(InfoArray::Integer)
        }
        InfoArray::Float(values) => {
            project_genotypes(values, alternate_count, alternate, None, context)
                .map(InfoArray::Float)
        }
        InfoArray::Character(values) => {
            project_genotypes(values, alternate_count, alternate, None, context)
                .map(InfoArray::Character)
        }
        InfoArray::String(values) => {
            project_genotypes(values, alternate_count, alternate, None, context)
                .map(InfoArray::String)
        }
    }
}

fn project_sample_genotypes(
    array: SampleArray,
    alternate_count: usize,
    alternate: usize,
    ploidy: Option<usize>,
    context: &str,
) -> Result<SampleArray> {
    match array {
        SampleArray::Integer(values) => {
            project_genotypes(values, alternate_count, alternate, ploidy, context)
                .map(SampleArray::Integer)
        }
        SampleArray::Float(values) => {
            project_genotypes(values, alternate_count, alternate, ploidy, context)
                .map(SampleArray::Float)
        }
        SampleArray::Character(values) => {
            project_genotypes(values, alternate_count, alternate, ploidy, context)
                .map(SampleArray::Character)
        }
        SampleArray::String(values) => {
            project_genotypes(values, alternate_count, alternate, ploidy, context)
                .map(SampleArray::String)
        }
    }
}

fn project<T: Clone>(
    values: Vec<Option<T>>,
    indexes: &[usize],
    expected: usize,
    context: &str,
) -> Result<Vec<Option<T>>> {
    if values.len() != expected {
        return Err(invalid(format!(
            "{context} has {} values, expected {expected}",
            values.len()
        )));
    }
    Ok(indexes.iter().map(|&index| values[index].clone()).collect())
}

fn project_genotypes<T: Clone>(
    values: Vec<Option<T>>,
    alternate_count: usize,
    alternate: usize,
    ploidy: Option<usize>,
    context: &str,
) -> Result<Vec<Option<T>>> {
    let allele_count = alternate_count + 1;
    let ploidy = match ploidy {
        Some(ploidy) => ploidy,
        None => infer_ploidy(allele_count, values.len()).ok_or_else(|| {
            invalid(format!(
                "{context} cardinality {} does not identify a ploidy for {allele_count} alleles",
                values.len()
            ))
        })?,
    };
    let expected = combinations(allele_count + ploidy - 1, ploidy)
        .ok_or_else(|| invalid(format!("{context} cardinality overflows")))?;
    if values.len() != expected {
        return Err(invalid(format!(
            "{context} has {} values, expected {expected} for ploidy {ploidy}",
            values.len()
        )));
    }
    let selected = alternate + 1;
    (0..=ploidy)
        .map(|copies| {
            genotype_index(ploidy, selected, copies)
                .and_then(|index| values.get(index).cloned())
                .ok_or_else(|| invalid(format!("{context} genotype index overflows")))
        })
        .collect()
}

fn genotype_index(ploidy: usize, selected: usize, copies: usize) -> Option<usize> {
    (0..ploidy).try_fold(0usize, |index, allele_index| {
        let allele = usize::from(allele_index >= ploidy - copies) * selected;
        combinations(allele + allele_index, allele_index + 1)
            .and_then(|offset| index.checked_add(offset))
    })
}

fn expected(number: info::Number, alternate_count: usize) -> Result<usize> {
    match number {
        info::Number::AlternateBases => Ok(alternate_count),
        info::Number::ReferenceAlternateBases => Ok(alternate_count + 1),
        _ => Err(invalid("invalid INFO projection cardinality")),
    }
}

fn expected_format(number: format::Number, alternate_count: usize) -> Result<usize> {
    match number {
        format::Number::AlternateBases => Ok(alternate_count),
        format::Number::ReferenceAlternateBases => Ok(alternate_count + 1),
        _ => Err(invalid("invalid FORMAT projection cardinality")),
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::io::Write as _};

    use super::*;

    #[test]
    fn remaps_typed_allele_fields_and_mixed_ploidy() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=FA,Number=A,Type=Integer,Description=\"A\">\n\
##FORMAT=<ID=FR,Number=R,Type=Integer,Description=\"R\">\n\
##FORMAT=<ID=FG,Number=G,Type=Integer,Description=\"G\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(
            b"chr1\t10\t.\tA\tC,G\t.\tPASS\tIA=10,20;IR=5,3,2;IG=0,10,20,30,40,50\tGT:FA:FR:FG\t1/2:11,22:7,4,3:0,10,20,30,40,50\t2:33,44:8,5,6:0,10,20"
                .as_slice(),
        )
        .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();

        let records = split(&header, &record, false).unwrap();
        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        assert_eq!(
            String::from_utf8(writer.into_inner()).unwrap(),
            "chr1\t10\t.\tA\tC\t.\tPASS\tIA=10;IR=5,3;IG=0,10,20\tGT:FA:FR:FG\t1/0:11:7,4:0,10,20\t0:33:8,5:0,10\n\
chr1\t10\t.\tA\tG\t.\tPASS\tIA=20;IR=5,2;IG=0,30,50\tGT:FA:FR:FG\t0/1:22:7,3:0,30,50\t1:44:8,6:0,20\n"
        );
    }

    #[test]
    fn malformed_allele_cardinality_fails() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(b"chr1\t10\t.\tA\tC,G\t.\tPASS\tIA=10".as_slice()).unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let error = split(&header, &record, false).unwrap_err().to_string();
        assert!(
            error.contains("INFO/IA has 1 values, expected 2"),
            "{error}"
        );
    }

    #[test]
    fn preserves_ad_sums_when_requested() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        let raw =
            vcf::Record::try_from(b"chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:10,3,2".as_slice())
                .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let records = split(&header, &record, true).unwrap();
        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        let output = String::from_utf8(writer.into_inner()).unwrap();
        assert!(output.contains("GT:AD\t1/0:12,3"), "{output}");
        assert!(output.contains("GT:AD\t0/1:13,2"), "{output}");
    }
}
