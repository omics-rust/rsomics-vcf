use noodles_core::Position;
use noodles_vcf::variant::{RecordBuf, record_buf::AlternateBases};
use rsomics_common::{Result, RsomicsError};

pub(super) fn atomize(record: RecordBuf) -> Result<(Vec<RecordBuf>, bool)> {
    let alternates = record.alternate_bases().as_ref();
    if alternates.len() != 1 {
        return Ok((vec![record], false));
    }
    let reference = record.reference_bases().as_bytes();
    let alternate = alternates[0].as_bytes();
    if reference.len() != alternate.len()
        || reference.len() <= 1
        || !is_sequence(reference)
        || !is_sequence(alternate)
    {
        return Ok((vec![record], false));
    }
    let start = record.variant_start().map(usize::from).ok_or_else(|| {
        RsomicsError::InvalidInput("atomizing a record without a position".to_owned())
    })?;
    let mut records = Vec::new();
    for (offset, (&reference, &alternate)) in reference.iter().zip(alternate).enumerate() {
        if reference.eq_ignore_ascii_case(&alternate) {
            continue;
        }
        let position = start.checked_add(offset).ok_or_else(|| {
            RsomicsError::InvalidInput("atomized position exceeds usize".to_owned())
        })?;
        let mut atom = record.clone();
        *atom.variant_start_mut() = Some(Position::try_from(position).map_err(|_| {
            RsomicsError::InvalidInput(
                "atomized position exceeds the VCF coordinate range".to_owned(),
            )
        })?);
        *atom.reference_bases_mut() = char::from(reference.to_ascii_uppercase()).to_string();
        *atom.alternate_bases_mut() =
            AlternateBases::from(vec![char::from(alternate.to_ascii_uppercase()).to_string()]);
        records.push(atom);
    }
    if records.is_empty() {
        Ok((vec![record], false))
    } else {
        Ok((records, true))
    }
}

fn is_sequence(allele: &[u8]) -> bool {
    allele
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::io::Write as _};

    use super::*;

    #[test]
    fn decomposes_mnvs_and_skips_matching_bases() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        let raw =
            vcf::Record::try_from(b"chr1\t20\t.\tACGT\tAGGA\t.\tPASS\tDP=7\tGT\t1/1".as_slice())
                .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let (records, changed) = atomize(record).unwrap();
        assert!(changed);

        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        assert_eq!(
            String::from_utf8(writer.into_inner()).unwrap(),
            "chr1\t21\t.\tC\tG\t.\tPASS\tDP=7\tGT\t1/1\n\
chr1\t23\t.\tT\tA\t.\tPASS\tDP=7\tGT\t1/1\n"
        );
    }
}
