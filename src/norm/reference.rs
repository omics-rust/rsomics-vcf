use std::path::Path;

use noodles_core::Position;
use noodles_vcf::variant::{RecordBuf, record_buf::AlternateBases};
use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::IndexedFasta;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Outcome {
    Changed,
    Unchanged,
    Unsupported,
}

pub(crate) struct ReferenceNormalizer {
    reference: IndexedFasta,
}

impl ReferenceNormalizer {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        IndexedFasta::open(path).map(|reference| Self { reference })
    }

    pub(crate) fn normalize(&mut self, record: &mut RecordBuf) -> Result<Outcome> {
        let Some(position) = record.variant_start() else {
            return Ok(Outcome::Unsupported);
        };
        let mut alleles = Vec::with_capacity(record.alternate_bases().as_ref().len() + 1);
        alleles.push(record.reference_bases().as_bytes().to_ascii_uppercase());
        for alternate in record.alternate_bases().as_ref() {
            alleles.push(alternate.as_bytes().to_ascii_uppercase());
        }
        if alleles.len() == 1 || alleles.iter().any(|allele| !is_sequence(allele)) {
            return Ok(Outcome::Unsupported);
        }

        let original_position = usize::from(position) - 1;
        let mut position = original_position;
        let reference_end = position
            .checked_add(alleles[0].len())
            .ok_or_else(|| invalid(record, "REF range overflows"))?;
        let expected = self.reference.fetch(
            record.reference_sequence_name().as_bytes(),
            position..reference_end,
        )?;
        if !alleles[0].eq_ignore_ascii_case(expected) {
            return Err(invalid(
                record,
                &format!(
                    "REF {} does not match indexed reference {}",
                    String::from_utf8_lossy(&alleles[0]),
                    String::from_utf8_lossy(expected)
                ),
            ));
        }

        loop {
            let last = alleles[0].last().copied().unwrap();
            if alleles
                .iter()
                .skip(1)
                .any(|allele| allele.last().copied() != Some(last))
            {
                break;
            }
            let minimum = alleles.iter().map(Vec::len).min().unwrap();
            if minimum <= 1 && position == 0 {
                break;
            }
            for allele in &mut alleles {
                allele.pop();
            }
            if alleles.iter().any(Vec::is_empty) {
                let previous = self.reference.fetch(
                    record.reference_sequence_name().as_bytes(),
                    position - 1..position,
                )?[0]
                    .to_ascii_uppercase();
                for allele in &mut alleles {
                    allele.insert(0, previous);
                }
                position -= 1;
            }
        }

        loop {
            let minimum = alleles.iter().map(Vec::len).min().unwrap();
            if minimum <= 1 {
                break;
            }
            let first = alleles[0][0];
            if alleles.iter().skip(1).any(|allele| allele[0] != first) {
                break;
            }
            for allele in &mut alleles {
                allele.remove(0);
            }
            position += 1;
        }

        let changed = position != original_position
            || alleles[0] != record.reference_bases().as_bytes()
            || alleles[1..].iter().map(Vec::as_slice).ne(record
                .alternate_bases()
                .as_ref()
                .iter()
                .map(String::as_bytes));
        if !changed {
            return Ok(Outcome::Unchanged);
        }
        *record.variant_start_mut() = Some(Position::try_from(position + 1).map_err(|_| {
            invalid(
                record,
                "normalized position exceeds the VCF coordinate range",
            )
        })?);
        *record.reference_bases_mut() = String::from_utf8(alleles.remove(0)).unwrap();
        *record.alternate_bases_mut() = AlternateBases::from(
            alleles
                .into_iter()
                .map(|allele| String::from_utf8(allele).unwrap())
                .collect::<Vec<_>>(),
        );
        Ok(Outcome::Changed)
    }
}

fn is_sequence(allele: &[u8]) -> bool {
    !allele.is_empty()
        && allele
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn invalid(record: &RecordBuf, message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "normalizing {}:{}: {message}",
        record.reference_sequence_name(),
        record.variant_start().map_or(0, usize::from)
    ))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use noodles_vcf as vcf;

    use super::*;

    fn record(header: &vcf::Header, line: &[u8]) -> RecordBuf {
        let raw = vcf::Record::try_from(line).unwrap();
        RecordBuf::try_from_variant_record(header, &raw).unwrap()
    }

    #[test]
    fn left_aligns_indels_and_trims_substitutions_with_shared_reference_access() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        fs::write(&path, b">chr1\nAAAAAACGTACGT\n").unwrap();
        fs::write(path.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let mut normalizer = ReferenceNormalizer::open(&path).unwrap();

        for (line, position, reference, alternate) in [
            (b"chr1\t4\t.\tA\tAA\t.\tPASS\t.".as_slice(), 1, "A", "AA"),
            (b"chr1\t4\t.\tAA\tA\t.\tPASS\t.".as_slice(), 1, "AA", "A"),
            (b"chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.".as_slice(), 11, "C", "G"),
        ] {
            let mut record = record(&header, line);
            assert_eq!(normalizer.normalize(&mut record).unwrap(), Outcome::Changed);
            assert_eq!(record.variant_start().map(usize::from), Some(position));
            assert_eq!(record.reference_bases(), reference);
            assert_eq!(record.alternate_bases().as_ref(), [alternate]);
        }
    }

    #[test]
    fn reference_mismatch_and_invalid_range_fail_with_record_context() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reference.fa");
        fs::write(&path, b">chr1\nACGT\n").unwrap();
        fs::write(path.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let mut normalizer = ReferenceNormalizer::open(&path).unwrap();

        let mut mismatch = record(&header, b"chr1\t2\t.\tT\tA\t.\tPASS\t.");
        let error = normalizer.normalize(&mut mismatch).unwrap_err().to_string();
        assert!(error.contains("chr1:2"), "{error}");
        assert!(error.contains("does not match"), "{error}");

        let mut outside = record(&header, b"chr1\t4\t.\tTT\tT\t.\tPASS\t.");
        let error = normalizer.normalize(&mut outside).unwrap_err().to_string();
        assert!(error.contains("chr1"), "{error}");
        assert!(error.contains("length 4"), "{error}");
    }
}
