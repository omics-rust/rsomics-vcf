mod fields;

use noodles_core::Position;
use noodles_vcf::{
    Header,
    header::record::value::{
        Map,
        map::{Info, info},
    },
    variant::{
        RecordBuf,
        record_buf::{AlternateBases, info::field::Value as InfoValue},
    },
};
use rsomics_common::{Result, RsomicsError};

use fields::{extend_allele_fields, remap_genotypes};

#[derive(Debug, Eq, Ord, PartialEq, PartialOrd)]
struct Atom {
    begin: usize,
    reference: String,
    alternate: String,
    end: usize,
    prefix_len: usize,
    source_alternate: usize,
}

struct OutputAtom {
    atom: Atom,
    mapping: Vec<usize>,
    star: bool,
}

pub(super) struct Origin<'a> {
    pub(super) record: &'a RecordBuf,
    pub(super) alternate: usize,
}

pub(super) fn prepare_header(header: &mut Header, tag: &str) -> Result<()> {
    if !valid_tag(tag) {
        return Err(RsomicsError::ConfigError(format!(
            "invalid --old-rec-tag ID: {tag}"
        )));
    }
    if let Some(schema) = header.infos().get(tag) {
        if schema.number() != info::Number::Count(1) || schema.ty() != info::Type::String {
            return Err(RsomicsError::ConfigError(format!(
                "--old-rec-tag {tag} requires INFO/{tag} Number=1,Type=String"
            )));
        }
    } else {
        header.infos_mut().insert(
            tag.to_owned(),
            Map::<Info>::new(
                info::Number::Count(1),
                info::Type::String,
                "Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX",
            ),
        );
    }
    Ok(())
}

pub(super) fn atomize(
    header: &Header,
    mut record: RecordBuf,
    star_allele: bool,
    old_record_tag: Option<&str>,
    origin: Option<Origin<'_>>,
) -> Result<(Vec<RecordBuf>, bool)> {
    let alternates = record.alternate_bases().as_ref();
    let reference = record.reference_bases().as_bytes();
    if !is_sequence(reference)
        || alternates
            .iter()
            .any(|alternate| !is_sequence(alternate.as_bytes()))
    {
        return Ok((vec![record], false));
    }
    let start = record.variant_start().map(usize::from).ok_or_else(|| {
        RsomicsError::InvalidInput("atomizing a record without a position".to_owned())
    })?;
    let mut atoms = Vec::new();
    for (alternate, sequence) in alternates.iter().enumerate() {
        atoms.extend(decompose(reference, sequence.as_bytes(), alternate + 1)?);
    }
    if atoms.len() == 1
        && alternates.len() == 1
        && atoms[0].begin == 0
        && atoms[0]
            .reference
            .as_bytes()
            .eq_ignore_ascii_case(reference)
        && atoms[0]
            .alternate
            .as_bytes()
            .eq_ignore_ascii_case(alternates[0].as_bytes())
    {
        if let Some(tag) = old_record_tag.filter(|tag| record.info().get(*tag).is_none()) {
            let source = origin.as_ref().map_or(&record, |origin| origin.record);
            let alternate = origin.as_ref().map(|origin| origin.alternate);
            let value = origin_value(source, &[0, 1], alternate)?;
            insert_origin(&mut record, tag, value);
        }
        return Ok((vec![record], false));
    }
    atoms.sort();
    let atoms = output_atoms(atoms, alternates.len());
    let mut records = Vec::new();
    for output_atom in atoms {
        let position = start.checked_add(output_atom.atom.begin).ok_or_else(|| {
            RsomicsError::InvalidInput("atomized position exceeds usize".to_owned())
        })?;
        let mut output = super::split::split_one(
            header,
            &record,
            output_atom.atom.source_alternate - 1,
            false,
        )?;
        *output.variant_start_mut() = Some(Position::try_from(position).map_err(|_| {
            RsomicsError::InvalidInput(
                "atomized position exceeds the VCF coordinate range".to_owned(),
            )
        })?);
        *output.reference_bases_mut() = output_atom.atom.reference;
        let mut alternates = vec![output_atom.atom.alternate];
        let has_star = output_atom.star && star_allele;
        if has_star {
            alternates.push("*".to_owned());
        }
        *output.alternate_bases_mut() = AlternateBases::from(alternates);
        remap_genotypes(&record, &mut output, &output_atom.mapping, star_allele)?;
        extend_allele_fields(header, &mut output, has_star)?;
        if let Some(tag) = old_record_tag {
            let source = origin.as_ref().map_or(&record, |origin| origin.record);
            let alternate = origin.as_ref().map(|origin| origin.alternate);
            annotate_origin(source, &mut output, tag, &output_atom.mapping, alternate)?;
        }
        records.push(output);
    }
    Ok((records, true))
}

fn annotate_origin(
    source: &RecordBuf,
    output: &mut RecordBuf,
    tag: &str,
    mapping: &[usize],
    source_alternate: Option<usize>,
) -> Result<()> {
    if output.info().get(tag).is_some() {
        return Ok(());
    }
    let value = origin_value(source, mapping, source_alternate)?;
    insert_origin(output, tag, value);
    Ok(())
}

fn origin_value(
    source: &RecordBuf,
    mapping: &[usize],
    source_alternate: Option<usize>,
) -> Result<String> {
    let position = source
        .variant_start()
        .map(usize::from)
        .ok_or_else(|| invalid("atomizing a record without a position"))?;
    let alternates = source.alternate_bases().as_ref().join(",");
    let used = source_alternate.map_or_else(
        || {
            mapping
                .iter()
                .enumerate()
                .filter(|(_, value)| **value == 1)
                .map(|(index, _)| (index + 1).to_string())
                .collect::<Vec<_>>()
                .join(",")
        },
        |alternate| alternate.to_string(),
    );
    Ok(format!(
        "{}|{position}|{}|{alternates}|{used}",
        source.reference_sequence_name(),
        source.reference_bases()
    ))
}

fn insert_origin(output: &mut RecordBuf, tag: &str, value: String) {
    output
        .info_mut()
        .insert(tag.to_owned(), Some(InfoValue::String(value)));
}

fn valid_tag(tag: &str) -> bool {
    if tag == "1000G" {
        return true;
    }
    let mut bytes = tag.bytes();
    bytes
        .next()
        .is_some_and(|byte| byte.is_ascii_alphabetic() || byte == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.'))
}

fn decompose(reference: &[u8], alternate: &[u8], source_alternate: usize) -> Result<Vec<Atom>> {
    let mut reference_len = reference.len();
    let mut alternate_len = alternate.len();
    while reference_len > 1
        && alternate_len > 1
        && reference[reference_len - 1].eq_ignore_ascii_case(&alternate[alternate_len - 1])
    {
        reference_len -= 1;
        alternate_len -= 1;
    }

    let mut atoms: Vec<Atom> = Vec::new();
    let mut current: Option<usize> = None;
    for offset in 0..reference_len.max(alternate_len) {
        let reference_base =
            (offset < reference_len).then(|| reference[offset].to_ascii_uppercase());
        let alternate_base =
            (offset < alternate_len).then(|| alternate[offset].to_ascii_uppercase());
        if reference_base != alternate_base {
            if reference_base.is_none() || alternate_base.is_none() {
                let atom = current
                    .and_then(|index| atoms.get_mut(index))
                    .ok_or_else(|| invalid("complex allele has no anchoring base"))?;
                if let Some(base) = alternate_base {
                    atom.alternate.push(char::from(base));
                }
                if let Some(base) = reference_base {
                    atom.reference.push(char::from(base));
                    atom.end += 1;
                }
                continue;
            }

            atoms.push(Atom {
                begin: offset,
                reference: char::from(reference_base.unwrap()).to_string(),
                alternate: char::from(alternate_base.unwrap()).to_string(),
                end: offset,
                prefix_len: 0,
                source_alternate,
            });
            current = Some(atoms.len() - 1);
            if reference_len != alternate_len
                && (offset + 1 >= reference_len || offset + 1 >= alternate_len)
            {
                let base = reference_base.unwrap();
                atoms.push(Atom {
                    begin: offset,
                    reference: char::from(base).to_string(),
                    alternate: char::from(base).to_string(),
                    end: offset,
                    prefix_len: 1,
                    source_alternate,
                });
                current = Some(atoms.len() - 1);
            }
        } else if offset + 1 >= reference_len || offset + 1 >= alternate_len {
            let base = reference_base.ok_or_else(|| invalid("complex allele is empty"))?;
            atoms.push(Atom {
                begin: offset,
                reference: char::from(base).to_string(),
                alternate: char::from(base).to_string(),
                end: offset,
                prefix_len: 0,
                source_alternate,
            });
            current = Some(atoms.len() - 1);
        }
    }
    Ok(atoms)
}

fn overlaps(left: &Atom, right: &Atom) -> bool {
    left.begin <= right.end && right.begin <= left.end
}

fn same_variant(left: &Atom, right: &Atom) -> bool {
    left.begin == right.begin
        && left.reference.eq_ignore_ascii_case(&right.reference)
        && left.alternate.eq_ignore_ascii_case(&right.alternate)
}

fn output_atoms(atoms: Vec<Atom>, alternate_count: usize) -> Vec<OutputAtom> {
    let mut outputs = Vec::new();
    for atom in &atoms {
        if outputs
            .last()
            .is_some_and(|output: &OutputAtom| same_variant(&output.atom, atom))
        {
            continue;
        }
        let mut mapping = vec![0; alternate_count];
        let mut star = false;
        for candidate in &atoms {
            if same_variant(atom, candidate) {
                mapping[candidate.source_alternate - 1] = 1;
            } else if overlaps(atom, candidate) {
                mapping[candidate.source_alternate - 1] = overlap_mapping(atom, candidate);
                star = true;
            }
        }
        outputs.push(OutputAtom {
            atom: Atom {
                begin: atom.begin,
                reference: atom.reference.clone(),
                alternate: atom.alternate.clone(),
                end: atom.end,
                prefix_len: atom.prefix_len,
                source_alternate: atom.source_alternate,
            },
            mapping,
            star,
        });
    }
    outputs
}

fn overlap_mapping(left: &Atom, right: &Atom) -> usize {
    if left.begin != right.begin {
        return 2;
    }
    if (left.prefix_len > 0 && left.prefix_len >= right.reference.len())
        || (right.prefix_len > 0 && right.prefix_len >= left.reference.len())
    {
        return 1;
    }
    if !left.reference.eq_ignore_ascii_case(&right.reference) {
        return 2;
    }
    if (left.prefix_len > 0 && left.prefix_len >= right.alternate.len())
        || (right.prefix_len > 0 && right.prefix_len >= left.alternate.len())
    {
        return 1;
    }
    usize::from(!left.alternate.eq_ignore_ascii_case(&right.alternate)) * 2
}

fn is_sequence(allele: &[u8]) -> bool {
    allele
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
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
        let (records, changed) = atomize(&header, record, true, None, None).unwrap();
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

    #[test]
    fn decomposes_complex_indels_and_marks_overlapping_atoms() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap();
        let raw =
            vcf::Record::try_from(b"chr1\t20\tb2\tAC\tGTG\t50\tPASS\tDP=7\tGT\t1/1".as_slice())
                .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let (records, changed) = atomize(&header, record, true, None, None).unwrap();
        assert!(changed);

        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        assert_eq!(
            String::from_utf8(writer.into_inner()).unwrap(),
            "chr1\t20\tb2\tA\tG\t50\tPASS\tDP=7\tGT\t1/1\n\
chr1\t21\tb2\tC\tCG,*\t50\tPASS\tDP=7\tGT\t1/1\n\
chr1\t21\tb2\tC\tT,*\t50\tPASS\tDP=7\tGT\t1/1\n"
        );
    }

    #[test]
    fn extends_allele_indexed_fields_for_spanning_deletions() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(
            b"chr1\t20\tb2\tAC\tGTG\t50\tPASS\tIA=8;IR=10,5\tGT:AD:PL\t1/1:10,5:0,10,20\t1:20,7:5,15"
                .as_slice(),
        )
        .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let (records, changed) = atomize(&header, record, true, None, None).unwrap();
        assert!(changed);

        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        assert_eq!(
            String::from_utf8(writer.into_inner()).unwrap(),
            "chr1\t20\tb2\tA\tG\t50\tPASS\tIA=8;IR=10,5\tGT:AD:PL\t1/1:10,5:0,10,20\t1:20,7:5,15\n\
chr1\t21\tb2\tC\tCG,*\t50\tPASS\tIA=8,.;IR=10,5,.\tGT:AD:PL\t1/1:10,5,.:0,10,20,.,.,.\t1:20,7,.:5,15,.\n\
chr1\t21\tb2\tC\tT,*\t50\tPASS\tIA=8,.;IR=10,5,.\tGT:AD:PL\t1/1:10,5,.:0,10,20,.,.,.\t1:20,7,.:5,15,.\n"
        );
    }

    #[test]
    fn remaps_multiallelic_genotypes_across_overlapping_atoms() {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(
            b"chr1\t50\tm1\tCC\tC,GG\t50\tPASS\tIA=3,8;IR=10,4,6\tGT:AD\t1/2:10,4,6\t0|2:20,1,7\t1:30,2,8"
                .as_slice(),
        )
        .unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        let (records, changed) = atomize(&header, record, true, None, None).unwrap();
        assert!(changed);

        let mut writer = vcf::io::Writer::new(Vec::new());
        for record in &records {
            writer.write_variant_record(&header, record).unwrap();
        }
        assert_eq!(
            String::from_utf8(writer.into_inner()).unwrap(),
            "chr1\t50\tm1\tC\tG,*\t50\tPASS\tIA=8,.;IR=10,6,.\tGT:AD\t2/1:10,6,.\t0|1:20,7,.\t2:30,8,.\n\
chr1\t50\tm1\tCC\tC,*\t50\tPASS\tIA=3,.;IR=10,4,.\tGT:AD\t1/2:10,4,.\t0|2:20,1,.\t1:30,2,.\n\
chr1\t51\tm1\tC\tG,*\t50\tPASS\tIA=8,.;IR=10,6,.\tGT:AD\t2/1:10,6,.\t0|1:20,7,.\t2:30,8,.\n"
        );
    }
}
