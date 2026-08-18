use std::collections::HashMap;
use std::io::{self, Read, Write};

use noodles_vcf::{self as vcf, header::StringMaps, variant::RecordBuf};
use rsomics_common::{Context, Result, RsomicsError};

use crate::format::bgzf::ValidatedReader;
use crate::format::{HeaderMode, OutputFormat, ParallelWriter, VariantWriter, Writer};

use super::header::HeaderText;
use super::{Encoding, Options, Summary, edit_header};

pub(super) fn rewrite_raw<R: Read, W: Write>(
    input: R,
    output: W,
    options: &Options,
) -> Result<Summary> {
    let reader = noodles_bcf::io::Reader::from(input);
    let mut writer = Writer::new(output, OutputFormat::BcfRaw);
    rewrite(reader, &mut writer, options, Encoding::RawBcf)
}

pub(super) fn rewrite_bgzf<R: Read, W: Write>(
    input: R,
    output: W,
    options: &Options,
) -> Result<Summary> {
    let reader = noodles_bcf::io::Reader::new(ValidatedReader::new(input));
    let mut writer = Writer::new(output, OutputFormat::Bcf);
    rewrite(reader, &mut writer, options, Encoding::BgzfBcf)
}

pub(super) fn rewrite_bgzf_parallel<R, W>(
    input: R,
    output: W,
    options: &Options,
    workers: usize,
) -> Result<Summary>
where
    R: Read,
    W: Write + Send + 'static,
{
    let reader = noodles_bcf::io::Reader::new(ValidatedReader::new(input));
    let mut writer = ParallelWriter::new(output, OutputFormat::Bcf, workers)?;
    rewrite(reader, &mut writer, options, Encoding::BgzfBcf)
}

fn rewrite<R: Read>(
    mut reader: noodles_bcf::io::Reader<R>,
    writer: &mut impl VariantWriter,
    options: &Options,
    encoding: Encoding,
) -> Result<Summary> {
    let (original_text, original) = read_header(&mut reader)?;
    let (edited_text, summary) = edit_header(original_text, options, encoding)?;
    let mut edited = edited_text.parse_noodles()?;
    preserve_indices(&original, &mut edited)?;
    writer.write_header(&edited, HeaderMode::Full)?;

    let mut raw = noodles_bcf::Record::default();
    let mut number = 0u64;
    loop {
        let next = number
            .checked_add(1)
            .ok_or_else(|| invalid("BCF record count exceeds u64"))?;
        let read = reader
            .read_record(&mut raw)
            .map_err(|error| input_error(format!("reading BCF record {next}"), error))?;
        if read == 0 {
            break;
        }
        number = next;
        let record = RecordBuf::try_from_variant_record(&original, &raw)
            .map_err(|error| invalid(format!("decoding BCF record {number}: {error}")))?;
        writer
            .write_record(&edited, &record, number)
            .rs_with_context(|| format!("reheadering BCF record {number}"))?;
    }
    writer.finish()?;
    Ok(summary)
}

fn read_header<R: Read>(
    reader: &mut noodles_bcf::io::Reader<R>,
) -> Result<(HeaderText, vcf::Header)> {
    let raw = {
        let mut reader = reader.header_reader();
        let magic = reader
            .read_magic_number()
            .map_err(|error| input_error("reading BCF magic number", error))?;
        if magic != *b"BCF" {
            return Err(invalid("invalid BCF magic number"));
        }
        let version = reader
            .read_format_version()
            .map_err(|error| input_error("reading BCF version", error))?;
        if version != (2, 2) {
            return Err(invalid(format!(
                "unsupported BCF version {}.{}",
                version.0, version.1
            )));
        }
        let mut raw_reader = reader
            .raw_vcf_header_reader()
            .map_err(|error| input_error("reading BCF header length", error))?;
        let mut raw = Vec::new();
        raw_reader
            .read_to_end(&mut raw)
            .map_err(|error| input_error("reading BCF header", error))?;
        raw_reader
            .discard_to_end()
            .map_err(|error| input_error("discarding BCF header padding", error))?;
        raw
    };
    let text = HeaderText::parse(&raw)?;
    let mut header = text.parse_noodles()?;
    *header.string_maps_mut() = StringMaps::try_from(&header)
        .map_err(|error| invalid(format!("reading BCF dictionaries: {error}")))?;
    Ok((text, header))
}

fn preserve_indices(original: &vcf::Header, edited: &mut vcf::Header) -> Result<()> {
    let mut next_contig = original
        .contigs()
        .keys()
        .filter_map(|name| original.string_maps().contigs().get_index_of(name))
        .max()
        .map_or(0, |index| index + 1);
    for (name, contig) in edited.contigs_mut() {
        let index = original
            .string_maps()
            .contigs()
            .get_index_of(name)
            .unwrap_or_else(|| {
                let index = next_contig;
                next_contig += 1;
                index
            });
        *contig.idx_mut() = Some(index);
    }

    let mut next_string = original_string_indices(original)
        .values()
        .copied()
        .max()
        .map_or(1, |index| index + 1);
    let mut assigned = HashMap::new();
    for (name, info) in edited.infos_mut() {
        let index = assign_string_index(original, &mut assigned, &mut next_string, name);
        *info.idx_mut() = Some(index);
    }
    for (name, filter) in edited.filters_mut() {
        let index = assign_string_index(original, &mut assigned, &mut next_string, name);
        *filter.idx_mut() = Some(index);
    }
    for (name, format) in edited.formats_mut() {
        let index = assign_string_index(original, &mut assigned, &mut next_string, name);
        *format.idx_mut() = Some(index);
    }

    let maps = StringMaps::try_from(&*edited)
        .map_err(|error| invalid(format!("building edited BCF dictionaries: {error}")))?;
    for (name, expected) in assigned {
        if maps.strings().get_index_of(&name) != Some(expected) {
            return Err(invalid(format!(
                "edited BCF string dictionary moved {name} from index {expected}"
            )));
        }
    }
    for name in edited.contigs().keys() {
        if let Some(expected) = original.string_maps().contigs().get_index_of(name)
            && maps.contigs().get_index_of(name) != Some(expected)
        {
            return Err(invalid(format!(
                "edited BCF contig dictionary moved {name} from index {expected}"
            )));
        }
    }
    *edited.string_maps_mut() = maps;
    Ok(())
}

fn original_string_indices(header: &vcf::Header) -> HashMap<String, usize> {
    let mut indices = HashMap::from([(String::from("PASS"), 0)]);
    for name in header
        .infos()
        .keys()
        .chain(header.filters().keys())
        .chain(header.formats().keys())
    {
        if let Some(index) = header.string_maps().strings().get_index_of(name) {
            indices.insert(name.clone(), index);
        }
    }
    indices
}

fn assign_string_index(
    original: &vcf::Header,
    assigned: &mut HashMap<String, usize>,
    next: &mut usize,
    name: &str,
) -> usize {
    if let Some(index) = assigned.get(name) {
        return *index;
    }
    let index = original
        .string_maps()
        .strings()
        .get_index_of(name)
        .unwrap_or_else(|| {
            let index = *next;
            *next += 1;
            index
        });
    assigned.insert(name.to_owned(), index);
    index
}

fn input_error(context: impl std::fmt::Display, error: io::Error) -> RsomicsError {
    let message = format!("{context}: {error}");
    if matches!(
        error.kind(),
        io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof
    ) {
        invalid(message)
    } else {
        RsomicsError::Io(io::Error::new(error.kind(), message))
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, header::StringMaps};

    use super::preserve_indices;

    fn header(source: &str) -> vcf::Header {
        let mut header: vcf::Header = source.parse().unwrap();
        *header.string_maps_mut() = StringMaps::try_from(&header).unwrap();
        header
    }

    #[test]
    fn retained_indices_stay_fixed_and_new_ids_follow_the_original_maps() {
        let original = header(
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100,IDX=2>\n\
##FILTER=<ID=q10,Description=\"low\",IDX=4>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\",IDX=2>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\",IDX=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
        );
        let mut edited = header(
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=200>\n\
##contig=<ID=chr2,length=300>\n\
##FILTER=<ID=q10,Description=\"low\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=XX,Number=1,Type=Integer,Description=\"new\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN1\n",
        );
        preserve_indices(&original, &mut edited).unwrap();

        assert_eq!(edited.contigs()["chr1"].idx(), Some(2));
        assert_eq!(edited.contigs()["chr2"].idx(), Some(3));
        assert_eq!(edited.formats()["GT"].idx(), Some(1));
        assert_eq!(edited.infos()["DP"].idx(), Some(2));
        assert_eq!(edited.filters()["q10"].idx(), Some(4));
        assert_eq!(edited.infos()["XX"].idx(), Some(5));
        assert_eq!(edited.string_maps().contigs().get_index(2), Some("chr1"));
        assert_eq!(edited.string_maps().contigs().get_index(3), Some("chr2"));
        assert_eq!(edited.string_maps().strings().get_index(1), Some("GT"));
        assert_eq!(edited.string_maps().strings().get_index(2), Some("DP"));
        assert_eq!(edited.string_maps().strings().get_index(4), Some("q10"));
        assert_eq!(edited.string_maps().strings().get_index(5), Some("XX"));
    }
}
