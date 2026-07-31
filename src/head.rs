use std::io::{BufWriter, Write};
use std::path::Path;

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::format::{HeaderTypes, Reader, reformat_record, trim_line_ending};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Options {
    pub header_lines: Option<usize>,
    pub records: usize,
    pub from_chrom: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub header_lines: usize,
    pub records: usize,
}

pub fn write(input: &Path, options: Options, mut output: impl Write) -> Result<Summary> {
    let mut output = BufWriter::with_capacity(1 << 20, &mut output);
    let mut reader = Reader::open(input)?;
    let (header, header_bytes, types) = reader.read_header()?;
    let header_lines = write_header(
        &mut output,
        &header_bytes,
        options.header_lines,
        options.from_chrom,
    )?;

    let records = if options.records == 0 {
        0
    } else if reader.is_text() {
        write_text_records(&mut reader, &types, options.records, input, &mut output)?
    } else {
        write_bcf_records(
            &mut reader,
            &header,
            &types,
            options.records,
            input,
            &mut output,
        )?
    };

    output.flush().map_err(RsomicsError::Io)?;
    Ok(Summary {
        header_lines,
        records,
    })
}

fn write_text_records(
    reader: &mut Reader,
    header: &HeaderTypes,
    limit: usize,
    input: &Path,
    output: &mut impl Write,
) -> Result<usize> {
    let mut raw = Vec::with_capacity(4096);
    let mut rendered = Vec::with_capacity(4096);
    let mut records = 0;
    while records < limit {
        let number = records + 1;
        if reader.read_text_record(&mut raw, number)? == 0 {
            break;
        }
        reformat_record(&raw, header, &mut rendered).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "{}: parsing variant record {number}: {error}",
                input.display()
            ))
        })?;
        output.write_all(&rendered).map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
        records += 1;
    }
    Ok(records)
}

fn write_bcf_records(
    reader: &mut Reader,
    header: &vcf::Header,
    types: &HeaderTypes,
    limit: usize,
    input: &Path,
    output: &mut impl Write,
) -> Result<usize> {
    let mut record = bcf::Record::default();
    let mut raw = Vec::with_capacity(4096);
    let mut rendered = Vec::with_capacity(4096);
    let mut records = 0;
    while records < limit {
        let number = records + 1;
        if reader.read_bcf_record(&mut record, number)? == 0 {
            break;
        }
        raw.clear();
        vcf::io::Writer::new(&mut raw)
            .write_variant_record(header, &record)
            .map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "{}: writing variant record {number}: {}",
                    input.display(),
                    error_chain(&error)
                ))
            })?;
        trim_line_ending(&mut raw);
        reformat_record(&raw, types, &mut rendered).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "{}: parsing variant record {number}: {error}",
                input.display()
            ))
        })?;
        output.write_all(&rendered).map_err(RsomicsError::Io)?;
        output.write_all(b"\n").map_err(RsomicsError::Io)?;
        records += 1;
    }
    Ok(records)
}

fn error_chain(error: &(dyn std::error::Error + 'static)) -> String {
    let mut message = error.to_string();
    let mut source = error.source();
    while let Some(error) = source {
        let detail = error.to_string();
        if !detail.is_empty() && !message.ends_with(&detail) {
            message.push_str(": ");
            message.push_str(&detail);
        }
        source = error.source();
    }
    message
}

fn write_header(
    output: &mut impl Write,
    header: &[u8],
    limit: Option<usize>,
    include_chrom: bool,
) -> Result<usize> {
    let lines: Vec<_> = header.split_inclusive(|byte| *byte == b'\n').collect();
    let take = limit.unwrap_or(lines.len()).min(lines.len());
    let mut written = 0;
    let mut chrom_written = false;

    for line in &lines[..take] {
        output.write_all(line).map_err(RsomicsError::Io)?;
        chrom_written |= line.starts_with(b"#CHROM\t") || line.starts_with(b"#CHROM\n");
        written += 1;
    }

    if include_chrom && !chrom_written {
        let line = lines
            .iter()
            .find(|line| line.starts_with(b"#CHROM\t") || line.starts_with(b"#CHROM\n"))
            .ok_or_else(|| {
                RsomicsError::InvalidInput("VCF header is missing the #CHROM line".to_owned())
            })?;
        output.write_all(line).map_err(RsomicsError::Io)?;
        written += 1;
    }

    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &[u8] =
        b"##fileformat=VCFv4.3\n##source=test\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn header_limit_counts_complete_lines() {
        let mut output = Vec::new();
        assert_eq!(
            write_header(&mut output, HEADER, Some(1), false).unwrap(),
            1
        );
        assert_eq!(output, b"##fileformat=VCFv4.3\n");
    }

    #[test]
    fn samples_mode_adds_chrom_line_after_limited_metadata() {
        let mut output = Vec::new();
        assert_eq!(write_header(&mut output, HEADER, Some(1), true).unwrap(), 2);
        assert_eq!(
            output,
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        );
    }

    #[test]
    fn samples_zero_writes_only_chrom_line() {
        let mut output = Vec::new();
        assert_eq!(write_header(&mut output, HEADER, Some(0), true).unwrap(), 1);
        assert_eq!(output, b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    }
}
