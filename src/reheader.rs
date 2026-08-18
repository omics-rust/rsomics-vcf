use std::fs::File;
use std::io::{self, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use serde::Serialize;

mod fai;
mod header;
mod samples;
mod vcf;

use header::HeaderText;
pub(crate) use samples::SampleSource;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum Encoding {
    PlainVcf,
    BgzfVcf,
    RawBcf,
}

pub(crate) struct Options {
    pub(crate) header: Option<PathBuf>,
    pub(crate) fai: Option<PathBuf>,
    pub(crate) samples: Option<SampleSource>,
}

#[derive(Debug, Serialize)]
pub(crate) struct Summary {
    encoding: Encoding,
    header_replaced: bool,
    fai_applied: bool,
    samples_renamed: bool,
    contigs_before: usize,
    contigs_after: usize,
    samples_before: usize,
    samples_after: usize,
}

pub(crate) fn write<W: Write>(input: &Path, options: &Options, output: W) -> Result<Summary> {
    let source: Box<dyn Read> = if input == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(input)
                .rs_with_context(|| format!("opening variant input {}", input.display()))?,
        )
    };
    let (prefix, source) = prefix(source)?;
    let encoding = detect(&prefix)?;
    let reader = BufReader::new(Cursor::new(prefix).chain(source));
    match encoding {
        Encoding::PlainVcf => vcf::rewrite_plain(reader, output, options),
        Encoding::BgzfVcf => Err(invalid("BGZF VCF reheader is not available")),
        Encoding::RawBcf => Err(invalid("BCF reheader is not available")),
    }
}

fn edit_header(
    original: HeaderText,
    options: &Options,
    encoding: Encoding,
) -> Result<(HeaderText, Summary)> {
    let contigs_before = original.contig_count();
    let samples_before = original.sample_names().len();
    let mut edited = if let Some(path) = &options.header {
        HeaderText::replace_from(path, samples_before)?
    } else {
        original
    };
    if let Some(path) = &options.fai {
        edited.apply_fai(path)?;
    }
    if let Some(source) = &options.samples {
        edited.apply_samples(source)?;
    }
    let summary = Summary {
        encoding,
        header_replaced: options.header.is_some(),
        fai_applied: options.fai.is_some(),
        samples_renamed: options.samples.is_some(),
        contigs_before,
        contigs_after: edited.contig_count(),
        samples_before,
        samples_after: edited.sample_names().len(),
    };
    Ok((edited, summary))
}

fn prefix(mut source: Box<dyn Read>) -> Result<(Vec<u8>, Box<dyn Read>)> {
    let mut buffer = [0; 18];
    let mut length = 0;
    while length < buffer.len() {
        match source.read(&mut buffer[length..]) {
            Ok(0) => break,
            Ok(read) => length += read,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(RsomicsError::Io(error)),
        }
    }
    Ok((buffer[..length].to_vec(), source))
}

fn detect(prefix: &[u8]) -> Result<Encoding> {
    if prefix.starts_with(b"BCF") {
        return Ok(Encoding::RawBcf);
    }
    if prefix.starts_with(&[0x1f, 0x8b]) {
        if prefix.len() >= 18
            && prefix[..4] == [0x1f, 0x8b, 0x08, 0x04]
            && prefix[10..16] == [0x06, 0x00, 0x42, 0x43, 0x02, 0x00]
        {
            return Ok(Encoding::BgzfVcf);
        }
        return Err(invalid(
            "ordinary gzip VCF is unsupported; convert it to plain VCF or BGZF VCF",
        ));
    }
    Ok(Encoding::PlainVcf)
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
