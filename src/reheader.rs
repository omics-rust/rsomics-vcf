use std::fs::File;
use std::io::{self, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use serde::Serialize;

use crate::format::bgzf::{Frame, FrameReader};

mod bcf;
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
    BgzfBcf,
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
    let source = open(input)?;
    let (prefix, source) = prefix(source)?;
    let encoding = detect(&prefix)?;
    let reader = BufReader::new(Cursor::new(prefix).chain(source));
    match encoding {
        Encoding::PlainVcf => vcf::rewrite_plain(reader, output, options),
        Encoding::BgzfVcf => route_bgzf(reader, output, options),
        Encoding::RawBcf => bcf::rewrite_raw(reader, output, options),
        Encoding::BgzfBcf => unreachable!("BGZF input is classified after inflation"),
    }
}

pub(crate) fn write_parallel<W>(
    input: &Path,
    options: &Options,
    output: W,
    workers: usize,
) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    let source = open(input)?;
    let (prefix, source) = prefix(source)?;
    let encoding = detect(&prefix)?;
    if encoding != Encoding::BgzfVcf {
        return Err(RsomicsError::ConfigError(
            "--threads is available only for BGZF BCF input".to_owned(),
        ));
    }
    let reader = BufReader::new(Cursor::new(prefix).chain(source));
    let (is_bcf, reader) = inspect_bgzf(reader)?;
    if !is_bcf {
        return Err(RsomicsError::ConfigError(
            "--threads is available only for BGZF BCF input".to_owned(),
        ));
    }
    bcf::rewrite_bgzf_parallel(reader, output, options, workers)
}

fn open(input: &Path) -> Result<Box<dyn Read>> {
    if input == Path::new("-") {
        Ok(Box::new(io::stdin()))
    } else {
        File::open(input)
            .map(|file| Box::new(file) as Box<dyn Read>)
            .rs_with_context(|| format!("opening variant input {}", input.display()))
    }
}

fn route_bgzf<R: Read, W: Write>(input: R, output: W, options: &Options) -> Result<Summary> {
    let (is_bcf, input) = inspect_bgzf(input)?;
    if is_bcf {
        bcf::rewrite_bgzf(input, output, options)
    } else {
        vcf::rewrite_bgzf(input, output, options)
    }
}

fn inspect_bgzf<R: Read>(input: R) -> Result<(bool, impl Read)> {
    let mut frames = FrameReader::new(input);
    let mut raw_prefix = Vec::new();
    loop {
        match frames.next().map_err(map_frame_error)? {
            Some(Frame::Data(raw)) => {
                let inflated = vcf::inflate_frame(&raw)?;
                raw_prefix.extend_from_slice(&raw);
                if !inflated.is_empty() {
                    return Ok((
                        inflated.starts_with(b"BCF"),
                        Cursor::new(raw_prefix).chain(frames.into_inner()),
                    ));
                }
            }
            Some(Frame::Eof) => return Err(invalid("BGZF stream contains no variant data")),
            None => return Err(invalid("canonical BGZF EOF block is missing")),
        }
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

fn map_frame_error(error: io::Error) -> RsomicsError {
    if error.kind() == io::ErrorKind::InvalidData {
        invalid(format!("reading BGZF stream: {error}"))
    } else {
        RsomicsError::Io(error)
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
