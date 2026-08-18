use std::io::{self, BufRead, Cursor, Read, Write};

use noodles_bgzf as bgzf;
use rsomics_common::{Result, RsomicsError};

use crate::format::bgzf::{Frame, FrameReader};

use super::header::HeaderText;
use super::{Encoding, Options, Summary, edit_header};

pub(super) fn rewrite_plain<R: BufRead, W: Write>(
    mut input: R,
    mut output: W,
    options: &Options,
) -> Result<Summary> {
    let mut raw_header = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        let read = input
            .read_until(b'\n', &mut line)
            .map_err(RsomicsError::Io)?;
        if read == 0 {
            break;
        }
        let chrom = line.starts_with(b"#CHROM\t") || line == b"#CHROM\n";
        raw_header.extend_from_slice(&line);
        if chrom {
            break;
        }
    }
    let original = HeaderText::parse(&raw_header)?;
    let (edited, summary) = edit_header(original, options, Encoding::PlainVcf)?;
    output
        .write_all(&edited.render())
        .map_err(RsomicsError::Io)?;
    std::io::copy(&mut input, &mut output).map_err(RsomicsError::Io)?;
    output.flush().map_err(RsomicsError::Io)?;
    Ok(summary)
}

pub(super) fn rewrite_bgzf<R: Read, W: Write>(
    input: R,
    output: W,
    options: &Options,
) -> Result<Summary> {
    let mut frames = FrameReader::new(input);
    let mut inflated = Vec::new();
    let mut scanned = 0;
    let (header_end, eof_consumed) = loop {
        match frames.next().map_err(map_frame_error)? {
            Some(Frame::Data(raw)) => {
                inflated.extend_from_slice(&inflate_frame(&raw)?);
                if inflated.starts_with(b"BCF") {
                    return Err(invalid("BGZF BCF reheader is not available"));
                }
                if let Some(end) = scan_header(&inflated, &mut scanned)? {
                    break (end, false);
                }
            }
            Some(Frame::Eof) => {
                if is_chrom_line(&inflated[scanned..]) {
                    break (inflated.len(), true);
                }
                return Err(invalid("BGZF VCF header is missing the #CHROM line"));
            }
            None => return Err(invalid("canonical BGZF EOF block is missing")),
        }
    };

    let original = HeaderText::parse(&inflated[..header_end])?;
    let (edited, summary) = edit_header(original, options, Encoding::BgzfVcf)?;
    let mut writer = bgzf::io::Writer::new(output);
    writer
        .write_all(&edited.render())
        .and_then(|()| writer.write_all(&inflated[header_end..]))
        .map_err(RsomicsError::Io)?;
    if eof_consumed {
        let mut output = writer.finish().map_err(RsomicsError::Io)?;
        output.flush().map_err(RsomicsError::Io)?;
    } else {
        writer.flush().map_err(RsomicsError::Io)?;
        let mut output = writer.into_inner();
        frames
            .copy_through_eof(&mut output)
            .map_err(map_frame_error)?;
        output.flush().map_err(RsomicsError::Io)?;
    }
    Ok(summary)
}

pub(super) fn inflate_frame(raw: &[u8]) -> Result<Vec<u8>> {
    let mut reader = bgzf::io::Reader::new(Cursor::new(raw));
    let mut inflated = Vec::new();
    reader
        .read_to_end(&mut inflated)
        .map_err(|error| invalid(format!("inflating BGZF header frame: {error}")))?;
    Ok(inflated)
}

fn scan_header(source: &[u8], scanned: &mut usize) -> Result<Option<usize>> {
    while let Some(relative_end) = source[*scanned..].iter().position(|byte| *byte == b'\n') {
        let end = *scanned + relative_end + 1;
        let mut line = &source[*scanned..end - 1];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        if is_chrom_line(line) {
            return Ok(Some(end));
        }
        if !line.starts_with(b"##") {
            return Err(invalid("BGZF VCF contains a non-header line before #CHROM"));
        }
        *scanned = end;
    }
    Ok(None)
}

fn is_chrom_line(line: &[u8]) -> bool {
    line.starts_with(b"#CHROM\t")
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
