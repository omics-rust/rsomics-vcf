use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};

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
