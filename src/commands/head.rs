use std::io;
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::head;

#[derive(Debug, Args)]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Print at most this many header lines
    #[arg(short = 'H', long, value_name = "INT")]
    headers: Option<usize>,

    /// Print at most this many variant records
    #[arg(short = 'n', long, value_name = "INT", conflicts_with = "samples")]
    records: Option<usize>,

    /// Start at the #CHROM line and print this many records
    #[arg(short = 's', long, value_name = "INT")]
    samples: Option<usize>,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json {
        return Err(RsomicsError::ConfigError(
            "--json cannot be combined with VCF stream output".to_owned(),
        ));
    }

    let from_chrom = arguments.samples.is_some();
    let records = arguments.samples.or(arguments.records).unwrap_or(0);
    let header_lines = if from_chrom && arguments.headers.is_none() {
        Some(0)
    } else {
        arguments.headers
    };
    let summary = head::write(
        &arguments.input,
        head::Options {
            header_lines,
            records,
            from_chrom,
        },
        io::stdout().lock(),
    )?;

    Ok(CommandOutput::Head { summary })
}
