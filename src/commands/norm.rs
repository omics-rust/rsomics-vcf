use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Result, RsomicsError, write_atomic};

use crate::cli::CommandOutput;
use crate::commands::variant::OutputType;
use crate::norm::{self, Options};

#[derive(Debug, Args)]
#[command(after_help = "\
Output types:
  v  uncompressed VCF
  z  BGZF-compressed VCF
  b  BGZF-compressed BCF
  u  uncompressed BCF

Example:
  rsomics-vcf norm -f reference.fa calls.vcf.gz -O z -o normalized.vcf.gz")]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Indexed FASTA reference
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE")]
    reference: PathBuf,

    /// Write normalized variants to this file instead of standard output
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        default_value = "-",
        hide_default_value = true
    )]
    output: PathBuf,

    /// Output encoding: v, z, b, or u
    #[arg(short = 'O', long, value_name = "TYPE", default_value = "v")]
    output_type: OutputType,

    /// Coordinate window retained for records moved by realignment
    #[arg(short = 'w', long, value_name = "INT", default_value_t = 1000)]
    site_window: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because variant data otherwise uses standard output"
                .to_owned(),
        ));
    }
    let options = Options {
        reference: arguments.reference,
        output_format: arguments.output_type.into(),
        site_window: arguments.site_window,
    };
    let summary = if arguments.output == Path::new("-") {
        norm::write(&arguments.input, &options, io::stdout().lock())?
    } else {
        write_atomic(&arguments.output, |output| {
            norm::write(&arguments.input, &options, output)
        })?
    };
    Ok(CommandOutput::Norm { summary })
}
