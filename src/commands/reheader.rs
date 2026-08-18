use std::io;
use std::path::{Path, PathBuf};

use clap::{ArgGroup, Args};
use rsomics_common::{AtomicFile, Result, RsomicsError, reject_output_alias};

use crate::cli::CommandOutput;
use crate::reheader::{self, Options, SampleSource};

#[derive(Debug, Args)]
#[command(group(
    ArgGroup::new("edit")
        .required(true)
        .multiple(true)
        .args(["header", "fai", "samples_list", "samples_file"])
))]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Replace the complete VCF header from FILE
    #[arg(short = 'H', long, value_name = "FILE")]
    header: Option<PathBuf>,

    /// Synchronize contig names and lengths from a FASTA index
    #[arg(short = 'f', long, value_name = "FILE")]
    fai: Option<PathBuf>,

    /// Replace samples positionally from a comma-separated list
    #[arg(
        short = 'n',
        long,
        value_name = "LIST",
        conflicts_with = "samples_file"
    )]
    samples_list: Option<String>,

    /// Read positional names or old-to-new sample pairs from FILE
    #[arg(short = 'N', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,

    /// Write variants to FILE instead of standard output
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        default_value = "-",
        hide_default_value = true
    )]
    output: PathBuf,

    /// Use INT BGZF compression workers for BGZF BCF
    #[arg(long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because variant data otherwise uses standard output"
                .to_owned(),
        ));
    }
    let samples = arguments
        .samples_list
        .map(SampleSource::List)
        .or_else(|| arguments.samples_file.clone().map(SampleSource::File));
    let mut inputs = vec![arguments.input.as_path()];
    inputs.extend(arguments.header.as_deref());
    inputs.extend(arguments.fai.as_deref());
    inputs.extend(arguments.samples_file.as_deref());
    reject_output_alias(&arguments.output, inputs)?;

    let options = Options {
        header: arguments.header,
        fai: arguments.fai,
        samples,
    };
    let summary = if arguments.output == Path::new("-") {
        if arguments.threads == 0 {
            reheader::write(&arguments.input, &options, io::stdout().lock())?
        } else {
            reheader::write_parallel(&arguments.input, &options, io::stdout(), arguments.threads)?
        }
    } else {
        let mut transaction = AtomicFile::new(&arguments.output)?;
        let summary = if arguments.threads == 0 {
            reheader::write(&arguments.input, &options, transaction.file_mut())?
        } else {
            reheader::write_parallel(
                &arguments.input,
                &options,
                transaction.reopen()?,
                arguments.threads,
            )?
        };
        transaction.commit()?;
        summary
    };
    Ok(CommandOutput::Reheader { summary })
}
