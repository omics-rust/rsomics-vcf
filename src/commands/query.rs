use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{Context, Result, RsomicsError, write_atomic};

use crate::cli::CommandOutput;
use crate::query::{self, HeaderMode};

#[derive(Debug, Args)]
#[command(after_help = "\
Supported fields:
  %CHROM %POS %POS0 %END %END0 %ID %REF %ALT %FIRST_ALT
  %QUAL %FILTER %TYPE %INFO %INFO/TAG %FORMAT %LINE
  [ %SAMPLE %GT %TGT %IUPACGT %TAG ] and {0}-based subscripts

Examples:
  rsomics-vcf query calls.bcf -f '%CHROM\\t%POS\\t%INFO/DP\\n'
  rsomics-vcf query calls.vcf.gz -s S1,S2 -f '%POS[\\t%SAMPLE=%GT]\\n'")]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Fields and literals to print for each variant
    #[arg(short = 'f', long, value_name = "FORMAT")]
    format: String,

    /// Write results to this file instead of standard output
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        default_value = "-",
        hide_default_value = true
    )]
    output: PathBuf,

    /// Comma-separated samples to include; prefix the first name with ^ to exclude
    #[arg(
        short = 's',
        long,
        value_name = "LIST",
        conflicts_with = "samples_file"
    )]
    samples: Option<String>,

    /// File containing one sample name per line
    #[arg(short = 'S', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,

    /// Print a header; repeat to omit column indices
    #[arg(short = 'H', long, action = clap::ArgAction::Count)]
    print_header: u8,

    /// Do not append a newline when FORMAT contains none
    #[arg(short = 'N', long)]
    disable_automatic_newline: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because query data otherwise uses standard output".to_owned(),
        ));
    }

    let options = query::Options {
        format: arguments.format,
        samples: read_samples(arguments.samples, arguments.samples_file.as_deref())?,
        header: match arguments.print_header {
            0 => HeaderMode::None,
            1 => HeaderMode::Indexed,
            _ => HeaderMode::Plain,
        },
        automatic_newline: !arguments.disable_automatic_newline,
    };

    let summary = if arguments.output == Path::new("-") {
        query::write(&arguments.input, &options, io::stdout().lock())?
    } else {
        write_atomic(&arguments.output, |output| {
            query::write(&arguments.input, &options, output)
        })?
    };
    Ok(CommandOutput::Query { summary })
}

fn read_samples(list: Option<String>, file: Option<&Path>) -> Result<Option<Vec<String>>> {
    if let Some(list) = list {
        let samples: Vec<_> = list.split(',').map(str::to_owned).collect();
        if samples.iter().any(String::is_empty) {
            return Err(RsomicsError::InvalidInput(
                "sample list contains an empty name".to_owned(),
            ));
        }
        return Ok(Some(samples));
    }
    let Some(file) = file else {
        return Ok(None);
    };
    let content = std::fs::read_to_string(file)
        .rs_with_context(|| format!("reading sample list {}", file.display()))?;
    let samples: Vec<_> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
        .collect();
    if samples.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "sample list is empty: {}",
            file.display()
        )));
    }
    Ok(Some(samples))
}
