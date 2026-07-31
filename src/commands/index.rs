use std::io::{self, Write};
use std::path::{Path, PathBuf};

use clap::{ArgGroup, Args};
use rsomics_common::{Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::index::{self, BuildOptions, IndexKind, InspectMode, Outcome};

#[derive(Debug, Args)]
#[command(
    group(
        ArgGroup::new("index_kind")
            .args(["csi", "tbi"])
            .multiple(false)
    ),
    after_help = "\
Examples:
  rsomics-vcf index calls.vcf.gz
  rsomics-vcf index --tbi calls.vcf.gz
  rsomics-vcf index --min-shift 18 calls.bcf
  rsomics-vcf index --stats calls.vcf.gz
  rsomics-vcf index --nrecords calls.bcf.csi"
)]
pub(crate) struct Arguments {
    /// Input BGZF-compressed VCF or BCF file, or an existing .csi/.tbi for statistics
    #[arg(value_name = "VARIANT")]
    input: PathBuf,

    /// Write a CSI index; this is the default
    #[arg(short = 'c', long)]
    csi: bool,

    /// Write a TBI index; VCF only
    #[arg(short = 't', long, conflicts_with = "min_shift")]
    tbi: bool,

    /// Overwrite an existing index
    #[arg(short = 'f', long, conflicts_with_all = ["stats", "nrecords"])]
    force: bool,

    /// Smallest CSI bin is 2^INT bases
    #[arg(
        short = 'm',
        long,
        value_name = "INT",
        default_value_t = 14,
        value_parser = clap::value_parser!(u8).range(1..=30),
        conflicts_with_all = ["stats", "nrecords"]
    )]
    min_shift: u8,

    /// Write the index to FILE
    #[arg(
        short = 'o',
        long,
        value_name = "FILE",
        conflicts_with_all = ["stats", "nrecords"]
    )]
    output: Option<PathBuf>,

    /// Use INT BGZF decompression workers
    #[arg(
        long,
        value_name = "INT",
        default_value_t = 0,
        conflicts_with_all = ["stats", "nrecords"]
    )]
    threads: usize,

    /// Print per-contig record counts from an existing index
    #[arg(
        short = 's',
        long,
        conflicts_with_all = ["nrecords", "csi", "tbi", "output"]
    )]
    stats: bool,

    /// Print the total record count from an existing index
    #[arg(
        short = 'n',
        long,
        conflicts_with_all = ["stats", "csi", "tbi", "output", "all"]
    )]
    nrecords: bool,

    /// Include contigs with zero records in --stats
    #[arg(short = 'a', long, requires = "stats")]
    all: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    let outcome = if arguments.stats || arguments.nrecords {
        let mode = if arguments.stats {
            InspectMode::PerContig {
                include_zero: arguments.all,
            }
        } else {
            InspectMode::Total
        };
        let report = index::inspect(&arguments.input, mode)?;
        if !json {
            write_report(&report, mode)?;
        }
        Outcome::Inspect(report)
    } else {
        if arguments.input == Path::new("-") && arguments.output.is_none() {
            return Err(RsomicsError::ConfigError(
                "--output is required when indexing standard input".to_owned(),
            ));
        }
        let kind = if arguments.tbi {
            IndexKind::Tbi
        } else {
            IndexKind::Csi
        };
        let output = arguments
            .output
            .unwrap_or_else(|| index::default_output_path(&arguments.input, kind));
        if output == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "index output must be a named file".to_owned(),
            ));
        }
        let summary = index::create(
            &arguments.input,
            &output,
            BuildOptions {
                kind,
                min_shift: arguments.min_shift,
                threads: arguments.threads,
                force: arguments.force,
            },
        )?;
        Outcome::Build(summary)
    };

    Ok(CommandOutput::Index { outcome })
}

fn write_report(report: &index::InspectReport, mode: InspectMode) -> Result<()> {
    let mut output = io::stdout().lock();
    match mode {
        InspectMode::Total => {
            writeln!(output, "{}", report.total).map_err(RsomicsError::Io)?;
        }
        InspectMode::PerContig { .. } => {
            for contig in &report.contigs {
                writeln!(
                    output,
                    "{}\t{}\t{}",
                    contig.name,
                    contig
                        .length
                        .map_or_else(|| ".".to_owned(), |length| length.to_string()),
                    contig
                        .records
                        .map_or_else(|| ".".to_owned(), |records| records.to_string())
                )
                .map_err(RsomicsError::Io)?;
            }
        }
    }
    output.flush().map_err(RsomicsError::Io)
}
