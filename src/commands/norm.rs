use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError, write_atomic};

use crate::cli::CommandOutput;
use crate::commands::variant::OutputType;
use crate::norm::{self, MismatchPolicy, Options};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum CheckReference {
    Exit,
    Warn,
    Skip,
}

#[derive(Clone, Copy, Debug)]
enum AtomOverlaps {
    Star,
    Missing,
}

impl ValueEnum for AtomOverlaps {
    fn value_variants<'a>() -> &'a [Self] {
        &[Self::Star, Self::Missing]
    }

    fn to_possible_value(&self) -> Option<clap::builder::PossibleValue> {
        Some(match self {
            Self::Star => clap::builder::PossibleValue::new("*"),
            Self::Missing => clap::builder::PossibleValue::new("."),
        })
    }
}

impl From<CheckReference> for MismatchPolicy {
    fn from(value: CheckReference) -> Self {
        match value {
            CheckReference::Exit => Self::Exit,
            CheckReference::Warn => Self::Warn,
            CheckReference::Skip => Self::Skip,
        }
    }
}

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

    /// Indexed FASTA reference for left alignment and REF validation
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE")]
    reference: Option<PathBuf>,

    /// Split multiallelic records into biallelic records
    #[arg(short = 'm', long)]
    split_multiallelic: bool,

    /// Decompose complex variants into primitive records
    #[arg(short = 'a', long)]
    atomize: bool,

    /// Allele used for conflicts between overlapping atoms
    #[arg(long, value_name = "'*'|'.'", requires = "atomize")]
    atom_overlaps: Option<AtomOverlaps>,

    /// Preserve FORMAT/AD depth sums while splitting
    #[arg(long, value_name = "TAG", requires = "split_multiallelic")]
    keep_sum: Option<String>,

    /// REF mismatch behavior: exit, warn, or skip
    #[arg(long, value_name = "MODE", requires = "reference")]
    check_ref: Option<CheckReference>,

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
    if arguments.reference.is_none() && !arguments.split_multiallelic && !arguments.atomize {
        return Err(RsomicsError::ConfigError(
            "norm requires --fasta-ref, --split-multiallelic, --atomize, or a combination"
                .to_owned(),
        ));
    }
    if arguments.keep_sum.as_deref().is_some_and(|tag| tag != "AD") {
        return Err(RsomicsError::ConfigError(
            "--keep-sum currently accepts only AD".to_owned(),
        ));
    }
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because variant data otherwise uses standard output"
                .to_owned(),
        ));
    }
    let options = Options {
        reference: arguments.reference,
        split_multiallelic: arguments.split_multiallelic,
        mismatch_policy: arguments.check_ref.unwrap_or(CheckReference::Exit).into(),
        atomize: arguments.atomize,
        atom_overlaps_star: !matches!(arguments.atom_overlaps, Some(AtomOverlaps::Missing)),
        keep_sum_ad: arguments.keep_sum.is_some(),
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
