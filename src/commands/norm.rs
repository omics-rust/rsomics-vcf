use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Result, RsomicsError, write_atomic};

use crate::cli::CommandOutput;
use crate::commands::variant::{OutputType, Overlap, read_regions, read_targets};
use crate::filter::Logic;
use crate::norm::{self, MismatchPolicy, Options};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum CheckReference {
    Exit,
    Warn,
    Skip,
    Fix,
}

#[derive(Clone, Copy, Debug)]
enum AtomOverlaps {
    Star,
    Missing,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum RemoveDuplicates {
    Snps,
    Indels,
    Both,
    All,
    Exact,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SplitOverlaps {
    Reference,
    Missing,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum JoinMultiallelic {
    Snps,
    Indels,
    Both,
    Any,
}

impl From<JoinMultiallelic> for norm::JoinPolicy {
    fn from(value: JoinMultiallelic) -> Self {
        match value {
            JoinMultiallelic::Snps => Self::Snps,
            JoinMultiallelic::Indels => Self::Indels,
            JoinMultiallelic::Both => Self::Both,
            JoinMultiallelic::Any => Self::Any,
        }
    }
}

impl From<RemoveDuplicates> for norm::DuplicatePolicy {
    fn from(value: RemoveDuplicates) -> Self {
        match value {
            RemoveDuplicates::Snps => Self::Snps,
            RemoveDuplicates::Indels => Self::Indels,
            RemoveDuplicates::Both => Self::Both,
            RemoveDuplicates::All => Self::All,
            RemoveDuplicates::Exact => Self::Exact,
        }
    }
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
            CheckReference::Fix => Self::Fix,
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

    /// Join records at the same position into multiallelic records
    #[arg(long, value_name = "MODE", conflicts_with = "split_multiallelic")]
    join_multiallelic: Option<JoinMultiallelic>,

    /// Reset prior filters when a later joined record is PASS
    #[arg(long, requires = "join_multiallelic")]
    strict_filter: bool,

    /// Replacement for non-selected ALT alleles while splitting
    #[arg(long, value_name = "MODE", requires = "split_multiallelic")]
    split_overlaps: Option<SplitOverlaps>,

    /// Decompose complex variants into primitive records
    #[arg(short = 'a', long)]
    atomize: bool,

    /// Allele used for conflicts between overlapping atoms
    #[arg(long, value_name = "'*'|'.'", requires = "atomize")]
    atom_overlaps: Option<AtomOverlaps>,

    /// INFO tag recording each atom's original variant and ALT index
    #[arg(long, value_name = "TAG", requires = "atomize")]
    old_rec_tag: Option<String>,

    /// Remove later records matching this duplicate policy
    #[arg(long, value_name = "POLICY")]
    remove_duplicates: Option<RemoveDuplicates>,

    /// Preserve FORMAT/AD depth sums while splitting
    #[arg(long, value_name = "TAG", requires = "split_multiallelic")]
    keep_sum: Option<String>,

    /// REF mismatch behavior: exit, warn, skip, or fix
    #[arg(long, value_name = "MODE", requires = "reference")]
    check_ref: Option<CheckReference>,

    /// Do not normalize records for which EXPR is true
    #[arg(short = 'e', long, value_name = "EXPR", conflicts_with = "include")]
    exclude: Option<String>,

    /// Normalize only records for which EXPR is true
    #[arg(short = 'i', long, value_name = "EXPR")]
    include: Option<String>,

    /// Indexed genomic regions, separated by commas
    #[arg(
        short = 'r',
        long,
        value_name = "REGIONS",
        conflicts_with = "regions_file"
    )]
    regions: Option<String>,

    /// File containing one indexed genomic region per line
    #[arg(short = 'R', long, value_name = "FILE")]
    regions_file: Option<PathBuf>,

    /// Indexed-region overlap rule
    #[arg(long, value_name = "MODE", default_value = "record")]
    regions_overlap: Overlap,

    /// Streaming target regions, separated by commas; prefix with ^ to exclude
    #[arg(
        short = 't',
        long,
        value_name = "REGIONS",
        conflicts_with = "targets_file"
    )]
    targets: Option<String>,

    /// File containing one streaming target region per line; prefix the path with ^ to exclude
    #[arg(short = 'T', long, value_name = "FILE")]
    targets_file: Option<PathBuf>,

    /// Target overlap rule
    #[arg(long, value_name = "MODE", default_value = "pos")]
    targets_overlap: Overlap,

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
    if arguments.reference.is_none()
        && !arguments.split_multiallelic
        && arguments.join_multiallelic.is_none()
        && !arguments.atomize
        && arguments.remove_duplicates.is_none()
    {
        return Err(RsomicsError::ConfigError(
            "norm requires --fasta-ref, --split-multiallelic, --join-multiallelic, --atomize, --remove-duplicates, or a combination".to_owned(),
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
    let (expression, expression_logic) = if let Some(expression) = arguments.include {
        (Some(expression), Logic::Include)
    } else {
        (arguments.exclude, Logic::Exclude)
    };
    let options = Options {
        reference: arguments.reference,
        expression,
        expression_logic,
        regions: read_regions(
            arguments.regions,
            arguments.regions_file.as_deref(),
            arguments.regions_overlap.into(),
        )?,
        targets: read_targets(
            arguments.targets,
            arguments.targets_file.as_deref(),
            arguments.targets_overlap.into(),
        )?,
        split_multiallelic: arguments.split_multiallelic,
        join_multiallelic: arguments.join_multiallelic.map(Into::into),
        strict_filter: arguments.strict_filter,
        split_overlaps_missing: matches!(arguments.split_overlaps, Some(SplitOverlaps::Missing)),
        mismatch_policy: arguments.check_ref.unwrap_or(CheckReference::Exit).into(),
        atomize: arguments.atomize,
        atom_overlaps_star: !matches!(arguments.atom_overlaps, Some(AtomOverlaps::Missing)),
        old_record_tag: arguments.old_rec_tag,
        duplicate_policy: arguments.remove_duplicates.map(Into::into),
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
