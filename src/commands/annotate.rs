use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{AtomicFile, Result, RsomicsError, reject_output_alias};

use crate::{
    annotate::{
        self, ColumnSpec, HeaderOptions, MarkSites, Options, OverlapFractions, PairLogic,
        SampleRequest, SourceOptions,
    },
    cli::CommandOutput,
    commands::variant::{OutputType, Overlap, read_regions},
    filter::Logic,
};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Pairing {
    Snps,
    Indels,
    Both,
    All,
    Some,
    Exact,
    Id,
}

impl From<Pairing> for PairLogic {
    fn from(value: Pairing) -> Self {
        match value {
            Pairing::Snps => Self::Snps,
            Pairing::Indels => Self::Indels,
            Pairing::Both => Self::Both,
            Pairing::All => Self::All,
            Pairing::Some => Self::Some,
            Pairing::Exact => Self::Exact,
            Pairing::Id => Self::Id,
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

Column plans:
  VCF and BCF sources infer their coordinate and allele fields.
  Tabular sources require CHROM,POS or CHROM,FROM,TO before transfer columns.
  Transfer ID, QUAL, FILTER, INFO/TAG, or FORMAT/TAG; DEST:=SOURCE renames.
  Prefixes: . writes missing source values; + fills missing targets; .+ writes missing source values
  while filling missing targets; = appends; .= appends missing values; - updates existing targets.

Examples:
  rsomics-vcf annotate calls.vcf.gz -a db.vcf.gz -c ID,INFO/AF
  rsomics-vcf annotate calls.bcf -a depths.vcf.gz -c FORMAT/DP -s S1,S2 -O z -o annotated.vcf.gz
  rsomics-vcf annotate calls.vcf.gz -x INFO/OLD --rename-annotations names.tsv")]
pub(crate) struct Arguments {
    /// Target VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Sorted VCF, BCF, BED, or tab-delimited annotation source
    #[arg(short = 'a', long, value_name = "FILE")]
    annotations: Option<PathBuf>,

    /// Annotation match and transfer columns
    #[arg(
        short = 'c',
        long,
        value_name = "LIST",
        conflicts_with = "columns_file"
    )]
    columns: Option<String>,

    /// Read one annotation column name per line
    #[arg(short = 'C', long, value_name = "FILE")]
    columns_file: Option<PathBuf>,

    /// Append one INFO, FORMAT, FILTER, or contig definition
    #[arg(short = 'H', long, value_name = "LINE", action = clap::ArgAction::Append)]
    header_line: Vec<String>,

    /// Append supported header definitions from FILE
    #[arg(long, value_name = "FILE", conflicts_with = "header_line")]
    header_lines: Option<PathBuf>,

    /// Set ID from a site-format expression; prefix with + to fill missing IDs only
    #[arg(short = 'I', long, value_name = "FORMAT")]
    set_id: Option<String>,

    /// Remove fixed fields or selected INFO, FORMAT, and FILTER entries
    #[arg(short = 'x', long, value_name = "LIST")]
    remove: Option<String>,

    /// Rename chromosomes from two-column FILE
    #[arg(long, value_name = "FILE")]
    rename_chromosomes: Option<PathBuf>,

    /// Rename INFO, FORMAT, and FILTER identifiers from FILE
    #[arg(long, value_name = "FILE")]
    rename_annotations: Option<PathBuf>,

    /// Annotate only records for which EXPR is true
    #[arg(short = 'i', long, value_name = "EXPR", conflicts_with = "exclude")]
    include: Option<String>,

    /// Annotate only records for which EXPR is false
    #[arg(short = 'e', long, value_name = "EXPR")]
    exclude: Option<String>,

    /// Keep expression-rejected records unchanged instead of dropping them
    #[arg(short = 'k', long)]
    keep_sites: bool,

    /// Add INFO/TAG to matched (+TAG) or unmatched (-TAG) records
    #[arg(
        short = 'm',
        long,
        value_name = "+TAG|-TAG",
        allow_hyphen_values = true
    )]
    mark_sites: Option<String>,

    /// Required overlap fraction as ANN:VCF
    #[arg(long, value_name = "ANN:VCF")]
    min_overlap: Option<String>,

    /// Variant pairing rule
    #[arg(long, value_name = "MODE")]
    pair_logic: Option<Pairing>,

    /// Comma-separated annotation samples; prefix with ^ to exclude
    #[arg(
        short = 's',
        long,
        value_name = "LIST",
        conflicts_with = "samples_file"
    )]
    samples: Option<String>,

    /// Annotation samples from FILE; prefix the path with ^ to exclude
    #[arg(short = 'S', long, value_name = "FILE")]
    samples_file: Option<PathBuf>,

    /// Indexed target regions, separated by commas
    #[arg(
        short = 'r',
        long,
        value_name = "REGIONS",
        conflicts_with = "regions_file"
    )]
    regions: Option<String>,

    /// File containing one indexed target region per line
    #[arg(short = 'R', long, value_name = "FILE")]
    regions_file: Option<PathBuf>,

    /// Indexed-region overlap rule
    #[arg(long, value_name = "MODE", default_value = "record")]
    regions_overlap: Overlap,

    /// Write annotated variants to FILE instead of standard output
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

    /// Use INT BGZF compression workers; 0 selects serial output
    #[arg(long, value_name = "INT", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn execute(mut arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because variant data otherwise uses standard output"
                .to_owned(),
        ));
    }
    let columns = match (arguments.columns.take(), arguments.columns_file.as_deref()) {
        (Some(columns), None) => Some(ColumnSpec::parse(&columns)?),
        (None, Some(path)) => Some(ColumnSpec::from_file(path)?),
        (None, None) => None,
        (Some(_), Some(_)) => unreachable!(),
    };
    let source = match (arguments.annotations.as_ref(), columns) {
        (Some(path), Some(columns)) => Some(SourceOptions {
            path: path.clone(),
            columns,
            samples: read_samples(
                arguments.samples.as_deref(),
                arguments.samples_file.as_deref(),
            )?,
            pair_logic: arguments.pair_logic.unwrap_or(Pairing::Some).into(),
            min_overlap: parse_overlap(arguments.min_overlap.as_deref().unwrap_or("0:0"))?,
        }),
        (Some(_), None) => {
            return Err(RsomicsError::ConfigError(
                "--annotations requires --columns or --columns-file".to_owned(),
            ));
        }
        (None, Some(_)) => {
            return Err(RsomicsError::ConfigError(
                "--columns and --columns-file require --annotations".to_owned(),
            ));
        }
        (None, None) => None,
    };
    if source.is_none() && (arguments.samples.is_some() || arguments.samples_file.is_some()) {
        return Err(RsomicsError::ConfigError(
            "--samples and --samples-file require --annotations".to_owned(),
        ));
    }
    if source.is_none() && (arguments.pair_logic.is_some() || arguments.min_overlap.is_some()) {
        return Err(RsomicsError::ConfigError(
            "--pair-logic and --min-overlap require --annotations".to_owned(),
        ));
    }
    let mut appended = arguments.header_line;
    if let Some(path) = &arguments.header_lines {
        let content = std::fs::read_to_string(path).map_err(|error| {
            RsomicsError::Io(io::Error::new(
                error.kind(),
                format!(
                    "reading annotation header lines {}: {error}",
                    path.display()
                ),
            ))
        })?;
        appended.extend(
            content
                .lines()
                .map(str::trim)
                .filter(|line| {
                    !line.is_empty() && (!line.starts_with('#') || line.starts_with("##"))
                })
                .map(str::to_owned),
        );
    }
    let (expression, expression_logic) = if let Some(expression) = arguments.include {
        (Some(expression), Logic::Include)
    } else {
        (arguments.exclude, Logic::Exclude)
    };
    let options = Options {
        source,
        header: HeaderOptions {
            appended,
            remove: arguments.remove,
            rename_chromosomes: arguments.rename_chromosomes,
            rename_annotations: arguments.rename_annotations,
        },
        set_id: arguments.set_id,
        expression,
        expression_logic,
        keep_sites: arguments.keep_sites,
        mark_sites: arguments
            .mark_sites
            .as_deref()
            .map(parse_mark)
            .transpose()?,
        regions: read_regions(
            arguments.regions,
            arguments.regions_file.as_deref(),
            arguments.regions_overlap.into(),
        )?,
        output_format: arguments.output_type.into(),
    };

    let mut inputs = vec![arguments.input.as_path()];
    if let Some(path) = arguments.annotations.as_deref() {
        inputs.push(path);
    }
    reject_output_alias(&arguments.output, inputs)?;
    let summary = if arguments.output == Path::new("-") {
        if arguments.threads == 0 {
            annotate::write(&arguments.input, &options, io::stdout().lock())?
        } else {
            annotate::write_parallel(&arguments.input, &options, io::stdout(), arguments.threads)?
        }
    } else {
        let mut transaction = AtomicFile::new(&arguments.output)?;
        let summary = if arguments.threads == 0 {
            annotate::write(&arguments.input, &options, transaction.file_mut())?
        } else {
            annotate::write_parallel(
                &arguments.input,
                &options,
                transaction.reopen()?,
                arguments.threads,
            )?
        };
        transaction.commit()?;
        summary
    };
    Ok(CommandOutput::Annotate { summary })
}

fn read_samples(list: Option<&str>, file: Option<&Path>) -> Result<Option<SampleRequest>> {
    match (list, file) {
        (Some(list), None) => Ok(Some(SampleRequest::List(list.to_owned()))),
        (None, Some(path)) => Ok(Some(SampleRequest::File(path.to_path_buf()))),
        (None, None) => Ok(None),
        (Some(_), Some(_)) => unreachable!(),
    }
}

fn parse_overlap(source: &str) -> Result<OverlapFractions> {
    let (annotation, target) = source
        .split_once(':')
        .ok_or_else(|| RsomicsError::ConfigError("--min-overlap requires ANN:VCF".to_owned()))?;
    let fractions = OverlapFractions {
        annotation: annotation.parse().map_err(|_| {
            RsomicsError::ConfigError(format!("invalid annotation overlap fraction: {annotation}"))
        })?,
        target: target.parse().map_err(|_| {
            RsomicsError::ConfigError(format!("invalid target overlap fraction: {target}"))
        })?,
    };
    fractions.validate()?;
    Ok(fractions)
}

fn parse_mark(source: &str) -> Result<MarkSites> {
    let (present, tag) = if let Some(tag) = source.strip_prefix('+') {
        (true, tag)
    } else if let Some(tag) = source.strip_prefix('-') {
        (false, tag)
    } else {
        return Err(RsomicsError::ConfigError(
            "--mark-sites requires +TAG or -TAG".to_owned(),
        ));
    };
    if tag.is_empty() {
        return Err(RsomicsError::ConfigError(
            "--mark-sites tag must not be empty".to_owned(),
        ));
    }
    Ok(if present {
        MarkSites::Present(tag.to_owned())
    } else {
        MarkSites::Absent(tag.to_owned())
    })
}
