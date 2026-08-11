use std::io;
use std::path::{Path, PathBuf};

use clap::{Args, ValueEnum};
use rsomics_common::{Context, Result, RsomicsError, write_atomic};

use crate::cli::CommandOutput;
use crate::view::{
    self, HeaderMode, IdSelection, Options, OutputFormat, OverlapMode, RegionSelection, RegionSet,
    SampleSelection, TypeSelection,
};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OutputType {
    V,
    Z,
    B,
    U,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Overlap {
    Pos,
    Record,
    Variant,
}

#[derive(Debug, Args)]
#[command(after_help = "\
Output types:
  v  uncompressed VCF
  z  BGZF-compressed VCF
  b  BGZF-compressed BCF
  u  uncompressed BCF

Variant types:
  snps, indels, mnps, ref, bnd, other, overlap

Examples:
  rsomics-vcf view calls.bcf -O z -o calls.vcf.gz
  rsomics-vcf view calls.vcf.gz -s tumor,normal -r chr1:1-100000")]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Write variants to this file instead of standard output
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

    /// Print only the header
    #[arg(long, conflicts_with = "no_header")]
    header_only: bool,

    /// Omit the VCF header; unavailable for BCF output
    #[arg(short = 'H', long)]
    no_header: bool,

    /// Comma-separated samples to include; prefix with ^ to exclude
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

    /// Ignore sample names that are absent from the input header
    #[arg(long)]
    force_samples: bool,

    /// Remove FORMAT and sample columns
    #[arg(short = 'G', long)]
    drop_genotypes: bool,

    /// Preserve existing INFO/AC and INFO/AN after sample selection
    #[arg(short = 'I', long)]
    no_update: bool,

    /// Keep records whose FILTER contains any listed value
    #[arg(short = 'f', long, value_name = "LIST", value_delimiter = ',')]
    apply_filters: Vec<String>,

    /// Keep records with a non-missing ID
    #[arg(short = 'k', long, conflicts_with = "novel")]
    known: bool,

    /// Keep records with a missing ID
    #[arg(short = 'n', long)]
    novel: bool,

    /// Keep records containing any listed variant type
    #[arg(
        short = 'v',
        long,
        value_name = "LIST",
        conflicts_with = "exclude_types"
    )]
    types: Option<String>,

    /// Remove records containing any listed variant type
    #[arg(short = 'V', long, value_name = "LIST")]
    exclude_types: Option<String>,

    /// Minimum number of REF plus ALT alleles
    #[arg(short = 'm', long, value_name = "INT")]
    min_alleles: Option<usize>,

    /// Maximum number of REF plus ALT alleles
    #[arg(short = 'M', long, value_name = "INT")]
    max_alleles: Option<usize>,

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
    #[arg(long, value_name = "MODE", default_value = "record")]
    targets_overlap: Overlap,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<CommandOutput> {
    if json && arguments.output == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "--json requires --output because variant data otherwise uses standard output"
                .to_owned(),
        ));
    }

    let mut samples = read_samples(arguments.samples, arguments.samples_file.as_deref())?;
    if let Some(selection) = &mut samples {
        selection.force = arguments.force_samples;
    } else if arguments.force_samples {
        return Err(RsomicsError::ConfigError(
            "--force-samples requires --samples or --samples-file".to_owned(),
        ));
    }
    let options = Options {
        output_format: match arguments.output_type {
            OutputType::V => OutputFormat::Vcf,
            OutputType::Z => OutputFormat::VcfBgzf,
            OutputType::B => OutputFormat::Bcf,
            OutputType::U => OutputFormat::BcfRaw,
        },
        header: if arguments.header_only {
            HeaderMode::HeaderOnly
        } else if arguments.no_header {
            HeaderMode::None
        } else {
            HeaderMode::Full
        },
        samples,
        drop_genotypes: arguments.drop_genotypes,
        update_info: !arguments.no_update,
        filters: arguments.apply_filters,
        ids: arguments
            .known
            .then_some(IdSelection::Known)
            .or(arguments.novel.then_some(IdSelection::Novel)),
        types: parse_types(arguments.types, arguments.exclude_types)?,
        min_alleles: arguments.min_alleles,
        max_alleles: arguments.max_alleles,
        targets: read_targets(
            arguments.targets,
            arguments.targets_file.as_deref(),
            arguments.targets_overlap.into(),
        )?,
        regions: read_regions(
            arguments.regions,
            arguments.regions_file.as_deref(),
            arguments.regions_overlap.into(),
        )?,
    };

    let summary = if arguments.output == Path::new("-") {
        view::write(&arguments.input, &options, io::stdout().lock())?
    } else {
        write_atomic(&arguments.output, |output| {
            view::write(&arguments.input, &options, output)
        })?
    };
    Ok(CommandOutput::View { summary })
}

fn parse_types(include: Option<String>, exclude: Option<String>) -> Result<Option<TypeSelection>> {
    let Some((values, include)) = include
        .map(|values| (values, true))
        .or_else(|| exclude.map(|values| (values, false)))
    else {
        return Ok(None);
    };
    if include {
        TypeSelection::include(&values).map(Some)
    } else {
        TypeSelection::exclude(&values).map(Some)
    }
}

fn read_samples(list: Option<String>, file: Option<&Path>) -> Result<Option<SampleSelection>> {
    let Some(raw) = read_list(list, file, "sample")? else {
        return Ok(None);
    };
    let (exclude, names) = if let Some(first) = raw.first().and_then(|name| name.strip_prefix('^'))
    {
        let mut names = raw.clone();
        names[0] = first.to_owned();
        (true, names)
    } else {
        (false, raw)
    };
    if names.iter().any(String::is_empty) {
        return Err(RsomicsError::InvalidInput(
            "sample list contains an empty name".to_owned(),
        ));
    }
    Ok(Some(SampleSelection {
        names,
        exclude,
        force: false,
    }))
}

fn read_regions(
    list: Option<String>,
    file: Option<&Path>,
    overlap: OverlapMode,
) -> Result<Option<RegionSet>> {
    let Some(values) = read_list(list, file, "region")? else {
        return Ok(None);
    };
    RegionSet::parse(values, overlap).map(Some)
}

fn read_targets(
    list: Option<String>,
    file: Option<&Path>,
    overlap: OverlapMode,
) -> Result<Option<RegionSelection>> {
    let list_supplied = list.is_some();
    let mut exclude = false;
    let mut target_file = file.map(Path::to_path_buf);
    if let Some(path) = file
        && let Some(value) = path.to_str().and_then(|value| value.strip_prefix('^'))
    {
        if value.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "target file path is empty after ^".to_owned(),
            ));
        }
        target_file = Some(PathBuf::from(value));
        exclude = true;
    }

    let Some(mut values) = read_list(list, target_file.as_deref(), "target region")? else {
        return Ok(None);
    };
    if list_supplied && let Some(value) = values.first().and_then(|value| value.strip_prefix('^')) {
        let value = value.to_owned();
        if value.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "target region list is empty after ^".to_owned(),
            ));
        }
        values[0] = value;
        exclude = true;
    }
    RegionSelection::parse(values, overlap, exclude).map(Some)
}

fn read_list(list: Option<String>, file: Option<&Path>, kind: &str) -> Result<Option<Vec<String>>> {
    if let Some(list) = list {
        return Ok(Some(list.split(',').map(str::to_owned).collect()));
    }
    let Some(file) = file else {
        return Ok(None);
    };
    let content = std::fs::read_to_string(file)
        .rs_with_context(|| format!("reading {kind} list {}", file.display()))?;
    let values: Vec<_> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    if values.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "{kind} list is empty: {}",
            file.display()
        )));
    }
    Ok(Some(values))
}

impl From<Overlap> for OverlapMode {
    fn from(value: Overlap) -> Self {
        match value {
            Overlap::Pos => Self::Position,
            Overlap::Record => Self::Record,
            Overlap::Variant => Self::Variant,
        }
    }
}
