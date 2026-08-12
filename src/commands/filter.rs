use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{AtomicFile, Result, RsomicsError};

use crate::cli::CommandOutput;
use crate::commands::variant::{OutputType, Overlap, read_mask, read_regions, read_targets};
use crate::filter::{
    self, AnnotationMode, GapOptions, GenotypeReplacement, Logic, Mask, Options, SnpGap,
    StreamOptions,
};
use crate::variant_type;

#[derive(Debug, Args)]
#[command(after_help = "\
Output types:
  v  uncompressed VCF
  z  BGZF-compressed VCF
  b  BGZF-compressed BCF
  u  uncompressed BCF

Examples:
  rsomics-vcf filter calls.vcf.gz -i 'QUAL >= 20'
  rsomics-vcf filter calls.bcf -e 'FMT/DP < 10' -s LowDepth -O z -o filtered.vcf.gz
  rsomics-vcf filter calls.vcf.gz -g 3:indel,mnp -G 5 -r chr1:1-100000")]
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

    /// Exclude records or samples for which EXPR is true
    #[arg(short = 'e', long, value_name = "EXPR", conflicts_with = "include")]
    exclude: Option<String>,

    /// Include only records or samples for which EXPR is true
    #[arg(short = 'i', long, value_name = "EXPR")]
    include: Option<String>,

    /// Annotate failures with FILTER instead of removing them; + chooses a unique name
    #[arg(short = 's', long, value_name = "FILTER")]
    soft_filter: Option<String>,

    /// FILTER annotation mode: + adds to existing values and x resets passing sites
    #[arg(short = 'm', long, value_name = "+|x|+x")]
    mode: Option<String>,

    /// Set genotypes of failed samples to missing (.) or reference (0)
    #[arg(short = 'S', long = "set-GTs", value_name = ".|0")]
    set_genotypes: Option<String>,

    /// Filter SNPs within INT bases of selected variant types
    #[arg(short = 'g', long = "SnpGap", value_name = "INT[:TYPE]")]
    snp_gap: Option<String>,

    /// Filter clustered indels separated by at most INT bases
    #[arg(short = 'G', long = "IndelGap", value_name = "INT")]
    indel_gap: Option<usize>,

    /// Soft-filter comma-separated regions; prefix with ^ to negate
    #[arg(long, value_name = "REGIONS", conflicts_with = "mask_file")]
    mask: Option<String>,

    /// Soft-filter regions from FILE; prefix the path with ^ to negate
    #[arg(short = 'M', long, value_name = "FILE")]
    mask_file: Option<PathBuf>,

    /// Mask overlap rule
    #[arg(long, value_name = "MODE", default_value = "record")]
    mask_overlap: Overlap,

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

    /// Use INT BGZF compression workers; 0 selects serial output
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

    let (expression, logic) = if let Some(expression) = arguments.include {
        (Some(expression), Logic::Include)
    } else {
        (arguments.exclude, Logic::Exclude)
    };
    let mask = read_mask(
        arguments.mask,
        arguments.mask_file.as_deref(),
        arguments.mask_overlap.into(),
    )?
    .map(|(regions, negate)| Mask::new(regions, negate));
    let options = StreamOptions {
        output_format: arguments.output_type.into(),
        expression,
        filter: Options {
            logic,
            mask,
            soft_filter: arguments.soft_filter,
            mode: parse_mode(arguments.mode.as_deref())?,
            set_genotypes: parse_genotype_replacement(arguments.set_genotypes.as_deref())?,
        },
        gaps: GapOptions {
            snp_gap: arguments
                .snp_gap
                .as_deref()
                .map(parse_snp_gap)
                .transpose()?,
            indel_gap: arguments.indel_gap,
            soft: false,
        },
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
    };

    let summary = if arguments.output == Path::new("-") {
        if arguments.threads == 0 {
            filter::write(&arguments.input, &options, io::stdout().lock())?
        } else {
            filter::write_parallel(&arguments.input, &options, io::stdout(), arguments.threads)?
        }
    } else {
        let mut transaction = AtomicFile::new(&arguments.output)?;
        let summary = if arguments.threads == 0 {
            filter::write(&arguments.input, &options, transaction.file_mut())?
        } else {
            filter::write_parallel(
                &arguments.input,
                &options,
                transaction.reopen()?,
                arguments.threads,
            )?
        };
        transaction.commit()?;
        summary
    };
    Ok(CommandOutput::Filter { summary })
}

fn parse_mode(value: Option<&str>) -> Result<AnnotationMode> {
    match value {
        None => Ok(AnnotationMode::Replace),
        Some("+") => Ok(AnnotationMode::Add),
        Some("x") => Ok(AnnotationMode::ResetPassed),
        Some("+x" | "x+") => Ok(AnnotationMode::AddAndResetPassed),
        Some(value) => Err(RsomicsError::ConfigError(format!(
            "invalid filter annotation mode: {value}"
        ))),
    }
}

fn parse_genotype_replacement(value: Option<&str>) -> Result<Option<GenotypeReplacement>> {
    match value {
        None => Ok(None),
        Some(".") => Ok(Some(GenotypeReplacement::Missing)),
        Some("0") => Ok(Some(GenotypeReplacement::Reference)),
        Some(value) => Err(RsomicsError::ConfigError(format!(
            "invalid failed-genotype replacement: {value}"
        ))),
    }
}

fn parse_snp_gap(value: &str) -> Result<SnpGap> {
    let (distance, types) = value
        .split_once(':')
        .map_or((value, "indel"), |(distance, types)| (distance, types));
    let distance = distance
        .parse::<usize>()
        .map_err(|_| RsomicsError::ConfigError(format!("invalid SnpGap distance: {distance}")))?;
    let types = variant_type::parse_mask(types).map_err(RsomicsError::ConfigError)?;
    if types & (variant_type::SNP | variant_type::REF) != 0 {
        return Err(RsomicsError::ConfigError(
            "SnpGap types are indel, mnp, bnd, other, or overlap".to_owned(),
        ));
    }
    Ok(SnpGap { distance, types })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_filter_modes_genotype_replacements_and_gap_types() {
        assert_eq!(
            parse_mode(Some("+x")).unwrap(),
            AnnotationMode::AddAndResetPassed
        );
        assert_eq!(
            parse_genotype_replacement(Some(".")).unwrap(),
            Some(GenotypeReplacement::Missing)
        );
        let gap = parse_snp_gap("3:indel,mnp").unwrap();
        assert_eq!(gap.distance, 3);
        assert_eq!(gap.types, variant_type::INDEL | variant_type::MNP);
        assert!(parse_snp_gap("3:snp").is_err());
    }
}
