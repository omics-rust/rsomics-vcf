use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use rsomics_common::{AtomicFile, Result, RsomicsError, reject_output_alias};

use crate::cli::CommandOutput;
use crate::commands::variant::OutputType;
use crate::setgt::{self, Options, Query, Replacement, Target};

#[derive(Debug, Args)]
#[command(after_help = "\
Targets:
  .                genotypes with one or more missing alleles
  ./x              partially missing genotypes
  ./.              completely missing genotypes
  a                 all genotypes
  q                 samples selected by --include or --exclude
  b:TAG<CMP>VALUE  heterozygotes selected by a two-tailed binomial probability
  r:FLOAT           random fraction, alone or after one principal target

New genotypes:
  .  0  0p  m  mp  M  Mp  X  p  u  i  c:GT

Output types:
  v  uncompressed VCF
  z  BGZF-compressed VCF
  b  BGZF-compressed BCF
  u  uncompressed BCF

Examples:
  rsomics-vcf setgt -t . -n 0 calls.vcf.gz
  rsomics-vcf setgt -t q -i 'FMT/DP < 10' -n . -O z -o edited.vcf.gz calls.bcf")]
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

    /// Target genotype selector; repeat once only to add r:FLOAT
    #[arg(
        short = 't',
        long = "target-gt",
        value_name = "TARGET",
        required = true
    )]
    targets: Vec<String>,

    /// Replacement genotype form
    #[arg(short = 'n', long = "new-gt", value_name = "TYPE", required = true)]
    replacement: String,

    /// Exclude samples for which EXPR is true; requires target q
    #[arg(short = 'e', long, value_name = "EXPR", conflicts_with = "include")]
    exclude: Option<String>,

    /// Include samples for which EXPR is true; requires target q
    #[arg(short = 'i', long, value_name = "EXPR")]
    include: Option<String>,

    /// Signed seed for a random target
    #[arg(short = 's', long, value_name = "INT")]
    seed: Option<i64>,

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
    if arguments.output != Path::new("-") {
        reject_output_alias(&arguments.output, [arguments.input.as_path()])?;
    }

    let target = Target::parse(&arguments.targets)?;
    let replacement = Replacement::parse(&arguments.replacement)?;
    let query = match (arguments.include, arguments.exclude) {
        (Some(source), None) => Some(Query::Include(source)),
        (None, Some(source)) => Some(Query::Exclude(source)),
        (None, None) => None,
        (Some(_), Some(_)) => unreachable!(),
    };
    let is_query = matches!(target.principal, setgt::Principal::Query);
    if is_query != query.is_some() {
        return Err(RsomicsError::ConfigError(if is_query {
            "target q requires exactly one --include or --exclude expression".to_owned()
        } else {
            "--include and --exclude require target q".to_owned()
        }));
    }
    if arguments.seed.is_some() && target.random_fraction.is_none() {
        return Err(RsomicsError::ConfigError(
            "--seed requires a random target r:FLOAT".to_owned(),
        ));
    }
    let options = Options {
        output_format: arguments.output_type.into(),
        target,
        replacement,
        query,
        seed: arguments.seed.unwrap_or(0),
    };

    let summary = if arguments.output == Path::new("-") {
        if arguments.threads == 0 {
            setgt::write(&arguments.input, &options, io::stdout().lock())?
        } else {
            setgt::write_parallel(&arguments.input, &options, io::stdout(), arguments.threads)?
        }
    } else {
        let mut transaction = AtomicFile::new(&arguments.output)?;
        let summary = if arguments.threads == 0 {
            setgt::write(&arguments.input, &options, transaction.file_mut())?
        } else {
            setgt::write_parallel(
                &arguments.input,
                &options,
                transaction.reopen()?,
                arguments.threads,
            )?
        };
        transaction.commit()?;
        summary
    };
    Ok(CommandOutput::Setgt { summary })
}
