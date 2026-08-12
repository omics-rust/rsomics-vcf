use std::process;

use clap::{Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, Validation, run_validation};
use serde::Serialize;

use crate::{commands, head};

const META: ToolMeta = ToolMeta {
    name: "rsomics-vcf",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-vcf",
    version,
    about = "VCF and BCF workflows",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Filter VCF or BCF records with typed expressions and genomic masks
    Filter(commands::filter::Arguments),

    /// Print VCF headers and the first variant records
    Head(commands::head::Arguments),

    /// Create or inspect random-access indexes
    Index(commands::index::Arguments),

    #[cfg(feature = "norm-preview")]
    /// Normalize variants against an indexed reference
    Norm(commands::norm::Arguments),

    /// Extract variant fields with a typed format string
    Query(commands::query::Arguments),

    /// Validate VCF or BCF structure and typed values
    Validate(commands::validate::Arguments),

    /// Convert, subset, and select VCF or BCF records
    View(commands::view::Arguments),
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub(crate) enum CommandOutput {
    Filter {
        summary: crate::filter::Summary,
    },
    Head {
        summary: head::Summary,
    },
    Index {
        outcome: crate::index::Outcome,
    },
    #[cfg(feature = "norm-preview")]
    Norm {
        summary: crate::norm::Summary,
    },
    Query {
        summary: crate::query::Summary,
    },
    Validate {
        report: crate::validate::Report,
    },
    View {
        summary: crate::view::Summary,
    },
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_validation(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<Validation<CommandOutput>> {
    match cli.command {
        Command::Filter(arguments) => {
            commands::filter::execute(arguments, cli.output.json).map(Validation::Valid)
        }
        Command::Head(arguments) => {
            commands::head::execute(arguments, cli.output.json).map(Validation::Valid)
        }
        Command::Index(arguments) => {
            commands::index::execute(arguments, cli.output.json).map(Validation::Valid)
        }
        #[cfg(feature = "norm-preview")]
        Command::Norm(arguments) => {
            commands::norm::execute(arguments, cli.output.json).map(Validation::Valid)
        }
        Command::Query(arguments) => {
            commands::query::execute(arguments, cli.output.json).map(Validation::Valid)
        }
        Command::Validate(arguments) => commands::validate::execute(arguments, cli.output.json),
        Command::View(arguments) => {
            commands::view::execute(arguments, cli.output.json).map(Validation::Valid)
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn head_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "head", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input VCF or BCF file"), "{help}");
        assert!(help.contains("-H, --headers <INT>"), "{help}");
        assert!(help.contains("-s, --samples <INT>"), "{help}");
    }

    #[test]
    fn filter_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "filter", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input VCF or BCF file"), "{help}");
        assert!(help.contains("-i, --include <EXPR>"), "{help}");
        assert!(help.contains("--threads <INT>"), "{help}");
        assert!(help.contains("-O, --output-type <TYPE>"), "{help}");
    }

    #[cfg(feature = "norm-preview")]
    #[test]
    fn norm_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "norm", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input VCF or BCF file"), "{help}");
        assert!(help.contains("-f, --fasta-ref <FILE>"), "{help}");
        assert!(help.contains("-m, --split-multiallelic"), "{help}");
        assert!(help.contains("--strict-filter"), "{help}");
        assert!(help.contains("--split-overlaps <MODE>"), "{help}");
        assert!(help.contains("-a, --atomize"), "{help}");
        assert!(help.contains("--atom-overlaps <'*'|'.'>"), "{help}");
        assert!(help.contains("--old-rec-tag <TAG>"), "{help}");
        assert!(help.contains("--remove-duplicates <POLICY>"), "{help}");
        assert!(help.contains("--keep-sum <TAG>"), "{help}");
        assert!(help.contains("--check-ref <MODE>"), "{help}");
        assert!(help.contains("-i, --include <EXPR>"), "{help}");
        assert!(help.contains("-e, --exclude <EXPR>"), "{help}");
        assert!(help.contains("-r, --regions <REGIONS>"), "{help}");
        assert!(help.contains("--regions-overlap <MODE>"), "{help}");
        assert!(help.contains("-t, --targets <REGIONS>"), "{help}");
        assert!(help.contains("--targets-overlap <MODE>"), "{help}");
        assert!(help.contains("-w, --site-window <INT>"), "{help}");
        assert!(help.contains("-O, --output-type <TYPE>"), "{help}");
    }

    #[test]
    fn query_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "query", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Fields and literals"), "{help}");
        assert!(help.contains("-f, --format <FORMAT>"), "{help}");
        assert!(help.contains("-S, --samples-file <FILE>"), "{help}");
    }

    #[test]
    fn index_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "index", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(
            help.contains("Input BGZF-compressed VCF or BCF file"),
            "{help}"
        );
        assert!(help.contains("-m, --min-shift <INT>"), "{help}");
        assert!(help.contains("-s, --stats"), "{help}");
        assert!(help.contains("-n, --nrecords"), "{help}");
    }

    #[test]
    fn validate_help_uses_family_layout() {
        let error =
            rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "validate", "--help"])
                .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input VCF or BCF file"), "{help}");
        assert!(help.contains("--max-diagnostics <INT>"), "{help}");
        assert!(help.contains("--require-evidence"), "{help}");
    }

    #[test]
    fn view_help_uses_family_layout() {
        let error = rsomics_help::try_parse_from::<Cli, _, _>(["rsomics-vcf", "view", "--help"])
            .unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input VCF or BCF file"), "{help}");
        assert!(help.contains("-O, --output-type <TYPE>"), "{help}");
        assert!(help.contains("-s, --samples <LIST>"), "{help}");
        assert!(help.contains("-r, --regions <REGIONS>"), "{help}");
    }
}
