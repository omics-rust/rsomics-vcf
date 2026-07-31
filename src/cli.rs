use std::process;

use clap::{Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, ToolMeta, run as run_tool};
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
    /// Print VCF headers and the first variant records
    Head(commands::head::Arguments),
}

#[derive(Debug, Serialize)]
#[serde(tag = "command", rename_all = "kebab-case")]
pub(crate) enum CommandOutput {
    Head { summary: head::Summary },
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<CommandOutput> {
    match cli.command {
        Command::Head(arguments) => commands::head::execute(arguments, cli.output.json),
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
}
