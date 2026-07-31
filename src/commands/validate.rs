use std::io::{self, Write};
use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, RsomicsError, Validation};

use crate::cli::CommandOutput;
use crate::validate::{self, Report};
use crate::validation::{Diagnostic, Severity};

#[derive(Debug, Args)]
#[command(after_help = "\
Examples:
  rsomics-vcf validate calls.vcf.gz
  rsomics-vcf validate calls.bcf --require-evidence
  rsomics-vcf validate calls.vcf --json")]
pub(crate) struct Arguments {
    /// Input VCF or BCF file; omit or use - for standard input
    #[arg(value_name = "VARIANT", default_value = "-")]
    input: PathBuf,

    /// Keep at most this many diagnostics in the report
    #[arg(long, value_name = "INT", default_value_t = 100)]
    max_diagnostics: usize,

    /// Require GT, INFO/AF, or both INFO/AC and INFO/AN on every record
    #[arg(long)]
    require_evidence: bool,
}

pub(crate) fn execute(arguments: Arguments, json: bool) -> Result<Validation<CommandOutput>> {
    let report = validate::check(
        &arguments.input,
        validate::Options {
            max_diagnostics: arguments.max_diagnostics,
            require_evidence: arguments.require_evidence,
        },
    )?;
    let valid = report.is_valid();

    if !json && valid {
        write_valid_report(&report)?;
    }

    let message = (!valid).then(|| {
        if json {
            summary(&report, "validation failed")
        } else {
            invalid_message(&report)
        }
    });
    let output = CommandOutput::Validate { report };
    if let Some(message) = message {
        Ok(Validation::Invalid {
            report: output,
            message,
        })
    } else {
        Ok(Validation::Valid(output))
    }
}

fn write_valid_report(report: &Report) -> Result<()> {
    let mut stderr = io::stderr().lock();
    for diagnostic in &report.diagnostics {
        writeln!(stderr, "{}", format_diagnostic(diagnostic)).map_err(RsomicsError::Io)?;
    }
    stderr.flush().map_err(RsomicsError::Io)?;

    let mut stdout = io::stdout().lock();
    writeln!(stdout, "{}", summary(report, "valid")).map_err(RsomicsError::Io)?;
    stdout.flush().map_err(RsomicsError::Io)
}

fn invalid_message(report: &Report) -> String {
    let mut message = summary(report, "validation failed");
    for diagnostic in &report.diagnostics {
        message.push_str("\n  ");
        message.push_str(&format_diagnostic(diagnostic));
    }
    if report.truncated {
        message.push_str("\n  additional diagnostics were omitted");
    }
    message
}

fn summary(report: &Report, status: &str) -> String {
    format!(
        "{status}: {} record{}, {} error{}, {} warning{}",
        report.records,
        plural(report.records),
        report.errors,
        plural(report.errors),
        report.warnings,
        plural(report.warnings)
    )
}

fn format_diagnostic(diagnostic: &Diagnostic) -> String {
    let severity = match diagnostic.severity {
        Severity::Error => "error",
        Severity::Warning => "warning",
    };
    let mut output = format!("{severity}[{}]", diagnostic.code);
    if let Some(line) = diagnostic.line {
        output.push_str(&format!(" line {line}"));
    }
    if let Some(field) = &diagnostic.field {
        output.push_str(&format!(", field {field}"));
    }
    output.push_str(": ");
    output.push_str(&diagnostic.message);
    output
}

fn plural(count: usize) -> &'static str {
    if count == 1 { "" } else { "s" }
}
