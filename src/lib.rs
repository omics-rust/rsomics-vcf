#![deny(unsafe_code)]

pub mod head;
pub mod query;

mod cli;
mod commands;
mod format;
mod query_bcf;
mod query_format;
mod variant_type;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
