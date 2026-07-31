#![deny(unsafe_code)]

pub mod head;
pub mod query;
pub mod validate;

mod cli;
mod commands;
mod format;
mod query_bcf;
mod query_format;
mod validation;
mod variant_type;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
