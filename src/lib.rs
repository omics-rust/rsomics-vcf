#![deny(unsafe_code)]

pub mod head;
pub mod index;
pub mod query;
pub mod validate;
pub mod view;

mod cli;
mod commands;
#[cfg(test)]
mod expression;
mod format;
mod query_bcf;
mod query_format;
mod validation;
mod variant_type;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
