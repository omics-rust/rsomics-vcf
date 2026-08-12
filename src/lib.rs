#![deny(unsafe_code)]

pub mod head;
pub mod index;
pub mod query;
pub mod validate;
pub mod view;

mod cli;
mod commands;
mod expression;
mod filter;
mod format;
#[cfg(feature = "norm-preview")]
mod norm;
mod query_bcf;
mod query_format;
mod regions;
mod validation;
mod variant_type;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
