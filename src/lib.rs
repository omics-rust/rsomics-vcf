#![deny(unsafe_code)]

#[allow(dead_code)]
mod annotate;
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
mod norm;
mod query_bcf;
mod query_format;
mod regions;
#[allow(dead_code)]
mod reheader;
mod validation;
mod variant_type;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
