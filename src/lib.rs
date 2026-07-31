#![deny(unsafe_code)]

pub mod head;

mod cli;
mod commands;
mod format;

#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
