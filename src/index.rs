mod bcf_record;
mod build;
mod csi;
mod stats;
mod vcf;

use std::path::{Path, PathBuf};

use rsomics_common::Result;
use serde::Serialize;

/// On-disk random-access index type.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum IndexKind {
    Csi,
    Tbi,
}

/// Variant stream format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum VariantFormat {
    Vcf,
    Bcf,
}

/// Index construction settings.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BuildOptions {
    pub kind: IndexKind,
    pub min_shift: u8,
    pub threads: usize,
    pub force: bool,
}

/// Completed index construction.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct BuildSummary {
    pub input: PathBuf,
    pub output: PathBuf,
    pub format: VariantFormat,
    pub kind: IndexKind,
    pub min_shift: u8,
    pub depth: u8,
    pub records: u64,
    pub reference_sequences: usize,
}

/// Existing-index inspection mode.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InspectMode {
    Total,
    PerContig { include_zero: bool },
}

/// Per-contig index metadata.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ContigStat {
    pub name: String,
    pub length: Option<usize>,
    pub records: Option<u64>,
}

/// Existing-index record counts.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct InspectReport {
    pub index: PathBuf,
    pub kind: IndexKind,
    pub total: u64,
    pub contigs: Vec<ContigStat>,
}

/// Result of the `index` command.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
#[serde(tag = "mode", rename_all = "kebab-case")]
pub enum Outcome {
    Build(BuildSummary),
    Inspect(InspectReport),
}

/// Returns the default sibling index path.
#[must_use]
pub fn default_output_path(input: &Path, kind: IndexKind) -> PathBuf {
    let extension = match kind {
        IndexKind::Csi => "csi",
        IndexKind::Tbi => "tbi",
    };
    PathBuf::from(format!("{}.{}", input.display(), extension))
}

/// Creates a random-access index without exposing a partial destination.
pub fn create(input: &Path, output: &Path, options: BuildOptions) -> Result<BuildSummary> {
    build::create(input, output, options)
}

/// Reads count metadata from an existing CSI or TBI.
pub fn inspect(input: &Path, mode: InspectMode) -> Result<InspectReport> {
    stats::inspect(input, mode)
}
