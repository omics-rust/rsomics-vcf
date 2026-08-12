use std::path::{Path, PathBuf};

use clap::ValueEnum;
use rsomics_common::{Context, Result, RsomicsError};

use crate::format::OutputFormat;
use crate::regions::{OverlapMode, RegionSelection, RegionSet};

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum OutputType {
    V,
    Z,
    B,
    U,
}

impl From<OutputType> for OutputFormat {
    fn from(value: OutputType) -> Self {
        match value {
            OutputType::V => Self::Vcf,
            OutputType::Z => Self::VcfBgzf,
            OutputType::B => Self::Bcf,
            OutputType::U => Self::BcfRaw,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum Overlap {
    Pos,
    Record,
    Variant,
}

impl From<Overlap> for OverlapMode {
    fn from(value: Overlap) -> Self {
        match value {
            Overlap::Pos => Self::Position,
            Overlap::Record => Self::Record,
            Overlap::Variant => Self::Variant,
        }
    }
}

pub(crate) fn read_regions(
    list: Option<String>,
    file: Option<&Path>,
    overlap: OverlapMode,
) -> Result<Option<RegionSet>> {
    let Some(values) = read_list(list, file, "region")? else {
        return Ok(None);
    };
    RegionSet::parse(values, overlap).map(Some)
}

pub(crate) fn read_targets(
    list: Option<String>,
    file: Option<&Path>,
    overlap: OverlapMode,
) -> Result<Option<RegionSelection>> {
    let Some((values, exclude)) = read_excludable_list(list, file, "target region")? else {
        return Ok(None);
    };
    RegionSelection::parse(values, overlap, exclude).map(Some)
}

pub(crate) fn read_mask(
    list: Option<String>,
    file: Option<&Path>,
    overlap: OverlapMode,
) -> Result<Option<(RegionSet, bool)>> {
    let Some((values, negate)) = read_excludable_list(list, file, "mask region")? else {
        return Ok(None);
    };
    RegionSet::parse(values, overlap).map(|regions| Some((regions, negate)))
}

pub(crate) fn read_list(
    list: Option<String>,
    file: Option<&Path>,
    kind: &str,
) -> Result<Option<Vec<String>>> {
    if let Some(list) = list {
        return Ok(Some(list.split(',').map(str::to_owned).collect()));
    }
    let Some(file) = file else {
        return Ok(None);
    };
    let content = std::fs::read_to_string(file)
        .rs_with_context(|| format!("reading {kind} list {}", file.display()))?;
    let values: Vec<_> = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    if values.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "{kind} list is empty: {}",
            file.display()
        )));
    }
    Ok(Some(values))
}

fn read_excludable_list(
    list: Option<String>,
    file: Option<&Path>,
    kind: &str,
) -> Result<Option<(Vec<String>, bool)>> {
    let list_supplied = list.is_some();
    let mut exclude = false;
    let mut selected_file = file.map(Path::to_path_buf);
    if let Some(path) = file
        && let Some(value) = path.to_str().and_then(|value| value.strip_prefix('^'))
    {
        if value.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "{kind} file path is empty after ^"
            )));
        }
        selected_file = Some(PathBuf::from(value));
        exclude = true;
    }

    let Some(mut values) = read_list(list, selected_file.as_deref(), kind)? else {
        return Ok(None);
    };
    if list_supplied && let Some(value) = values.first().and_then(|value| value.strip_prefix('^')) {
        let value = value.to_owned();
        if value.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "{kind} list is empty after ^"
            )));
        }
        values[0] = value;
        exclude = true;
    }
    Ok(Some((values, exclude)))
}
