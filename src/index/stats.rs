use std::collections::HashMap;
use std::path::{Path, PathBuf};

use noodles_csi::{
    self as csi, BinningIndex,
    binning_index::{
        self,
        index::{ReferenceSequence, reference_sequence::Metadata},
    },
};
use noodles_tabix as tabix;
use noodles_vcf as vcf;
use rsomics_common::{Context, Result, RsomicsError};

use super::{ContigStat, IndexKind, InspectMode, InspectReport, default_output_path};

enum Index {
    Csi(csi::Index),
    Tbi(tabix::Index),
}

struct Resolved {
    variant: Option<PathBuf>,
    index: PathBuf,
    kind: IndexKind,
}

pub(super) fn inspect(input: &Path, mode: InspectMode) -> Result<InspectReport> {
    let resolved = resolve(input)?;
    let loaded = match resolved.kind {
        IndexKind::Csi => csi::fs::read(&resolved.index)
            .map(Index::Csi)
            .rs_with_context(|| format!("reading CSI index {}", resolved.index.display()))?,
        IndexKind::Tbi => tabix::fs::read(&resolved.index)
            .map(Index::Tbi)
            .rs_with_context(|| format!("reading TBI index {}", resolved.index.display()))?,
    };
    let variant_header = resolved
        .variant
        .as_deref()
        .filter(|path| path.exists())
        .map(read_variant_header)
        .transpose()?;
    let lengths = variant_header
        .as_ref()
        .map(contig_lengths)
        .unwrap_or_default();

    let (names, counts) = metadata(&loaded);
    let mut contigs = Vec::with_capacity(counts.len());
    let mut total = 0u64;
    for (i, records) in counts.into_iter().enumerate() {
        if let Some(records) = records {
            total = total.checked_add(records).ok_or_else(|| {
                RsomicsError::InvalidInput("index record count overflow".to_owned())
            })?;
        }
        let name = names
            .get(i)
            .cloned()
            .or_else(|| {
                variant_header
                    .as_ref()
                    .and_then(|header| header.contigs().get_index(i).map(|(name, _)| name.clone()))
            })
            .unwrap_or_else(|| "n/a".to_owned());
        let include = match mode {
            InspectMode::Total => false,
            InspectMode::PerContig { include_zero } => {
                include_zero || records.is_some_and(|records| records > 0)
            }
        };
        if include {
            contigs.push(ContigStat {
                length: lengths.get(&name).copied().flatten(),
                name,
                records,
            });
        }
    }

    Ok(InspectReport {
        index: resolved.index,
        kind: resolved.kind,
        total,
        contigs,
    })
}

fn resolve(input: &Path) -> Result<Resolved> {
    let text = input.to_string_lossy();
    if let Some((variant, index)) = text.split_once("##idx##") {
        let index = PathBuf::from(index);
        return Ok(Resolved {
            variant: Some(PathBuf::from(variant)),
            kind: kind_from_path(&index)?,
            index,
        });
    }
    if let Some(kind) = optional_kind_from_path(input) {
        let variant = text
            .strip_suffix(".csi")
            .or_else(|| text.strip_suffix(".tbi"))
            .map(PathBuf::from)
            .filter(|path| path.exists());
        return Ok(Resolved {
            variant,
            index: input.to_path_buf(),
            kind,
        });
    }

    let csi = default_output_path(input, IndexKind::Csi);
    if csi.exists() {
        return Ok(Resolved {
            variant: Some(input.to_path_buf()),
            index: csi,
            kind: IndexKind::Csi,
        });
    }
    let tbi = default_output_path(input, IndexKind::Tbi);
    if tbi.exists() {
        return Ok(Resolved {
            variant: Some(input.to_path_buf()),
            index: tbi,
            kind: IndexKind::Tbi,
        });
    }
    Err(RsomicsError::InvalidInput(format!(
        "no CSI or TBI index found for {}",
        input.display()
    )))
}

fn kind_from_path(path: &Path) -> Result<IndexKind> {
    optional_kind_from_path(path).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "index path must end in .csi or .tbi: {}",
            path.display()
        ))
    })
}

fn optional_kind_from_path(path: &Path) -> Option<IndexKind> {
    let extension = path.extension()?.to_str()?;
    if extension.eq_ignore_ascii_case("csi") {
        Some(IndexKind::Csi)
    } else if extension.eq_ignore_ascii_case("tbi") {
        Some(IndexKind::Tbi)
    } else {
        None
    }
}

fn read_variant_header(path: &Path) -> Result<vcf::Header> {
    let mut reader = crate::format::Reader::open(path)?;
    reader.read_header().map(|(header, _, _)| header)
}

fn contig_lengths(header: &vcf::Header) -> HashMap<String, Option<usize>> {
    header
        .contigs()
        .iter()
        .map(|(name, contig)| (name.clone(), contig.length()))
        .collect()
}

fn metadata(index: &Index) -> (Vec<String>, Vec<Option<u64>>) {
    match index {
        Index::Csi(index) => {
            let names = index
                .header()
                .map(|header| {
                    header
                        .reference_sequence_names()
                        .iter()
                        .map(|name| String::from_utf8_lossy(name).into_owned())
                        .collect()
                })
                .unwrap_or_default();
            let counts = index
                .reference_sequences()
                .iter()
                .map(mapped_records)
                .collect();
            (names, counts)
        }
        Index::Tbi(index) => {
            let names = index
                .header()
                .map(|header| {
                    header
                        .reference_sequence_names()
                        .iter()
                        .map(|name| String::from_utf8_lossy(name).into_owned())
                        .collect()
                })
                .unwrap_or_default();
            let counts = index
                .reference_sequences()
                .iter()
                .map(mapped_records)
                .collect();
            (names, counts)
        }
    }
}

fn mapped_records<I>(reference: &ReferenceSequence<I>) -> Option<u64>
where
    I: Default,
    ReferenceSequence<I>: binning_index::ReferenceSequence,
{
    use binning_index::ReferenceSequence as _;
    reference.metadata().map(Metadata::mapped_record_count)
}
