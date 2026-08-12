use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use noodles_vcf::variant::RecordBuf;
use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Direction {
    Left,
    Right,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Strand {
    Forward,
    Reverse,
    Unknown,
}

#[derive(Clone, Debug)]
struct Transcript {
    start: usize,
    end: usize,
    strand: Strand,
}

#[derive(Debug)]
struct TranscriptNode {
    transcript: Transcript,
    max_end: usize,
    left: Option<Box<Self>>,
    right: Option<Box<Self>>,
}

impl TranscriptNode {
    fn build(transcripts: &[Transcript]) -> Option<Box<Self>> {
        if transcripts.is_empty() {
            return None;
        }
        let middle = transcripts.len() / 2;
        let left = Self::build(&transcripts[..middle]);
        let right = Self::build(&transcripts[middle + 1..]);
        let max_end = left
            .as_ref()
            .map_or(transcripts[middle].end, |node| node.max_end)
            .max(right.as_ref().map_or(0, |node| node.max_end));
        Some(Box::new(Self {
            transcript: transcripts[middle].clone(),
            max_end,
            left,
            right,
        }))
    }

    fn has_reverse(&self, start: usize, end: usize, forward: &mut bool) -> bool {
        if self.max_end < start {
            return false;
        }
        if let Some(left) = &self.left
            && left.max_end >= start
            && left.has_reverse(start, end, forward)
        {
            return true;
        }
        if self.transcript.start <= end && self.transcript.end >= start {
            match self.transcript.strand {
                Strand::Forward => *forward = true,
                Strand::Reverse => return true,
                Strand::Unknown => {}
            }
        }
        if self.transcript.start <= end
            && let Some(right) = &self.right
            && right.max_end >= start
            && right.has_reverse(start, end, forward)
        {
            return true;
        }
        false
    }
}

#[derive(Debug)]
struct Gene {
    sequence: String,
}

#[derive(Debug)]
struct PendingTranscript {
    line: usize,
    sequence: String,
    parent: String,
    start: usize,
    end: usize,
    strand: Strand,
}

#[derive(Debug)]
pub(super) struct TranscriptIndex {
    sequences: HashMap<String, Box<TranscriptNode>>,
}

impl TranscriptIndex {
    pub(super) fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .rs_with_context(|| format!("opening GFF3 annotation {}", path.display()))?;
        let mut source = BufReader::new(file);
        let compressed = source
            .fill_buf()
            .rs_with_context(|| format!("reading GFF3 annotation {}", path.display()))?
            .get(..2)
            .is_some_and(|magic| magic == [0x1f, 0x8b]);
        let reader: Box<dyn BufRead> = if compressed {
            Box::new(BufReader::new(MultiGzDecoder::new(source)))
        } else {
            Box::new(source)
        };
        Self::read(path, reader)
    }

    fn read(path: &Path, mut reader: Box<dyn BufRead>) -> Result<Self> {
        let mut genes = HashMap::new();
        let mut transcripts = Vec::new();
        let mut transcript_ids = HashSet::new();
        let mut text = String::new();
        let mut line_number = 0;
        loop {
            text.clear();
            if reader
                .read_line(&mut text)
                .rs_with_context(|| format!("reading GFF3 annotation {}", path.display()))?
                == 0
            {
                break;
            }
            line_number += 1;
            let line = text.trim_end_matches(['\r', '\n']);
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 9 {
                return Err(invalid(
                    path,
                    line_number,
                    "expected nine tab-separated fields",
                ));
            }
            let start = fields[3]
                .parse::<usize>()
                .map_err(|_| invalid(path, line_number, "start is not a positive integer"))?;
            let end = fields[4]
                .parse::<usize>()
                .map_err(|_| invalid(path, line_number, "end is not a positive integer"))?;
            if start == 0 || end < start {
                return Err(invalid(path, line_number, "invalid one-based interval"));
            }
            let strand = match fields[6] {
                "+" => Strand::Forward,
                "-" => Strand::Reverse,
                "." | "?" => Strand::Unknown,
                _ => return Err(invalid(path, line_number, "strand must be +, -, ., or ?")),
            };
            let attributes = attributes(path, line_number, fields[8])?;
            if fields[2] == "gene"
                || attributes
                    .get("ID")
                    .is_some_and(|id| id.starts_with("gene:"))
            {
                let Some(id) = attributes.get("ID") else {
                    continue;
                };
                let id = strip_prefix(id, "gene:").to_owned();
                if genes
                    .insert(
                        id,
                        Gene {
                            sequence: fields[0].to_owned(),
                        },
                    )
                    .is_some()
                {
                    return Err(invalid(path, line_number, "duplicate gene ID"));
                }
                continue;
            }
            if matches!(
                fields[2],
                "CDS" | "exon" | "three_prime_UTR" | "five_prime_UTR"
            ) {
                continue;
            }
            let supported = attributes
                .get("biotype")
                .or_else(|| attributes.get("gene_biotype"))
                .is_some_and(|value| is_supported(value))
                || is_supported(fields[2]);
            if !supported {
                continue;
            }
            let id = attributes
                .get("ID")
                .ok_or_else(|| invalid(path, line_number, "transcript is missing ID"))?;
            let id = strip_prefix(id, "transcript:").to_owned();
            if !transcript_ids.insert(id) {
                return Err(invalid(path, line_number, "duplicate transcript ID"));
            }
            let parent = attributes
                .get("Parent")
                .ok_or_else(|| invalid(path, line_number, "transcript is missing Parent"))?;
            transcripts.push(PendingTranscript {
                line: line_number,
                sequence: fields[0].to_owned(),
                parent: strip_prefix(parent, "gene:").to_owned(),
                start: start - 1,
                end: end - 1,
                strand,
            });
        }

        let mut sequences: HashMap<String, Vec<Transcript>> = HashMap::new();
        for transcript in transcripts {
            let gene = genes.get(&transcript.parent).ok_or_else(|| {
                invalid(
                    path,
                    transcript.line,
                    &format!("transcript parent {} is absent", transcript.parent),
                )
            })?;
            if gene.sequence != transcript.sequence {
                return Err(invalid(
                    path,
                    transcript.line,
                    "transcript and parent gene use different sequences",
                ));
            }
            sequences
                .entry(gene.sequence.clone())
                .or_default()
                .push(Transcript {
                    start: transcript.start,
                    end: transcript.end,
                    strand: transcript.strand,
                });
        }
        if sequences.values().all(Vec::is_empty) {
            return Err(RsomicsError::InvalidInput(format!(
                "GFF3 annotation {} contains no usable transcripts",
                path.display()
            )));
        }
        for transcripts in sequences.values_mut() {
            transcripts.sort_unstable_by_key(|transcript| transcript.start);
        }
        Ok(Self {
            sequences: sequences
                .into_iter()
                .map(|(sequence, transcripts)| {
                    (sequence, TranscriptNode::build(&transcripts).unwrap())
                })
                .collect(),
        })
    }

    pub(super) fn direction(&self, record: &RecordBuf) -> Direction {
        let Some(start) = record
            .variant_start()
            .map(|position| usize::from(position) - 1)
        else {
            return Direction::Left;
        };
        let end = start.saturating_add(record.reference_bases().len());
        let Some(transcripts) = self.sequences.get(record.reference_sequence_name()) else {
            return Direction::Left;
        };
        let mut forward = false;
        if transcripts.has_reverse(start, end, &mut forward) {
            return Direction::Left;
        }
        if forward {
            Direction::Right
        } else {
            Direction::Left
        }
    }
}

fn attributes<'a>(path: &Path, line: usize, field: &'a str) -> Result<HashMap<&'a str, &'a str>> {
    let mut attributes = HashMap::new();
    for item in field.split(';').filter(|item| !item.is_empty()) {
        let (key, value) = item
            .split_once('=')
            .ok_or_else(|| invalid(path, line, "attribute is missing ="))?;
        if key.is_empty() || value.is_empty() {
            return Err(invalid(
                path,
                line,
                "attribute key and value must be non-empty",
            ));
        }
        if attributes.insert(key, value).is_some() {
            return Err(invalid(path, line, "duplicate attribute key"));
        }
    }
    Ok(attributes)
}

fn strip_prefix<'a>(value: &'a str, prefix: &str) -> &'a str {
    value.strip_prefix(prefix).unwrap_or(value)
}

fn is_supported(value: &str) -> bool {
    if value.eq_ignore_ascii_case("mRNA")
        || value.eq_ignore_ascii_case("Mt_rRNA")
        || value.eq_ignore_ascii_case("Mt_tRNA")
        || value.eq_ignore_ascii_case("3prime_overlapping_ncRNA")
        || value.eq_ignore_ascii_case("3_prime_overlapping_ncRNA")
    {
        return true;
    }
    matches!(
        value,
        "protein_coding"
            | "polymorphic_pseudogene"
            | "IG_C"
            | "IG_D"
            | "IG_J"
            | "IG_LV"
            | "IG_V"
            | "TR_C"
            | "TR_D"
            | "TR_J"
            | "TR_V"
            | "NMD"
            | "nonsense_mediated_decay"
            | "non_stop_decay"
            | "lincRNA"
            | "lncRNA"
            | "miRNA"
            | "misc_RNA"
            | "rRNA"
            | "snRNA"
            | "snoRNA"
            | "processed_transcript"
            | "antisense"
            | "macro_lncRNA"
            | "ribozyme"
            | "sRNA"
            | "scRNA"
            | "scaRNA"
            | "sense_intronic"
            | "sense_overlapping"
            | "pseudogene"
            | "processed_pseudogene"
            | "artifact"
            | "IG_pseudogene"
            | "IG_C_pseudogene"
            | "IG_J_pseudogene"
            | "IG_V_pseudogene"
            | "TR_V_pseudogene"
            | "TR_J_pseudogene"
            | "MT_tRNA_pseudogene"
            | "misc_RNA_pseudogene"
            | "miRNA_pseudogene"
            | "retained_intron"
            | "retrotransposed"
            | "tRNA_pseudogene"
            | "Trna_pseudogene"
            | "transcribed_processed_pseudogene"
            | "transcribed_unprocessed_pseudogene"
            | "transcribed_unitary_pseudogene"
            | "translated_unprocessed_pseudogene"
            | "translated_processed_pseudogene"
            | "known_ncRNA"
            | "known_ncrna"
            | "unitary_pseudogene"
            | "unprocessed_pseudogene"
            | "LRG_gene"
            | "disrupted_domain"
            | "vaultRNA"
            | "bidirectional_promoter_lncRNA"
            | "ambiguous_orf"
    )
}

fn invalid(path: &Path, line: usize, message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(format!("GFF3 {} line {line}: {message}", path.display()))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use noodles_vcf as vcf;

    use super::*;

    fn record(line: &[u8]) -> vcf::variant::RecordBuf {
        let header: vcf::Header = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=300>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .parse()
            .unwrap();
        let raw = vcf::Record::try_from(line).unwrap();
        vcf::variant::RecordBuf::try_from_variant_record(&header, &raw).unwrap()
    }

    #[test]
    fn forward_only_overlap_selects_right_alignment() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("annotation.gff3");
        fs::write(
            &path,
            b"##gff-version 3\n\
chr1\t.\tgene\t1\t300\t.\t+\t.\tID=gene:g1;biotype=protein_coding\n\
chr1\t.\tmRNA\t10\t200\t.\t+\t.\tID=transcript:t1;Parent=gene:g1;biotype=protein_coding\n",
        )
        .unwrap();

        let index = TranscriptIndex::open(&path).unwrap();
        assert_eq!(
            index.direction(&record(b"chr1\t20\t.\tAA\tA\t.\tPASS\t.")),
            Direction::Right
        );
        assert_eq!(
            index.direction(&record(b"chr1\t250\t.\tAA\tA\t.\tPASS\t.")),
            Direction::Left
        );
    }

    #[test]
    fn reverse_or_conflicting_overlap_selects_left_alignment() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("annotation.gff3");
        fs::write(
            &path,
            b"##gff-version 3\n\
chr1\t.\tgene\t1\t300\t.\t+\t.\tID=g1;biotype=protein_coding\n\
chr1\t.\tmRNA\t10\t100\t.\t+\t.\tID=t1;Parent=g1;biotype=protein_coding\n\
chr1\t.\tgene\t1\t300\t.\t-\t.\tID=g2;biotype=antisense\n\
chr1\t.\ttranscript\t50\t200\t.\t-\t.\tID=t2;Parent=g2;biotype=antisense\n",
        )
        .unwrap();

        let index = TranscriptIndex::open(&path).unwrap();
        assert_eq!(
            index.direction(&record(b"chr1\t20\t.\tAA\tA\t.\tPASS\t.")),
            Direction::Right
        );
        assert_eq!(
            index.direction(&record(b"chr1\t60\t.\tAA\tA\t.\tPASS\t.")),
            Direction::Left
        );
        assert_eq!(
            index.direction(&record(b"chr1\t150\t.\tAA\tA\t.\tPASS\t.")),
            Direction::Left
        );
    }

    #[test]
    fn malformed_or_empty_annotations_fail() {
        let directory = tempfile::tempdir().unwrap();
        let malformed = directory.path().join("malformed.gff3");
        fs::write(
            &malformed,
            b"chr1\t.\tmRNA\t10\t5\t.\t+\t.\tID=t1;Parent=g1\n",
        )
        .unwrap();
        assert!(TranscriptIndex::open(&malformed).is_err());

        let empty = directory.path().join("empty.gff3");
        fs::write(&empty, b"##gff-version 3\n").unwrap();
        let error = TranscriptIndex::open(&empty).unwrap_err().to_string();
        assert!(error.contains("no usable transcripts"), "{error}");
    }
}
