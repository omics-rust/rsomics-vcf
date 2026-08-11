use std::collections::{HashSet, VecDeque};

use noodles_vcf::{
    self as vcf,
    header::record::value::{Map, map::Filter as HeaderFilter},
    variant::{
        RecordBuf,
        record_buf::{
            info::field::{Value as InfoValue, value::Array},
            samples::sample::Value as SampleValue,
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use crate::{filter::Disposition, variant_type};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SnpGap {
    pub distance: usize,
    pub types: u32,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct Options {
    pub snp_gap: Option<SnpGap>,
    pub indel_gap: Option<usize>,
    pub soft: bool,
}

pub(crate) struct GapOutput {
    record: RecordBuf,
    disposition: Disposition,
}

impl GapOutput {
    pub(crate) fn record(&self) -> &RecordBuf {
        &self.record
    }

    pub(crate) fn disposition(&self) -> Disposition {
        self.disposition
    }
}

struct Entry {
    id: u64,
    mask: u32,
    end: usize,
    failed: bool,
    record: RecordBuf,
}

struct IndelCluster {
    until: usize,
    members: Vec<u64>,
}

pub(crate) struct GapBuffer {
    options: Options,
    entries: VecDeque<Entry>,
    cluster: Option<IndelCluster>,
    next_id: u64,
    reference: Option<String>,
    position: Option<usize>,
    completed_references: HashSet<String>,
}

impl GapBuffer {
    pub(crate) fn new(header: &mut vcf::Header, options: Options) -> Result<Self> {
        if options.snp_gap.is_some_and(|gap| gap.types == 0) {
            return Err(RsomicsError::ConfigError(
                "SnpGap requires at least one target variant type".to_owned(),
            ));
        }
        if let Some(gap) = options.snp_gap {
            header
                .filters_mut()
                .entry("SnpGap".to_owned())
                .or_insert_with(|| {
                    Map::<HeaderFilter>::new(format!(
                        "SNP within {} bp of {}",
                        gap.distance,
                        type_names(gap.types)
                    ))
                });
        }
        if let Some(gap) = options.indel_gap {
            header
                .filters_mut()
                .entry("IndelGap".to_owned())
                .or_insert_with(|| {
                    Map::<HeaderFilter>::new(format!("Indel within {gap} bp of an indel"))
                });
        }
        Ok(Self {
            options,
            entries: VecDeque::new(),
            cluster: None,
            next_id: 0,
            reference: None,
            position: None,
            completed_references: HashSet::new(),
        })
    }

    pub(crate) fn push(&mut self, record: RecordBuf) -> Result<Vec<GapOutput>> {
        let reference = record.reference_sequence_name();
        let position = record
            .variant_start()
            .map(usize::from)
            .ok_or_else(|| RsomicsError::InvalidInput("gap filters require POS".to_owned()))?;
        let mut output = Vec::new();
        if self
            .reference
            .as_deref()
            .is_some_and(|name| name != reference)
        {
            if self.completed_references.contains(reference) {
                return Err(RsomicsError::InvalidInput(format!(
                    "reference sequence {reference} is not contiguous"
                )));
            }
            self.finalize_cluster()?;
            output.extend(self.drain_all());
            self.completed_references
                .insert(self.reference.take().unwrap());
            self.position = None;
        }
        if self.reference.is_none() {
            self.reference = Some(reference.to_owned());
        }
        if self.position.is_some_and(|previous| position < previous) {
            return Err(RsomicsError::InvalidInput(format!(
                "gap filter input is not sorted at {reference}:{position}"
            )));
        }
        self.position = Some(position);

        if self
            .cluster
            .as_ref()
            .is_some_and(|cluster| position > cluster.until)
        {
            self.finalize_cluster()?;
        }

        let mask = variant_type::record_mask(&record);
        let end = position
            .checked_add(record.reference_bases().len().max(1) - 1)
            .ok_or_else(|| RsomicsError::InvalidInput("variant end overflows usize".to_owned()))?;
        let id = self.next_id;
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or_else(|| RsomicsError::InvalidInput("variant count exceeds u64".to_owned()))?;
        self.entries.push_back(Entry {
            id,
            mask,
            end,
            failed: false,
            record,
        });

        if let Some(distance) = self
            .options
            .indel_gap
            .filter(|_| mask & variant_type::INDEL != 0)
        {
            let until = end.checked_add(distance).ok_or_else(|| {
                RsomicsError::InvalidInput("IndelGap coordinate overflows usize".to_owned())
            })?;
            match &mut self.cluster {
                Some(cluster) => {
                    cluster.until = until;
                    cluster.members.push(id);
                }
                None => {
                    self.cluster = Some(IndelCluster {
                        until,
                        members: vec![id],
                    });
                }
            }
        }

        self.apply_snp_gap(position, mask)?;
        output.extend(self.drain_ready(position)?);
        Ok(output)
    }

    pub(crate) fn finish(&mut self) -> Result<Vec<GapOutput>> {
        self.finalize_cluster()?;
        self.reference = None;
        self.position = None;
        Ok(self.drain_all())
    }

    #[cfg(test)]
    pub(crate) fn buffered_len(&self) -> usize {
        self.entries.len()
    }

    fn apply_snp_gap(&mut self, position: usize, mask: u32) -> Result<()> {
        let Some(gap) = self.options.snp_gap else {
            return Ok(());
        };
        let mut mark_new = false;
        for entry in &mut self.entries {
            let in_range = entry.end.checked_add(gap.distance).ok_or_else(|| {
                RsomicsError::InvalidInput("SnpGap coordinate overflows usize".to_owned())
            })? >= position;
            if !in_range {
                continue;
            }
            if mask & gap.types != 0 && entry.mask & variant_type::SNP != 0 {
                mark(entry, "SnpGap");
            } else if mask & variant_type::SNP != 0 && entry.mask & gap.types != 0 {
                mark_new = true;
            }
        }
        if mark_new {
            let entry = self.entries.back_mut().unwrap();
            mark(entry, "SnpGap");
        }
        Ok(())
    }

    fn finalize_cluster(&mut self) -> Result<()> {
        let Some(cluster) = self.cluster.take() else {
            return Ok(());
        };
        let members: Vec<_> = cluster
            .members
            .iter()
            .map(|id| {
                self.entries
                    .iter()
                    .position(|entry| entry.id == *id)
                    .ok_or_else(|| {
                        RsomicsError::InvalidInput(
                            "IndelGap state lost a pending record".to_owned(),
                        )
                    })
            })
            .collect::<Result<_>>()?;
        let first = members[0];
        let mut best_quality = first;
        let mut maximum_quality = self.entries[first]
            .record
            .quality_score()
            .unwrap_or(f32::NAN);
        for &index in &members[1..] {
            let quality = self.entries[index]
                .record
                .quality_score()
                .unwrap_or(f32::NAN);
            if maximum_quality < quality {
                maximum_quality = quality;
                best_quality = index;
            }
        }
        let best = if maximum_quality > 0.0 {
            best_quality
        } else {
            let mut best = first;
            let mut maximum_ac = first_alt_count(&self.entries[first].record)?.unwrap_or(0);
            for &index in &members[1..] {
                if let Some(ac) = first_alt_count(&self.entries[index].record)?
                    && maximum_ac < ac
                {
                    maximum_ac = ac;
                    best = index;
                }
            }
            best
        };
        for index in members {
            if index != best {
                mark(&mut self.entries[index], "IndelGap");
            }
        }
        Ok(())
    }

    fn drain_ready(&mut self, position: usize) -> Result<Vec<GapOutput>> {
        let mut output = Vec::new();
        loop {
            let Some(entry) = self.entries.front() else {
                break;
            };
            let snp_safe = match self.options.snp_gap {
                Some(gap) => {
                    entry.end.checked_add(gap.distance).ok_or_else(|| {
                        RsomicsError::InvalidInput("SnpGap coordinate overflows usize".to_owned())
                    })? < position
                }
                None => true,
            };
            let indel_safe = self
                .cluster
                .as_ref()
                .is_none_or(|cluster| entry.id < cluster.members[0]);
            if !snp_safe || !indel_safe {
                break;
            }
            output.push(self.pop_front());
        }
        Ok(output)
    }

    fn drain_all(&mut self) -> Vec<GapOutput> {
        let mut output = Vec::with_capacity(self.entries.len());
        while !self.entries.is_empty() {
            output.push(self.pop_front());
        }
        output
    }

    fn pop_front(&mut self) -> GapOutput {
        let entry = self.entries.pop_front().unwrap();
        let disposition = if entry.failed && !self.options.soft {
            Disposition::Drop
        } else {
            Disposition::Keep
        };
        GapOutput {
            record: entry.record,
            disposition,
        }
    }
}

fn mark(entry: &mut Entry, filter: &str) {
    entry.failed = true;
    let filters = entry.record.filters_mut().as_mut();
    filters.shift_remove("PASS");
    filters.insert(filter.to_owned());
}

fn first_alt_count(record: &RecordBuf) -> Result<Option<i32>> {
    if let (Some(Some(InfoValue::Integer(an))), Some(Some(InfoValue::Array(Array::Integer(ac))))) =
        (record.info().get("AN"), record.info().get("AC"))
        && *an >= 0
    {
        if ac.len() != record.alternate_bases().as_ref().len() {
            return Ok(None);
        }
        let Some(values) = ac.iter().copied().collect::<Option<Vec<_>>>() else {
            return Ok(None);
        };
        let total = values.iter().try_fold(0i32, |total, value| {
            total
                .checked_add(*value)
                .ok_or_else(|| RsomicsError::InvalidInput("INFO/AC count overflow".to_owned()))
        })?;
        if total > *an {
            return Err(RsomicsError::InvalidInput(
                "INFO/AC exceeds INFO/AN".to_owned(),
            ));
        }
        return Ok(values.first().copied());
    }

    let Some(genotypes) = record.samples().select("GT") else {
        return Ok(None);
    };
    let mut count = 0i32;
    for sample in 0..record.samples().values().count() {
        let Some(Some(value)) = genotypes.get(sample) else {
            continue;
        };
        let SampleValue::Genotype(genotype) = value else {
            return Err(RsomicsError::InvalidInput(
                "FORMAT/GT is not encoded as a genotype".to_owned(),
            ));
        };
        for allele in genotype.as_ref() {
            if let Some(position) = allele.position()
                && position > record.alternate_bases().as_ref().len()
            {
                return Err(RsomicsError::InvalidInput(format!(
                    "genotype allele index {position} exceeds {} ALT alleles",
                    record.alternate_bases().as_ref().len()
                )));
            }
            if allele.position() == Some(1) {
                count = count.checked_add(1).ok_or_else(|| {
                    RsomicsError::InvalidInput("allele count exceeds int32".to_owned())
                })?;
            }
        }
    }
    Ok(Some(count))
}

fn type_names(mask: u32) -> String {
    let mut names = Vec::new();
    for (kind, name) in [
        (variant_type::INDEL, "indel"),
        (variant_type::MNP, "mnp"),
        (variant_type::BND, "bnd"),
        (variant_type::OTHER, "other"),
        (variant_type::OVERLAP, "overlap"),
    ] {
        if mask & kind != 0 {
            names.push(name);
        }
    }
    names.join(",")
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use super::*;
    use crate::{filter::Disposition, variant_type};

    fn header() -> vcf::Header {
        "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"allele number\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap()
    }

    fn record(header: &vcf::Header, line: &str) -> RecordBuf {
        let record = vcf::Record::try_from(line.as_bytes()).unwrap();
        RecordBuf::try_from_variant_record(header, &record).unwrap()
    }

    fn collect(
        buffer: &mut GapBuffer,
        records: impl IntoIterator<Item = RecordBuf>,
    ) -> Vec<GapOutput> {
        let mut output = Vec::new();
        for record in records {
            output.extend(buffer.push(record).unwrap());
        }
        output.extend(buffer.finish().unwrap());
        output
    }

    fn positions(records: &[GapOutput], disposition: Disposition) -> Vec<usize> {
        records
            .iter()
            .filter(|record| record.disposition() == disposition)
            .map(|record| usize::from(record.record().variant_start().unwrap()))
            .collect()
    }

    #[test]
    fn snp_gap_marks_snps_on_both_sides_of_indels() {
        let mut header = header();
        let records = [
            "chr1\t1\ts1\tA\tG\t10\tPASS\t.\tGT\t0/1",
            "chr1\t2\td1\tAT\tA\t10\tPASS\t.\tGT\t0/1",
            "chr1\t7\ts2\tA\tC\t10\tPASS\t.\tGT\t0/1",
            "chr1\t8\ti1\tA\tAT\t10\tPASS\t.\tGT\t0/1",
            "chr1\t12\ts3\tA\tT\t10\tPASS\t.\tGT\t0/1",
        ]
        .map(|line| record(&header, line));
        let mut buffer = GapBuffer::new(
            &mut header,
            Options {
                snp_gap: Some(SnpGap {
                    distance: 3,
                    types: variant_type::INDEL,
                }),
                ..Options::default()
            },
        )
        .unwrap();
        let output = collect(&mut buffer, records);
        assert_eq!(positions(&output, Disposition::Drop), [1, 7]);
        assert!(header.filters().contains_key("SnpGap"));
    }

    #[test]
    fn indel_gap_prefers_positive_quality_then_first_alt_count_then_input_order() {
        let mut header = header();
        let records = [
            "chr1\t1\tq1\tA\tAT\t10\tPASS\tAC=5;AN=10\tGT\t0/1",
            "chr1\t3\tq2\tA\tAT\t20\tPASS\tAC=1;AN=10\tGT\t0/1",
            "chr1\t7\ta1\tA\tAT\t.\tPASS\tAC=5;AN=10\tGT\t0/1",
            "chr1\t9\ta2\tA\tAT\t.\tPASS\tAC=2;AN=10\tGT\t0/1",
            "chr1\t13\tf1\tA\tAT\t.\tPASS\t.\tGT\t0/1",
            "chr1\t15\tf2\tA\tAT\t.\tPASS\t.\tGT\t0/1",
        ]
        .map(|line| record(&header, line));
        let mut buffer = GapBuffer::new(
            &mut header,
            Options {
                indel_gap: Some(2),
                ..Options::default()
            },
        )
        .unwrap();
        let output = collect(&mut buffer, records);
        assert_eq!(positions(&output, Disposition::Drop), [1, 9, 15]);
        assert!(header.filters().contains_key("IndelGap"));
    }

    #[test]
    fn soft_gap_filters_annotate_without_dropping_and_state_is_bounded() {
        let mut header = header();
        let mut buffer = GapBuffer::new(
            &mut header,
            Options {
                snp_gap: Some(SnpGap {
                    distance: 2,
                    types: variant_type::INDEL,
                }),
                indel_gap: Some(2),
                soft: true,
            },
        )
        .unwrap();
        let mut output = Vec::new();
        for position in (1..=100).step_by(10) {
            output.extend(
                buffer
                    .push(record(
                        &header,
                        &format!("chr1\t{position}\t.\tA\tG\t10\tPASS\t.\tGT\t0/1"),
                    ))
                    .unwrap(),
            );
            assert!(buffer.buffered_len() <= 2);
        }
        output.extend(buffer.finish().unwrap());
        assert!(
            output
                .iter()
                .all(|record| record.disposition() == Disposition::Keep)
        );
    }

    #[test]
    fn chromosome_changes_flush_pending_clusters() {
        let mut header = header();
        let mut buffer = GapBuffer::new(
            &mut header,
            Options {
                indel_gap: Some(5),
                ..Options::default()
            },
        )
        .unwrap();
        assert!(
            buffer
                .push(record(&header, "chr1\t10\ta\tA\tAT\t10\tPASS\t.\tGT\t0/1"))
                .unwrap()
                .is_empty()
        );
        let flushed = buffer
            .push(record(&header, "chr2\t1\tb\tA\tAT\t10\tPASS\t.\tGT\t0/1"))
            .unwrap();
        assert_eq!(positions(&flushed, Disposition::Keep), [10]);
    }
}
