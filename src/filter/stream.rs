use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::expression::Compiled;
use crate::format::{
    HeaderMode, OutputFormat, ParallelWriter, Reader, RecordScratch, VariantWriter, Writer,
    reformat_record,
};
use crate::regions::{IndexedRecords, RegionSelection, RegionSet};

use super::{Options, Processor, gaps};

#[derive(Clone, Debug, Default)]
pub(crate) struct StreamOptions {
    pub output_format: OutputFormat,
    pub expression: Option<String>,
    pub filter: Options,
    pub gaps: gaps::Options,
    pub regions: Option<RegionSet>,
    pub targets: Option<RegionSelection>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct Summary {
    pub read: u64,
    pub written: u64,
    pub output_format: OutputFormat,
}

pub(crate) fn write(input: &Path, options: &StreamOptions, output: impl Write) -> Result<Summary> {
    let writer = Writer::new(output, options.output_format);
    write_with_writer(input, options, writer)
}

pub(crate) fn write_parallel<W>(
    input: &Path,
    options: &StreamOptions,
    output: W,
    workers: usize,
) -> Result<Summary>
where
    W: Write + Send + 'static,
{
    let writer = ParallelWriter::new(output, options.output_format, workers)?;
    write_with_writer(input, options, writer)
}

fn write_with_writer(
    input: &Path,
    options: &StreamOptions,
    mut writer: impl VariantWriter,
) -> Result<Summary> {
    if let Some(regions) = &options.regions {
        if input == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "indexed regions require a named input".to_owned(),
            ));
        }
        return write_indexed(input, regions, options, &mut writer);
    }

    let mut reader = Reader::open(input)?;
    let (mut header, _, schema) = reader.read_header()?;
    if reader.is_text()
        && writer.supports_vcf_records()
        && options.expression.is_some()
        && options.filter.mask.is_none()
        && options.filter.soft_filter.is_none()
        && options.filter.set_genotypes.is_none()
        && options.gaps.snp_gap.is_none()
        && options.gaps.indel_gap.is_none()
        && options.targets.is_none()
    {
        let expression = options.expression.as_deref().unwrap();
        let compiled = Compiled::bind(expression, &header).map_err(|error| {
            RsomicsError::ConfigError(format!("invalid filter expression: {error}"))
        })?;
        if let Some(predicate) = compiled.raw() {
            return write_raw_text(
                input,
                &mut reader,
                &header,
                &schema,
                options,
                predicate,
                &mut writer,
            );
        }
    }
    let mut processor = Processor::bind(
        &mut header,
        options.expression.as_deref(),
        options.filter.clone(),
        options.gaps,
    )?;
    writer.write_header(&header, HeaderMode::Full)?;

    let mut scratch = RecordScratch::default();
    let mut read = 0;
    let mut written = 0;
    {
        let mut emit = |record| {
            written += 1;
            writer.write_record(&header, &record, written)
        };
        loop {
            let number = read + 1;
            let Some(record) = reader.read_record(&header, &mut scratch, number)? else {
                break;
            };
            read += 1;
            if options
                .targets
                .as_ref()
                .is_some_and(|targets| !targets.keeps(&record))
            {
                continue;
            }
            processor.push(&header, record, &mut emit)?;
        }
        processor.finish(&mut emit)?;
    }
    writer.finish()?;
    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}

fn write_raw_text(
    input: &Path,
    reader: &mut Reader,
    header: &noodles_vcf::Header,
    schema: &crate::format::HeaderTypes,
    options: &StreamOptions,
    predicate: crate::expression::RawPredicate,
    writer: &mut impl VariantWriter,
) -> Result<Summary> {
    writer.write_header(header, HeaderMode::Full)?;
    let mut raw = Vec::new();
    let mut canonical = Vec::new();
    let mut read = 0;
    let mut written = 0;
    loop {
        let number = read + 1;
        if reader.read_text_record(
            &mut raw,
            usize::try_from(number).map_err(|_| {
                RsomicsError::InvalidInput("variant record count exceeds usize".to_owned())
            })?,
        )? == 0
        {
            break;
        }
        read += 1;
        reformat_record(&raw, schema, &mut canonical).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "{}: parsing variant record {number}: {error}",
                input.display()
            ))
        })?;
        normalize_missing_filter(&mut canonical);
        let passes = predicate.evaluate(&canonical).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "{}: evaluating filter expression at record {number}: scalar numeric field has multiple or invalid values",
                input.display()
            ))
        })?;
        if options.filter.logic.accepts(passes) {
            written += 1;
            writer.write_vcf_record(&canonical, written)?;
        }
    }
    writer.finish()?;
    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}

fn normalize_missing_filter(record: &mut Vec<u8>) {
    let mut tabs = record
        .iter()
        .enumerate()
        .filter(|(_, byte)| **byte == b'\t');
    let Some(start) = tabs.nth(5).map(|(index, _)| index + 1) else {
        return;
    };
    let Some(end) = tabs.next().map(|(index, _)| index) else {
        return;
    };
    if &record[start..end] == b"." {
        record.splice(start..end, *b"PASS");
    }
}

fn write_indexed(
    input: &Path,
    regions: &RegionSet,
    options: &StreamOptions,
    writer: &mut impl VariantWriter,
) -> Result<Summary> {
    let mut reader = IndexedRecords::open(input, regions)?;
    let mut header = reader.header().clone();
    let mut processor = Processor::bind(
        &mut header,
        options.expression.as_deref(),
        options.filter.clone(),
        options.gaps,
    )?;
    writer.write_header(&header, HeaderMode::Full)?;

    let mut written = 0;
    let read = {
        let mut emit = |record| {
            written += 1;
            writer.write_record(&header, &record, written)
        };
        let read = reader.visit(|_, record, _| {
            if options
                .targets
                .as_ref()
                .is_some_and(|targets| !targets.keeps(&record))
            {
                return Ok(());
            }
            processor.push(&header, record, &mut emit)
        })?;
        processor.finish(&mut emit)?;
        read
    };
    writer.finish()?;
    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::Read;

    use crate::format::OutputFormat;
    use crate::index::{self, BuildOptions, IndexKind};
    use crate::regions::{OverlapMode, RegionSelection, RegionSet};
    use crate::variant_type;
    use noodles_bgzf as bgzf;
    use rsomics_common::AtomicFile;

    use super::*;

    #[test]
    fn targets_are_selected_before_expression_and_gap_filters() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("calls.vcf");
        fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tC\t10\tPASS\t.\n\
chr1\t12\t.\tA\tAT\t1\tPASS\t.\n\
chr1\t20\t.\tG\tT\t10\tPASS\t.\n",
        )
        .unwrap();
        let options = StreamOptions {
            expression: Some("QUAL >= 10".to_owned()),
            gaps: gaps::Options {
                snp_gap: Some(gaps::SnpGap {
                    distance: 3,
                    types: variant_type::INDEL,
                }),
                ..gaps::Options::default()
            },
            targets: Some(
                RegionSelection::parse(["chr1:12-12".to_owned()], OverlapMode::Position, true)
                    .unwrap(),
            ),
            ..StreamOptions::default()
        };
        let mut output = Vec::new();
        let summary = write(&input, &options, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        let records: Vec<_> = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect();

        assert_eq!(summary.read, 3);
        assert_eq!(summary.written, 2, "{output}");
        assert_eq!(summary.output_format, OutputFormat::Vcf);
        assert_eq!(records.len(), 2);
        assert!(records[0].starts_with("chr1\t10\t"));
        assert!(records[1].starts_with("chr1\t20\t"));
        assert!(output.contains("##FILTER=<ID=SnpGap"));
    }

    #[test]
    fn indexed_regions_deduplicate_before_target_filtering() {
        let directory = tempfile::tempdir().unwrap();
        let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/index.vcf");
        let input = directory.path().join("calls.vcf.gz");
        let mut writer = bgzf::io::Writer::new(File::create(&input).unwrap());
        writer.write_all(&fs::read(source).unwrap()).unwrap();
        writer.try_finish().unwrap();
        index::create(
            &input,
            &index::default_output_path(&input, IndexKind::Tbi),
            BuildOptions {
                kind: IndexKind::Tbi,
                min_shift: 14,
                threads: 1,
                force: false,
            },
        )
        .unwrap();

        let options = StreamOptions {
            expression: Some("QUAL >= 60".to_owned()),
            regions: Some(
                RegionSet::parse(
                    ["chr1:70000-70000".to_owned(), "chr1:69990-70010".to_owned()],
                    OverlapMode::Record,
                )
                .unwrap(),
            ),
            targets: Some(
                RegionSelection::parse(
                    ["chr1:70000-70000".to_owned()],
                    OverlapMode::Position,
                    true,
                )
                .unwrap(),
            ),
            ..StreamOptions::default()
        };
        let mut output = Vec::new();
        let summary = write(&input, &options, &mut output).unwrap();
        let output = String::from_utf8(output).unwrap();

        assert_eq!(summary.read, 2, "{output}");
        assert_eq!(summary.written, 1, "{output}");
        assert_eq!(output.matches("\tdel1\t").count(), 1, "{output}");
        assert!(!output.contains("\tsnv4\t"), "{output}");
    }

    #[test]
    fn parallel_filter_commits_owned_bgzf_output() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("calls.vcf");
        let output = directory.path().join("filtered.vcf.gz");
        fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tC\t10\tPASS\t.\n\
chr1\t20\t.\tG\tT\t5\tPASS\t.\n",
        )
        .unwrap();
        let options = StreamOptions {
            output_format: OutputFormat::VcfBgzf,
            expression: Some("QUAL >= 10".to_owned()),
            ..StreamOptions::default()
        };
        let transaction = AtomicFile::new(&output).unwrap();
        let summary = write_parallel(&input, &options, transaction.reopen().unwrap(), 2).unwrap();
        transaction.commit().unwrap();

        let mut decoded = String::new();
        bgzf::io::Reader::new(File::open(output).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(summary.read, 2);
        assert_eq!(summary.written, 1);
        assert!(decoded.contains("chr1\t10\t"), "{decoded}");
        assert!(!decoded.contains("chr1\t20\t"), "{decoded}");
    }
}
