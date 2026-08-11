use std::io::Write;
use std::path::Path;

use rsomics_common::Result;
use serde::Serialize;

use crate::format::{HeaderMode, OutputFormat, Reader, RecordScratch, Writer};
use crate::regions::RegionSelection;

use super::{Options, Processor, gaps};

#[derive(Clone, Debug, Default)]
pub(crate) struct StreamOptions {
    pub output_format: OutputFormat,
    pub expression: Option<String>,
    pub filter: Options,
    pub gaps: gaps::Options,
    pub targets: Option<RegionSelection>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct Summary {
    pub read: u64,
    pub written: u64,
    pub output_format: OutputFormat,
}

pub(crate) fn write(input: &Path, options: &StreamOptions, output: impl Write) -> Result<Summary> {
    let mut reader = Reader::open(input)?;
    let (mut header, _, _) = reader.read_header()?;
    let mut processor = Processor::bind(
        &mut header,
        options.expression.as_deref(),
        options.filter.clone(),
        options.gaps,
    )?;
    let mut writer = Writer::new(output, options.output_format);
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

#[cfg(test)]
mod tests {
    use std::fs;

    use crate::format::OutputFormat;
    use crate::regions::{OverlapMode, RegionSelection};
    use crate::variant_type;

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
}
