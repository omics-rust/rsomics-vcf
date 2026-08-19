use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::format::{
    HeaderMode, OutputFormat, ParallelWriter, Reader, RecordScratch, VariantWriter, Writer,
};

use super::{Options, Program};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct Summary {
    pub(crate) read: u64,
    pub(crate) changed_records: u64,
    pub(crate) changed_genotypes: u64,
    pub(crate) changed_alleles: u64,
    pub(crate) output_format: OutputFormat,
}

pub(crate) fn write(input: &Path, options: &Options, output: impl Write) -> Result<Summary> {
    write_with_writer(input, options, Writer::new(output, options.output_format))
}

pub(crate) fn write_parallel<W>(
    input: &Path,
    options: &Options,
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
    options: &Options,
    mut writer: impl VariantWriter,
) -> Result<Summary> {
    let mut reader = Reader::open(input)?;
    let (header, _, _) = reader.read_header()?;
    let mut program = Program::bind(
        &header,
        options.target.clone(),
        options.replacement.clone(),
        options.query.clone(),
        options.seed,
    )?;
    writer.write_header(&header, HeaderMode::Full)?;

    let mut scratch = RecordScratch::default();
    let mut summary = Summary {
        read: 0,
        changed_records: 0,
        changed_genotypes: 0,
        changed_alleles: 0,
        output_format: options.output_format,
    };
    loop {
        let number = checked_add(summary.read, 1, "record count")?;
        let Some(mut record) = reader.read_record(&header, &mut scratch, number)? else {
            break;
        };
        summary.read = number;
        let change = program.apply(&header, &mut record, number)?;
        if change.genotypes > 0 {
            summary.changed_records =
                checked_add(summary.changed_records, 1, "changed record count")?;
        }
        summary.changed_genotypes = checked_add(
            summary.changed_genotypes,
            change.genotypes,
            "changed genotype count",
        )?;
        summary.changed_alleles = checked_add(
            summary.changed_alleles,
            change.alleles,
            "changed allele count",
        )?;
        writer.write_record(&header, &record, number)?;
    }
    writer.finish()?;
    Ok(summary)
}

fn checked_add(left: u64, right: u64, name: &str) -> Result<u64> {
    left.checked_add(right)
        .ok_or_else(|| RsomicsError::InvalidInput(format!("{name} exceeds u64")))
}
