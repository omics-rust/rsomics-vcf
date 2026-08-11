use std::io::Write;
use std::path::Path;

use rsomics_common::Result;

use crate::format::Writer;
use crate::regions::{IndexedRecords, RegionSet};

use super::{HeaderMode, Options, Summary, samples, selection};

pub(super) fn write_indexed(
    input: &Path,
    options: &Options,
    regions: &RegionSet,
    output: impl Write,
) -> Result<Summary> {
    let mut reader = IndexedRecords::open(input, regions)?;
    let projection = samples::Projection::new(reader.header(), options)?;
    let output_header = projection.header(reader.header(), options);
    let mut writer = Writer::new(output, options.output_format);
    writer.write_header(&output_header, options.header)?;

    if options.header == HeaderMode::HeaderOnly {
        writer.finish()?;
        return Ok(Summary {
            read: 0,
            written: 0,
            output_format: options.output_format,
        });
    }

    let mut written = 0;
    let read = reader.visit(|_, record, number| {
        if !selection::keep(&record, options)
            || options
                .targets
                .as_ref()
                .is_some_and(|targets| !targets.keeps(&record))
        {
            return Ok(());
        }
        let record = projection.apply(record, options)?;
        writer.write_record(&output_header, &record, number)?;
        written += 1;
        Ok(())
    })?;
    writer.finish()?;

    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}
