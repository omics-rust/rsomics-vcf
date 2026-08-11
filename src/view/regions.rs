use std::io::Write;
use std::path::Path;

use noodles_util::variant::io::indexed_reader;
use noodles_vcf::variant::RecordBuf;
use rsomics_common::{Context, Result};

use crate::format::Writer;
use crate::regions::{RegionSet, overlaps};

use super::{HeaderMode, Options, Summary, samples, selection};

pub(super) fn write_indexed(
    input: &Path,
    options: &Options,
    regions: &RegionSet,
    output: impl Write,
) -> Result<Summary> {
    let mut reader = indexed_reader::Builder::default()
        .build_from_path(input)
        .rs_with_context(|| format!("opening indexed variant input {}", input.display()))?;
    let header = reader
        .read_header()
        .rs_with_context(|| format!("reading variant header {}", input.display()))?;
    let query_regions = regions.merged(&header)?;
    let projection = samples::Projection::new(&header, options)?;
    let output_header = projection.header(&header, options);
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

    let mut read = 0;
    let mut written = 0;
    for (region_index, region) in query_regions.iter().enumerate() {
        let records = reader
            .query(&header, region)
            .rs_with_context(|| format!("querying region {region}"))?;
        for result in records {
            let number = read + 1;
            let record =
                result.rs_with_context(|| format!("reading indexed variant record {number}"))?;
            let record = RecordBuf::try_from_variant_record(&header, record.as_ref())
                .map_err(|error| super::invalid(input, number, "decoding indexed input", error))?;
            read += 1;

            if !overlaps(&record, region.interval(), regions.overlap())
                || query_regions[..region_index].iter().any(|previous| {
                    previous.name() == region.name()
                        && overlaps(&record, previous.interval(), regions.overlap())
                })
                || !selection::keep(&record, options)
                || options
                    .targets
                    .as_ref()
                    .is_some_and(|targets| !targets.matches(&record))
            {
                continue;
            }

            let record = projection.apply(record, options)?;
            writer.write_record(&output_header, &record, number)?;
            written += 1;
        }
    }
    writer.finish()?;

    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}
