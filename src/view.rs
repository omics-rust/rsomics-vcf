mod regions;
mod samples;
mod selection;

use std::io::Write;
use std::path::Path;

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::RecordBuf};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::format::{Reader, Writer, reformat_record};

pub use crate::format::{HeaderMode, OutputFormat};
pub use crate::regions::{OverlapMode, RegionSelection, RegionSet};
pub use samples::SampleSelection;
pub use selection::TypeSelection;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IdSelection {
    Known,
    Novel,
}

#[derive(Clone, Debug)]
pub struct Options {
    pub output_format: OutputFormat,
    pub header: HeaderMode,
    pub samples: Option<SampleSelection>,
    pub drop_genotypes: bool,
    pub update_info: bool,
    pub filters: Vec<String>,
    pub ids: Option<IdSelection>,
    pub types: Option<TypeSelection>,
    pub min_alleles: Option<usize>,
    pub max_alleles: Option<usize>,
    pub targets: Option<RegionSelection>,
    pub regions: Option<RegionSet>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            output_format: OutputFormat::default(),
            header: HeaderMode::default(),
            samples: None,
            drop_genotypes: false,
            update_info: true,
            filters: Vec::new(),
            ids: None,
            types: None,
            min_alleles: None,
            max_alleles: None,
            targets: None,
            regions: None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub read: u64,
    pub written: u64,
    pub output_format: OutputFormat,
}

pub fn write(input: &Path, options: &Options, output: impl Write) -> Result<Summary> {
    validate_options(input, options)?;
    if let Some(regions) = &options.regions {
        return write_regions(input, options, regions, output);
    }

    let mut reader = Reader::open(input)?;
    let (header, _, schema) = reader.read_header()?;
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
    if reader.is_text()
        && matches!(
            options.output_format,
            OutputFormat::Vcf | OutputFormat::VcfBgzf
        )
        && can_use_text_path(options)
    {
        let mut input_record = Vec::new();
        let mut canonical = Vec::new();
        let mut records = 0;
        loop {
            let number = records + 1;
            if reader.read_text_record(&mut input_record, usize_number(number)?)? == 0 {
                break;
            }
            reformat_record(&input_record, &schema, &mut canonical)
                .map_err(|error| invalid(input, number, "parsing VCF", error))?;
            writer.write_vcf_record(&canonical, number)?;
            records += 1;
        }
        writer.finish()?;
        return Ok(Summary {
            read: records,
            written: records,
            output_format: options.output_format,
        });
    }

    let mut bcf_record = bcf::Record::default();
    let mut text_record = Vec::new();
    let mut read = 0;
    let mut written = 0;
    loop {
        let number = read + 1;
        let record = if reader.is_text() {
            if reader.read_text_record(&mut text_record, usize_number(number)?)? == 0 {
                break;
            }
            let record = vcf::Record::try_from(text_record.as_slice())
                .map_err(|error| invalid(input, number, "parsing VCF", error))?;
            RecordBuf::try_from_variant_record(&header, &record)
                .map_err(|error| invalid(input, number, "decoding VCF", error))?
        } else {
            if reader.read_bcf_record(&mut bcf_record, usize_number(number)?)? == 0 {
                break;
            }
            RecordBuf::try_from_variant_record(&header, &bcf_record)
                .map_err(|error| invalid(input, number, "decoding BCF", error))?
        };
        read += 1;

        if !selection::keep(&record, options)
            || options
                .targets
                .as_ref()
                .is_some_and(|targets| !targets.keeps(&record))
        {
            continue;
        }

        let record = projection.apply(record, options)?;
        writer.write_record(&output_header, &record, number)?;
        written += 1;
    }

    writer.finish()?;
    Ok(Summary {
        read,
        written,
        output_format: options.output_format,
    })
}

fn can_use_text_path(options: &Options) -> bool {
    options.samples.is_none()
        && !options.drop_genotypes
        && options.filters.is_empty()
        && options.ids.is_none()
        && options.types.is_none()
        && options.min_alleles.is_none()
        && options.max_alleles.is_none()
        && options.targets.is_none()
}

fn write_regions(
    input: &Path,
    options: &Options,
    regions: &RegionSet,
    output: impl Write,
) -> Result<Summary> {
    regions::write_indexed(input, options, regions, output)
}

fn validate_options(input: &Path, options: &Options) -> Result<()> {
    if options.header == HeaderMode::None
        && matches!(
            options.output_format,
            OutputFormat::Bcf | OutputFormat::BcfRaw
        )
    {
        return Err(RsomicsError::ConfigError(
            "BCF output requires a header".to_owned(),
        ));
    }
    if options.regions.is_some() && input == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "indexed regions require a named input".to_owned(),
        ));
    }
    if options
        .min_alleles
        .zip(options.max_alleles)
        .is_some_and(|(minimum, maximum)| minimum > maximum)
    {
        return Err(RsomicsError::ConfigError(
            "--min-alleles cannot exceed --max-alleles".to_owned(),
        ));
    }
    if options.min_alleles == Some(0) || options.max_alleles == Some(0) {
        return Err(RsomicsError::ConfigError(
            "allele limits must be at least 1".to_owned(),
        ));
    }
    Ok(())
}

fn invalid(input: &Path, number: u64, action: &str, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "{}: {action} variant record {number}: {error}",
        input.display()
    ))
}

fn usize_number(number: u64) -> Result<usize> {
    usize::try_from(number)
        .map_err(|_| RsomicsError::InvalidInput("variant record count exceeds usize".to_owned()))
}
