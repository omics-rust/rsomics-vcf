use std::io::{BufWriter, Write};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_common::{Result, RsomicsError};

use super::{HeaderMode, OutputFormat};

pub(crate) enum Writer<W>
where
    W: Write,
{
    Vcf(vcf::io::Writer<BufWriter<W>>),
    VcfBgzf(vcf::io::Writer<bgzf::io::Writer<W>>),
    Bcf(bcf::io::Writer<bgzf::io::Writer<W>>),
    BcfRaw(bcf::io::Writer<BufWriter<W>>),
}

impl<W> Writer<W>
where
    W: Write,
{
    pub(crate) fn new(output: W, format: OutputFormat) -> Self {
        match format {
            OutputFormat::Vcf => Self::Vcf(vcf::io::Writer::new(BufWriter::new(output))),
            OutputFormat::VcfBgzf => {
                Self::VcfBgzf(vcf::io::Writer::new(bgzf::io::Writer::new(output)))
            }
            OutputFormat::Bcf => Self::Bcf(bcf::io::Writer::new(output)),
            OutputFormat::BcfRaw => Self::BcfRaw(bcf::io::Writer::from(BufWriter::new(output))),
        }
    }

    pub(crate) fn write_header(&mut self, header: &vcf::Header, mode: HeaderMode) -> Result<()> {
        if mode == HeaderMode::None {
            return Ok(());
        }
        match self {
            Self::Vcf(writer) => writer.write_header(header),
            Self::VcfBgzf(writer) => writer.write_header(header),
            Self::Bcf(writer) => writer.write_header(header),
            Self::BcfRaw(writer) => writer.write_header(header),
        }
        .map_err(|error| map_write_error(error, "writing variant header"))
    }

    pub(crate) fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
        number: u64,
    ) -> Result<()> {
        match self {
            Self::Vcf(writer) => writer.write_variant_record(header, record),
            Self::VcfBgzf(writer) => writer.write_variant_record(header, record),
            Self::Bcf(writer) => writer.write_variant_record(header, record),
            Self::BcfRaw(writer) => writer.write_variant_record(header, record),
        }
        .map_err(|error| map_write_error(error, &format!("writing variant record {number}")))
    }

    pub(crate) fn write_vcf_record(&mut self, record: &[u8], number: u64) -> Result<()> {
        let result = match self {
            Self::Vcf(writer) => writer.get_mut().write_all(record),
            Self::VcfBgzf(writer) => writer.get_mut().write_all(record),
            Self::Bcf(_) | Self::BcfRaw(_) => unreachable!(),
        }
        .and_then(|()| match self {
            Self::Vcf(writer) => writer.get_mut().write_all(b"\n"),
            Self::VcfBgzf(writer) => writer.get_mut().write_all(b"\n"),
            Self::Bcf(_) | Self::BcfRaw(_) => unreachable!(),
        });
        result.map_err(|error| {
            RsomicsError::Io(std::io::Error::new(
                error.kind(),
                format!("writing variant record {number}: {error}"),
            ))
        })
    }

    pub(crate) fn finish(&mut self) -> Result<()> {
        match self {
            Self::Vcf(writer) => writer.get_mut().flush(),
            Self::VcfBgzf(writer) => writer.get_mut().try_finish(),
            Self::Bcf(writer) => writer.try_finish(),
            Self::BcfRaw(writer) => writer.get_mut().flush(),
        }
        .map_err(RsomicsError::Io)
    }
}

fn map_write_error(error: std::io::Error, context: &str) -> RsomicsError {
    let message = format!("{context}: {error}");
    if error.kind() == std::io::ErrorKind::InvalidInput {
        RsomicsError::InvalidInput(message)
    } else {
        RsomicsError::Io(std::io::Error::new(error.kind(), message))
    }
}
