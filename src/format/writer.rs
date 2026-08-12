use std::io::{BufWriter, Write};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rayon::{ThreadPool, ThreadPoolBuilder};
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

pub(crate) enum ParallelWriter<W>
where
    W: Write + Send + 'static,
{
    VcfBgzf(vcf::io::Writer<BoundedBgzf<W>>),
    Bcf(bcf::io::Writer<BoundedBgzf<W>>),
}

pub(crate) trait VariantWriter {
    fn write_header(&mut self, header: &vcf::Header, mode: HeaderMode) -> Result<()>;

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
        number: u64,
    ) -> Result<()>;

    fn supports_vcf_records(&self) -> bool;

    fn write_vcf_record(&mut self, record: &[u8], number: u64) -> Result<()>;

    fn finish(&mut self) -> Result<()>;
}

pub(crate) struct BoundedBgzf<W>
where
    W: Write + Send + 'static,
{
    writer: bgzf::io::MultithreadedWriter<W>,
    pool: ThreadPool,
    finished: bool,
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

impl<W> VariantWriter for Writer<W>
where
    W: Write,
{
    fn write_header(&mut self, header: &vcf::Header, mode: HeaderMode) -> Result<()> {
        Self::write_header(self, header, mode)
    }

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
        number: u64,
    ) -> Result<()> {
        Self::write_record(self, header, record, number)
    }

    fn supports_vcf_records(&self) -> bool {
        matches!(self, Self::Vcf(_) | Self::VcfBgzf(_))
    }

    fn write_vcf_record(&mut self, record: &[u8], number: u64) -> Result<()> {
        Self::write_vcf_record(self, record, number)
    }

    fn finish(&mut self) -> Result<()> {
        Self::finish(self)
    }
}

impl<W> ParallelWriter<W>
where
    W: Write + Send + 'static,
{
    pub(crate) fn new(output: W, format: OutputFormat, workers: usize) -> Result<Self> {
        if !(1..=256).contains(&workers) {
            return Err(RsomicsError::ConfigError(format!(
                "BGZF compression workers must be between 1 and 256, got {workers}"
            )));
        }
        match format {
            OutputFormat::VcfBgzf => Ok(Self::VcfBgzf(vcf::io::Writer::new(BoundedBgzf::new(
                output, workers,
            )?))),
            OutputFormat::Bcf => Ok(Self::Bcf(bcf::io::Writer::from(BoundedBgzf::new(
                output, workers,
            )?))),
            OutputFormat::Vcf | OutputFormat::BcfRaw => Err(RsomicsError::ConfigError(
                "compression workers require BGZF VCF or BCF output".to_owned(),
            )),
        }
    }

    #[cfg(test)]
    fn worker_count(&self) -> usize {
        match self {
            Self::VcfBgzf(writer) => writer.get_ref().worker_count(),
            Self::Bcf(writer) => writer.get_ref().worker_count(),
        }
    }
}

impl<W> VariantWriter for ParallelWriter<W>
where
    W: Write + Send + 'static,
{
    fn write_header(&mut self, header: &vcf::Header, mode: HeaderMode) -> Result<()> {
        if mode == HeaderMode::None {
            return Ok(());
        }
        match self {
            Self::VcfBgzf(writer) => writer.write_header(header),
            Self::Bcf(writer) => writer.write_header(header),
        }
        .map_err(|error| map_write_error(error, "writing variant header"))
    }

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
        number: u64,
    ) -> Result<()> {
        match self {
            Self::VcfBgzf(writer) => writer.write_variant_record(header, record),
            Self::Bcf(writer) => writer.write_variant_record(header, record),
        }
        .map_err(|error| map_write_error(error, &format!("writing variant record {number}")))
    }

    fn supports_vcf_records(&self) -> bool {
        matches!(self, Self::VcfBgzf(_))
    }

    fn write_vcf_record(&mut self, record: &[u8], number: u64) -> Result<()> {
        let result = match self {
            Self::VcfBgzf(writer) => writer.get_mut().write_all(record),
            Self::Bcf(_) => unreachable!(),
        }
        .and_then(|()| match self {
            Self::VcfBgzf(writer) => writer.get_mut().write_all(b"\n"),
            Self::Bcf(_) => unreachable!(),
        });
        result.map_err(|error| {
            RsomicsError::Io(std::io::Error::new(
                error.kind(),
                format!("writing variant record {number}: {error}"),
            ))
        })
    }

    fn finish(&mut self) -> Result<()> {
        match self {
            Self::VcfBgzf(writer) => writer.get_mut().finish(),
            Self::Bcf(writer) => writer.get_mut().finish(),
        }
        .map_err(RsomicsError::Io)
    }
}

impl<W> BoundedBgzf<W>
where
    W: Write + Send + 'static,
{
    fn new(output: W, workers: usize) -> Result<Self> {
        let pool = ThreadPoolBuilder::new()
            .num_threads(workers)
            .thread_name(|index| format!("rsomics-bgzf-{index}"))
            .build()
            .map_err(|error| {
                RsomicsError::ConfigError(format!("creating BGZF compression workers: {error}"))
            })?;
        let writer = pool.install(|| bgzf::io::MultithreadedWriter::new(output));
        Ok(Self {
            writer,
            pool,
            finished: false,
        })
    }

    #[cfg(test)]
    fn worker_count(&self) -> usize {
        self.pool.current_num_threads()
    }

    fn finish(&mut self) -> std::io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        self.pool.install(|| self.writer.flush())?;
        self.writer.finish().map(drop)
    }
}

impl<W> Write for BoundedBgzf<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.pool.install(|| self.writer.write(buffer))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.pool.install(|| self.writer.flush())
    }
}

impl<W> Drop for BoundedBgzf<W>
where
    W: Write + Send + 'static,
{
    fn drop(&mut self) {
        let _ = self.finish();
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

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read};
    use std::sync::{Arc, Mutex};

    use super::*;

    #[derive(Clone, Default)]
    struct SharedOutput(Arc<Mutex<Vec<u8>>>);

    impl Write for SharedOutput {
        fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(buffer);
            Ok(buffer.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn parallel_bgzf_uses_the_requested_bounded_pool() {
        let output = SharedOutput::default();
        let bytes = output.0.clone();
        let mut writer = ParallelWriter::new(output, OutputFormat::VcfBgzf, 2).unwrap();
        writer
            .write_header(&vcf::Header::default(), HeaderMode::Full)
            .unwrap();
        assert_eq!(writer.worker_count(), 2);
        writer.finish().unwrap();

        let compressed = bytes.lock().unwrap().clone();
        let mut reader = bgzf::io::Reader::new(Cursor::new(compressed));
        let mut decoded = String::new();
        reader.read_to_string(&mut decoded).unwrap();
        assert!(decoded.starts_with("##fileformat=VCF"), "{decoded}");
    }

    #[test]
    fn parallel_bcf_finishes_and_decodes() {
        let output = SharedOutput::default();
        let bytes = output.0.clone();
        let mut writer = ParallelWriter::new(output, OutputFormat::Bcf, 1).unwrap();
        writer
            .write_header(&vcf::Header::default(), HeaderMode::Full)
            .unwrap();
        writer.finish().unwrap();

        let compressed = bytes.lock().unwrap().clone();
        let mut reader = bcf::io::Reader::new(Cursor::new(compressed));
        reader.read_header().unwrap();
    }

    #[test]
    fn parallel_bgzf_rejects_unbounded_or_uncompressed_requests() {
        for workers in [0, 257] {
            assert!(
                ParallelWriter::new(SharedOutput::default(), OutputFormat::VcfBgzf, workers)
                    .is_err()
            );
        }
        assert!(ParallelWriter::new(SharedOutput::default(), OutputFormat::Vcf, 2).is_err());
    }
}
