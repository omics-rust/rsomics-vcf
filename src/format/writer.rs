use std::borrow::Cow;
use std::io::{BufWriter, Write};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{
    self as vcf,
    header::record::value::map::{format::Type as FormatType, info::Type},
    variant::{
        RecordBuf,
        io::Write as _,
        record_buf::{
            Samples,
            info::field::{Value, value::Array},
            samples::sample::Value as SampleValue,
        },
    },
};
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
        let result = match self {
            Self::Vcf(writer) => writer.write_variant_record(header, record),
            Self::VcfBgzf(writer) => writer.write_variant_record(header, record),
            Self::Bcf(writer) => {
                let record = bcf_record(header, record)?;
                writer.write_variant_record(header, record.as_ref())
            }
            Self::BcfRaw(writer) => {
                let record = bcf_record(header, record)?;
                writer.write_variant_record(header, record.as_ref())
            }
        };
        result.map_err(|error| map_write_error(error, &format!("writing variant record {number}")))
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
        let result = match self {
            Self::VcfBgzf(writer) => writer.write_variant_record(header, record),
            Self::Bcf(writer) => {
                let record = bcf_record(header, record)?;
                writer.write_variant_record(header, record.as_ref())
            }
        };
        result.map_err(|error| map_write_error(error, &format!("writing variant record {number}")))
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

fn bcf_record<'a>(header: &vcf::Header, record: &'a RecordBuf) -> Result<Cow<'a, RecordBuf>> {
    let missing_info = record
        .info()
        .as_ref()
        .iter()
        .filter_map(|(key, value)| value.is_none().then_some(key.clone()))
        .collect::<Vec<_>>();
    let genotype_index = record.samples().keys().as_ref().get_index_of("GT");
    let missing_genotypes = genotype_index.is_some_and(|index| {
        record
            .samples()
            .values()
            .any(|sample| sample.values().get(index).is_none_or(Option::is_none))
    });
    let missing_text = record
        .samples()
        .keys()
        .as_ref()
        .iter()
        .enumerate()
        .filter_map(|(index, key)| {
            let ty = header.formats().get(key)?.ty();
            if matches!(ty, FormatType::String | FormatType::Character)
                && key != "GT"
                && record
                    .samples()
                    .values()
                    .all(|sample| sample.values().get(index).is_none_or(Option::is_none))
            {
                Some((index, ty))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if missing_info.is_empty() && !missing_genotypes && missing_text.is_empty() {
        return Ok(Cow::Borrowed(record));
    }

    let mut record = record.clone();
    for key in missing_info {
        let definition = header.infos().get(&key).ok_or_else(|| {
            RsomicsError::InvalidInput(format!("INFO/{key} has no header definition"))
        })?;
        let value = match definition.ty() {
            Type::Integer => Value::Array(Array::Integer(vec![None])),
            Type::Float => Value::Array(Array::Float(vec![None])),
            Type::Character => Value::Array(Array::Character(vec![None])),
            Type::String => Value::Array(Array::String(vec![None])),
            Type::Flag => {
                return Err(RsomicsError::InvalidInput(format!(
                    "INFO/{key} flag has an explicit missing value"
                )));
            }
        };
        record.info_mut().insert(key, Some(value));
    }
    if genotype_index.is_some() || !missing_text.is_empty() {
        let (keys, mut values) = record.samples().clone().into();
        for row in &mut values {
            row.resize(keys.as_ref().len(), None);
            if let Some(index) = genotype_index
                && row[index].is_none()
            {
                row[index] = Some(SampleValue::Genotype(
                    ".".parse().expect("missing genotype is valid"),
                ));
            }
            for &(index, ty) in &missing_text {
                row[index] = Some(match ty {
                    FormatType::Character => SampleValue::Character('.'),
                    FormatType::String => SampleValue::String(".".to_owned()),
                    _ => unreachable!("missing text field has a text type"),
                });
            }
        }
        *record.samples_mut() = Samples::new(keys, values);
    }
    Ok(Cow::Owned(record))
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
    fn bcf_encodes_explicitly_missing_info_values() {
        let raw = "##fileformat=VCFv4.3\n\
##INFO=<ID=I,Number=1,Type=Integer,Description=\"I\">\n\
##INFO=<ID=F,Number=1,Type=Float,Description=\"F\">\n\
##INFO=<ID=C,Number=1,Type=Character,Description=\"C\">\n\
##INFO=<ID=S,Number=1,Type=String,Description=\"S\">\n\
##contig=<ID=chr1,length=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header: vcf::Header = raw.parse().unwrap();
        let record = vcf::Record::try_from(b"chr1\t1\t.\tA\tC\t.\tPASS\t.".as_slice()).unwrap();
        let mut record =
            vcf::variant::RecordBuf::try_from_variant_record(&header, &record).unwrap();
        for key in ["I", "F", "C", "S"] {
            record.info_mut().insert(key.to_owned(), None);
        }
        for parallel in [false, true] {
            let output = SharedOutput::default();
            let bytes = output.0.clone();
            let mut writer: Box<dyn VariantWriter> = if parallel {
                Box::new(ParallelWriter::new(output, OutputFormat::Bcf, 1).unwrap())
            } else {
                Box::new(Writer::new(output, OutputFormat::Bcf))
            };
            writer.write_header(&header, HeaderMode::Full).unwrap();
            writer.write_record(&header, &record, 1).unwrap();
            writer.finish().unwrap();

            let compressed = bytes.lock().unwrap().clone();
            let mut reader = bcf::io::Reader::new(Cursor::new(compressed));
            let header = reader.read_header().unwrap();
            let record = reader.record_bufs(&header).next().unwrap().unwrap();
            assert_eq!(record.info().as_ref().len(), 4, "parallel={parallel}");
        }
    }

    #[test]
    fn bcf_encodes_missing_genotypes_and_text_fields() {
        let raw = "##fileformat=VCFv4.3\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=C,Number=1,Type=Character,Description=\"C\">\n\
##FORMAT=<ID=S,Number=1,Type=String,Description=\"S\">\n\
##contig=<ID=chr1,length=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n";
        let header: vcf::Header = raw.parse().unwrap();
        let record = vcf::Record::try_from(
            b"chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:C:S\t0/1:.:.\t.:.:.".as_slice(),
        )
        .unwrap();
        let record = vcf::variant::RecordBuf::try_from_variant_record(&header, &record).unwrap();

        let output = SharedOutput::default();
        let bytes = output.0.clone();
        let mut writer = Writer::new(output, OutputFormat::Bcf);
        writer.write_header(&header, HeaderMode::Full).unwrap();
        writer.write_record(&header, &record, 1).unwrap();
        writer.finish().unwrap();

        let compressed = bytes.lock().unwrap().clone();
        let mut reader = bcf::io::Reader::new(Cursor::new(compressed));
        let header = reader.read_header().unwrap();
        let record = reader.record_bufs(&header).next().unwrap().unwrap();
        let sample = record.samples().get(&header, "B").unwrap();
        for key in ["GT", "C", "S"] {
            assert!(sample.get(key).is_some(), "{key}");
        }
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
