mod columns;
mod edit;
mod header;
mod matching;
mod set_id;
mod source;

use std::io::Write;
use std::path::{Path, PathBuf};

use noodles_vcf::{
    Header,
    header::record::value::{
        Map,
        map::{Info, info},
    },
    variant::{RecordBuf, record_buf::info::field::Value},
};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

pub(crate) use columns::ColumnSpec;
pub(crate) use header::HeaderOptions;
pub(crate) use matching::{OverlapFractions, PairLogic};

use crate::{
    expression::Compiled,
    filter::Logic,
    format::{
        HeaderMode, HeaderTypes, OutputFormat, ParallelWriter, Reader, RecordScratch,
        VariantWriter, Writer,
    },
    regions::{IndexedRecords, RegionSet},
};
use edit::{Editor, SampleSelection};
use header::HeaderPlan;
use set_id::IdPlan;
use source::AnnotationSource;

#[derive(Clone, Debug)]
pub(crate) struct Options {
    pub(crate) source: Option<SourceOptions>,
    pub(crate) header: HeaderOptions,
    pub(crate) set_id: Option<String>,
    pub(crate) expression: Option<String>,
    pub(crate) expression_logic: Logic,
    pub(crate) keep_sites: bool,
    pub(crate) mark_sites: Option<MarkSites>,
    pub(crate) regions: Option<RegionSet>,
    pub(crate) output_format: OutputFormat,
}

#[derive(Clone, Debug)]
pub(crate) struct SourceOptions {
    pub(crate) path: PathBuf,
    pub(crate) columns: ColumnSpec,
    pub(crate) samples: Option<SampleRequest>,
    pub(crate) pair_logic: PairLogic,
    pub(crate) min_overlap: OverlapFractions,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum SampleRequest {
    List(String),
    File(PathBuf),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum MarkSites {
    Present(String),
    Absent(String),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub(crate) struct Summary {
    pub(crate) read: u64,
    pub(crate) written: u64,
    pub(crate) annotated: u64,
    pub(crate) unchanged: u64,
    pub(crate) filtered: u64,
    pub(crate) output_format: OutputFormat,
}

struct SourcePlan {
    source: AnnotationSource,
    editor: Editor,
    pair_logic: PairLogic,
    min_overlap: OverlapFractions,
}

struct Processor {
    input_header: Header,
    output_header: Header,
    header: HeaderPlan,
    id: Option<IdPlan>,
    expression: Option<Compiled>,
    expression_logic: Logic,
    keep_sites: bool,
    mark_sites: Option<MarkSites>,
    source: Option<SourcePlan>,
    summary: Summary,
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
    write_with_writer(
        input,
        options,
        ParallelWriter::new(output, options.output_format, workers)?,
    )
}

fn write_with_writer(
    input: &Path,
    options: &Options,
    mut writer: impl VariantWriter,
) -> Result<Summary> {
    validate_options(options)?;
    if let Some(regions) = &options.regions {
        if input == Path::new("-") {
            return Err(RsomicsError::ConfigError(
                "indexed regions require a named input".to_owned(),
            ));
        }
        let mut reader = IndexedRecords::open(input, regions)?;
        let mut processor = Processor::bind(reader.header(), options)?;
        writer.write_header(&processor.output_header, HeaderMode::Full)?;
        let read =
            reader.visit(|_, record, number| processor.process(record, number, &mut writer))?;
        processor.summary.read = read;
        writer.finish()?;
        return Ok(processor.summary);
    }

    let mut reader = Reader::open(input)?;
    let (input_header, _, _) = reader.read_header()?;
    let mut processor = Processor::bind(&input_header, options)?;
    writer.write_header(&processor.output_header, HeaderMode::Full)?;
    let mut scratch = RecordScratch::default();
    loop {
        let number = processor
            .summary
            .read
            .checked_add(1)
            .ok_or_else(|| invalid("variant record count exceeds u64"))?;
        let Some(record) = reader.read_record(&input_header, &mut scratch, number)? else {
            break;
        };
        processor.summary.read = number;
        processor.process(record, number, &mut writer)?;
    }
    writer.finish()?;
    Ok(processor.summary)
}

impl Processor {
    fn bind(input_header: &Header, options: &Options) -> Result<Self> {
        let header = HeaderPlan::bind(input_header, options.header.clone())?;
        let mut output_header = input_header.clone();
        header.prepare_with_retained_removals(
            &mut output_header,
            options.expression.is_some() && options.keep_sites,
        )?;

        let source = options
            .source
            .as_ref()
            .map(|options| -> Result<SourcePlan> {
                let source =
                    AnnotationSource::open(&options.path, input_header, options.columns.clone())?;
                let mut editor = Editor::bind(
                    &mut output_header,
                    source.header(),
                    source.columns().clone(),
                )?;
                if let Some(request) = &options.samples {
                    let source_header = source.header().ok_or_else(|| {
                        invalid("annotation sample selection requires a VCF or BCF source")
                    })?;
                    let selection = match request {
                        SampleRequest::List(list) => {
                            SampleSelection::bind(source_header, &output_header, Some(list), None)?
                        }
                        SampleRequest::File(path) => {
                            SampleSelection::bind(source_header, &output_header, None, Some(path))?
                        }
                    };
                    editor.set_samples(selection);
                }
                Ok(SourcePlan {
                    source,
                    editor,
                    pair_logic: options.pair_logic,
                    min_overlap: options.min_overlap,
                })
            })
            .transpose()?;

        if let Some(mark) = &options.mark_sites {
            prepare_mark(&mut output_header, mark.tag())?;
        }
        let expression = options
            .expression
            .as_deref()
            .map(|source| {
                Compiled::bind(source, input_header).map_err(|error| {
                    RsomicsError::ConfigError(format!("invalid annotate expression: {error}"))
                })
            })
            .transpose()?;
        let schema = header_types(&output_header)?;
        let id = IdPlan::bind(options.set_id.as_deref(), &schema)?;

        Ok(Self {
            input_header: input_header.clone(),
            output_header,
            header,
            id,
            expression,
            expression_logic: options.expression_logic,
            keep_sites: options.keep_sites,
            mark_sites: options.mark_sites.clone(),
            source,
            summary: Summary {
                read: 0,
                written: 0,
                annotated: 0,
                unchanged: 0,
                filtered: 0,
                output_format: options.output_format,
            },
        })
    }

    fn process(
        &mut self,
        mut record: RecordBuf,
        number: u64,
        writer: &mut impl VariantWriter,
    ) -> Result<()> {
        if let Some(expression) = &self.expression {
            let passes = expression
                .evaluate(&self.input_header, &record)
                .map_err(|error| {
                    invalid(format!(
                        "evaluating annotate expression at record {number}: {error}"
                    ))
                })?
                .site_passes();
            if !self.expression_logic.accepts(passes) {
                if !self.keep_sites {
                    self.summary.filtered += 1;
                    return Ok(());
                }
                if self.header.apply_renames(&mut record)? {
                    self.summary.annotated += 1;
                } else {
                    self.summary.unchanged += 1;
                }
                return self.write_record(record, writer);
            }
        }

        let header = &self.header;
        let id = &mut self.id;
        let output_header = &self.output_header;
        let source = &mut self.source;
        let mut changed;
        let matched = if let Some(source) = source {
            let current = if source.editor.requires_all_matches() {
                source
                    .source
                    .matches(&record, source.pair_logic, source.min_overlap)?
            } else {
                source
                    .source
                    .first_match(&record, source.pair_logic, source.min_overlap)?
                    .into_iter()
                    .collect()
            };
            changed = header.apply(&mut record)?;
            if let Some(id) = id {
                changed |= id.apply(output_header, &mut record)?;
            }
            if !current.is_empty() {
                changed |=
                    source
                        .editor
                        .apply_info_matches(output_header, &current, &mut record)?;
                changed |= source
                    .editor
                    .apply_samples(output_header, &current[0], &mut record)?;
                true
            } else {
                false
            }
        } else {
            changed = header.apply(&mut record)?;
            if let Some(id) = id {
                changed |= id.apply(output_header, &mut record)?;
            }
            false
        };
        if let Some(mark) = &self.mark_sites
            && mark.applies(matched)
        {
            changed |= set_mark(&mut record, mark.tag());
        }

        if changed {
            self.summary.annotated += 1;
        } else {
            self.summary.unchanged += 1;
        }
        self.write_record(record, writer)
    }

    fn write_record(&mut self, record: RecordBuf, writer: &mut impl VariantWriter) -> Result<()> {
        self.summary.written = self
            .summary
            .written
            .checked_add(1)
            .ok_or_else(|| invalid("written record count exceeds u64"))?;
        writer.write_record(&self.output_header, &record, self.summary.written)
    }
}

impl MarkSites {
    fn tag(&self) -> &str {
        match self {
            Self::Present(tag) | Self::Absent(tag) => tag,
        }
    }

    fn applies(&self, matched: bool) -> bool {
        match self {
            Self::Present(_) => matched,
            Self::Absent(_) => !matched,
        }
    }
}

fn validate_options(options: &Options) -> Result<()> {
    let header = &options.header;
    let has_header_action = !header.appended.is_empty()
        || header.remove.is_some()
        || header.rename_chromosomes.is_some()
        || header.rename_annotations.is_some();
    if options.source.is_none()
        && !has_header_action
        && options.set_id.is_none()
        && options.mark_sites.is_none()
    {
        return Err(RsomicsError::ConfigError(
            "annotate requires an annotation source, header edit, removal, rename, or --set-id"
                .to_owned(),
        ));
    }
    if options.mark_sites.is_some() && options.source.is_none() {
        return Err(RsomicsError::ConfigError(
            "--mark-sites requires --annotations".to_owned(),
        ));
    }
    if let Some(source) = &options.source {
        source.min_overlap.validate()?;
        if source.columns.transfers().is_empty() && options.mark_sites.is_none() {
            return Err(RsomicsError::ConfigError(
                "annotation columns must transfer a field unless --mark-sites is used".to_owned(),
            ));
        }
    }
    if options.keep_sites && options.expression.is_none() {
        return Err(RsomicsError::ConfigError(
            "--keep-sites requires --include or --exclude".to_owned(),
        ));
    }
    Ok(())
}

fn prepare_mark(header: &mut Header, tag: &str) -> Result<()> {
    if !valid_tag(tag) {
        return Err(RsomicsError::ConfigError(format!(
            "invalid mark-sites tag: {tag:?}"
        )));
    }
    if let Some(existing) = header.infos().get(tag) {
        if existing.number() != info::Number::Count(0) || existing.ty() != info::Type::Flag {
            return Err(RsomicsError::ConfigError(format!(
                "INFO/{tag} must have Number=0,Type=Flag for --mark-sites"
            )));
        }
    } else {
        header.infos_mut().insert(
            tag.to_owned(),
            Map::<Info>::new(
                info::Number::Count(0),
                info::Type::Flag,
                format!("Sites listed in {tag}"),
            ),
        );
    }
    Ok(())
}

fn set_mark(record: &mut RecordBuf, tag: &str) -> bool {
    !matches!(
        record.info_mut().insert(tag.to_owned(), Some(Value::Flag)),
        Some(Some(Value::Flag))
    )
}

fn header_types(header: &Header) -> Result<HeaderTypes> {
    let mut raw = Vec::new();
    noodles_vcf::io::Writer::new(&mut raw)
        .write_header(header)
        .map_err(RsomicsError::Io)?;
    HeaderTypes::parse(&raw).map_err(|error| {
        RsomicsError::InvalidInput(format!("building annotation output schema: {error}"))
    })
}

fn valid_tag(tag: &str) -> bool {
    !tag.is_empty()
        && tag
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    struct FailingOutput;

    #[derive(Default)]
    struct FailingFlush(Vec<u8>);

    impl Write for FailingOutput {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed output"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed output"))
        }
    }

    impl Write for FailingFlush {
        fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
            self.0.extend_from_slice(buffer);
            Ok(buffer.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "closed output"))
        }
    }

    fn options() -> Options {
        Options {
            source: None,
            header: HeaderOptions {
                remove: Some("ID".to_owned()),
                ..HeaderOptions::default()
            },
            set_id: None,
            expression: None,
            expression_logic: Logic::Include,
            keep_sites: false,
            mark_sites: None,
            regions: None,
            output_format: OutputFormat::Vcf,
        }
    }

    fn input(directory: &Path) -> PathBuf {
        let input = directory.join("input.vcf");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\told\tA\tC\t.\tPASS\t.\n",
        )
        .unwrap();
        input
    }

    #[test]
    fn propagates_output_write_failures() {
        let directory = tempfile::tempdir().unwrap();
        assert!(write(&input(directory.path()), &options(), FailingOutput).is_err());
    }

    #[test]
    fn propagates_output_finish_failures() {
        let directory = tempfile::tempdir().unwrap();
        assert!(
            write(
                &input(directory.path()),
                &options(),
                FailingFlush::default()
            )
            .is_err()
        );
    }
}
