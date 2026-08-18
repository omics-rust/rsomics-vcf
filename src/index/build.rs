use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

use bstr::BString;
use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{
    self as noodles_csi,
    binning_index::{
        self,
        index::{
            header::{Builder as HeaderBuilder, ReferenceSequenceNames},
            reference_sequence::{bin::Chunk, index::BinnedIndex},
        },
    },
};
use noodles_tabix as tabix;
use noodles_vcf;
use rayon::ThreadPoolBuilder;
use rsomics_common::{Context, Result, RsomicsError, write_atomic};

use crate::format::bgzf::EOF_BLOCK;

use super::bcf_record;
use super::csi::{LinearIndex, apply_loffsets, settings};
use super::vcf::{parse_interval, trim_line_ending};
use super::{BuildOptions, BuildSummary, IndexKind, VariantFormat};

type Source = Box<dyn Read + Send>;

enum Index {
    Csi(noodles_csi::Index),
    Tbi(tabix::Index),
}

struct Built {
    index: Index,
    format: VariantFormat,
    min_shift: u8,
    depth: u8,
    records: u64,
    reference_sequences: usize,
}

pub(super) fn create(input: &Path, output: &Path, options: BuildOptions) -> Result<BuildSummary> {
    if input != Path::new("-") && input == output {
        return Err(RsomicsError::ConfigError(
            "input and index output must be different files".to_owned(),
        ));
    }
    if !options.force && output.exists() {
        return Err(RsomicsError::InvalidInput(format!(
            "index already exists: {} (use --force to replace it)",
            output.display()
        )));
    }
    if !(1..=30).contains(&options.min_shift) {
        return Err(RsomicsError::ConfigError(
            "CSI min-shift must be in 1..=30".to_owned(),
        ));
    }

    let source = open_source(input)?;
    let built = if options.threads == 0 {
        build(
            bgzf::io::Reader::new(source),
            options.kind,
            options.min_shift,
        )
    } else {
        let pool = ThreadPoolBuilder::new()
            .num_threads(options.threads)
            .build()
            .map_err(|error| RsomicsError::ConfigError(error.to_string()))?;
        pool.install(|| {
            build(
                bgzf::io::MultithreadedReader::new(source),
                options.kind,
                options.min_shift,
            )
        })
    }
    .map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "{}: indexing variant stream: {error}",
            source_name(input)
        ))
    })?;

    write_atomic(output, |file| write_index(file, &built.index))?;

    Ok(BuildSummary {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        format: built.format,
        kind: options.kind,
        min_shift: built.min_shift,
        depth: built.depth,
        records: built.records,
        reference_sequences: built.reference_sequences,
    })
}

fn open_source(input: &Path) -> Result<Source> {
    if input == Path::new("-") {
        return Ok(Box::new(io::stdin()));
    }

    let mut file = File::open(input)
        .rs_with_context(|| format!("opening variant input {}", input.display()))?;
    let length = file
        .metadata()
        .rs_with_context(|| format!("reading variant metadata {}", input.display()))?
        .len();
    if length < EOF_BLOCK.len() as u64 {
        return Err(RsomicsError::InvalidInput(format!(
            "{}: input is not a complete BGZF stream",
            input.display()
        )));
    }
    file.seek(SeekFrom::End(-(EOF_BLOCK.len() as i64)))
        .rs_with_context(|| format!("checking BGZF terminator {}", input.display()))?;
    let mut eof = [0; EOF_BLOCK.len()];
    file.read_exact(&mut eof)
        .rs_with_context(|| format!("checking BGZF terminator {}", input.display()))?;
    if eof != EOF_BLOCK {
        return Err(RsomicsError::InvalidInput(format!(
            "{}: BGZF end-of-file marker is missing",
            input.display()
        )));
    }
    file.rewind()
        .rs_with_context(|| format!("rewinding variant input {}", input.display()))?;
    Ok(Box::new(file))
}

fn build<R>(mut reader: R, kind: IndexKind, min_shift: u8) -> io::Result<Built>
where
    R: bgzf::io::BufRead,
{
    let format = if reader.fill_buf()?.starts_with(b"BCF") {
        VariantFormat::Bcf
    } else {
        VariantFormat::Vcf
    };

    match (format, kind) {
        (VariantFormat::Vcf, _) => build_vcf(reader, kind, min_shift),
        (VariantFormat::Bcf, IndexKind::Csi) => build_bcf(reader, min_shift),
        (VariantFormat::Bcf, IndexKind::Tbi) => Err(invalid("TBI does not support BCF")),
    }
}

fn build_vcf<R>(reader: R, kind: IndexKind, requested_min_shift: u8) -> io::Result<Built>
where
    R: bgzf::io::BufRead,
{
    let mut reader = noodles_vcf::io::Reader::new(reader);
    let header = reader.read_header()?;
    let max_length = header
        .contigs()
        .values()
        .filter_map(|contig| contig.length())
        .max();
    let (min_shift, depth) = match kind {
        IndexKind::Csi => {
            let initial_depth = 31u8.saturating_sub(requested_min_shift).div_ceil(3);
            settings(requested_min_shift, initial_depth, max_length)
        }
        IndexKind::Tbi => (14, 5),
    };

    let mut reader = reader.into_inner();
    let mut indexer = VcfIndexer::new(kind, min_shift, depth);
    let mut names = Vec::<BString>::new();
    let mut name_ids = HashMap::<Vec<u8>, usize>::new();
    let mut line = Vec::new();
    let mut start_offset = reader.virtual_position();
    let mut records = 0u64;
    let mut order = CoordinateOrder::default();

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        trim_line_ending(&mut line);
        records += 1;
        let end_offset = reader.virtual_position();
        let interval = parse_interval(&line).map_err(|error| {
            io::Error::new(error.kind(), format!("VCF record {records}: {error}"))
        })?;
        let reference_sequence_id = match name_ids.get(interval.name).copied() {
            Some(id) => id,
            None => {
                let id = names.len();
                name_ids.insert(interval.name.to_vec(), id);
                names.push(BString::from(interval.name));
                indexer.add_reference_sequence();
                id
            }
        };
        order.check(reference_sequence_id, interval.start)?;
        indexer.add_record(
            reference_sequence_id,
            std::str::from_utf8(interval.name)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?,
            interval.start,
            interval.end,
            Chunk::new(start_offset, end_offset),
        )?;
        start_offset = end_offset;
    }

    let reference_sequences = names.len();
    let index = indexer.finish(names)?;
    Ok(Built {
        index,
        format: VariantFormat::Vcf,
        min_shift,
        depth,
        records,
        reference_sequences,
    })
}

fn build_bcf<R>(reader: R, requested_min_shift: u8) -> io::Result<Built>
where
    R: bgzf::io::BufRead,
{
    let mut bcf_reader = bcf::io::Reader::from(reader);
    let header = bcf_reader.read_header()?;
    let max_length = header
        .contigs()
        .values()
        .filter_map(|contig| contig.length())
        .max()
        .or(Some(i32::MAX as usize));
    let (min_shift, depth) = settings(requested_min_shift, 0, max_length);
    let reference_sequence_count = header.contigs().len();
    let mut indexer = binning_index::Indexer::<BinnedIndex>::new(min_shift, depth);
    let mut linear = (0..reference_sequence_count)
        .map(|_| LinearIndex::new(min_shift, depth))
        .collect::<Vec<_>>();
    let mut reader = bcf_reader.into_inner();
    let mut shared = Vec::new();
    let mut samples = Vec::new();
    let mut start_offset = reader.virtual_position();
    let mut records = 0u64;
    let mut order = CoordinateOrder::default();

    while let Some(record) =
        bcf_record::read(&mut reader, &mut shared, &mut samples, header.string_maps()).map_err(
            |error| io::Error::new(error.kind(), format!("BCF record {}: {error}", records + 1)),
        )?
    {
        records += 1;
        let end_offset = reader.virtual_position();
        let reference_sequence_id = usize::try_from(record.reference_sequence_id)
            .map_err(|_| invalid("BCF reference sequence ID is negative"))?;
        let start_number = if record.position == u32::MAX {
            1
        } else {
            usize::try_from(record.position)
                .map_err(|_| invalid("BCF position is outside the supported range"))?
                .checked_add(1)
                .ok_or_else(|| invalid("BCF position is outside the supported range"))?
        };
        let span = usize::try_from(record.span)
            .map_err(|_| invalid("BCF record span must be positive"))?;
        if span == 0 {
            return Err(invalid("BCF record span must be positive"));
        }
        let start = Position::try_from(start_number)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        let end = start
            .checked_add(span - 1)
            .ok_or_else(|| invalid("BCF record end overflows the supported range"))?;
        order.check(reference_sequence_id, start)?;
        let chunk = Chunk::new(start_offset, end_offset);
        let linear_index = linear
            .get_mut(reference_sequence_id)
            .ok_or_else(|| invalid("BCF reference sequence ID is outside the header"))?;
        linear_index.insert(start, end, start_offset);
        indexer
            .add_record(Some((reference_sequence_id, start, end, true)), chunk)
            .map_err(|error| {
                io::Error::new(error.kind(), format!("BCF record {records}: {error}"))
            })?;
        start_offset = end_offset;
    }

    let index = apply_loffsets(indexer.build(reference_sequence_count), &mut linear, None);
    Ok(Built {
        index: Index::Csi(index),
        format: VariantFormat::Bcf,
        min_shift,
        depth,
        records,
        reference_sequences: reference_sequence_count,
    })
}

#[derive(Default)]
struct CoordinateOrder {
    reference_sequence_id: Option<usize>,
    start: Option<Position>,
}

impl CoordinateOrder {
    fn check(&mut self, reference_sequence_id: usize, start: Position) -> io::Result<()> {
        if let Some(current) = self.reference_sequence_id {
            if reference_sequence_id < current {
                return Err(invalid("reference sequence blocks are not contiguous"));
            }
            if reference_sequence_id == current
                && self.start.is_some_and(|previous| start < previous)
            {
                return Err(invalid("positions are not sorted"));
            }
        }
        if self.reference_sequence_id != Some(reference_sequence_id) {
            self.reference_sequence_id = Some(reference_sequence_id);
            self.start = None;
        }
        self.start = Some(start);
        Ok(())
    }
}

enum VcfIndexer {
    Csi {
        indexer: binning_index::Indexer<BinnedIndex>,
        linear: Vec<LinearIndex>,
        min_shift: u8,
        depth: u8,
    },
    Tbi(tabix::index::Indexer),
}

impl VcfIndexer {
    fn new(kind: IndexKind, min_shift: u8, depth: u8) -> Self {
        match kind {
            IndexKind::Csi => Self::Csi {
                indexer: binning_index::Indexer::new(min_shift, depth),
                linear: Vec::new(),
                min_shift,
                depth,
            },
            IndexKind::Tbi => {
                let mut indexer = tabix::index::Indexer::default();
                indexer.set_header(HeaderBuilder::vcf().build());
                Self::Tbi(indexer)
            }
        }
    }

    fn add_reference_sequence(&mut self) {
        if let Self::Csi {
            linear,
            min_shift,
            depth,
            ..
        } = self
        {
            linear.push(LinearIndex::new(*min_shift, *depth));
        }
    }

    fn add_record(
        &mut self,
        reference_sequence_id: usize,
        name: &str,
        start: Position,
        end: Position,
        chunk: Chunk,
    ) -> io::Result<()> {
        match self {
            Self::Csi {
                indexer, linear, ..
            } => {
                linear[reference_sequence_id].insert(start, end, chunk.start());
                indexer.add_record(Some((reference_sequence_id, start, end, true)), chunk)
            }
            Self::Tbi(indexer) => indexer.add_record(name, start, end, chunk),
        }
    }

    fn finish(self, names: Vec<BString>) -> io::Result<Index> {
        match self {
            Self::Csi {
                indexer,
                mut linear,
                ..
            } => {
                let reference_sequence_count = names.len();
                let header = HeaderBuilder::vcf()
                    .set_reference_sequence_names(
                        names.into_iter().collect::<ReferenceSequenceNames>(),
                    )
                    .build();
                Ok(Index::Csi(apply_loffsets(
                    indexer.build(reference_sequence_count),
                    &mut linear,
                    Some(header),
                )))
            }
            Self::Tbi(indexer) => Ok(Index::Tbi(indexer.build())),
        }
    }
}

fn write_index(output: &mut File, index: &Index) -> Result<()> {
    match index {
        Index::Csi(index) => {
            let mut writer = noodles_csi::io::Writer::new(output);
            writer.write_index(index).map_err(RsomicsError::Io)?;
            writer.into_inner().try_finish().map_err(RsomicsError::Io)?;
        }
        Index::Tbi(index) => {
            let mut writer = tabix::io::Writer::new(output);
            writer.write_index(index).map_err(RsomicsError::Io)?;
            writer.try_finish().map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}

fn source_name(path: &Path) -> String {
    if path == Path::new("-") {
        "standard input".to_owned()
    } else {
        path.display().to_string()
    }
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}
