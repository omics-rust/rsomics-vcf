use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use noodles_vcf::{
    self as vcf,
    variant::{RecordBuf, record::info::field::key, record_buf::info::field::Value},
};
use rsomics_common::{Context, Result, RsomicsError};

use super::columns::{BoundColumns, Column, ColumnSpec, MatchField, MatchLayout, SourceKind};
use crate::format::{Reader, RecordScratch};

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Payload {
    Variant(Box<RecordBuf>),
    Tabular(Vec<Option<Vec<u8>>>),
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct AnnotationRecord {
    pub(crate) contig: usize,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) serial: u64,
    pub(crate) reference: Option<String>,
    pub(crate) alternates: Vec<String>,
    pub(crate) id: Option<String>,
    pub(crate) info_end: Option<usize>,
    pub(crate) payload: Payload,
    pub(crate) zero_based: bool,
    contig_name: String,
}

enum SourceReader {
    Variant {
        reader: Box<Reader>,
        header: Box<vcf::Header>,
        scratch: Box<RecordScratch>,
        columns: BoundColumns,
    },
    Tabular {
        reader: Box<dyn BufRead>,
        line: Vec<u8>,
        columns: BoundColumns,
        bed: bool,
        line_number: u64,
    },
}

pub(crate) struct AnnotationSource {
    reader: SourceReader,
    pub(crate) next: Option<AnnotationRecord>,
    pub(crate) active: VecDeque<AnnotationRecord>,
    pub(super) contigs: HashMap<String, usize>,
    last_coordinate: Option<(usize, usize, u64)>,
    pub(super) last_target_coordinate: Option<(usize, usize)>,
}

impl std::fmt::Debug for AnnotationSource {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("AnnotationSource")
            .field("next", &self.next)
            .field("active", &self.active)
            .field("last_coordinate", &self.last_coordinate)
            .finish_non_exhaustive()
    }
}

impl AnnotationSource {
    pub(crate) fn open(path: &Path, target: &vcf::Header, spec: ColumnSpec) -> Result<Self> {
        if path == Path::new("-") {
            return Err(invalid(
                "annotation source cannot use standard input while the target stream is active",
            ));
        }

        let reader = if spec.match_layout().is_some() {
            let bed = is_bed(path);
            if bed && spec.match_layout() != Some(MatchLayout::Interval) {
                return Err(invalid(
                    "BED annotation columns require CHROM, FROM, and TO",
                ));
            }
            let columns = BoundColumns::bind(spec, SourceKind::Tabular, target, None)?;
            SourceReader::Tabular {
                reader: open_tabular(path)?,
                line: Vec::new(),
                columns,
                bed,
                line_number: 0,
            }
        } else {
            let mut reader = Reader::open(path)?;
            let (header, _, _) = reader.read_header()?;
            let columns = BoundColumns::bind(spec, SourceKind::Variant, target, Some(&header))?;
            SourceReader::Variant {
                reader: Box::new(reader),
                header: Box::new(header),
                scratch: Box::new(RecordScratch::default()),
                columns,
            }
        };
        let contigs = target
            .contigs()
            .keys()
            .enumerate()
            .map(|(index, name)| (name.clone(), index))
            .collect();
        let mut source = Self {
            reader,
            next: None,
            active: VecDeque::new(),
            contigs,
            last_coordinate: None,
            last_target_coordinate: None,
        };
        source.next = source.read_next()?;
        Ok(source)
    }

    pub(crate) fn read_record(&mut self) -> Result<Option<AnnotationRecord>> {
        let Some(record) = self.next.take() else {
            return Ok(None);
        };
        self.next = self.read_next()?;
        Ok(Some(record))
    }

    pub(crate) fn columns(&self) -> &BoundColumns {
        match &self.reader {
            SourceReader::Variant { columns, .. } | SourceReader::Tabular { columns, .. } => {
                columns
            }
        }
    }

    fn read_next(&mut self) -> Result<Option<AnnotationRecord>> {
        let serial = match self.last_coordinate {
            Some((_, _, serial)) => serial
                .checked_add(1)
                .ok_or_else(|| invalid("annotation record count exceeds u64"))?,
            None => 1,
        };
        let record = match &mut self.reader {
            SourceReader::Variant {
                reader,
                header,
                scratch,
                ..
            } => match reader.read_record(header, scratch, serial)? {
                Some(record) => Some(variant_record(record, serial)?),
                None => None,
            },
            SourceReader::Tabular {
                reader,
                line,
                columns,
                bed,
                line_number,
            } => read_tabular_record(reader.as_mut(), line, columns, *bed, line_number, serial)?,
        };
        let Some(mut record) = record else {
            return Ok(None);
        };
        record.contig = *self.contigs.get(record_contig(&record)).ok_or_else(|| {
            invalid(format!(
                "annotation record {} uses unknown contig {:?}",
                record.serial,
                record_contig(&record)
            ))
        })?;
        self.check_order(&record)?;
        self.last_coordinate = Some((record.contig, record.start, record.serial));
        Ok(Some(record))
    }

    fn check_order(&self, record: &AnnotationRecord) -> Result<()> {
        if let Some((contig, start, _)) = self.last_coordinate
            && (record.contig, record.start) < (contig, start)
        {
            return Err(invalid(format!(
                "annotation record {} is out of coordinate order",
                record.serial
            )));
        }
        Ok(())
    }
}

impl AnnotationRecord {
    pub(crate) fn inclusive_start(&self) -> usize {
        if self.zero_based {
            self.start + 1
        } else {
            self.start
        }
    }

    pub(crate) fn inclusive_end(&self) -> usize {
        self.end
    }
}

fn open_tabular(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .rs_with_context(|| format!("opening annotation source {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let compressed = reader
        .fill_buf()
        .rs_with_context(|| format!("reading annotation source {}", path.display()))?
        .get(..2)
        .is_some_and(|magic| magic == [0x1f, 0x8b]);
    if compressed {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

fn is_bed(path: &Path) -> bool {
    let path = path.to_string_lossy().to_ascii_lowercase();
    [".bed", ".bed.gz", ".bed.bgz", ".bed.bgzf"]
        .iter()
        .any(|suffix| path.ends_with(suffix))
}

fn read_tabular_record(
    reader: &mut dyn BufRead,
    line: &mut Vec<u8>,
    columns: &BoundColumns,
    bed: bool,
    line_number: &mut u64,
    serial: u64,
) -> Result<Option<AnnotationRecord>> {
    loop {
        line.clear();
        let read = reader.read_until(b'\n', line).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading annotation line {}: {error}",
                line_number.saturating_add(1)
            ))
        })?;
        if read == 0 {
            return Ok(None);
        }
        *line_number = line_number.saturating_add(1);
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.is_empty() || line.first() == Some(&b'#') {
            continue;
        }
        return parse_tabular_line(line, columns, bed, *line_number, serial).map(Some);
    }
}

fn parse_tabular_line(
    line: &[u8],
    columns: &BoundColumns,
    bed: bool,
    line_number: u64,
    serial: u64,
) -> Result<AnnotationRecord> {
    let fields: Vec<_> = line.split(|byte| *byte == b'\t').collect();
    let expected = columns.spec().fields().len();
    if fields.len() < expected {
        return Err(invalid(format!(
            "annotation line {line_number} is missing column {}",
            fields.len() + 1
        )));
    }

    let mut contig = None;
    let mut start = None;
    let mut end = None;
    let mut reference = None;
    let mut alternates = Vec::new();
    let mut id = None;
    let mut info_end = None;

    for (index, column) in columns.spec().fields().iter().enumerate() {
        let Column::Match(field) = column else {
            continue;
        };
        let value = fields[index];
        match field {
            MatchField::Chrom => contig = Some(text(value, line_number, "CHROM")?.to_owned()),
            MatchField::Pos => {
                let value = integer(value, line_number, "POS")?;
                start = Some(value);
                end = Some(value);
            }
            MatchField::From => start = Some(integer(value, line_number, "FROM")?),
            MatchField::To => end = Some(integer(value, line_number, "TO")?),
            MatchField::Ref => reference = optional_text(value, line_number, "REF")?,
            MatchField::Alt => {
                if value != b"." {
                    alternates = text(value, line_number, "ALT")?
                        .split(',')
                        .map(str::to_owned)
                        .collect();
                }
            }
            MatchField::Id => id = optional_text(value, line_number, "ID")?,
            MatchField::End => {
                if value != b"." {
                    info_end = Some(integer(value, line_number, "END")?);
                }
            }
            MatchField::Ignore => {}
        }
    }

    let start = start.expect("validated coordinate layout has a start");
    let end = end.expect("validated coordinate layout has an end");
    if !bed && (start == 0 || end == 0) {
        return Err(invalid(format!(
            "annotation line {line_number} uses zero in one-based coordinates"
        )));
    }
    if end < start {
        return Err(invalid(format!(
            "annotation line {line_number} end {end} precedes start {start}"
        )));
    }
    if bed && columns.spec().match_layout() == Some(MatchLayout::Interval) && end == start {
        return Err(invalid(format!(
            "annotation line {line_number} has an empty BED interval"
        )));
    }

    let contig_name = contig.expect("validated coordinate layout has CHROM");
    Ok(AnnotationRecord {
        contig: 0,
        start,
        end,
        serial,
        reference,
        alternates,
        id,
        info_end,
        payload: Payload::Tabular(
            fields
                .into_iter()
                .map(|value| (value != b".").then(|| value.to_vec()))
                .collect(),
        ),
        zero_based: bed,
        contig_name,
    })
}

pub(super) fn variant_record(record: RecordBuf, serial: u64) -> Result<AnnotationRecord> {
    let start = record.variant_start().map_or(0, usize::from);
    let info_end = match record.info().get(key::END_POSITION) {
        None => None,
        Some(Some(Value::Integer(value))) => Some(usize::try_from(*value).map_err(|_| {
            invalid(format!(
                "annotation record {serial} has invalid END value {value}"
            ))
        })?),
        Some(_) => {
            return Err(invalid(format!(
                "annotation record {serial} has invalid END value"
            )));
        }
    };
    let end = match info_end {
        Some(end) => end,
        None => start
            .checked_add(record.reference_bases().len().max(1) - 1)
            .ok_or_else(|| invalid(format!("annotation record {serial} span exceeds usize")))?,
    };
    if end < start {
        return Err(invalid(format!(
            "annotation record {serial} end {end} precedes start {start}"
        )));
    }
    let contig_name = record.reference_sequence_name().to_owned();
    let reference = Some(record.reference_bases().to_owned());
    let alternates = record.alternate_bases().as_ref().to_vec();
    let id = (!record.ids().as_ref().is_empty()).then(|| {
        record
            .ids()
            .as_ref()
            .iter()
            .cloned()
            .collect::<Vec<_>>()
            .join(";")
    });
    Ok(AnnotationRecord {
        contig: 0,
        start,
        end,
        serial,
        reference,
        alternates,
        id,
        info_end,
        payload: Payload::Variant(Box::new(record)),
        zero_based: false,
        contig_name,
    })
}

fn record_contig(record: &AnnotationRecord) -> &str {
    &record.contig_name
}

fn integer(value: &[u8], line: u64, field: &str) -> Result<usize> {
    text(value, line, field)?.parse().map_err(|_| {
        invalid(format!(
            "annotation line {line} has invalid {field} integer {:?}",
            String::from_utf8_lossy(value)
        ))
    })
}

fn optional_text(value: &[u8], line: u64, field: &str) -> Result<Option<String>> {
    if value == b"." {
        Ok(None)
    } else {
        text(value, line, field).map(|value| Some(value.to_owned()))
    }
}

fn text<'a>(value: &'a [u8], line: u64, field: &str) -> Result<&'a str> {
    std::str::from_utf8(value).map_err(|error| {
        invalid(format!(
            "annotation line {line} has invalid UTF-8 in {field}: {error}"
        ))
    })
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use std::fs::{self, File};
    use std::io::Write;
    use std::path::{Path, PathBuf};

    use flate2::Compression;
    use flate2::write::GzEncoder;
    use noodles_bgzf as bgzf;
    use noodles_vcf as vcf;

    use super::*;
    use crate::annotate::columns::ColumnSpec;
    use crate::format::{HeaderMode, OutputFormat, Reader, RecordScratch, Writer};

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn target_header() -> vcf::Header {
        HEADER.parse().unwrap()
    }

    fn write(path: &Path, content: &[u8]) {
        fs::write(path, content).unwrap();
    }

    fn read_one(path: &Path, columns: &str) -> AnnotationRecord {
        let spec = ColumnSpec::parse(columns).unwrap();
        let mut source = AnnotationSource::open(path, &target_header(), spec).unwrap();
        let record = source.read_record().unwrap().unwrap();
        assert!(source.read_record().unwrap().is_none());
        record
    }

    #[test]
    fn reads_plain_gzip_and_bgzf_tabular_sources() {
        let directory = tempfile::tempdir().unwrap();
        let plain = directory.path().join("regions.tsv");
        let gzip = directory.path().join("regions.tsv.gz");
        let bgzf = directory.path().join("regions.tsv.bgz");
        let content = b"# source\r\n\r\nchr1\t10\t20\tDB1\r\n";
        write(&plain, content);

        let mut writer = GzEncoder::new(File::create(&gzip).unwrap(), Compression::default());
        writer.write_all(content).unwrap();
        writer.finish().unwrap();

        let mut writer = bgzf::io::Writer::new(File::create(&bgzf).unwrap());
        writer.write_all(content).unwrap();
        writer.try_finish().unwrap();

        for path in [&plain, &gzip, &bgzf] {
            let record = read_one(path, "CHROM,FROM,TO,INFO/DB");
            assert_eq!((record.contig, record.start, record.end), (0, 10, 20));
            let Payload::Tabular(fields) = record.payload else {
                panic!("expected tabular payload")
            };
            assert_eq!(fields[3].as_deref(), Some(b"DB1".as_slice()));
        }
    }

    #[test]
    fn bed_and_tab_coordinates_remain_distinct() {
        let directory = tempfile::tempdir().unwrap();
        let bed = directory.path().join("regions.BED.GZ");
        let tab = directory.path().join("regions.tsv");
        let content = b"chr1\t9\t20\tDB\n";

        let mut writer = bgzf::io::Writer::new(File::create(&bed).unwrap());
        writer.write_all(content).unwrap();
        writer.try_finish().unwrap();
        write(&tab, b"chr1\t10\t20\tDB\n");

        let bed = read_one(&bed, "CHROM,FROM,TO,INFO/DB");
        let tab = read_one(&tab, "CHROM,FROM,TO,INFO/DB");
        assert_eq!((bed.start, bed.end), (9, 20));
        assert_eq!((tab.start, tab.end), (10, 20));
        assert!(bed.zero_based);
        assert!(!tab.zero_based);
        assert_eq!(bed.inclusive_start(), tab.inclusive_start());
    }

    #[test]
    fn reads_committed_plain_bed_asset_in_target_order() {
        let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/DB").unwrap();
        let mut source =
            AnnotationSource::open(&fixture("regions.bed"), &target_header(), spec).unwrap();
        let first = source.read_record().unwrap().unwrap();
        let second = source.read_record().unwrap().unwrap();

        assert_eq!((first.contig, first.start, first.end), (0, 50, 250));
        assert_eq!((second.contig, second.start, second.end), (1, 200, 300));
        assert!(first.zero_based && second.zero_based);
        assert!(source.read_record().unwrap().is_none());
    }

    #[test]
    fn extracts_tabular_match_fields_and_long_payloads() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("annotations.tsv");
        let long = "x".repeat(256 * 1024);
        write(
            &path,
            format!("chr2\t7\tA\tC,G\trs7\t12\t{long}\n").as_bytes(),
        );

        let record = read_one(&path, "CHROM,POS,REF,ALT,~ID,~INFO/END,INFO/NOTE");
        assert_eq!((record.contig, record.start, record.end), (1, 7, 7));
        assert_eq!(record.reference.as_deref(), Some("A"));
        assert_eq!(record.alternates, ["C", "G"]);
        assert_eq!(record.id.as_deref(), Some("rs7"));
        assert_eq!(record.info_end, Some(12));
        let Payload::Tabular(fields) = record.payload else {
            panic!("expected tabular payload")
        };
        assert_eq!(fields[6].as_ref().unwrap().len(), long.len());
    }

    #[test]
    fn rejects_malformed_and_unsorted_tabular_records() {
        let cases = [
            ("missing.tsv", "chr1\t10\n", "missing column 3"),
            ("integer.tsv", "chr1\tx\t20\tDB\n", "invalid FROM"),
            ("zero.tsv", "chr1\t0\t20\tDB\n", "one-based"),
            ("range.tsv", "chr1\t20\t10\tDB\n", "precedes start"),
            ("contig.tsv", "chr9\t10\t20\tDB\n", "unknown contig"),
            (
                "order.tsv",
                "chr1\t20\t30\tA\nchr1\t10\t15\tB\n",
                "out of coordinate order",
            ),
            (
                "contig-order.tsv",
                "chr2\t10\t20\tA\nchr1\t30\t40\tB\n",
                "out of coordinate order",
            ),
        ];
        let directory = tempfile::tempdir().unwrap();
        for (name, content, expected) in cases {
            let path = directory.path().join(name);
            write(&path, content.as_bytes());
            let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/DB").unwrap();
            let result = AnnotationSource::open(&path, &target_header(), spec)
                .and_then(|mut source| source.read_record());
            let error = result.unwrap_err();
            assert!(error.to_string().contains(expected), "{name}: {error}");
        }
    }

    #[test]
    fn rejects_empty_bed_intervals_and_standard_input() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("empty.bed");
        write(&path, b"chr1\t10\t10\tDB\n");
        let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/DB").unwrap();
        let error = AnnotationSource::open(&path, &target_header(), spec.clone()).unwrap_err();
        assert!(error.to_string().contains("empty BED interval"), "{error}");
        let error = AnnotationSource::open(Path::new("-"), &target_header(), spec).unwrap_err();
        assert!(error.to_string().contains("standard input"), "{error}");
    }

    #[test]
    fn reads_every_variant_encoding() {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("input.vcf");
        write(
            &input,
            format!("{HEADER}chr1\t10\trs1\tAC\tA,G\t.\tPASS\tEND=14\n").as_bytes(),
        );

        for (format, extension) in [
            (OutputFormat::Vcf, "plain.vcf"),
            (OutputFormat::VcfBgzf, "compressed.vcf.gz"),
            (OutputFormat::BcfRaw, "raw.bcf"),
            (OutputFormat::Bcf, "compressed.bcf"),
        ] {
            let path = directory.path().join(extension);
            transcode(&input, &path, format);
            let record = read_one(&path, "ID");
            assert_eq!((record.contig, record.start, record.end), (0, 10, 14));
            assert_eq!(record.reference.as_deref(), Some("AC"));
            assert_eq!(record.alternates, ["A", "G"]);
            assert_eq!(record.id.as_deref(), Some("rs1"));
            assert_eq!(record.info_end, Some(14));
            assert!(matches!(record.payload, Payload::Variant(_)));
        }
    }

    #[test]
    fn rejects_variant_coordinate_regressions_and_unknown_contigs() {
        let directory = tempfile::tempdir().unwrap();
        for (name, records, expected) in [
            (
                "order.vcf",
                "chr1\t20\t.\tA\tC\t.\tPASS\t.\nchr1\t10\t.\tA\tC\t.\tPASS\t.\n",
                "out of coordinate order",
            ),
            (
                "contig.vcf",
                "chr9\t10\t.\tA\tC\t.\tPASS\t.\n",
                "unknown contig",
            ),
            (
                "end.vcf",
                "chr1\t10\t.\tA\tC\t.\tPASS\tEND=-1\n",
                "invalid END",
            ),
        ] {
            let path = directory.path().join(name);
            write(&path, format!("{HEADER}{records}").as_bytes());
            let spec = ColumnSpec::parse("ID").unwrap();
            let result =
                AnnotationSource::open(&path, &target_header(), spec).and_then(|mut source| {
                    while source.read_record()?.is_some() {}
                    Ok(())
                });
            let error = result.unwrap_err();
            assert!(error.to_string().contains(expected), "{name}: {error}");
        }
    }

    #[test]
    fn truncated_compressed_sources_fail_loud() {
        let directory = tempfile::tempdir().unwrap();
        let tab = directory.path().join("truncated.tsv.gz");
        let mut writer = GzEncoder::new(File::create(&tab).unwrap(), Compression::default());
        writer.write_all(b"chr1\t10\t20\tDB\n").unwrap();
        writer.finish().unwrap();
        truncate(&tab, 4);
        let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/DB").unwrap();
        let result = AnnotationSource::open(&tab, &target_header(), spec)
            .and_then(|mut source| source.read_record());
        assert!(result.is_err());

        let input = directory.path().join("input.vcf");
        write(
            &input,
            format!("{HEADER}chr1\t10\t.\tA\tC\t.\tPASS\t.\n").as_bytes(),
        );
        for (format, name) in [
            (OutputFormat::VcfBgzf, "truncated.vcf.gz"),
            (OutputFormat::BcfRaw, "truncated.raw.bcf"),
            (OutputFormat::Bcf, "truncated.bcf"),
        ] {
            let path = directory.path().join(name);
            transcode(&input, &path, format);
            truncate(&path, 4);
            let spec = ColumnSpec::parse("ID").unwrap();
            let result = AnnotationSource::open(&path, &target_header(), spec)
                .and_then(|mut source| source.read_record());
            assert!(result.is_err(), "{name}");
        }
    }

    fn transcode(input: &Path, output: &Path, format: OutputFormat) {
        let mut reader = Reader::open(input).unwrap();
        let (header, _, _) = reader.read_header().unwrap();
        let mut scratch = RecordScratch::default();
        let mut records = Vec::new();
        let mut number = 1;
        while let Some(record) = reader.read_record(&header, &mut scratch, number).unwrap() {
            records.push(record);
            number += 1;
        }
        let mut writer = Writer::new(File::create(output).unwrap(), format);
        writer.write_header(&header, HeaderMode::Full).unwrap();
        for (index, record) in records.iter().enumerate() {
            writer
                .write_record(&header, record, index as u64 + 1)
                .unwrap();
        }
        writer.finish().unwrap();
    }

    fn truncate(path: &Path, amount: u64) {
        let length = fs::metadata(path).unwrap().len();
        File::options()
            .write(true)
            .open(path)
            .unwrap()
            .set_len(length - amount)
            .unwrap();
    }

    fn fixture(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/upstream/bcftools-annotate")
            .join(name)
    }
}
