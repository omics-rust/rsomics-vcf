use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;
use noodles_bcf as bcf;
use noodles_vcf as vcf;
use rsomics_common::{Context, Result, RsomicsError};

use super::HeaderTypes;

type Source = Box<dyn Read>;
type Buffered = BufReader<Source>;
type Compressed = MultiGzDecoder<Buffered>;

enum Inner {
    Bcf(bcf::io::Reader<Compressed>),
    BcfRaw(bcf::io::Reader<Buffered>),
    Vcf(Buffered),
    VcfGz(BufReader<Compressed>),
}

pub(crate) struct Reader {
    inner: Inner,
    source: String,
}

impl Reader {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let source = if path == Path::new("-") {
            "standard input".to_owned()
        } else {
            path.display().to_string()
        };
        let input: Source = if path == Path::new("-") {
            Box::new(io::stdin())
        } else {
            Box::new(
                File::open(path)
                    .rs_with_context(|| format!("opening variant input {}", path.display()))?,
            )
        };
        let mut input = BufReader::new(input);
        let compression = is_gzip(&mut input)
            .map_err(|error| invalid(&source, "detecting compression", error))?;
        let bcf = is_bcf(&mut input, compression)
            .map_err(|error| invalid(&source, "detecting variant format", error))?;
        let inner = match (bcf, compression) {
            (true, true) => Inner::Bcf(bcf::io::Reader::from(MultiGzDecoder::new(input))),
            (true, false) => Inner::BcfRaw(bcf::io::Reader::from(input)),
            (false, true) => Inner::VcfGz(BufReader::new(MultiGzDecoder::new(input))),
            (false, false) => Inner::Vcf(input),
        };

        Ok(Self { inner, source })
    }

    pub(crate) fn read_header(&mut self) -> Result<(vcf::Header, Vec<u8>, HeaderTypes)> {
        let raw = self.read_raw_header()?;
        let header = self.parse_header(&raw)?;
        let output = canonical_header(&raw)
            .map_err(|error| invalid(&self.source, "formatting VCF header", error))?;
        let types = HeaderTypes::parse(&output)
            .map_err(|error| invalid(&self.source, "reading VCF schema", error))?;

        Ok((header, output, types))
    }

    pub(crate) fn read_raw_header(&mut self) -> Result<Vec<u8>> {
        match &mut self.inner {
            Inner::Bcf(reader) => read_bcf_header(reader),
            Inner::BcfRaw(reader) => read_bcf_header(reader),
            Inner::Vcf(reader) => read_vcf_header(reader),
            Inner::VcfGz(reader) => read_vcf_header(reader),
        }
        .map_err(|error| invalid(&self.source, "reading VCF header", error))
    }

    pub(crate) fn parse_header(&self, raw: &[u8]) -> Result<vcf::Header> {
        let text = std::str::from_utf8(raw).map_err(|error| {
            invalid(
                &self.source,
                "decoding VCF header",
                io::Error::new(io::ErrorKind::InvalidData, error),
            )
        })?;
        let mut header: vcf::Header = text.parse().map_err(|error| {
            invalid(
                &self.source,
                "parsing VCF header",
                io::Error::new(io::ErrorKind::InvalidData, error),
            )
        })?;
        if !self.is_text() {
            *header.string_maps_mut() = text.parse().map_err(|error| {
                invalid(
                    &self.source,
                    "reading BCF string maps",
                    io::Error::new(io::ErrorKind::InvalidData, error),
                )
            })?;
        }
        Ok(header)
    }

    pub(crate) fn is_text(&self) -> bool {
        matches!(self.inner, Inner::Vcf(_) | Inner::VcfGz(_))
    }

    pub(crate) fn read_text_record(
        &mut self,
        record: &mut Vec<u8>,
        number: usize,
    ) -> Result<usize> {
        self.read_text_record_with_termination(record, number)
            .map(|(read, _)| read)
    }

    pub(crate) fn read_text_record_with_termination(
        &mut self,
        record: &mut Vec<u8>,
        number: usize,
    ) -> Result<(usize, bool)> {
        let result = match &mut self.inner {
            Inner::Vcf(reader) => read_line_with_termination(reader, record),
            Inner::VcfGz(reader) => read_line_with_termination(reader, record),
            _ => unreachable!(),
        };
        result.map_err(|error| {
            invalid(
                &self.source,
                &format!("reading variant record {number}"),
                error,
            )
        })
    }

    pub(crate) fn read_bcf_record(
        &mut self,
        record: &mut bcf::Record,
        number: usize,
    ) -> Result<usize> {
        let result = match &mut self.inner {
            Inner::Bcf(reader) => reader.read_record(record),
            Inner::BcfRaw(reader) => reader.read_record(record),
            _ => unreachable!(),
        };
        result.map_err(|error| {
            invalid(
                &self.source,
                &format!("reading variant record {number}"),
                error,
            )
        })
    }
}

fn is_gzip(reader: &mut impl BufRead) -> io::Result<bool> {
    Ok(reader
        .fill_buf()?
        .get(..2)
        .is_some_and(|magic| magic == [0x1f, 0x8b]))
}

fn is_bcf(reader: &mut impl BufRead, compressed: bool) -> io::Result<bool> {
    let source = reader.fill_buf()?;
    if compressed {
        let mut decoder = MultiGzDecoder::new(source);
        let mut magic = [0; 3];
        decoder.read_exact(&mut magic)?;
        Ok(magic == *b"BCF")
    } else {
        Ok(source.get(..3).is_some_and(|magic| magic == b"BCF"))
    }
}

fn read_vcf_header<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: BufRead,
{
    let mut raw = Vec::new();
    loop {
        let source = reader.fill_buf()?;
        if source.first() != Some(&b'#') {
            break;
        }
        reader.read_until(b'\n', &mut raw)?;
    }
    Ok(raw)
}

fn read_bcf_header<R>(reader: &mut bcf::io::Reader<R>) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut reader = reader.header_reader();
    if reader.read_magic_number()? != *b"BCF" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF magic number",
        ));
    }
    let version = reader.read_format_version()?;
    if version != (2, 2) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported BCF version {}.{}", version.0, version.1),
        ));
    }

    let mut raw_reader = reader.raw_vcf_header_reader()?;
    let mut raw = Vec::new();
    raw_reader.read_to_end(&mut raw)?;
    raw_reader.discard_to_end()?;
    Ok(raw)
}

fn read_line_with_termination(
    reader: &mut impl BufRead,
    record: &mut Vec<u8>,
) -> io::Result<(usize, bool)> {
    record.clear();
    let read = reader.read_until(b'\n', record)?;
    let terminated = record.last() == Some(&b'\n');
    if terminated {
        record.pop();
        if record.last() == Some(&b'\r') {
            record.pop();
        }
    }
    Ok((read, terminated))
}

fn canonical_header(raw: &[u8]) -> io::Result<Vec<u8>> {
    let mut lines = Vec::new();
    for raw_line in raw.split(|byte| *byte == b'\n') {
        let line = raw_line.strip_suffix(b"\r").unwrap_or(raw_line);
        if line.is_empty() || is_pass_filter(line) {
            continue;
        }
        lines.push(strip_idx(line));
    }
    if lines
        .first()
        .is_none_or(|line| !line.starts_with(b"##fileformat="))
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "VCF header must start with ##fileformat",
        ));
    }

    let mut output = Vec::with_capacity(raw.len() + 64);
    write_line(&mut output, &lines[0]);
    write_line(
        &mut output,
        br#"##FILTER=<ID=PASS,Description="All filters passed">"#,
    );
    for line in &lines[1..] {
        write_line(&mut output, line);
    }
    Ok(output)
}

fn is_pass_filter(line: &[u8]) -> bool {
    let Some(body) = line
        .strip_prefix(b"##FILTER=<")
        .and_then(|line| line.strip_suffix(b">"))
    else {
        return false;
    };
    fields(body).any(|field| trim(field) == b"ID=PASS")
}

fn strip_idx(line: &[u8]) -> Vec<u8> {
    let Some(start) = line.windows(2).position(|window| window == b"=<") else {
        return line.to_vec();
    };
    let Some(body) = line
        .get(start + 2..)
        .and_then(|line| line.strip_suffix(b">"))
    else {
        return line.to_vec();
    };
    let kept: Vec<_> = fields(body)
        .filter(|field| !trim(field).starts_with(b"IDX="))
        .collect();
    if kept.len() == fields(body).count() {
        return line.to_vec();
    }

    let mut output = Vec::with_capacity(line.len());
    output.extend_from_slice(&line[..start + 2]);
    for (index, field) in kept.into_iter().enumerate() {
        if index > 0 {
            output.push(b',');
        }
        output.extend_from_slice(field);
    }
    output.push(b'>');
    output
}

fn fields(body: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut quoted = false;
    let mut escaped = false;
    body.split(move |byte| {
        if escaped {
            escaped = false;
            return false;
        }
        match *byte {
            b'\\' if quoted => {
                escaped = true;
                false
            }
            b'"' => {
                quoted = !quoted;
                false
            }
            b',' if !quoted => true,
            _ => false,
        }
    })
}

fn trim(mut field: &[u8]) -> &[u8] {
    while field.first().is_some_and(u8::is_ascii_whitespace) {
        field = &field[1..];
    }
    while field.last().is_some_and(u8::is_ascii_whitespace) {
        field = &field[..field.len() - 1];
    }
    field
}

fn write_line(output: &mut Vec<u8>, line: &[u8]) {
    output.extend_from_slice(line);
    output.push(b'\n');
}

fn invalid(source: &str, operation: &str, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{source}: {operation}: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_preserves_order_and_inserts_pass() {
        let raw = b"##fileformat=VCFv4.3\n##source=x\n##contig=<ID=chr1,length=9>\n#CHROM\tPOS\n";
        let actual = canonical_header(raw).unwrap();
        assert_eq!(
            actual,
            b"##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##source=x\n##contig=<ID=chr1,length=9>\n#CHROM\tPOS\n"
        );
    }

    #[test]
    fn header_replaces_pass_and_removes_bcf_indices() {
        let raw = b"##fileformat=VCFv4.3\n##FILTER=<Description=\"custom\",ID=PASS,IDX=0>\n##INFO=<ID=X,Number=1,Type=String,Description=\"a,b\",IDX=2>\n#CHROM\tPOS\n";
        let actual = canonical_header(raw).unwrap();
        assert_eq!(
            actual,
            b"##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##INFO=<ID=X,Number=1,Type=String,Description=\"a,b\">\n#CHROM\tPOS\n"
        );
    }
}
