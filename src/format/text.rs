use std::io;

use super::value::ValueType;

#[derive(Default)]
struct TypeMap(Vec<(Vec<u8>, ValueType)>);

impl TypeMap {
    fn insert(&mut self, id: &[u8], value_type: ValueType) {
        self.0.push((id.to_vec(), value_type));
    }

    fn get(&self, id: &[u8]) -> Option<ValueType> {
        self.0
            .iter()
            .find(|(key, _)| key == id)
            .map(|(_, value_type)| *value_type)
    }
}

#[derive(Default)]
pub(crate) struct HeaderTypes {
    info: TypeMap,
    format: TypeMap,
    samples: usize,
}

impl HeaderTypes {
    pub(super) fn parse(header: &[u8]) -> io::Result<Self> {
        let mut types = Self::default();
        let mut chrom = false;
        for line in header.split(|byte| *byte == b'\n') {
            if let Some(body) = line
                .strip_prefix(b"##INFO=<")
                .and_then(|line| line.strip_suffix(b">"))
            {
                absorb(body, &mut types.info);
            } else if let Some(body) = line
                .strip_prefix(b"##FORMAT=<")
                .and_then(|line| line.strip_suffix(b">"))
            {
                absorb(body, &mut types.format);
            } else if line.starts_with(b"#CHROM\t") {
                let columns: Vec<_> = line.split(|byte| *byte == b'\t').collect();
                if columns.len() < 8 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "#CHROM line has fewer than 8 columns",
                    ));
                }
                let fixed: [&[u8]; 8] = [
                    b"#CHROM", b"POS", b"ID", b"REF", b"ALT", b"QUAL", b"FILTER", b"INFO",
                ];
                if columns[..8] != fixed {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "#CHROM line has invalid fixed columns",
                    ));
                }
                if columns.len() == 9 || (columns.len() > 9 && columns[8] != b"FORMAT") {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "#CHROM sample columns require FORMAT and at least one sample",
                    ));
                }
                types.samples = columns.len().saturating_sub(9);
                chrom = true;
            }
        }
        if !chrom {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF header is missing the #CHROM line",
            ));
        }
        Ok(types)
    }

    pub(super) fn info(&self, id: &[u8]) -> Option<ValueType> {
        self.info.get(id)
    }

    pub(super) fn format(&self, id: &[u8]) -> Option<ValueType> {
        self.format.get(id)
    }

    pub(super) fn samples(&self) -> usize {
        self.samples
    }
}

fn absorb(body: &[u8], target: &mut TypeMap) {
    let mut id = None;
    let mut value_type = None;
    for field in fields(body) {
        if let Some(value) = field.strip_prefix(b"ID=") {
            id = Some(value);
        } else if let Some(value) = field.strip_prefix(b"Type=") {
            value_type = ValueType::parse(value);
        }
    }
    if let (Some(id), Some(value_type)) = (id, value_type) {
        target.insert(id, value_type);
    }
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
