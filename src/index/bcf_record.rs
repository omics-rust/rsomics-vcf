use std::io::{self, Read};

use noodles_vcf::header::string_maps::{StringMap, StringMaps};

pub(super) struct Site {
    pub reference_sequence_id: i32,
    pub position: u32,
    pub span: i32,
}

pub(super) fn read(
    reader: &mut impl Read,
    shared: &mut Vec<u8>,
    samples: &mut Vec<u8>,
    maps: &StringMaps,
) -> io::Result<Option<Site>> {
    let mut lengths = [0; 8];
    if !read_exact_or_eof(reader, &mut lengths)? {
        return Ok(None);
    }

    let shared_length = usize::try_from(u32::from_le_bytes(lengths[..4].try_into().unwrap()))
        .map_err(invalid_value)?;
    let sample_length = usize::try_from(u32::from_le_bytes(lengths[4..].try_into().unwrap()))
        .map_err(invalid_value)?;
    if shared_length < 24 {
        return Err(invalid("shared block is shorter than 24 bytes"));
    }

    resize(shared, shared_length)?;
    resize(samples, sample_length)?;
    reader.read_exact(shared)?;
    reader.read_exact(samples)?;

    let reference_sequence_id = i32::from_le_bytes(shared[..4].try_into().unwrap());
    let position = u32::from_le_bytes(shared[4..8].try_into().unwrap());
    let span = i32::from_le_bytes(shared[8..12].try_into().unwrap());
    let info_count = usize::from(u16::from_le_bytes(shared[16..18].try_into().unwrap()));
    let allele_count = usize::from(u16::from_le_bytes(shared[18..20].try_into().unwrap()));
    let sample_count = usize::try_from(u32::from_le_bytes([shared[20], shared[21], shared[22], 0]))
        .map_err(invalid_value)?;
    let format_count = usize::from(shared[23]);

    let reference_index =
        usize::try_from(reference_sequence_id).map_err(|_| invalid("contig ID is negative"))?;
    if maps.contigs().get_index(reference_index).is_none() {
        return Err(invalid("contig ID is outside the header dictionary"));
    }
    if span <= 0 {
        return Err(invalid("record span must be positive"));
    }
    if allele_count == 0 {
        return Err(invalid("record has no REF allele"));
    }

    validate_shared(&shared[24..], allele_count, info_count, maps.strings())?;
    validate_samples(samples, sample_count, format_count, maps.strings())?;

    Ok(Some(Site {
        reference_sequence_id,
        position,
        span,
    }))
}

fn validate_shared(
    mut input: &[u8],
    allele_count: usize,
    info_count: usize,
    strings: &StringMap,
) -> io::Result<()> {
    consume_string(&mut input, "ID")?;
    for _ in 0..allele_count {
        consume_string(&mut input, "REF/ALT")?;
    }

    let filter = read_type(&mut input)?;
    if filter.length > 0 {
        if !filter.kind.is_integer() {
            return Err(invalid("FILTER is not integer encoded"));
        }
        for _ in 0..filter.length {
            validate_dictionary_id(read_integer(&mut input, filter.kind)?, strings, "FILTER")?;
        }
    }

    for _ in 0..info_count {
        let key = read_typed_integer(&mut input)?;
        validate_dictionary_id(key, strings, "INFO")?;
        consume_value(&mut input, 1, "INFO")?;
    }

    Ok(())
}

fn validate_samples(
    mut input: &[u8],
    sample_count: usize,
    format_count: usize,
    strings: &StringMap,
) -> io::Result<()> {
    for _ in 0..format_count {
        let key = read_typed_integer(&mut input)?;
        validate_dictionary_id(key, strings, "FORMAT")?;
        consume_value(&mut input, sample_count, "FORMAT")?;
    }
    Ok(())
}

fn consume_string(input: &mut &[u8], field: &str) -> io::Result<()> {
    let encoded = read_type(input)?;
    if encoded.kind != Kind::String {
        return Err(invalid(format!("{field} is not string encoded")));
    }
    take(input, encoded.length)?;
    Ok(())
}

fn consume_value(input: &mut &[u8], multiplier: usize, field: &str) -> io::Result<()> {
    let encoded = read_type(input)?;
    if !encoded.kind.is_value() || encoded.kind == Kind::Null && encoded.length > 0 {
        return Err(invalid(format!("{field} has an invalid value type")));
    }
    let length = encoded
        .length
        .checked_mul(multiplier)
        .and_then(|length| length.checked_mul(encoded.kind.width()))
        .ok_or_else(|| invalid(format!("{field} value length overflows")))?;
    take(input, length)?;
    Ok(())
}

fn validate_dictionary_id(id: i32, strings: &StringMap, field: &str) -> io::Result<()> {
    let index = usize::try_from(id).map_err(|_| invalid(format!("{field} ID is negative")))?;
    if strings.get_index(index).is_none() {
        return Err(invalid(format!(
            "{field} ID is outside the header dictionary"
        )));
    }
    Ok(())
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum Kind {
    Null,
    Int8,
    Int16,
    Int32,
    Float,
    String,
}

impl Kind {
    fn parse(value: u8) -> io::Result<Self> {
        match value {
            0 => Ok(Self::Null),
            1 => Ok(Self::Int8),
            2 => Ok(Self::Int16),
            3 => Ok(Self::Int32),
            5 => Ok(Self::Float),
            7 => Ok(Self::String),
            _ => Err(invalid("value uses an invalid BCF type")),
        }
    }

    fn width(self) -> usize {
        match self {
            Self::Null => 0,
            Self::Int8 | Self::String => 1,
            Self::Int16 => 2,
            Self::Int32 | Self::Float => 4,
        }
    }

    fn is_integer(self) -> bool {
        matches!(self, Self::Int8 | Self::Int16 | Self::Int32)
    }

    fn is_value(self) -> bool {
        matches!(
            self,
            Self::Null | Self::Int8 | Self::Int16 | Self::Int32 | Self::Float | Self::String
        )
    }
}

struct EncodedType {
    length: usize,
    kind: Kind,
}

fn read_type(input: &mut &[u8]) -> io::Result<EncodedType> {
    let descriptor = take(input, 1)?[0];
    let kind = Kind::parse(descriptor & 0x0f)?;
    let mut length = usize::from(descriptor >> 4);
    if length == 15 {
        length = usize::try_from(read_typed_integer(input)?)
            .map_err(|_| invalid("value length is negative"))?;
    }
    Ok(EncodedType { length, kind })
}

fn read_typed_integer(input: &mut &[u8]) -> io::Result<i32> {
    let encoded = read_type(input)?;
    if encoded.length != 1 || !encoded.kind.is_integer() {
        return Err(invalid("dictionary ID is not an integer scalar"));
    }
    read_integer(input, encoded.kind)
}

fn read_integer(input: &mut &[u8], kind: Kind) -> io::Result<i32> {
    match kind {
        Kind::Int8 => Ok(i32::from(i8::from_le_bytes(
            take(input, 1)?.try_into().unwrap(),
        ))),
        Kind::Int16 => Ok(i32::from(i16::from_le_bytes(
            take(input, 2)?.try_into().unwrap(),
        ))),
        Kind::Int32 => Ok(i32::from_le_bytes(take(input, 4)?.try_into().unwrap())),
        _ => Err(invalid("value is not integer encoded")),
    }
}

fn take<'a>(input: &mut &'a [u8], length: usize) -> io::Result<&'a [u8]> {
    if input.len() < length {
        return Err(invalid("record block is malformed or too short"));
    }
    let (value, remaining) = input.split_at(length);
    *input = remaining;
    Ok(value)
}

fn resize(buffer: &mut Vec<u8>, length: usize) -> io::Result<()> {
    if length > buffer.len() {
        buffer
            .try_reserve(length - buffer.len())
            .map_err(invalid_value)?;
    }
    buffer.resize(length, 0);
    Ok(())
}

fn read_exact_or_eof(reader: &mut impl Read, buffer: &mut [u8]) -> io::Result<bool> {
    let Some((first, rest)) = buffer.split_first_mut() else {
        return Ok(true);
    };
    loop {
        match reader.read(std::slice::from_mut(first)) {
            Ok(0) => return Ok(false),
            Ok(_) => break,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    reader.read_exact(rest)?;
    Ok(true)
}

fn invalid_value(error: impl std::error::Error + Send + Sync + 'static) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error)
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}
