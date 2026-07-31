use std::io;

use noodles_core::Position;

pub(super) struct Interval<'a> {
    pub name: &'a [u8],
    pub start: Position,
    pub end: Position,
}

pub(super) fn parse_interval(line: &[u8]) -> io::Result<Interval<'_>> {
    let mut fields = line.split(|byte| *byte == b'\t');
    let name = next_field(&mut fields, "CHROM")?;
    if name.is_empty() {
        return Err(invalid("CHROM is empty"));
    }
    let raw_position = parse_u64(next_field(&mut fields, "POS")?, "POS")?;
    let start_number = usize::try_from(raw_position.max(1))
        .map_err(|_| invalid("POS is outside the supported range"))?;
    let start = Position::try_from(start_number)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    next_field(&mut fields, "ID")?;
    let reference = next_field(&mut fields, "REF")?;
    if reference.is_empty() || reference == b"." {
        return Err(invalid("REF is empty"));
    }
    let alternates = next_field(&mut fields, "ALT")?;
    next_field(&mut fields, "QUAL")?;
    next_field(&mut fields, "FILTER")?;
    let info = next_field(&mut fields, "INFO")?;
    let format = fields.next();
    let samples = fields;

    let mut span = reference.len();
    if let Some(value) = info_value(info, b"SVLEN") {
        span = span.max(max_reference_svlen(alternates, value)?);
    }
    if has_reference_block(alternates)
        && let Some(format) = format
        && let Some(index) = format
            .split(|byte| *byte == b':')
            .position(|key| key == b"LEN")
    {
        span = span.max(max_sample_length(samples, index)?);
    }
    let mut end_number = start_number
        .checked_add(span.saturating_sub(1))
        .ok_or_else(|| invalid("variant end overflows the supported range"))?;
    if let Some(value) = info_value(info, b"END") {
        let end = parse_u64(value, "INFO/END")?;
        if end > raw_position {
            end_number = end_number.max(
                usize::try_from(end)
                    .map_err(|_| invalid("INFO/END is outside the supported range"))?,
            );
        }
    }
    let end = Position::try_from(end_number.max(start_number))
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;

    Ok(Interval { name, start, end })
}

pub(super) fn trim_line_ending(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

fn next_field<'a>(fields: &mut impl Iterator<Item = &'a [u8]>, name: &str) -> io::Result<&'a [u8]> {
    fields
        .next()
        .ok_or_else(|| invalid(format!("record is missing {name}")))
}

fn info_value<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    info.split(|byte| *byte == b';').find_map(|field| {
        let separator = field.iter().position(|byte| *byte == b'=')?;
        let candidate = &field[..separator];
        let value = &field[separator + 1..];
        (candidate == key).then_some(value)
    })
}

fn max_reference_svlen(alternates: &[u8], values: &[u8]) -> io::Result<usize> {
    let mut max_length = 0usize;
    for (alternate, value) in alternates
        .split(|byte| *byte == b',')
        .zip(values.split(|byte| *byte == b','))
    {
        if uses_svlen(alternate) && value != b"." {
            let value = parse_i64(value, "INFO/SVLEN")?;
            let length = usize::try_from(value.unsigned_abs())
                .map_err(|_| invalid("INFO/SVLEN is outside the supported range"))?;
            max_length = max_length.max(length);
        }
    }
    Ok(max_length)
}

fn uses_svlen(alternate: &[u8]) -> bool {
    alternate.len() >= 5
        && alternate.last() == Some(&b'>')
        && alternate
            .get(4)
            .is_some_and(|byte| matches!(byte, b'>' | b':'))
        && matches!(
            alternate.get(..4),
            Some(b"<CNV" | b"<DEL" | b"<DUP" | b"<INV")
        )
}

fn has_reference_block(alternates: &[u8]) -> bool {
    alternates
        .split(|byte| *byte == b',')
        .any(|alternate| matches!(alternate, b"<*>" | b"<NON_REF>"))
}

fn max_sample_length<'a>(
    samples: impl Iterator<Item = &'a [u8]>,
    index: usize,
) -> io::Result<usize> {
    let mut max_length = 0usize;
    for sample in samples {
        let Some(value) = sample.split(|byte| *byte == b':').nth(index) else {
            continue;
        };
        if value != b"." && !value.is_empty() {
            let length = parse_u64(value, "FORMAT/LEN")?;
            max_length = max_length.max(
                usize::try_from(length)
                    .map_err(|_| invalid("FORMAT/LEN is outside the supported range"))?,
            );
        }
    }
    Ok(max_length)
}

fn parse_u64(value: &[u8], field: &str) -> io::Result<u64> {
    let value = std::str::from_utf8(value)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    value
        .parse()
        .map_err(|_| invalid(format!("{field} is not an unsigned integer")))
}

fn parse_i64(value: &[u8], field: &str) -> io::Result<i64> {
    let value = std::str::from_utf8(value)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    value
        .parse()
        .map_err(|_| invalid(format!("{field} is not an integer")))
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_reference_spans() {
        let interval =
            parse_interval(b"chr1\t100\t.\tN\t<DEL>,A\t.\tPASS\tEND=500;SVLEN=-350,1").unwrap();
        assert_eq!(interval.name, b"chr1");
        assert_eq!(usize::from(interval.start), 100);
        assert_eq!(usize::from(interval.end), 500);
    }

    #[test]
    fn parses_reference_block_lengths() {
        let interval =
            parse_interval(b"chr1\t100\t.\tN\t<NON_REF>\t.\tPASS\t.\tGT:LEN\t0/0:81\t0/0:20")
                .unwrap();
        assert_eq!(usize::from(interval.end), 180);
    }
}
