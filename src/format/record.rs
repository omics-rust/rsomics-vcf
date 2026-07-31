use super::HeaderTypes;
use super::value::{ValueType, reformat_into};

pub(crate) fn reformat_record(
    line: &[u8],
    header: &HeaderTypes,
    output: &mut Vec<u8>,
) -> Result<(), &'static str> {
    output.clear();
    let mut fields = line.splitn(9, |byte| *byte == b'\t');
    let (
        Some(chrom),
        Some(position),
        Some(ids),
        Some(reference),
        Some(alternates),
        Some(quality),
        Some(filters),
        Some(info),
    ) = (
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
    )
    else {
        return Err("expected at least 8 tab-delimited fields");
    };
    if chrom.is_empty() || reference.is_empty() || alternates.is_empty() {
        return Err("CHROM, REF, and ALT must not be empty");
    }

    output.extend_from_slice(chrom);
    output.push(b'\t');
    write_position(output, position)?;
    output.push(b'\t');
    output.extend_from_slice(ids);
    output.push(b'\t');
    output.extend_from_slice(reference);
    output.push(b'\t');
    output.extend_from_slice(alternates);
    output.push(b'\t');
    write_quality(output, quality)?;
    output.push(b'\t');
    output.extend_from_slice(filters);
    output.push(b'\t');
    write_info(output, info, header)?;

    match fields.next() {
        Some(samples) => write_samples(output, samples, header)?,
        None if header.samples() > 0 => {
            return Err("record is missing FORMAT and sample columns");
        }
        None => {}
    }
    Ok(())
}

fn write_position(output: &mut Vec<u8>, raw: &[u8]) -> Result<(), &'static str> {
    let valid = match unsigned_digits(raw) {
        Some(digits) if digits.len() <= 9 => digits.iter().any(|digit| *digit != b'0'),
        _ => std::str::from_utf8(raw)
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
            .is_some_and(|value| (1..=i32::MAX as u64).contains(&value)),
    };
    if !valid {
        return Err("POS must be an integer between 1 and 2147483647");
    }
    reformat_into(output, raw, ValueType::Integer);
    Ok(())
}

fn write_quality(output: &mut Vec<u8>, raw: &[u8]) -> Result<(), &'static str> {
    if raw == b"." {
        output.push(b'.');
    } else if valid_float(raw) {
        reformat_into(output, raw, ValueType::Float);
    } else {
        return Err("QUAL is not a floating-point value");
    }
    Ok(())
}

fn write_info(output: &mut Vec<u8>, raw: &[u8], header: &HeaderTypes) -> Result<(), &'static str> {
    if raw == b"." {
        output.push(b'.');
        return Ok(());
    }
    for (index, field) in raw.split(|byte| *byte == b';').enumerate() {
        if index > 0 {
            output.push(b';');
        }
        match split_once(field, b'=') {
            Some((key, value)) => {
                output.extend_from_slice(key);
                output.push(b'=');
                match header.info(key) {
                    Some(value_type) => write_typed(output, value, value_type)?,
                    None => output.extend_from_slice(value),
                }
            }
            None => output.extend_from_slice(field),
        }
    }
    Ok(())
}

fn write_samples(
    output: &mut Vec<u8>,
    raw: &[u8],
    header: &HeaderTypes,
) -> Result<(), &'static str> {
    if header.samples() == 0 {
        return Err("record has FORMAT data but header has no samples");
    }
    let mut fields = raw.split(|byte| *byte == b'\t');
    let format = fields.next().unwrap_or_default();
    if format.is_empty() {
        return Err("FORMAT must not be empty");
    }
    output.push(b'\t');
    output.extend_from_slice(format);

    let mut samples = 0;
    for sample in fields {
        if samples == header.samples() {
            return Err("record has more samples than the header");
        }
        output.push(b'\t');
        write_sample(output, sample, format, header)?;
        samples += 1;
    }
    if samples != header.samples() {
        return Err("record has fewer samples than the header");
    }
    Ok(())
}

fn write_sample(
    output: &mut Vec<u8>,
    raw: &[u8],
    format: &[u8],
    header: &HeaderTypes,
) -> Result<(), &'static str> {
    let mut values = raw.split(|byte| *byte == b':');
    for (index, key) in format.split(|byte| *byte == b':').enumerate() {
        if index > 0 {
            output.push(b':');
        }
        let value = values.next().unwrap_or(b".");
        match header.format(key) {
            Some(value_type) => write_typed(output, value, value_type)?,
            None => output.extend_from_slice(value),
        }
    }
    if values.next().is_some() {
        return Err("sample has more values than FORMAT declares");
    }
    Ok(())
}

fn write_typed(
    output: &mut Vec<u8>,
    raw: &[u8],
    value_type: ValueType,
) -> Result<(), &'static str> {
    match value_type {
        ValueType::Integer | ValueType::Float => {
            for (index, value) in raw.split(|byte| *byte == b',').enumerate() {
                if index > 0 {
                    output.push(b',');
                }
                if value == b"." {
                    output.push(b'.');
                } else if valid_numeric(value, value_type) {
                    reformat_into(output, value, value_type);
                } else {
                    return Err("typed numeric value is invalid or out of range");
                }
            }
        }
        ValueType::Flag | ValueType::Character | ValueType::String => {
            reformat_into(output, raw, value_type);
        }
    }
    Ok(())
}

const INT32_MIN_VALID: i64 = i32::MIN as i64 + 8;
const INT32_MAX_VALID: i64 = i32::MAX as i64;

fn valid_numeric(value: &[u8], value_type: ValueType) -> bool {
    match value_type {
        ValueType::Integer => match digits(value) {
            Some(digits) if digits.len() <= 9 => true,
            Some(_) => std::str::from_utf8(value)
                .ok()
                .and_then(|value| value.parse::<i64>().ok())
                .is_some_and(|value| (INT32_MIN_VALID..=INT32_MAX_VALID).contains(&value)),
            None => false,
        },
        _ => valid_float(value),
    }
}

fn digits(raw: &[u8]) -> Option<&[u8]> {
    let digits = match raw.split_first() {
        Some((b'-', rest)) => rest,
        _ => raw,
    };
    (!digits.is_empty() && digits.iter().all(u8::is_ascii_digit)).then_some(digits)
}

fn unsigned_digits(raw: &[u8]) -> Option<&[u8]> {
    (!raw.is_empty() && raw.iter().all(u8::is_ascii_digit)).then_some(raw)
}

fn valid_float(raw: &[u8]) -> bool {
    plain_decimal(raw)
        || std::str::from_utf8(raw)
            .ok()
            .and_then(|value| value.parse::<f64>().ok())
            .is_some()
}

fn plain_decimal(raw: &[u8]) -> bool {
    let digits = match raw.split_first() {
        Some((b'-', rest)) => rest,
        _ => raw,
    };
    let mut decimal = false;
    let mut digit = false;
    for byte in digits {
        match byte {
            b'0'..=b'9' => digit = true,
            b'.' if !decimal => decimal = true,
            _ => return false,
        }
    }
    digit
}

fn split_once(field: &[u8], separator: u8) -> Option<(&[u8], &[u8])> {
    let position = field.iter().position(|byte| *byte == separator)?;
    Some((&field[..position], &field[position + 1..]))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> HeaderTypes {
        HeaderTypes::parse(
            b"##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
        )
        .unwrap()
    }

    #[test]
    fn record_is_reformatted_without_allocation_per_field() {
        let mut output = Vec::new();
        reformat_record(
            b"chr1\t005\t.\tA\tC\t50.50\tPASS\tDP=007\tGT:DP\t0/1:020",
            &header(),
            &mut output,
        )
        .unwrap();
        assert_eq!(output, b"chr1\t5\t.\tA\tC\t50.5\tPASS\tDP=7\tGT:DP\t0/1:20");
    }

    #[test]
    fn malformed_position_and_sample_cardinality_fail() {
        let mut output = Vec::new();
        assert!(
            reformat_record(
                b"chr1\tbad\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
                &header(),
                &mut output
            )
            .is_err()
        );
        assert!(
            reformat_record(
                b"chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0",
                &header(),
                &mut output
            )
            .is_err()
        );
    }

    #[test]
    fn malformed_typed_values_fail() {
        let mut output = Vec::new();
        assert!(
            reformat_record(
                b"chr1\t1\t.\tA\tC\tbad\tPASS\t.\tGT\t0/1",
                &header(),
                &mut output
            )
            .is_err()
        );
        assert!(
            reformat_record(
                b"chr1\t1\t.\tA\tC\t.\tPASS\tDP=bad\tGT\t0/1",
                &header(),
                &mut output
            )
            .is_err()
        );
    }
}
