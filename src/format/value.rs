#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ValueType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl ValueType {
    pub(crate) fn parse(value: &[u8]) -> Option<Self> {
        Some(match value {
            b"Integer" => Self::Integer,
            b"Float" => Self::Float,
            b"Flag" => Self::Flag,
            b"Character" => Self::Character,
            b"String" => Self::String,
            _ => return None,
        })
    }
}

pub(crate) fn reformat_into(output: &mut Vec<u8>, raw: &[u8], value_type: ValueType) {
    match value_type {
        ValueType::String | ValueType::Character | ValueType::Flag => {
            output.extend_from_slice(raw);
        }
        ValueType::Integer => each_element(output, raw, write_integer),
        ValueType::Float => each_element(output, raw, write_float),
    }
}

pub(crate) fn write_f32(output: &mut Vec<u8>, value: f32) {
    write_g6(output, f64::from(value));
}

fn each_element(output: &mut Vec<u8>, raw: &[u8], render: fn(&mut Vec<u8>, &[u8])) {
    for (index, value) in raw.split(|byte| *byte == b',').enumerate() {
        if index > 0 {
            output.push(b',');
        }
        if value == b"." {
            output.push(b'.');
        } else {
            render(output, value);
        }
    }
}

const INT32_MIN_VALID: i64 = i32::MIN as i64 + 8;
const INT32_MAX_VALID: i64 = i32::MAX as i64;

fn write_integer(output: &mut Vec<u8>, value: &[u8]) {
    if is_canonical_integer(value) && integer_digit_len(value) <= 9 {
        output.extend_from_slice(value);
        return;
    }
    match std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse::<i64>().ok())
    {
        Some(value) if (INT32_MIN_VALID..=INT32_MAX_VALID).contains(&value) => {
            write_i64(output, value);
        }
        _ => output.push(b'.'),
    }
}

fn integer_digit_len(value: &[u8]) -> usize {
    value
        .first()
        .map_or(value.len(), |byte| value.len() - usize::from(*byte == b'-'))
}

fn write_float(output: &mut Vec<u8>, value: &[u8]) {
    match std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse::<f64>().ok())
    {
        Some(value) => write_g6(output, value as f32 as f64),
        None => output.extend_from_slice(value),
    }
}

fn is_canonical_integer(value: &[u8]) -> bool {
    let digits = match value.first() {
        Some(b'-') => &value[1..],
        _ => value,
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return false;
    }
    if digits.len() > 1 && digits[0] == b'0' {
        return false;
    }
    !(value[0] == b'-' && digits == b"0")
}

fn write_i64(output: &mut Vec<u8>, value: i64) {
    let mut buffer = [0; 20];
    let mut number = if value < 0 {
        (value as i128).unsigned_abs()
    } else {
        value as u128
    };
    let mut start = buffer.len();
    loop {
        start -= 1;
        buffer[start] = b'0' + (number % 10) as u8;
        number /= 10;
        if number == 0 {
            break;
        }
    }
    if value < 0 {
        output.push(b'-');
    }
    output.extend_from_slice(&buffer[start..]);
}

fn write_g6(output: &mut Vec<u8>, value: f64) {
    if value == 0.0 {
        output.extend_from_slice(if value.is_sign_negative() {
            b"-0"
        } else {
            b"0"
        });
        return;
    }
    if value.is_nan() {
        output.extend_from_slice(b"nan");
        return;
    }
    if value.is_infinite() {
        output.extend_from_slice(if value < 0.0 { b"-inf" } else { b"inf" });
        return;
    }
    let absolute = if value < 0.0 {
        output.push(b'-');
        -value
    } else {
        value
    };
    if !(0.0001..=999_999.0).contains(&absolute) {
        output.extend_from_slice(format_g6(absolute).as_bytes());
        return;
    }

    let (scale, mut exponent) = if absolute < 0.001 {
        (1e9, -4)
    } else if absolute < 0.01 {
        (1e8, -3)
    } else if absolute < 0.1 {
        (1e7, -2)
    } else if absolute < 1.0 {
        (1e6, -1)
    } else if absolute < 10.0 {
        (1e5, 0)
    } else if absolute < 100.0 {
        (1e4, 1)
    } else if absolute < 1000.0 {
        (1e3, 2)
    } else if absolute < 10_000.0 {
        (1e2, 3)
    } else if absolute < 100_000.0 {
        (1e1, 4)
    } else {
        (1.0, 5)
    };
    let mut rounded = (absolute * scale).round_ties_even() as u64;
    if rounded >= 1_000_000 {
        rounded /= 10;
        exponent += 1;
    }

    let mut digits = [0; 6];
    for digit in digits.iter_mut().rev() {
        *digit = b'0' + (rounded % 10) as u8;
        rounded /= 10;
    }
    let mut buffer = [0; 16];
    let length = if exponent < 0 {
        buffer[0] = b'0';
        buffer[1] = b'.';
        let zeros = usize::try_from(-exponent - 1).unwrap();
        buffer[2..2 + zeros].fill(b'0');
        buffer[2 + zeros..8 + zeros].copy_from_slice(&digits);
        8 + zeros
    } else if exponent >= 5 {
        buffer[..6].copy_from_slice(&digits);
        6
    } else {
        let integer_digits = usize::try_from(exponent + 1).unwrap();
        buffer[..integer_digits].copy_from_slice(&digits[..integer_digits]);
        buffer[integer_digits] = b'.';
        buffer[integer_digits + 1..7].copy_from_slice(&digits[integer_digits..]);
        7
    };
    push_stripped(output, &buffer[..length]);
}

fn push_stripped(output: &mut Vec<u8>, value: &[u8]) {
    let mut end = value.len();
    if value.contains(&b'.') {
        while end > 0 && value[end - 1] == b'0' {
            end -= 1;
        }
        if end > 0 && value[end - 1] == b'.' {
            end -= 1;
        }
    }
    output.extend_from_slice(&value[..end]);
}

fn format_g6(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    if value == 0.0 {
        return "0".to_owned();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if (-4..6).contains(&exponent) {
        let decimals = usize::try_from((5 - exponent).max(0)).unwrap();
        strip_zeros(format!("{value:.decimals$}"))
    } else {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            strip_zeros(mantissa.to_owned()),
            exponent.abs()
        )
    }
}

fn strip_zeros(mut value: String) -> String {
    if value.contains('.') {
        while value.ends_with('0') {
            value.pop();
        }
        if value.ends_with('.') {
            value.pop();
        }
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render(value: &[u8], value_type: ValueType) -> String {
        let mut output = Vec::new();
        reformat_into(&mut output, value, value_type);
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn integers_use_int32_and_missing_sentinels() {
        assert_eq!(render(b"007", ValueType::Integer), "7");
        assert_eq!(render(b"-0", ValueType::Integer), "0");
        assert_eq!(render(b"2147483647", ValueType::Integer), "2147483647");
        assert_eq!(render(b"2147483648", ValueType::Integer), ".");
        assert_eq!(render(b"-2147483648", ValueType::Integer), ".");
    }

    #[test]
    fn floats_round_through_f32_and_g6() {
        assert_eq!(render(b"0.250", ValueType::Float), "0.25");
        assert_eq!(render(b"1.0e-2", ValueType::Float), "0.01");
        assert_eq!(render(b"1e308", ValueType::Float), "inf");
        assert_eq!(render(b"0.00005", ValueType::Float), "5e-05");
    }
}
