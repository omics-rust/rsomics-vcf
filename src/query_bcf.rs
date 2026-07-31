use noodles_bcf as bcf;
use noodles_vcf::{
    self as vcf,
    variant::record::{
        AlternateBases as _, Filters as _,
        info::field::{Value as InfoValue, value::Array as InfoArray},
    },
};
use rsomics_common::Result;

use crate::format::{HeaderTypes, ValueType, write_f32};
use crate::query::{invalid, write_integer};
use crate::query_format::{Field, Token};

pub(crate) fn supports(tokens: &[Token], schema: &HeaderTypes) -> bool {
    tokens.iter().all(|token| match token {
        Token::Literal(_) => true,
        Token::Samples(_) => false,
        Token::Field { field, .. } => match field {
            Field::Chrom
            | Field::Pos
            | Field::Pos0
            | Field::End
            | Field::End0
            | Field::Id
            | Field::Ref
            | Field::Alt
            | Field::FirstAlt
            | Field::Qual
            | Field::Filter => true,
            Field::InfoTag(name) | Field::Named(name) => matches!(
                schema.info(name),
                Some(ValueType::Integer | ValueType::Float | ValueType::Flag)
            ),
            Field::Fixed(name) => matches!(
                name.as_slice(),
                b"CHROM"
                    | b"POS"
                    | b"POS0"
                    | b"END"
                    | b"END0"
                    | b"ID"
                    | b"REF"
                    | b"ALT"
                    | b"FIRST_ALT"
                    | b"QUAL"
                    | b"FILTER"
            ),
            Field::Type
            | Field::Info
            | Field::Format
            | Field::Line
            | Field::Gt
            | Field::Tgt
            | Field::IupacGt
            | Field::Sample
            | Field::Vkey => false,
        },
    })
}

pub(crate) fn render(
    tokens: &[Token],
    header: &vcf::Header,
    record: &bcf::Record,
    output: &mut Vec<u8>,
) -> Result<()> {
    for token in tokens {
        match token {
            Token::Literal(value) => output.extend_from_slice(value),
            Token::Samples(_) => {
                return Err(invalid("sample loop requires the typed BCF fallback"));
            }
            Token::Field { field, subscript } => {
                write_field(header, record, field, *subscript, output)?;
            }
        }
    }
    Ok(())
}

fn write_field(
    header: &vcf::Header,
    record: &bcf::Record,
    field: &Field,
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) -> Result<()> {
    match field {
        Field::Chrom => write_chrom(header, record, output)?,
        Field::Pos => write_integer(output, position(record)?),
        Field::Pos0 => write_integer(output, position(record)? - 1),
        Field::End => write_integer(output, record.end()?.get() as i64),
        Field::End0 => write_integer(output, record.end()?.get() as i64 - 1),
        Field::Id => write_id(record, output),
        Field::Ref => output.extend_from_slice(record.reference_bases().as_ref()),
        Field::Alt => write_alternates(record, subscript, output)?,
        Field::FirstAlt => write_alternates(record, Some(0), output)?,
        Field::Qual => write_quality(record, output)?,
        Field::Filter => write_filters(header, record, output)?,
        Field::InfoTag(name) | Field::Named(name) => {
            write_info(header, record, name, subscript, output)?;
        }
        Field::Fixed(name) => write_fixed(header, record, name, subscript, output)?,
        Field::Type
        | Field::Info
        | Field::Format
        | Field::Line
        | Field::Gt
        | Field::Tgt
        | Field::IupacGt
        | Field::Sample
        | Field::Vkey => return Err(invalid("query field requires the typed BCF fallback")),
    }
    Ok(())
}

fn write_fixed(
    header: &vcf::Header,
    record: &bcf::Record,
    name: &[u8],
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) -> Result<()> {
    match name {
        b"CHROM" => write_chrom(header, record, output)?,
        b"POS" => write_integer(output, position(record)?),
        b"POS0" => write_integer(output, position(record)? - 1),
        b"END" => write_integer(output, record.end()?.get() as i64),
        b"END0" => write_integer(output, record.end()?.get() as i64 - 1),
        b"ID" => write_id(record, output),
        b"REF" => output.extend_from_slice(record.reference_bases().as_ref()),
        b"ALT" => write_alternates(record, subscript, output)?,
        b"FIRST_ALT" => write_alternates(record, Some(0), output)?,
        b"QUAL" => write_quality(record, output)?,
        b"FILTER" => write_filters(header, record, output)?,
        _ => return Err(invalid("fixed query field requires the typed BCF fallback")),
    }
    Ok(())
}

fn write_chrom(header: &vcf::Header, record: &bcf::Record, output: &mut Vec<u8>) -> Result<()> {
    output.extend_from_slice(
        record
            .reference_sequence_name(header.string_maps())?
            .as_bytes(),
    );
    Ok(())
}

fn write_id(record: &bcf::Record, output: &mut Vec<u8>) {
    let ids = record.ids();
    if ids.as_ref().is_empty() {
        output.push(b'.');
    } else {
        output.extend_from_slice(ids.as_ref());
    }
}

fn write_quality(record: &bcf::Record, output: &mut Vec<u8>) -> Result<()> {
    match record.quality_score()? {
        Some(value) => write_f32(output, value),
        None => output.push(b'.'),
    }
    Ok(())
}

fn position(record: &bcf::Record) -> Result<i64> {
    record
        .variant_start()
        .transpose()?
        .map(|position| position.get() as i64)
        .ok_or_else(|| invalid("BCF record is missing POS"))
}

fn write_alternates(
    record: &bcf::Record,
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) -> Result<()> {
    let alternates = record.alternate_bases();
    let mut values = alternates.iter();
    if let Some(index) = subscript {
        match values.nth(index).transpose()? {
            Some(value) => output.extend_from_slice(value.as_bytes()),
            None => output.push(b'.'),
        }
        return Ok(());
    }

    let mut count = 0;
    for result in values {
        if count > 0 {
            output.push(b',');
        }
        output.extend_from_slice(result?.as_bytes());
        count += 1;
    }
    if count == 0 {
        output.push(b'.');
    }
    Ok(())
}

fn write_filters(header: &vcf::Header, record: &bcf::Record, output: &mut Vec<u8>) -> Result<()> {
    let filters = record.filters();
    let mut count = 0;
    for result in filters.iter(header) {
        if count > 0 {
            output.push(b';');
        }
        output.extend_from_slice(result?.as_bytes());
        count += 1;
    }
    if count == 0 {
        output.push(b'.');
    }
    Ok(())
}

fn write_info(
    header: &vcf::Header,
    record: &bcf::Record,
    name: &[u8],
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) -> Result<()> {
    let name =
        std::str::from_utf8(name).map_err(|_| invalid("query INFO tag is not valid UTF-8"))?;
    let info = record.info();
    let Some(value) = info.get(header, name) else {
        output.push(b'.');
        return Ok(());
    };
    match value? {
        Some(value) => write_info_value(value, subscript, output),
        None => {
            output.push(b'1');
            Ok(())
        }
    }
}

fn write_info_value(
    value: InfoValue<'_>,
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) -> Result<()> {
    match value {
        InfoValue::Integer(value) => {
            write_integer(output, i64::from(value));
            Ok(())
        }
        InfoValue::Float(value) => {
            write_f32(output, value);
            Ok(())
        }
        InfoValue::Flag => {
            output.push(b'1');
            Ok(())
        }
        InfoValue::Array(InfoArray::Integer(values)) => {
            write_array(values.iter(), subscript, output, |output, value| {
                write_integer(output, i64::from(value))
            })
        }
        InfoValue::Array(InfoArray::Float(values)) => {
            write_array(values.iter(), subscript, output, write_f32)
        }
        InfoValue::Character(_)
        | InfoValue::String(_)
        | InfoValue::Array(InfoArray::Character(_) | InfoArray::String(_)) => {
            Err(invalid("string INFO value requires the typed BCF fallback"))
        }
    }
}

fn write_array<T>(
    mut values: impl Iterator<Item = std::io::Result<Option<T>>>,
    subscript: Option<usize>,
    output: &mut Vec<u8>,
    write: impl Fn(&mut Vec<u8>, T),
) -> Result<()>
where
    T: Copy,
{
    if let Some(index) = subscript {
        let first = values.next().transpose()?;
        let selected = if index == 0 {
            first
        } else {
            match values.next().transpose()? {
                None => first,
                second if index == 1 => second,
                Some(_) => values.nth(index - 2).transpose()?,
            }
        };
        match selected {
            Some(Some(value)) => write(output, value),
            Some(None) | None => output.push(b'.'),
        }
        return Ok(());
    }

    let mut count = 0;
    for result in values {
        if count > 0 {
            output.push(b',');
        }
        match result? {
            Some(value) => write(output, value),
            None => output.push(b'.'),
        }
        count += 1;
    }
    if count == 0 {
        output.push(b'.');
    }
    Ok(())
}
