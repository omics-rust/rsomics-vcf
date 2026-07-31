use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::Path;

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::format::{HeaderTypes, Reader, ValueType, reformat_record, trim_line_ending};
use crate::query_bcf;
use crate::query_format::{self, Field, Token};
use crate::variant_type;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum HeaderMode {
    #[default]
    None,
    Indexed,
    Plain,
}

#[derive(Debug)]
pub struct Options {
    pub format: String,
    pub samples: Option<Vec<String>>,
    pub header: HeaderMode,
    pub automatic_newline: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct Summary {
    pub records: usize,
}

pub fn write(input: &Path, options: &Options, mut output: impl Write) -> Result<Summary> {
    let mut reader = Reader::open(input)?;
    let (header, _, schema) = reader.read_header()?;
    let tokens = query_format::parse(&options.format)?;
    validate(&tokens, &schema, false)?;
    let selected = select_samples(&schema, options.samples.as_deref())?;
    let automatic_newline = options.automatic_newline && !query_format::contains_newline(&tokens);
    let context = Context {
        tokens: &tokens,
        schema: &schema,
        selected: &selected,
        automatic_newline,
    };

    let mut output = BufWriter::with_capacity(1 << 20, &mut output);
    write_header(&context, options.header, &mut output)?;

    let mut scratch = Scratch::default();
    let mut bcf_record = bcf::Record::default();
    let text = reader.is_text();
    let direct_bcf = !text && query_bcf::supports(&tokens, &schema);
    let mut records = 0;
    loop {
        let number = records + 1;
        let read = if text {
            reader.read_text_record(&mut scratch.input, number)?
        } else {
            let read = reader.read_bcf_record(&mut bcf_record, number)?;
            if read > 0 && direct_bcf {
                scratch.rendered.clear();
                query_bcf::render(context.tokens, &header, &bcf_record, &mut scratch.rendered)
                    .map_err(|error| {
                        RsomicsError::InvalidInput(format!(
                            "{}: projecting variant record {number}: {error}",
                            input.display()
                        ))
                    })?;
                if context.automatic_newline {
                    scratch.rendered.push(b'\n');
                }
                output
                    .write_all(&scratch.rendered)
                    .map_err(RsomicsError::Io)?;
                records += 1;
                continue;
            }
            scratch.input.clear();
            if read > 0 {
                vcf::io::Writer::new(&mut scratch.input)
                    .write_variant_record(&header, &bcf_record)
                    .map_err(|error| {
                        RsomicsError::InvalidInput(format!(
                            "{}: decoding variant record {number}: {error}",
                            input.display()
                        ))
                    })?;
                trim_line_ending(&mut scratch.input);
            }
            read
        };
        if read == 0 {
            break;
        }

        reformat_record(&scratch.input, &schema, &mut scratch.canonical).map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "{}: parsing variant record {number}: {error}",
                input.display()
            ))
        })?;
        render_record(&context, &mut scratch, &mut output)?;
        records += 1;
    }

    output.flush().map_err(RsomicsError::Io)?;
    Ok(Summary { records })
}

#[derive(Default)]
struct Scratch {
    input: Vec<u8>,
    canonical: Vec<u8>,
    columns: Vec<Range<usize>>,
    format_keys: Vec<Range<usize>>,
    sample_values: Vec<Range<usize>>,
    rendered: Vec<u8>,
}

struct Context<'a> {
    tokens: &'a [Token],
    schema: &'a HeaderTypes,
    selected: &'a [usize],
    automatic_newline: bool,
}

fn render_record(
    context: &Context<'_>,
    scratch: &mut Scratch,
    output: &mut impl Write,
) -> Result<()> {
    scratch.columns.clear();
    split_spans(
        &scratch.canonical,
        0..scratch.canonical.len(),
        b'\t',
        &mut scratch.columns,
    );
    scratch.format_keys.clear();
    if let Some(format) = scratch.columns.get(8) {
        split_spans(
            &scratch.canonical,
            format.clone(),
            b':',
            &mut scratch.format_keys,
        );
    }

    scratch.rendered.clear();
    render_tokens(
        context,
        context.tokens,
        &scratch.canonical,
        &scratch.columns,
        &scratch.format_keys,
        None,
        &mut scratch.sample_values,
        &mut scratch.rendered,
    )?;
    if context.automatic_newline {
        scratch.rendered.push(b'\n');
    }
    output
        .write_all(&scratch.rendered)
        .map_err(RsomicsError::Io)
}

#[allow(clippy::too_many_arguments)]
fn render_tokens(
    context: &Context<'_>,
    tokens: &[Token],
    line: &[u8],
    columns: &[Range<usize>],
    format_keys: &[Range<usize>],
    sample: Option<usize>,
    sample_values: &mut Vec<Range<usize>>,
    output: &mut Vec<u8>,
) -> Result<()> {
    for token in tokens {
        match token {
            Token::Literal(value) => output.extend_from_slice(value),
            Token::Field { field, subscript } => write_field(
                context,
                field,
                *subscript,
                line,
                columns,
                format_keys,
                sample,
                sample_values,
                output,
            )?,
            Token::Samples(inner) => {
                for &index in context.selected {
                    sample_values.clear();
                    if let Some(value) = columns.get(9 + index) {
                        split_spans(line, value.clone(), b':', sample_values);
                    }
                    render_tokens(
                        context,
                        inner,
                        line,
                        columns,
                        format_keys,
                        Some(index),
                        sample_values,
                        output,
                    )?;
                }
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_field(
    context: &Context<'_>,
    field: &Field,
    subscript: Option<usize>,
    line: &[u8],
    columns: &[Range<usize>],
    format_keys: &[Range<usize>],
    sample: Option<usize>,
    sample_values: &[Range<usize>],
    output: &mut Vec<u8>,
) -> Result<()> {
    let column = |index| resolve(line, columns, index);
    match field {
        Field::Chrom => output.extend_from_slice(column(0)),
        Field::Pos => output.extend_from_slice(column(1)),
        Field::Pos0 => write_integer(output, parse_integer(column(1))? - 1),
        Field::End => write_integer(output, record_end(column(1), column(3), column(7))?),
        Field::End0 => write_integer(output, record_end(column(1), column(3), column(7))? - 1),
        Field::Id => output.extend_from_slice(column(2)),
        Field::Ref => output.extend_from_slice(column(3)),
        Field::Alt => write_subscript(output, column(4), subscript),
        Field::FirstAlt => write_subscript(output, column(4), Some(0)),
        Field::Qual => output.extend_from_slice(column(5)),
        Field::Filter => output.extend_from_slice(column(6)),
        Field::Type => variant_type::write(output, column(3), column(4)),
        Field::Info => output.extend_from_slice(column(7)),
        Field::InfoTag(name) => write_info(context.schema, column(7), name, subscript, output),
        Field::Format => write_format(context.selected, line, columns, output),
        Field::Line => write_line(context.selected, line, columns, output),
        Field::Gt => {
            output.extend_from_slice(sample_value(b"GT", line, format_keys, sample_values))
        }
        Field::Tgt => write_translated_genotype(
            sample_value(b"GT", line, format_keys, sample_values),
            column(3),
            column(4),
            output,
        ),
        Field::IupacGt => write_iupac_genotype(
            sample_value(b"GT", line, format_keys, sample_values),
            column(3),
            column(4),
            output,
        ),
        Field::Sample => {
            let index = sample.ok_or_else(|| invalid("%SAMPLE must be inside [ ]"))?;
            output.extend_from_slice(
                context
                    .schema
                    .sample_name(index)
                    .ok_or_else(|| invalid("sample index is outside the VCF header"))?,
            );
        }
        Field::Named(name) => {
            if sample.is_some() && context.schema.format(name).is_some() {
                write_subscript(
                    output,
                    sample_value(name, line, format_keys, sample_values),
                    subscript,
                );
            } else if sample.is_some() {
                write_fixed(name, subscript, context.selected, line, columns, output)?;
            } else {
                write_info(context.schema, column(7), name, subscript, output);
            }
        }
        Field::Fixed(name) => {
            write_fixed(name, subscript, context.selected, line, columns, output)?
        }
        Field::Vkey => {
            return Err(invalid("%VKX is not enabled in the current query engine"));
        }
    }
    Ok(())
}

fn write_info(
    schema: &HeaderTypes,
    info: &[u8],
    name: &[u8],
    subscript: Option<usize>,
    output: &mut Vec<u8>,
) {
    for field in info.split(|value| *value == b';') {
        match split_once(field, b'=') {
            Some((key, value)) if key == name => {
                write_info_subscript(output, value, subscript);
                return;
            }
            None if field == name && schema.info(name) == Some(ValueType::Flag) => {
                output.push(b'1');
                return;
            }
            _ => {}
        }
    }
    output.push(b'.');
}

fn write_format(selected: &[usize], line: &[u8], columns: &[Range<usize>], output: &mut Vec<u8>) {
    let Some(format) = columns.get(8) else {
        output.push(b'.');
        return;
    };
    if selected.is_empty() {
        output.extend_from_slice(b"\t.");
        return;
    }
    output.extend_from_slice(&line[format.clone()]);
    for &sample in selected {
        output.push(b'\t');
        output.extend_from_slice(resolve(line, columns, 9 + sample));
    }
}

fn write_line(selected: &[usize], line: &[u8], columns: &[Range<usize>], output: &mut Vec<u8>) {
    let end = if columns.len() > 8 && !selected.is_empty() {
        9
    } else {
        columns.len().min(8)
    };
    for index in 0..end {
        if index > 0 {
            output.push(b'\t');
        }
        output.extend_from_slice(resolve(line, columns, index));
    }
    if columns.len() > 8 {
        for &sample in selected {
            output.push(b'\t');
            output.extend_from_slice(resolve(line, columns, 9 + sample));
        }
    }
}

fn write_fixed(
    name: &[u8],
    subscript: Option<usize>,
    selected: &[usize],
    line: &[u8],
    columns: &[Range<usize>],
    output: &mut Vec<u8>,
) -> Result<()> {
    match name {
        b"CHROM" => output.extend_from_slice(resolve(line, columns, 0)),
        b"POS" => output.extend_from_slice(resolve(line, columns, 1)),
        b"POS0" => write_integer(output, parse_integer(resolve(line, columns, 1))? - 1),
        b"END" => write_integer(
            output,
            record_end(
                resolve(line, columns, 1),
                resolve(line, columns, 3),
                resolve(line, columns, 7),
            )?,
        ),
        b"END0" => write_integer(
            output,
            record_end(
                resolve(line, columns, 1),
                resolve(line, columns, 3),
                resolve(line, columns, 7),
            )? - 1,
        ),
        b"ID" => output.extend_from_slice(resolve(line, columns, 2)),
        b"REF" => output.extend_from_slice(resolve(line, columns, 3)),
        b"ALT" => {
            write_subscript(output, resolve(line, columns, 4), subscript);
        }
        b"FIRST_ALT" => write_subscript(output, resolve(line, columns, 4), Some(0)),
        b"QUAL" => output.extend_from_slice(resolve(line, columns, 5)),
        b"FILTER" => output.extend_from_slice(resolve(line, columns, 6)),
        b"TYPE" => {
            variant_type::write(output, resolve(line, columns, 3), resolve(line, columns, 4))
        }
        b"LINE" => write_line(selected, line, columns, output),
        _ => return Err(invalid("unknown fixed VCF field")),
    }
    Ok(())
}

fn write_header(context: &Context<'_>, mode: HeaderMode, output: &mut impl Write) -> Result<()> {
    if mode == HeaderMode::None {
        return Ok(());
    }
    let mut row = vec![b'#'];
    let mut column = 1;
    header_tokens(context, context.tokens, None, mode, &mut column, &mut row);
    if context.automatic_newline {
        row.push(b'\n');
    }
    output.write_all(&row).map_err(RsomicsError::Io)
}

fn header_tokens(
    context: &Context<'_>,
    tokens: &[Token],
    sample: Option<usize>,
    mode: HeaderMode,
    column: &mut usize,
    output: &mut Vec<u8>,
) {
    for token in tokens {
        match token {
            Token::Literal(value) => output.extend_from_slice(value),
            Token::Field { field, .. } => {
                if mode == HeaderMode::Indexed {
                    output.extend_from_slice(format!("[{column}]").as_bytes());
                }
                if let Some(index) = sample {
                    output.extend_from_slice(context.schema.sample_name(index).unwrap_or(b"."));
                    output.push(b':');
                }
                output.extend_from_slice(field_label(field));
                *column += 1;
            }
            Token::Samples(inner) => {
                for &index in context.selected {
                    header_tokens(context, inner, Some(index), mode, column, output);
                }
            }
        }
    }
}

fn field_label(field: &Field) -> &[u8] {
    match field {
        Field::Chrom => b"CHROM",
        Field::Pos => b"POS",
        Field::Pos0 => b"POS0",
        Field::End => b"END",
        Field::End0 => b"END0",
        Field::Id => b"ID",
        Field::Ref => b"REF",
        Field::Alt => b"ALT",
        Field::FirstAlt => b"FIRST_ALT",
        Field::Qual => b"QUAL",
        Field::Filter => b"FILTER",
        Field::Type => b"TYPE",
        Field::Info => b"INFO",
        Field::Format => b"FORMAT",
        Field::Line => b"LINE",
        Field::Gt | Field::Tgt | Field::IupacGt => b"GT",
        Field::Sample => b"SAMPLE",
        Field::Vkey => b"VKX",
        Field::InfoTag(name) | Field::Named(name) | Field::Fixed(name) => name,
    }
}

fn validate(tokens: &[Token], schema: &HeaderTypes, in_samples: bool) -> Result<()> {
    for token in tokens {
        match token {
            Token::Literal(_) => {}
            Token::Samples(inner) => validate(inner, schema, true)?,
            Token::Field { field, .. } => match field {
                Field::InfoTag(name) if schema.info(name).is_none() => {
                    return Err(no_tag("INFO", name));
                }
                Field::Named(name) if in_samples => {
                    if schema.format(name).is_none() && !is_fixed(name) {
                        return Err(no_tag("FORMAT", name));
                    }
                }
                Field::Named(name) if schema.info(name).is_none() => {
                    return Err(no_tag("INFO", name));
                }
                Field::Fixed(name) if !is_fixed(name) => {
                    return Err(invalid("unknown fixed VCF field"));
                }
                Field::Vkey => {
                    return Err(invalid("%VKX is not enabled in the current query engine"));
                }
                _ => {}
            },
        }
    }
    Ok(())
}

fn select_samples(schema: &HeaderTypes, wanted: Option<&[String]>) -> Result<Vec<usize>> {
    let Some(wanted) = wanted else {
        return Ok((0..schema.samples()).collect());
    };
    if let Some(first) = wanted.first().and_then(|name| name.strip_prefix('^')) {
        let mut excluded = vec![first.as_bytes()];
        excluded.extend(wanted[1..].iter().map(String::as_bytes));
        for name in &excluded {
            find_sample(schema, name)?;
        }
        return Ok((0..schema.samples())
            .filter(|index| {
                schema
                    .sample_name(*index)
                    .is_some_and(|name| !excluded.contains(&name))
            })
            .collect());
    }
    wanted
        .iter()
        .map(|name| find_sample(schema, name.as_bytes()))
        .collect()
}

fn find_sample(schema: &HeaderTypes, name: &[u8]) -> Result<usize> {
    (0..schema.samples())
        .find(|index| schema.sample_name(*index) == Some(name))
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "sample is not present in the VCF header: {}",
                String::from_utf8_lossy(name)
            ))
        })
}

fn sample_value<'a>(
    key: &[u8],
    line: &'a [u8],
    format_keys: &[Range<usize>],
    sample_values: &[Range<usize>],
) -> &'a [u8] {
    let Some(index) = format_keys
        .iter()
        .position(|range| &line[range.clone()] == key)
    else {
        return b".";
    };
    sample_values
        .get(index)
        .map_or(b".", |range| &line[range.clone()])
}

fn write_translated_genotype(
    genotype: &[u8],
    reference: &[u8],
    alternates: &[u8],
    output: &mut Vec<u8>,
) {
    let mut start = 0;
    for (index, value) in genotype.iter().enumerate() {
        if matches!(value, b'/' | b'|') {
            write_allele(output, &genotype[start..index], reference, alternates);
            output.push(*value);
            start = index + 1;
        }
    }
    write_allele(output, &genotype[start..], reference, alternates);
}

fn write_iupac_genotype(
    genotype: &[u8],
    reference: &[u8],
    alternates: &[u8],
    output: &mut Vec<u8>,
) {
    let alleles: Vec<_> = genotype
        .split(|value| matches!(value, b'/' | b'|'))
        .map(|value| allele(value, reference, alternates))
        .collect();
    if alleles.len() == 2
        && alleles
            .iter()
            .all(|value| value.is_some_and(|value| value.len() == 1))
    {
        let first = alleles[0].unwrap()[0].to_ascii_uppercase();
        let second = alleles[1].unwrap()[0].to_ascii_uppercase();
        if let Some(code) = iupac(first, second) {
            output.push(code);
            return;
        }
    }
    write_translated_genotype(genotype, reference, alternates, output);
}

fn iupac(first: u8, second: u8) -> Option<u8> {
    Some(match (first.min(second), first.max(second)) {
        (b'A', b'A') => b'A',
        (b'A', b'C') => b'M',
        (b'A', b'G') => b'R',
        (b'A', b'T') => b'W',
        (b'C', b'C') => b'C',
        (b'C', b'G') => b'S',
        (b'C', b'T') => b'Y',
        (b'G', b'G') => b'G',
        (b'G', b'T') => b'K',
        (b'T', b'T') => b'T',
        _ => return None,
    })
}

fn write_allele(output: &mut Vec<u8>, index: &[u8], reference: &[u8], alternates: &[u8]) {
    match allele(index, reference, alternates) {
        Some(value) => output.extend_from_slice(value),
        None => output.push(b'.'),
    }
}

fn allele<'a>(index: &[u8], reference: &'a [u8], alternates: &'a [u8]) -> Option<&'a [u8]> {
    let index = std::str::from_utf8(index).ok()?.parse::<usize>().ok()?;
    if index == 0 {
        Some(reference)
    } else {
        alternates.split(|value| *value == b',').nth(index - 1)
    }
}

fn record_end(position: &[u8], reference: &[u8], info: &[u8]) -> Result<i64> {
    if let Some(end) = info
        .split(|value| *value == b';')
        .filter_map(|field| split_once(field, b'='))
        .find_map(|(key, value)| (key == b"END").then_some(value))
    {
        return parse_integer(end);
    }
    parse_integer(position)?
        .checked_add(reference.len() as i64 - 1)
        .ok_or_else(|| invalid("POS plus REF length overflows i64"))
}

fn parse_integer(value: &[u8]) -> Result<i64> {
    std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "expected an integer, found {}",
                String::from_utf8_lossy(value)
            ))
        })
}

pub(crate) fn write_integer(output: &mut Vec<u8>, value: i64) {
    let mut buffer = [0; 20];
    let mut number = if value < 0 {
        (value as i128).unsigned_abs()
    } else {
        value as u128
    };
    let mut index = buffer.len();
    loop {
        index -= 1;
        buffer[index] = b'0' + (number % 10) as u8;
        number /= 10;
        if number == 0 {
            break;
        }
    }
    if value < 0 {
        output.push(b'-');
    }
    output.extend_from_slice(&buffer[index..]);
}

fn write_subscript(output: &mut Vec<u8>, value: &[u8], subscript: Option<usize>) {
    output.extend_from_slice(match subscript {
        Some(index) => value
            .split(|value| *value == b',')
            .nth(index)
            .unwrap_or(b"."),
        None => value,
    });
}

fn write_info_subscript(output: &mut Vec<u8>, value: &[u8], subscript: Option<usize>) {
    let Some(index) = subscript else {
        output.extend_from_slice(value);
        return;
    };
    let mut values = value.split(|value| *value == b',');
    let first = values.next().unwrap_or(b".");
    let selected = if index == 0 {
        first
    } else {
        match values.next() {
            None => first,
            Some(second) if index == 1 => second,
            Some(_) => values.nth(index - 2).unwrap_or(b"."),
        }
    };
    output.extend_from_slice(selected);
}

fn split_spans(line: &[u8], within: Range<usize>, separator: u8, output: &mut Vec<Range<usize>>) {
    let mut start = within.start;
    for index in within.clone() {
        if line[index] == separator {
            output.push(start..index);
            start = index + 1;
        }
    }
    output.push(start..within.end);
}

fn resolve<'a>(line: &'a [u8], columns: &[Range<usize>], index: usize) -> &'a [u8] {
    columns
        .get(index)
        .map_or(b".", |range| &line[range.clone()])
}

fn split_once(field: &[u8], separator: u8) -> Option<(&[u8], &[u8])> {
    let position = field.iter().position(|value| *value == separator)?;
    Some((&field[..position], &field[position + 1..]))
}

fn is_fixed(name: &[u8]) -> bool {
    matches!(
        name,
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
            | b"TYPE"
            | b"LINE"
    )
}

fn no_tag(kind: &str, name: &[u8]) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "no such {kind} tag defined in the VCF header: {}",
        String::from_utf8_lossy(name)
    ))
}

pub(crate) fn invalid(message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(message.to_owned())
}
