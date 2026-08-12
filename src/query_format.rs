use std::{cell::RefCell, ops::Range};

use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_common::{Result, RsomicsError};

use crate::format::{HeaderTypes, trim_line_ending};

#[derive(Debug)]
pub(crate) struct SiteFormat {
    tokens: Vec<Token>,
    schema: HeaderTypes,
    scratch: RefCell<SiteScratch>,
}

#[derive(Debug, Default)]
struct SiteScratch {
    line: Vec<u8>,
    columns: Vec<Range<usize>>,
    sample_values: Vec<Range<usize>>,
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Field {
    Alt,
    Chrom,
    End,
    End0,
    Filter,
    FirstAlt,
    Format,
    Gt,
    Id,
    Info,
    InfoTag(Vec<u8>),
    IupacGt,
    Line,
    Named(Vec<u8>),
    Pos,
    Pos0,
    Qual,
    Ref,
    Sample,
    Tgt,
    Type,
    Vkey,
    Fixed(Vec<u8>),
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Token {
    Literal(Vec<u8>),
    Field {
        field: Field,
        subscript: Option<usize>,
    },
    Samples(Vec<Token>),
}

pub(crate) fn parse(format: &str) -> Result<Vec<Token>> {
    let mut position = 0;
    let tokens = parse_sequence(format.as_bytes(), &mut position, false)?;
    if position != format.len() {
        return Err(invalid("unmatched ] in query format"));
    }
    Ok(tokens)
}

pub(crate) fn contains_newline(tokens: &[Token]) -> bool {
    tokens.iter().any(|token| match token {
        Token::Literal(value) => value.contains(&b'\n'),
        Token::Samples(inner) => contains_newline(inner),
        Token::Field { .. } => false,
    })
}

impl SiteFormat {
    pub(crate) fn bind(source: &str, schema: &HeaderTypes) -> Result<Self> {
        let tokens = parse(source)?;
        if tokens.is_empty() {
            return Err(invalid("site format must not be empty"));
        }
        validate(&tokens, schema, false)?;
        validate_site(&tokens)?;
        Ok(Self {
            tokens,
            schema: schema.clone(),
            scratch: RefCell::new(SiteScratch::default()),
        })
    }

    pub(crate) fn render(
        &self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
        output: &mut Vec<u8>,
    ) -> Result<()> {
        let mut scratch = self.scratch.borrow_mut();
        scratch.line.clear();
        vcf::io::Writer::new(&mut scratch.line)
            .write_variant_record(header, record)
            .map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "serializing variant record for site format: {error}"
                ))
            })?;
        trim_line_ending(&mut scratch.line);
        output.clear();
        let SiteScratch {
            line,
            columns,
            sample_values,
        } = &mut *scratch;
        crate::query::render_site(
            &self.tokens,
            &self.schema,
            line,
            columns,
            sample_values,
            output,
        )
    }
}

pub(crate) fn validate(tokens: &[Token], schema: &HeaderTypes, in_samples: bool) -> Result<()> {
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

fn validate_site(tokens: &[Token]) -> Result<()> {
    for token in tokens {
        match token {
            Token::Literal(value) if value.iter().any(u8::is_ascii_whitespace) => {
                return Err(invalid("site format literals must not contain whitespace"));
            }
            Token::Literal(_) => {}
            Token::Samples(_) => return Err(invalid("site formats do not support sample loops")),
            Token::Field { field, .. } => match field {
                Field::Format
                | Field::Gt
                | Field::IupacGt
                | Field::Line
                | Field::Sample
                | Field::Tgt => {
                    return Err(invalid(
                        "site formats do not support sample or record fields",
                    ));
                }
                Field::Fixed(name) if name == b"FORMAT" || name == b"LINE" => {
                    return Err(invalid(
                        "site formats do not support sample or record fields",
                    ));
                }
                _ => {}
            },
        }
    }
    Ok(())
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

fn parse_sequence(input: &[u8], position: &mut usize, in_samples: bool) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let mut literal = Vec::new();

    while *position < input.len() {
        match input[*position] {
            b']' if in_samples => {
                *position += 1;
                flush(&mut literal, &mut tokens);
                return Ok(tokens);
            }
            b']' => return Err(invalid("unmatched ] in query format")),
            b'[' if in_samples => return Err(invalid("nested sample loops are not supported")),
            b'[' => {
                flush(&mut literal, &mut tokens);
                *position += 1;
                tokens.push(Token::Samples(parse_sequence(input, position, true)?));
            }
            b'\\' => {
                *position += 1;
                match input.get(*position) {
                    Some(b'n') => literal.push(b'\n'),
                    Some(b't') => literal.push(b'\t'),
                    Some(b'r') => literal.push(b'\r'),
                    Some(b'\\') => literal.push(b'\\'),
                    Some(value) => literal.push(*value),
                    None => literal.push(b'\\'),
                }
                *position += 1;
            }
            b'%' => {
                flush(&mut literal, &mut tokens);
                *position += 1;
                tokens.push(parse_field(input, position, in_samples)?);
            }
            value => {
                literal.push(value);
                *position += 1;
            }
        }
    }

    if in_samples {
        return Err(invalid("unterminated sample loop in query format"));
    }
    flush(&mut literal, &mut tokens);
    Ok(tokens)
}

fn parse_field(input: &[u8], position: &mut usize, in_samples: bool) -> Result<Token> {
    let start = *position;
    while input
        .get(*position)
        .is_some_and(|value| value.is_ascii_alphanumeric() || matches!(value, b'_' | b'.' | b'/'))
    {
        *position += 1;
    }
    let name = &input[start..*position];
    if name.is_empty() {
        return Err(invalid("a query field must follow %"));
    }
    if input.get(*position) == Some(&b'(') {
        return Err(invalid(
            "query functions require the expression engine and are not enabled",
        ));
    }

    let subscript = parse_subscript(input, position)?;
    let field = resolve(name, in_samples)?;
    Ok(Token::Field { field, subscript })
}

fn parse_subscript(input: &[u8], position: &mut usize) -> Result<Option<usize>> {
    if input.get(*position) != Some(&b'{') {
        return Ok(None);
    }
    let start = *position + 1;
    let Some(length) = input[start..].iter().position(|value| *value == b'}') else {
        return Err(invalid("unterminated query subscript"));
    };
    let end = start + length;
    let value = std::str::from_utf8(&input[start..end])
        .ok()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| invalid("query subscript must be a non-negative integer"))?;
    *position = end + 1;
    Ok(Some(value))
}

fn resolve(name: &[u8], in_samples: bool) -> Result<Field> {
    if in_samples {
        return Ok(match name {
            b"GT" => Field::Gt,
            b"TGT" => Field::Tgt,
            b"IUPACGT" => Field::IupacGt,
            b"SAMPLE" => Field::Sample,
            b"INFO" => {
                return Err(invalid(
                    "%INFO inside [ ] must name an INFO tag, such as %INFO/DP",
                ));
            }
            _ if name.starts_with(b"INFO/") && name.len() > 5 => Field::InfoTag(name[5..].to_vec()),
            _ if name.starts_with(b"/") && name.len() > 1 => Field::Fixed(name[1..].to_vec()),
            _ if name.contains(&b'/') => {
                return Err(invalid("unsupported query field path"));
            }
            _ => Field::Named(name.to_vec()),
        });
    }

    let field = match name {
        b"CHROM" => Field::Chrom,
        b"POS" => Field::Pos,
        b"POS0" => Field::Pos0,
        b"END" => Field::End,
        b"END0" => Field::End0,
        b"ID" => Field::Id,
        b"REF" => Field::Ref,
        b"ALT" => Field::Alt,
        b"FIRST_ALT" => Field::FirstAlt,
        b"QUAL" => Field::Qual,
        b"FILTER" => Field::Filter,
        b"TYPE" => Field::Type,
        b"INFO" => Field::Info,
        b"FORMAT" => Field::Format,
        b"LINE" => Field::Line,
        b"VKX" => Field::Vkey,
        _ if name.starts_with(b"INFO/") && name.len() > 5 => Field::InfoTag(name[5..].to_vec()),
        _ if name.starts_with(b"/") && name.len() > 1 => Field::Fixed(name[1..].to_vec()),
        _ if name.contains(&b'/') => {
            return Err(invalid("unsupported query field path"));
        }
        _ => Field::Named(name.to_vec()),
    };
    Ok(field)
}

fn flush(literal: &mut Vec<u8>, tokens: &mut Vec<Token>) {
    if !literal.is_empty() {
        tokens.push(Token::Literal(std::mem::take(literal)));
    }
}

fn invalid(message: &str) -> RsomicsError {
    RsomicsError::InvalidInput(message.to_owned())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{self as vcf, variant::RecordBuf};

    use crate::format::HeaderTypes;

    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn site() -> (vcf::Header, HeaderTypes, RecordBuf) {
        let header: vcf::Header = HEADER.parse().unwrap();
        let schema = HeaderTypes::parse(HEADER.as_bytes()).unwrap();
        let mut reader =
            vcf::io::Reader::new(b"chr1\t10\trs1\tA\tC,G\t12\tPASS\tDP=7\n".as_slice());
        let raw = reader.records().next().unwrap().unwrap();
        let record = RecordBuf::try_from_variant_record(&header, &raw).unwrap();
        (header, schema, record)
    }

    #[test]
    fn parses_site_and_sample_fields() {
        let tokens = parse(r"%CHROM\t%INFO/DP[\t%SAMPLE=%GT]\n").unwrap();
        assert!(contains_newline(&tokens));
        assert_eq!(
            tokens,
            vec![
                Token::Field {
                    field: Field::Chrom,
                    subscript: None,
                },
                Token::Literal(b"\t".to_vec()),
                Token::Field {
                    field: Field::InfoTag(b"DP".to_vec()),
                    subscript: None,
                },
                Token::Samples(vec![
                    Token::Literal(b"\t".to_vec()),
                    Token::Field {
                        field: Field::Sample,
                        subscript: None,
                    },
                    Token::Literal(b"=".to_vec()),
                    Token::Field {
                        field: Field::Gt,
                        subscript: None,
                    },
                ]),
                Token::Literal(b"\n".to_vec()),
            ]
        );
    }

    #[test]
    fn newline_anywhere_disables_automatic_newline() {
        let tokens = parse(r"%CHROM\n%POS").unwrap();
        assert!(contains_newline(&tokens));
    }

    #[test]
    fn rejects_incomplete_and_expression_syntax() {
        assert!(parse("%").is_err());
        assert!(parse("[%GT").is_err());
        assert!(parse("%SUM(FMT/AD)").is_err());
    }

    #[test]
    fn site_format_renders_fixed_info_and_escaped_literals() {
        let (header, schema, record) = site();
        let format = SiteFormat::bind(
            r"%CHROM\_%POS\_%ID\_%REF\_%ALT{1}\_%FIRST_ALT\_%QUAL\_%FILTER\_%INFO/DP\_%TYPE",
            &schema,
        )
        .unwrap();
        let mut output = Vec::new();

        format.render(&header, &record, &mut output).unwrap();

        assert_eq!(output, b"chr1_10_rs1_A_G_C_12_PASS_7_SNP");
    }

    #[test]
    fn site_format_rejects_sample_fields_and_invalid_literals() {
        let (_, schema, _) = site();

        for value in ["[%GT]", "%FORMAT", "%GT", "x\\ty", "x\\ny", "", "x y"] {
            assert!(SiteFormat::bind(value, &schema).is_err(), "{value:?}");
        }
    }
}
