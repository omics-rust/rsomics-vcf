use rsomics_common::{Result, RsomicsError};

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
    use super::*;

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
}
