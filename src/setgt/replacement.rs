use noodles_vcf::variant::record::samples::series::value::genotype::Phasing;
use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) enum Replacement {
    Missing,
    Reference { phased: bool },
    Minor { phased: bool },
    Major { phased: bool },
    Depth,
    Phase,
    Unphase,
    Invert,
    Custom(Template),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Template {
    pub terms: Vec<Term>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Term {
    pub value: TemplateValue,
    pub phasing: Phasing,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) enum TemplateValue {
    Position(usize),
    Missing,
    Minor,
    Major,
    Depth,
}

impl Replacement {
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "." => Ok(Self::Missing),
            "0" => Ok(Self::Reference { phased: false }),
            "0p" => Ok(Self::Reference { phased: true }),
            "m" => Ok(Self::Minor { phased: false }),
            "mp" => Ok(Self::Minor { phased: true }),
            "M" => Ok(Self::Major { phased: false }),
            "Mp" => Ok(Self::Major { phased: true }),
            "X" => Ok(Self::Depth),
            "p" => Ok(Self::Phase),
            "u" => Ok(Self::Unphase),
            "i" => Ok(Self::Invert),
            value if value.starts_with("c:") => parse_template(value).map(Self::Custom),
            _ => Err(config(format!("invalid setGT replacement: {value}"))),
        }
    }
}

fn parse_template(value: &str) -> Result<Template> {
    let expression = value
        .strip_prefix("c:")
        .ok_or_else(|| config(format!("invalid setGT replacement: {value}")))?;
    let mut raw_terms = Vec::new();
    let mut separators = Vec::new();
    let mut start = 0;

    for (index, character) in expression.char_indices() {
        if matches!(character, '/' | '|') {
            if start == index {
                return Err(config(format!("invalid setGT replacement: {value}")));
            }
            raw_terms.push(&expression[start..index]);
            separators.push(character);
            start = index + character.len_utf8();
        }
    }
    if start == expression.len() {
        return Err(config(format!("invalid setGT replacement: {value}")));
    }
    raw_terms.push(&expression[start..]);

    let first_phasing = separators
        .first()
        .map_or(Phasing::Phased, |separator| phasing(*separator));
    let mut terms = Vec::with_capacity(raw_terms.len());
    for (index, raw) in raw_terms.into_iter().enumerate() {
        let term_phasing = if index == 0 {
            first_phasing
        } else {
            phasing(separators[index - 1])
        };
        terms.push(Term {
            value: parse_term(raw, value)?,
            phasing: term_phasing,
        });
    }
    Ok(Template { terms })
}

fn parse_term(term: &str, replacement: &str) -> Result<TemplateValue> {
    match term {
        "." => Ok(TemplateValue::Missing),
        "m" => Ok(TemplateValue::Minor),
        "M" => Ok(TemplateValue::Major),
        "X" => Ok(TemplateValue::Depth),
        value if value.bytes().all(|byte| byte.is_ascii_digit()) => value
            .parse::<usize>()
            .map(TemplateValue::Position)
            .map_err(|_| config(format!("invalid setGT replacement: {replacement}"))),
        _ => Err(config(format!("invalid setGT replacement: {replacement}"))),
    }
}

fn phasing(separator: char) -> Phasing {
    match separator {
        '|' => Phasing::Phased,
        '/' => Phasing::Unphased,
        _ => unreachable!(),
    }
}

fn config(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_every_fixed_replacement() {
        for (value, expected) in [
            (".", Replacement::Missing),
            ("0", Replacement::Reference { phased: false }),
            ("0p", Replacement::Reference { phased: true }),
            ("m", Replacement::Minor { phased: false }),
            ("mp", Replacement::Minor { phased: true }),
            ("M", Replacement::Major { phased: false }),
            ("Mp", Replacement::Major { phased: true }),
            ("X", Replacement::Depth),
            ("p", Replacement::Phase),
            ("u", Replacement::Unphase),
            ("i", Replacement::Invert),
        ] {
            assert_eq!(Replacement::parse(value).unwrap(), expected);
        }
    }

    #[test]
    fn custom_templates_preserve_terms_and_mixed_separators() {
        assert_eq!(
            Replacement::parse("c:0/2|.").unwrap(),
            Replacement::Custom(Template {
                terms: vec![
                    Term {
                        value: TemplateValue::Position(0),
                        phasing: Phasing::Unphased,
                    },
                    Term {
                        value: TemplateValue::Position(2),
                        phasing: Phasing::Unphased,
                    },
                    Term {
                        value: TemplateValue::Missing,
                        phasing: Phasing::Phased,
                    },
                ],
            })
        );
        assert_eq!(
            Replacement::parse("c:m|M/X").unwrap(),
            Replacement::Custom(Template {
                terms: vec![
                    Term {
                        value: TemplateValue::Minor,
                        phasing: Phasing::Phased,
                    },
                    Term {
                        value: TemplateValue::Major,
                        phasing: Phasing::Phased,
                    },
                    Term {
                        value: TemplateValue::Depth,
                        phasing: Phasing::Unphased,
                    },
                ],
            })
        );
        assert_eq!(
            Replacement::parse("c:7").unwrap(),
            Replacement::Custom(Template {
                terms: vec![Term {
                    value: TemplateValue::Position(7),
                    phasing: Phasing::Phased,
                }],
            })
        );
    }

    #[test]
    fn rejects_empty_overflowing_and_malformed_custom_terms() {
        for value in [
            "c:",
            "c:/0",
            "c:0/",
            "c:0//1",
            "c:0||1",
            "c:0/|1",
            "c:+1",
            "c:-1",
            "c:x",
            "c:184467440737095516160000",
            "c:c:0",
        ] {
            assert!(Replacement::parse(value).is_err(), "{value}");
        }
    }

    #[test]
    fn rejects_accidental_character_mask_combinations() {
        for value in [
            "", "0u", "0i", "mu", "Mi", "Xp", "Xu", "pi", "ui", "..", "cc:0",
        ] {
            assert!(Replacement::parse(value).is_err(), "{value}");
        }
    }
}
