use noodles_vcf::{
    Header,
    header::record::value::map::format,
    variant::{
        RecordBuf,
        record::samples::series::value::genotype::Phasing,
        record_buf::samples::sample::{
            Value as SampleValue,
            value::{Array as SampleArray, Genotype, genotype::Allele},
        },
    },
};
use rsomics_common::{Result, RsomicsError};

use crate::genotype::allele_counts;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Replacement {
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
pub(crate) struct Template {
    pub(crate) terms: Vec<Term>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Term {
    pub(crate) value: TemplateValue,
    pub(crate) phasing: Phasing,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum TemplateValue {
    Position(usize),
    Missing,
    Minor,
    Major,
    Depth,
}

impl Replacement {
    pub(crate) fn parse(value: &str) -> Result<Self> {
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

    pub(crate) fn validate(&self, header: &Header) -> Result<()> {
        if !self.uses_depth() {
            return Ok(());
        }
        let schema = header
            .formats()
            .get("AD")
            .ok_or_else(|| config("setGT replacement X requires FORMAT/AD in the header"))?;
        if schema.number() != format::Number::ReferenceAlternateBases
            || schema.ty() != format::Type::Integer
        {
            return Err(config(
                "setGT replacement X requires FORMAT/AD Number=R,Type=Integer",
            ));
        }
        Ok(())
    }

    pub(crate) fn prepare(&self, record: &RecordBuf, selected: &[bool]) -> Result<Prepared> {
        let sample_count = record.samples().values().count();
        if selected.len() != sample_count {
            return Err(invalid(format!(
                "genotype selection has {} samples but the record has {sample_count}",
                selected.len()
            )));
        }
        let any_selected = selected.iter().any(|selected| *selected);
        let allele_count = record.alternate_bases().as_ref().len() + 1;
        let (major, minor) = if any_selected && (self.uses_major() || self.uses_minor()) {
            ranked_alleles(record, self.uses_minor())?
        } else {
            (None, None)
        };
        let depths = if any_selected && self.uses_depth() {
            depth_alleles(record, selected, allele_count)?
        } else {
            vec![None; sample_count]
        };
        Ok(Prepared {
            replacement: self.clone(),
            allele_count,
            major,
            minor,
            depths,
        })
    }

    fn uses_major(&self) -> bool {
        matches!(self, Self::Major { .. })
            || matches!(self, Self::Custom(template) if template.terms.iter().any(|term| term.value == TemplateValue::Major))
    }

    fn uses_minor(&self) -> bool {
        matches!(self, Self::Minor { .. })
            || matches!(self, Self::Custom(template) if template.terms.iter().any(|term| term.value == TemplateValue::Minor))
    }

    fn uses_depth(&self) -> bool {
        matches!(self, Self::Depth)
            || matches!(self, Self::Custom(template) if template.terms.iter().any(|term| term.value == TemplateValue::Depth))
    }
}

pub(crate) struct Prepared {
    replacement: Replacement,
    allele_count: usize,
    major: Option<usize>,
    minor: Option<usize>,
    depths: Vec<Option<usize>>,
}

impl Prepared {
    pub(crate) fn resolve(&self, sample: usize, genotype: &Genotype) -> Result<Genotype> {
        match &self.replacement {
            Replacement::Missing => Ok(filled(genotype.as_ref().len(), None, Phasing::Unphased)),
            Replacement::Reference { phased } => {
                Ok(filled(genotype.as_ref().len(), Some(0), phase(*phased)))
            }
            Replacement::Minor { phased } => Ok(filled(
                genotype.as_ref().len(),
                Some(
                    self.minor
                        .ok_or_else(|| invalid("minor allele is unavailable"))?,
                ),
                phase(*phased),
            )),
            Replacement::Major { phased } => Ok(filled(
                genotype.as_ref().len(),
                Some(
                    self.major
                        .ok_or_else(|| invalid("major allele is unavailable"))?,
                ),
                phase(*phased),
            )),
            Replacement::Depth => Ok(filled(
                genotype.as_ref().len(),
                self.depth(sample)?,
                Phasing::Unphased,
            )),
            Replacement::Phase => Ok(genotype
                .as_ref()
                .iter()
                .map(|allele| Allele::new(allele.position(), Phasing::Phased))
                .collect()),
            Replacement::Unphase => {
                let mut positions = genotype
                    .as_ref()
                    .iter()
                    .map(|allele| allele.position())
                    .collect::<Vec<_>>();
                positions.sort_by_key(|position| (position.is_some(), position.unwrap_or(0)));
                Ok(positions
                    .into_iter()
                    .map(|position| Allele::new(position, Phasing::Unphased))
                    .collect())
            }
            Replacement::Invert if genotype.as_ref().len() == 2 => {
                let separator = genotype.as_ref()[1].phasing();
                Ok(genotype
                    .as_ref()
                    .iter()
                    .rev()
                    .map(|allele| Allele::new(allele.position(), separator))
                    .collect())
            }
            Replacement::Invert => Ok(genotype.clone()),
            Replacement::Custom(template) => template
                .terms
                .iter()
                .map(|term| {
                    self.term(sample, &term.value)
                        .map(|position| Allele::new(position, term.phasing))
                })
                .collect(),
        }
    }

    fn term(&self, sample: usize, value: &TemplateValue) -> Result<Option<usize>> {
        match value {
            TemplateValue::Position(position) => {
                Ok((*position < self.allele_count).then_some(*position))
            }
            TemplateValue::Missing => Ok(None),
            TemplateValue::Minor => self
                .minor
                .map(Some)
                .ok_or_else(|| invalid("minor allele is unavailable")),
            TemplateValue::Major => self
                .major
                .map(Some)
                .ok_or_else(|| invalid("major allele is unavailable")),
            TemplateValue::Depth => self.depth(sample),
        }
    }

    fn depth(&self, sample: usize) -> Result<Option<usize>> {
        self.depths
            .get(sample)
            .copied()
            .ok_or_else(|| invalid(format!("sample {} has no depth resolution", sample + 1)))
    }
}

fn ranked_alleles(record: &RecordBuf, needs_minor: bool) -> Result<(Option<usize>, Option<usize>)> {
    let counts = allele_counts(record)?;
    if counts.total == 0 {
        return Err(invalid(
            "major or minor allele cannot be resolved without called alleles",
        ));
    }
    let mut ranked = (0..counts.counts.len()).collect::<Vec<_>>();
    ranked.sort_by(|left, right| {
        counts.counts[*right]
            .cmp(&counts.counts[*left])
            .then_with(|| left.cmp(right))
    });
    if needs_minor && ranked.len() < 2 {
        return Err(invalid(
            "minor allele requires at least two declared alleles",
        ));
    }
    Ok((ranked.first().copied(), ranked.get(1).copied()))
}

fn depth_alleles(
    record: &RecordBuf,
    selected: &[bool],
    allele_count: usize,
) -> Result<Vec<Option<usize>>> {
    let values = record
        .samples()
        .select("AD")
        .ok_or_else(|| invalid("record has no FORMAT/AD field"))?;
    selected
        .iter()
        .enumerate()
        .map(|(sample, selected)| {
            if !selected {
                return Ok(None);
            }
            let Some(value) = values.get(sample) else {
                return Err(invalid(format!(
                    "sample {} has no FORMAT/AD value",
                    sample + 1
                )));
            };
            let Some(value) = value else {
                return Ok(None);
            };
            let SampleValue::Array(SampleArray::Integer(depths)) = value else {
                return Err(invalid(format!(
                    "sample {} FORMAT/AD is not encoded as an integer array",
                    sample + 1
                )));
            };
            if depths.len() != allele_count {
                return Err(invalid(format!(
                    "sample {} FORMAT/AD has {} values, expected {allele_count}",
                    sample + 1,
                    depths.len()
                )));
            }
            if depths.iter().flatten().any(|depth| *depth < 0) {
                return Err(invalid(format!(
                    "sample {} FORMAT/AD contains a negative depth",
                    sample + 1
                )));
            }
            let mut maximum = None;
            for (allele, depth) in depths.iter().enumerate() {
                let Some(depth) = depth else {
                    continue;
                };
                if maximum.is_none_or(|(_, current)| depth > current) {
                    maximum = Some((allele, depth));
                }
            }
            Ok(maximum.map(|(allele, _)| allele))
        })
        .collect()
}

fn filled(ploidy: usize, position: Option<usize>, phasing: Phasing) -> Genotype {
    (0..ploidy)
        .map(|_| Allele::new(position, phasing))
        .collect()
}

fn phase(phased: bool) -> Phasing {
    if phased {
        Phasing::Phased
    } else {
        Phasing::Unphased
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

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
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
