use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Target {
    pub principal: Principal,
    pub random_fraction: Option<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) enum Principal {
    AnyMissing,
    PartialMissing,
    CompleteMissing,
    All,
    Query,
    Binomial(Binomial),
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Binomial {
    pub tag: String,
    pub comparison: Comparison,
    pub threshold: f64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Comparison {
    Less,
    LessEqual,
    Equal,
    GreaterEqual,
    Greater,
}

impl Target {
    pub fn parse(values: &[String]) -> Result<Self> {
        match values {
            [value] if value.starts_with("r:") => Ok(Self {
                principal: Principal::All,
                random_fraction: Some(parse_fraction(value)?),
            }),
            [principal] => Ok(Self {
                principal: parse_principal(principal)?,
                random_fraction: None,
            }),
            [principal, random] if random.starts_with("r:") => Ok(Self {
                principal: parse_principal(principal)?,
                random_fraction: Some(parse_fraction(random)?),
            }),
            _ => Err(config("invalid setGT target composition")),
        }
    }
}

fn parse_principal(value: &str) -> Result<Principal> {
    match value {
        "." => Ok(Principal::AnyMissing),
        "./x" => Ok(Principal::PartialMissing),
        "./." => Ok(Principal::CompleteMissing),
        "a" => Ok(Principal::All),
        "q" => Ok(Principal::Query),
        value if value.starts_with("b:") => parse_binomial(value).map(Principal::Binomial),
        _ => Err(config(format!("invalid setGT target: {value}"))),
    }
}

fn parse_fraction(value: &str) -> Result<f64> {
    let fraction = value
        .strip_prefix("r:")
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value > 0.0 && *value < 1.0)
        .ok_or_else(|| config(format!("invalid setGT random fraction: {value}")))?;
    Ok(fraction)
}

fn parse_binomial(value: &str) -> Result<Binomial> {
    let expression = value
        .strip_prefix("b:")
        .ok_or_else(|| config(format!("invalid setGT binomial target: {value}")))?;
    let operator_start = expression
        .find(['<', '=', '>'])
        .ok_or_else(|| config(format!("invalid setGT binomial target: {value}")))?;
    let operator_end = expression[operator_start..]
        .find(|character: char| !matches!(character, '<' | '=' | '>'))
        .map_or(expression.len(), |offset| operator_start + offset);
    let tag = &expression[..operator_start];
    let operator = &expression[operator_start..operator_end];
    let threshold = &expression[operator_end..];

    if !valid_tag(tag) || threshold.contains(['<', '=', '>']) {
        return Err(config(format!("invalid setGT binomial target: {value}")));
    }

    let comparison = match operator {
        "<" => Comparison::Less,
        "<=" => Comparison::LessEqual,
        "=" | "==" => Comparison::Equal,
        ">=" => Comparison::GreaterEqual,
        ">" => Comparison::Greater,
        _ => return Err(config(format!("invalid setGT binomial target: {value}"))),
    };
    let threshold = threshold
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| config(format!("invalid setGT binomial target: {value}")))?;

    Ok(Binomial {
        tag: tag.to_owned(),
        comparison,
        threshold,
    })
}

fn valid_tag(value: &str) -> bool {
    let mut characters = value.chars();
    characters
        .next()
        .is_some_and(|character| character.is_ascii_alphabetic() || character == '_')
        && characters
            .all(|character| character.is_ascii_alphanumeric() || matches!(character, '_' | '.'))
}

fn config(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(values: &[&str]) -> rsomics_common::Result<Target> {
        Target::parse(
            &values
                .iter()
                .map(|value| (*value).to_owned())
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn parses_every_principal_and_random_composition() {
        for (value, expected) in [
            (".", Principal::AnyMissing),
            ("./x", Principal::PartialMissing),
            ("./.", Principal::CompleteMissing),
            ("a", Principal::All),
            ("q", Principal::Query),
        ] {
            assert_eq!(parse(&[value]).unwrap().principal, expected);
        }

        let target = parse(&["./x", "r:0.25"]).unwrap();
        assert_eq!(target.principal, Principal::PartialMissing);
        assert_eq!(target.random_fraction, Some(0.25));

        let target = parse(&["r:0.75"]).unwrap();
        assert_eq!(target.principal, Principal::All);
        assert_eq!(target.random_fraction, Some(0.75));
    }

    #[test]
    fn parses_every_binomial_comparison() {
        for (operator, expected) in [
            ("<", Comparison::Less),
            ("<=", Comparison::LessEqual),
            ("=", Comparison::Equal),
            ("==", Comparison::Equal),
            (">=", Comparison::GreaterEqual),
            (">", Comparison::Greater),
        ] {
            let target = parse(&[&format!("b:AD{operator}0.125")]).unwrap();
            let Principal::Binomial(binomial) = target.principal else {
                panic!("expected binomial target")
            };
            assert_eq!(binomial.tag, "AD");
            assert_eq!(binomial.comparison, expected);
            assert_eq!(binomial.threshold, 0.125);
        }
    }

    #[test]
    fn rejects_empty_duplicate_conflicting_and_undocumented_targets() {
        for values in [
            vec![],
            vec![""],
            vec!["a", "."],
            vec!["r:0.5", "a"],
            vec!["a", "r:0.5", "r:0.25"],
            vec!["r:0.5", "r:0.25"],
            vec!["ar:0.5"],
            vec!["A"],
        ] {
            assert!(parse(&values).is_err(), "{values:?}");
        }
    }

    #[test]
    fn random_fraction_must_be_finite_and_inside_the_unit_interval() {
        for value in [
            "r:", "r:0", "r:-0", "r:-0.1", "r:1", "r:1.1", "r:NaN", "r:inf", "r:-inf",
        ] {
            assert!(parse(&[value]).is_err(), "{value}");
        }
    }

    #[test]
    fn rejects_malformed_binomial_targets() {
        for value in [
            "b:",
            "b:AD",
            "b:<0.1",
            "b:1AD<0.1",
            "b:AD<",
            "b:AD!=0.1",
            "b:AD<>0.1",
            "b:AD<0.1>0",
            "b:AD===0.1",
            "b:AD<NaN",
            "b:AD<inf",
        ] {
            assert!(parse(&[value]).is_err(), "{value}");
        }
    }
}
