use std::fmt;

use noodles_vcf::{
    self as vcf,
    header::record::value::map::{format, info},
};

use super::syntax::{
    BinaryOperator, Expression, Field, Namespace, Selector, Subscript, UnaryOperator, Value,
};

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum BoundExpression {
    Value(BoundValue),
    Unary {
        operator: UnaryOperator,
        expression: Box<Self>,
    },
    Binary {
        operator: BinaryOperator,
        left: Box<Self>,
        right: Box<Self>,
    },
    Function {
        name: String,
        arguments: Vec<Self>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum BoundValue {
    Number(f64),
    String(String),
    Missing,
    File(std::path::PathBuf),
    Field(BoundField),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct BoundField {
    pub kind: FieldKind,
    pub value_type: ValueType,
    pub cardinality: Cardinality,
    pub subscript: Option<Subscript>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum FieldKind {
    Fixed(FixedField),
    Info(String),
    Format(String),
    Calculated(CalculatedField),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum FixedField {
    Chrom,
    Position,
    Id,
    Reference,
    Alternate,
    Quality,
    Filter,
    Type,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum CalculatedField {
    AlternateCount,
    SampleCount,
    AlleleCount,
    MinorAlleleCount,
    AlleleFrequency,
    MinorAlleleFrequency,
    TotalAlleleCount,
    MissingSampleCount,
    MissingSampleFraction,
    IndelLength,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ValueType {
    Integer,
    Float,
    Flag,
    Character,
    String,
    Genotype,
    VariantType,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Cardinality {
    Count(usize),
    AlternateBases,
    ReferenceAlternateBases,
    Genotypes,
    LocalAlternateBases,
    LocalReferenceAlternateBases,
    LocalGenotypes,
    Ploidy,
    BaseModifications,
    Unknown,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct BindError(String);

impl BindError {
    fn new(message: impl Into<String>) -> Self {
        Self(message.into())
    }
}

impl fmt::Display for BindError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for BindError {}

pub(crate) fn bind(
    expression: Expression,
    header: &vcf::Header,
) -> Result<BoundExpression, BindError> {
    match expression {
        Expression::Value(value) => bind_value(value, header).map(BoundExpression::Value),
        Expression::Unary {
            operator,
            expression,
        } => Ok(BoundExpression::Unary {
            operator,
            expression: Box::new(bind(*expression, header)?),
        }),
        Expression::Binary {
            operator,
            left,
            right,
        } => Ok(BoundExpression::Binary {
            operator,
            left: Box::new(bind(*left, header)?),
            right: Box::new(bind(*right, header)?),
        }),
        Expression::Function { name, arguments } => Ok(BoundExpression::Function {
            name: name.to_ascii_uppercase(),
            arguments: arguments
                .into_iter()
                .map(|argument| bind(argument, header))
                .collect::<Result<_, _>>()?,
        }),
    }
}

fn bind_value(value: Value, header: &vcf::Header) -> Result<BoundValue, BindError> {
    match value {
        Value::Number(value) => Ok(BoundValue::Number(value)),
        Value::String(value) => Ok(BoundValue::String(value)),
        Value::Missing => Ok(BoundValue::Missing),
        Value::File(path) => Ok(BoundValue::File(path)),
        Value::Field(field) => bind_field(field, header).map(BoundValue::Field),
    }
}

fn bind_field(field: Field, header: &vcf::Header) -> Result<BoundField, BindError> {
    let name = field.name.as_str();
    let upper = name.to_ascii_uppercase();
    let mut bound = match field.namespace {
        Some(Namespace::Info) => bind_info(name, header)?,
        Some(Namespace::Format) => bind_format(name, header)?,
        None => {
            if let Some(fixed) = fixed_field(&upper) {
                fixed
            } else if let Some(calculated) = calculated_field(&upper, header) {
                calculated
            } else {
                let info = header.infos().get(name);
                let format = header.formats().get(name);
                match (info, format) {
                    (Some(_), Some(_)) => {
                        return Err(BindError::new(format!(
                            "ambiguous field {name}: both INFO/{name} and FORMAT/{name} are declared"
                        )));
                    }
                    (Some(_), None) => bind_info(name, header)?,
                    (None, Some(_)) => bind_format(name, header)?,
                    (None, None) => {
                        return Err(BindError::new(format!("undefined field {name}")));
                    }
                }
            }
        }
    };
    validate_subscript(&mut bound, field.subscript)?;
    Ok(bound)
}

fn bind_info(name: &str, header: &vcf::Header) -> Result<BoundField, BindError> {
    let definition = header
        .infos()
        .get(name)
        .ok_or_else(|| BindError::new(format!("undefined INFO field {name}")))?;
    Ok(BoundField {
        kind: FieldKind::Info(name.to_owned()),
        value_type: match definition.ty() {
            info::Type::Integer => ValueType::Integer,
            info::Type::Float => ValueType::Float,
            info::Type::Flag => ValueType::Flag,
            info::Type::Character => ValueType::Character,
            info::Type::String => ValueType::String,
        },
        cardinality: match definition.number() {
            info::Number::Count(value) => Cardinality::Count(value),
            info::Number::AlternateBases => Cardinality::AlternateBases,
            info::Number::ReferenceAlternateBases => Cardinality::ReferenceAlternateBases,
            info::Number::Samples => Cardinality::Genotypes,
            info::Number::Unknown => Cardinality::Unknown,
        },
        subscript: None,
    })
}

fn bind_format(name: &str, header: &vcf::Header) -> Result<BoundField, BindError> {
    let definition = header
        .formats()
        .get(name)
        .ok_or_else(|| BindError::new(format!("undefined FORMAT field {name}")))?;
    Ok(BoundField {
        kind: FieldKind::Format(name.to_owned()),
        value_type: if name == "GT" {
            ValueType::Genotype
        } else {
            match definition.ty() {
                format::Type::Integer => ValueType::Integer,
                format::Type::Float => ValueType::Float,
                format::Type::Character => ValueType::Character,
                format::Type::String => ValueType::String,
            }
        },
        cardinality: match definition.number() {
            format::Number::Count(value) => Cardinality::Count(value),
            format::Number::AlternateBases => Cardinality::AlternateBases,
            format::Number::ReferenceAlternateBases => Cardinality::ReferenceAlternateBases,
            format::Number::Samples => Cardinality::Genotypes,
            format::Number::LocalAlternateBases => Cardinality::LocalAlternateBases,
            format::Number::LocalReferenceAlternateBases => {
                Cardinality::LocalReferenceAlternateBases
            }
            format::Number::LocalSamples => Cardinality::LocalGenotypes,
            format::Number::Ploidy => Cardinality::Ploidy,
            format::Number::BaseModifications => Cardinality::BaseModifications,
            format::Number::Unknown => Cardinality::Unknown,
        },
        subscript: None,
    })
}

fn fixed_field(name: &str) -> Option<BoundField> {
    let (field, value_type, cardinality) = match name {
        "CHROM" => (FixedField::Chrom, ValueType::String, Cardinality::Count(1)),
        "POS" => (
            FixedField::Position,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "ID" => (FixedField::Id, ValueType::String, Cardinality::Unknown),
        "REF" => (
            FixedField::Reference,
            ValueType::String,
            Cardinality::Count(1),
        ),
        "ALT" => (
            FixedField::Alternate,
            ValueType::String,
            Cardinality::AlternateBases,
        ),
        "QUAL" => (FixedField::Quality, ValueType::Float, Cardinality::Count(1)),
        "FILTER" => (FixedField::Filter, ValueType::String, Cardinality::Unknown),
        "TYPE" => (
            FixedField::Type,
            ValueType::VariantType,
            Cardinality::AlternateBases,
        ),
        _ => return None,
    };
    Some(BoundField {
        kind: FieldKind::Fixed(field),
        value_type,
        cardinality,
        subscript: None,
    })
}

fn calculated_field(name: &str, header: &vcf::Header) -> Option<BoundField> {
    let (field, value_type, cardinality) = match name {
        "N_ALT" => (
            CalculatedField::AlternateCount,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "N_SAMPLES" => (
            CalculatedField::SampleCount,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "AC" if !header.infos().contains_key("AC") => (
            CalculatedField::AlleleCount,
            ValueType::Integer,
            Cardinality::AlternateBases,
        ),
        "MAC" if !header.infos().contains_key("MAC") => (
            CalculatedField::MinorAlleleCount,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "AF" if !header.infos().contains_key("AF") => (
            CalculatedField::AlleleFrequency,
            ValueType::Float,
            Cardinality::AlternateBases,
        ),
        "MAF" if !header.infos().contains_key("MAF") => (
            CalculatedField::MinorAlleleFrequency,
            ValueType::Float,
            Cardinality::Count(1),
        ),
        "AN" if !header.infos().contains_key("AN") => (
            CalculatedField::TotalAlleleCount,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "N_MISSING" => (
            CalculatedField::MissingSampleCount,
            ValueType::Integer,
            Cardinality::Count(1),
        ),
        "F_MISSING" => (
            CalculatedField::MissingSampleFraction,
            ValueType::Float,
            Cardinality::Count(1),
        ),
        "ILEN" => (
            CalculatedField::IndelLength,
            ValueType::Integer,
            Cardinality::AlternateBases,
        ),
        _ => return None,
    };
    Some(BoundField {
        kind: FieldKind::Calculated(field),
        value_type,
        cardinality,
        subscript: None,
    })
}

fn validate_subscript(
    field: &mut BoundField,
    subscript: Option<Subscript>,
) -> Result<(), BindError> {
    let Some(subscript) = subscript else {
        return Ok(());
    };
    if matches!(subscript, Subscript::SampleValues { .. })
        && !matches!(field.kind, FieldKind::Format(_))
    {
        return Err(BindError::new(
            "sample:value subscripts require a FORMAT field",
        ));
    }
    for selector in selectors(&subscript) {
        if matches!(selector, Selector::Genotype)
            && !matches!(
                field.cardinality,
                Cardinality::ReferenceAlternateBases | Cardinality::AlternateBases
            )
        {
            return Err(BindError::new(
                "GT-selected subscripts require Number=R or Number=A",
            ));
        }
        if matches!(selector, Selector::File(_)) && !matches!(field.kind, FieldKind::Format(_)) {
            return Err(BindError::new(
                "file-selected samples require a FORMAT field",
            ));
        }
    }
    field.subscript = Some(subscript);
    Ok(())
}

fn selectors(subscript: &Subscript) -> impl Iterator<Item = &Selector> {
    let (first, second) = match subscript {
        Subscript::Values(values) => (values, None),
        Subscript::SampleValues { samples, values } => (samples, Some(values)),
    };
    std::iter::once(first).chain(second)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expression::syntax::parse;

    fn header() -> vcf::Header {
        "##fileformat=VCFv4.3\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"site depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"frequency\">\n\
##INFO=<ID=R,Number=R,Type=Integer,Description=\"alleles\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"sample depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"allele depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            .parse()
            .unwrap()
    }

    fn bound_field(source: &str) -> BoundField {
        let BoundExpression::Value(BoundValue::Field(field)) =
            bind(parse(source).unwrap(), &header()).unwrap()
        else {
            panic!("expected field");
        };
        field
    }

    #[test]
    fn header_binding_retains_namespace_type_number_and_subscript() {
        assert_eq!(
            bound_field("INFO/AF[0]"),
            BoundField {
                kind: FieldKind::Info("AF".to_owned()),
                value_type: ValueType::Float,
                cardinality: Cardinality::AlternateBases,
                subscript: Some(Subscript::Values(Selector::Indices(vec![
                    super::super::syntax::IndexSelector {
                        start: 0,
                        end: None,
                    }
                ]))),
            }
        );
        let ad = bound_field("FMT/AD[0:GT]");
        assert_eq!(ad.kind, FieldKind::Format("AD".to_owned()));
        assert_eq!(ad.value_type, ValueType::Integer);
        assert_eq!(ad.cardinality, Cardinality::ReferenceAlternateBases);
    }

    #[test]
    fn fixed_and_calculated_fields_have_stable_types() {
        assert_eq!(
            bound_field("QUAL").kind,
            FieldKind::Fixed(FixedField::Quality)
        );
        assert_eq!(bound_field("TYPE").value_type, ValueType::VariantType);
        assert_eq!(
            bound_field("N_MISSING").kind,
            FieldKind::Calculated(CalculatedField::MissingSampleCount)
        );
        assert_eq!(bound_field("AF").kind, FieldKind::Info("AF".to_owned()));
        assert_eq!(
            bound_field("AC").kind,
            FieldKind::Calculated(CalculatedField::AlleleCount)
        );
    }

    #[test]
    fn ambiguous_undefined_and_invalid_subscripts_fail_during_binding() {
        for source in [
            "DP > 1",
            "INFO/MISSING > 1",
            "FORMAT/MISSING > 1",
            "QUAL[0:1] > 1",
            "INFO/R[0:GT] > 1",
            "INFO/R[@samples] > 1",
            "FMT/DP[GT] > 1",
        ] {
            assert!(bind(parse(source).unwrap(), &header()).is_err(), "{source}");
        }
    }

    #[test]
    fn function_names_are_bound_case_insensitively() {
        let BoundExpression::Function { name, arguments } =
            bind(parse("smpl_mean(FMT/AD)").unwrap(), &header()).unwrap()
        else {
            panic!("expected function");
        };
        assert_eq!(name, "SMPL_MEAN");
        assert_eq!(arguments.len(), 1);
    }
}
