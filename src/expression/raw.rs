use super::bind::{BoundExpression, BoundValue, Cardinality, FieldKind, FixedField, ValueType};
use super::syntax::{BinaryOperator, UnaryOperator};

#[derive(Clone, Debug)]
pub(crate) struct Predicate(Truth);

#[derive(Clone, Debug)]
enum Truth {
    Compare(BinaryOperator, Number, Number),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
}

#[derive(Clone, Debug)]
enum Number {
    Literal(f64),
    Position,
    Quality,
    Info(String),
    Negate(Box<Self>),
    Arithmetic(BinaryOperator, Box<Self>, Box<Self>),
}

#[derive(Clone, Copy)]
enum Scalar {
    Value(f64),
    Missing,
    Unsupported,
}

struct Record<'a> {
    position: &'a [u8],
    quality: &'a [u8],
    info: &'a [u8],
}

impl Predicate {
    pub(crate) fn compile(expression: &BoundExpression) -> Option<Self> {
        compile_truth(expression).map(Self)
    }

    pub(crate) fn evaluate(&self, line: &[u8]) -> Option<bool> {
        let record = Record::parse(line)?;
        self.0.evaluate(&record)
    }
}

impl Truth {
    fn evaluate(&self, record: &Record<'_>) -> Option<bool> {
        match self {
            Self::Compare(operator, left, right) => {
                let left = left.evaluate(record);
                let right = right.evaluate(record);
                compare(*operator, left, right)
            }
            Self::And(left, right) => {
                let left = left.evaluate(record)?;
                let right = right.evaluate(record)?;
                Some(left && right)
            }
            Self::Or(left, right) => {
                let left = left.evaluate(record)?;
                let right = right.evaluate(record)?;
                Some(left || right)
            }
        }
    }
}

impl Number {
    fn evaluate(&self, record: &Record<'_>) -> Scalar {
        match self {
            Self::Literal(value) => Scalar::Value(*value),
            Self::Position => parse_number(record.position),
            Self::Quality => parse_number(record.quality),
            Self::Info(name) => record.info(name),
            Self::Negate(value) => match value.evaluate(record) {
                Scalar::Value(value) => Scalar::Value(-value),
                other => other,
            },
            Self::Arithmetic(operator, left, right) => {
                let left = left.evaluate(record);
                let right = right.evaluate(record);
                let (Scalar::Value(left), Scalar::Value(right)) = (left, right) else {
                    return merge_scalar(left, right);
                };
                Scalar::Value(match operator {
                    BinaryOperator::Add => left + right,
                    BinaryOperator::Subtract => left - right,
                    BinaryOperator::Multiply => left * right,
                    BinaryOperator::Divide => left / right,
                    _ => return Scalar::Unsupported,
                })
            }
        }
    }
}

impl<'a> Record<'a> {
    fn parse(line: &'a [u8]) -> Option<Self> {
        let mut fields = line.splitn(9, |byte| *byte == b'\t');
        fields.next()?;
        let position = fields.next()?;
        fields.next()?;
        fields.next()?;
        fields.next()?;
        let quality = fields.next()?;
        fields.next()?;
        let info = fields.next()?;
        Some(Self {
            position,
            quality,
            info,
        })
    }

    fn info(&self, requested: &str) -> Scalar {
        let mut found = None;
        for field in self.info.split(|byte| *byte == b';') {
            let Some((name, value)) = split_once(field, b'=') else {
                continue;
            };
            if name == requested.as_bytes() {
                if found.is_some() || value.contains(&b',') {
                    return Scalar::Unsupported;
                }
                found = Some(parse_number(value));
            }
        }
        found.unwrap_or(Scalar::Missing)
    }
}

fn compile_truth(expression: &BoundExpression) -> Option<Truth> {
    let BoundExpression::Binary {
        operator,
        left,
        right,
    } = expression
    else {
        return None;
    };
    match operator {
        BinaryOperator::Less
        | BinaryOperator::LessEqual
        | BinaryOperator::Greater
        | BinaryOperator::GreaterEqual => Some(Truth::Compare(
            *operator,
            compile_number(left)?,
            compile_number(right)?,
        )),
        BinaryOperator::SiteAnd => Some(Truth::And(
            Box::new(compile_truth(left)?),
            Box::new(compile_truth(right)?),
        )),
        BinaryOperator::SiteOr => Some(Truth::Or(
            Box::new(compile_truth(left)?),
            Box::new(compile_truth(right)?),
        )),
        _ => None,
    }
}

fn compile_number(expression: &BoundExpression) -> Option<Number> {
    match expression {
        BoundExpression::Value(BoundValue::Number(value)) => Some(Number::Literal(*value)),
        BoundExpression::Value(BoundValue::Field(field))
            if field.subscript.is_none() && field.cardinality == Cardinality::Count(1) =>
        {
            match (&field.kind, field.value_type) {
                (FieldKind::Fixed(FixedField::Position), ValueType::Integer) => {
                    Some(Number::Position)
                }
                (FieldKind::Fixed(FixedField::Quality), ValueType::Float) => Some(Number::Quality),
                (FieldKind::Info(name), ValueType::Integer | ValueType::Float) => {
                    Some(Number::Info(name.clone()))
                }
                _ => None,
            }
        }
        BoundExpression::Unary {
            operator: UnaryOperator::Negate,
            expression,
        } => Some(Number::Negate(Box::new(compile_number(expression)?))),
        BoundExpression::Binary {
            operator,
            left,
            right,
        } if matches!(
            operator,
            BinaryOperator::Add
                | BinaryOperator::Subtract
                | BinaryOperator::Multiply
                | BinaryOperator::Divide
        ) =>
        {
            Some(Number::Arithmetic(
                *operator,
                Box::new(compile_number(left)?),
                Box::new(compile_number(right)?),
            ))
        }
        _ => None,
    }
}

fn compare(operator: BinaryOperator, left: Scalar, right: Scalar) -> Option<bool> {
    match (left, right) {
        (Scalar::Unsupported, _) | (_, Scalar::Unsupported) => None,
        (Scalar::Missing, _) | (_, Scalar::Missing) => Some(false),
        (Scalar::Value(left), Scalar::Value(right)) => Some(match operator {
            BinaryOperator::Less => left < right,
            BinaryOperator::LessEqual => left <= right,
            BinaryOperator::Greater => left > right,
            BinaryOperator::GreaterEqual => left >= right,
            _ => return None,
        }),
    }
}

fn merge_scalar(left: Scalar, right: Scalar) -> Scalar {
    if matches!(left, Scalar::Unsupported) || matches!(right, Scalar::Unsupported) {
        Scalar::Unsupported
    } else {
        Scalar::Missing
    }
}

fn parse_number(value: &[u8]) -> Scalar {
    if value == b"." {
        return Scalar::Missing;
    }
    std::str::from_utf8(value)
        .ok()
        .and_then(|value| value.parse().ok())
        .map_or(Scalar::Unsupported, Scalar::Value)
}

fn split_once(value: &[u8], separator: u8) -> Option<(&[u8], &[u8])> {
    let index = value.iter().position(|byte| *byte == separator)?;
    Some((&value[..index], &value[index + 1..]))
}
