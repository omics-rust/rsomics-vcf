use noodles_vcf::variant::RecordBuf;

use crate::expression::{
    bind::{BoundExpression, BoundSubscript, BoundValue, Cardinality, FieldKind, ValueSelector},
    value::{self, Atom, SampleValues, Values},
};

use super::{
    super::{EvaluateError, sample_width},
    numeric_count,
};

pub(super) fn evaluate<'a>(
    mut arguments: Vec<Values<'a>>,
    expressions: &[BoundExpression],
    record: &RecordBuf,
) -> Result<Values<'a>, EvaluateError> {
    match arguments.len() {
        1 => one(arguments.pop().expect("single FISHER argument")),
        2 => {
            let right = arguments.pop().expect("second FISHER argument");
            let left = arguments.pop().expect("first FISHER argument");
            two(left, right, expressions, record)
        }
        _ => Err(EvaluateError::new("FISHER requires one or two arguments")),
    }
}

fn one<'a>(values: Values<'a>) -> Result<Values<'a>, EvaluateError> {
    match values {
        Values::Site(values) => {
            let value = match values.as_slice() {
                [n11, n12, n21, n22] => probability(n11, n12, n21, n22)?,
                _ => Atom::Missing,
            };
            Ok(Values::Site(vec![value]))
        }
        Values::Samples(samples) => {
            let values = samples
                .values
                .iter()
                .zip(&samples.selected)
                .map(|(counts, selected)| {
                    if !selected {
                        return Ok(vec![Atom::Missing]);
                    }
                    match counts.as_slice() {
                        [n11, n12, n21, n22] => {
                            probability(n11, n12, n21, n22).map(|value| vec![value])
                        }
                        _ => Ok(vec![Atom::Missing]),
                    }
                })
                .collect::<Result<_, _>>()?;
            Ok(Values::Samples(SampleValues {
                values,
                selected: samples.selected,
            }))
        }
    }
}

fn two<'a>(
    left: Values<'a>,
    right: Values<'a>,
    expressions: &[BoundExpression],
    record: &RecordBuf,
) -> Result<Values<'a>, EvaluateError> {
    match (left, right) {
        (Values::Site(left), Values::Site(right)) => {
            let value = match (left.as_slice(), right.as_slice()) {
                ([n11, n21, ..], [n12, n22, ..]) => probability(n11, n12, n21, n22)?,
                _ => Atom::Missing,
            };
            Ok(Values::Site(vec![value]))
        }
        (Values::Samples(left), Values::Samples(right)) => {
            if left.values.len() != right.values.len() {
                return Err(EvaluateError::new("FISHER sample counts differ"));
            }
            let explicit = expressions.len() == 2
                && expressions.iter().all(explicit_reference_pair)
                && sample_width(&left.values) == 2
                && sample_width(&right.values) == 2;
            if !explicit && !expressions.iter().all(reference_format_field) {
                return Err(EvaluateError::new(
                    "two-argument FORMAT FISHER requires Number=R fields",
                ));
            }
            samples(left, right, explicit, record)
        }
        _ => Err(EvaluateError::new(
            "FISHER arguments must both be INFO or both be FORMAT values",
        )),
    }
}

fn samples<'a>(
    left: SampleValues<'a>,
    right: SampleValues<'a>,
    explicit: bool,
    record: &RecordBuf,
) -> Result<Values<'a>, EvaluateError> {
    let selected = left
        .selected
        .iter()
        .zip(&right.selected)
        .map(|(left, right)| *left && *right)
        .collect::<Vec<_>>()
        .into_boxed_slice();
    let genotypes = if explicit {
        None
    } else {
        Some(
            value::diploid_genotype_indices(record)
                .map_err(|error| EvaluateError::new(error.to_string()))?,
        )
    };
    if let Some(genotypes) = &genotypes
        && genotypes.len() != left.values.len()
    {
        return Err(EvaluateError::new(
            "FISHER genotype and FORMAT sample counts differ",
        ));
    }
    let values = left
        .values
        .iter()
        .zip(&right.values)
        .zip(&selected)
        .enumerate()
        .map(|(index, ((left, right), selected))| {
            if !selected {
                return Ok(vec![Atom::Missing]);
            }
            let atoms = if explicit {
                match (left.as_slice(), right.as_slice()) {
                    ([n11, n12], [n21, n22]) => Some((n11, n12, n21, n22)),
                    _ => None,
                }
            } else {
                let Some([first, second]) = genotypes.as_ref().expect("genotypes")[index] else {
                    return Ok(vec![Atom::Missing]);
                };
                match (
                    left.get(first),
                    left.get(second),
                    right.get(first),
                    right.get(second),
                ) {
                    (Some(n11), Some(n12), Some(n21), Some(n22)) => Some((n11, n12, n21, n22)),
                    _ => None,
                }
            };
            match atoms {
                Some((n11, n12, n21, n22)) => {
                    probability(n11, n12, n21, n22).map(|value| vec![value])
                }
                None => Ok(vec![Atom::Missing]),
            }
        })
        .collect::<Result<_, _>>()?;
    Ok(Values::Samples(SampleValues { values, selected }))
}

fn explicit_reference_pair(expression: &BoundExpression) -> bool {
    let Some(field) = reference_field(expression) else {
        return false;
    };
    matches!(
        field.subscript,
        Some(BoundSubscript::SampleValues {
            values: ValueSelector::Indices(_),
            ..
        })
    )
}

fn reference_format_field(expression: &BoundExpression) -> bool {
    reference_field(expression).is_some()
}

fn reference_field(expression: &BoundExpression) -> Option<&crate::expression::bind::BoundField> {
    let BoundExpression::Value(BoundValue::Field(field)) = expression else {
        return None;
    };
    (matches!(field.kind, FieldKind::Format(_))
        && field.cardinality == Cardinality::ReferenceAlternateBases)
        .then_some(field)
}

fn probability<'a>(
    n11: &Atom<'_>,
    n12: &Atom<'_>,
    n21: &Atom<'_>,
    n22: &Atom<'_>,
) -> Result<Atom<'a>, EvaluateError> {
    let (Some(n11), Some(n12), Some(n21), Some(n22)) = (
        numeric_count(n11, "FISHER")?,
        numeric_count(n12, "FISHER")?,
        numeric_count(n21, "FISHER")?,
        numeric_count(n22, "FISHER")?,
    ) else {
        return Ok(Atom::Missing);
    };
    Ok(Atom::Number(fisher_two_sided(
        i64::from(n11),
        i64::from(n12),
        i64::from(n21),
        i64::from(n22),
    )))
}

fn fisher_two_sided(n11: i64, n12: i64, n21: i64, n22: i64) -> f64 {
    let row = n11 + n12;
    let column = n11 + n21;
    let total = n11 + n12 + n21 + n22;
    let maximum = row.min(column);
    let minimum = (row + column - total).max(0);
    if minimum == maximum {
        return 1.0;
    }
    let mut accumulator = Hypergeometric::default();
    let observed = accumulator.probability(n11, row, column, total);
    if observed == 0.0 {
        return 0.0;
    }

    let mut probability = accumulator.probability(minimum, 0, 0, 0);
    let mut left = 0.0;
    let mut value = minimum + 1;
    while probability < 0.999_999_99 * observed && value <= maximum {
        left += probability;
        probability = accumulator.probability(value, 0, 0, 0);
        value += 1;
    }
    if probability < 1.000_000_01 * observed {
        left += probability;
    }

    let mut probability = accumulator.probability(maximum, 0, 0, 0);
    let mut right = 0.0;
    let mut value = maximum - 1;
    while probability < 0.999_999_99 * observed && value >= 0 {
        right += probability;
        probability = accumulator.probability(value, 0, 0, 0);
        value -= 1;
    }
    if probability < 1.000_000_01 * observed {
        right += probability;
    }
    (left + right).min(1.0)
}

#[derive(Default)]
struct Hypergeometric {
    n11: i64,
    row: i64,
    column: i64,
    total: i64,
    value: f64,
}

impl Hypergeometric {
    fn probability(&mut self, n11: i64, row: i64, column: i64, total: i64) -> f64 {
        if row != 0 || column != 0 || total != 0 {
            self.n11 = n11;
            self.row = row;
            self.column = column;
            self.total = total;
        } else if n11 % 11 != 0 && n11 + self.total - self.row - self.column != 0 {
            if n11 == self.n11 + 1 {
                self.value *= (self.row - self.n11) as f64 / n11 as f64
                    * (self.column - self.n11) as f64
                    / (n11 + self.total - self.row - self.column) as f64;
                self.n11 = n11;
                return self.value;
            }
            if n11 == self.n11 - 1 {
                self.value *= self.n11 as f64 / (self.row - n11) as f64
                    * (self.n11 + self.total - self.row - self.column) as f64
                    / (self.column - n11) as f64;
                self.n11 = n11;
                return self.value;
            }
            self.n11 = n11;
        } else {
            self.n11 = n11;
        }
        self.value = hypergeometric(self.n11, self.row, self.column, self.total);
        self.value
    }
}

fn hypergeometric(n11: i64, row: i64, column: i64, total: i64) -> f64 {
    (log_binomial(row, n11) + log_binomial(total - row, column - n11) - log_binomial(total, column))
        .exp()
}

fn log_binomial(total: i64, selected: i64) -> f64 {
    if selected == 0 || total == selected {
        0.0
    } else {
        log_gamma(total as f64 + 1.0)
            - log_gamma(selected as f64 + 1.0)
            - log_gamma((total - selected) as f64 + 1.0)
    }
}

fn log_gamma(value: f64) -> f64 {
    const G: f64 = 607.0 / 128.0;
    const COEFFICIENTS: [f64; 15] = [
        0.999_999_999_999_997_1,
        57.156_235_665_862_92,
        -59.597_960_355_475_49,
        14.136_097_974_741_746,
        -0.491_913_816_097_620_2,
        3.399_464_998_481_189e-5,
        4.652_362_892_704_858e-5,
        -9.837_447_530_487_956e-5,
        1.580_887_032_249_125e-4,
        -2.102_644_417_241_049e-4,
        2.174_396_181_152_126e-4,
        -1.643_181_065_367_639e-4,
        8.441_822_398_385_275e-5,
        -2.619_083_840_158_141e-5,
        3.689_918_265_953_162e-6,
    ];
    let shifted = value - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (shifted + index as f64);
    }
    let scaled = shifted + G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (shifted + 0.5) * scaled.ln() - scaled + series.ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_htslib_reference_probabilities() {
        assert_close(fisher_two_sided(2, 1, 0, 31), 0.005_347_593_583);
        assert_close(fisher_two_sided(2, 1, 0, 1), 1.0);
        assert_close(fisher_two_sided(3, 15, 37, 45), 0.033_161_943_699);
        assert_close(fisher_two_sided(12, 5, 29, 2), 0.080_268_552_074);
        assert_eq!(fisher_two_sided(781, 23_171, 4_963, 2_455_001), 0.0);
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
    }
}
