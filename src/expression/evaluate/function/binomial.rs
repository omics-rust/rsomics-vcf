use noodles_vcf::variant::RecordBuf;

use crate::expression::value::{self, Atom, SampleValues, Values};

use super::super::{EvaluateError, sample_width};

pub(super) fn evaluate<'a>(
    mut arguments: Vec<Values<'a>>,
    record: &RecordBuf,
) -> Result<Values<'a>, EvaluateError> {
    match arguments.len() {
        1 => one(arguments.pop().expect("single BINOM argument"), record),
        2 => {
            let right = arguments.pop().expect("second BINOM argument");
            let left = arguments.pop().expect("first BINOM argument");
            two(left, right)
        }
        _ => Err(EvaluateError::new("BINOM requires one or two arguments")),
    }
}

fn one<'a>(values: Values<'a>, record: &RecordBuf) -> Result<Values<'a>, EvaluateError> {
    match values {
        Values::Site(values) => {
            let value = match values.as_slice() {
                [left, right] => probability(left, right)?,
                _ => Atom::Missing,
            };
            Ok(Values::Site(vec![value]))
        }
        Values::Samples(samples) => {
            let genotypes = value::diploid_genotype_indices(record)
                .map_err(|error| EvaluateError::new(error.to_string()))?;
            if samples.values.len() != genotypes.len() {
                return Err(EvaluateError::new(
                    "BINOM genotype and FORMAT sample counts differ",
                ));
            }
            let values = samples
                .values
                .iter()
                .zip(&samples.selected)
                .zip(genotypes)
                .map(|((counts, selected), genotype)| {
                    if !selected {
                        return Ok(vec![Atom::Missing]);
                    }
                    let Some([left, right]) = genotype else {
                        return Ok(vec![Atom::Missing]);
                    };
                    match (counts.get(left), counts.get(right)) {
                        (Some(left), Some(right)) => {
                            probability(left, right).map(|value| vec![value])
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

fn two<'a>(left: Values<'a>, right: Values<'a>) -> Result<Values<'a>, EvaluateError> {
    match (left, right) {
        (Values::Site(left), Values::Site(right)) => {
            let value = match (left.as_slice(), right.as_slice()) {
                ([left], [right]) => probability(left, right)?,
                _ => Atom::Missing,
            };
            Ok(Values::Site(vec![value]))
        }
        (Values::Samples(left), Values::Samples(right)) => {
            if left.values.len() != right.values.len() {
                return Err(EvaluateError::new("BINOM sample counts differ"));
            }
            if sample_width(&left.values) != 1 || sample_width(&right.values) != 1 {
                return Err(EvaluateError::new(
                    "BINOM requires one value per explicit FORMAT argument",
                ));
            }
            let selected = left
                .selected
                .iter()
                .zip(&right.selected)
                .map(|(left, right)| *left && *right)
                .collect::<Vec<_>>()
                .into_boxed_slice();
            let values = left
                .values
                .iter()
                .zip(&right.values)
                .zip(&selected)
                .map(|((left, right), selected)| {
                    if !selected {
                        Ok(vec![Atom::Missing])
                    } else {
                        match (left.first(), right.first()) {
                            (Some(left), Some(right)) => {
                                probability(left, right).map(|value| vec![value])
                            }
                            _ => Ok(vec![Atom::Missing]),
                        }
                    }
                })
                .collect::<Result<_, _>>()?;
            Ok(Values::Samples(SampleValues { values, selected }))
        }
        _ => Err(EvaluateError::new(
            "BINOM arguments must both be INFO or both be FORMAT values",
        )),
    }
}

fn probability<'a>(left: &Atom<'_>, right: &Atom<'_>) -> Result<Atom<'a>, EvaluateError> {
    let (Some(left), Some(right)) = (count(left)?, count(right)?) else {
        return Ok(Atom::Missing);
    };
    Ok(binomial_two_sided(left, right)
        .map(Atom::Number)
        .unwrap_or(Atom::Missing))
}

fn count(value: &Atom<'_>) -> Result<Option<i32>, EvaluateError> {
    match value {
        Atom::Absent | Atom::Missing => Ok(None),
        Atom::Number(value) if value.is_finite() && *value >= 0.0 && *value <= i32::MAX as f64 => {
            Ok(Some(*value as i32))
        }
        Atom::Flag => Ok(Some(1)),
        Atom::Number(_) => Err(EvaluateError::new(
            "BINOM counts must be finite nonnegative 32-bit values",
        )),
        _ => Err(EvaluateError::new("BINOM received a nonnumeric value")),
    }
}

fn binomial_two_sided(left: i32, right: i32) -> Option<f64> {
    if left == 0 && right == 0 {
        None
    } else if left == right {
        Some(1.0)
    } else {
        let probability = if left > right {
            2.0 * regularized_beta(left as f64, right as f64 + 1.0, 0.5)
        } else {
            2.0 * regularized_beta(right as f64, left as f64 + 1.0, 0.5)
        };
        Some(probability.min(1.0))
    }
}

fn regularized_beta(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        beta_fraction(a, b, x)
    } else {
        1.0 - beta_fraction(b, a, 1.0 - x)
    }
}

fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let mut fraction = 1.0;
    let mut c = fraction;
    let mut d = 0.0;
    for index in 1..200 {
        let middle = index >> 1;
        let middle = middle as f64;
        let coefficient = if index & 1 == 1 {
            -(a + middle) * (a + b + middle) * x / ((a + 2.0 * middle) * (a + 2.0 * middle + 1.0))
        } else {
            middle * (b - middle) * x / ((a + 2.0 * middle - 1.0) * (a + 2.0 * middle))
        };
        d = 1.0 + coefficient * d;
        if d < 1e-290 {
            d = 1e-290;
        }
        c = 1.0 + coefficient / c;
        if c < 1e-290 {
            c = 1e-290;
        }
        d = 1.0 / d;
        let delta = c * d;
        fraction *= delta;
        if (delta - 1.0).abs() < 1e-14 {
            break;
        }
    }
    (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp()
        / a
        / fraction
}

fn log_gamma(value: f64) -> f64 {
    let mut sum = 0.0;
    sum += 0.165_947_018_740_846_2e-6 / (value + 7.0);
    sum += 0.993_493_711_393_074_8e-5 / (value + 6.0);
    sum -= 0.138_571_033_129_652_6 / (value + 5.0);
    sum += 12.507_343_240_090_56 / (value + 4.0);
    sum -= 176.615_029_149_838_6 / (value + 3.0);
    sum += 771.323_428_775_767_4 / (value + 2.0);
    sum -= 1_259.139_216_722_289 / (value + 1.0);
    sum += 676.520_368_121_883_5 / value;
    sum += 0.999_999_999_999_518_3;
    sum.ln() - 5.581_061_466_795_328 - value + (value - 0.5) * (value + 6.5).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_reference_two_sided_probabilities() {
        assert_close(binomial_two_sided(3, 5), 0.7265625);
        assert_close(binomial_two_sided(0, 3), 0.25);
        assert_close(binomial_two_sided(2, 0), 0.5);
        assert_close(binomial_two_sided(4, 4), 1.0);
        assert_eq!(binomial_two_sided(0, 0), None);
    }

    fn assert_close(actual: Option<f64>, expected: f64) {
        assert!((actual.expect("probability") - expected).abs() < 1e-12);
    }
}
