use std::collections::HashMap;

use super::Diagnostics;

pub(crate) fn inspect(columns: &[&str], line: usize, diagnostics: &mut Diagnostics) {
    let alternates: Vec<_> = columns[4].split(',').collect();
    let info = parse_info(columns[7]);

    inspect_structural_variants(&alternates, &info, line, diagnostics);
    inspect_repeats(&alternates, &info, line, diagnostics);
    inspect_info_intervals(&alternates, &info, line, diagnostics);

    if columns.len() > 8 {
        inspect_format(&alternates, &info, &columns[8..], line, diagnostics);
    }
}

fn parse_info(value: &str) -> HashMap<&str, &str> {
    value
        .split(';')
        .filter_map(|item| item.split_once('='))
        .collect()
}

fn inspect_structural_variants(
    alternates: &[&str],
    info: &HashMap<&str, &str>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let mut symbolic = false;
    for alternate in alternates {
        let Some(inner) = symbolic_inner(alternate) else {
            continue;
        };
        let kind = inner.split(':').next().unwrap_or(inner);
        symbolic |= matches!(kind, "DEL" | "INS" | "DUP" | "INV" | "CNV");
        if inner.contains(':') && !matches!(kind, "DEL" | "INS" | "DUP" | "INV" | "CNV") {
            diagnostics.error(
                "record.symbolic-alt",
                Some(line),
                Some("ALT"),
                format!("symbolic ALT type {kind:?} is not defined by VCF 4.4"),
            );
        }
    }

    if symbolic && !info.contains_key("SVLEN") {
        diagnostics.error(
            "record.svlen-required",
            Some(line),
            Some("INFO/SVLEN"),
            "symbolic structural variants require INFO/SVLEN",
        );
    }

    let claims = info.get("SVCLAIM").map(|value| values(value));
    if alternates
        .iter()
        .any(|alternate| symbolic_kind(alternate).is_some_and(|kind| matches!(kind, "DEL" | "DUP")))
        && claims.is_none()
    {
        diagnostics.error(
            "record.svclaim-required",
            Some(line),
            Some("INFO/SVCLAIM"),
            "DEL and DUP alleles require INFO/SVCLAIM",
        );
    }
    let Some(claims) = claims else {
        return;
    };
    if claims.len() != alternates.len() {
        diagnostics.error(
            "record.svclaim-cardinality",
            Some(line),
            Some("INFO/SVCLAIM"),
            format!(
                "SVCLAIM has {} values; ALT has {} alleles",
                claims.len(),
                alternates.len()
            ),
        );
        return;
    }

    for (alternate, claim) in alternates.iter().zip(claims) {
        let valid = match symbolic_kind(alternate) {
            Some("DEL" | "DUP") => matches!(claim, "D" | "J" | "DJ"),
            Some("CNV") => matches!(claim, "D" | "."),
            Some("INV" | "INS") => matches!(claim, "J" | "DJ" | "."),
            Some(_) => claim == ".",
            None if is_breakend(alternate) => matches!(claim, "J" | "."),
            None => claim == ".",
        };
        if !valid {
            diagnostics.error(
                "record.svclaim",
                Some(line),
                Some("INFO/SVCLAIM"),
                format!("SVCLAIM value {claim:?} is invalid for ALT allele {alternate}"),
            );
        }
    }
}

fn inspect_repeats(
    alternates: &[&str],
    info: &HashMap<&str, &str>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let tandem_repeats = alternates
        .iter()
        .filter(|alternate| symbolic_inner(alternate) == Some("CNV:TR"))
        .count();
    if tandem_repeats == 0 {
        return;
    }

    let repeat_count = info.get("RN").map_or(tandem_repeats, |raw| {
        values(raw)
            .into_iter()
            .filter_map(|value| value.parse::<usize>().ok())
            .sum()
    });
    let rus = info.get("RUS").map(|raw| values(raw));
    let rul = info.get("RUL").map(|raw| values(raw));

    if rus.as_ref().is_none_or(|values| all_missing(values))
        && rul.as_ref().is_none_or(|values| all_missing(values))
    {
        diagnostics.error(
            "record.repeat-unit-required",
            Some(line),
            Some("INFO"),
            "CNV:TR alleles require INFO/RUS or INFO/RUL",
        );
    }
    if rus.is_some() && rul.is_some() {
        diagnostics.error(
            "record.repeat-unit-redundant",
            Some(line),
            Some("INFO"),
            "RUS and RUL must not be present together",
        );
    }

    for key in ["RUS", "RUL", "RUC", "RB"] {
        if let Some(raw) = info.get(key) {
            let count = values(raw).len();
            if count != repeat_count {
                diagnostics.error(
                    "record.repeat-cardinality",
                    Some(line),
                    Some(key),
                    format!("{key} has {count} values; RN requires {repeat_count}"),
                );
            }
        }
    }

    if let (Some(rus), Some(rul)) = (&rus, &rul) {
        for (sequence, length) in rus.iter().zip(rul) {
            if *sequence != "."
                && length
                    .parse::<usize>()
                    .is_ok_and(|length| length != sequence.len())
            {
                diagnostics.error(
                    "record.repeat-unit-length",
                    Some(line),
                    Some("INFO/RUL"),
                    format!("RUL value {length} does not match RUS value {sequence:?}"),
                );
            }
        }
    }

    inspect_repeat_interval("CIRUC", "RUC", info, line, diagnostics);
    inspect_repeat_interval("CIRB", "RB", info, line, diagnostics);
    inspect_repeat_bases(info, rus.as_deref(), rul.as_deref(), line, diagnostics);
    inspect_repeat_units(info, line, diagnostics);
}

fn inspect_repeat_interval(
    key: &str,
    parent: &str,
    info: &HashMap<&str, &str>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let Some(raw) = info.get(key) else {
        return;
    };
    let interval = values(raw);
    let Some(parent_raw) = info.get(parent) else {
        diagnostics.error(
            "record.repeat-interval-parent",
            Some(line),
            Some(key),
            format!("{key} requires {parent}"),
        );
        return;
    };
    let parent_values = values(parent_raw);
    if interval.len() != parent_values.len() * 2 {
        diagnostics.error(
            "record.repeat-interval-cardinality",
            Some(line),
            Some(key),
            format!(
                "{key} has {} values; {parent} requires {}",
                interval.len(),
                parent_values.len() * 2
            ),
        );
        return;
    }
    for (parent_value, bounds) in parent_values.iter().zip(interval.chunks_exact(2)) {
        if *parent_value == "." && bounds != [".", "."] {
            diagnostics.error(
                "record.repeat-interval-missing",
                Some(line),
                Some(key),
                format!("{key} bounds must be missing when {parent} is missing"),
            );
        }
    }
    inspect_interval(key, &interval, line, diagnostics);
}

fn inspect_repeat_bases(
    info: &HashMap<&str, &str>,
    rus: Option<&[&str]>,
    rul: Option<&[&str]>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let (Some(rb), Some(ruc)) = (info.get("RB"), info.get("RUC")) else {
        return;
    };
    let lengths = rul
        .map(|values| {
            values
                .iter()
                .map(|value| value.parse::<f64>().ok())
                .collect::<Vec<_>>()
        })
        .or_else(|| {
            rus.map(|values| {
                values
                    .iter()
                    .map(|value| (value != &".").then_some(value.len() as f64))
                    .collect()
            })
        });
    let Some(lengths) = lengths else {
        return;
    };

    for ((rb, ruc), length) in values(rb).iter().zip(values(ruc)).zip(lengths) {
        let (Ok(rb), Ok(ruc), Some(length)) = (rb.parse::<f64>(), ruc.parse::<f64>(), length)
        else {
            continue;
        };
        let expected = ruc * length;
        if (rb - expected).abs() / rb.abs().max(1.0) > 0.05 {
            diagnostics.error(
                "record.repeat-bases",
                Some(line),
                Some("INFO/RB"),
                format!("RB value {rb} differs from RUC × repeat-unit length {expected}"),
            );
        }
    }
}

fn inspect_repeat_units(info: &HashMap<&str, &str>, line: usize, diagnostics: &mut Diagnostics) {
    let Some(rub) = info.get("RUB") else {
        return;
    };
    let Some(ruc) = info.get("RUC") else {
        diagnostics.error(
            "record.rub-parent",
            Some(line),
            Some("INFO/RUB"),
            "RUB requires RUC",
        );
        return;
    };
    let mut expected = 0usize;
    let mut valid = true;
    for count in values(ruc) {
        if count == "." {
            continue;
        }
        match count.parse::<usize>() {
            Ok(count) => expected = expected.saturating_add(count),
            Err(_) => valid = false,
        }
    }
    if !valid {
        diagnostics.error(
            "record.rub-ruc",
            Some(line),
            Some("INFO/RUC"),
            "RUC values must be integers when RUB is present",
        );
    }
    let actual = values(rub).len();
    if actual != expected {
        diagnostics.error(
            "record.rub-cardinality",
            Some(line),
            Some("INFO/RUB"),
            format!("RUB has {actual} values; RUC requires {expected}"),
        );
    }
}

fn inspect_info_intervals(
    alternates: &[&str],
    info: &HashMap<&str, &str>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let Some(raw) = info.get("CICN") else {
        return;
    };
    let interval = values(raw);
    if interval.len() != alternates.len() * 2 {
        diagnostics.error(
            "record.cicn-cardinality",
            Some(line),
            Some("INFO/CICN"),
            format!(
                "CICN has {} values; ALT requires {}",
                interval.len(),
                alternates.len() * 2
            ),
        );
        return;
    }
    inspect_interval("INFO/CICN", &interval, line, diagnostics);
}

fn inspect_interval(key: &str, values: &[&str], line: usize, diagnostics: &mut Diagnostics) {
    for bounds in values.chunks_exact(2) {
        let lower = missing_as_zero(bounds[0]);
        let upper = missing_as_zero(bounds[1]);
        if lower.is_some_and(|value| value > 0.0) || upper.is_some_and(|value| value < 0.0) {
            diagnostics.error(
                "record.confidence-interval",
                Some(line),
                Some(key),
                format!("{key} bounds must contain zero"),
            );
        }
    }
}

fn inspect_format(
    alternates: &[&str],
    info: &HashMap<&str, &str>,
    columns: &[&str],
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let keys: Vec<_> = columns[0].split(':').collect();
    if keys.contains(&"CICN") && !keys.contains(&"CN") {
        diagnostics.error(
            "record.format-cicn-parent",
            Some(line),
            Some("FORMAT/CICN"),
            "FORMAT/CICN requires FORMAT/CN",
        );
    }
    if keys.contains(&"CN") {
        inspect_copy_number_svlen(alternates, info, line, diagnostics);
    }

    for sample in &columns[1..] {
        let fields: HashMap<_, _> = keys.iter().copied().zip(sample.split(':')).collect();
        if let Some(raw) = fields.get("CICN") {
            let interval = values(raw);
            if !interval.len().is_multiple_of(2) {
                diagnostics.error(
                    "record.format-cicn-cardinality",
                    Some(line),
                    Some("FORMAT/CICN"),
                    "FORMAT/CICN must contain pairs of bounds",
                );
            } else {
                inspect_interval("FORMAT/CICN", &interval, line, diagnostics);
            }
        }
        inspect_phasing(&fields, line, diagnostics);
    }
}

fn inspect_copy_number_svlen(
    alternates: &[&str],
    info: &HashMap<&str, &str>,
    line: usize,
    diagnostics: &mut Diagnostics,
) {
    let Some(raw) = info.get("SVLEN") else {
        return;
    };
    let lengths = values(raw);
    if lengths.len() != alternates.len() {
        return;
    }
    let mut expected = None;
    for (alternate, length) in alternates.iter().zip(lengths) {
        if !symbolic_kind(alternate).is_some_and(|kind| matches!(kind, "DEL" | "DUP" | "CNV")) {
            continue;
        }
        if let Some(expected) = expected {
            if length != expected {
                diagnostics.error(
                    "record.copy-number-svlen",
                    Some(line),
                    Some("INFO/SVLEN"),
                    "DEL, DUP, and CNV alleles must have the same SVLEN when FORMAT/CN is present",
                );
                return;
            }
        } else {
            expected = Some(length);
        }
    }
}

fn inspect_phasing(fields: &HashMap<&str, &str>, line: usize, diagnostics: &mut Diagnostics) {
    let Some(psl) = fields.get("PSL") else {
        return;
    };
    let psl = values(psl);
    if let Some(gt) = fields.get("GT") {
        for (phase, value) in phases(gt).into_iter().zip(&psl) {
            if phase == '/' && *value != "." {
                diagnostics.error(
                    "record.psl-unphased",
                    Some(line),
                    Some("FORMAT/PSL"),
                    "PSL must be missing for unphased genotype alleles",
                );
            }
        }
    }
    for key in ["PSO", "PSQ"] {
        let Some(raw) = fields.get(key) else {
            continue;
        };
        for (psl, value) in psl.iter().zip(values(raw)) {
            if *psl == "." && value != "." {
                diagnostics.error(
                    "record.phase-value",
                    Some(line),
                    Some(key),
                    format!("{key} must be missing when PSL is missing"),
                );
            }
        }
    }
}

fn phases(genotype: &str) -> Vec<char> {
    let genotype = genotype.strip_suffix(['/', '|']).unwrap_or(genotype);
    let first = if genotype.starts_with(['/', '|']) {
        genotype.as_bytes()[0] as char
    } else if genotype.get(1..).is_some_and(|value| value.contains('/')) {
        '/'
    } else {
        '|'
    };
    std::iter::once(first)
        .chain(
            genotype
                .chars()
                .filter(|character| matches!(character, '/' | '|'))
                .skip(usize::from(genotype.starts_with(['/', '|']))),
        )
        .collect()
}

fn values(value: &str) -> Vec<&str> {
    value.split(',').collect()
}

fn all_missing(values: &[&str]) -> bool {
    values.iter().all(|value| *value == ".")
}

fn missing_as_zero(value: &str) -> Option<f64> {
    if value == "." {
        Some(0.0)
    } else {
        value.parse().ok()
    }
}

fn symbolic_inner(alternate: &str) -> Option<&str> {
    alternate
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
}

fn symbolic_kind(alternate: &str) -> Option<&str> {
    symbolic_inner(alternate).map(|value| value.split(':').next().unwrap_or(value))
}

fn is_breakend(alternate: &str) -> bool {
    alternate.contains(['[', ']'])
}
