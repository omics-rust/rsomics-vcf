use std::collections::HashSet;

use noodles_vcf as vcf;

use super::{
    Diagnostics, Field, Number, Schema, Type, format_definition, info_definition,
    inspect_v44_record,
};

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct VariantKey {
    pub(crate) position: u64,
    pub(crate) reference: String,
    pub(crate) alternate: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct RecordPosition<'a> {
    pub(crate) chrom: &'a str,
    pub(crate) position: u64,
    pub(crate) variants: Vec<VariantKey>,
}

pub(crate) fn inspect_record<'a>(
    line: &'a [u8],
    line_number: usize,
    schema: &Schema,
    header: Option<&vcf::Header>,
    require_evidence: bool,
    diagnostics: &mut Diagnostics,
) -> Option<RecordPosition<'a>> {
    let Ok(line) = std::str::from_utf8(line) else {
        diagnostics.error(
            "record.utf8",
            Some(line_number),
            None,
            "record is not valid UTF-8",
        );
        return None;
    };
    if line.is_empty() {
        diagnostics.error(
            "record.empty",
            Some(line_number),
            None,
            "empty lines are not valid variant records",
        );
        return None;
    }

    let columns: Vec<_> = line.split('\t').collect();
    let expected = if schema.samples.is_empty() {
        8
    } else {
        9 + schema.samples.len()
    };
    if columns.len() != expected {
        diagnostics.error(
            "record.columns",
            Some(line_number),
            None,
            format!(
                "record has {} columns; header requires {expected}",
                columns.len()
            ),
        );
        return None;
    }

    let chrom = columns[0];
    if !valid_chrom(chrom, schema.version) {
        diagnostics.error(
            "record.chrom",
            Some(line_number),
            Some("CHROM"),
            "CHROM is not a valid contig identifier",
        );
    }
    let position = parse_position(columns[1], line_number, diagnostics);
    inspect_ids(columns[2], line_number, diagnostics);
    if !valid_base_sequence(columns[3]) {
        diagnostics.error(
            "record.ref",
            Some(line_number),
            Some("REF"),
            "REF must contain only A, C, G, T, or N bases",
        );
    }
    let alternates = inspect_alternates(columns[4], schema.version, line_number, diagnostics);
    inspect_quality(columns[5], line_number, diagnostics);
    inspect_filters(columns[6], line_number, diagnostics);
    let evidence = inspect_info(
        columns[7],
        line_number,
        alternates,
        schema,
        position,
        chrom,
        diagnostics,
    );

    if !schema.samples.is_empty() {
        let has_genotype =
            inspect_samples(&columns[8..], line_number, alternates, schema, diagnostics);
        if require_evidence && !has_genotype && !evidence {
            diagnostics.error(
                "record.evidence",
                Some(line_number),
                None,
                "record has no GT samples, INFO/AF, or INFO/AC with INFO/AN",
            );
        }
    } else if require_evidence && !evidence {
        diagnostics.error(
            "record.evidence",
            Some(line_number),
            None,
            "record has no INFO/AF or INFO/AC with INFO/AN",
        );
    }

    if schema.version.is_some_and(|version| version >= (4, 4)) {
        inspect_v44_record(&columns, line_number, diagnostics);
    }

    if let Some(header) = header {
        let mut typed = vcf::variant::RecordBuf::default();
        let mut reader = vcf::io::Reader::new(line.as_bytes());
        if let Err(error) = reader.read_record_buf(header, &mut typed) {
            diagnostics.error(
                "record.typed-value",
                Some(line_number),
                None,
                error_chain(&error),
            );
        }
    }

    position.map(|position| RecordPosition {
        chrom: chrom
            .strip_prefix('<')
            .and_then(|chrom| chrom.strip_suffix('>'))
            .unwrap_or(chrom),
        position,
        variants: variant_keys(position, columns[3], columns[4]),
    })
}

fn parse_position(value: &str, line_number: usize, diagnostics: &mut Diagnostics) -> Option<u64> {
    match value.parse::<u64>() {
        Ok(position) if position <= i32::MAX as u64 => Some(position),
        _ => {
            diagnostics.error(
                "record.pos",
                Some(line_number),
                Some("POS"),
                "POS must be an integer between 0 and 2147483647",
            );
            None
        }
    }
}

fn inspect_ids(value: &str, line_number: usize, diagnostics: &mut Diagnostics) {
    if value == "." {
        return;
    }
    let mut seen = HashSet::new();
    for id in value.split(';') {
        if id.is_empty() || id.bytes().any(|byte| byte.is_ascii_whitespace()) {
            diagnostics.error(
                "record.id",
                Some(line_number),
                Some("ID"),
                "ID values must be nonempty and contain no whitespace",
            );
            return;
        }
        if !seen.insert(id) {
            diagnostics.error(
                "record.duplicate-id",
                Some(line_number),
                Some("ID"),
                format!("ID {id} is duplicated"),
            );
        }
    }
}

fn inspect_alternates(
    value: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> usize {
    if value == "." {
        return 0;
    }
    let mut seen = HashSet::new();
    let mut count = 0;
    for alternate in value.split(',') {
        count += 1;
        if !valid_alternate(alternate, version) {
            diagnostics.error(
                "record.alt",
                Some(line_number),
                Some("ALT"),
                format!("ALT allele {alternate:?} is invalid"),
            );
        } else if (valid_base_sequence(alternate) || version.is_none_or(|version| version < (4, 4)))
            && !seen.insert(alternate)
        {
            diagnostics.error(
                "record.duplicate-alt",
                Some(line_number),
                Some("ALT"),
                format!("ALT allele {alternate} is duplicated"),
            );
        }
    }
    count
}

fn inspect_quality(value: &str, line_number: usize, diagnostics: &mut Diagnostics) {
    if value == "." {
        return;
    }
    match parse_float(value) {
        Some(value) if value.is_nan() || value >= 0.0 => {}
        _ => diagnostics.error(
            "record.qual",
            Some(line_number),
            Some("QUAL"),
            "QUAL must be a nonnegative number or missing",
        ),
    }
}

fn inspect_filters(value: &str, line_number: usize, diagnostics: &mut Diagnostics) {
    if matches!(value, "." | "PASS") {
        return;
    }
    let mut seen = HashSet::new();
    for filter in value.split(';') {
        if filter.is_empty()
            || matches!(filter, "." | "0")
            || filter.bytes().any(|byte| byte.is_ascii_whitespace())
        {
            diagnostics.error(
                "record.filter",
                Some(line_number),
                Some("FILTER"),
                "FILTER must be missing, PASS, or a semicolon-delimited list of filter IDs",
            );
            return;
        }
        if !seen.insert(filter) {
            diagnostics.error(
                "record.duplicate-filter",
                Some(line_number),
                Some("FILTER"),
                format!("FILTER value {filter} is duplicated"),
            );
        }
    }
}

fn inspect_info(
    value: &str,
    line_number: usize,
    alternates: usize,
    schema: &Schema,
    position: Option<u64>,
    chrom: &str,
    diagnostics: &mut Diagnostics,
) -> bool {
    if value == "." {
        return false;
    }
    if value.is_empty() {
        diagnostics.error(
            "record.info",
            Some(line_number),
            Some("INFO"),
            "INFO must be missing or contain at least one field",
        );
        return false;
    }

    let mut seen = HashSet::new();
    let mut has_af = false;
    let mut has_ac = false;
    let mut has_an = false;
    for item in value.split(';') {
        let (key, raw) = item
            .split_once('=')
            .map_or((item, None), |(key, value)| (key, Some(value)));
        if !valid_info_key(key, schema.version) {
            diagnostics.error(
                "record.info-key",
                Some(line_number),
                Some("INFO"),
                format!("INFO key {key:?} is invalid"),
            );
            continue;
        }
        if !seen.insert(key) {
            diagnostics.error(
                "record.duplicate-info",
                Some(line_number),
                Some("INFO"),
                format!("INFO key {key} is duplicated"),
            );
            continue;
        }
        has_af |= key == "AF";
        has_ac |= key == "AC";
        has_an |= key == "AN";

        let field = schema
            .info
            .get(key)
            .copied()
            .or_else(|| info_definition(key, schema.version));
        if let Some(field) = field {
            inspect_value(
                key,
                raw,
                field,
                alternates,
                None,
                None,
                true,
                schema.version.is_some_and(|version| version >= (4, 5)),
                line_number,
                diagnostics,
            );
        }

        if key == "END"
            && let Some(raw) = raw
            && let Ok(end) = raw.parse::<u64>()
            && let Some(Some(length)) = schema.contigs.get(chrom)
            && end > *length
        {
            diagnostics.error(
                "record.contig-bound",
                Some(line_number),
                Some("INFO/END"),
                format!("INFO/END {end} exceeds contig {chrom} length {length}"),
            );
        }
    }

    if let Some(position) = position
        && let Some(Some(length)) = schema.contigs.get(chrom)
        && position > *length + 1
    {
        diagnostics.error(
            "record.contig-bound",
            Some(line_number),
            Some("POS"),
            format!("POS {position} exceeds contig {chrom} length {length}"),
        );
    }
    has_af || (has_ac && has_an)
}

fn inspect_samples(
    columns: &[&str],
    line_number: usize,
    alternates: usize,
    schema: &Schema,
    diagnostics: &mut Diagnostics,
) -> bool {
    let format = columns[0];
    if format.is_empty() || format == "." {
        diagnostics.error(
            "record.format",
            Some(line_number),
            Some("FORMAT"),
            "FORMAT must contain at least one key when samples are present",
        );
        return false;
    }
    let keys: Vec<_> = format.split(':').collect();
    let mut seen = HashSet::new();
    for key in &keys {
        if !valid_format_key(key, schema.version) || !seen.insert(*key) {
            diagnostics.error(
                "record.format-key",
                Some(line_number),
                Some("FORMAT"),
                format!("FORMAT key {key:?} is invalid or duplicated"),
            );
        }
    }
    if let Some(index) = keys.iter().position(|key| *key == "GT")
        && index != 0
    {
        diagnostics.error(
            "record.gt-order",
            Some(line_number),
            Some("FORMAT/GT"),
            "GT must be the first FORMAT field",
        );
    }

    for (sample_index, sample) in columns[1..].iter().enumerate() {
        if *sample == "." {
            continue;
        }
        if sample.is_empty() && schema.version.is_some_and(|version| version >= (4, 5)) {
            continue;
        }
        if sample.is_empty() {
            diagnostics.error(
                "record.sample",
                Some(line_number),
                Some("sample"),
                format!("sample {} is empty", schema.samples[sample_index]),
            );
            continue;
        }
        let values: Vec<_> = sample.split(':').collect();
        if values.len() > keys.len() {
            diagnostics.error(
                "record.sample-width",
                Some(line_number),
                Some("sample"),
                format!(
                    "sample {} has {} values but FORMAT declares {}",
                    schema.samples[sample_index],
                    values.len(),
                    keys.len()
                ),
            );
            continue;
        }
        let ploidy = if keys[0] == "GT" {
            inspect_genotype(
                values.first().copied().unwrap_or("."),
                alternates,
                schema.version,
                line_number,
                &schema.samples[sample_index],
                diagnostics,
            )
        } else {
            2
        };
        let local_alternates = keys
            .iter()
            .position(|key| *key == "LAA")
            .and_then(|index| values.get(index).copied())
            .and_then(|value| {
                inspect_local_alternates(
                    value,
                    alternates,
                    line_number,
                    &schema.samples[sample_index],
                    diagnostics,
                )
            });
        for (key, raw) in keys.iter().zip(values) {
            if *key == "GT" {
                continue;
            }
            let field = schema
                .format
                .get(*key)
                .copied()
                .or_else(|| format_definition(key, schema.version));
            if let Some(field) = field {
                inspect_value(
                    key,
                    Some(raw),
                    field,
                    alternates,
                    local_alternates,
                    Some(ploidy),
                    false,
                    schema.version.is_some_and(|version| version >= (4, 5)),
                    line_number,
                    diagnostics,
                );
            }
        }
    }

    keys.contains(&"GT")
}

fn inspect_genotype(
    value: &str,
    alternates: usize,
    version: Option<(u8, u8)>,
    line_number: usize,
    sample: &str,
    diagnostics: &mut Diagnostics,
) -> usize {
    if value == "." {
        return 1;
    }
    let value = value.strip_prefix(['/', '|']).unwrap_or(value);
    let value = if version.is_some_and(|version| version >= (4, 4)) {
        value.strip_suffix(['/', '|']).unwrap_or(value)
    } else {
        value
    };
    let alleles: Vec<_> = value.split(['/', '|']).collect();
    if alleles.is_empty() || alleles.iter().any(|allele| allele.is_empty()) {
        diagnostics.error(
            "record.genotype",
            Some(line_number),
            Some("FORMAT/GT"),
            format!("sample {sample} has an empty genotype allele"),
        );
        return alleles.len().max(1);
    }
    for allele in &alleles {
        if *allele == "." {
            continue;
        }
        match allele.parse::<usize>() {
            Ok(index) if index <= alternates => {}
            _ => diagnostics.error(
                "record.genotype-allele",
                Some(line_number),
                Some("FORMAT/GT"),
                format!("sample {sample} genotype allele {allele:?} is outside 0..={alternates}"),
            ),
        }
    }
    alleles.len()
}

fn inspect_local_alternates(
    value: &str,
    alternates: usize,
    line_number: usize,
    sample: &str,
    diagnostics: &mut Diagnostics,
) -> Option<usize> {
    if value == "." {
        return None;
    }
    if value.is_empty() {
        return Some(0);
    }
    let mut seen = HashSet::new();
    let values: Vec<_> = value.split(',').collect();
    for value in &values {
        match value.parse::<usize>() {
            Ok(index) if (1..=alternates).contains(&index) && seen.insert(index) => {}
            _ => diagnostics.error(
                "record.local-allele",
                Some(line_number),
                Some("FORMAT/LAA"),
                format!("sample {sample} has invalid local ALT index {value:?}"),
            ),
        }
    }
    Some(values.len())
}

#[allow(clippy::too_many_arguments)]
fn inspect_value(
    key: &str,
    raw: Option<&str>,
    field: Field,
    alternates: usize,
    local_alternates: Option<usize>,
    ploidy: Option<usize>,
    info: bool,
    allow_empty: bool,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) {
    if field.ty == Type::Flag {
        if raw.is_some_and(|value| !matches!(value, "0" | "1")) {
            diagnostics.error(
                "record.flag-value",
                Some(line_number),
                Some(key),
                format!("{key} is a flag and must not have a value"),
            );
        }
        return;
    }
    let Some(raw) = raw else {
        diagnostics.error(
            "record.missing-value",
            Some(line_number),
            Some(key),
            format!("{key} requires a value"),
        );
        return;
    };
    if raw == "." {
        return;
    }
    if raw.is_empty() && !allow_empty && field.number != Number::Count(0) {
        diagnostics.error(
            "record.empty-value",
            Some(line_number),
            Some(key),
            format!("{key} has an empty value"),
        );
        return;
    }

    let values: Vec<_> = if raw.is_empty() {
        Vec::new()
    } else {
        split_values(raw).collect()
    };
    let expected = match field.number {
        Number::Count(count) => Some(count),
        Number::Alternate => Some(alternates),
        Number::Alleles => Some(alternates + 1),
        Number::Genotypes if info => None,
        Number::Genotypes => ploidy.and_then(|ploidy| genotype_count(alternates + 1, ploidy)),
        Number::LocalAlternate => local_alternates,
        Number::LocalAlleles => local_alternates.and_then(|count| count.checked_add(1)),
        Number::LocalGenotypes => local_alternates.and_then(|count| {
            ploidy.and_then(|ploidy| genotype_count(count.saturating_add(1), ploidy))
        }),
        Number::BaseModifications => None,
        Number::Ploidy => ploidy,
        Number::Unknown => None,
    };
    if let Some(expected) = expected
        && values.len() != expected
    {
        diagnostics.error(
            "record.cardinality",
            Some(line_number),
            Some(key),
            format!(
                "{key} has {} values; Number requires {expected}",
                values.len()
            ),
        );
    }

    for value in values {
        if value == "." {
            continue;
        }
        let valid = match field.ty {
            Type::Integer => valid_integer(value),
            Type::Float => parse_float(value).is_some(),
            Type::Character => value.chars().count() == 1,
            Type::String => !value.bytes().any(|byte| byte.is_ascii_whitespace()),
            Type::Flag => unreachable!(),
        };
        if !valid {
            diagnostics.error(
                "record.value-type",
                Some(line_number),
                Some(key),
                format!("{key} value {value:?} does not match {:?}", field.ty),
            );
        }
        if matches!(key, "AC" | "AN" | "DP" | "END" | "MQ0" | "NS")
            && value.parse::<i64>().is_ok_and(|value| value < 0)
        {
            diagnostics.error(
                "record.value-range",
                Some(line_number),
                Some(key),
                format!("{key} must not be negative"),
            );
        }
        if matches!(key, "AF" | "GP")
            && parse_float(value).is_some_and(|value| !(0.0..=1.0).contains(&value))
        {
            diagnostics.error(
                "record.value-range",
                Some(line_number),
                Some(key),
                format!("{key} must be between 0 and 1"),
            );
        }
        if key == "CIGAR" && !valid_cigar(value) {
            diagnostics.error(
                "record.cigar",
                Some(line_number),
                Some(key),
                format!("CIGAR value {value:?} is invalid"),
            );
        }
    }
}

fn genotype_count(alleles: usize, ploidy: usize) -> Option<usize> {
    if ploidy == 0 {
        return Some(0);
    }
    let choose = ploidy.checked_sub(1)?.checked_add(alleles)?;
    binomial(choose, ploidy)
}

fn binomial(n: usize, k: usize) -> Option<usize> {
    let k = k.min(n.checked_sub(k)?);
    let mut value = 1usize;
    for index in 0..k {
        value = value.checked_mul(n - index)?.checked_div(index + 1)?;
    }
    Some(value)
}

fn valid_info_key(key: &str, version: Option<(u8, u8)>) -> bool {
    if version.is_some_and(|version| version >= (4, 3)) {
        key == "1000G" || valid_encoded_id(key, version, true)
    } else {
        valid_encoded_id(key, version, false)
    }
}

fn valid_format_key(key: &str, version: Option<(u8, u8)>) -> bool {
    valid_encoded_id(key, version, true)
}

fn valid_encoded_id(key: &str, version: Option<(u8, u8)>, require_alpha: bool) -> bool {
    let bytes = key.as_bytes();
    if bytes.is_empty() || require_alpha && !bytes[0].is_ascii_alphabetic() && bytes[0] != b'_' {
        return false;
    }
    let encoded = version.is_some_and(|version| version >= (4, 4));
    let underscored = version.is_some_and(|version| version >= (4, 3));
    let mut index = 0;
    while index < bytes.len() {
        let byte = bytes[index];
        if byte.is_ascii_alphanumeric()
            || byte == b'.'
            || byte == b'_' && (index == 0 || underscored || !require_alpha)
        {
            index += 1;
        } else if encoded
            && byte == b'%'
            && bytes
                .get(index + 1..index + 3)
                .is_some_and(|code| code.iter().all(u8::is_ascii_hexdigit))
        {
            index += 3;
        } else {
            return false;
        }
    }
    true
}

fn valid_chrom(value: &str, version: Option<(u8, u8)>) -> bool {
    if let Some(inner) = value
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
    {
        return !inner.is_empty()
            && !inner
                .bytes()
                .any(|byte| byte.is_ascii_whitespace() || matches!(byte, b'<' | b'>' | b','));
    }
    let _ = version;
    valid_contig_name(value)
}

fn valid_contig_name(value: &str) -> bool {
    let mut bytes = value.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };
    valid_contig_byte(first, false) && bytes.all(|byte| valid_contig_byte(byte, true))
}

fn valid_contig_byte(byte: u8, tail: bool) -> bool {
    byte.is_ascii_alphanumeric()
        || matches!(
            byte,
            b'!' | b'#'
                | b'$'
                | b'%'
                | b'&'
                | b'+'
                | b'.'
                | b'/'
                | b':'
                | b';'
                | b'?'
                | b'@'
                | b'^'
                | b'_'
                | b'|'
                | b'~'
                | b'-'
        )
        || tail && matches!(byte, b'*' | b'=')
}

fn valid_base_sequence(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| matches!(byte.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

fn valid_alternate(value: &str, version: Option<(u8, u8)>) -> bool {
    if value == "*" || valid_base_sequence(value) {
        return true;
    }
    if value.len() > 1
        && ((value.starts_with('.') && valid_base_sequence(&value[1..]))
            || (value.ends_with('.') && valid_base_sequence(&value[..value.len() - 1])))
    {
        return true;
    }
    if let Some(inner) = value
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
    {
        return !inner.is_empty()
            && !inner
                .bytes()
                .any(|byte| byte.is_ascii_whitespace() || matches!(byte, b'<' | b'>' | b','));
    }
    valid_breakend(value, version)
}

fn valid_breakend(value: &str, version: Option<(u8, u8)>) -> bool {
    let brackets: Vec<_> = value
        .char_indices()
        .filter(|(_, character)| matches!(character, '[' | ']'))
        .collect();
    if brackets.len() != 2 || brackets[0].1 != brackets[1].1 {
        return false;
    }
    let first = brackets[0].0;
    let second = brackets[1].0;
    let remote = &value[first + 1..second];
    let Some((chrom, position)) = remote.rsplit_once(':') else {
        return false;
    };
    if !valid_chrom(chrom, version) || !position.parse::<u64>().is_ok_and(|position| position > 0) {
        return false;
    }
    if first == 0 {
        valid_base_sequence(&value[second + 1..])
    } else {
        second + 1 == value.len() && valid_base_sequence(&value[..first])
    }
}

fn valid_integer(value: &str) -> bool {
    value
        .parse::<i32>()
        .is_ok_and(|value| value >= i32::MIN + 8)
}

fn split_values(value: &str) -> impl Iterator<Item = &str> {
    let mut quoted = false;
    let mut escaped = false;
    value.split(move |character| {
        if escaped {
            escaped = false;
            return false;
        }
        match character {
            '\\' if quoted => {
                escaped = true;
                false
            }
            '"' => {
                quoted = !quoted;
                false
            }
            ',' if !quoted => true,
            _ => false,
        }
    })
}

fn valid_cigar(value: &str) -> bool {
    let mut bytes = value.bytes().peekable();
    let mut operations = 0;
    while bytes.peek().is_some() {
        let mut length = 0u64;
        while let Some(byte) = bytes.peek().copied().filter(u8::is_ascii_digit) {
            bytes.next();
            length = match length
                .checked_mul(10)
                .and_then(|length| length.checked_add(u64::from(byte - b'0')))
            {
                Some(length) => length,
                None => return false,
            };
        }
        if length == 0
            || !bytes.next().is_some_and(|byte| {
                matches!(
                    byte,
                    b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X'
                )
            })
        {
            return false;
        }
        operations += 1;
    }
    operations > 0
}

fn variant_keys(position: u64, reference: &str, alternates: &str) -> Vec<VariantKey> {
    alternates
        .split(',')
        .filter(|alternate| valid_base_sequence(alternate))
        .map(|alternate| {
            let (position, reference, alternate) =
                normalize_variant(position, reference, alternate);
            VariantKey {
                position,
                reference: reference.to_owned(),
                alternate: alternate.to_owned(),
            }
        })
        .collect()
}

fn normalize_variant<'a>(
    mut position: u64,
    mut reference: &'a str,
    mut alternate: &'a str,
) -> (u64, &'a str, &'a str) {
    if !valid_base_sequence(reference) || !valid_base_sequence(alternate) {
        return (position, reference, alternate);
    }
    while reference.len() > 1
        && alternate.len() > 1
        && reference.as_bytes().last() == alternate.as_bytes().last()
    {
        reference = &reference[..reference.len() - 1];
        alternate = &alternate[..alternate.len() - 1];
    }
    while reference.len() > 1
        && alternate.len() > 1
        && reference.as_bytes().first() == alternate.as_bytes().first()
    {
        reference = &reference[1..];
        alternate = &alternate[1..];
        position += 1;
    }
    (position, reference, alternate)
}

fn parse_float(value: &str) -> Option<f64> {
    let (negative, token) = if let Some(token) = value.strip_prefix('-') {
        (true, token)
    } else {
        (false, value.strip_prefix('+').unwrap_or(value))
    };
    if token.eq_ignore_ascii_case("inf") || token.eq_ignore_ascii_case("infinity") {
        Some(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        })
    } else if token.eq_ignore_ascii_case("nan") {
        Some(f64::NAN)
    } else {
        value.parse().ok()
    }
}

fn error_chain(error: &(dyn std::error::Error + 'static)) -> String {
    let mut message = error.to_string();
    let mut source = error.source();
    while let Some(error) = source {
        let detail = error.to_string();
        if !detail.is_empty() && !message.ends_with(&detail) {
            message.push_str(": ");
            message.push_str(&detail);
        }
        source = error.source();
    }
    message
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_cardinality_uses_ploidy_and_allele_count() {
        assert_eq!(genotype_count(2, 2), Some(3));
        assert_eq!(genotype_count(3, 2), Some(6));
        assert_eq!(genotype_count(3, 3), Some(10));
    }
}
