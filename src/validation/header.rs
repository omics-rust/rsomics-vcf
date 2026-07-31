use std::collections::{HashMap, HashSet};

use super::{Diagnostics, Field, Number, Schema, Type, format_definition, info_definition};

pub(crate) fn inspect_header(raw: &[u8], diagnostics: &mut Diagnostics) -> Schema {
    let mut schema = Schema::default();
    let mut chrom_line = None;
    let mut info_lines = HashMap::new();
    let mut format_lines = HashMap::new();
    let mut filter_lines = HashMap::new();
    let mut alternate_lines = HashMap::new();
    let mut alternates = HashSet::new();
    let mut reference = false;

    let lines = raw.split_inclusive(|byte| *byte == b'\n');
    for (index, raw_line) in lines.enumerate() {
        let line_number = index + 1;
        schema.header_lines = line_number;
        let line = trim_ending(raw_line);
        let Ok(line) = std::str::from_utf8(line) else {
            diagnostics.error(
                "header.utf8",
                Some(line_number),
                None,
                "header line is not valid UTF-8",
            );
            continue;
        };

        if line_number == 1 {
            schema.version = parse_file_format(line, diagnostics);
        } else if line.starts_with("##") {
            if chrom_line.is_some() {
                diagnostics.error(
                    "header.meta-after-columns",
                    Some(line_number),
                    None,
                    "metadata must precede the #CHROM line",
                );
            }
            match line
                .strip_prefix("##")
                .and_then(|line| line.split_once('='))
            {
                Some(("INFO", value)) => {
                    if let Some((id, field)) =
                        parse_field_map(value, true, schema.version, line_number, diagnostics)
                    {
                        insert_unique(
                            &mut schema.info,
                            &mut info_lines,
                            id,
                            field,
                            "INFO",
                            line_number,
                            diagnostics,
                        );
                    }
                }
                Some(("FORMAT", value)) => {
                    if let Some((id, field)) =
                        parse_field_map(value, false, schema.version, line_number, diagnostics)
                    {
                        insert_unique(
                            &mut schema.format,
                            &mut format_lines,
                            id,
                            field,
                            "FORMAT",
                            line_number,
                            diagnostics,
                        );
                    }
                }
                Some(("FILTER", value)) => {
                    if let Some(id) =
                        parse_id_map(value, "FILTER", schema.version, line_number, diagnostics)
                    {
                        insert_name(
                            &mut schema.filters,
                            &mut filter_lines,
                            id,
                            "FILTER",
                            line_number,
                            diagnostics,
                        );
                    }
                }
                Some(("ALT", value)) => {
                    if let Some(id) =
                        parse_alternate(value, schema.version, line_number, diagnostics)
                    {
                        insert_name(
                            &mut alternates,
                            &mut alternate_lines,
                            id,
                            "ALT",
                            line_number,
                            diagnostics,
                        );
                    }
                }
                Some(("contig", value)) => {
                    if let Some((id, length)) =
                        parse_contig(value, schema.version, line_number, diagnostics)
                    {
                        schema
                            .contigs
                            .entry(id)
                            .and_modify(|current| {
                                if current.is_none() {
                                    *current = length;
                                }
                            })
                            .or_insert(length);
                    }
                }
                Some(("META", value)) => {
                    inspect_meta(value, line_number, diagnostics);
                }
                Some(("PEDIGREE", value)) => {
                    inspect_pedigree(value, schema.version, line_number, diagnostics);
                }
                Some(("SAMPLE", value)) => {
                    inspect_sample(value, schema.version, line_number, diagnostics);
                }
                Some(("assembly" | "pedigreeDB", value)) => {
                    if !valid_url(value) {
                        diagnostics.error(
                            "header.url",
                            Some(line_number),
                            None,
                            format!("metadata URL {value:?} is invalid"),
                        );
                    }
                }
                Some(("reference", value)) => {
                    if value.is_empty() {
                        diagnostics.error(
                            "header.meta-value",
                            Some(line_number),
                            Some("reference"),
                            "reference metadata must not be empty",
                        );
                    } else {
                        reference = true;
                    }
                }
                Some(("", _)) => diagnostics.error(
                    "header.meta-key",
                    Some(line_number),
                    None,
                    "metadata key must not be empty",
                ),
                Some((key, "")) => diagnostics.error(
                    "header.meta-value",
                    Some(line_number),
                    Some(key),
                    "metadata value must not be empty",
                ),
                Some(_) => {}
                None => diagnostics.error(
                    "header.meta-syntax",
                    Some(line_number),
                    None,
                    "metadata line must contain a key and value",
                ),
            }
        } else if line.starts_with("#CHROM") {
            if let Some(first) = chrom_line.replace(line_number) {
                diagnostics.error(
                    "header.duplicate-columns",
                    Some(line_number),
                    None,
                    format!("#CHROM line duplicates line {first}"),
                );
            }
            parse_columns(line, line_number, &mut schema, diagnostics);
        } else if line.starts_with('#') {
            diagnostics.error(
                "header.columns",
                Some(line_number),
                None,
                "header line is neither metadata nor the #CHROM line",
            );
        } else if !line.is_empty() {
            diagnostics.error(
                "header.record-before-columns",
                Some(line_number),
                None,
                "variant data begins before the #CHROM line",
            );
        }
    }

    if raw.is_empty() {
        diagnostics.error("header.missing", None, None, "input has no VCF header");
    } else if chrom_line.is_none() {
        diagnostics.error(
            "header.missing-columns",
            None,
            None,
            "VCF header is missing the #CHROM line",
        );
    }
    if !raw.is_empty() && !raw.ends_with(b"\n") {
        diagnostics.error(
            "header.newline",
            Some(schema.header_lines),
            None,
            "the final VCF header line is not newline-terminated",
        );
    }
    if !raw.is_empty() && !reference {
        diagnostics.warning(
            "header.reference",
            None,
            None,
            "header has no reference metadata",
        );
    }

    schema
}

fn parse_file_format(line: &str, diagnostics: &mut Diagnostics) -> Option<(u8, u8)> {
    let Some(version) = line.strip_prefix("##fileformat=VCFv") else {
        diagnostics.error(
            "header.fileformat",
            Some(1),
            None,
            "first line must be ##fileformat=VCFv4.x",
        );
        return None;
    };
    let Some((major, minor)) = version.split_once('.') else {
        diagnostics.error(
            "header.fileformat",
            Some(1),
            None,
            "fileformat version must contain major and minor numbers",
        );
        return None;
    };
    let parsed = major.parse::<u8>().ok().zip(minor.parse::<u8>().ok());
    match parsed {
        Some((4, minor @ 1..=5)) => Some((4, minor)),
        _ => {
            diagnostics.error(
                "header.version",
                Some(1),
                None,
                format!("unsupported VCF version {version}; expected VCFv4.1 through VCFv4.5"),
            );
            None
        }
    }
}

fn parse_columns(
    line: &str,
    line_number: usize,
    schema: &mut Schema,
    diagnostics: &mut Diagnostics,
) {
    let columns: Vec<_> = line.split('\t').collect();
    const FIXED: [&str; 8] = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ];
    if columns.len() < FIXED.len() || columns[..FIXED.len()] != FIXED {
        diagnostics.error(
            "header.fixed-columns",
            Some(line_number),
            None,
            "#CHROM line must begin with CHROM, POS, ID, REF, ALT, QUAL, FILTER, and INFO",
        );
        return;
    }
    match columns.len() {
        8 => {}
        9 => diagnostics.error(
            "header.samples",
            Some(line_number),
            Some("FORMAT"),
            "FORMAT requires at least one sample column",
        ),
        _ if columns[8] != "FORMAT" => diagnostics.error(
            "header.format-column",
            Some(line_number),
            Some("FORMAT"),
            "sample columns must be preceded by FORMAT",
        ),
        _ => {
            let mut names = HashSet::new();
            for sample in &columns[9..] {
                if sample.is_empty() {
                    diagnostics.error(
                        "header.sample-name",
                        Some(line_number),
                        Some("sample"),
                        "sample name must not be empty",
                    );
                } else if !names.insert(*sample) {
                    diagnostics.error(
                        "header.duplicate-sample",
                        Some(line_number),
                        Some("sample"),
                        format!("sample name {sample} is duplicated"),
                    );
                } else {
                    schema.samples.push((*sample).to_owned());
                }
            }
        }
    }
}

fn parse_field_map(
    value: &str,
    info: bool,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<(String, Field)> {
    let fields = parse_map(
        value,
        if info { "INFO" } else { "FORMAT" },
        line_number,
        diagnostics,
    )?;
    validate_field_order(
        value,
        &["ID", "Number", "Type", "Description"],
        if info { "INFO" } else { "FORMAT" },
        line_number,
        diagnostics,
    )?;
    let id = required(&fields, "ID", line_number, diagnostics)?;
    validate_id(
        id,
        if info { "INFO" } else { "FORMAT" },
        line_number,
        diagnostics,
    )?;
    let number = required(&fields, "Number", line_number, diagnostics)
        .and_then(|value| parse_number(value, info, version, line_number, diagnostics))?;
    let ty = required(&fields, "Type", line_number, diagnostics)
        .and_then(|value| parse_type(value, info, line_number, diagnostics))?;
    required(&fields, "Description", line_number, diagnostics)?;

    if info && ty != Type::Flag && number == Number::Count(0) {
        diagnostics.error(
            "header.number-zero",
            Some(line_number),
            Some("Number"),
            "only INFO fields with Type=Flag may use Number=0",
        );
        return None;
    }

    let field = Field { number, ty };
    let expected = if info {
        info_definition(id, version)
    } else {
        format_definition(id, version)
    };
    if let Some(expected) = expected
        && field != expected
    {
        diagnostics.error(
            "header.reserved-definition",
            Some(line_number),
            Some(id),
            format!(
                "reserved {} field {id} requires Number={:?} and Type={:?}",
                if info { "INFO" } else { "FORMAT" },
                expected.number,
                expected.ty
            ),
        );
        return None;
    }

    Some((id.to_owned(), field))
}

fn parse_id_map(
    value: &str,
    kind: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<String> {
    let fields = parse_map(value, kind, line_number, diagnostics)?;
    let keys = map_keys(value)?;
    match keys.as_slice() {
        ["ID"] if version.is_some_and(|version| version >= (4, 4)) => {}
        ["ID", "Description"] => {
            require_quoted_description(value, kind, line_number, diagnostics)?;
        }
        _ => {
            diagnostics.error(
                "header.field-order",
                Some(line_number),
                Some(kind),
                format!("{kind} metadata fields are not in specification order"),
            );
            return None;
        }
    }
    let id = required(&fields, "ID", line_number, diagnostics)?;
    validate_id(id, kind, line_number, diagnostics)?;
    Some(id.to_owned())
}

fn parse_alternate(
    value: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<String> {
    let fields = parse_map(value, "ALT", line_number, diagnostics)?;
    let keys = map_keys(value)?;
    let legacy = match keys.as_slice() {
        ["ID"] if version.is_some_and(|version| version >= (4, 4)) => false,
        ["ID", "Description"] => false,
        ["ID", "Number", "Type", "Description"]
            if version.is_some_and(|version| version <= (4, 2))
                && fields.get("Number") == Some(&"1")
                && fields.get("Type") == Some(&"String") =>
        {
            true
        }
        _ => {
            diagnostics.error(
                "header.alt-fields",
                Some(line_number),
                Some("ALT"),
                "ALT metadata requires ID and quoted Description in specification order",
            );
            return None;
        }
    };
    if fields.contains_key("Description") {
        require_quoted_description(value, "ALT", line_number, diagnostics)?;
    }
    let id = required(&fields, "ID", line_number, diagnostics)?;
    validate_alt_id(id, legacy, line_number, diagnostics)?;
    Some(id.to_owned())
}

fn parse_contig(
    value: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<(String, Option<u64>)> {
    let fields = parse_map(value, "contig", line_number, diagnostics)?;
    let id = required(&fields, "ID", line_number, diagnostics)?;
    if !valid_contig_id(id, version) {
        diagnostics.error(
            "header.id",
            Some(line_number),
            Some("ID"),
            format!("contig ID {id:?} is invalid"),
        );
        return None;
    }
    let length = fields.get("length").and_then(|value| {
        value
            .parse::<u64>()
            .ok()
            .filter(|value| *value > 0)
            .or_else(|| {
                diagnostics.error(
                    "header.contig-length",
                    Some(line_number),
                    Some("length"),
                    "contig length must be a positive integer",
                );
                None
            })
    });
    Some((id.to_owned(), length))
}

fn inspect_meta(value: &str, line_number: usize, diagnostics: &mut Diagnostics) {
    let Some(fields) = parse_map(value, "META", line_number, diagnostics) else {
        return;
    };
    if map_keys(value).as_deref() != Some(&["ID", "Number", "Type", "Values"]) {
        diagnostics.error(
            "header.meta-fields",
            Some(line_number),
            Some("META"),
            "META fields must be ID, Number, Type, and Values in specification order",
        );
        return;
    }
    let Some(id) = required(&fields, "ID", line_number, diagnostics) else {
        return;
    };
    if validate_id(id, "META", line_number, diagnostics).is_none() {
        return;
    }
    if fields.get("Number") != Some(&".") {
        diagnostics.error(
            "header.meta-number",
            Some(line_number),
            Some("Number"),
            "META Number must be .",
        );
    }
    if fields.get("Type") != Some(&"String") {
        diagnostics.error(
            "header.meta-type",
            Some(line_number),
            Some("Type"),
            "META Type must be String",
        );
    }
    if !fields.get("Values").is_some_and(|value| {
        value
            .strip_prefix('[')
            .and_then(|value| value.strip_suffix(']'))
            .is_some_and(|value| {
                !value.is_empty() && value.split(',').all(|value| !value.trim().is_empty())
            })
    }) {
        diagnostics.error(
            "header.meta-values",
            Some(line_number),
            Some("Values"),
            "META Values must be a nonempty square-bracketed list",
        );
    }
}

fn inspect_pedigree(
    value: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) {
    let Some(fields) = parse_map(value, "PEDIGREE", line_number, diagnostics) else {
        return;
    };
    if version.is_some_and(|version| version >= (4, 3))
        && required(&fields, "ID", line_number, diagnostics).is_none()
    {
        return;
    }
    for (key, value) in fields {
        if !valid_genome_id(value) {
            diagnostics.error(
                "header.pedigree-value",
                Some(line_number),
                Some(key),
                format!("PEDIGREE value {value:?} is invalid"),
            );
        }
    }
}

fn inspect_sample(
    value: &str,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) {
    let Some(sample) = parse_map(value, "SAMPLE", line_number, diagnostics) else {
        return;
    };
    let Some(id) = required(&sample, "ID", line_number, diagnostics) else {
        return;
    };
    if !valid_genome_id(id) {
        diagnostics.error(
            "header.sample-id",
            Some(line_number),
            Some("ID"),
            format!("SAMPLE ID {id:?} is invalid"),
        );
    }
    if version.is_some_and(|version| version <= (4, 2)) {
        for key in ["Genomes", "Mixture", "Description"] {
            required(&sample, key, line_number, diagnostics);
        }
        let Some(body) = value
            .strip_prefix('<')
            .and_then(|value| value.strip_suffix('>'))
        else {
            return;
        };
        for field in fields(body) {
            if let Some((key, raw)) = field.split_once('=')
                && key != "Description"
                && raw.starts_with('"')
            {
                diagnostics.error(
                    "header.sample-value",
                    Some(line_number),
                    Some(key),
                    "only SAMPLE Description may be quoted in VCF 4.2",
                );
            }
        }
    }
}

fn valid_genome_id(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
}

fn valid_url(value: &str) -> bool {
    let value = value
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
        .unwrap_or(value);
    let Some((scheme, rest)) = value.split_once("://") else {
        return false;
    };
    if !matches!(scheme, "ftp" | "http" | "https") {
        return false;
    }
    let authority = rest.split('/').next().unwrap_or_default();
    let host_port = authority.rsplit('@').next().unwrap_or_default();
    let host = host_port
        .rsplit_once(':')
        .map_or(host_port, |(host, _)| host);
    !host.is_empty()
        && (host.contains('.') || host.bytes().any(|byte| byte.is_ascii_alphabetic()))
        && host
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'-'))
        && !host.split('.').any(str::is_empty)
}

fn parse_map<'a>(
    value: &'a str,
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<HashMap<&'a str, &'a str>> {
    let Some(body) = value
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
    else {
        diagnostics.error(
            "header.map-syntax",
            Some(line_number),
            Some(kind),
            format!("{kind} metadata must use angle-bracket fields"),
        );
        return None;
    };
    let mut output = HashMap::new();
    for field in fields(body) {
        let Some((key, value)) = field.split_once('=') else {
            diagnostics.error(
                "header.map-field",
                Some(line_number),
                Some(kind),
                format!("{kind} metadata field {field:?} has no value"),
            );
            return None;
        };
        if key.is_empty() || output.insert(key, unquote(value)).is_some() {
            diagnostics.error(
                "header.map-field",
                Some(line_number),
                Some(kind),
                format!("{kind} metadata has an empty or duplicate field {key:?}"),
            );
            return None;
        }
    }
    Some(output)
}

fn validate_field_order(
    value: &str,
    expected: &[&str],
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<()> {
    let keys = map_keys(value)?;
    if keys.len() < expected.len() || keys[..expected.len()] != expected[..] {
        diagnostics.error(
            "header.field-order",
            Some(line_number),
            Some(kind),
            format!("{kind} metadata fields are not in specification order"),
        );
        return None;
    }
    require_quoted_description(value, kind, line_number, diagnostics)
}

fn require_quoted_description(
    value: &str,
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<()> {
    let body = value.strip_prefix('<')?.strip_suffix('>')?;
    let raw = fields(body)
        .filter_map(|field| field.split_once('='))
        .find(|(key, _)| *key == "Description")
        .map(|(_, value)| value)?;
    if !valid_quoted_value(raw) {
        diagnostics.error(
            "header.description",
            Some(line_number),
            Some(kind),
            format!("{kind} Description must be enclosed in double quotes"),
        );
        None
    } else {
        Some(())
    }
}

fn valid_quoted_value(value: &str) -> bool {
    let Some(inner) = value
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
    else {
        return false;
    };
    let mut escaped = false;
    for character in inner.chars() {
        if escaped {
            escaped = false;
        } else if character == '\\' {
            escaped = true;
        } else if character == '"' {
            return false;
        }
    }
    !escaped
}

fn map_keys(value: &str) -> Option<Vec<&str>> {
    let body = value.strip_prefix('<')?.strip_suffix('>')?;
    fields(body)
        .map(|field| field.split_once('=').map(|(key, _)| key))
        .collect()
}

fn fields(body: &str) -> impl Iterator<Item = &str> {
    let mut quoted = false;
    let mut escaped = false;
    let mut value_start = false;
    let mut bracketed = false;
    body.split(move |character| {
        if escaped {
            escaped = false;
            return false;
        }
        match character {
            '\\' if quoted => {
                escaped = true;
                false
            }
            '"' if quoted => {
                quoted = false;
                false
            }
            '"' if value_start => {
                quoted = true;
                value_start = false;
                false
            }
            '=' if !quoted => {
                value_start = true;
                false
            }
            '[' if value_start => {
                bracketed = true;
                value_start = false;
                false
            }
            ']' if bracketed => {
                bracketed = false;
                false
            }
            ',' if !quoted && !bracketed => {
                value_start = false;
                true
            }
            _ => {
                value_start = false;
                false
            }
        }
    })
}

fn unquote(value: &str) -> &str {
    value
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
        .unwrap_or(value)
}

fn required<'a>(
    fields: &'a HashMap<&str, &'a str>,
    key: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<&'a str> {
    fields.get(key).copied().or_else(|| {
        diagnostics.error(
            "header.required-field",
            Some(line_number),
            Some(key),
            format!("structured metadata is missing {key}"),
        );
        None
    })
}

fn parse_number(
    value: &str,
    info: bool,
    version: Option<(u8, u8)>,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<Number> {
    let number = match value {
        "A" => Number::Alternate,
        "R" => Number::Alleles,
        "G" => Number::Genotypes,
        "LA" if !info && version.is_some_and(|version| version >= (4, 5)) => Number::LocalAlternate,
        "LR" if !info && version.is_some_and(|version| version >= (4, 5)) => Number::LocalAlleles,
        "LG" if !info && version.is_some_and(|version| version >= (4, 5)) => Number::LocalGenotypes,
        "P" if !info && version.is_some_and(|version| version >= (4, 4)) => Number::Ploidy,
        "M" if !info && version.is_some_and(|version| version >= (4, 5)) => {
            Number::BaseModifications
        }
        "." => Number::Unknown,
        _ => match value.parse::<usize>() {
            Ok(value) => Number::Count(value),
            Err(_) => {
                diagnostics.error(
                    "header.number",
                    Some(line_number),
                    Some("Number"),
                    format!("invalid Number value {value:?}"),
                );
                return None;
            }
        },
    };
    Some(number)
}

fn parse_type(
    value: &str,
    info: bool,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<Type> {
    let ty = match value {
        "Integer" => Type::Integer,
        "Float" => Type::Float,
        "Flag" if info => Type::Flag,
        "Character" => Type::Character,
        "String" => Type::String,
        _ => {
            diagnostics.error(
                "header.type",
                Some(line_number),
                Some("Type"),
                format!("invalid Type value {value:?}"),
            );
            return None;
        }
    };
    Some(ty)
}

fn validate_id(
    id: &str,
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<()> {
    let valid = match kind {
        "INFO" | "FORMAT" => {
            (kind == "INFO" && id == "1000G")
                || id
                    .bytes()
                    .next()
                    .is_some_and(|byte| byte.is_ascii_alphabetic() || byte == b'_')
                    && id
                        .bytes()
                        .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.'))
        }
        _ => {
            !id.is_empty()
                && !id
                    .bytes()
                    .any(|byte| byte.is_ascii_whitespace() || matches!(byte, b',' | b'<' | b'>'))
        }
    };
    if !valid {
        diagnostics.error(
            "header.id",
            Some(line_number),
            Some("ID"),
            format!("{kind} ID {id:?} is invalid"),
        );
        None
    } else {
        Some(())
    }
}

fn validate_alt_id(
    id: &str,
    legacy: bool,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) -> Option<()> {
    let supported_prefix = id
        .split(':')
        .next()
        .is_some_and(|prefix| matches!(prefix, "DEL" | "INS" | "DUP" | "INV" | "CNV" | "BND"));
    let valid = (!legacy && !id.contains(':') || supported_prefix)
        && !id.is_empty()
        && !id
            .bytes()
            .any(|byte| byte.is_ascii_whitespace() || matches!(byte, b',' | b'<' | b'>'));
    if valid {
        Some(())
    } else {
        diagnostics.error(
            "header.alt-id",
            Some(line_number),
            Some("ALT"),
            format!("ALT ID {id:?} is invalid"),
        );
        None
    }
}

fn valid_contig_id(id: &str, _version: Option<(u8, u8)>) -> bool {
    let mut bytes = id.bytes();
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

fn insert_unique(
    target: &mut HashMap<String, Field>,
    lines: &mut HashMap<String, usize>,
    id: String,
    field: Field,
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) {
    if let Some(first) = lines.insert(id.clone(), line_number) {
        diagnostics.error(
            "header.duplicate-definition",
            Some(line_number),
            Some(kind),
            format!("{kind} {id} duplicates the definition on line {first}"),
        );
    } else {
        target.insert(id, field);
    }
}

fn insert_name(
    target: &mut HashSet<String>,
    lines: &mut HashMap<String, usize>,
    id: String,
    kind: &str,
    line_number: usize,
    diagnostics: &mut Diagnostics,
) {
    if let Some(first) = lines.insert(id.clone(), line_number) {
        diagnostics.error(
            "header.duplicate-definition",
            Some(line_number),
            Some(kind),
            format!("{kind} {id} duplicates the definition on line {first}"),
        );
    } else {
        target.insert(id);
    }
}

fn trim_ending(mut line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\n') {
        line = &line[..line.len() - 1];
    }
    if line.last() == Some(&b'\r') {
        line = &line[..line.len() - 1];
    }
    line
}
