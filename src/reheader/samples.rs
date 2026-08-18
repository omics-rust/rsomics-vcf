use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;

use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum SampleSource {
    List(String),
    File(PathBuf),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) enum SampleEdit {
    Positional(Vec<String>),
    Pairs(Vec<(String, String)>),
}

impl SampleEdit {
    pub(super) fn read(source: &SampleSource) -> Result<Self> {
        match source {
            SampleSource::List(source) => Self::parse_list(source),
            SampleSource::File(path) => {
                let source = fs::read(path)
                    .rs_with_context(|| format!("reading sample names from {}", path.display()))?;
                Self::parse_file(&source)
            }
        }
    }

    fn parse_list(source: &str) -> Result<Self> {
        let names: Vec<_> = source.split(',').map(str::to_owned).collect();
        for name in &names {
            validate_name(name)?;
        }
        Ok(Self::Positional(names))
    }

    fn parse_file(source: &[u8]) -> Result<Self> {
        let source = std::str::from_utf8(source)
            .map_err(|error| invalid(format!("sample file is not valid UTF-8: {error}")))?;
        let mut rows = Vec::new();
        for (index, line) in source.lines().enumerate() {
            if line.trim_ascii().is_empty() {
                continue;
            }
            let fields = fields(line)
                .map_err(|error| invalid(format!("sample file line {}: {error}", index + 1)))?;
            if !(1..=2).contains(&fields.len()) {
                return Err(invalid(format!(
                    "sample file line {} must contain one or two fields",
                    index + 1
                )));
            }
            for name in &fields {
                validate_name(name)
                    .map_err(|error| invalid(format!("sample file line {}: {error}", index + 1)))?;
            }
            rows.push(fields);
        }
        let width = rows
            .first()
            .map(Vec::len)
            .ok_or_else(|| invalid("sample file contains no names"))?;
        if rows.iter().any(|row| row.len() != width) {
            return Err(invalid(
                "sample file cannot mix positional names and old-to-new pairs",
            ));
        }
        match width {
            1 => Ok(Self::Positional(
                rows.into_iter().map(|mut row| row.remove(0)).collect(),
            )),
            2 => Ok(Self::Pairs(
                rows.into_iter()
                    .map(|mut row| (row.remove(0), row.remove(0)))
                    .collect(),
            )),
            _ => unreachable!("sample rows have one or two fields"),
        }
    }

    pub(super) fn apply(&self, current: &[String]) -> Result<Vec<String>> {
        if current.is_empty() {
            return Err(invalid("cannot rename samples in a sites-only VCF header"));
        }
        let names = match self {
            Self::Positional(names) => {
                if names.len() != current.len() {
                    return Err(invalid(format!(
                        "sample edit must contain {} names, found {}",
                        current.len(),
                        names.len()
                    )));
                }
                names.clone()
            }
            Self::Pairs(pairs) => {
                let mut names = current.to_vec();
                let mut sources = HashSet::with_capacity(pairs.len());
                for (source, target) in pairs {
                    if !sources.insert(source) {
                        return Err(invalid(format!(
                            "duplicate sample mapping source: {source}"
                        )));
                    }
                    let index =
                        current
                            .iter()
                            .position(|name| name == source)
                            .ok_or_else(|| {
                                invalid(format!("unknown sample mapping source: {source}"))
                            })?;
                    names[index] = target.clone();
                }
                names
            }
        };
        let mut unique = HashSet::with_capacity(names.len());
        for name in &names {
            validate_name(name)?;
            if !unique.insert(name) {
                return Err(invalid(format!("duplicate final sample name: {name}")));
            }
        }
        Ok(names)
    }
}

fn fields(line: &str) -> Result<Vec<String>> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut escaped = false;
    let mut active = false;
    for character in line.chars() {
        if escaped {
            field.push(character);
            active = true;
            escaped = false;
        } else if character == '\\' {
            escaped = true;
            active = true;
        } else if character.is_ascii_whitespace() {
            if active {
                fields.push(std::mem::take(&mut field));
                active = false;
            }
        } else {
            field.push(character);
            active = true;
        }
    }
    if escaped {
        return Err(invalid("sample name ends with an incomplete escape"));
    }
    if active {
        fields.push(field);
    }
    Ok(fields)
}

fn validate_name(name: &str) -> Result<()> {
    if name.is_empty() {
        return Err(invalid("sample names cannot be empty"));
    }
    if name
        .chars()
        .any(|character| character == '\t' || character == '\r' || character == '\n')
    {
        return Err(invalid(
            "sample names cannot contain tab, carriage return, or line feed",
        ));
    }
    if let Some(character) = name.chars().find(|character| character.is_control()) {
        return Err(invalid(format!(
            "sample name contains control character U+{:04X}",
            u32::from(character)
        )));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[cfg(test)]
mod tests {
    use super::SampleEdit;

    fn current() -> Vec<String> {
        ["S1", "S2"].map(str::to_owned).to_vec()
    }

    #[test]
    fn list_and_one_column_file_replace_positionally() {
        let list = SampleEdit::parse_list("Tumor,Normal").unwrap();
        assert_eq!(list.apply(&current()).unwrap(), ["Tumor", "Normal"]);

        let file = SampleEdit::parse_file(b"Tumor\n\nNormal\n").unwrap();
        assert_eq!(file.apply(&current()).unwrap(), ["Tumor", "Normal"]);
    }

    #[test]
    fn escaped_pairs_rename_known_samples() {
        let edit = SampleEdit::parse_file(b"S1\tTumor\\ One\nS2 Normal\\ Two\n").unwrap();
        assert_eq!(edit.apply(&current()).unwrap(), ["Tumor One", "Normal Two"]);
    }

    #[test]
    fn a_partial_pair_map_retains_unmentioned_samples() {
        let edit = SampleEdit::parse_file(b"S1\tTumor\n").unwrap();
        assert_eq!(edit.apply(&current()).unwrap(), ["Tumor", "S2"]);
    }

    #[test]
    fn mappings_reject_unknown_duplicate_and_conflicting_names() {
        for source in [
            b"missing\tN1\n".as_slice(),
            b"S1\tN1\nS1\tN2\n",
            b"S1\tN\nS2\tN\n",
            b"S1\tS2\n",
            b"S1\tN1\textra\nS2\tN2\n",
            b"S1\tN1\nN2\n",
        ] {
            assert!(
                SampleEdit::parse_file(source)
                    .and_then(|edit| edit.apply(&current()))
                    .is_err(),
                "{source:?}"
            );
        }
    }

    #[test]
    fn positional_edits_require_exactly_unique_nonempty_names() {
        for list in ["", "N1", "N1,N1", "N1,", ",N2", "N1,\t"] {
            assert!(
                SampleEdit::parse_list(list)
                    .and_then(|edit| edit.apply(&current()))
                    .is_err(),
                "{list:?}"
            );
        }
    }

    #[test]
    fn file_grammar_rejects_invalid_encoding_and_escapes() {
        for source in [
            b"S1\\".as_slice(),
            b"\n\n",
            b"S1\tN1\textra\n",
            b"S1\tN\\\t1\n",
            &[0xff, b'\n'],
        ] {
            assert!(SampleEdit::parse_file(source).is_err(), "{source:?}");
        }
    }

    #[test]
    fn sites_only_headers_reject_sample_edits() {
        for edit in [
            SampleEdit::parse_list("S1").unwrap(),
            SampleEdit::parse_file(b"S1\n").unwrap(),
            SampleEdit::parse_file(b"old\tnew\n").unwrap(),
        ] {
            assert!(edit.apply(&[]).is_err());
        }
    }
}
