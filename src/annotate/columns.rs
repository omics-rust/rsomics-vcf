use std::collections::HashSet;
use std::path::Path;

use noodles_vcf::Header;
use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) enum MatchField {
    Chrom,
    Pos,
    From,
    To,
    Ref,
    Alt,
    Id,
    End,
    Ignore,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MatchLayout {
    Position,
    Interval,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SourceKind {
    Tabular,
    Variant,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum SourceField {
    Id,
    Qual,
    Filter,
    Info(String),
    Format(String),
    AllInfo,
    AllFormat,
    Tabular(usize),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Destination {
    Id,
    Qual,
    Filter,
    Info(String),
    Format(String),
    AllInfo,
    AllFormat,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum WriteMode {
    Replace,
    ReplaceMissing,
    Add,
    AddMissing,
    Append,
    AppendMissing,
    ReplaceExisting,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Transfer {
    pub(crate) source: SourceField,
    pub(crate) destination: Destination,
    pub(crate) mode: WriteMode,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Column {
    Match(MatchField),
    Transfer(Transfer),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ColumnSpec {
    fields: Vec<Column>,
    transfers: Vec<Transfer>,
    match_layout: Option<MatchLayout>,
    matches_id: bool,
    matches_end: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct BoundColumns {
    spec: ColumnSpec,
    source_kind: SourceKind,
}

impl ColumnSpec {
    pub(crate) fn parse(source: &str) -> Result<Self> {
        let raw = split_checked_csv(source)?;
        let mut fields = Vec::with_capacity(raw.len());
        let mut transfers = Vec::new();
        let mut destinations = HashSet::new();

        for value in raw {
            let column = parse_column(value)?;
            if let Column::Transfer(transfer) = &column {
                let key = destination_key(&transfer.destination);
                if !destinations.insert(key) {
                    return Err(invalid(format!(
                        "annotation destination is assigned more than once: {value}"
                    )));
                }
                transfers.push(transfer.clone());
            }
            fields.push(column);
        }

        let match_layout = validate_layout(&fields)?;
        let matches_id = fields
            .iter()
            .any(|field| matches!(field, Column::Match(MatchField::Id)));
        let matches_end = fields
            .iter()
            .any(|field| matches!(field, Column::Match(MatchField::End)));

        Ok(Self {
            fields,
            transfers,
            match_layout,
            matches_id,
            matches_end,
        })
    }

    pub(crate) fn from_file(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .rs_with_context(|| format!("reading annotation columns {}", path.display()))?;
        let fields: Vec<_> = content
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .collect();
        if fields.is_empty() {
            return Err(invalid(format!(
                "annotation columns file is empty: {}",
                path.display()
            )));
        }
        Self::parse(&fields.join(","))
    }

    pub(crate) fn fields(&self) -> &[Column] {
        &self.fields
    }

    pub(crate) fn transfers(&self) -> &[Transfer] {
        &self.transfers
    }

    pub(crate) fn match_layout(&self) -> Option<MatchLayout> {
        self.match_layout
    }

    pub(crate) fn matches_id(&self) -> bool {
        self.matches_id
    }

    pub(crate) fn matches_end(&self) -> bool {
        self.matches_end
    }
}

impl BoundColumns {
    pub(crate) fn bind(
        mut spec: ColumnSpec,
        source_kind: SourceKind,
        _target: &Header,
        source: Option<&Header>,
    ) -> Result<Self> {
        match source_kind {
            SourceKind::Tabular if spec.match_layout.is_none() => {
                return Err(invalid(
                    "tabular annotation columns require CHROM with POS or FROM and TO",
                ));
            }
            SourceKind::Variant if source.is_none() => {
                return Err(invalid(
                    "variant annotation columns require a source header",
                ));
            }
            SourceKind::Tabular | SourceKind::Variant => {}
        }
        if source_kind == SourceKind::Tabular {
            let reference = spec
                .fields
                .iter()
                .any(|field| matches!(field, Column::Match(MatchField::Ref)));
            let alternate = spec
                .fields
                .iter()
                .any(|field| matches!(field, Column::Match(MatchField::Alt)));
            if reference != alternate {
                return Err(invalid(
                    "tabular allele matching requires both REF and ALT columns",
                ));
            }
            if reference && spec.match_layout == Some(MatchLayout::Interval) {
                return Err(invalid(
                    "tabular REF and ALT matching requires CHROM and POS rather than FROM and TO",
                ));
            }
        }
        if source_kind == SourceKind::Tabular {
            for (index, field) in spec.fields.iter_mut().enumerate() {
                if let Column::Transfer(transfer) = field {
                    transfer.source = SourceField::Tabular(index);
                }
            }
            spec.transfers = spec
                .fields
                .iter()
                .filter_map(|field| match field {
                    Column::Transfer(transfer) => Some(transfer.clone()),
                    Column::Match(_) => None,
                })
                .collect();
        }
        Ok(Self { spec, source_kind })
    }

    pub(crate) fn spec(&self) -> &ColumnSpec {
        &self.spec
    }

    pub(crate) fn source_kind(&self) -> SourceKind {
        self.source_kind
    }
}

fn split_checked_csv(source: &str) -> Result<Vec<&str>> {
    if source.trim().is_empty() {
        return Err(invalid("annotation columns must not be empty"));
    }
    source
        .split(',')
        .map(str::trim)
        .map(|field| {
            if field.is_empty() {
                Err(invalid("annotation columns contain an empty field"))
            } else {
                Ok(field)
            }
        })
        .collect()
}

fn parse_column(source: &str) -> Result<Column> {
    if source == "-" {
        return Ok(Column::Match(MatchField::Ignore));
    }
    if let Some(field) = parse_match(source)? {
        return Ok(Column::Match(field));
    }

    let (mode, body) = parse_mode(source)?;
    let (destination, source_field) = match body.split_once(":=") {
        Some((destination, source_field)) => {
            if destination.is_empty() || source_field.is_empty() || source_field.contains(":=") {
                return Err(invalid(format!(
                    "invalid annotation rename column: {source}"
                )));
            }
            let source_field = parse_source(source_field)?;
            (
                parse_renamed_destination(destination, &source_field)?,
                source_field,
            )
        }
        None => {
            let destination = parse_destination(body)?;
            let source_field = source_from_destination(&destination);
            (destination, source_field)
        }
    };
    validate_transfer(&destination, &source_field, source)?;
    Ok(Column::Transfer(Transfer {
        source: source_field,
        destination,
        mode,
    }))
}

fn parse_match(source: &str) -> Result<Option<MatchField>> {
    let field = match source {
        "CHROM" => Some(MatchField::Chrom),
        "POS" => Some(MatchField::Pos),
        "FROM" | "BEG" => Some(MatchField::From),
        "TO" | "END" => Some(MatchField::To),
        "REF" => Some(MatchField::Ref),
        "ALT" => Some(MatchField::Alt),
        "~ID" => Some(MatchField::Id),
        "~INFO/END" => Some(MatchField::End),
        value if value.starts_with('~') => {
            return Err(invalid(format!(
                "unsupported annotation match column: {value}"
            )));
        }
        _ => None,
    };
    Ok(field)
}

fn parse_mode(source: &str) -> Result<(WriteMode, &str)> {
    let (mode, body) = if let Some(body) = source.strip_prefix(".+") {
        (WriteMode::AddMissing, body)
    } else if let Some(body) = source.strip_prefix(".=") {
        (WriteMode::AppendMissing, body)
    } else if let Some(body) = source.strip_prefix('.') {
        (WriteMode::ReplaceMissing, body)
    } else if let Some(body) = source.strip_prefix('+') {
        (WriteMode::Add, body)
    } else if let Some(body) = source.strip_prefix('=') {
        (WriteMode::Append, body)
    } else if let Some(body) = source.strip_prefix('-') {
        (WriteMode::ReplaceExisting, body)
    } else {
        (WriteMode::Replace, source)
    };
    if body.is_empty()
        || body
            .as_bytes()
            .first()
            .is_some_and(|byte| matches!(byte, b'.' | b'+' | b'=' | b'-'))
    {
        return Err(invalid(format!("invalid annotation write mode: {source}")));
    }
    Ok((mode, body))
}

fn parse_destination(source: &str) -> Result<Destination> {
    match source {
        "ID" => Ok(Destination::Id),
        "QUAL" => Ok(Destination::Qual),
        "FILTER" => Ok(Destination::Filter),
        "INFO" => Ok(Destination::AllInfo),
        "FMT" | "FORMAT" => Ok(Destination::AllFormat),
        _ => parse_tag(source, true).map(|(format, tag)| {
            if format {
                Destination::Format(tag)
            } else {
                Destination::Info(tag)
            }
        }),
    }
}

fn parse_renamed_destination(source: &str, source_field: &SourceField) -> Result<Destination> {
    let reserved = matches!(source, "ID" | "QUAL" | "FILTER" | "INFO" | "FMT" | "FORMAT");
    if !reserved && !source.contains('/') && matches!(source_field, SourceField::Format(_)) {
        let (_, tag) = parse_tag(source, true)?;
        return Ok(Destination::Format(tag));
    }
    parse_destination(source)
}

fn parse_source(source: &str) -> Result<SourceField> {
    match source {
        "ID" => Ok(SourceField::Id),
        "QUAL" => Ok(SourceField::Qual),
        "FILTER" => Ok(SourceField::Filter),
        "INFO" => Ok(SourceField::AllInfo),
        "FMT" | "FORMAT" => Ok(SourceField::AllFormat),
        _ => parse_tag(source, false).map(|(format, tag)| {
            if format {
                SourceField::Format(tag)
            } else {
                SourceField::Info(tag)
            }
        }),
    }
}

fn parse_tag(source: &str, destination: bool) -> Result<(bool, String)> {
    let (format, tag) = if let Some(tag) = source.strip_prefix("INFO/") {
        (false, tag)
    } else if let Some(tag) = source
        .strip_prefix("FMT/")
        .or_else(|| source.strip_prefix("FORMAT/"))
    {
        (true, tag)
    } else if source.contains('/') {
        return Err(invalid(format!(
            "unsupported annotation {} field: {source}",
            if destination { "destination" } else { "source" }
        )));
    } else {
        (false, source)
    };
    if !valid_tag(tag) {
        return Err(invalid(format!(
            "invalid annotation {} tag: {source}",
            if destination { "destination" } else { "source" }
        )));
    }
    Ok((format, tag.to_owned()))
}

fn valid_tag(tag: &str) -> bool {
    !tag.is_empty()
        && tag
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'.' | b'-'))
}

fn source_from_destination(destination: &Destination) -> SourceField {
    match destination {
        Destination::Id => SourceField::Id,
        Destination::Qual => SourceField::Qual,
        Destination::Filter => SourceField::Filter,
        Destination::Info(tag) => SourceField::Info(tag.clone()),
        Destination::Format(tag) => SourceField::Format(tag.clone()),
        Destination::AllInfo => SourceField::AllInfo,
        Destination::AllFormat => SourceField::AllFormat,
    }
}

fn validate_transfer(destination: &Destination, source: &SourceField, raw: &str) -> Result<()> {
    let destination_is_format =
        matches!(destination, Destination::Format(_) | Destination::AllFormat);
    let source_is_format = matches!(source, SourceField::Format(_) | SourceField::AllFormat);
    if destination_is_format != source_is_format {
        return Err(invalid(format!(
            "annotation INFO and FORMAT fields cannot be renamed across namespaces: {raw}"
        )));
    }
    if matches!(destination, Destination::AllInfo | Destination::AllFormat)
        && !matches!(source, SourceField::AllInfo | SourceField::AllFormat)
    {
        return Err(invalid(format!(
            "whole annotation namespaces cannot be rename destinations: {raw}"
        )));
    }
    Ok(())
}

fn validate_layout(fields: &[Column]) -> Result<Option<MatchLayout>> {
    let mut chrom = 0;
    let mut pos = 0;
    let mut from = 0;
    let mut to = 0;
    let mut unique = HashSet::new();
    for field in fields {
        let Column::Match(field) = field else {
            continue;
        };
        if *field != MatchField::Ignore && !unique.insert(field.clone()) {
            return Err(invalid(format!(
                "annotation match column is assigned more than once: {field:?}"
            )));
        }
        match field {
            MatchField::Chrom => chrom += 1,
            MatchField::Pos => pos += 1,
            MatchField::From => from += 1,
            MatchField::To => to += 1,
            MatchField::Ref
            | MatchField::Alt
            | MatchField::Id
            | MatchField::End
            | MatchField::Ignore => {}
        }
    }

    let has_coordinate = chrom + pos + from + to > 0;
    if !has_coordinate {
        return Ok(None);
    }
    if chrom != 1 {
        return Err(invalid(
            "annotation coordinates require exactly one CHROM column",
        ));
    }
    match (pos, from, to) {
        (1, 0, 0) => Ok(Some(MatchLayout::Position)),
        (0, 1, 1) => Ok(Some(MatchLayout::Interval)),
        _ => Err(invalid(
            "annotation coordinates require POS or the pair FROM and TO",
        )),
    }
}

fn destination_key(destination: &Destination) -> String {
    match destination {
        Destination::Id => "ID".to_owned(),
        Destination::Qual => "QUAL".to_owned(),
        Destination::Filter => "FILTER".to_owned(),
        Destination::Info(tag) => format!("INFO/{tag}"),
        Destination::Format(tag) => format!("FORMAT/{tag}"),
        Destination::AllInfo => "INFO".to_owned(),
        Destination::AllFormat => "FORMAT".to_owned(),
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::ConfigError(message.into())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use noodles_vcf as vcf;

    use super::*;

    #[test]
    fn parses_match_and_transfer_columns() {
        let plan =
            ColumnSpec::parse("CHROM,FROM,TO,REF,ALT,+INFO/DB,FMT/NEW:=FMT/OLD,~ID").unwrap();

        assert_eq!(plan.match_layout(), Some(MatchLayout::Interval));
        assert_eq!(plan.transfers().len(), 2);
        assert!(plan.matches_id());
        assert_eq!(plan.fields().len(), 8);
    }

    #[test]
    fn parses_every_write_mode() {
        let plan =
            ColumnSpec::parse("CHROM,POS,INFO/A,.INFO/B,+INFO/C,.+INFO/D,=INFO/E,.=INFO/F,-INFO/G")
                .unwrap();
        let modes: Vec<_> = plan
            .transfers()
            .iter()
            .map(|transfer| transfer.mode)
            .collect();

        assert_eq!(
            modes,
            [
                WriteMode::Replace,
                WriteMode::ReplaceMissing,
                WriteMode::Add,
                WriteMode::AddMissing,
                WriteMode::Append,
                WriteMode::AppendMissing,
                WriteMode::ReplaceExisting,
            ]
        );
    }

    #[test]
    fn parses_variant_source_fields_without_coordinates() {
        let plan = ColumnSpec::parse("ID,QUAL,FILTER,INFO,FMT/NEW:=FMT/OLD").unwrap();

        assert_eq!(plan.match_layout(), None);
        assert_eq!(plan.transfers().len(), 5);
    }

    #[test]
    fn infers_bare_rename_destination_namespace() {
        let plan = ColumnSpec::parse("XX:=FORMAT/X-X,AA:=INFO/A-A").unwrap();

        assert!(matches!(
            &plan.transfers()[0],
            Transfer {
                source: SourceField::Format(source),
                destination: Destination::Format(destination),
                ..
            } if source == "X-X" && destination == "XX"
        ));
        assert!(matches!(
            &plan.transfers()[1],
            Transfer {
                source: SourceField::Info(source),
                destination: Destination::Info(destination),
                ..
            } if source == "A-A" && destination == "AA"
        ));
    }

    #[test]
    fn accepts_ignored_and_end_match_columns() {
        let plan = ColumnSpec::parse("CHROM,POS,-,REF,ALT,~INFO/END,INFO/DB").unwrap();

        assert!(plan.matches_end());
        assert!(matches!(
            plan.fields()[2],
            Column::Match(MatchField::Ignore)
        ));
    }

    #[test]
    fn rejects_ambiguous_coordinate_and_mode_grammar() {
        for value in [
            "POS,FROM,TO,DB",
            "CHROM,FROM,DB",
            "CHROM,POS,++DB",
            "CHROM,POS,INFO/A,INFO/A",
            "CHROM,POS,~QUAL,INFO/A",
            "CHROM,POS,FMT/A:=INFO/B",
        ] {
            assert!(ColumnSpec::parse(value).is_err(), "{value}");
        }
    }

    #[test]
    fn reads_columns_file_in_order() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("columns.txt");
        fs::write(&path, "# coordinates\nCHROM\nPOS\n\nINFO/DB\n").unwrap();

        let plan = ColumnSpec::from_file(&path).unwrap();

        assert_eq!(plan.match_layout(), Some(MatchLayout::Position));
        assert_eq!(plan.transfers().len(), 1);
    }

    #[test]
    fn binds_only_to_the_matching_source_kind() {
        let header = vcf::Header::default();
        let tabular = ColumnSpec::parse("CHROM,POS,INFO/DB").unwrap();
        assert!(BoundColumns::bind(tabular, SourceKind::Tabular, &header, None).is_ok());

        let variant = ColumnSpec::parse("INFO/DB").unwrap();
        assert!(BoundColumns::bind(variant.clone(), SourceKind::Variant, &header, None).is_err());
        let bound =
            BoundColumns::bind(variant, SourceKind::Variant, &header, Some(&header)).unwrap();
        assert_eq!(bound.source_kind(), SourceKind::Variant);
        assert_eq!(bound.spec().transfers().len(), 1);
    }

    #[test]
    fn binds_tabular_transfer_sources_to_column_positions() {
        let header = vcf::Header::default();
        let spec = ColumnSpec::parse("CHROM,FROM,TO,INFO/DB,FMT/X").unwrap();
        let bound = BoundColumns::bind(spec, SourceKind::Tabular, &header, None).unwrap();

        assert_eq!(bound.spec().transfers()[0].source, SourceField::Tabular(3));
        assert_eq!(bound.spec().transfers()[1].source, SourceField::Tabular(4));
    }

    #[test]
    fn rejects_incomplete_or_interval_tabular_allele_keys() {
        let header = vcf::Header::default();
        for raw in [
            "CHROM,POS,REF,INFO/X",
            "CHROM,POS,ALT,INFO/X",
            "CHROM,FROM,TO,REF,ALT,INFO/X",
        ] {
            let spec = ColumnSpec::parse(raw).unwrap();
            assert!(
                BoundColumns::bind(spec, SourceKind::Tabular, &header, None).is_err(),
                "{raw}"
            );
        }
    }
}
