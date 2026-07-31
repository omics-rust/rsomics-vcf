mod definitions;
mod header;
mod record;
mod v44;

pub(crate) use definitions::{format_definition, info_definition};
pub(crate) use header::inspect_header;
pub(crate) use record::{RecordPosition, inspect_record};
pub(crate) use v44::inspect as inspect_v44_record;

use serde::Serialize;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Diagnostic {
    pub severity: Severity,
    pub code: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub line: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub field: Option<String>,
    pub message: String,
}

pub(crate) struct Diagnostics {
    limit: usize,
    errors: usize,
    warnings: usize,
    items: Vec<Diagnostic>,
}

impl Diagnostics {
    pub(crate) fn new(limit: usize) -> Self {
        Self {
            limit,
            errors: 0,
            warnings: 0,
            items: Vec::with_capacity(limit.min(128)),
        }
    }

    pub(crate) fn error(
        &mut self,
        code: &str,
        line: Option<usize>,
        field: Option<&str>,
        message: impl Into<String>,
    ) {
        self.push(Severity::Error, code, line, field, message);
    }

    pub(crate) fn warning(
        &mut self,
        code: &str,
        line: Option<usize>,
        field: Option<&str>,
        message: impl Into<String>,
    ) {
        self.push(Severity::Warning, code, line, field, message);
    }

    fn push(
        &mut self,
        severity: Severity,
        code: &str,
        line: Option<usize>,
        field: Option<&str>,
        message: impl Into<String>,
    ) {
        match severity {
            Severity::Error => self.errors += 1,
            Severity::Warning => self.warnings += 1,
        }
        if self.items.len() < self.limit {
            self.items.push(Diagnostic {
                severity,
                code: code.to_owned(),
                line,
                field: field.map(str::to_owned),
                message: message.into(),
            });
        }
    }

    pub(crate) fn errors(&self) -> usize {
        self.errors
    }

    pub(crate) fn warnings(&self) -> usize {
        self.warnings
    }

    pub(crate) fn into_items(self) -> (Vec<Diagnostic>, bool) {
        let truncated = self.errors + self.warnings > self.items.len();
        (self.items, truncated)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Number {
    Count(usize),
    Alternate,
    Alleles,
    Genotypes,
    LocalAlternate,
    LocalAlleles,
    LocalGenotypes,
    Ploidy,
    BaseModifications,
    Unknown,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Type {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Field {
    pub(crate) number: Number,
    pub(crate) ty: Type,
}

impl Field {
    pub(crate) const fn new(number: Number, ty: Type) -> Self {
        Self { number, ty }
    }
}

#[derive(Default)]
pub(crate) struct Schema {
    pub(crate) version: Option<(u8, u8)>,
    pub(crate) header_lines: usize,
    pub(crate) samples: Vec<String>,
    pub(crate) info: std::collections::HashMap<String, Field>,
    pub(crate) format: std::collections::HashMap<String, Field>,
    pub(crate) filters: std::collections::HashSet<String>,
    pub(crate) contigs: std::collections::HashMap<String, Option<u64>>,
}
