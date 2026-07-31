use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_common::Result;
use serde::Serialize;

use crate::format::{Reader, trim_line_ending};
use crate::validation::{Diagnostics, RecordPosition, inspect_header, inspect_record};

pub use crate::validation::{Diagnostic, Severity};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum InputFormat {
    Vcf,
    Bcf,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options {
    pub max_diagnostics: usize,
    pub require_evidence: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_diagnostics: 100,
            require_evidence: false,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Report {
    pub format: InputFormat,
    pub version: Option<String>,
    pub records: usize,
    pub errors: usize,
    pub warnings: usize,
    pub diagnostics: Vec<Diagnostic>,
    pub truncated: bool,
}

impl Report {
    pub fn is_valid(&self) -> bool {
        self.errors == 0
    }
}

pub fn check(input: &Path, options: Options) -> Result<Report> {
    let mut reader = Reader::open(input)?;
    let format = if reader.is_text() {
        InputFormat::Vcf
    } else {
        InputFormat::Bcf
    };
    let mut diagnostics = Diagnostics::new(options.max_diagnostics);
    let raw_header = reader.read_raw_header()?;
    let schema = inspect_header(&raw_header, &mut diagnostics);
    let header = if reader.is_text() {
        None
    } else {
        match reader.parse_header(&raw_header) {
            Ok(header) => Some(header),
            Err(error) => {
                diagnostics.error("header.parse", None, None, error.to_string());
                None
            }
        }
    };

    let records = if reader.is_text() {
        inspect_text_records(
            &mut reader,
            &schema,
            options.require_evidence,
            &mut diagnostics,
        )
    } else {
        header.as_ref().map_or(0, |header| {
            inspect_bcf_records(
                &mut reader,
                header,
                &schema,
                options.require_evidence,
                &mut diagnostics,
            )
        })
    };

    let errors = diagnostics.errors();
    let warnings = diagnostics.warnings();
    let (diagnostics, truncated) = diagnostics.into_items();
    Ok(Report {
        format,
        version: schema
            .version
            .map(|(major, minor)| format!("{major}.{minor}")),
        records,
        errors,
        warnings,
        diagnostics,
        truncated,
    })
}

fn inspect_text_records(
    reader: &mut Reader,
    schema: &crate::validation::Schema,
    require_evidence: bool,
    diagnostics: &mut Diagnostics,
) -> usize {
    let mut raw = Vec::with_capacity(4096);
    let mut records = 0;
    let mut lines = 0;
    let mut order = RecordOrder::default();
    loop {
        let number = lines + 1;
        let (read, terminated) = match reader.read_text_record_with_termination(&mut raw, number) {
            Ok(result) => result,
            Err(error) => {
                diagnostics.error(
                    "input.read",
                    Some(schema.header_lines + number),
                    None,
                    error.to_string(),
                );
                break;
            }
        };
        if read == 0 {
            break;
        }
        lines += 1;
        let line = schema.header_lines + lines;
        if !terminated && schema.version.is_none_or(|version| version < (4, 5)) {
            diagnostics.error(
                "input.newline",
                Some(line),
                None,
                "the final VCF record is not newline-terminated",
            );
        }
        if raw.is_empty() && schema.version.is_some_and(|version| version >= (4, 4)) {
            continue;
        }
        let position = inspect_record(&raw, line, schema, None, require_evidence, diagnostics);
        order.inspect(position, line, diagnostics);
        records += 1;
    }
    records
}

fn inspect_bcf_records(
    reader: &mut Reader,
    header: &vcf::Header,
    schema: &crate::validation::Schema,
    require_evidence: bool,
    diagnostics: &mut Diagnostics,
) -> usize {
    let mut record = bcf::Record::default();
    let mut raw = Vec::with_capacity(4096);
    let mut records = 0;
    let mut order = RecordOrder::default();
    loop {
        let number = records + 1;
        let read = match reader.read_bcf_record(&mut record, number) {
            Ok(read) => read,
            Err(error) => {
                diagnostics.error(
                    "bcf.record",
                    Some(schema.header_lines + number),
                    None,
                    error.to_string(),
                );
                break;
            }
        };
        if read == 0 {
            break;
        }

        raw.clear();
        match vcf::io::Writer::new(&mut raw).write_variant_record(header, &record) {
            Ok(()) => {
                trim_line_ending(&mut raw);
                let position = inspect_record(
                    &raw,
                    schema.header_lines + number,
                    schema,
                    Some(header),
                    require_evidence,
                    diagnostics,
                );
                order.inspect(position, schema.header_lines + number, diagnostics);
            }
            Err(error) => diagnostics.error(
                "bcf.decode",
                Some(schema.header_lines + number),
                None,
                error_chain(&error),
            ),
        }
        records += 1;
    }
    records
}

#[derive(Default)]
struct RecordOrder {
    current_chrom: Option<String>,
    previous_position: u64,
    finished_chroms: HashSet<String>,
    variants: BTreeMap<u64, HashMap<(String, String), usize>>,
}

impl RecordOrder {
    fn inspect(
        &mut self,
        record: Option<RecordPosition<'_>>,
        line: usize,
        diagnostics: &mut Diagnostics,
    ) {
        let Some(record) = record else {
            return;
        };
        let same_chrom = self.current_chrom.as_deref() == Some(record.chrom);
        if same_chrom {
            if record.position < self.previous_position {
                diagnostics.error(
                    "record.position-order",
                    Some(line),
                    Some("POS"),
                    format!(
                        "position {} follows position {} on contig {}",
                        record.position, self.previous_position, record.chrom
                    ),
                );
            } else {
                while self
                    .variants
                    .first_key_value()
                    .is_some_and(|(position, _)| *position < record.position)
                {
                    self.variants.pop_first();
                }
            }
        } else {
            if let Some(previous) = self.current_chrom.replace(record.chrom.to_owned()) {
                self.finished_chroms.insert(previous);
            }
            if self.finished_chroms.contains(record.chrom) {
                diagnostics.error(
                    "record.contig-order",
                    Some(line),
                    Some("CHROM"),
                    format!("contig {} is not a contiguous block", record.chrom),
                );
            }
            self.variants.clear();
        }
        for variant in record.variants {
            let alleles = (variant.reference, variant.alternate);
            if let Some(first) = self
                .variants
                .entry(variant.position)
                .or_default()
                .insert(alleles, line)
            {
                diagnostics.error(
                    "record.duplicate-variant",
                    Some(line),
                    None,
                    format!("variant duplicates the normalized allele on line {first}"),
                );
            }
        }
        self.previous_position = record.position;
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
    use std::io::Write;

    use super::*;

    fn check_text(text: &str, options: Options) -> Report {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        check(file.path(), options).unwrap()
    }

    #[test]
    fn valid_typed_record_passes() {
        let report = check_text(
            "##fileformat=VCFv4.5\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele counts\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
1\t1\t.\tA\tC,G\t30\tPASS\tAC=1,1\tGT:PL\t1|2:4,3,2,1,0,5\n",
            Options::default(),
        );
        assert!(report.is_valid(), "{:?}", report.diagnostics);
        assert_eq!(report.records, 1);
    }

    #[test]
    fn invalid_cardinality_and_genotype_are_reported() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele counts\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
1\t1\t.\tA\tC,G\t30\tPASS\tAC=1\tGT:PL\t3/1:1,2,3\n",
            Options::default(),
        );
        assert!(!report.is_valid());
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.cardinality")
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.genotype-allele")
        );
    }

    #[test]
    fn evidence_is_an_explicit_policy() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\tC\t.\tPASS\tDP=3\n",
            Options {
                require_evidence: true,
                ..Options::default()
            },
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.evidence")
        );
    }

    #[test]
    fn v44_cross_field_invariants_are_enforced() {
        let report = check_text(
            "##fileformat=VCFv4.4\n\
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Length\">\n\
##INFO=<ID=SVCLAIM,Number=A,Type=String,Description=\"Claim\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
1\t1\t.\tA\t<DEL>\t.\tPASS\tSVLEN=1;SVCLAIM=JD\tGT\t0/1\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.svclaim")
        );
    }

    #[test]
    fn invalid_float_is_not_accepted_for_oracle_compatibility() {
        let report = check_text(
            "##fileformat=VCFv4.4\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT:GL\t0/1:-0.13,-0r.58,-3.62\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.value-type")
        );
    }

    #[test]
    fn missing_alt_has_zero_alternate_values() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele counts\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\t.\t.\tPASS\tAC=1\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.cardinality")
        );
    }

    #[test]
    fn records_are_sorted_in_contiguous_contig_blocks() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\t.\tA\tC\t.\tPASS\t.\n\
1\t1\t.\tA\tC\t.\tPASS\t.\n\
2\t1\t.\tA\tC\t.\tPASS\t.\n\
1\t3\t.\tA\tC\t.\tPASS\t.\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.position-order")
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.contig-order")
        );
    }

    #[test]
    fn normalized_duplicate_variants_are_rejected() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tTAT\tTGT\t.\tPASS\t.\n\
1\t101\t.\tA\tG\t.\tPASS\t.\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.duplicate-variant")
        );
    }

    #[test]
    fn modern_info_keys_require_a_valid_prefix() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\tC\t.\tPASS\t42=1;1000G\n",
            Options::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|item| item.code == "record.info-key")
        );
        assert!(
            !report
                .diagnostics
                .iter()
                .any(|item| { item.code == "record.info-key" && item.message.contains("1000G") })
        );
    }

    #[test]
    fn non_structural_symbolic_alleles_do_not_require_svlen() {
        let report = check_text(
            "##fileformat=VCFv4.4\n\
##ALT=<ID=R,Description=\"IUPAC ambiguity code\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\t<R>\t.\tPASS\t.\n",
            Options::default(),
        );
        assert!(report.is_valid(), "{:?}", report.diagnostics);
    }

    #[test]
    fn float_keywords_are_case_insensitive() {
        let report = check_text(
            "##fileformat=VCFv4.3\n\
##INFO=<ID=X,Number=1,Type=Float,Description=\"Value\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\tC\t.\tPASS\tX=-INFINITY\n",
            Options::default(),
        );
        assert!(report.is_valid(), "{:?}", report.diagnostics);
    }
}
