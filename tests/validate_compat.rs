use std::fs;
use std::path::{Path, PathBuf};

use rsomics_vcf::validate::{self, Options, Severity};

#[test]
#[ignore = "release oracle: requires the pinned hts-specs repository"]
fn hts_specs_vcf_corpus() {
    let root = std::env::var_os("HTS_SPECS_DIR")
        .map(PathBuf::from)
        .expect("HTS_SPECS_DIR must point to the pinned hts-specs repository");
    let root = root.join("test/vcf");
    check_corpus(&root, &["4.2", "4.3", "4.5"], 450, known_hts_fixture_defect);
}

#[test]
#[ignore = "release oracle: requires the pinned EBI vcf-validator repository"]
fn ebi_vcf_corpus() {
    let root = std::env::var_os("VCF_VALIDATOR_DIR")
        .map(PathBuf::from)
        .expect("VCF_VALIDATOR_DIR must point to the pinned EBI vcf-validator repository");
    let root = root.join("test/input_files");
    check_corpus(
        &root,
        &["v4.1", "v4.2", "v4.3", "v4.4"],
        950,
        known_ebi_fixture_defect,
    );
}

fn check_corpus(
    root: &Path,
    versions: &[&str],
    minimum: usize,
    known_defect: fn(&Path, &validate::Report) -> bool,
) {
    let mut mismatches = Vec::new();
    let mut checked = 0;

    for version in versions {
        let version = root.join(version);
        for (directory, expected) in [("passed", true), ("failed", false)] {
            let directory = version.join(directory);
            if !directory.is_dir() {
                continue;
            }
            for path in files(&directory) {
                if path.extension().and_then(|value| value.to_str()) != Some("vcf") {
                    continue;
                }
                let report = validate::check(
                    &path,
                    Options {
                        max_diagnostics: 8,
                        require_evidence: false,
                    },
                )
                .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
                checked += 1;
                if report.is_valid() != expected && !known_defect(&path, &report) {
                    mismatches.push(format!(
                        "{}: expected {}, errors={}, diagnostics={:?}",
                        path.strip_prefix(root).unwrap().display(),
                        if expected { "valid" } else { "invalid" },
                        report.errors,
                        report.diagnostics
                    ));
                }
            }
        }
    }

    assert!(checked >= minimum, "only checked {checked} fixtures");
    assert!(
        mismatches.is_empty(),
        "{} of {checked} fixtures disagreed:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}

fn known_ebi_fixture_defect(path: &Path, report: &validate::Report) -> bool {
    if [
        "v4.3/failed/failed_body_chrom_003.vcf",
        "v4.3/failed/failed_meta_contig_003.vcf",
        "v4.4/failed/failed_body_chrom_003.vcf",
        "v4.4/failed/failed_meta_contig_003.vcf",
    ]
    .iter()
    .any(|suffix| path.ends_with(suffix))
    {
        return report.errors == 0;
    }

    if known_missing_alt_fixture_defect(path, report) {
        return true;
    }

    let defect = if path.ends_with("v4.1/passed/passed_meta_pedigree.vcf")
        || path.ends_with("v4.2/passed/passed_meta_pedigree.vcf")
    {
        Some(("record.empty", 5, None))
    } else if path.ends_with("v4.3/passed/passed_body_format.vcf") {
        Some(("record.format-key", 7, Some("FORMAT")))
    } else if path.ends_with("v4.3/passed/passed_body_samples.vcf") {
        Some(("record.empty", 15, None))
    } else if path.ends_with("v4.3/passed/passed_meta_pedigree.vcf") {
        Some(("record.empty", 8, None))
    } else if path.ends_with("v4.4/passed/passed_body_alt.vcf") {
        Some(("record.value-type", 12, Some("GL")))
    } else {
        None
    };
    let Some((code, line, field)) = defect else {
        return false;
    };
    report.errors == 1
        && report.diagnostics.iter().any(|item| {
            item.code == code && item.line == Some(line) && item.field.as_deref() == field
        })
}

fn known_hts_fixture_defect(path: &Path, report: &validate::Report) -> bool {
    if path.ends_with("4.5/passed/zero_length_LAA.vcf") {
        return report.errors == 1
            && report.diagnostics.iter().any(|item| {
                item.severity == Severity::Error
                    && item.code == "record.position-order"
                    && item.line == Some(8)
                    && item.field.as_deref() == Some("POS")
            });
    }

    ([
        "4.2/failed/failed_body_chrom_001.vcf",
        "4.3/failed/failed_body_chrom_001.vcf",
        "4.3/failed/failed_body_chrom_004.vcf",
        "4.3/failed/failed_meta_contig_003.vcf",
    ]
    .iter()
    .any(|suffix| path.ends_with(suffix))
        && report.errors == 0)
        || known_missing_alt_fixture_defect(path, report)
}

fn known_missing_alt_fixture_defect(path: &Path, report: &validate::Report) -> bool {
    let complex = [
        "4.2/passed/complexfile_passed_000.vcf",
        "4.3/passed/complexfile_passed_000.vcf",
        "v4.1/passed/complexfile_passed_000.vcf",
        "v4.2/passed/complexfile_passed_000.vcf",
        "v4.3/passed/complexfile_passed_000.vcf",
        "v4.4/passed/complexfile_passed_000.vcf",
    ]
    .iter()
    .any(|suffix| path.ends_with(suffix));
    if complex {
        return report.errors == 102
            && report
                .diagnostics
                .iter()
                .filter(|item| item.severity == Severity::Error)
                .all(|item| {
                    item.line == Some(55)
                        && matches!(
                            item.code.as_str(),
                            "record.cardinality" | "record.genotype-allele"
                        )
                });
    }

    let alternate = [
        "4.2/passed/passed_body_alt.vcf",
        "4.3/passed/passed_body_alt.vcf",
        "v4.1/passed/passed_body_alt.vcf",
        "v4.2/passed/passed_body_alt.vcf",
        "v4.3/passed/passed_body_alt.vcf",
        "v4.4/passed/passed_body_alt.vcf",
    ]
    .iter()
    .any(|suffix| path.ends_with(suffix));
    if !alternate {
        return false;
    }

    let version_44 = path.ends_with("v4.4/passed/passed_body_alt.vcf");
    report.errors == if version_44 { 6 } else { 5 }
        && report
            .diagnostics
            .iter()
            .filter(|item| item.severity == Severity::Error)
            .all(|item| {
                (item.line == Some(22)
                    && matches!(
                        item.code.as_str(),
                        "record.cardinality" | "record.genotype-allele"
                    ))
                    || (version_44
                        && item.line == Some(12)
                        && item.code == "record.value-type"
                        && item.field.as_deref() == Some("GL"))
            })
}

fn files(root: &Path) -> Vec<PathBuf> {
    let mut pending = vec![root.to_owned()];
    let mut output = Vec::new();
    while let Some(directory) = pending.pop() {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else {
                output.push(path);
            }
        }
    }
    output.sort();
    output
}
