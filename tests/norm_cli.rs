#![cfg(feature = "norm-preview")]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture(directory: &Path) -> (PathBuf, PathBuf) {
    let reference = directory.join("reference.fa");
    let input = directory.join("input.vcf");
    fs::write(&reference, b">chr1\nAAAAAACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.\n",
    )
    .unwrap();
    (reference, input)
}

#[test]
fn public_command_normalizes_and_reports_json_separately() {
    let directory = tempfile::tempdir().unwrap();
    let (reference, input) = fixture(directory.path());
    let output_path = directory.path().join("normalized.vcf");
    let output = Command::new(binary())
        .args([
            "--json",
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "norm");
    assert_eq!(envelope["result"]["summary"]["read"], 2);
    assert_eq!(envelope["result"]["summary"]["changed"], 2);
    let normalized = fs::read_to_string(output_path).unwrap();
    assert!(normalized.contains("chr1\t1\t.\tA\tAA"), "{normalized}");
    assert!(normalized.contains("chr1\t11\t.\tC\tG"), "{normalized}");
}

#[test]
fn failed_normalization_does_not_replace_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let (reference, input) = fixture(directory.path());
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tT\tA\t.\tPASS\t.\n",
    )
    .unwrap();
    let output_path = directory.path().join("normalized.vcf");
    fs::write(&output_path, b"existing").unwrap();

    let output = Command::new(binary())
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(output_path).unwrap(), b"existing");
}
