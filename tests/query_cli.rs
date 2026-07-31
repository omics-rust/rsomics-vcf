use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

fn command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

#[test]
fn query_writes_the_selected_fields() {
    let output = command()
        .args(["query", "-f", r"%CHROM\t%POS[\t%SAMPLE=%GT]\n", "-s", "S2"])
        .arg(fixture())
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t10\tS2=0/0\nchr1\t20\tS2=0/1\nchr1\t30\tS2=0/1\n"
    );
}

#[test]
fn failed_query_does_not_replace_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("result.tsv");
    fs::write(&output_path, b"existing\n").unwrap();

    let output = command()
        .args(["query", "-f", "%INFO/UNDECLARED", "-o"])
        .arg(&output_path)
        .arg(fixture())
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert_eq!(fs::read(output_path).unwrap(), b"existing\n");
}

#[test]
fn json_requires_named_query_output() {
    let output = command()
        .args(["query", "--json", "-f", "%CHROM"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("--json requires --output"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn json_summary_is_separate_from_named_query_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("result.tsv");
    let output = command()
        .args(["query", "--json", "-f", "%CHROM", "-o"])
        .arg(&output_path)
        .arg(fixture())
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read(output_path).unwrap(), b"chr1\nchr1\nchr1\n");
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "query");
    assert_eq!(envelope["result"]["summary"]["records"], 3);
}
