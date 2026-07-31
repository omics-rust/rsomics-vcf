use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

fn invalid_fixture() -> tempfile::NamedTempFile {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    input
        .write_all(
            b"##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\tbad\t.\tA\tC\t.\tPASS\t.\n",
        )
        .unwrap();
    input
}

#[test]
fn valid_input_prints_diagnostics_and_summary() {
    let output = Command::new(binary())
        .arg("validate")
        .arg(fixture())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "valid: 3 records, 0 errors, 1 warning\n"
    );
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("warning[header.reference]"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn invalid_input_exits_one_with_line_and_field_context() {
    let input = invalid_fixture();
    let output = Command::new(binary())
        .arg("validate")
        .arg(input.path())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        stderr.contains("validation failed: 1 record, 1 error"),
        "{stderr}"
    );
    assert!(
        stderr.contains("error[record.pos] line 3, field POS"),
        "{stderr}"
    );
}

#[test]
fn json_failure_keeps_the_structured_report() {
    let input = invalid_fixture();
    let output = Command::new(binary())
        .args(["validate", "--json"])
        .arg(input.path())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());

    let envelope: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
    assert_eq!(envelope["status"], "error");
    assert_eq!(envelope["exit_code"], 1);
    assert_eq!(envelope["report"]["command"], "validate");
    assert_eq!(envelope["report"]["report"]["records"], 1);
    assert_eq!(envelope["report"]["report"]["errors"], 1);
    assert_eq!(
        envelope["report"]["report"]["diagnostics"][1]["field"],
        "POS"
    );
}

#[test]
fn json_success_uses_the_common_success_envelope() {
    let output = Command::new(binary())
        .args(["validate", "--json"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());

    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["status"], "ok");
    assert_eq!(envelope["result"]["command"], "validate");
    assert_eq!(envelope["result"]["report"]["records"], 3);
}
