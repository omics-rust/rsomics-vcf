use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/view.vcf")
}

fn body(region: &str) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-vcf"))
        .args(["view", "--targets", region, "--no-header"])
        .arg(fixture())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn single_coordinates_are_exact_and_trailing_hyphens_are_open_ended() {
    let exact = body("chr1:20");
    assert_eq!(exact.lines().count(), 1, "{exact}");
    assert!(exact.starts_with("chr1\t20\t"), "{exact}");

    let open = body("chr1:20-");
    assert_eq!(open.lines().count(), 3, "{open}");
    assert!(open.starts_with("chr1\t20\t"), "{open}");
    assert!(open.contains("chr1\t30\t"), "{open}");
}
