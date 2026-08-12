use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

fn run(arguments: &[&str]) -> std::process::Output {
    Command::new(binary()).args(arguments).output().unwrap()
}

#[test]
fn filters_typed_expressions_from_the_public_command() {
    let output = run(&[
        "filter",
        "--include",
        "QUAL >= 10",
        fixture().to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    let records: Vec<_> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(records.len(), 1, "{output}");
    assert!(records[0].starts_with("chr1\t10\t"), "{output}");
}

#[test]
fn compressed_filter_output_is_transactional_and_separate_from_json() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("filtered.vcf.gz");
    let output = run(&[
        "filter",
        "--json",
        "--include",
        "QUAL >= 10",
        "--output-type",
        "z",
        "--threads",
        "2",
        "--output",
        output_path.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "filter");
    assert_eq!(envelope["result"]["summary"]["read"], 3);
    assert_eq!(envelope["result"]["summary"]["written"], 1);

    let mut decoded = String::new();
    bgzf::io::Reader::new(File::open(&output_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert!(decoded.contains("chr1\t10\t"), "{decoded}");
    assert!(!decoded.contains("chr1\t30\t"), "{decoded}");

    fs::write(&output_path, b"existing").unwrap();
    let failed = run(&[
        "filter",
        "--include",
        "INFO/UNDECLARED > 0",
        "--output",
        output_path.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);
    assert!(!failed.status.success());
    assert_eq!(fs::read(output_path).unwrap(), b"existing");
}

#[test]
fn filter_rejects_parallel_uncompressed_output() {
    let output = run(&[
        "filter",
        "--include",
        "QUAL >= 10",
        "--threads",
        "2",
        fixture().to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("compression workers require BGZF VCF or BCF output")
    );
}

#[test]
fn soft_filter_modes_and_failed_genotypes_reach_the_public_command() {
    let output = run(&[
        "filter",
        "--include",
        "FMT/DP >= 10",
        "--soft-filter",
        "LowDepth",
        "--mode",
        "+x",
        "--set-GTs",
        ".",
        fixture().to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("##FILTER=<ID=LowDepth"), "{output}");
    assert!(
        output.contains("chr1\t10\t.\tA\tC\t50.5\tPASS\tDP=7\tGT:DP\t0/1:20\t./.:."),
        "{output}"
    );
    assert!(
        output.contains("chr1\t20\trs2\tG\tT\t.\tLowDepth\t.\tGT:DP\t./.:5\t./.:8"),
        "{output}"
    );
}

#[test]
fn masks_and_gap_filters_are_exposed_without_separate_crates() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tC\t10\tPASS\t.\n\
chr1\t12\t.\tA\tAT\t10\tPASS\t.\n\
chr1\t20\t.\tG\tT\t10\tPASS\t.\n",
    )
    .unwrap();
    let output = run(&[
        "filter",
        "--mask",
        "chr1:20-20",
        "--soft-filter",
        "Masked",
        "--SnpGap",
        "3:indel",
        input.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("##FILTER=<ID=Masked"), "{output}");
    assert!(output.contains("##FILTER=<ID=SnpGap"), "{output}");
    assert!(output.contains("chr1\t10\t.\tA\tC\t10\tSnpGap"), "{output}");
    assert!(output.contains("chr1\t20\t.\tG\tT\t10\tMasked"), "{output}");
}

#[test]
fn every_filter_output_encoding_round_trips() {
    let directory = tempfile::tempdir().unwrap();
    for (kind, extension) in [
        ("v", "vcf"),
        ("z", "vcf.gz"),
        ("b", "bcf"),
        ("u", "raw.bcf"),
    ] {
        let encoded = directory.path().join(extension);
        let decoded = directory.path().join(format!("{extension}.decoded.vcf"));
        let mut arguments = vec![
            "filter",
            "--include",
            "QUAL >= 10",
            "--output-type",
            kind,
            "--output",
            encoded.to_str().unwrap(),
        ];
        if matches!(kind, "z" | "b") {
            arguments.extend(["--threads", "2"]);
        }
        let input = fixture();
        arguments.push(input.to_str().unwrap());
        let filtered = run(&arguments);
        assert!(
            filtered.status.success(),
            "{}: {}",
            kind,
            String::from_utf8_lossy(&filtered.stderr)
        );
        let viewed = run(&[
            "view",
            "--output",
            decoded.to_str().unwrap(),
            encoded.to_str().unwrap(),
        ]);
        assert!(
            viewed.status.success(),
            "{}: {}",
            kind,
            String::from_utf8_lossy(&viewed.stderr)
        );
        let decoded = fs::read_to_string(decoded).unwrap();
        assert_eq!(
            decoded
                .lines()
                .filter(|line| !line.starts_with('#'))
                .count(),
            1,
            "{kind}: {decoded}"
        );
    }
}
