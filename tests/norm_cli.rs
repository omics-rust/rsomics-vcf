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

#[test]
fn public_command_splits_typed_multiallelic_records_without_a_reference() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\tAF=0.25,0.5\tGT\t1/2\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args(["norm", "--split-multiallelic", input.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(
        output.contains("A\tC\t.\tPASS\tAF=0.25\tGT\t1/0"),
        "{output}"
    );
    assert!(
        output.contains("A\tG\t.\tPASS\tAF=0.5\tGT\t0/1"),
        "{output}"
    );
}

#[test]
fn public_command_can_use_missing_for_split_overlaps() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT\t1|2\t./2\t0/2\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--split-multiallelic",
            "--split-overlaps",
            "missing",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("GT\t1|.\t./.\t0/."), "{output}");
    assert!(output.contains("GT\t.|1\t./1\t0/1"), "{output}");
}

#[test]
fn public_command_preserves_ad_sums_while_splitting() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:10,3,2\t0/2:10,.,2\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--split-multiallelic",
            "--keep-sum",
            "AD",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("GT:AD\t1/0:12,3\t0/0:12,."), "{output}");
    assert!(output.contains("GT:AD\t0/1:13,2\t0/1:10,2"), "{output}");
}

#[test]
fn reference_mismatch_warn_and_skip_are_observable() {
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tT\tA\t.\tPASS\t.\n\
chr1\t3\t.\tG\tA\t.\tPASS\t.\n",
    )
    .unwrap();

    let warn = Command::new(binary())
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "warn",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(warn.status.success());
    assert!(String::from_utf8_lossy(&warn.stderr).contains("REF_MISMATCH\tchr1\t2"));
    assert_eq!(
        String::from_utf8_lossy(&warn.stdout)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .count(),
        2
    );

    let skip = Command::new(binary())
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "skip",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(skip.status.success());
    assert_eq!(
        String::from_utf8_lossy(&skip.stdout)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .count(),
        1
    );
}

#[test]
fn public_command_atomizes_mnvs_without_a_reference() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t20\t.\tACGT\tAGGA\t.\tPASS\tDP=7\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args(["norm", "--atomize", input.to_str().unwrap()])
        .output()
        .unwrap();
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
    assert_eq!(records.len(), 2, "{output}");
    assert!(records[0].starts_with("chr1\t21\t.\tC\tG"), "{output}");
    assert!(records[1].starts_with("chr1\t23\t.\tT\tA"), "{output}");
}

#[test]
fn public_command_can_use_missing_for_overlapping_atoms() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t50\t.\tCC\tC,GG\t.\tPASS\t.\tGT\t1/2\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--atomize",
            "--atom-overlaps",
            ".",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("C\tG\t.\tPASS\t.\tGT\t./1"), "{output}");
    assert!(output.contains("CC\tC\t.\tPASS\t.\tGT\t1/."), "{output}");
    assert!(!output.contains("*"), "{output}");
}

#[test]
fn public_command_can_trace_atoms_to_the_original_record() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t50\t.\tCC\tC,GG\t.\tPASS\t.\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--atomize",
            "--old-rec-tag",
            "ORIG",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(
        output.contains("##INFO=<ID=ORIG,Number=1,Type=String,Description=\"Original variant."),
        "{output}"
    );
    assert!(
        output.contains("C\tG,*\t.\tPASS\tORIG=chr1|50|CC|C,GG|2"),
        "{output}"
    );
    assert!(
        output.contains("CC\tC,*\t.\tPASS\tORIG=chr1|50|CC|C,GG|1"),
        "{output}"
    );
}

#[test]
fn public_command_removes_duplicates_by_explicit_policy() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.4\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SV length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ts1\tA\tC\t.\tPASS\t.\n\
chr1\t10\ts2\tA\tG\t.\tPASS\t.\n\
chr1\t10\ti1\tA\tAT\t.\tPASS\t.\n\
chr1\t10\ti2\tA\tAG\t.\tPASS\t.\n\
chr1\t10\to1\tA\t<DEL>\t.\tPASS\tSVLEN=-5\n\
chr1\t20\te1\tA\tC\t.\tPASS\t.\n\
chr1\t20\te2\ta\tc\t.\tPASS\t.\n\
chr1\t20\te3\tA\tG\t.\tPASS\t.\n\
chr1\t30\tx1\tN\t<DEL>\t.\tPASS\tSVLEN=-10\n\
chr1\t30\tx2\tN\t<DEL>\t.\tPASS\tSVLEN=-20\n\
chr1\t30\tx3\tN\t<DEL>\t.\tPASS\tSVLEN=-10\n\
chr1\t40\tm1\tA\tC,G\t.\tPASS\t.\n\
chr1\t40\tm2\tA\tG,C\t.\tPASS\t.\n\
chr1\t40\tm3\tA\tC,C\t.\tPASS\t.\n",
    )
    .unwrap();

    for (policy, expected) in [
        (
            "snps",
            vec!["s1", "i1", "i2", "o1", "e1", "x1", "x2", "x3", "m1"],
        ),
        (
            "indels",
            vec![
                "s1", "s2", "i1", "o1", "e1", "e2", "e3", "x1", "x2", "x3", "m1", "m2", "m3",
            ],
        ),
        ("both", vec!["s1", "i1", "o1", "e1", "x1", "x2", "x3", "m1"]),
        ("all", vec!["s1", "e1", "x1", "m1"]),
        (
            "exact",
            vec!["s1", "s2", "i1", "i2", "o1", "e1", "e3", "x1", "x2", "m1"],
        ),
    ] {
        let output = Command::new(binary())
            .args([
                "norm",
                "--remove-duplicates",
                policy,
                input.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{policy}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = String::from_utf8(output.stdout).unwrap();
        let ids: Vec<_> = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| line.split('\t').nth(2))
            .collect();
        assert_eq!(ids, expected, "{policy}");
    }
}
