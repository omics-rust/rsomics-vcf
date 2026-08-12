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
fn public_command_joins_biallelic_snps_and_indels() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ts\tA\tC\t10\tPASS\t.\n\
chr1\t10\td\tAT\tA\t20\tPASS\t.\n\
chr1\t20\tx\tG\tA\t.\tPASS\t.\n\
chr1\t20\ty\tG\tT\t.\tPASS\t.\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--join-multiallelic",
            "any",
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
        output.contains("chr1\t10\ts;d\tAT\tCT,A\t20\tPASS\t."),
        "{output}"
    );
    assert!(
        output.contains("chr1\t20\tx;y\tG\tA,T\t.\tPASS\t."),
        "{output}"
    );
}

#[test]
fn strict_filter_uses_upstream_join_precedence() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FILTER=<ID=q10,Description=\"q10\">\n\
##FILTER=<ID=s20,Description=\"s20\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t.\tq10\t.\n\
chr1\t10\tb\tA\tG\t.\tPASS\t.\n\
chr1\t20\ta\tA\tC\t.\tPASS\t.\n\
chr1\t20\tb\tA\tG\t.\tq10\t.\n\
chr1\t30\ta\tA\tC\t.\tq10\t.\n\
chr1\t30\tb\tA\tG\t.\tPASS\t.\n\
chr1\t30\tc\tA\tT\t.\ts20\t.\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--join-multiallelic",
            "any",
            "--strict-filter",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let filters = String::from_utf8(output.stdout)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(6).unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(filters, ["PASS", "q10", "s20"]);
}

#[test]
fn classified_join_modes_keep_unselected_types_separate() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ti1\tA\tAT\t.\tPASS\t.\n\
chr1\t10\ts1\tA\tC\t.\tPASS\t.\n\
chr1\t10\ti2\tA\tAG\t.\tPASS\t.\n\
chr1\t10\ts2\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();
    for (mode, expected, joined) in [
        ("snps", vec!["s1;s2", "i1", "i2"], 1),
        ("indels", vec!["s1", "s2", "i1;i2"], 1),
        ("both", vec!["s1;s2", "i1;i2"], 2),
    ] {
        let output_path = directory.path().join(format!("{mode}.vcf"));
        let output = Command::new(binary())
            .args([
                "--json",
                "norm",
                "--join-multiallelic",
                mode,
                "--output",
                output_path.to_str().unwrap(),
                input.to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{mode}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let summary: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(summary["result"]["summary"]["joined"], joined, "{mode}");
        let output = fs::read_to_string(output_path).unwrap();
        let ids = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').nth(2).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(ids, expected, "{mode}");
    }
}

#[test]
fn join_rejects_incompatible_reference_alleles() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t.\tPASS\t.\n\
chr1\t10\tt\tT\tG\t.\tPASS\t.\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--join-multiallelic",
            "any",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("cannot join incompatible REF alleles")
    );
}

#[test]
fn join_rejects_invalid_allele_field_cardinality() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Count\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t.\tPASS\tAC=1,2\n\
chr1\t10\tt\tA\tG\t.\tPASS\tAC=3\n",
    )
    .unwrap();
    let output = Command::new(binary())
        .args([
            "norm",
            "--join-multiallelic",
            "any",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("INFO/AC has 2 values, expected 1"));
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
fn reference_fix_swaps_alleles_genotypes_and_allele_count() {
    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t2\t.\tT\tA,C\t.\tPASS\tAC=4,5\tGT\t0/2\t1|0\n",
    )
    .unwrap();

    let output = Command::new(binary())
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "fix",
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stdout)
            .contains("chr1\t2\t.\tC\tA,T\t.\tPASS\tAC=4,2\tGT\t2/0\t1|2\n")
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
