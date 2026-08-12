#![cfg(feature = "norm-preview")]

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

fn run(command: &mut Command) -> Output {
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn body(output: Output) -> Vec<u8> {
    output
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.starts_with(b"#") && !line.is_empty())
        .flat_map(|line| line.iter().copied().chain([b'\n']))
        .collect()
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_realignments_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nARAAAACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tA\tC\t.\tPASS\t.\n\
chr1\t4\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t4\t.\tAA\tA\t.\tPASS\t.\n\
chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn typed_multiallelic_split_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=FA,Number=A,Type=Integer,Description=\"A\">\n\
##FORMAT=<ID=FR,Number=R,Type=Integer,Description=\"R\">\n\
##FORMAT=<ID=FG,Number=G,Type=Integer,Description=\"G\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\tIA=10,20;IR=5,3,2;IG=0,10,20,30,40,50\tGT:FA:FR:FG\t1/2:11,22:7,4,3:0,10,20,30,40,50\t2:33,44:8,5,6:0,10,20\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args(["norm", "--split-multiallelic", input.to_str().unwrap()])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn split_and_reference_realign_compose_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nAAAAAACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tAA,AAA\t.\tPASS\tAF=0.25,0.5\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--split-multiallelic",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "-m",
        "-any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_warn_and_skip_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

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

    for (ours_policy, oracle_policy) in [("warn", "w"), ("skip", "x")] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            ours_policy,
            input.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            oracle_policy,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{ours_policy}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn biallelic_mnv_atomization_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t10\t.\tACGT\tTGCA\t.\tPASS\tDP=5\tGT\t0/1\n\
chr1\t20\t.\tACGT\tAGGA\t.\tPASS\tDP=7\tGT\t1/1\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args(["norm", "--atomize", input.to_str().unwrap()])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn split_and_atomize_compose_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tAC\tGT,AT\t.\tPASS\tAF=0.25,0.5\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--split-multiallelic",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}
