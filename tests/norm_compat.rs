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
    fs::write(&reference, b">chr1\nAAAAAACGTACGT\n").unwrap();
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
