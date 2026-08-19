use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/setgt.vcf")
}

fn run(arguments: &[&str]) -> std::process::Output {
    Command::new(binary()).args(arguments).output().unwrap()
}

fn run_stdin(arguments: &[&str], input: &[u8]) -> std::process::Output {
    let mut child = Command::new(binary())
        .args(arguments)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn edits_missing_genotypes_from_the_public_command() {
    let output = run(&[
        "setgt",
        "--target-gt",
        ".",
        "--new-gt",
        "0",
        fixture().to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("AC=1;AN=6\tGT:AD:DP:GQ\t0/1:4,4:8:30\t0/0:3,.:3:10\t0/0:.,.:.:."));
    assert!(output.contains("AC=1;AN=5\tGT:AD:DP:GQ\t0/1:.,.:5:20\t0/0:4,0:4:30\t0:.,.:.:."));
}

#[test]
fn json_summary_is_separate_from_transactional_variant_output() {
    let directory = tempfile::tempdir().unwrap();
    let destination = directory.path().join("set.vcf.gz");
    let output = run(&[
        "setgt",
        "--json",
        "-t",
        ".",
        "-n",
        "0",
        "-O",
        "z",
        "--threads",
        "2",
        "-o",
        destination.to_str().unwrap(),
        fixture().to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "setgt");
    assert_eq!(envelope["result"]["summary"]["read"], 3);
    assert_eq!(envelope["result"]["summary"]["changed_records"], 2);
    assert_eq!(envelope["result"]["summary"]["changed_genotypes"], 3);
    assert_eq!(envelope["result"]["summary"]["changed_alleles"], 4);

    let mut decoded = String::new();
    bgzf::io::Reader::new(File::open(destination).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert!(decoded.contains("chr1\t30"));
}

#[test]
fn argument_relationships_fail_before_variant_output() {
    let input = fixture();
    let input = input.to_str().unwrap();
    for arguments in [
        vec!["setgt", "-n", "0", input],
        vec!["setgt", "-t", "a", input],
        vec!["setgt", "-t", "q", "-n", "0", input],
        vec!["setgt", "-t", "a", "-n", "0", "-i", "QUAL>0", input],
        vec!["setgt", "-t", "a", "-n", "0", "--seed", "7", input],
        vec!["setgt", "-t", "a", "-t", ".", "-n", "0", input],
        vec!["setgt", "-t", "r:0.5", "-t", "r:0.2", "-n", "0", input],
        vec!["setgt", "-t", "a", "-n", "0u", input],
        vec![
            "setgt", "-t", "a", "-n", "0", "-i", "QUAL>0", "-e", "QUAL<0", input,
        ],
        vec!["setgt", "--json", "-t", "a", "-n", "0", input],
        vec!["setgt", "-t", "a", "-n", "0", "--threads", "2", input],
    ] {
        let output = run(&arguments);
        assert!(!output.status.success(), "{arguments:?}");
        assert!(output.stdout.is_empty(), "{arguments:?}");
    }
}

#[test]
fn named_output_survives_a_late_record_failure() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("late.vcf");
    let destination = directory.path().join("output.vcf");
    let mut source = fs::read_to_string(fixture()).unwrap();
    source.push_str("chr1\tBAD\n");
    fs::write(&input, source).unwrap();
    fs::write(&destination, b"existing").unwrap();

    let output = run(&[
        "setgt",
        "-t",
        "a",
        "-n",
        "0",
        "-o",
        destination.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(fs::read(destination).unwrap(), b"existing");
}

#[test]
fn named_output_rejects_the_input_path() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    let source = fs::read(fixture()).unwrap();
    fs::write(&input, &source).unwrap();

    let output = run(&[
        "setgt",
        "-t",
        "a",
        "-n",
        "0",
        "-o",
        input.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(fs::read(input).unwrap(), source);
}

#[test]
fn standard_input_and_output_share_the_typed_stream() {
    let input = fs::read(fixture()).unwrap();
    let output = run_stdin(&["setgt", "-t", "./x", "-n", "0"], &input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("0/0:3,.:3:10"));
    assert!(output.contains("./.:.,.:.:."));
}

#[test]
fn sites_only_is_a_no_op_and_sample_headers_require_gt() {
    let directory = tempfile::tempdir().unwrap();
    let sites = directory.path().join("sites.vcf");
    fs::write(
        &sites,
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t10\tPASS\t.\n",
    )
    .unwrap();
    let output = run(&["setgt", "-t", "a", "-n", "0", sites.to_str().unwrap()]);
    assert!(output.status.success());
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .contains("chr1\t1")
    );

    let samples = directory.path().join("samples.vcf");
    fs::write(
        &samples,
        "##fileformat=VCFv4.3\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t1\t.\tA\tC\t10\tPASS\t.\tDP\t3\n",
    )
    .unwrap();
    let output = run(&["setgt", "-t", "a", "-n", "0", samples.to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}

#[test]
fn every_output_encoding_round_trips() {
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
            "setgt",
            "-t",
            ".",
            "-n",
            "0",
            "-O",
            kind,
            "-o",
            encoded.to_str().unwrap(),
        ];
        if matches!(kind, "z" | "b") {
            arguments.extend(["--threads", "2"]);
        }
        let input = fixture();
        arguments.push(input.to_str().unwrap());
        let output = run(&arguments);
        assert!(
            output.status.success(),
            "{kind}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = run(&[
            "view",
            "-o",
            decoded.to_str().unwrap(),
            encoded.to_str().unwrap(),
        ]);
        assert!(
            output.status.success(),
            "{kind}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let decoded = fs::read_to_string(decoded).unwrap();
        assert!(decoded.contains("0/0:3,.:3:10"), "{kind}: {decoded}");
    }
}
