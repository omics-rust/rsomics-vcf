use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

fn query_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/query.vcf")
}

fn run(mut command: Command) -> Output {
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

fn bcftools(arguments: &[&str]) -> Output {
    let mut command = Command::new("bcftools");
    command.args(arguments);
    run(command)
}

fn ours(arguments: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-vcf"));
    command.args(arguments);
    run(command)
}

fn assert_bcftools_1_24() {
    let output = bcftools(&["--version"]);
    let version = String::from_utf8(output.stdout).unwrap();
    assert!(version.starts_with("bcftools 1.24\n"), "{version}");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn head_matches_bcftools_1_24_for_vcf_bgzf_and_bcf() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let vcf = fixture();
    let gz = directory.path().join("head.vcf.gz");
    let bcf = directory.path().join("head.bcf");
    let vcf_text = vcf.to_str().unwrap();
    let gz_text = gz.to_str().unwrap();
    let bcf_text = bcf.to_str().unwrap();

    bcftools(&["view", "--no-version", "-Oz", "-o", gz_text, vcf_text]);
    bcftools(&["view", "--no-version", "-Ob", "-o", bcf_text, vcf_text]);

    for input in [vcf_text, gz_text, bcf_text] {
        for (our_options, oracle_options) in [
            (vec![], vec![]),
            (vec!["-H", "2"], vec!["-h", "2"]),
            (vec!["-n", "2"], vec!["-n", "2"]),
            (vec!["-s", "0"], vec!["-s", "0"]),
            (vec!["-s", "2"], vec!["-s", "2"]),
            (vec!["-H", "2", "-s", "1"], vec!["-h", "2", "-s", "1"]),
        ] {
            let mut our_arguments = vec!["head"];
            our_arguments.extend(our_options);
            our_arguments.push(input);
            let mut oracle_arguments = vec!["head"];
            oracle_arguments.extend(oracle_options);
            oracle_arguments.push(input);

            assert_eq!(
                ours(&our_arguments).stdout,
                bcftools(&oracle_arguments).stdout,
                "{input} {our_arguments:?}"
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn head_stdin_matches_bcftools_1_24() {
    assert_bcftools_1_24();
    let input = fs::read(fixture()).unwrap();

    let mut our_command = Command::new(env!("CARGO_BIN_EXE_rsomics-vcf"));
    our_command
        .args(["head", "-n", "2"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped());
    let mut our_child = our_command.spawn().unwrap();
    our_child.stdin.take().unwrap().write_all(&input).unwrap();
    let our_output = our_child.wait_with_output().unwrap();
    assert!(our_output.status.success());

    let mut oracle_command = Command::new("bcftools");
    oracle_command
        .args(["head", "-n", "2"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped());
    let mut oracle_child = oracle_command.spawn().unwrap();
    oracle_child
        .stdin
        .take()
        .unwrap()
        .write_all(&input)
        .unwrap();
    let oracle_output = oracle_child.wait_with_output().unwrap();
    assert!(oracle_output.status.success());

    assert_eq!(our_output.stdout, oracle_output.stdout);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn head_canonicalizes_bcf_numeric_values() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let bcf = directory.path().join("query.bcf");
    let input = query_fixture();
    bcftools(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);

    let arguments = ["head", "-n", "3", bcf.to_str().unwrap()];
    assert_eq!(ours(&arguments).stdout, bcftools(&arguments).stdout);
}
