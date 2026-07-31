use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/view.vcf")
}

fn index_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/index.vcf")
}

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

fn assert_bcftools_1_24() {
    let output = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&output.stdout).starts_with("bcftools 1.24\n"));
}

fn bgzip(source: &[u8], destination: &Path) {
    let mut writer = bgzf::io::Writer::new(File::create(destination).unwrap());
    writer.write_all(source).unwrap();
    writer.try_finish().unwrap();
}

fn body_ours(input: &Path, arguments: &[&str]) -> Vec<u8> {
    let mut command = Command::new(binary());
    command.arg("view").arg(input);
    command.args(arguments).arg("--no-header");
    run(&mut command).stdout
}

fn body_bcftools(input: &Path, arguments: &[&str]) -> Vec<u8> {
    let mut command = Command::new("bcftools");
    command
        .args(["view", "--no-version", "-H"])
        .args(arguments)
        .arg(input);
    run(&mut command).stdout
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn selection_and_sample_projection_match_bcftools() {
    assert_bcftools_1_24();
    let input = fixture();
    for (ours, oracle) in [
        (vec!["-s", "S2,S1"], vec!["-s", "S2,S1"]),
        (vec!["-s", "^S3"], vec!["-s", "^S3"]),
        (vec!["-v", "snps,mnps"], vec!["-v", "snps,mnps"]),
        (vec!["-V", "snps,mnps"], vec!["-V", "snps,mnps"]),
        (vec!["-f", "PASS"], vec!["-f", "PASS"]),
        (vec!["--known"], vec!["--known"]),
        (vec!["--novel"], vec!["--novel"]),
        (vec!["-m", "2", "-M", "2"], vec!["-m", "2", "-M", "2"]),
        (vec!["-s", "S2,S1", "-G"], vec!["-s", "S2,S1", "-G"]),
    ] {
        assert_eq!(
            body_ours(&input, &ours),
            body_bcftools(&input, &oracle),
            "{ours:?}"
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn indexed_overlap_modes_match_bcftools() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.vcf.gz");
    bgzip(&fs::read(index_fixture()).unwrap(), &input);
    run(Command::new(binary()).args(["index", input.to_str().unwrap()]));

    for (ours_mode, oracle_mode) in [("pos", "0"), ("record", "1"), ("variant", "2")] {
        let ours = body_ours(
            &input,
            &["-r", "chr1:70000-70000", "--regions-overlap", ours_mode],
        );
        let oracle = body_bcftools(
            &input,
            &["-r", "chr1:70000-70000", "--regions-overlap", oracle_mode],
        );
        assert_eq!(ours, oracle, "{ours_mode}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn literal_variant_overlap_matches_bcftools() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.vcf.gz");
    bgzip(&fs::read(fixture()).unwrap(), &input);
    run(Command::new(binary()).args(["index", input.to_str().unwrap()]));

    for region in ["chr1:20-20", "chr1:21-21", "chr1:25-25", "chr1:26-26"] {
        let ours = body_ours(&input, &["-r", region, "--regions-overlap", "variant"]);
        let oracle = body_bcftools(&input, &["-r", region, "--regions-overlap", "2"]);
        assert_eq!(ours, oracle, "{region}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn every_output_encoding_round_trips_through_bcftools() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = fixture();

    for (kind, extension) in [("z", "vcf.gz"), ("b", "bcf"), ("u", "raw.bcf")] {
        let output = directory.path().join(extension);
        run(Command::new(binary()).args([
            "view",
            input.to_str().unwrap(),
            "-O",
            kind,
            "-o",
            output.to_str().unwrap(),
        ]));
        assert_eq!(
            body_bcftools(&output, &[]),
            body_bcftools(&input, &[]),
            "{kind}"
        );
    }
}
