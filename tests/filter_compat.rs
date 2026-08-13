use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
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

fn body(output: Output) -> Vec<u8> {
    output
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.starts_with(b"#") && !line.is_empty())
        .flat_map(|line| line.iter().copied().chain(*b"\n"))
        .collect()
}

fn body_ours(input: &Path, arguments: &[&str]) -> Vec<u8> {
    let mut command = Command::new(binary());
    command.arg("filter").args(arguments).arg(input);
    body(run(&mut command))
}

fn body_bcftools(input: &Path, arguments: &[&str]) -> Vec<u8> {
    let mut command = Command::new("bcftools");
    command
        .args(["filter", "--no-version"])
        .args(arguments)
        .arg(input);
    body(run(&mut command))
}

fn bgzip(source: &[u8], destination: &Path) {
    let mut writer = bgzf::io::Writer::new(File::create(destination).unwrap());
    writer.write_all(source).unwrap();
    writer.try_finish().unwrap();
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn expressions_soft_filters_and_failed_genotypes_match_bcftools() {
    assert_bcftools_1_24();
    let input = fixture("head.vcf");
    for arguments in [
        vec!["-i", "QUAL >= 10"],
        vec!["-e", "QUAL < 10"],
        vec!["-i", "INFO/DP >= 7 && QUAL >= 30"],
        vec!["-i", "POS + 1 > 20"],
        vec!["-e", "QUAL / 2 >= 10"],
        vec!["-e", "QUAL < 10", "-s", "LowQual"],
        vec!["-i", "FMT/DP >= 10", "-s", "LowDepth", "-S", "."],
        vec!["-i", "FMT/DP >= 10", "-s", "LowDepth", "-S", "0"],
    ] {
        assert_eq!(
            body_ours(&input, &arguments),
            body_bcftools(&input, &arguments),
            "{arguments:?}"
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn masks_gap_filters_and_streaming_targets_match_bcftools() {
    assert_bcftools_1_24();
    let input = fixture("filter-gap.vcf");
    for (ours, oracle) in [
        (
            vec!["--mask", "chr1:20-20", "-s", "Masked"],
            vec!["--mask", "chr1:20-20", "-s", "Masked"],
        ),
        (vec!["-g", "3:indel"], vec!["-g", "3:indel"]),
        (vec!["-G", "3"], vec!["-G", "3"]),
        (
            vec!["-i", "QUAL >= 10", "-t", "^chr1:20-20"],
            vec![
                "-i",
                "QUAL >= 10",
                "-t",
                "^chr1:20-20",
                "--targets-overlap",
                "0",
            ],
        ),
        (
            vec!["-i", "QUAL >= 10", "-t", "chr1:10"],
            vec![
                "-i",
                "QUAL >= 10",
                "-t",
                "chr1:10",
                "--targets-overlap",
                "0",
            ],
        ),
        (
            vec!["-i", "QUAL >= 10", "-t", "chr1:10-"],
            vec![
                "-i",
                "QUAL >= 10",
                "-t",
                "chr1:10-",
                "--targets-overlap",
                "0",
            ],
        ),
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
fn indexed_regions_and_single_coordinates_match_bcftools() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.vcf.gz");
    bgzip(&fs::read(fixture("index.vcf")).unwrap(), &input);
    run(Command::new(binary()).args(["index", input.to_str().unwrap()]));

    for (ours, oracle) in [
        (
            vec!["-i", "QUAL >= 0", "-r", "chr1:70000"],
            vec![
                "-i",
                "QUAL >= 0",
                "-r",
                "chr1:70000",
                "--regions-overlap",
                "1",
            ],
        ),
        (
            vec!["-i", "QUAL >= 0", "-r", "chr1:70000-70000,chr1:69990-70010"],
            vec![
                "-i",
                "QUAL >= 0",
                "-r",
                "chr1:70000-70000,chr1:69990-70010",
                "--regions-overlap",
                "1",
            ],
        ),
    ] {
        assert_eq!(
            body_ours(&input, &ours),
            body_bcftools(&input, &oracle),
            "{ours:?}"
        );
    }
}
