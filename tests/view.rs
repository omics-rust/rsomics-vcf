use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture(directory: &Path) -> PathBuf {
    let path = directory.join("calls.vcf");
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/view.vcf"),
        &path,
    )
    .unwrap();
    path
}

fn run(arguments: &[&str]) -> std::process::Output {
    let output = Command::new(binary()).args(arguments).output().unwrap();
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn bgzip(source: &[u8], destination: &Path) {
    let mut writer = bgzf::io::Writer::new(File::create(destination).unwrap());
    writer.write_all(source).unwrap();
    writer.try_finish().unwrap();
}

#[test]
fn projects_samples_and_recalculates_ac_an() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(directory.path());
    let output = run(&["view", input.to_str().unwrap(), "--samples", "S2,S1"]);
    let output = String::from_utf8(output.stdout).unwrap();

    assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1"));
    assert!(output.contains("chr1\t10\trs1\tA\tG\t50\tPASS\tAC=3;AN=4\tGT\t1/1\t0/1"));
    assert!(output.contains("chr1\t20\t.\tAT\tA\t8\tq10\tAC=1;AN=4\tGT\t0/1\t0/0"));
}

#[test]
fn sample_selection_rejects_duplicates_and_can_drop_missing_names() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(directory.path());
    let input = input.to_str().unwrap();

    let duplicate = Command::new(binary())
        .args(["view", input, "-s", "S1,S1"])
        .output()
        .unwrap();
    assert!(!duplicate.status.success());
    assert!(String::from_utf8_lossy(&duplicate.stderr).contains("duplicate"));

    let output = run(&["view", input, "-s", "missing", "--force-samples"]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(!output.contains("##FORMAT="));
    assert!(
        output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .all(|line| line.split('\t').count() == 8)
    );
}

#[test]
fn filters_records_with_complete_type_categories() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(directory.path());
    let input = input.to_str().unwrap();

    let output = run(&["view", input, "--types", "indels", "--no-header"]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("chr1\t20\t.\tAT\tA"));
    assert!(output.contains("chr1\t25\tins\tA\tAT"));
    assert_eq!(output.lines().count(), 2);

    let output = run(&[
        "view",
        input,
        "--exclude-types",
        "snps,mnps,indels",
        "--no-header",
    ]);
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "chr2\t15\tref\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\t0/0\n"
    );

    let output = run(&[
        "view",
        input,
        "--known",
        "--apply-filters",
        "PASS",
        "--min-alleles",
        "2",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("\trs1\t"));
    assert!(output.contains("\tmnp\t"));
    assert!(!output.contains("\tref\t"));
}

#[test]
fn converts_between_all_output_encodings() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(directory.path());

    for (kind, extension) in [("z", "vcf.gz"), ("b", "bcf"), ("u", "raw.bcf")] {
        let encoded = directory.path().join(extension);
        let decoded = directory.path().join(format!("{extension}.vcf"));
        run(&[
            "view",
            input.to_str().unwrap(),
            "-O",
            kind,
            "-o",
            encoded.to_str().unwrap(),
        ]);
        run(&[
            "view",
            encoded.to_str().unwrap(),
            "-O",
            "v",
            "-o",
            decoded.to_str().unwrap(),
        ]);
        let output = fs::read_to_string(decoded).unwrap();
        assert_eq!(
            output.lines().filter(|line| !line.starts_with('#')).count(),
            5
        );
        assert!(output.contains("chr1\t10\trs1\tA\tG"));
    }
}

#[test]
fn targets_stream_and_header_modes_are_explicit() {
    let directory = tempfile::tempdir().unwrap();
    let input = fixture(directory.path());
    let input = input.to_str().unwrap();

    let output = run(&[
        "view",
        input,
        "--targets",
        "chr1:15-24",
        "--targets-overlap",
        "record",
        "--no-header",
    ]);
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "chr1\t20\t.\tAT\tA\t8\tq10\t.\tGT\t0/0\t0/1\t./.\n"
    );

    let output = run(&[
        "view",
        input,
        "--targets",
        "^chr1:15-24",
        "--targets-overlap",
        "record",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert_eq!(output.lines().count(), 4);
    assert!(!output.contains("chr1\t20\t"));

    let targets = directory.path().join("targets.txt");
    fs::write(&targets, "chr1:15-24\n").unwrap();
    let output = run(&[
        "view",
        input,
        "--targets-file",
        &format!("^{}", targets.display()),
        "--targets-overlap",
        "record",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert_eq!(output.lines().count(), 4);
    assert!(!output.contains("chr1\t20\t"));

    let output = run(&["view", input, "--header-only"]);
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .all(|line| line.starts_with('#'))
    );

    let result = Command::new(binary())
        .args(["view", input, "-O", "b", "--no-header"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("requires a header"));
}

#[test]
fn malformed_input_does_not_replace_named_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("bad.vcf");
    let output = directory.path().join("output.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\tbad\t.\tA\tG\t.\tPASS\t.\n",
    )
    .unwrap();
    fs::write(&output, b"existing").unwrap();

    let result = Command::new(binary())
        .args([
            "view",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(output).unwrap(), b"existing");
}

#[test]
fn indexed_regions_include_spanning_records_without_duplicates() {
    let directory = tempfile::tempdir().unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/index.vcf");
    let input = directory.path().join("calls.vcf.gz");
    bgzip(&fs::read(source).unwrap(), &input);
    let input = input.to_str().unwrap();
    run(&["index", input]);

    let output = run(&[
        "view",
        input,
        "--regions",
        "chr1:70000-70000,chr1:69990-70010",
        "--regions-overlap",
        "record",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert_eq!(output.lines().count(), 2);
    assert_eq!(output.matches("\tdel1\t").count(), 1);
    assert_eq!(output.matches("\tsnv4\t").count(), 1);

    let output = run(&[
        "view",
        input,
        "--regions",
        "chr1:30000-30000,chr1:70000-70000",
        "--regions-overlap",
        "record",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert_eq!(output.matches("\tdel1\t").count(), 1);

    let output = run(&[
        "view",
        input,
        "--regions",
        "chr1:70000-70000",
        "--regions-overlap",
        "pos",
        "--no-header",
    ]);
    let output = String::from_utf8(output.stdout).unwrap();
    assert_eq!(output.lines().count(), 1);
    assert!(output.contains("\tsnv4\t"));
}
