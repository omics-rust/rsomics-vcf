use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use flate2::Compression;
use flate2::write::GzEncoder;

const HEADER: &[u8] = b"##fileformat=VCFv4.3\n\
##source=input\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

struct Fixture {
    _directory: tempfile::TempDir,
    input: PathBuf,
    replacement: PathBuf,
    fai: PathBuf,
    samples: PathBuf,
    body: Vec<u8>,
}

impl Fixture {
    fn new() -> Self {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("input.vcf");
        let replacement = directory.path().join("replacement.vcfh");
        let fai = directory.path().join("reference.fai");
        let samples = directory.path().join("samples.tsv");
        let body = b"chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\r\n".to_vec();
        fs::write(&input, [HEADER, body.as_slice()].concat()).unwrap();
        fs::write(
            &replacement,
            b"##fileformat=VCFv4.3\r\n\
##source=replacement\r\n\
##contig=<ID=chr1,length=20,assembly=test>\r\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\tR2\r\n",
        )
        .unwrap();
        fs::write(&fai, b"chr1\t1000\t0\t0\t0\nchr2\t2000\t0\t0\t0\n").unwrap();
        fs::write(&samples, b"R1\tTumor\nR2\tNormal\n").unwrap();
        Self {
            _directory: directory,
            input,
            replacement,
            fai,
            samples,
            body,
        }
    }
}

fn command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn success(output: Output) -> Output {
    assert!(
        output.status.success(),
        "stdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn with_stdin(arguments: &[&str], input: &[u8]) -> Output {
    let mut child = command()
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
fn help_uses_the_family_surface_and_requires_an_edit() {
    let output = success(command().args(["reheader", "--help"]).output().unwrap());
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("-H, --header <FILE>"), "{help}");
    assert!(help.contains("-f, --fai <FILE>"), "{help}");
    assert!(help.contains("-n, --samples-list <LIST>"), "{help}");
    assert!(help.contains("-N, --samples-file <FILE>"), "{help}");
    assert!(help.contains("--threads <INT>"), "{help}");

    let fixture = Fixture::new();
    let output = command()
        .arg("reheader")
        .arg(&fixture.input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("required"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn composes_header_fai_and_samples_without_touching_the_body() {
    let fixture = Fixture::new();
    let output = success(
        command()
            .args(["reheader", "-H"])
            .arg(&fixture.replacement)
            .arg("-f")
            .arg(&fixture.fai)
            .args(["-n", "Tumor,Normal"])
            .arg(&fixture.input)
            .output()
            .unwrap(),
    );
    let expected_header = b"##fileformat=VCFv4.3\n\
##source=replacement\n\
##contig=<ID=chr1,length=1000,assembly=test>\n\
##contig=<ID=chr2,length=2000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal\n";
    assert_eq!(
        output.stdout,
        [expected_header.as_slice(), fixture.body.as_slice()].concat()
    );
}

#[test]
fn reads_plain_vcf_from_standard_input() {
    let fixture = Fixture::new();
    let input = fs::read(&fixture.input).unwrap();
    let output = success(with_stdin(&["reheader", "-n", "N1,N2"], &input));
    assert!(output.stdout.starts_with(b"##fileformat=VCFv4.3\n"));
    assert!(output.stdout.windows(6).any(|window| window == b"N1\tN2\n"));
    assert!(output.stdout.ends_with(&fixture.body));
}

#[test]
fn json_summary_is_separate_from_named_variant_output() {
    let fixture = Fixture::new();
    let output_path = fixture._directory.path().join("output.vcf");
    let output = success(
        command()
            .args(["reheader", "--json", "-n", "N1,N2", "-o"])
            .arg(&output_path)
            .arg(&fixture.input)
            .output()
            .unwrap(),
    );
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "reheader");
    assert_eq!(envelope["result"]["summary"]["encoding"], "plain-vcf");
    assert_eq!(envelope["result"]["summary"]["samples_before"], 2);
    assert_eq!(envelope["result"]["summary"]["samples_after"], 2);
    assert!(fs::read(output_path).unwrap().ends_with(&fixture.body));
}

#[test]
fn json_rejects_variant_output_on_standard_output() {
    let fixture = Fixture::new();
    let output = command()
        .args(["reheader", "--json", "-n", "N1,N2"])
        .arg(&fixture.input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(!output.stdout.starts_with(b"##fileformat"));
    let diagnostic = [output.stdout, output.stderr].concat();
    assert!(
        String::from_utf8_lossy(&diagnostic).contains("--json requires --output"),
        "{}",
        String::from_utf8_lossy(&diagnostic)
    );
}

#[test]
fn output_cannot_alias_any_input() {
    let fixture = Fixture::new();
    let cases: [(&Path, &[&str]); 4] = [
        (&fixture.input, &["-n", "N1,N2"]),
        (&fixture.replacement, &["-H"]),
        (&fixture.fai, &["-f"]),
        (&fixture.samples, &["-N"]),
    ];
    for (output_path, arguments) in cases {
        let before = fs::read(output_path).unwrap();
        let mut process = command();
        process.arg("reheader").args(arguments);
        if arguments.len() == 1 {
            process.arg(output_path);
        }
        let output = process
            .arg("-o")
            .arg(output_path)
            .arg(&fixture.input)
            .output()
            .unwrap();
        assert!(!output.status.success(), "{}", output_path.display());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("also an input path"),
            "{}: {}",
            output_path.display(),
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(fs::read(output_path).unwrap(), before);
    }
}

#[test]
fn a_failed_edit_does_not_replace_an_existing_destination() {
    let fixture = Fixture::new();
    let output_path = fixture._directory.path().join("existing.vcf");
    fs::write(&output_path, b"keep").unwrap();
    let output = command()
        .args(["reheader", "-n", "only-one", "-o"])
        .arg(&output_path)
        .arg(&fixture.input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("2 names"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read(output_path).unwrap(), b"keep");
}

#[test]
fn rejects_ordinary_gzip_with_a_conversion_diagnostic() {
    let fixture = Fixture::new();
    let gzip = fixture._directory.path().join("input.vcf.gz");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder
        .write_all(&fs::read(&fixture.input).unwrap())
        .unwrap();
    fs::write(&gzip, encoder.finish().unwrap()).unwrap();
    let output = command()
        .args(["reheader", "-n", "N1,N2"])
        .arg(gzip)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let error = String::from_utf8_lossy(&output.stderr);
    assert!(error.contains("ordinary gzip"), "{error}");
    assert!(error.contains("convert"), "{error}");
}

#[test]
fn nonzero_threads_are_rejected_for_plain_vcf() {
    let fixture = Fixture::new();
    let output = command()
        .args(["reheader", "-n", "N1,N2", "--threads", "1"])
        .arg(&fixture.input)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("BGZF BCF"),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
