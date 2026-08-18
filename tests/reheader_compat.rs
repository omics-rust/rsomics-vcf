use std::ffi::OsStr;
use std::fs;
use std::io::{Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use flate2::read::MultiGzDecoder;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn oracle() -> PathBuf {
    std::env::var_os("RSOMICS_BCFTOOLS")
        .map(PathBuf::from)
        .expect("RSOMICS_BCFTOOLS must point to bcftools 1.24")
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bcftools-reheader")
        .join(name)
}

fn run(command: &mut Command) -> Output {
    command.output().unwrap()
}

fn require_success(output: Output, context: impl std::fmt::Display) -> Output {
    assert!(
        output.status.success(),
        "{context}: status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn assert_oracle_version() {
    let output = require_success(run(Command::new(oracle()).arg("--version")), "version");
    assert!(
        String::from_utf8_lossy(&output.stdout).starts_with("bcftools 1.24\n"),
        "{}",
        String::from_utf8_lossy(&output.stdout)
    );
}

#[derive(Clone, Copy, Debug)]
enum TestEncoding {
    PlainVcf,
    BgzfVcf,
    RawBcf,
    BgzfBcf,
}

impl TestEncoding {
    fn bcftools_flag(self) -> &'static str {
        match self {
            Self::PlainVcf => "-Ov",
            Self::BgzfVcf => "-Oz",
            Self::RawBcf => "-Ou",
            Self::BgzfBcf => "-Ob",
        }
    }

    fn extension(self) -> &'static str {
        match self {
            Self::PlainVcf => "vcf",
            Self::BgzfVcf => "vcf.gz",
            Self::RawBcf => "raw.bcf",
            Self::BgzfBcf => "bcf",
        }
    }

    fn is_vcf(self) -> bool {
        matches!(self, Self::PlainVcf | Self::BgzfVcf)
    }
}

#[derive(Clone, Copy, Debug)]
enum EditCase {
    Header,
    Fai,
    PositionalSamples,
    PairSamples,
    Composed,
}

impl EditCase {
    fn name(self) -> &'static str {
        match self {
            Self::Header => "header",
            Self::Fai => "fai",
            Self::PositionalSamples => "positional",
            Self::PairSamples => "pairs",
            Self::Composed => "composed",
        }
    }
}

fn create_input(encoding: TestEncoding, destination: &Path) {
    require_success(
        run(Command::new(oracle()).args([
            OsStr::new("view"),
            OsStr::new("--no-version"),
            OsStr::new(encoding.bcftools_flag()),
            OsStr::new("-o"),
            destination.as_os_str(),
            fixture("input.vcf").as_os_str(),
        ])),
        format_args!("creating {encoding:?} input"),
    );
}

fn add_edit(command: &mut Command, edit: EditCase, ours: bool) {
    let header_flag = if ours { "-H" } else { "-h" };
    match edit {
        EditCase::Header => {
            command.arg(header_flag).arg(fixture("replacement.vcfh"));
        }
        EditCase::Fai => {
            command.arg("-f").arg(fixture("reference.fai"));
        }
        EditCase::PositionalSamples => {
            command.args(["-n", "Tumor,Normal"]);
        }
        EditCase::PairSamples => {
            command.arg("-N").arg(fixture("samples.txt"));
        }
        EditCase::Composed => {
            command
                .arg(header_flag)
                .arg(fixture("replacement.vcfh"))
                .arg("-f")
                .arg(fixture("reference.fai"))
                .arg("-N")
                .arg(fixture("samples.txt"));
        }
    }
}

fn decoded_vcf(source: &[u8], encoding: TestEncoding) -> Vec<u8> {
    if matches!(encoding, TestEncoding::PlainVcf) {
        return source.to_vec();
    }
    let mut decoded = Vec::new();
    MultiGzDecoder::new(Cursor::new(source))
        .read_to_end(&mut decoded)
        .unwrap();
    decoded
}

#[derive(Debug)]
struct DecodedBcf {
    header: noodles_vcf::Header,
    records: Vec<u8>,
}

fn canonical_bcf(path: &Path) -> DecodedBcf {
    let output = require_success(
        run(Command::new(oracle()).args([
            OsStr::new("view"),
            OsStr::new("--no-version"),
            OsStr::new("-Ov"),
            path.as_os_str(),
        ])),
        format_args!("decoding {}", path.display()),
    );
    let mut header = Vec::new();
    let mut records = Vec::new();
    for line in output.stdout.split(|byte| *byte == b'\n') {
        if line.is_empty() || line.starts_with(b"##bcftools_") {
            continue;
        }
        if line.starts_with(b"#") {
            header.extend_from_slice(line);
            header.push(b'\n');
        } else {
            records.extend_from_slice(line);
            records.push(b'\n');
        }
    }
    DecodedBcf {
        header: String::from_utf8(header).unwrap().parse().unwrap(),
        records,
    }
}

fn assert_equivalent(encoding: TestEncoding, edit: EditCase) {
    let directory = tempfile::tempdir().unwrap();
    let input = directory
        .path()
        .join(format!("input.{}", encoding.extension()));
    let ours = directory
        .path()
        .join(format!("{}.ours.{}", edit.name(), encoding.extension()));
    let expected =
        directory
            .path()
            .join(format!("{}.bcftools.{}", edit.name(), encoding.extension()));
    create_input(encoding, &input);

    let mut ours_command = Command::new(binary());
    ours_command.arg("reheader");
    add_edit(&mut ours_command, edit, true);
    ours_command.arg("-o").arg(&ours).arg(&input);
    require_success(
        run(&mut ours_command),
        format_args!("rsomics {encoding:?} {edit:?}"),
    );

    let mut oracle_command = Command::new(oracle());
    oracle_command.arg("reheader");
    add_edit(&mut oracle_command, edit, false);
    oracle_command.arg("-o").arg(&expected).arg(&input);
    require_success(
        run(&mut oracle_command),
        format_args!("bcftools {encoding:?} {edit:?}"),
    );

    if encoding.is_vcf() {
        assert_eq!(
            decoded_vcf(&fs::read(&ours).unwrap(), encoding),
            decoded_vcf(&fs::read(&expected).unwrap(), encoding),
            "{encoding:?} {edit:?}"
        );
    } else {
        let ours = canonical_bcf(&ours);
        let expected = canonical_bcf(&expected);
        assert_eq!(ours.header, expected.header, "{encoding:?} {edit:?} header");
        assert_eq!(
            ours.records, expected.records,
            "{encoding:?} {edit:?} records"
        );
    }
}

fn run_stdin(command: &mut Command, input: &[u8]) -> Output {
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

fn assert_stdin_equivalent(encoding: TestEncoding, edit: EditCase) {
    let directory = tempfile::tempdir().unwrap();
    let input = directory
        .path()
        .join(format!("stdin.{}", encoding.extension()));
    create_input(encoding, &input);
    let input = fs::read(input).unwrap();

    let mut ours = Command::new(binary());
    ours.arg("reheader");
    add_edit(&mut ours, edit, true);
    ours.arg("-");
    let ours = require_success(run_stdin(&mut ours, &input), "rsomics stdin");

    let mut expected = Command::new(oracle());
    expected.arg("reheader");
    add_edit(&mut expected, edit, false);
    expected.arg("-");
    let expected = require_success(run_stdin(&mut expected, &input), "bcftools stdin");

    assert_eq!(
        decoded_vcf(&ours.stdout, encoding),
        decoded_vcf(&expected.stdout, encoding),
        "stdin {encoding:?} {edit:?}"
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn declared_success_matrix_matches_bcftools_1_24() {
    assert_oracle_version();
    for encoding in [
        TestEncoding::PlainVcf,
        TestEncoding::BgzfVcf,
        TestEncoding::RawBcf,
        TestEncoding::BgzfBcf,
    ] {
        for edit in [
            EditCase::Header,
            EditCase::Fai,
            EditCase::PositionalSamples,
            EditCase::PairSamples,
            EditCase::Composed,
        ] {
            assert_equivalent(encoding, edit);
        }
    }
    assert_stdin_equivalent(TestEncoding::PlainVcf, EditCase::Composed);
    assert_stdin_equivalent(TestEncoding::BgzfVcf, EditCase::Composed);
}

fn assert_controlled_failure(output: Output, case: &str) {
    assert!(!output.status.success(), "{case}");
    assert!(
        output.status.code().is_none_or(|code| code < 128),
        "{case}: {:?}",
        output.status
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("panicked"), "{case}: {stderr}");
    assert!(!stderr.contains("not available"), "{case}: {stderr}");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn fail_loud_divergences_are_controlled() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let plain = directory.path().join("input.vcf");
    let bgzf = directory.path().join("input.vcf.gz");
    let bcf = directory.path().join("input.bcf");
    create_input(TestEncoding::PlainVcf, &plain);
    create_input(TestEncoding::BgzfVcf, &bgzf);
    create_input(TestEncoding::BgzfBcf, &bcf);

    assert_controlled_failure(
        run(Command::new(binary()).args([
            OsStr::new("reheader"),
            OsStr::new("-n"),
            OsStr::new("OnlyOne"),
            plain.as_os_str(),
        ])),
        "sample count",
    );

    let unknown = directory.path().join("unknown.tsv");
    fs::write(&unknown, b"missing\tN1\n").unwrap();
    assert_controlled_failure(
        run(Command::new(binary()).args([
            OsStr::new("reheader"),
            OsStr::new("-N"),
            unknown.as_os_str(),
            plain.as_os_str(),
        ])),
        "unknown sample",
    );

    let duplicate = directory.path().join("duplicate.tsv");
    fs::write(&duplicate, b"S1\tN\nS2\tN\n").unwrap();
    assert_controlled_failure(
        run(Command::new(binary()).args([
            OsStr::new("reheader"),
            OsStr::new("-N"),
            duplicate.as_os_str(),
            plain.as_os_str(),
        ])),
        "duplicate final sample",
    );

    let missing_contig = directory.path().join("missing-contig.fai");
    fs::write(&missing_contig, b"chr2\t200\n").unwrap();
    assert_controlled_failure(
        run(Command::new(binary()).args([
            OsStr::new("reheader"),
            OsStr::new("-f"),
            missing_contig.as_os_str(),
            bcf.as_os_str(),
        ])),
        "used BCF contig removed",
    );

    let truncated = directory.path().join("truncated.vcf.gz");
    let mut bytes = fs::read(bgzf).unwrap();
    bytes.truncate(bytes.len() - 10);
    fs::write(&truncated, bytes).unwrap();
    assert_controlled_failure(
        run(Command::new(binary()).args([
            OsStr::new("reheader"),
            OsStr::new("-n"),
            OsStr::new("Tumor,Normal"),
            truncated.as_os_str(),
        ])),
        "truncated BGZF",
    );
}
