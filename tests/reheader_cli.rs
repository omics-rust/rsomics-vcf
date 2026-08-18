use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use flate2::Compression;
use flate2::write::GzEncoder;

const HEADER: &[u8] = b"##fileformat=VCFv4.3\n\
##source=input\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

const BCF_SOURCE: &[u8] = b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100,IDX=2>\n\
##FILTER=<ID=q10,Description=\"low quality\",IDX=4>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\",IDX=2>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\",IDX=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t1\t.\tA\tC\t.\tq10\tDP=7\tGT\t0/1\t0/0\n";

struct Fixture {
    _directory: tempfile::TempDir,
    input: PathBuf,
    replacement: PathBuf,
    fai: PathBuf,
    samples: PathBuf,
    body: Vec<u8>,
}

struct BcfFixture {
    _directory: tempfile::TempDir,
    input: PathBuf,
    replacement: PathBuf,
    replacement_without_dp: PathBuf,
    fai: PathBuf,
    fai_without_used_contig: PathBuf,
    samples: PathBuf,
    expected_body: Vec<u8>,
}

impl BcfFixture {
    fn new(encoding: &str) -> Self {
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("source.vcf");
        let input = directory.path().join(format!("input-{encoding}.bcf"));
        let replacement = directory.path().join("replacement.vcfh");
        let replacement_without_dp = directory.path().join("without-dp.vcfh");
        let fai = directory.path().join("reference.fai");
        let fai_without_used_contig = directory.path().join("without-chr1.fai");
        let samples = directory.path().join("samples.tsv");
        fs::write(&source, BCF_SOURCE).unwrap();
        success(
            command()
                .args(["view", "-O", encoding, "-o"])
                .arg(&input)
                .arg(&source)
                .output()
                .unwrap(),
        );
        fs::write(
            &replacement,
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=200>\n\
##contig=<ID=chr2,length=300>\n\
##FILTER=<ID=q10,Description=\"low quality\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n\
##INFO=<ID=XX,Number=1,Type=Integer,Description=\"new\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\tR2\n",
        )
        .unwrap();
        fs::write(
            &replacement_without_dp,
            b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=200>\n\
##FILTER=<ID=q10,Description=\"low quality\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1\tR2\n",
        )
        .unwrap();
        fs::write(&fai, b"chr1\t1000\nchr2\t2000\n").unwrap();
        fs::write(&fai_without_used_contig, b"chr2\t2000\n").unwrap();
        fs::write(&samples, b"R1\tN1\nR2\tN2\n").unwrap();
        Self {
            _directory: directory,
            input,
            replacement,
            replacement_without_dp,
            fai,
            fai_without_used_contig,
            samples,
            expected_body: body(BCF_SOURCE).to_vec(),
        }
    }
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

fn bgzf_bytes(chunks: &[&[u8]]) -> Vec<u8> {
    let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
    for chunk in chunks {
        writer.write_all(chunk).unwrap();
        writer.flush().unwrap();
    }
    writer.finish().unwrap()
}

fn raw_frames(source: &[u8]) -> Vec<&[u8]> {
    let mut frames = Vec::new();
    let mut offset = 0;
    while offset < source.len() {
        assert!(source.len() - offset >= 18);
        let length = usize::from(u16::from_le_bytes([
            source[offset + 16],
            source[offset + 17],
        ])) + 1;
        assert!(source.len() - offset >= length);
        frames.push(&source[offset..offset + length]);
        offset += length;
    }
    frames
}

fn inflate(source: &[u8]) -> Vec<u8> {
    let mut reader = noodles_bgzf::io::Reader::new(source);
    let mut output = Vec::new();
    reader.read_to_end(&mut output).unwrap();
    output
}

fn body(source: &[u8]) -> &[u8] {
    let mut offset = 0;
    for line in source.split_inclusive(|byte| *byte == b'\n') {
        offset += line.len();
        if line.starts_with(b"#CHROM\t") {
            return &source[offset..];
        }
    }
    panic!("VCF header has no #CHROM line")
}

fn bcf_header(source: &[u8]) -> String {
    let raw = if source.starts_with(&[0x1f, 0x8b]) {
        inflate(source)
    } else {
        source.to_vec()
    };
    assert_eq!(&raw[..5], b"BCF\x02\x02");
    let length = u32::from_le_bytes(raw[5..9].try_into().unwrap()) as usize;
    let text = &raw[9..9 + length];
    let text = text.strip_suffix(&[0]).unwrap_or(text);
    String::from_utf8(text.to_vec()).unwrap()
}

fn decode_bcf(source: &[u8]) -> Vec<u8> {
    let file = tempfile::NamedTempFile::new().unwrap();
    fs::write(file.path(), source).unwrap();
    success(
        command()
            .args(["view", "-O", "v"])
            .arg(file.path())
            .output()
            .unwrap(),
    )
    .stdout
}

fn definition<'a>(header: &'a str, kind: &str, id: &str) -> &'a str {
    let prefix = format!("##{kind}=<");
    let id = format!("ID={id}");
    header
        .lines()
        .find(|line| line.starts_with(&prefix) && line.split([',', '<', '>']).any(|v| v == id))
        .unwrap_or_else(|| panic!("missing {kind}/{id} in {header}"))
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
fn nonzero_threads_are_rejected_for_vcf() {
    let fixture = Fixture::new();
    let compressed = fixture._directory.path().join("input.vcf.gz");
    fs::write(
        &compressed,
        bgzf_bytes(&[&fs::read(&fixture.input).unwrap()]),
    )
    .unwrap();
    for input in [&fixture.input, &compressed] {
        let output = command()
            .args(["reheader", "-n", "N1,N2", "--threads", "1"])
            .arg(input)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("BGZF BCF"),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn bgzf_rewrites_only_the_prefix_frames() {
    let fixture = Fixture::new();
    let mut later_body = Vec::new();
    for position in 2..20_000 {
        writeln!(
            later_body,
            "chr1\t{position}\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0"
        )
        .unwrap();
    }
    let first = [HEADER, fixture.body.as_slice()].concat();
    let input = bgzf_bytes(&[&first, &later_body]);
    let input_path = fixture._directory.path().join("input.vcf.gz");
    fs::write(&input_path, &input).unwrap();
    let original_frames = raw_frames(&input);
    let unchanged_tail = original_frames[1..].concat();

    let output = success(
        command()
            .args(["reheader", "-n", "N1,N2"])
            .arg(input_path)
            .output()
            .unwrap(),
    );
    let inflated = inflate(&output.stdout);
    assert!(inflated.windows(6).any(|window| window == b"N1\tN2\n"));
    assert_eq!(body(&inflated), [fixture.body, later_body].concat());
    assert!(output.stdout.ends_with(&unchanged_tail));
}

#[test]
fn bgzf_handles_a_header_spanning_multiple_frames() {
    let fixture = Fixture::new();
    let mut header = b"##fileformat=VCFv4.3\n".to_vec();
    for index in 0..6000 {
        writeln!(header, "##meta{index}=abcdefghijklmnop").unwrap();
    }
    header.extend_from_slice(
        b"##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
    );
    let input = bgzf_bytes(&[&[header.as_slice(), fixture.body.as_slice()].concat()]);
    let input_path = fixture._directory.path().join("large-header.vcf.gz");
    fs::write(&input_path, input).unwrap();
    let output = success(
        command()
            .args(["reheader", "-n", "N1,N2"])
            .arg(input_path)
            .output()
            .unwrap(),
    );
    let inflated = inflate(&output.stdout);
    assert!(inflated.windows(6).any(|window| window == b"N1\tN2\n"));
    assert_eq!(body(&inflated), fixture.body);
}

#[test]
fn bgzf_header_only_stdin_emits_one_complete_eof() {
    let input = bgzf_bytes(&[HEADER]);
    let output = success(with_stdin(&["reheader", "-n", "N1,N2"], &input));
    let frames = raw_frames(&output.stdout);
    assert_eq!(
        frames
            .iter()
            .filter(|frame| **frame == crate_eof_block())
            .count(),
        1
    );
    assert!(inflate(&output.stdout).ends_with(b"N1\tN2\n"));
}

fn crate_eof_block() -> &'static [u8] {
    &[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 27, 0, 3, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
    ]
}

#[test]
fn malformed_bgzf_tails_do_not_replace_an_existing_destination() {
    let fixture = Fixture::new();
    let first = [HEADER, fixture.body.as_slice()].concat();
    let valid = bgzf_bytes(&[
        first.as_slice(),
        b"chr1\t2\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n",
    ]);
    let first_frame = raw_frames(&valid)[0].len();
    let mut cases = Vec::new();

    let mut missing_eof = valid.clone();
    missing_eof.truncate(missing_eof.len() - crate_eof_block().len());
    cases.push(missing_eof);
    cases.push(valid[..first_frame + 19].to_vec());

    let mut invalid_size = valid.clone();
    invalid_size[first_frame + 16..first_frame + 18].copy_from_slice(&0u16.to_le_bytes());
    cases.push(invalid_size);

    let mut trailing = valid;
    trailing.extend_from_slice(b"trailing");
    cases.push(trailing);

    for (index, case) in cases.into_iter().enumerate() {
        let input = fixture
            ._directory
            .path()
            .join(format!("broken-{index}.vcf.gz"));
        let output_path = fixture
            ._directory
            .path()
            .join(format!("output-{index}.vcf.gz"));
        fs::write(&input, case).unwrap();
        fs::write(&output_path, b"keep").unwrap();
        let output = command()
            .args(["reheader", "-n", "N1,N2", "-o"])
            .arg(&output_path)
            .arg(&input)
            .output()
            .unwrap();
        assert!(!output.status.success(), "case={index}");
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(error.contains("BGZF"), "case={index}: {error}");
        assert!(!error.contains("not available"), "case={index}: {error}");
        assert_eq!(fs::read(output_path).unwrap(), b"keep", "case={index}");
    }
}

#[test]
fn raw_and_bgzf_bcf_preserve_records_encoding_and_dictionary_indices() {
    for encoding in ["u", "b"] {
        let fixture = BcfFixture::new(encoding);
        let output_path = fixture
            ._directory
            .path()
            .join(format!("output-{encoding}.bcf"));
        success(
            command()
                .args(["reheader", "-H"])
                .arg(&fixture.replacement)
                .arg("-f")
                .arg(&fixture.fai)
                .arg("-N")
                .arg(&fixture.samples)
                .arg("-o")
                .arg(&output_path)
                .arg(&fixture.input)
                .output()
                .unwrap(),
        );
        let output = fs::read(&output_path).unwrap();
        assert_eq!(output.starts_with(&[0x1f, 0x8b]), encoding == "b");
        assert_eq!(body(&decode_bcf(&output)), fixture.expected_body);

        let header = bcf_header(&output);
        assert!(definition(&header, "contig", "chr1").contains("IDX=2"));
        assert!(definition(&header, "contig", "chr2").contains("IDX=3"));
        assert!(definition(&header, "FORMAT", "GT").contains("IDX=1"));
        assert!(definition(&header, "INFO", "DP").contains("IDX=2"));
        assert!(definition(&header, "FILTER", "q10").contains("IDX=4"));
        assert!(definition(&header, "INFO", "XX").contains("IDX=5"));
        assert!(header.contains("\tN1\tN2\n"), "{header}");
    }
}

#[test]
fn removing_used_bcf_definitions_fails_without_replacing_output() {
    let fixture = BcfFixture::new("b");
    for (index, arguments) in [
        ("contig", &fixture.fai_without_used_contig),
        ("info", &fixture.replacement_without_dp),
    ] {
        let output_path = fixture
            ._directory
            .path()
            .join(format!("existing-{index}.bcf"));
        fs::write(&output_path, b"keep").unwrap();
        let option = if index == "contig" { "-f" } else { "-H" };
        let output = command()
            .args(["reheader", option])
            .arg(arguments)
            .arg("-o")
            .arg(&output_path)
            .arg(&fixture.input)
            .output()
            .unwrap();
        assert!(!output.status.success(), "{index}");
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(!error.contains("panicked"), "{error}");
        assert!(!error.contains("not available"), "{error}");
        assert!(error.contains("record"), "{error}");
        assert_eq!(fs::read(output_path).unwrap(), b"keep", "{index}");
    }
}

#[test]
fn truncated_bcf_and_missing_bgzf_eof_fail_transactionally() {
    for encoding in ["u", "b"] {
        let fixture = BcfFixture::new(encoding);
        let broken = fixture
            ._directory
            .path()
            .join(format!("broken-{encoding}.bcf"));
        let mut bytes = fs::read(&fixture.input).unwrap();
        if encoding == "b" {
            bytes.truncate(bytes.len() - crate_eof_block().len());
        } else {
            bytes.truncate(bytes.len() - 3);
        }
        fs::write(&broken, bytes).unwrap();
        let output_path = fixture
            ._directory
            .path()
            .join(format!("existing-{encoding}.bcf"));
        fs::write(&output_path, b"keep").unwrap();
        let output = command()
            .args(["reheader", "-n", "N1,N2", "-o"])
            .arg(&output_path)
            .arg(&broken)
            .output()
            .unwrap();
        assert!(!output.status.success(), "{encoding}");
        let error = String::from_utf8_lossy(&output.stderr);
        assert!(!error.contains("not available"), "{encoding}: {error}");
        assert_eq!(fs::read(output_path).unwrap(), b"keep", "{encoding}");
    }
}

#[test]
fn bcf_threads_and_standard_input_follow_the_encoding_contract() {
    let raw = BcfFixture::new("u");
    let rejected = command()
        .args(["reheader", "-n", "N1,N2", "--threads", "1"])
        .arg(&raw.input)
        .output()
        .unwrap();
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("BGZF BCF"));

    let compressed = BcfFixture::new("b");
    let parallel = compressed._directory.path().join("parallel.bcf");
    success(
        command()
            .args(["reheader", "-n", "N1,N2", "--threads", "2", "-o"])
            .arg(&parallel)
            .arg(&compressed.input)
            .output()
            .unwrap(),
    );
    let parallel = fs::read(parallel).unwrap();
    assert!(parallel.starts_with(&[0x1f, 0x8b]));
    assert_eq!(body(&decode_bcf(&parallel)), compressed.expected_body);

    for fixture in [raw, compressed] {
        let input = fs::read(&fixture.input).unwrap();
        let output = success(with_stdin(&["reheader", "-n", "N1,N2"], &input));
        assert_eq!(
            output.stdout.starts_with(&[0x1f, 0x8b]),
            input.starts_with(&[0x1f, 0x8b])
        );
        assert_eq!(body(&decode_bcf(&output.stdout)), fixture.expected_body);
    }
}
