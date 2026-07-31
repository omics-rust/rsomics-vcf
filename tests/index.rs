use std::fs::{self, File};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/index.vcf")
}

fn bgzip(source: &[u8], destination: &Path) {
    let mut writer = bgzf::io::Writer::new(File::create(destination).unwrap());
    writer.write_all(source).unwrap();
    writer.try_finish().unwrap();
}

fn write_bcf(input: &Path, output: &Path) {
    let mut reader = vcf::io::Reader::new(BufReader::new(File::open(input).unwrap()));
    let header = reader.read_header().unwrap();
    let mut writer = bcf::io::Writer::new(File::create(output).unwrap());
    writer.write_header(&header).unwrap();

    for record in reader.records() {
        writer
            .write_variant_record(&header, &record.unwrap())
            .unwrap();
    }

    writer.try_finish().unwrap();
}

fn corrupt_bcf_reference_type(input: &Path, output: &Path) {
    let mut data = Vec::new();
    bgzf::io::Reader::new(File::open(input).unwrap())
        .read_to_end(&mut data)
        .unwrap();
    let header_length =
        usize::try_from(u32::from_le_bytes(data[5..9].try_into().unwrap())).unwrap();
    let record_start = 9 + header_length;
    let id_start = record_start + 8 + 24;
    let id_length = usize::from(data[id_start] >> 4);
    data[id_start + 1 + id_length] = 0x11;
    bgzip(&data, output);
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

#[test]
fn builds_csi_tbi_and_reports_counts() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.vcf.gz");
    bgzip(&fs::read(fixture()).unwrap(), &input);
    let input = input.to_str().unwrap();

    run(&["index", "--threads", "2", input]);
    assert!(PathBuf::from(format!("{input}.csi")).exists());
    assert_eq!(run(&["index", "--nrecords", input]).stdout, b"14\n");
    assert_eq!(
        run(&["index", "--stats", input]).stdout,
        b"chr1\t248956422\t9\nchr2\t242193529\t5\n"
    );
    assert_eq!(
        run(&["index", "--stats", "--all", input]).stdout,
        b"chr1\t248956422\t9\nchr2\t242193529\t5\n"
    );

    fs::remove_file(format!("{input}.csi")).unwrap();
    run(&["index", "--tbi", input]);
    assert!(PathBuf::from(format!("{input}.tbi")).exists());
    assert_eq!(run(&["index", "--nrecords", input]).stdout, b"14\n");
}

#[test]
fn builds_bcf_csi_and_reports_empty_contigs() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("calls.bcf");
    write_bcf(&fixture(), &input);
    let input = input.to_str().unwrap();

    run(&["index", input]);
    assert!(PathBuf::from(format!("{input}.csi")).exists());
    assert_eq!(run(&["index", "--nrecords", input]).stdout, b"14\n");
    assert_eq!(
        run(&["index", "--stats", input]).stdout,
        b"chr1\t248956422\t9\nchr2\t242193529\t5\n"
    );
    assert_eq!(
        run(&["index", "--stats", "--all", input]).stdout,
        b"chr1\t248956422\t9\nchr2\t242193529\t5\nchr3\t198295559\t.\n"
    );

    let result = Command::new(binary())
        .args(["index", "--tbi", input])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("BCF"));
}

#[test]
fn rejects_malformed_bcf_without_replacing_an_index() {
    let directory = tempfile::tempdir().unwrap();
    let valid = directory.path().join("valid.bcf");
    let input = directory.path().join("malformed.bcf");
    let output = directory.path().join("malformed.bcf.csi");
    write_bcf(&fixture(), &valid);
    corrupt_bcf_reference_type(&valid, &input);
    fs::write(&output, b"existing").unwrap();

    let result = Command::new(binary())
        .args([
            "index",
            "--force",
            "--output",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("REF/ALT"));
    assert_eq!(fs::read(output).unwrap(), b"existing");
}

#[test]
fn overwrite_is_explicit_and_transactional() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("bad.vcf.gz");
    let output = directory.path().join("bad.vcf.gz.csi");
    bgzip(
        b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\tbad\t.\tA\tT\t.\tPASS\t.\n",
        &input,
    );
    fs::write(&output, b"existing").unwrap();

    let result = Command::new(binary())
        .args([
            "index",
            "--force",
            "--output",
            output.to_str().unwrap(),
            input.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("VCF record 1"));
    assert_eq!(fs::read(output).unwrap(), b"existing");
}

#[test]
fn standard_input_requires_and_uses_named_output() {
    let input = fs::read(fixture()).unwrap();
    let directory = tempfile::tempdir().unwrap();
    let compressed = directory.path().join("input.vcf.gz");
    let output = directory.path().join("stdin.csi");
    bgzip(&input, &compressed);
    let compressed = fs::read(compressed).unwrap();

    let missing = Command::new(binary())
        .args(["index", "-"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(!missing.status.success());

    let mut child = Command::new(binary())
        .args(["index", "--output", output.to_str().unwrap(), "-"])
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&compressed).unwrap();
    assert!(child.wait().unwrap().success());
    assert!(output.exists());
}

#[test]
fn rejects_incomplete_and_unsorted_inputs() {
    let directory = tempfile::tempdir().unwrap();
    let incomplete = directory.path().join("incomplete.vcf.gz");
    let mut compressed = Vec::new();
    {
        let mut writer = bgzf::io::Writer::new(&mut compressed);
        writer.write_all(&fs::read(fixture()).unwrap()).unwrap();
        writer.try_finish().unwrap();
    }
    compressed.truncate(compressed.len() - 29);
    fs::write(&incomplete, compressed).unwrap();
    let result = Command::new(binary())
        .args(["index", incomplete.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("end-of-file marker"));

    let unsorted = directory.path().join("unsorted.vcf.gz");
    bgzip(
        b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t20\t.\tA\tT\t.\tPASS\t.\nchr1\t10\t.\tA\tT\t.\tPASS\t.\n",
        &unsorted,
    );
    let result = Command::new(binary())
        .args(["index", unsorted.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("positions are not sorted"));
    assert!(!PathBuf::from(format!("{}.csi", unsorted.display())).exists());
}
