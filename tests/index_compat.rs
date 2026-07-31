use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn fixture() -> PathBuf {
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

fn query(path: &Path, region: &str) -> Vec<u8> {
    run(Command::new("bcftools")
        .args(["view", "-H", "-r", region])
        .arg(path))
    .stdout
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn vcf_indexes_match_bcftools_region_queries_and_stats() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let ours = directory.path().join("ours.vcf.gz");
    let oracle = directory.path().join("oracle.vcf.gz");
    bgzip(&fs::read(fixture()).unwrap(), &ours);
    fs::copy(&ours, &oracle).unwrap();

    run(Command::new(binary()).args(["index", ours.to_str().unwrap()]));
    run(Command::new("bcftools").args(["index", "-c"]).arg(&oracle));
    for region in [
        "chr1:20000-20000",
        "chr1:75000-75000",
        "chr1:151000-151000",
        "chr1:196000-196000",
        "chr1:205000-205000",
        "chr2:100000-100000",
        "chr2:130000-130000",
    ] {
        assert_eq!(query(&ours, region), query(&oracle, region), "{region}");
    }
    assert_eq!(
        run(Command::new(binary()).args(["index", "--stats", ours.to_str().unwrap()])).stdout,
        run(Command::new("bcftools")
            .args(["index", "--stats"])
            .arg(&oracle))
        .stdout
    );

    fs::remove_file(format!("{}.csi", ours.display())).unwrap();
    run(Command::new(binary()).args(["index", "--tbi", ours.to_str().unwrap()]));
    for region in ["chr1:151000-151000", "chr2:100000-100000"] {
        assert_eq!(query(&ours, region), query(&oracle, region), "{region}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn custom_csi_and_bcf_match_bcftools_region_queries() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let ours_vcf = directory.path().join("ours.vcf.gz");
    let oracle_vcf = directory.path().join("oracle.vcf.gz");
    bgzip(&fs::read(fixture()).unwrap(), &ours_vcf);
    fs::copy(&ours_vcf, &oracle_vcf).unwrap();

    run(Command::new(binary()).args(["index", "--min-shift", "18", ours_vcf.to_str().unwrap()]));
    run(Command::new("bcftools")
        .args(["index", "--min-shift", "18"])
        .arg(&oracle_vcf));
    for region in ["chr1:151000-151000", "chr2:100000-100000"] {
        assert_eq!(query(&ours_vcf, region), query(&oracle_vcf, region));
    }

    let ours_bcf = directory.path().join("ours.bcf");
    let oracle_bcf = directory.path().join("oracle.bcf");
    run(Command::new("bcftools")
        .args(["view", "-Ob", "-o", ours_bcf.to_str().unwrap()])
        .arg(&oracle_vcf));
    fs::copy(&ours_bcf, &oracle_bcf).unwrap();
    run(Command::new(binary()).args(["index", ours_bcf.to_str().unwrap()]));
    run(Command::new("bcftools")
        .args(["index", "-c"])
        .arg(&oracle_bcf));
    for region in ["chr1:151000-151000", "chr2:100000-100000"] {
        assert_eq!(query(&ours_bcf, region), query(&oracle_bcf, region));
    }
    assert_eq!(
        run(Command::new(binary()).args(["index", "--stats", "--all", ours_bcf.to_str().unwrap()]))
            .stdout,
        run(Command::new("bcftools")
            .args(["index", "--stats", "--all"])
            .arg(&oracle_bcf))
        .stdout
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn malformed_bcf_is_rejected_like_bcftools() {
    assert_bcftools_1_24();
    let directory = tempfile::tempdir().unwrap();
    let vcf = directory.path().join("input.vcf.gz");
    let valid = directory.path().join("valid.bcf");
    let malformed = directory.path().join("malformed.bcf");
    bgzip(&fs::read(fixture()).unwrap(), &vcf);
    run(Command::new("bcftools")
        .args(["view", "-Ob", "-o", valid.to_str().unwrap()])
        .arg(&vcf));
    corrupt_bcf_reference_type(&valid, &malformed);

    let ours = Command::new(binary())
        .args(["index", malformed.to_str().unwrap()])
        .output()
        .unwrap();
    let oracle = Command::new("bcftools")
        .args(["index", "-c"])
        .arg(&malformed)
        .output()
        .unwrap();
    assert!(!ours.status.success());
    assert!(!oracle.status.success());
}
