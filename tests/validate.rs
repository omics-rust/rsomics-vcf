use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_vcf::validate::{self, InputFormat, Options};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

#[test]
fn vcf_bgzf_and_bcf_share_the_validation_contract() {
    let directory = tempfile::tempdir().unwrap();
    let bgzf_vcf = directory.path().join("input.vcf.gz");
    let raw_bcf = directory.path().join("input.bcf");
    let bgzf_bcf = directory.path().join("input.bcf.gz");

    let mut vcf_writer = bgzf::io::Writer::new(File::create(&bgzf_vcf).unwrap());
    vcf_writer
        .write_all(&std::fs::read(fixture()).unwrap())
        .unwrap();
    vcf_writer.try_finish().unwrap();

    let mut raw_writer = bcf::io::Writer::from(File::create(&raw_bcf).unwrap());
    write_bcf(&fixture(), &mut raw_writer);
    raw_writer.get_mut().flush().unwrap();
    drop(raw_writer);

    let mut bgzf_writer = bcf::io::Writer::new(File::create(&bgzf_bcf).unwrap());
    write_bcf(&fixture(), &mut bgzf_writer);
    bgzf_writer.try_finish().unwrap();

    for (path, format) in [
        (fixture(), InputFormat::Vcf),
        (bgzf_vcf, InputFormat::Vcf),
        (raw_bcf, InputFormat::Bcf),
        (bgzf_bcf, InputFormat::Bcf),
    ] {
        let report = validate::check(&path, Options::default()).unwrap();
        assert!(report.is_valid(), "{}: {:?}", path.display(), report);
        assert_eq!(report.format, format);
        assert_eq!(report.version.as_deref(), Some("4.3"));
        assert_eq!(report.records, 3);
    }
}

#[test]
fn malformed_bgzf_vcf_reports_the_record_context() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("invalid.vcf.gz");
    let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
    writer
        .write_all(
            b"##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\tbad\t.\tA\tC\t.\tPASS\t.\n",
        )
        .unwrap();
    writer.try_finish().unwrap();

    let report = validate::check(&path, Options::default()).unwrap();
    assert!(!report.is_valid());
    assert!(report.diagnostics.iter().any(|item| {
        item.code == "record.pos" && item.line == Some(3) && item.field.as_deref() == Some("POS")
    }));
}

#[test]
fn truncated_bcf_reports_the_record_context() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("truncated.bcf");
    let mut writer = bcf::io::Writer::from(File::create(&path).unwrap());
    write_bcf(&fixture(), &mut writer);
    writer.get_mut().flush().unwrap();
    drop(writer);

    let length = std::fs::metadata(&path).unwrap().len();
    File::options()
        .write(true)
        .open(&path)
        .unwrap()
        .set_len(length - 1)
        .unwrap();

    let report = validate::check(&path, Options::default()).unwrap();
    assert!(!report.is_valid());
    assert!(
        report
            .diagnostics
            .iter()
            .any(|item| item.code == "bcf.record" && item.line.is_some()),
        "{:?}",
        report.diagnostics
    );
}

fn write_bcf<W>(input: &Path, writer: &mut bcf::io::Writer<W>)
where
    W: Write,
{
    let mut reader = vcf::io::Reader::new(BufReader::new(File::open(input).unwrap()));
    let header = reader.read_header().unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.records() {
        writer
            .write_variant_record(&header, &record.unwrap())
            .unwrap();
    }
}
