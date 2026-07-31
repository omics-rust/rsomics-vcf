use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_vcf::head::{self, Options};

const HEADER: &str = "##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##source=rsomics-test\n\
##contig=<ID=chr1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

const RECORD_1: &str = "chr1\t10\t.\tA\tC\t50.5\tPASS\tDP=7\tGT:DP\t0/1:20\t0/0:.\n";
const RECORD_2: &str = "chr1\t20\trs2\tG\tT\t.\t.\t.\tGT:DP\t1/1:5\t0/1:8\n";

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(path: &Path, options: Options) -> (String, head::Summary) {
    let mut output = Vec::new();
    let summary = head::write(path, options, &mut output).unwrap();
    (String::from_utf8(output).unwrap(), summary)
}

#[test]
fn plain_vcf_matches_the_stable_head_contract() {
    let path = fixture("head.vcf");
    let (default, summary) = run(&path, Options::default());
    assert_eq!(default, HEADER);
    assert_eq!(summary.header_lines, 8);
    assert_eq!(summary.records, 0);

    let (records, summary) = run(
        &path,
        Options {
            records: 2,
            ..Options::default()
        },
    );
    assert_eq!(records, format!("{HEADER}{RECORD_1}{RECORD_2}"));
    assert_eq!(summary.records, 2);
}

#[test]
fn samples_mode_starts_at_chrom() {
    let (output, summary) = run(
        &fixture("head.vcf"),
        Options {
            header_lines: Some(0),
            records: 1,
            from_chrom: true,
        },
    );
    assert_eq!(
        output,
        format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n{RECORD_1}")
    );
    assert_eq!(summary.header_lines, 1);
    assert_eq!(summary.records, 1);
}

#[test]
fn bgzf_vcf_matches_plain_vcf() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("head.vcf.gz");
    let mut writer = bgzf::io::Writer::new(File::create(&path).unwrap());
    writer
        .write_all(&std::fs::read(fixture("head.vcf")).unwrap())
        .unwrap();
    writer.try_finish().unwrap();

    let (output, summary) = run(
        &path,
        Options {
            records: 2,
            ..Options::default()
        },
    );
    assert_eq!(output, format!("{HEADER}{RECORD_1}{RECORD_2}"));
    assert_eq!(summary.records, 2);
}

#[test]
fn bcf_is_decoded_to_ordered_vcf_without_internal_indices() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("head.bcf");
    write_bcf(&fixture("head.vcf"), &path);

    let (output, summary) = run(
        &path,
        Options {
            records: 2,
            ..Options::default()
        },
    );
    assert!(output.starts_with("##fileformat=VCFv4.3\n"));
    assert!(output.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">\n"));
    assert!(!output.contains("IDX="));
    assert!(output.ends_with(&format!("{RECORD_1}{RECORD_2}")));
    assert_eq!(summary.records, 2);
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

#[test]
fn malformed_record_fails_instead_of_returning_partial_success() {
    let mut output = Vec::new();
    let error = head::write(
        &fixture("malformed.vcf"),
        Options {
            records: 1,
            ..Options::default()
        },
        &mut output,
    )
    .unwrap_err();
    assert!(error.to_string().contains("variant record 1"), "{error}");
}
