use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};
use rsomics_vcf::query::{self, HeaderMode, Options};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/head.vcf")
}

fn query_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/query.vcf")
}

fn run(input: &Path, format: &str, samples: Option<Vec<String>>, header: HeaderMode) -> String {
    let mut output = Vec::new();
    query::write(
        input,
        &Options {
            format: format.to_owned(),
            samples,
            header,
            automatic_newline: true,
        },
        &mut output,
    )
    .unwrap();
    String::from_utf8(output).unwrap()
}

#[test]
fn renders_site_and_sample_fields() {
    assert_eq!(
        run(
            &fixture(),
            r"%CHROM\t%POS\t%INFO/DP[\t%SAMPLE=%GT]\n",
            None,
            HeaderMode::None,
        ),
        "chr1\t10\t7\tS1=0/1\tS2=0/0\n\
chr1\t20\t.\tS1=1/1\tS2=0/1\n\
chr1\t30\t3\tS1=./.\tS2=0/1\n"
    );
}

#[test]
fn sample_selection_applies_to_format_and_line() {
    assert_eq!(
        run(
            &fixture(),
            r"%FORMAT\n%LINE\n",
            Some(vec!["S2".to_owned()]),
            HeaderMode::None,
        ),
        "GT:DP\t0/0:.\n\
chr1\t10\t.\tA\tC\t50.5\tPASS\tDP=7\tGT:DP\t0/0:.\n\
GT:DP\t0/1:8\n\
chr1\t20\trs2\tG\tT\t.\t.\t.\tGT:DP\t0/1:8\n\
GT:DP\t0/1:3\n\
chr1\t30\t.\tC\tG\t9\tPASS\tDP=3\tGT:DP\t0/1:3\n"
    );
}

#[test]
fn numbered_header_and_newline_rule_match_query_contract() {
    assert_eq!(
        run(
            &fixture(),
            r"%CHROM\t%POS[\t%GT]",
            Some(vec!["S1".to_owned()]),
            HeaderMode::Indexed,
        ),
        "#[1]CHROM\t[2]POS\t[3]S1:GT\n\
chr1\t10\t0/1\n\
chr1\t20\t1/1\n\
chr1\t30\t./.\n"
    );
    assert_eq!(
        run(&fixture(), r"%CHROM\n%POS", None, HeaderMode::None,),
        "chr1\n10chr1\n20chr1\n30"
    );
}

#[test]
fn sample_fields_shadow_fixed_columns_until_forced() {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        input,
        "##fileformat=VCFv4.3\n\
##contig=<ID=1>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n\
##INFO=<ID=A.B,Number=1,Type=Integer,Description=\"Dot tag\">\n\
##FORMAT=<ID=ID,Number=1,Type=String,Description=\"ID\">\n\
##FORMAT=<ID=ALT,Number=1,Type=String,Description=\"ALT\">\n\
##FORMAT=<ID=POS0,Number=1,Type=String,Description=\"POS0\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
1\t5\tsite\tA\tC\t.\tPASS\tAF=0.25,0.75;A.B=4\tID:ALT:POS0:AD\tformat-id:format-alt:format-pos:7,3"
    )
    .unwrap();

    assert_eq!(
        run(
            input.path(),
            r"[%ID\t%/ID\t%ALT\t%/ALT\t%POS0\t%/POS0\t%AD{1}]\t%INFO/AF{0}\t%A.B\n",
            None,
            HeaderMode::None,
        ),
        "format-id\tsite\tformat-alt\tC\tformat-pos\t4\t3\t0.25\t4\n"
    );
}

#[test]
fn singleton_info_subscripts_reuse_the_only_value() {
    assert_eq!(
        run(
            &query_fixture(),
            r"%INFO/DP{9}\t%INFO/AF{9}[\t%DP{9}]\n",
            None,
            HeaderMode::None,
        ),
        "30\t0.25\t.\t.\t.\n\
14\t.\t.\t.\t.\n\
150\t.\t.\t.\t.\n"
    );
}

#[test]
fn exclusion_keeps_remaining_samples_in_header_order() {
    assert_eq!(
        run(
            &fixture(),
            r"[%SAMPLE=%GT\n]",
            Some(vec!["^S1".to_owned()]),
            HeaderMode::None,
        ),
        "S2=0/0\nS2=0/1\nS2=0/1\n"
    );
    assert_eq!(
        run(
            &fixture(),
            r"%FORMAT\n%LINE\n",
            Some(vec!["^S1".to_owned(), "S2".to_owned()]),
            HeaderMode::None,
        ),
        "\t.\n\
chr1\t10\t.\tA\tC\t50.5\tPASS\tDP=7\n\
\t.\n\
chr1\t20\trs2\tG\tT\t.\t.\t.\n\
\t.\n\
chr1\t30\t.\tC\tG\t9\tPASS\tDP=3\n"
    );
}

#[test]
fn bgzf_and_bcf_match_plain_vcf() {
    let directory = tempfile::tempdir().unwrap();
    let bgzf_path = directory.path().join("query.vcf.gz");
    let bcf_path = directory.path().join("query.bcf");

    let mut bgzf_writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
    bgzf_writer
        .write_all(&std::fs::read(fixture()).unwrap())
        .unwrap();
    bgzf_writer.try_finish().unwrap();

    let mut reader = vcf::io::Reader::new(BufReader::new(File::open(fixture()).unwrap()));
    let header = reader.read_header().unwrap();
    let mut bcf_writer = bcf::io::Writer::new(File::create(&bcf_path).unwrap());
    bcf_writer.write_header(&header).unwrap();
    for record in reader.records() {
        bcf_writer
            .write_variant_record(&header, &record.unwrap())
            .unwrap();
    }
    bcf_writer.try_finish().unwrap();

    let format = r"%CHROM\t%POS\t%TYPE[\t%SAMPLE=%IUPACGT]\n";
    let expected = run(&fixture(), format, None, HeaderMode::None);
    assert_eq!(run(&bgzf_path, format, None, HeaderMode::None), expected);
    assert_eq!(run(&bcf_path, format, None, HeaderMode::None), expected);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn supported_query_contract_matches_bcftools_1_24() {
    let version = Command::new("bcftools").arg("--version").output().unwrap();
    assert!(
        String::from_utf8(version.stdout)
            .unwrap()
            .starts_with("bcftools 1.24\n")
    );

    let directory = tempfile::tempdir().unwrap();
    let bgzf_path = directory.path().join("query.vcf.gz");
    let bcf_path = directory.path().join("query.bcf");
    for (output_type, output) in [("-Oz", &bgzf_path), ("-Ob", &bcf_path)] {
        let result = Command::new("bcftools")
            .args(["view", "--no-version", output_type, "-o"])
            .arg(output)
            .arg(query_fixture())
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    let cases = [
        (
            r"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\n",
            None,
            HeaderMode::None,
        ),
        (
            r"%POS0\t%END\t%TYPE\t%INFO\t%INFO/DP\n",
            None,
            HeaderMode::None,
        ),
        (
            r"%INFO/AF\t%INFO/AF{0}\t%INFO/AF{1}\t%INFO/DB\n",
            None,
            HeaderMode::None,
        ),
        (
            r"%FORMAT\n%LINE\n",
            Some(vec!["S2".to_owned()]),
            HeaderMode::None,
        ),
        (
            r"%CHROM\t%POS[\t%SAMPLE=%GT:%TGT:%IUPACGT:%DP]\n",
            None,
            HeaderMode::Indexed,
        ),
        (
            r"%CHROM\t%POS[\t%GT]",
            Some(vec!["S1".to_owned()]),
            HeaderMode::Plain,
        ),
        (r"%CHROM\n%POS", None, HeaderMode::None),
    ];

    for input in [query_fixture(), bgzf_path, bcf_path.clone()] {
        for (format, samples, header_mode) in &cases {
            let expected = run(&input, format, samples.clone(), *header_mode);
            let mut command = Command::new("bcftools");
            command.arg("query").arg("-f").arg(format);
            if let Some(samples) = samples {
                command.arg("-s").arg(samples.join(","));
            }
            match header_mode {
                HeaderMode::None => {}
                HeaderMode::Indexed => {
                    command.arg("-H");
                }
                HeaderMode::Plain => {
                    command.arg("-HH");
                }
            }
            command.arg(&input);
            let output = command.output().unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(
                expected.as_bytes(),
                output.stdout,
                "input={} format={format:?} samples={samples:?} header={header_mode:?}",
                input.display()
            );
        }
    }

    let format = r"%FORMAT\n%LINE\n";
    let samples = Some(vec!["^S1".to_owned(), "S2".to_owned(), "S3".to_owned()]);
    let expected = run(&query_fixture(), format, samples.clone(), HeaderMode::None);
    assert_eq!(run(&bcf_path, format, samples, HeaderMode::None), expected);
    let vcf_output = Command::new("bcftools")
        .args(["query", "-f", format, "-s", "^S1,S2,S3"])
        .arg(query_fixture())
        .output()
        .unwrap();
    assert!(vcf_output.status.success());
    assert_eq!(vcf_output.stdout, expected.as_bytes());
    let bcf_output = Command::new("bcftools")
        .args(["query", "-f", format, "-s", "^S1,S2,S3"])
        .arg(&bcf_path)
        .output()
        .unwrap();
    assert!(bcf_output.status.success());
    assert_ne!(bcf_output.stdout, expected.as_bytes());
}
