use std::ffi::OsStr;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn oracle() -> PathBuf {
    std::env::var_os("RSOMICS_BCFTOOLS")
        .map(PathBuf::from)
        .expect("RSOMICS_BCFTOOLS must point to bcftools 1.24")
}

fn companion(name: &str) -> PathBuf {
    let oracle = oracle();
    let parent = oracle.parent().unwrap_or(Path::new("."));
    for path in [
        parent.join(name),
        parent.join("htslib-1.24").join(name),
        parent.join("htslib").join(name),
    ] {
        if path.is_file() {
            return path;
        }
    }
    PathBuf::from(name)
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bcftools-annotate")
        .join(name)
}

fn run(command: &mut Command) -> Output {
    command.output().unwrap()
}

fn require_success(output: Output, context: &str) -> Output {
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
    let output = require_success(
        run(Command::new(oracle()).arg("--version")),
        "bcftools version",
    );
    assert!(
        String::from_utf8_lossy(&output.stdout).starts_with("bcftools 1.24\n"),
        "{}",
        String::from_utf8_lossy(&output.stdout)
    );
}

fn compressed_vcf(source: &Path, destination: &Path) {
    require_success(
        run(Command::new(oracle()).args([
            OsStr::new("view"),
            OsStr::new("--no-version"),
            OsStr::new("-Oz"),
            OsStr::new("-o"),
            destination.as_os_str(),
            source.as_os_str(),
        ])),
        "compressing VCF fixture",
    );
    require_success(
        run(Command::new(oracle()).args([
            OsStr::new("index"),
            OsStr::new("-f"),
            destination.as_os_str(),
        ])),
        "indexing VCF fixture",
    );
}

fn compressed_tabular(source: &Path, destination: &Path, end_column: usize) {
    let output = File::create(destination).unwrap();
    let status = Command::new(companion("bgzip"))
        .args([OsStr::new("-c"), source.as_os_str()])
        .stdout(Stdio::from(output))
        .status()
        .unwrap();
    assert!(status.success(), "compressing tabular fixture");

    let end_column = end_column.to_string();
    let mut command = Command::new(companion("tabix"));
    command.args(["-f", "-s", "1", "-b", "2", "-e", &end_column]);
    require_success(run(command.arg(destination)), "indexing tabular fixture");
}

fn compressed_bed(source: &Path, destination: &Path) {
    let output = File::create(destination).unwrap();
    let status = Command::new(companion("bgzip"))
        .args([OsStr::new("-c"), source.as_os_str()])
        .stdout(Stdio::from(output))
        .status()
        .unwrap();
    assert!(status.success(), "compressing BED fixture");
    require_success(
        run(Command::new(companion("tabix")).args([
            OsStr::new("-f"),
            OsStr::new("-p"),
            OsStr::new("bed"),
            destination.as_os_str(),
        ])),
        "indexing BED fixture",
    );
}

#[derive(Debug)]
struct Decoded {
    header: noodles_vcf::Header,
    records: Vec<u8>,
}

fn decoded(path: &Path) -> Decoded {
    let output = require_success(
        run(Command::new(oracle()).args([
            OsStr::new("view"),
            OsStr::new("--no-version"),
            OsStr::new("-Ov"),
            path.as_os_str(),
        ])),
        "decoding annotation output",
    );
    let mut raw_header = Vec::new();
    let mut records = Vec::new();
    for line in output.stdout.split(|byte| *byte == b'\n') {
        if line.starts_with(b"##bcftools_") || line.is_empty() {
            continue;
        }
        if line.starts_with(b"#") {
            raw_header.extend_from_slice(line);
            raw_header.push(b'\n');
        } else {
            records.extend_from_slice(line);
            records.push(b'\n');
        }
    }
    Decoded {
        header: String::from_utf8(raw_header).unwrap().parse().unwrap(),
        records,
    }
}

fn assert_equivalent(
    directory: &Path,
    name: &str,
    target: &Path,
    ours: &[&str],
    oracle_arguments: &[&str],
) {
    let ours_output = directory.join(format!("{name}.ours.bcf"));
    let oracle_output = directory.join(format!("{name}.oracle.bcf"));
    let mut ours_command = Command::new(binary());
    ours_command
        .arg("annotate")
        .args(ours)
        .args(["-O", "b", "-o"])
        .arg(&ours_output)
        .arg(target);
    let mut oracle_command = Command::new(oracle());
    oracle_command
        .args(["annotate", "--no-version"])
        .args(oracle_arguments)
        .args(["-Ob", "-o"])
        .arg(&oracle_output)
        .arg(target);
    let ours_result = run(&mut ours_command);
    let oracle_result = run(&mut oracle_command);
    assert_eq!(
        ours_result.status.success(),
        oracle_result.status.success(),
        "{name}\nours={}\noracle={}",
        String::from_utf8_lossy(&ours_result.stderr),
        String::from_utf8_lossy(&oracle_result.stderr)
    );
    if ours_result.status.success() {
        let ours = decoded(&ours_output);
        let oracle = decoded(&oracle_output);
        assert_eq!(ours.header, oracle.header, "{name} header");
        assert_eq!(
            ours.records,
            oracle.records,
            "{name}\nours={}\noracle={}",
            String::from_utf8_lossy(&ours.records),
            String::from_utf8_lossy(&oracle.records)
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn official_missing_value_write_modes_match_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let target = directory.path().join("target.vcf.gz");
    let annotations = directory.path().join("annotations.tsv.gz");
    compressed_vcf(&fixture("annotate.missing.vcf"), &target);
    compressed_tabular(&fixture("annotate.missing.tab"), &annotations, 2);

    for (name, prefix) in [
        ("replace", ""),
        ("replace-missing", "."),
        ("add", "+"),
        ("add-missing", ".+"),
        ("append", "="),
        ("append-missing", ".="),
    ] {
        let columns =
            format!("CHROM,POS,REF,ALT,{prefix}INFO/TSTR,{prefix}INFO/TFLT,{prefix}INFO/TINT");
        let arguments = ["-a", annotations.to_str().unwrap(), "-c", columns.as_str()];
        assert_equivalent(directory.path(), name, &target, &arguments, &arguments);
    }

    let columns = "CHROM,POS,REF,ALT,-INFO/TSTR,-INFO/TFLT,-INFO/TINT";
    let oracle_result = run(Command::new(oracle()).args([
        OsStr::new("annotate"),
        OsStr::new("--no-version"),
        OsStr::new("-a"),
        annotations.as_os_str(),
        OsStr::new("-c"),
        OsStr::new(columns),
        target.as_os_str(),
    ]));
    assert!(!oracle_result.status.success());
    assert!(
        String::from_utf8_lossy(&oracle_result.stderr)
            .contains("the -INFO/TAG feature has not been implemented yet")
    );

    let output = directory.path().join("replace-existing.ours.bcf");
    let ours_result = run(Command::new(binary()).args([
        OsStr::new("annotate"),
        OsStr::new("-a"),
        annotations.as_os_str(),
        OsStr::new("-c"),
        OsStr::new(columns),
        OsStr::new("-O"),
        OsStr::new("b"),
        OsStr::new("-o"),
        output.as_os_str(),
        target.as_os_str(),
    ]));
    require_success(ours_result, "replace-existing extension");
    let records = decoded(&output)
        .records
        .split(|byte| *byte == b'\n')
        .filter(|line| line.starts_with(b"chr1\t"))
        .map(|line| String::from_utf8_lossy(line).into_owned())
        .collect::<Vec<_>>();
    assert_eq!(
        records,
        [
            "chr1\t10\t.\tA\tT\t.\t.\tTFLT=0.1;TINT=1;TSTR=old",
            "chr1\t20\t.\tA\tT\t.\t.\tTFLT=0.9;TINT=9;TSTR=new",
            "chr1\t30\t.\tA\tT\t.\t.\tTFLT=0.9;TINT=9;TSTR=new",
            "chr1\t40\t.\tA\tT\t.\t.\t.",
            "chr1\t50\t.\tA\tT\t.\t.\t.",
        ]
    );

    let character_source = directory.path().join("character.vcf");
    let character_table = directory.path().join("character.tsv");
    fs::write(
        &character_source,
        b"##fileformat=VCFv4.3\n\
##INFO=<ID=C,Number=1,Type=Character,Description=\"Character\">\n\
##contig=<ID=chr1,length=3>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\tC\t.\tPASS\tC=x\n\
chr1\t2\t.\tA\tC\t.\tPASS\tC=.\n\
chr1\t3\t.\tA\tC\t.\tPASS\t.\n",
    )
    .unwrap();
    fs::write(
        &character_table,
        b"chr1\t1\tA\tC\t.\nchr1\t2\tA\tC\tz\nchr1\t3\tA\tC\t.\n",
    )
    .unwrap();
    let character_target = directory.path().join("character.vcf.gz");
    let character_annotations = directory.path().join("character.tsv.gz");
    compressed_vcf(&character_source, &character_target);
    compressed_tabular(&character_table, &character_annotations, 2);
    for (name, column) in [
        ("character-replace", "INFO/C"),
        ("character-replace-missing", ".INFO/C"),
    ] {
        let columns = format!("CHROM,POS,REF,ALT,{column}");
        let arguments = [
            "-a",
            character_annotations.to_str().unwrap(),
            "-c",
            columns.as_str(),
        ];
        assert_equivalent(
            directory.path(),
            name,
            &character_target,
            &arguments,
            &arguments,
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn typed_variant_samples_and_allele_cardinality_match_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let target_source = directory.path().join("target.vcf");
    let annotation_source = directory.path().join("source.vcf");
    fs::write(
        &target_source,
        b"##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t10\told\tA\tC,G\t20\tPASS\tIA=1,2;IR=3,4,5;IG=0,1,2,3,4,5\tGT:DP\t1/2:3\t0/1:4\n",
    )
    .unwrap();
    fs::write(
        &annotation_source,
        b"##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA\n\
chr1\t10\tsource\tA\tC,G\t30\tPASS\tIA=10,20;IR=30,40,50;IG=0,1,2,3,4,5\tGT:DP\t1/2:40\t1|2:30\n",
    )
    .unwrap();
    let target = directory.path().join("target.vcf.gz");
    let annotations = directory.path().join("source.vcf.gz");
    compressed_vcf(&target_source, &target);
    compressed_vcf(&annotation_source, &annotations);
    let arguments = [
        "-a",
        annotations.to_str().unwrap(),
        "-c",
        "INFO/IA,INFO/IR,INFO/IG,FORMAT/GT,FORMAT/DP",
        "-s",
        "A",
    ];
    assert_equivalent(
        directory.path(),
        "typed-samples",
        &target,
        &arguments,
        &arguments,
    );

    let ar_target = directory.path().join("annotate.AR.vcf.gz");
    let ar_annotations = directory.path().join("annotate.AR.tab.gz");
    compressed_vcf(&fixture("annotate.AR.vcf"), &ar_target);
    compressed_tabular(&fixture("annotate.AR.tab"), &ar_annotations, 2);
    for (name, prefix) in [("official-ar-replace", ""), ("official-ar-add", "+")] {
        let columns = format!(
            "CHROM,POS,REF,ALT,{prefix}INFO/IA,{prefix}INFO/FA,{prefix}INFO/IR,{prefix}INFO/FR"
        );
        let arguments = [
            "-a",
            ar_annotations.to_str().unwrap(),
            "-c",
            columns.as_str(),
        ];
        assert_equivalent(directory.path(), name, &ar_target, &arguments, &arguments);
    }

    let official_target = directory.path().join("annotate2.vcf.gz");
    let official_annotations = directory.path().join("annots2.vcf.gz");
    compressed_vcf(&fixture("annotate2.vcf"), &official_target);
    compressed_vcf(&fixture("annots2.vcf"), &official_annotations);
    let arguments = [
        "-a",
        official_annotations.to_str().unwrap(),
        "-c",
        "ID,QUAL,+FILTER,FORMAT/GT",
        "-s",
        "A",
    ];
    assert_equivalent(
        directory.path(),
        "official-variant-sample",
        &official_target,
        &arguments,
        &arguments,
    );

    let incompatible = [
        "annotate",
        "-a",
        official_annotations.to_str().unwrap(),
        "-c",
        "+INFO",
        official_target.to_str().unwrap(),
    ];
    let ours_result = run(Command::new(binary()).args(incompatible));
    assert!(!ours_result.status.success());
    assert!(String::from_utf8_lossy(&ours_result.stderr).contains("incompatible schema"));
    let oracle_result = run(Command::new(oracle()).args(incompatible));
    assert!(oracle_result.status.success());
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn edits_expressions_and_indexed_regions_match_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let source = directory.path().join("target.vcf");
    let chromosome_map = directory.path().join("chromosomes.tsv");
    let annotation_map = directory.path().join("annotations.tsv");
    fs::write(
        &source,
        b"##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##INFO=<ID=OLD,Number=1,Type=String,Description=\"Old\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t20\tPASS\tDP=3;OLD=x\n\
chr1\t20\tb\tA\tT\t5\tPASS\tDP=5;OLD=y\n",
    )
    .unwrap();
    fs::write(&chromosome_map, "chr1\tchrZ\n").unwrap();
    fs::write(&annotation_map, "INFO/OLD\tNEW\n").unwrap();
    let target = directory.path().join("target.vcf.gz");
    compressed_vcf(&source, &target);

    assert_equivalent(
        directory.path(),
        "expression-edits",
        &target,
        &[
            "-x",
            "INFO/DP",
            "--rename-chromosomes",
            chromosome_map.to_str().unwrap(),
            "--rename-annotations",
            annotation_map.to_str().unwrap(),
            "-i",
            "QUAL>=10",
            "-k",
        ],
        &[
            "-x",
            "INFO/DP",
            "--rename-chrs",
            chromosome_map.to_str().unwrap(),
            "--rename-annots",
            annotation_map.to_str().unwrap(),
            "-i",
            "QUAL>=10",
            "-k",
        ],
    );
    assert_equivalent(
        directory.path(),
        "set-id-regions",
        &target,
        &["-I", "%CHROM\\_%POS", "-r", "chr1:10-10"],
        &["-I", "%CHROM\\_%POS", "-r", "chr1:10-10"],
    );

    let rename_target = directory.path().join("annotate21.vcf.gz");
    compressed_vcf(&fixture("annotate21.vcf"), &rename_target);
    assert_equivalent(
        directory.path(),
        "official-annotation-renames",
        &rename_target,
        &[
            "--rename-annotations",
            fixture("annotate21.txt").to_str().unwrap(),
        ],
        &[
            "--rename-annots",
            fixture("annotate21.txt").to_str().unwrap(),
        ],
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn interval_mark_overlap_and_pair_logic_match_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let target_source = directory.path().join("target.vcf");
    let table = directory.path().join("annotations.tsv");
    let variants_source = directory.path().join("source.vcf");
    fs::write(
        &target_source,
        b"##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##INFO=<ID=ANN,Number=1,Type=String,Description=\"Annotation\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\tid1\tA\t<CNV>\t.\tPASS\tEND=19\n\
chr1\t30\tid2\tA\tC\t.\tPASS\t.\n",
    )
    .unwrap();
    fs::write(&table, "chr1\t10\t15\tone\nchr1\t30\t30\ttwo\n").unwrap();
    fs::write(
        &variants_source,
        b"##fileformat=VCFv4.3\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t30\tid2\tA\tG,C\t.\tPASS\t.\n",
    )
    .unwrap();
    let target = directory.path().join("target.vcf.gz");
    let annotations = directory.path().join("annotations.tsv.gz");
    let variants = directory.path().join("source.vcf.gz");
    compressed_vcf(&target_source, &target);
    compressed_tabular(&table, &annotations, 3);
    compressed_vcf(&variants_source, &variants);

    let interval = [
        "-a",
        annotations.to_str().unwrap(),
        "-c",
        "CHROM,FROM,TO,INFO/ANN",
        "--min-overlap",
        "0.5:0.5",
        "-m",
        "+HIT",
    ];
    assert_equivalent(
        directory.path(),
        "interval-overlap",
        &target,
        &interval,
        &interval,
    );
    for logic in ["some", "exact", "id"] {
        let arguments = [
            "-a",
            variants.to_str().unwrap(),
            "-c",
            "+ID",
            "--pair-logic",
            logic,
        ];
        assert_equivalent(
            directory.path(),
            &format!("pair-{logic}"),
            &target,
            &arguments,
            &arguments,
        );
    }

    let mark_target = directory.path().join("annots-mark.vcf.gz");
    let mark_source = directory.path().join("annots-mark.bed");
    let mark_annotations = directory.path().join("annots-mark.bed.gz");
    fs::write(&mark_source, b"chr1\t10611\t10671").unwrap();
    compressed_vcf(&fixture("annots-mark.vcf"), &mark_target);
    compressed_bed(&mark_source, &mark_annotations);
    let arguments = [
        "-a",
        mark_annotations.to_str().unwrap(),
        "-c",
        "CHROM,FROM,TO",
        "-m",
        "+TAG",
    ];
    assert_equivalent(
        directory.path(),
        "official-mark-sites",
        &mark_target,
        &arguments,
        &arguments,
    );
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn malformed_inputs_and_rename_collisions_fail_like_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let source = directory.path().join("target.vcf");
    let renames = directory.path().join("renames.tsv");
    fs::write(
        &source,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=A,Number=1,Type=Integer,Description=\"A\">\n\
##INFO=<ID=B,Number=1,Type=Integer,Description=\"B\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tA\tC\t.\tPASS\tA=1;B=2\n",
    )
    .unwrap();
    fs::write(&renames, "INFO/A\tB\n").unwrap();
    let target = directory.path().join("target.vcf.gz");
    compressed_vcf(&source, &target);
    assert_equivalent(
        directory.path(),
        "rename-collision",
        &target,
        &["--rename-annotations", renames.to_str().unwrap()],
        &["--rename-annots", renames.to_str().unwrap()],
    );
}
