use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use noodles_bgzf as bgzf;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn run(arguments: &[&str]) -> std::process::Output {
    Command::new(binary()).args(arguments).output().unwrap()
}

fn fixtures(directory: &Path) -> (PathBuf, PathBuf) {
    let target = directory.join("target.vcf");
    let source = directory.join("source.vcf");
    fs::write(
        &target,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t10\told\tA\tC,G\t20\tPASS\tDP=3\tGT:DP\t1/2:3\t0/1:4\n\
chr1\t20\tkeep\tA\tT\t5\tPASS\tDP=5\tGT:DP\t0/1:5\t0/0:6\n",
    )
    .unwrap();
    fs::write(
        &source,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA\n\
chr1\t10\tsource\tA\tG,C\t30\tPASS\tAF=0.2,0.8\tGT:DP\t1/2:40\t2|1:30\n",
    )
    .unwrap();
    (target, source)
}

#[test]
fn lifecycle_rejects_incomplete_or_conflicting_actions() {
    let directory = tempfile::tempdir().unwrap();
    let (target, source) = fixtures(directory.path());
    let target = target.to_str().unwrap();
    let source = source.to_str().unwrap();
    let cases = [
        vec!["annotate", target],
        vec!["annotate", "--annotations", source, target],
        vec!["annotate", "--columns", "INFO/AF", target],
        vec!["annotate", "--mark-sites", "+HIT", target],
        vec![
            "annotate",
            "--remove",
            "ID",
            "--pair-logic",
            "exact",
            target,
        ],
        vec![
            "annotate",
            "--remove",
            "ID",
            "--min-overlap",
            "0.5:0.5",
            target,
        ],
        vec![
            "annotate",
            "--annotations",
            source,
            "--columns",
            "CHROM,POS",
            target,
        ],
        vec![
            "annotate",
            "--remove",
            "ID",
            "--include",
            "QUAL>10",
            "--exclude",
            "QUAL<10",
            target,
        ],
        vec![
            "annotate",
            "--remove",
            "ID",
            "--regions",
            "chr1:1-10",
            "--regions-file",
            target,
            target,
        ],
        vec![
            "annotate",
            "--annotations",
            source,
            "--columns",
            "FORMAT/DP",
            "--samples",
            "A",
            "--samples-file",
            target,
            target,
        ],
        vec![
            "annotate",
            "--annotations",
            source,
            "--columns",
            "INFO/AF",
            "--columns-file",
            target,
            target,
        ],
        vec![
            "annotate",
            "--header-line",
            "##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">",
            "--header-lines",
            target,
            target,
        ],
    ];
    for arguments in cases {
        let output = run(&arguments);
        assert!(!output.status.success(), "{arguments:?}");
        assert!(output.stdout.is_empty(), "{arguments:?}");
    }
}

#[test]
fn transfers_typed_fields_samples_and_reports_json_separately() {
    let directory = tempfile::tempdir().unwrap();
    let (target, source) = fixtures(directory.path());
    let output_path = directory.path().join("annotated.vcf");
    let output = run(&[
        "--json",
        "annotate",
        "--annotations",
        source.to_str().unwrap(),
        "--columns",
        "INFO/AF,FORMAT/GT,FORMAT/DP",
        "--samples",
        "A",
        "--output",
        output_path.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["command"], "annotate");
    assert_eq!(envelope["result"]["summary"]["read"], 2);
    assert_eq!(envelope["result"]["summary"]["annotated"], 1);
    assert_eq!(envelope["result"]["summary"]["unchanged"], 1);

    let output = fs::read_to_string(output_path).unwrap();
    assert!(
        output.contains("##INFO=<ID=AF,Number=A,Type=Float"),
        "{output}"
    );
    assert!(
        output.contains("AF=0.8,0.2\tGT:DP\t1|2:30\t0/1:4"),
        "{output}"
    );
}

#[test]
fn chromosome_renames_do_not_change_source_matching_coordinates() {
    let directory = tempfile::tempdir().unwrap();
    let (target, source) = fixtures(directory.path());
    let renames = directory.path().join("chromosomes.tsv");
    fs::write(&renames, "chr1\tchrZ\n").unwrap();
    let output = run(&[
        "annotate",
        "--annotations",
        source.to_str().unwrap(),
        "--columns",
        "INFO/AF",
        "--rename-chromosomes",
        renames.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("##contig=<ID=chrZ,length=100>"), "{output}");
    assert!(
        output.contains("chrZ\t10\told\tA\tC,G\t20\tPASS\tDP=3;AF=0.8,0.2"),
        "{output}"
    );
}

#[test]
fn tabular_columns_header_files_and_id_edits_share_the_public_stream() {
    let directory = tempfile::tempdir().unwrap();
    let (target, _) = fixtures(directory.path());
    let source = directory.path().join("depths.tsv");
    let columns = directory.path().join("columns.txt");
    let header = directory.path().join("header.txt");
    fs::write(&source, "chr1\t10\t9\n").unwrap();
    fs::write(&columns, "CHROM\nPOS\nINFO/X\n").unwrap();
    fs::write(
        &header,
        "##INFO=<ID=X,Number=1,Type=Integer,Description=\"X\">\n",
    )
    .unwrap();
    let output = run(&[
        "annotate",
        "--annotations",
        source.to_str().unwrap(),
        "--columns-file",
        columns.to_str().unwrap(),
        "--header-lines",
        header.to_str().unwrap(),
        "--set-id",
        r"%CHROM\_%POS",
        target.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(
        output.contains("##INFO=<ID=X,Number=1,Type=Integer"),
        "{output}"
    );
    assert!(
        output.contains("chr1\t10\tchr1_10\tA\tC,G\t20\tPASS\tDP=3;X=9"),
        "{output}"
    );
    assert!(output.contains("chr1\t20\tchr1_20\tA\tT"), "{output}");
}

#[test]
fn expressions_keep_sites_and_mark_both_match_states() {
    let directory = tempfile::tempdir().unwrap();
    let (target, source) = fixtures(directory.path());
    for (mark, tag, first_marked) in [("+HIT", "HIT", true), ("-MISS", "MISS", false)] {
        let output = run(&[
            "annotate",
            "--annotations",
            source.to_str().unwrap(),
            "--columns",
            "INFO/AF",
            "--mark-sites",
            mark,
            target.to_str().unwrap(),
        ]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = String::from_utf8(output.stdout).unwrap();
        assert!(output.contains(&format!("##INFO=<ID={tag}")));
        let records = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), 2, "{output}");
        assert_eq!(
            records[0].split('\t').nth(7).unwrap().contains(tag),
            first_marked
        );
        assert_eq!(
            records[1].split('\t').nth(7).unwrap().contains(tag),
            !first_marked
        );
    }

    for keep in [false, true] {
        let mut arguments = vec!["annotate", "--remove", "ID", "--include", "QUAL >= 10"];
        if keep {
            arguments.push("--keep-sites");
        }
        arguments.push(target.to_str().unwrap());
        let output = run(&arguments);
        assert!(output.status.success());
        let output = String::from_utf8(output.stdout).unwrap();
        let records = output
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>();
        assert_eq!(records.len(), usize::from(keep) + 1, "{output}");
        assert!(records[0].starts_with("chr1\t10\t.\t"), "{output}");
        if keep {
            assert!(records[1].starts_with("chr1\t20\tkeep\t"), "{output}");
        }
    }
}

#[test]
fn kept_expression_failures_retain_removed_fields_but_receive_global_renames() {
    let directory = tempfile::tempdir().unwrap();
    let (target, _) = fixtures(directory.path());
    let renames = directory.path().join("chromosomes.tsv");
    fs::write(&renames, "chr1\tchrZ\n").unwrap();
    let output = run(&[
        "annotate",
        "--remove",
        "INFO/DP",
        "--rename-chromosomes",
        renames.to_str().unwrap(),
        "--include",
        "QUAL >= 10",
        "--keep-sites",
        target.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let output = String::from_utf8(output.stdout).unwrap();
    assert!(output.contains("##contig=<ID=chrZ,length=100>"), "{output}");
    assert!(
        output.contains("##INFO=<ID=DP,Number=1,Type=Integer"),
        "{output}"
    );
    assert!(
        output.contains("chrZ\t10\told\tA\tC,G\t20\tPASS\t."),
        "{output}"
    );
    assert!(
        output.contains("chrZ\t20\tkeep\tA\tT\t5\tPASS\tDP=5"),
        "{output}"
    );
}

#[test]
fn supports_stdin_regions_bgzf_workers_and_every_output_encoding() {
    let directory = tempfile::tempdir().unwrap();
    let (target, _) = fixtures(directory.path());

    let mut child = Command::new(binary())
        .args(["annotate", "--remove", "ID", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(&target).unwrap())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .contains("chr1\t10\t.\t")
    );

    let compressed = directory.path().join("target.vcf.gz");
    let viewed = run(&[
        "view",
        "--output-type",
        "z",
        "--output",
        compressed.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(
        viewed.status.success(),
        "{}",
        String::from_utf8_lossy(&viewed.stderr)
    );
    let indexed = run(&["index", compressed.to_str().unwrap()]);
    assert!(
        indexed.status.success(),
        "{}",
        String::from_utf8_lossy(&indexed.stderr)
    );
    let selected = run(&[
        "annotate",
        "--remove",
        "ID",
        "--regions",
        "chr1:20-20",
        compressed.to_str().unwrap(),
    ]);
    assert!(
        selected.status.success(),
        "{}",
        String::from_utf8_lossy(&selected.stderr)
    );
    let selected = String::from_utf8(selected.stdout).unwrap();
    assert!(!selected.contains("chr1\t10\t"), "{selected}");
    assert!(selected.contains("chr1\t20\t.\t"), "{selected}");

    for (kind, extension) in [
        ("v", "vcf"),
        ("z", "vcf.gz"),
        ("b", "bcf"),
        ("u", "raw.bcf"),
    ] {
        let encoded = directory.path().join(format!("output.{extension}"));
        let decoded = directory.path().join(format!("output.{extension}.vcf"));
        let mut arguments = vec![
            "annotate",
            "--remove",
            "ID",
            "--output-type",
            kind,
            "--output",
            encoded.to_str().unwrap(),
        ];
        if matches!(kind, "z" | "b") {
            arguments.extend(["--threads", "2"]);
        }
        arguments.push(target.to_str().unwrap());
        let annotated = run(&arguments);
        assert!(
            annotated.status.success(),
            "{kind}: {}",
            String::from_utf8_lossy(&annotated.stderr)
        );
        let viewed = run(&[
            "view",
            "--output",
            decoded.to_str().unwrap(),
            encoded.to_str().unwrap(),
        ]);
        assert!(
            viewed.status.success(),
            "{kind}: {}",
            String::from_utf8_lossy(&viewed.stderr)
        );
        assert_eq!(
            fs::read_to_string(decoded)
                .unwrap()
                .lines()
                .filter(|line| !line.starts_with('#'))
                .count(),
            2
        );
    }
}

#[test]
fn named_output_is_transactional_and_rejects_input_aliases() {
    let directory = tempfile::tempdir().unwrap();
    let (target, source) = fixtures(directory.path());
    let malformed = directory.path().join("malformed.vcf");
    let mut content = fs::read(&target).unwrap();
    content.extend_from_slice(b"chr1\t30\tbad\n");
    fs::write(&malformed, content).unwrap();
    let output_path = directory.path().join("output.vcf");
    fs::write(&output_path, b"existing").unwrap();
    let failed = run(&[
        "annotate",
        "--remove",
        "ID",
        "--output",
        output_path.to_str().unwrap(),
        malformed.to_str().unwrap(),
    ]);
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output_path).unwrap(), b"existing");

    for arguments in [
        vec![
            "annotate",
            "--remove",
            "ID",
            "--output",
            target.to_str().unwrap(),
            target.to_str().unwrap(),
        ],
        vec![
            "annotate",
            "--annotations",
            source.to_str().unwrap(),
            "--columns",
            "INFO/AF",
            "--output",
            source.to_str().unwrap(),
            target.to_str().unwrap(),
        ],
    ] {
        let rejected = run(&arguments);
        assert!(!rejected.status.success(), "{arguments:?}");
        assert!(rejected.stdout.is_empty(), "{arguments:?}");
    }
}

#[test]
fn annotate_help_is_complete_and_bounded() {
    let output = run(&["annotate", "--help"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let help = String::from_utf8(output.stdout).unwrap();
    for option in [
        "-a, --annotations",
        "-c, --columns",
        "-C, --columns-file",
        "-H, --header-line",
        "--header-lines",
        "-I, --set-id",
        "-x, --remove",
        "--rename-chromosomes",
        "--rename-annotations",
        "-i, --include",
        "-e, --exclude",
        "-k, --keep-sites",
        "-m, --mark-sites",
        "--min-overlap",
        "--pair-logic",
        "-s, --samples",
        "-S, --samples-file",
        "-r, --regions",
        "-R, --regions-file",
        "--regions-overlap",
        "-O, --output-type",
        "-o, --output",
        "--threads",
    ] {
        assert!(help.contains(option), "missing {option}: {help}");
    }
    for excluded in [
        "--force",
        "--merge-logic",
        "--single-overlaps",
        "--write-index",
        "--compression-level",
        "--no-version",
    ] {
        assert!(!help.contains(excluded), "unexpected {excluded}: {help}");
    }
}

#[test]
fn parallel_output_finishes_before_the_named_file_is_committed() {
    let directory = tempfile::tempdir().unwrap();
    let (target, _) = fixtures(directory.path());
    let output_path = directory.path().join("annotated.vcf.gz");
    let output = run(&[
        "annotate",
        "--remove",
        "ID",
        "--output-type",
        "z",
        "--threads",
        "2",
        "--output",
        output_path.to_str().unwrap(),
        target.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut decoded = String::new();
    bgzf::io::Reader::new(File::open(output_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert!(decoded.contains("chr1\t10\t.\t"), "{decoded}");
}
