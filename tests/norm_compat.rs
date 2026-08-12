#![cfg(feature = "norm-preview")]

use std::fmt::Write as _;
use std::fs;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::{Command, Output};

use noodles_bgzf as bgzf;

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

fn body(output: Output) -> Vec<u8> {
    output
        .stdout
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.starts_with(b"#") && !line.is_empty())
        .flat_map(|line| line.iter().copied().chain([b'\n']))
        .collect()
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_realignments_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nARAAAACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tA\tC\t.\tPASS\t.\n\
chr1\t4\t.\tA\tAA\t.\tPASS\t.\n\
chr1\t4\t.\tAA\tA\t.\tPASS\t.\n\
chr1\t9\t.\tTAC\tTAG\t.\tPASS\t.\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn typed_multiallelic_split_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=FA,Number=A,Type=Integer,Description=\"A\">\n\
##FORMAT=<ID=FR,Number=R,Type=Integer,Description=\"R\">\n\
##FORMAT=<ID=FG,Number=G,Type=Integer,Description=\"G\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\tIA=10,20;IR=5,3,2;IG=0,10,20,30,40,50\tGT:FA:FR:FG\t1/2:11,22:7,4,3:0,10,20,30,40,50\t2:33,44:8,5,6:0,10,20\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args(["norm", "--split-multiallelic", input.to_str().unwrap()])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--split-multiallelic",
        "--split-overlaps",
        "missing",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        "--multi-overlaps",
        ".",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn expression_selected_split_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\tselected\tA\tC,G\t.\tPASS\tDP=20\n\
chr1\t20\tunchanged\tA\tC,G\t.\tPASS\tDP=5\n",
    )
    .unwrap();

    for (argument, expression) in [("--include", "INFO/DP>=10"), ("--exclude", "INFO/DP<10")] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--split-multiallelic",
            argument,
            expression,
            input.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "-m",
            "-any",
            argument,
            expression,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{argument}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn streaming_target_modes_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tAT\tA,G\t.\tPASS\t.\n\
chr1\t20\tb\tA\tC,G\t.\tPASS\t.\n",
    )
    .unwrap();

    for (mode, oracle_mode) in [("pos", "0"), ("record", "1"), ("variant", "2")] {
        for targets in ["chr1:11", "^chr1:11"] {
            let ours = body(run(Command::new(PathBuf::from(env!(
                "CARGO_BIN_EXE_rsomics-vcf"
            )))
            .args([
                "norm",
                "--split-multiallelic",
                "--targets",
                targets,
                "--targets-overlap",
                mode,
                input.to_str().unwrap(),
            ])));
            let oracle = body(run(Command::new("bcftools").args([
                "norm",
                "--no-version",
                "-m",
                "-any",
                "--targets",
                targets,
                "--targets-overlap",
                oracle_mode,
                input.to_str().unwrap(),
            ])));
            assert_eq!(ours, oracle, "{mode} {targets}");
        }
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn typed_biallelic_join_any_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FILTER=<ID=q10,Description=\"Low quality\">\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Float,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=String,Description=\"G\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=FA,Number=A,Type=Integer,Description=\"A\">\n\
##FORMAT=<ID=FR,Number=R,Type=Float,Description=\"R\">\n\
##FORMAT=<ID=FG,Number=G,Type=String,Description=\"G\">\n\
##FORMAT=<ID=NG,Number=G,Type=Integer,Description=\"G\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\tc\tA\tC\t10\tq10\tIA=10;IR=5.0,3.0;IG=r,rc,cc;DP=7\tGT:FA:FR:FG\t0/1:11:7.0,4.0:r,rc,cc\t./1:12:8.0,5.0:R,RC,CC\n\
chr1\t10\tg\tA\tG\t20\tPASS\tIA=20;IR=6.0,2.0;IG=s,sg,gg;DP=9\tGT:FA:FR:FG\t1/0:21:9.0,3.0:s,sg,gg\t0/1:22:10.0,6.0:S,SG,GG\n\
chr1\t20\tt\tC\tT\t30\tPASS\tIA=30;IR=7.0,1.0;IG=x,xt,tt;DP=11\tGT:FA:FR:FG\t0/1:31:11.0,2.0:x,xt,tt\t0/0:32:12.0,1.0:X,XT,TT\n\
chr1\t30\tlo\ta\tc\t.\tPASS\t.\tGT\t0/1\t0/0\n\
chr1\t30\tup\tA\tG\t.\tPASS\t.\tGT\t1/0\t0/1\n\
chr1\t40\tb\tA\tA]chr1:80]\t.\tPASS\t.\tGT\t0/1\t0/0\n\
chr1\t40\td\tAT\tA\t.\tPASS\t.\tGT\t1/0\t0/1\n\
chr1\t50\th\tA\tC\t.\tPASS\t.\tGT:NG:FG\t1:10,11:h0,h1\t0:30,31:H0,H1\n\
chr1\t50\td\tA\tG\t.\tPASS\t.\tGT:NG:FG\t0/1:20,21,22:d0,d1,d2\t0/1:40,41,42:D0,D1,D2\n\
chr1\t60\td\tA\tC\t.\tPASS\t.\tGT:NG:FG\t0/1:10,11,12:d0,d1,d2\t0/0:30,31,32:D0,D1,D2\n\
chr1\t60\th\tA\tG\t.\tPASS\t.\tGT:NG:FG\t1:20,21:h0,h1\t0:40,41:H0,H1\n",
    )
    .unwrap();

    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "+any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(
        oracle,
        b"chr1\t10\tc;g\tA\tC,G\t20\tq10\tIA=10,20;IR=6,3,2;IG=r,rc,cc,sg,.,gg;DP=7\tGT:FA:FR:FG\t2/1:11,21:9,4,3:r,rc,cc,sg,.,gg\t2/1:12,22:10,5,6:R,RC,CC,SG,.,GG\n\
chr1\t20\tt\tC\tT\t30\tPASS\tIA=30;IR=7,1;IG=x,xt,tt;DP=11\tGT:FA:FR:FG\t0/1:31:11,2:x,xt,tt\t0/0:32:12,1:X,XT,TT\n\
chr1\t30\tlo;up\tA\tC,G\t.\tPASS\t.\tGT\t2/1\t0/2\n\
chr1\t40\tb;d\tAT\tA]chr1:80]T,A\t.\tPASS\t.\tGT\t2/1\t0/2\n\
chr1\t50\th;d\tA\tC,G\t.\tPASS\t.\tGT:NG:FG\t1/2:20,11,.,21,.,22:h0,h1,d1\t0/2:40,31,.,41,.,42:H0,H1,D1\n\
chr1\t60\td;h\tA\tC,G\t.\tPASS\t.\tGT:NG:FG\t2/1:20,11,21,.,.,.:d0,d1,d2,h1,.,.\t0/0:40,31,41,.,.,.:D0,D1,D2,H1,.,.\n"
    );
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--join-multiallelic",
        "any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);

    for output_type in ["v", "z", "b", "u"] {
        let ours = directory.path().join(format!("joined-ours.{output_type}"));
        let oracle = directory
            .path()
            .join(format!("joined-oracle.{output_type}"));
        run(
            Command::new(PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))).args([
                "norm",
                "--join-multiallelic",
                "any",
                "-O",
                output_type,
                "-o",
                ours.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
        );
        run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "-m",
            "+any",
            "-O",
            output_type,
            "-o",
            oracle.to_str().unwrap(),
            input.to_str().unwrap(),
        ]));
        let ours = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            ours.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            oracle.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{output_type}");
    }

    for input_type in ["v", "z", "b", "u"] {
        let encoded = directory.path().join(format!("joined-input.{input_type}"));
        run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-O",
            input_type,
            "-o",
            encoded.to_str().unwrap(),
            input.to_str().unwrap(),
        ]));
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--join-multiallelic",
            "any",
            encoded.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "-m",
            "+any",
            encoded.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{input_type}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn strict_filter_join_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FILTER=<ID=q10,Description=\"q10\">\n\
##FILTER=<ID=s20,Description=\"s20\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t.\tq10\t.\n\
chr1\t10\tb\tA\tG\t.\tPASS\t.\n\
chr1\t20\ta\tA\tC\t.\tPASS\t.\n\
chr1\t20\tb\tA\tG\t.\tq10\t.\n\
chr1\t30\ta\tA\tC\t.\tq10\t.\n\
chr1\t30\tb\tA\tG\t.\tPASS\t.\n\
chr1\t30\tc\tA\tT\t.\ts20\t.\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--join-multiallelic",
        "any",
        "--strict-filter",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "+any",
        "--strict-filter",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn expression_selected_join_matches_bcftools_records() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ta\tA\tC\t.\tPASS\tDP=10\n\
chr1\t10\tc\tA\tT\t.\tPASS\tDP=20\n\
chr1\t20\tb\tA\tG\t.\tPASS\tDP=0\n\
chr1\t30\td\tA\tC\t.\tPASS\tDP=10\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--join-multiallelic",
        "any",
        "--include",
        "INFO/DP>=10",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "+any",
        "--include",
        "INFO/DP>=10",
        input.to_str().unwrap(),
    ])));
    let mut ours_records = ours.split(|byte| *byte == b'\n').collect::<Vec<_>>();
    let mut oracle_records = oracle.split(|byte| *byte == b'\n').collect::<Vec<_>>();
    ours_records.sort_unstable();
    oracle_records.sort_unstable();
    assert_eq!(ours_records, oracle_records);
    let positions = ours
        .split(|byte| *byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            std::str::from_utf8(line)
                .unwrap()
                .split('\t')
                .nth(1)
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["10", "20", "30"]);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn classified_join_modes_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ti1\tA\tAT\t.\tPASS\t.\n\
chr1\t10\ts1\tA\tC\t.\tPASS\t.\n\
chr1\t10\to1\tA\t<DEL>\t.\tPASS\t.\n\
chr1\t10\tr1\tA\t.\t.\tPASS\t.\n\
chr1\t10\tm1\tAC\tGT\t.\tPASS\t.\n\
chr1\t10\tb1\tA\tA]chr1:80]\t.\tPASS\t.\n\
chr1\t10\tstar\tA\t*\t.\tPASS\t.\n\
chr1\t10\ts2\tA\tG\t.\tPASS\t.\n\
chr1\t10\ti2\tA\tAG\t.\tPASS\t.\n\
chr1\t10\tm2\tAC\tTT\t.\tPASS\t.\n\
chr1\t10\to2\tA\t<DUP>\t.\tPASS\t.\n\
chr1\t10\tr2\tA\t.\t.\tPASS\t.\n",
    )
    .unwrap();

    for (mode, expected_ids) in [
        (
            "snps",
            vec![
                "r1;r2;s1;s2",
                "m1",
                "m2",
                "i1",
                "i2",
                "o1",
                "o2",
                "b1",
                "star",
            ],
        ),
        (
            "indels",
            vec![
                "r1;r2", "s1", "s2", "m1", "m2", "i1;i2", "o1", "o2", "b1", "star",
            ],
        ),
        (
            "both",
            vec!["r1;r2;s1;s2", "m1;m2", "i1;i2", "o1;o2", "b1", "star"],
        ),
    ] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args(["norm", "--join-multiallelic", mode, input.to_str().unwrap()])));
        let join = format!("+{mode}");
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "-m",
            &join,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{mode}");
        let ids = ours
            .split(|byte| *byte == b'\n')
            .filter(|line| !line.is_empty())
            .map(|line| {
                std::str::from_utf8(line)
                    .unwrap()
                    .split('\t')
                    .nth(2)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        assert_eq!(ids, expected_ids, "{mode}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn split_ad_sum_preservation_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t1/2:10,3,2\t0/2:10,.,2\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--split-multiallelic",
        "--keep-sum",
        "AD",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        "--keep-sum",
        "AD",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn split_and_reference_realign_compose_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nAAAAAACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t13\t6\t13\t14\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=13>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tA\tAA,AAA\t.\tPASS\tAF=0.25,0.5\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--split-multiallelic",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "-m",
        "-any",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_warn_and_skip_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tT\tA\t.\tPASS\t.\n\
chr1\t3\t.\tG\tA\t.\tPASS\t.\n",
    )
    .unwrap();

    for (ours_policy, oracle_policy) in [("warn", "w"), ("skip", "x")] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            ours_policy,
            input.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            oracle_policy,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{ours_policy}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_fix_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t8\t6\t8\t9\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=8>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t2\tswap\tT\tA,C,G\t.\tPASS\tAC=4,5,6\tGT\t0/2\t2|0\t./2\t0\n\
chr1\t3\tset-snp\tT\tA,C\t.\tPASS\tAC=7,8\tGT\t0/1\t1|2\t./0\t2\n\
chr1\t5\tset-mnv\tAT\tTT,AG,<DEL>\t.\tPASS\tAC=1,2,3\tGT\t0/1\t1|2\t./0\t3\n",
    )
    .unwrap();

    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "s",
        input.to_str().unwrap(),
    ])));
    assert_eq!(
        oracle,
        b"chr1\t2\tswap\tC\tA,T,G\t.\tPASS\tAC=4,3,6\tGT\t2/0\t0|2\t./0\t2\n\
chr1\t3\tset-snp\tG\tA,C\t.\tPASS\tAC=7,8\tGT\t0/1\t1|2\t./0\t2\n\
chr1\t5\tset-mnv\tAC\tTC,AG,<DEL>\t.\tPASS\tAC=1,2,3\tGT\t0/1\t1|2\t./0\t3\n"
    );
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "fix",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);

    for output_type in ["v", "z", "b", "u"] {
        let ours = directory.path().join(format!("ours.{output_type}"));
        let oracle = directory.path().join(format!("oracle.{output_type}"));
        run(
            Command::new(PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))).args([
                "norm",
                "--fasta-ref",
                reference.to_str().unwrap(),
                "--check-ref",
                "fix",
                "-O",
                output_type,
                "-o",
                ours.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
        );
        run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "s",
            "-O",
            output_type,
            "-o",
            oracle.to_str().unwrap(),
            input.to_str().unwrap(),
        ]));
        let ours = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            ours.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            oracle.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{output_type}");
    }

    for input_type in ["v", "z", "b", "u"] {
        let encoded = directory.path().join(format!("input.{input_type}"));
        run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-O",
            input_type,
            "-o",
            encoded.to_str().unwrap(),
            input.to_str().unwrap(),
        ]));
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "fix",
            encoded.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--fasta-ref",
            reference.to_str().unwrap(),
            "--check-ref",
            "s",
            encoded.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{input_type}");
    }

    let compressed_reference = directory.path().join("reference.fa.gz");
    let mut writer = bgzf::io::Writer::new(fs::File::create(&compressed_reference).unwrap());
    writer.write_all(b">chr1\nACGTACGT\n").unwrap();
    writer.try_finish().unwrap();
    fs::write(
        compressed_reference.with_extension("gz.fai"),
        b"chr1\t8\t6\t8\t9\n",
    )
    .unwrap();
    fs::write(
        compressed_reference.with_extension("gz.gzi"),
        0_u64.to_le_bytes(),
    )
    .unwrap();
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        compressed_reference.to_str().unwrap(),
        "--check-ref",
        "fix",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        compressed_reference.to_str().unwrap(),
        "--check-ref",
        "s",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_fix_resolves_missing_and_ambiguous_bases_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGTACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t8\t6\t8\t9\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=8>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\tmissing\t.\tG\t.\tPASS\t.\n\
chr1\t2\tiupac\tY\tR\t.\tPASS\t.\n\
chr1\t3\tunknown\tGN\tNC\t.\tPASS\t.\n\
chr1\t4\tlower\tt\tr\t.\tPASS\t.\n",
    )
    .unwrap();

    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "s",
        input.to_str().unwrap(),
    ])));
    assert_eq!(
        oracle,
        b"chr1\t1\tmissing\tA\tG\t.\tPASS\t.\n\
chr1\t2\tiupac\tC\tA\t.\tPASS\t.\n\
chr1\t3\tunknown\tGT\tNC\t.\tPASS\t.\n\
chr1\t4\tlower\tt\ta\t.\tPASS\t.\n"
    );
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "fix",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_fix_split_atomization_and_origin_trace_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tT\tA,C\t.\tPASS\t.\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "fix",
        "--split-multiallelic",
        "--atomize",
        "--old-rec-tag",
        "ORIG",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "s",
        "-m",
        "-any",
        "--atomize",
        "--old-rec-tag",
        "ORIG",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn reference_fix_removes_alternates_that_become_the_reference_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let reference = directory.path().join("reference.fa");
    let input = directory.path().join("input.vcf");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t4\t6\t4\t5\n").unwrap();
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=4>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t2\t.\tT\tC,C,A\t.\tPASS\t.\tGT\t0/2\t2/1\t1/1\t3/0\n",
    )
    .unwrap();

    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "s",
        input.to_str().unwrap(),
    ])));
    assert_eq!(
        oracle,
        b"chr1\t2\t.\tC\tT,A\t.\tPASS\t.\tGT\t1/0\t0/0\t0/0\t2/1\n"
    );
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--check-ref",
        "fix",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn biallelic_mnv_atomization_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t10\t.\tACGT\tTGCA\t.\tPASS\tDP=5\tGT\t0/1\n\
chr1\t20\t.\tACGT\tAGGA\t.\tPASS\tDP=7\tGT\t1/1\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args(["norm", "--atomize", input.to_str().unwrap()])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn complex_and_multiallelic_atomization_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">\n\
##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">\n\
##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t10\tb1\tAC\tG\t50\tPASS\tIA=7;IR=9,4;IG=1,2,3\tGT:AD:PL\t0/1:9,4:0,10,20\t1:20,7:5,15\n\
chr1\t20\tb2\tAC\tGTG\t50\tPASS\tIA=8;IR=10,5;IG=1,2,3\tGT:AD:PL\t1/1:10,5:0,10,20\t1:20,7:5,15\n\
chr1\t30\tb3\tACG\tAT\t50\tPASS\tIA=9;IR=11,6;IG=1,2,3\tGT:AD:PL\t0/1:11,6:0,10,20\t0:20,7:5,15\n\
chr1\t40\tb4\tCA\tTCG\t50\tPASS\tIA=10;IR=12,7;IG=1,2,3\tGT:AD:PL\t0/1:12,7:0,10,20\t1:20,7:5,15\n\
chr1\t50\tm1\tCC\tC,GG\t50\tPASS\tIA=3,8;IR=10,4,6;IG=0,1,2,3,4,5\tGT:AD:PL\t1/2:10,4,6:0,10,20,30,40,50\t0|2:20,1,7:0,5,10,15,20,25\n",
    )
    .unwrap();

    for overlap in ["*", "."] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--atomize",
            "--atom-overlaps",
            overlap,
            input.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--atomize",
            "--atom-overlaps",
            overlap,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{overlap}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn atom_origin_trace_matches_bcftools_1_24_in_all_encodings() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t50\t.\tA\tC,G\t.\tPASS\t.\n",
    )
    .unwrap();

    let ours_text = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--atomize",
        "--split-multiallelic",
        "--old-rec-tag",
        "ORIG",
        input.to_str().unwrap(),
    ])));
    let oracle_text = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--atomize",
        "-m",
        "-any",
        "--old-rec-tag",
        "ORIG",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours_text, oracle_text);

    for output_type in ["v", "z", "b", "u"] {
        let ours = directory.path().join(format!("ours.{output_type}"));
        let oracle = directory.path().join(format!("oracle.{output_type}"));
        run(
            Command::new(PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))).args([
                "norm",
                "--atomize",
                "--split-multiallelic",
                "--old-rec-tag",
                "ORIG",
                "-O",
                output_type,
                "-o",
                ours.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
        );
        run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--atomize",
            "-m",
            "-any",
            "--old-rec-tag",
            "ORIG",
            "-O",
            output_type,
            "-o",
            oracle.to_str().unwrap(),
            input.to_str().unwrap(),
        ]));
        let ours = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            ours.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            oracle.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{output_type}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn duplicate_policies_match_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("duplicates.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.4\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"SV length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\ts1\tA\tC\t.\tPASS\t.\n\
chr1\t10\ts2\tA\tG\t.\tPASS\t.\n\
chr1\t10\ti1\tA\tAT\t.\tPASS\t.\n\
chr1\t10\ti2\tA\tAG\t.\tPASS\t.\n\
chr1\t10\to1\tA\t<DEL>\t.\tPASS\tSVLEN=-5\n\
chr1\t20\te1\tA\tC\t.\tPASS\t.\n\
chr1\t20\te2\ta\tc\t.\tPASS\t.\n\
chr1\t20\te3\tA\tG\t.\tPASS\t.\n\
chr1\t30\tx1\tN\t<DEL>\t.\tPASS\tSVLEN=-10\n\
chr1\t30\tx2\tN\t<DEL>\t.\tPASS\tSVLEN=-20\n\
chr1\t30\tx3\tN\t<DEL>\t.\tPASS\tSVLEN=-10\n\
chr1\t40\tm1\tA\tC,G\t.\tPASS\t.\n\
chr1\t40\tm2\tA\tG,C\t.\tPASS\t.\n\
chr1\t40\tm3\tA\tC,C\t.\tPASS\t.\n",
    )
    .unwrap();

    for policy in ["snps", "indels", "both", "all", "exact"] {
        let ours = body(run(Command::new(PathBuf::from(env!(
            "CARGO_BIN_EXE_rsomics-vcf"
        )))
        .args([
            "norm",
            "--remove-duplicates",
            policy,
            input.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--rm-dup",
            policy,
            input.to_str().unwrap(),
        ])));
        assert_eq!(ours, oracle, "{policy}");
    }

    for output_type in ["v", "z", "b", "u"] {
        let ours = directory.path().join(format!("ours.{output_type}"));
        run(
            Command::new(PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))).args([
                "norm",
                "--remove-duplicates",
                "exact",
                "-O",
                output_type,
                "-o",
                ours.to_str().unwrap(),
                input.to_str().unwrap(),
            ]),
        );
        let decoded = body(run(Command::new("bcftools").args([
            "view",
            "--no-version",
            "-Ov",
            ours.to_str().unwrap(),
        ])));
        let oracle = body(run(Command::new("bcftools").args([
            "norm",
            "--no-version",
            "--rm-dup",
            "exact",
            input.to_str().unwrap(),
        ])));
        assert_eq!(decoded, oracle, "{output_type}");
    }

    let reference = directory.path().join("reference.fa");
    let realign = directory.path().join("realign.vcf");
    fs::write(&reference, b">chr1\nAAAAAA\n").unwrap();
    fs::write(reference.with_extension("fa.fai"), b"chr1\t6\t6\t6\t7\n").unwrap();
    fs::write(
        &realign,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=6>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\tr1\tA\tAA\t.\tPASS\t.\n\
chr1\t5\tr2\tA\tAA\t.\tPASS\t.\n",
    )
    .unwrap();
    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--remove-duplicates",
        "exact",
        realign.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--fasta-ref",
        reference.to_str().unwrap(),
        "--rm-dup",
        "exact",
        realign.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn exhaustive_short_allele_atomization_matches_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input_path = directory.path().join("input.vcf");
    let sequences = short_sequences();
    let record_count = sequences.len() * (sequences.len() - 1);
    let mut input = format!(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length={}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        record_count * 4 + 4
    );
    let mut position = 1;
    for reference in &sequences {
        for alternate in &sequences {
            if reference == alternate {
                continue;
            }
            writeln!(
                input,
                "chr1\t{position}\t.\t{reference}\t{alternate}\t.\tPASS\t."
            )
            .unwrap();
            position += 4;
        }
    }
    fs::write(&input_path, input).unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args(["norm", "--atomize", input_path.to_str().unwrap()])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "--atomize",
        input_path.to_str().unwrap(),
    ])));
    let ours = String::from_utf8(ours).unwrap();
    let oracle = String::from_utf8(oracle).unwrap();
    let ours_lines: Vec<_> = ours.lines().collect();
    let oracle_lines: Vec<_> = oracle.lines().collect();
    for (index, (ours, oracle)) in ours_lines.iter().zip(&oracle_lines).enumerate() {
        assert_eq!(ours, oracle, "output line {}", index + 1);
    }
    assert_eq!(ours_lines.len(), oracle_lines.len());
}

fn short_sequences() -> Vec<String> {
    let mut sequences = Vec::new();
    for length in 1..=3 {
        let count = 4usize.pow(length);
        for value in 0..count {
            let mut value = value;
            let mut sequence = vec![b'A'; length as usize];
            for base in sequence.iter_mut().rev() {
                *base = b"ACGT"[value % 4];
                value /= 4;
            }
            sequences.push(String::from_utf8(sequence).unwrap());
        }
    }
    sequences
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn split_and_atomize_compose_like_bcftools_1_24() {
    let version = run(Command::new("bcftools").arg("--version"));
    assert!(String::from_utf8_lossy(&version.stdout).starts_with("bcftools 1.24\n"));

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.vcf");
    fs::write(
        &input,
        b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t10\t.\tAC\tGT,AT\t.\tPASS\tAF=0.25,0.5\n",
    )
    .unwrap();

    let ours = body(run(Command::new(PathBuf::from(env!(
        "CARGO_BIN_EXE_rsomics-vcf"
    )))
    .args([
        "norm",
        "--split-multiallelic",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    let oracle = body(run(Command::new("bcftools").args([
        "norm",
        "--no-version",
        "-m",
        "-any",
        "--atomize",
        input.to_str().unwrap(),
    ])));
    assert_eq!(ours, oracle);
}
