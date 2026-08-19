use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf"))
}

fn oracle() -> PathBuf {
    env::var_os("RSOMICS_BCFTOOLS")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("bcftools"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/upstream/bcftools-setgt")
        .join(name)
}

fn run(command: &mut Command) -> Output {
    command.output().unwrap()
}

fn success(command: &mut Command, context: &str) -> Output {
    let output = run(command);
    assert!(
        output.status.success(),
        "{context}\nstatus={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn assert_oracle_version() {
    let output = success(Command::new(oracle()).arg("--version"), "bcftools version");
    assert!(
        String::from_utf8_lossy(&output.stdout).starts_with("bcftools 1.24\n"),
        "{}",
        String::from_utf8_lossy(&output.stdout)
    );
}

#[derive(Clone, Copy)]
enum Encoding {
    Vcf,
    VcfGz,
    Bcf,
    RawBcf,
}

impl Encoding {
    const ALL: [Self; 4] = [Self::Vcf, Self::VcfGz, Self::Bcf, Self::RawBcf];

    fn rsomics(self) -> &'static str {
        match self {
            Self::Vcf => "v",
            Self::VcfGz => "z",
            Self::Bcf => "b",
            Self::RawBcf => "u",
        }
    }

    fn bcftools(self) -> &'static str {
        match self {
            Self::Vcf => "-Ov",
            Self::VcfGz => "-Oz",
            Self::Bcf => "-Ob",
            Self::RawBcf => "-Ou",
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Vcf => "vcf",
            Self::VcfGz => "vcf_gz",
            Self::Bcf => "bcf",
            Self::RawBcf => "raw_bcf",
        }
    }
}

#[derive(Debug)]
struct Normalized {
    metadata: Vec<String>,
    columns: String,
    records: Vec<Vec<String>>,
}

fn normalize(path: &Path) -> Normalized {
    let output = success(
        Command::new(oracle())
            .args(["view", "--no-version", "-Ov"])
            .arg(path),
        &format!("decode {}", path.display()),
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let mut metadata = Vec::new();
    let mut columns = None;
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with("##bcftools_") {
            continue;
        }
        if line.starts_with("##") {
            metadata.push(line.to_owned());
        } else if line.starts_with("#CHROM") {
            columns = Some(line.to_owned());
        } else if !line.is_empty() {
            records.push(line.split('\t').map(str::to_owned).collect());
        }
    }
    metadata.sort_unstable();
    Normalized {
        metadata,
        columns: columns.unwrap(),
        records,
    }
}

fn assert_equivalent(label: &str, ours: &Path, oracle_output: &Path) {
    let ours = normalize(ours);
    let oracle_output = normalize(oracle_output);
    assert_eq!(ours.metadata, oracle_output.metadata, "{label}: metadata");
    assert_eq!(ours.columns, oracle_output.columns, "{label}: columns");
    assert_eq!(
        ours.records.len(),
        oracle_output.records.len(),
        "{label}: record count"
    );
    for (record_index, (ours, oracle_record)) in
        ours.records.iter().zip(&oracle_output.records).enumerate()
    {
        assert_eq!(
            ours.len(),
            oracle_record.len(),
            "{label}: record {} field count",
            record_index + 1
        );
        assert_eq!(
            &ours[..9],
            &oracle_record[..9],
            "{label}: record {} site fields",
            record_index + 1
        );
        for (sample_index, (ours, oracle_sample)) in
            ours[9..].iter().zip(&oracle_record[9..]).enumerate()
        {
            assert_eq!(
                ours,
                oracle_sample,
                "{label}: record {} sample {}",
                record_index + 1,
                sample_index + 1
            );
        }
    }
}

fn safe_name(value: &str) -> String {
    value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() {
                character
            } else {
                '_'
            }
        })
        .collect()
}

fn output_paths(directory: &Path, name: &str, encoding: Encoding) -> (PathBuf, PathBuf) {
    let stem = format!("{}-{}", safe_name(name), encoding.name());
    (
        directory.join(format!("{stem}-ours.out")),
        directory.join(format!("{stem}-oracle.out")),
    )
}

fn ours_command(input: &Path, output: &Path, encoding: Encoding, arguments: &[&str]) -> Command {
    let mut command = Command::new(binary());
    command
        .arg("setgt")
        .args(arguments)
        .args(["-O", encoding.rsomics(), "-o"])
        .arg(output)
        .arg(input);
    command
}

fn oracle_command(input: &Path, output: &Path, encoding: Encoding, arguments: &[&str]) -> Command {
    let mut command = Command::new(oracle());
    command
        .arg("+setGT")
        .arg(input)
        .args(["--no-version", encoding.bcftools(), "-o"])
        .arg(output)
        .arg("--")
        .args(arguments);
    command
}

fn compare_case(
    directory: &Path,
    name: &str,
    input: &Path,
    encoding: Encoding,
    arguments: &[&str],
) {
    let (ours, expected) = output_paths(directory, name, encoding);
    success(
        &mut ours_command(input, &ours, encoding, arguments),
        &format!("{name}: rsomics"),
    );
    success(
        &mut oracle_command(input, &expected, encoding, arguments),
        &format!("{name}: bcftools"),
    );
    assert_equivalent(
        &format!("{name}: output={}", encoding.name()),
        &ours,
        &expected,
    );
}

fn encode(input: &Path, output: &Path, encoding: Encoding) {
    success(
        Command::new(oracle())
            .args(["view", "--no-version", encoding.bcftools(), "-o"])
            .arg(output)
            .arg(input),
        &format!("encode {}", encoding.name()),
    );
}

fn diploid_core(directory: &Path) -> PathBuf {
    let destination = directory.join("diploid.vcf");
    let source = fs::read_to_string(fixture("core.vcf")).unwrap();
    let filtered = source
        .lines()
        .filter(|line| !line.starts_with("1\t200\t"))
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(&destination, format!("{filtered}\n")).unwrap();
    destination
}

fn record(path: &Path, index: usize) -> Vec<String> {
    normalize(path).records[index].clone()
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn normal_targets_match_bcftools_1_24() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let core = fixture("core.vcf");
    let diploid = diploid_core(directory.path());
    for (name, input, arguments) in [
        ("any_missing", core.as_path(), &["-t", ".", "-n", "0"][..]),
        ("partial_missing", core.as_path(), &["-t", "./x", "-n", "."]),
        (
            "complete_missing",
            core.as_path(),
            &["-t", "./.", "-n", "0p"],
        ),
        ("all", core.as_path(), &["-t", "a", "-n", "0"]),
        (
            "query_include",
            core.as_path(),
            &["-t", "q", "-i", "FMT/DP>=10", "-n", "."],
        ),
        (
            "query_exclude",
            core.as_path(),
            &["-t", "q", "-e", "FMT/DP>=10", "-n", "."],
        ),
        (
            "binomial_lt",
            diploid.as_path(),
            &["-t", "b:AD<0.5", "-n", "0"],
        ),
        (
            "binomial_le",
            diploid.as_path(),
            &["-t", "b:AD<=0.5", "-n", "0"],
        ),
        (
            "binomial_eq",
            diploid.as_path(),
            &["-t", "b:AD=0.5", "-n", "0"],
        ),
        (
            "binomial_eq2",
            diploid.as_path(),
            &["-t", "b:AD==0.5", "-n", "0"],
        ),
        (
            "binomial_ge",
            diploid.as_path(),
            &["-t", "b:AD>=0.5", "-n", "0"],
        ),
        (
            "binomial_gt",
            diploid.as_path(),
            &["-t", "b:AD>0.5", "-n", "0"],
        ),
        (
            "random",
            core.as_path(),
            &["-t", "r:0.4", "-s", "7", "-n", "."],
        ),
        (
            "random_composed",
            core.as_path(),
            &["-t", ".", "-t", "r:0.6", "-s", "7", "-n", "0"],
        ),
    ] {
        compare_case(directory.path(), name, input, Encoding::Vcf, arguments);
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn replacement_families_match_in_every_output_encoding() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let core = fixture("core.vcf");
    let diploid = diploid_core(directory.path());
    for replacement in [
        ".", "0", "0p", "m", "mp", "M", "Mp", "X", "p", "u", "i", "c:0/1", "c:2|0", "c:m/M",
        "c:0/X", "c:./.",
    ] {
        let input = if replacement == "i" {
            diploid.as_path()
        } else {
            core.as_path()
        };
        for encoding in Encoding::ALL {
            compare_case(
                directory.path(),
                &format!("replacement_{replacement}"),
                input,
                encoding,
                &["-t", "a", "-n", replacement],
            );
        }
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn every_input_encoding_matches_a_pairwise_output() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let source = fixture("core.vcf");
    for (input_encoding, output_encoding) in [
        (Encoding::Vcf, Encoding::Bcf),
        (Encoding::VcfGz, Encoding::RawBcf),
        (Encoding::Bcf, Encoding::VcfGz),
        (Encoding::RawBcf, Encoding::Vcf),
    ] {
        let input = directory
            .path()
            .join(format!("input-{}.data", input_encoding.name()));
        encode(&source, &input, input_encoding);
        compare_case(
            directory.path(),
            &format!("input_{}", input_encoding.name()),
            &input,
            output_encoding,
            &["-t", ".", "-n", "0"],
        );
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn query_inversion_edits_the_selected_sample_instead_of_sample_zero() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("core.vcf");
    let (ours, upstream) = output_paths(directory.path(), "query_inversion", Encoding::Vcf);
    let arguments = ["-t", "q", "-n", "i", "-i", "FMT/DP=7"];
    success(
        &mut ours_command(&input, &ours, Encoding::Vcf, &arguments),
        "query inversion: rsomics",
    );
    success(
        &mut oracle_command(&input, &upstream, Encoding::Vcf, &arguments),
        "query inversion: bcftools",
    );
    let ours = record(&ours, 0);
    let upstream = record(&upstream, 0);
    assert!(upstream[9].starts_with("1/0:"));
    assert!(upstream[11].starts_with("./1:"));
    assert!(ours[9].starts_with("0/1:"));
    assert!(ours[11].starts_with("1/.:"));
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn missing_required_ad_fails_instead_of_becoming_an_upstream_no_op() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("malformed-ad.vcf");
    let (ours, upstream) = output_paths(directory.path(), "missing_ad", Encoding::Vcf);
    let arguments = ["-t", "a", "-n", "X"];
    let upstream_result = run(&mut oracle_command(
        &input,
        &upstream,
        Encoding::Vcf,
        &arguments,
    ));
    assert!(upstream_result.status.success());
    let upstream_record = record(&upstream, 0);
    assert_eq!(upstream_record[8], "GT");
    assert_eq!(upstream_record[9], "0/1");
    assert_eq!(upstream_record[10], "0/0");

    let ours_result = run(&mut ours_command(&input, &ours, Encoding::Vcf, &arguments));
    assert!(!ours_result.status.success());
    assert!(!ours.exists());
    assert!(String::from_utf8_lossy(&ours_result.stderr).contains("FORMAT/AD"));
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn existing_ac_an_are_reconciled_instead_of_left_stale() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/setgt.vcf");
    let (ours, upstream) = output_paths(directory.path(), "ac_an", Encoding::Vcf);
    let arguments = ["-t", ".", "-n", "0"];
    success(
        &mut ours_command(&input, &ours, Encoding::Vcf, &arguments),
        "AC/AN: rsomics",
    );
    success(
        &mut oracle_command(&input, &upstream, Encoding::Vcf, &arguments),
        "AC/AN: bcftools",
    );
    assert_eq!(record(&upstream, 0)[7], "AC=1;AN=3");
    assert_eq!(record(&ours, 0)[7], "AC=1;AN=6");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn ambiguous_replacements_fail_instead_of_using_upstream_mask_precedence() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("core.vcf");
    for (ambiguous, canonical) in [("0u", "u"), ("Xp", "X")] {
        let (upstream, expected) = output_paths(
            directory.path(),
            &format!("ambiguous_{ambiguous}"),
            Encoding::Vcf,
        );
        success(
            &mut oracle_command(
                &input,
                &upstream,
                Encoding::Vcf,
                &["-t", "a", "-n", ambiguous],
            ),
            &format!("ambiguous {ambiguous}: bcftools"),
        );
        success(
            &mut oracle_command(
                &input,
                &expected,
                Encoding::Vcf,
                &["-t", "a", "-n", canonical],
            ),
            &format!("canonical {canonical}: bcftools"),
        );
        assert_equivalent(
            &format!("bcftools precedence {ambiguous} -> {canonical}"),
            &upstream,
            &expected,
        );

        let ours = directory.path().join(format!("ours-{ambiguous}.vcf"));
        let result = run(&mut ours_command(
            &input,
            &ours,
            Encoding::Vcf,
            &["-t", "a", "-n", ambiguous],
        ));
        assert!(!result.status.success(), "{ambiguous}");
        assert!(!ours.exists(), "{ambiguous}");
    }
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn named_destination_survives_a_late_error_instead_of_being_truncated() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("late.vcf");
    let mut source = fs::read_to_string(fixture("core.vcf")).unwrap();
    source.push_str("1\tBAD\n");
    fs::write(&input, source).unwrap();

    let upstream = directory.path().join("upstream.vcf");
    fs::write(&upstream, b"existing").unwrap();
    let upstream_result = run(&mut oracle_command(
        &input,
        &upstream,
        Encoding::Vcf,
        &["-t", "a", "-n", "0"],
    ));
    assert!(!upstream_result.status.success());
    assert_eq!(fs::read(&upstream).unwrap(), b"");

    let ours = directory.path().join("ours.vcf");
    fs::write(&ours, b"existing").unwrap();
    let ours_result = run(&mut ours_command(
        &input,
        &ours,
        Encoding::Vcf,
        &["-t", "a", "-n", "0"],
    ));
    assert!(!ours_result.status.success());
    assert_eq!(fs::read(&ours).unwrap(), b"existing");
}

#[test]
#[ignore = "release oracle: requires bcftools 1.24"]
fn typed_ploidy_rules_do_not_follow_record_width_accidents() {
    assert_oracle_version();
    let directory = tempfile::tempdir().unwrap();
    let input = fixture("core.vcf");

    let (ours, upstream) = output_paths(directory.path(), "partial_width", Encoding::Vcf);
    let arguments = ["-t", "./x", "-n", "0"];
    success(
        &mut ours_command(&input, &ours, Encoding::Vcf, &arguments),
        "partial width: rsomics",
    );
    success(
        &mut oracle_command(&input, &upstream, Encoding::Vcf, &arguments),
        "partial width: bcftools",
    );
    assert!(record(&upstream, 0)[12].starts_with("0/0:"));
    assert!(record(&ours, 0)[12].starts_with("./.:"));

    let (ours, upstream) = output_paths(directory.path(), "binomial_width", Encoding::Vcf);
    let arguments = ["-t", "b:AD>=0.5", "-n", "0"];
    success(
        &mut ours_command(&input, &ours, Encoding::Vcf, &arguments),
        "binomial width: rsomics",
    );
    success(
        &mut oracle_command(&input, &upstream, Encoding::Vcf, &arguments),
        "binomial width: bcftools",
    );
    assert!(record(&upstream, 1)[13].starts_with("0/0/0:"));
    assert!(record(&ours, 1)[13].starts_with("0/1/2:"));

    let (ours, upstream) = output_paths(directory.path(), "invert_width", Encoding::Vcf);
    let arguments = ["-t", "a", "-n", "i"];
    success(
        &mut ours_command(&input, &ours, Encoding::Vcf, &arguments),
        "invert width: rsomics",
    );
    success(
        &mut oracle_command(&input, &upstream, Encoding::Vcf, &arguments),
        "invert width: bcftools",
    );
    assert!(record(&upstream, 1)[9].starts_with("1/2:"));
    assert!(record(&ours, 1)[9].starts_with("2/1:"));
}
