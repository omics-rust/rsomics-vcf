use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

fn run(binary: &Path, arguments: &[&str]) {
    let output = Command::new(binary)
        .args(arguments)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn benchmark(c: &mut Criterion) {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    {
        let mut writer = BufWriter::new(input.as_file_mut());
        writer
            .write_all(
                b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            )
            .unwrap();
        for position in 1..=250_000 {
            writeln!(writer, "chr1\t{position}\t.\tA\tC\t50\tPASS\tDP=20").unwrap();
        }
        writer.flush().unwrap();
    }

    let current = std::env::current_exe().unwrap();
    let release = current.parent().unwrap().parent().unwrap();
    let ours = release.join("rsomics-vcf");
    let bcftools = Path::new("bcftools");
    let input_path = input.path().to_str().unwrap();
    let expression = "INFO/DP >= 20 && QUAL >= 30";
    let bytes = input.as_file().metadata().unwrap().len();
    let mut group = c.benchmark_group("filter");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(10));
    group.throughput(Throughput::Bytes(bytes));
    group.bench_function(BenchmarkId::new("rsomics", 250_000), |b| {
        b.iter(|| run(&ours, &["filter", "-i", expression, input_path]));
    });
    group.bench_function(BenchmarkId::new("bcftools_1_24", 250_000), |b| {
        b.iter(|| {
            run(
                bcftools,
                &["filter", "--no-version", "-i", expression, input_path],
            )
        });
    });
    group.bench_function(BenchmarkId::new("rsomics_view", 250_000), |b| {
        b.iter(|| run(&ours, &["view", input_path]));
    });
    group.bench_function(BenchmarkId::new("bcftools_view_1_24", 250_000), |b| {
        b.iter(|| run(bcftools, &["view", "--no-version", input_path]));
    });
    group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
