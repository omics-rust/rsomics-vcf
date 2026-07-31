use std::io::{BufWriter, Write};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rsomics_vcf::head::{self, Options};

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
        for position in 1..=100_000 {
            writeln!(writer, "chr1\t{position}\t.\tA\tC\t50.0\tPASS\tDP=20").unwrap();
        }
        writer.flush().unwrap();
    }
    let bytes = input.as_file().metadata().unwrap().len();
    let mut group = c.benchmark_group("head");
    group.throughput(Throughput::Bytes(bytes));
    group.bench_with_input(
        BenchmarkId::new("vcf_records", 100_000),
        &input,
        |b, input| {
            b.iter(|| {
                head::write(
                    input.path(),
                    Options {
                        records: 100_000,
                        ..Options::default()
                    },
                    std::io::sink(),
                )
                .unwrap()
            });
        },
    );
    group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
