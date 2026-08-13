# Performance

## `annotate` 0.4 release gate

The release benchmark compares complete normalized VCF output against
bcftools 1.24 before any timed run. It uses three warmups followed by ten
measured pairs in alternating order. Wall time, user and system CPU time, and
peak resident memory come from macOS `/usr/bin/time -lp`.

The measured implementation revision is
`f5a76ed58dcc20496eec3b160d245af358c9c80b`. The worktree was clean. The
result ledger is retained at
`/Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111`.

### Decision

The release gate passes on strict resource-use advantage. The interval join
uses 26.02 times less peak memory than bcftools and the typed transfer uses
1.64 times less. Throughput is not a win: the interval join is 5.26% slower
and the typed transfer is 5.66 times slower. The 0.4 release can therefore be
described as bounded-memory annotation with verified compatibility, but not as
a generally faster replacement for bcftools. Typed transfer throughput remains
an explicit optimization target after release.

| Workload | Tool | Median wall | p99 wall | Median peak RSS | Result |
|---|---:|---:|---:|---:|---|
| 2,000,000 interval records | rsomics | 3.900 s | 4.000 s | 4,767,744 B | 5.26% slower, 96.16% less RSS |
| 2,000,000 interval records | bcftools | 3.705 s | 3.820 s | 124,059,648 B | reference |
| 300,000 typed records, 8 samples | rsomics | 8.825 s | 8.920 s | 5,046,272 B | 5.66 times slower, 38.89% less RSS |
| 300,000 typed records, 8 samples | bcftools | 1.560 s | 1.580 s | 8,257,536 B | reference |

### Workloads and equality

The interval workload contains 2,000,000 biallelic target SNVs and one
coordinate-sorted bounded interval per target. It transfers one integer
`INFO/SCORE` value from a BGZF and tabix-indexed tabular source. The target is
57,444,666 bytes, the source is 10,660,527 bytes, and each plain VCF output is
73,224,666 bytes.

The typed workload contains 300,000 triallelic records and eight samples. It
transfers `Number=A/R/G` integer INFO values plus `GT`, `AD`, and `PL`; the
source sample header is the reverse of the target header. The BGZF target is
1,176,699 bytes, the BGZF source is 1,881,545 bytes, and each plain VCF output
is 95,289,421 bytes.

The typed performance fixture deliberately keeps the ALT order identical.
bcftools 1.24 copies allele-indexed values and genotypes without remapping when
the source ALT order differs, while rsomics remaps them. That case cannot
produce semantically equal output and is covered by the compatibility tests
instead of being used as a favorable performance comparison.

Canonicalization removes only `##bcftools_` provenance lines and sorts metadata
lines. The complete remaining header and record body must compare byte for
byte. The canonical SHA-256 values match for both workloads:

| Workload | rsomics canonical output | bcftools canonical output |
|---|---|---|
| interval | `f9354d0fe551aaca2f9c3ac14d176638840ce453d43687834ee1dd7172e356a9` | `f9354d0fe551aaca2f9c3ac14d176638840ce453d43687834ee1dd7172e356a9` |
| typed | `0e2eba0b57db893a2fd564b31eee5391b1181b69bf3ab53fa2f53954061c3ee6` | `0e2eba0b57db893a2fd564b31eee5391b1181b69bf3ab53fa2f53954061c3ee6` |

### Environment and command

- Mac14,3 with Apple M2 and 8,589,934,592 bytes physical memory
- macOS 26.6.1 build 25G76, Darwin 25.6.0 arm64
- rustc 1.97.1, commit `8bab26f4f68e0e26f0bb7960be334d5b520ea452`
- rsomics-vcf 0.3.0 release binary built from the measured revision
- bcftools 1.24 with HTSlib 1.24

```console
env RSOMICS_RUSTC=/opt/homebrew/Cellar/rust/1.97.1/bin/rustc \
  benchmarks/annotate-vs-bcftools.sh generate \
  /Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf \
  /opt/homebrew/bin/bcftools \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111
```

The measured commands were:

```console
/Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf annotate \
  -a /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-source.tsv.gz \
  -c CHROM,FROM,TO,INFO/SCORE \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-target.vcf

/opt/homebrew/Cellar/bcftools/1.24/bin/bcftools annotate --no-version \
  -a /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-source.tsv.gz \
  -c CHROM,FROM,TO,INFO/SCORE \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-target.vcf

/Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf annotate \
  -a /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-source.vcf.gz \
  -c INFO/IA,INFO/IR,INFO/IG,FORMAT/GT,FORMAT/AD,FORMAT/PL \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-target.vcf.gz

/opt/homebrew/Cellar/bcftools/1.24/bin/bcftools annotate --no-version \
  -a /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-source.vcf.gz \
  -c INFO/IA,INFO/IR,INFO/IG,FORMAT/GT,FORMAT/AD,FORMAT/PL \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-target.vcf.gz
```

### Paired summaries

| Workload | Tool | Metric | n | Median | p99 | Mean | Sample standard deviation |
|---|---|---|---:|---:|---:|---:|---:|
| interval | rsomics | wall seconds | 10 | 3.900 | 4.000 | 3.913 | 0.039455 |
| interval | rsomics | user seconds | 10 | 3.790 | 3.840 | 3.793 | 0.030930 |
| interval | rsomics | system seconds | 10 | 0.065 | 0.070 | 0.060 | 0.012472 |
| interval | rsomics | peak RSS bytes | 10 | 4,767,744 | 4,800,512 | 4,769,382 | 23,743 |
| interval | bcftools | wall seconds | 10 | 3.705 | 3.820 | 3.720 | 0.045216 |
| interval | bcftools | user seconds | 10 | 3.575 | 3.650 | 3.586 | 0.030984 |
| interval | bcftools | system seconds | 10 | 0.080 | 0.080 | 0.073 | 0.011595 |
| interval | bcftools | peak RSS bytes | 10 | 124,059,648 | 126,222,336 | 123,825,357 | 1,807,331 |
| typed | rsomics | wall seconds | 10 | 8.825 | 8.920 | 8.834 | 0.049261 |
| typed | rsomics | user seconds | 10 | 8.595 | 8.650 | 8.599 | 0.038137 |
| typed | rsomics | system seconds | 10 | 0.130 | 0.130 | 0.128 | 0.004216 |
| typed | rsomics | peak RSS bytes | 10 | 5,046,272 | 5,062,656 | 5,047,910 | 16,293 |
| typed | bcftools | wall seconds | 10 | 1.560 | 1.580 | 1.560 | 0.011547 |
| typed | bcftools | user seconds | 10 | 1.510 | 1.530 | 1.512 | 0.009189 |
| typed | bcftools | system seconds | 10 | 0.020 | 0.020 | 0.020 | 0.000000 |
| typed | bcftools | peak RSS bytes | 10 | 8,257,536 | 8,290,304 | 8,254,259 | 28,691 |

### Raw paired distribution

Order is the execution order within each pair. RSS values are bytes.

```text
workload pair order rs_wall rs_user rs_sys rs_rss bcf_wall bcf_user bcf_sys bcf_rss
interval 1 rs>bcf 3.90 3.79 0.04 4767744 3.67 3.57 0.05 124174336
interval 2 bcf>rs 4.00 3.84 0.06 4784128 3.75 3.62 0.06 123944960
interval 3 rs>bcf 3.89 3.79 0.04 4767744 3.71 3.57 0.06 122322944
interval 4 bcf>rs 3.90 3.76 0.05 4800512 3.82 3.65 0.08 120520704
interval 5 rs>bcf 3.94 3.83 0.06 4718592 3.68 3.55 0.08 124796928
interval 6 bcf>rs 3.90 3.79 0.07 4767744 3.76 3.61 0.08 126205952
interval 7 rs>bcf 3.88 3.77 0.07 4751360 3.72 3.58 0.08 123305984
interval 8 bcf>rs 3.92 3.81 0.07 4800512 3.69 3.56 0.08 124616704
interval 9 rs>bcf 3.86 3.74 0.07 4767744 3.70 3.57 0.08 122142720
interval 10 bcf>rs 3.94 3.81 0.07 4767744 3.70 3.58 0.08 126222336
typed 1 rs>bcf 8.85 8.62 0.12 5062656 1.55 1.50 0.02 8257536
typed 2 bcf>rs 8.83 8.60 0.13 5062656 1.55 1.51 0.02 8273920
typed 3 rs>bcf 8.81 8.58 0.13 5046272 1.56 1.51 0.02 8290304
typed 4 bcf>rs 8.81 8.58 0.13 5013504 1.57 1.52 0.02 8290304
typed 5 rs>bcf 8.92 8.65 0.13 5062656 1.56 1.51 0.02 8257536
typed 6 bcf>rs 8.89 8.64 0.13 5062656 1.56 1.51 0.02 8208384
typed 7 rs>bcf 8.82 8.59 0.13 5046272 1.57 1.52 0.02 8241152
typed 8 bcf>rs 8.78 8.55 0.13 5046272 1.54 1.50 0.02 8208384
typed 9 rs>bcf 8.76 8.54 0.12 5046272 1.56 1.51 0.02 8257536
typed 10 bcf>rs 8.87 8.64 0.13 5029888 1.58 1.53 0.02 8257536
```

### Fingerprints

| Artifact | SHA-256 |
|---|---|
| benchmark harness | `d7c75e785f142764729bbf218d479e8f04ad1d1636a0c15d46e0aef662d06ddc` |
| rsomics-vcf binary | `2f1f17b751a5749f464d099faea157b66e6383ce3719628f879099d8557c0d6c` |
| bcftools binary | `33100a6b961c529e915394d53b4737a0f8dd7a164eac352afe4e74e1ced51f60` |
| rustc binary | `d69d40bfd2e11825feb3538512b6ffcd63de91c35ec36bb876849f0f9f8fe6bd` |
| interval target | `46e9f02024f8cc62ab2cf09a6025b4551e675182713fdd3d902b632bf543e15c` |
| interval source | `b3d2fa693fd8e55ed61ce45e89af79662c1bdd37fffb951901255dc13fa3dc10` |
| interval raw distribution | `2b8e1f544b17e40c09d1cde827a30d4d5690f718711d096f09c11ab511fb225f` |
| interval summary | `9ca4a2c193cf9880755ba24b53b65231c4fb8245f80678be424a70215d3c8d03` |
| typed target | `185e2082e08b815ada2ec06cb63bf9b087d837f6d5396880bed2ab24f10e7c53` |
| typed source | `c33c8b0741ca5138a3237181f152daab3109fdf352265b2f6112b82c12429448` |
| typed raw distribution | `84b7a47da9345863d0743a1fc86c82245f550b3f058bd6f19ecc867a34fce68d` |
| typed summary | `c96c68fce3f7b7c4192322568b4158d3717f21ad4da578d18e42963efc3800af` |
