# Performance

## `reheader` 0.5 release gate

The release benchmark compares complete plain VCF output and decompressed BGZF
VCF output against bcftools 1.24 before timing. It also requires the BGZF
record frames and canonical EOF written by rsomics to match the input byte for
byte. Three warmups precede ten measured pairs in alternating order. Wall time,
CPU time, and peak resident memory come from macOS `/usr/bin/time -lp`.

The measured implementation revision is
`b1e50bd0cc998b4c6da7768f61a9244f38a9a81e`. The worktree was clean. The
generated inputs, raw measurements, summaries, decisions, and artifact hashes
are retained under
`/Volumes/KIOXIA/Developments/tmp/rsomics-vcf-reheader-gate-20260818-b1e50bd`.

### Decision

Both paths pass the release gate. Plain VCF reheadering is 5.97 times faster
and uses 38.72% less peak memory than bcftools. BGZF reheadering is twice as
slow at this input size, but uses 33.88% less peak memory. The BGZF result is a
strict memory advantage, not a throughput claim.

| Input path | Tool | Median wall | p99 wall | Median peak RSS | Result |
|---|---|---:|---:|---:|---|
| 2,000,000-record plain VCF | rsomics | 0.515 s | 0.560 s | 3,915,776 B | 5.97 times faster, 38.72% less RSS |
| 2,000,000-record plain VCF | bcftools | 3.075 s | 3.370 s | 6,389,760 B | reference |
| 2,000,000-record BGZF VCF | rsomics | 0.040 s | 0.110 s | 4,636,672 B | 2.00 times slower, 33.88% less RSS |
| 2,000,000-record BGZF VCF | bcftools | 0.020 s | 0.030 s | 7,012,352 B | reference |

### Workload and equality

The fixture contains 2,000,000 biallelic SNVs and eight diploid samples. The
plain input is 188,889,310 bytes and the BGZF input is 4,458,189 bytes. Each
run replaces the complete header, synchronizes contigs from a FAI, and renames
all samples from a two-column map.

The record-body SHA-256 is
`269ad78c213a2c6ba641b4356ecdb8a8692656b48c91309484229bdc949db662`
for both input encodings and both tools. The input and rsomics BGZF raw-tail
SHA-256 is
`24f39c0939ffce222637b02569fb2f75d7f0619f80104572bec84af76bdadd8f`.
The complete plain outputs compare byte for byte; BGZF outputs compare after
decompression because block boundaries are not part of the VCF contract.

### Environment and command

- Mac14,3 with Apple M2 and 8,589,934,592 bytes physical memory
- macOS 26.6.1 build 25G76, Darwin 25.6.0 arm64
- rustc 1.97.1, commit `8bab26f4f68e0e26f0bb7960be334d5b520ea452`
- rsomics-vcf 0.4.0 release binary built from the measured revision
- bcftools 1.24 with HTSlib 1.24

```console
env RSOMICS_RUSTC=/opt/homebrew/Cellar/rust/1.97.1/bin/rustc \
  benchmarks/reheader-vs-bcftools.sh run \
  /Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf \
  /opt/homebrew/bin/bcftools \
  /Volumes/KIOXIA/Developments/tmp/rsomics-vcf-reheader-gate-20260818-b1e50bd \
  /Volumes/KIOXIA/Developments/tmp/rsomics-vcf-reheader-gate-20260818-b1e50bd/results
```

### Measured distributions

| Path | Tool | Metric | n | Median | p99 | Mean | Sample standard deviation |
|---|---|---|---:|---:|---:|---:|---:|
| plain | rsomics | wall seconds | 10 | 0.515 | 0.560 | 0.523 | 0.018886 |
| plain | rsomics | user seconds | 10 | 0.000 | 0.000 | 0.000 | 0.000000 |
| plain | rsomics | system seconds | 10 | 0.170 | 0.180 | 0.170 | 0.004714 |
| plain | rsomics | peak RSS bytes | 10 | 3,915,776 | 3,915,776 | 3,909,222 | 8,461 |
| plain | bcftools | wall seconds | 10 | 3.075 | 3.370 | 3.115 | 0.131593 |
| plain | bcftools | user seconds | 10 | 0.190 | 0.210 | 0.194 | 0.010750 |
| plain | bcftools | system seconds | 10 | 2.770 | 2.960 | 2.771 | 0.083593 |
| plain | bcftools | peak RSS bytes | 10 | 6,389,760 | 6,406,144 | 6,393,037 | 10,362 |
| BGZF | rsomics | wall seconds | 10 | 0.040 | 0.110 | 0.048 | 0.022010 |
| BGZF | rsomics | user seconds | 10 | 0.000 | 0.000 | 0.000 | 0.000000 |
| BGZF | rsomics | system seconds | 10 | 0.000 | 0.010 | 0.003 | 0.004830 |
| BGZF | rsomics | peak RSS bytes | 10 | 4,636,672 | 4,653,056 | 4,638,310 | 14,346 |
| BGZF | bcftools | wall seconds | 10 | 0.020 | 0.030 | 0.021 | 0.003162 |
| BGZF | bcftools | user seconds | 10 | 0.000 | 0.000 | 0.000 | 0.000000 |
| BGZF | bcftools | system seconds | 10 | 0.000 | 0.000 | 0.000 | 0.000000 |
| BGZF | bcftools | peak RSS bytes | 10 | 7,012,352 | 7,028,736 | 7,004,160 | 15,922 |

### Fingerprints

| Artifact | SHA-256 |
|---|---|
| benchmark harness | `ae317ecd5548920e83086227f58ff3f77cec5eea60a21687fd9f6b4d6016db06` |
| rsomics-vcf binary | `5d07c74749efe4120ac6be49bec5e74d0127309938996bb3624893fbe03c71a4` |
| bcftools binary | `33100a6b961c529e915394d53b4737a0f8dd7a164eac352afe4e74e1ced51f60` |
| rustc binary | `d69d40bfd2e11825feb3538512b6ffcd63de91c35ec36bb876849f0f9f8fe6bd` |
| plain input | `d18a222e7b72242955d1370c0b55e366eac4b732bc4ec1112cedaf5a0fed4e84` |
| BGZF input | `799007fe927c2eb5f8746627f6ce85034f1c794525377314e3bcd9613990bd4b` |
| raw distribution | `127648bca21912fd143d0ea44dd9dfb46c1b454f473b67be509f706ff8eb80a0` |
| summary | `76b0995ddc7f8af9975aed4e359e320fd2dee27e09b8ca56257d189d810f3a13` |
| decision | `a67a43ef8cf5a3b4ff628127cc6189019a89de84566df713576efad1cfaaed83` |
| equality ledger | `456a653d5634f1fe8456469fc0a528248dc1d9160408c34bddb7121bedbed139` |

## `annotate` 0.4 release gate

The release benchmark compares complete normalized VCF output against
bcftools 1.24 before any timed run. It uses three warmups followed by ten
measured pairs in alternating order. Wall time, user and system CPU time, and
peak resident memory come from macOS `/usr/bin/time -lp`.

The measured implementation revision is
`799a7681e64f2f869d00766ca21ae4d62085d129`. The worktree was clean. Inputs
and their generation manifest are retained under
`/Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111`;
the final result ledger is retained under
`/Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-799a768`.

### Decision

The release gate passes on strict resource-use advantage. The interval join
uses 26.00 times less peak memory than bcftools and the typed transfer uses
1.63 times less. Throughput is not a win: the interval join is 4.28% slower
and the typed transfer is 5.45 times slower. The 0.4 release can therefore be
described as bounded-memory annotation with verified compatibility, but not as
a generally faster replacement for bcftools. Typed transfer throughput remains
an explicit optimization target after release.

| Workload | Tool | Median wall | p99 wall | Median peak RSS | Result |
|---|---:|---:|---:|---:|---|
| 2,000,000 interval records | rsomics | 3.895 s | 3.980 s | 4,784,128 B | 4.28% slower, 96.15% less RSS |
| 2,000,000 interval records | bcftools | 3.735 s | 3.810 s | 124,395,520 B | reference |
| 300,000 typed records, 8 samples | rsomics | 8.495 s | 8.570 s | 5,070,848 B | 5.45 times slower, 38.59% less RSS |
| 300,000 typed records, 8 samples | bcftools | 1.560 s | 1.590 s | 8,257,536 B | reference |

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
  benchmarks/annotate-vs-bcftools.sh run \
  /Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf \
  /opt/homebrew/bin/bcftools \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-target.vcf \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/interval-source.tsv.gz \
  CHROM,FROM,TO,INFO/SCORE \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-799a768/results/interval

env RSOMICS_RUSTC=/opt/homebrew/Cellar/rust/1.97.1/bin/rustc \
  benchmarks/annotate-vs-bcftools.sh run \
  /Volumes/KIOXIA/Developments/cargo-target/rsomics-vcf/release/rsomics-vcf \
  /opt/homebrew/bin/bcftools \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-target.vcf.gz \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-0111/inputs/typed-source.vcf.gz \
  INFO/IA,INFO/IR,INFO/IG,FORMAT/GT,FORMAT/AD,FORMAT/PL \
  /Volumes/KIOXIA/Developments/tmp/rsomics-annotate-benchmark-20260813-799a768/results/typed
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
| interval | rsomics | wall seconds | 10 | 3.895 | 3.980 | 3.895 | 0.052757 |
| interval | rsomics | user seconds | 10 | 3.785 | 3.830 | 3.778 | 0.046380 |
| interval | rsomics | system seconds | 10 | 0.070 | 0.100 | 0.064 | 0.017127 |
| interval | rsomics | peak RSS bytes | 10 | 4,784,128 | 4,800,512 | 4,789,043 | 11,058 |
| interval | bcftools | wall seconds | 10 | 3.735 | 3.810 | 3.741 | 0.040947 |
| interval | bcftools | user seconds | 10 | 3.585 | 3.640 | 3.596 | 0.025906 |
| interval | bcftools | system seconds | 10 | 0.080 | 0.100 | 0.075 | 0.013540 |
| interval | bcftools | peak RSS bytes | 10 | 124,395,520 | 126,222,336 | 124,284,109 | 1,478,702 |
| typed | rsomics | wall seconds | 10 | 8.495 | 8.570 | 8.500 | 0.035901 |
| typed | rsomics | user seconds | 10 | 8.280 | 8.340 | 8.278 | 0.030840 |
| typed | rsomics | system seconds | 10 | 0.120 | 0.130 | 0.121 | 0.003162 |
| typed | rsomics | peak RSS bytes | 10 | 5,070,848 | 5,079,040 | 5,061,018 | 27,252 |
| typed | bcftools | wall seconds | 10 | 1.560 | 1.590 | 1.565 | 0.010801 |
| typed | bcftools | user seconds | 10 | 1.520 | 1.530 | 1.517 | 0.006749 |
| typed | bcftools | system seconds | 10 | 0.020 | 0.020 | 0.020 | 0.000000 |
| typed | bcftools | peak RSS bytes | 10 | 8,257,536 | 8,306,688 | 8,260,813 | 28,691 |

### Raw paired distribution

Order is the execution order within each pair. RSS values are bytes.

```text
workload pair order rs_wall rs_user rs_sys rs_rss bcf_wall bcf_user bcf_sys bcf_rss
interval 1 rs>bcf 3.98 3.83 0.10 4800512 3.77 3.58 0.10 125304832
interval 2 bcf>rs 3.93 3.81 0.05 4800512 3.81 3.64 0.07 121798656
interval 3 rs>bcf 3.92 3.81 0.05 4784128 3.79 3.62 0.08 122781696
interval 4 bcf>rs 3.90 3.79 0.05 4800512 3.71 3.59 0.05 123813888
interval 5 rs>bcf 3.80 3.71 0.04 4784128 3.71 3.58 0.06 126222336
interval 6 bcf>rs 3.87 3.76 0.07 4767744 3.76 3.63 0.07 122978304
interval 7 rs>bcf 3.89 3.78 0.07 4800512 3.70 3.58 0.08 124846080
interval 8 bcf>rs 3.83 3.69 0.07 4784128 3.69 3.56 0.08 126189568
interval 9 rs>bcf 3.94 3.82 0.07 4784128 3.72 3.58 0.08 123944960
interval 10 bcf>rs 3.89 3.78 0.07 4784128 3.75 3.60 0.08 124960768
typed 1 rs>bcf 8.49 8.28 0.12 5079040 1.56 1.52 0.02 8290304
typed 2 bcf>rs 8.55 8.31 0.12 5062656 1.56 1.52 0.02 8273920
typed 3 rs>bcf 8.48 8.26 0.12 5079040 1.59 1.53 0.02 8257536
typed 4 bcf>rs 8.57 8.34 0.13 5029888 1.56 1.51 0.02 8224768
typed 5 rs>bcf 8.46 8.24 0.12 4997120 1.57 1.51 0.02 8224768
typed 6 bcf>rs 8.51 8.29 0.12 5062656 1.56 1.51 0.02 8290304
typed 7 rs>bcf 8.50 8.28 0.12 5079040 1.57 1.52 0.02 8257536
typed 8 bcf>rs 8.48 8.26 0.12 5079040 1.55 1.51 0.02 8241152
typed 9 rs>bcf 8.50 8.28 0.12 5079040 1.57 1.52 0.02 8241152
typed 10 bcf>rs 8.46 8.24 0.12 5062656 1.56 1.52 0.02 8306688
```

### Fingerprints

| Artifact | SHA-256 |
|---|---|
| benchmark harness | `d7c75e785f142764729bbf218d479e8f04ad1d1636a0c15d46e0aef662d06ddc` |
| rsomics-vcf binary | `f8296129fb2c66beccfa855844d9a532a9dedada14fa5b8de165ef89adc224fd` |
| bcftools binary | `33100a6b961c529e915394d53b4737a0f8dd7a164eac352afe4e74e1ced51f60` |
| rustc binary | `d69d40bfd2e11825feb3538512b6ffcd63de91c35ec36bb876849f0f9f8fe6bd` |
| interval target | `46e9f02024f8cc62ab2cf09a6025b4551e675182713fdd3d902b632bf543e15c` |
| interval source | `b3d2fa693fd8e55ed61ce45e89af79662c1bdd37fffb951901255dc13fa3dc10` |
| interval raw distribution | `26a0c93d3347e773dd26237d84e4011397070806c0b283cf40f7659d99948574` |
| interval summary | `92075c45b8f0f3817c049feb6e089c4544783e7a613cc0cee77c321937dc7042` |
| typed target | `185e2082e08b815ada2ec06cb63bf9b087d837f6d5396880bed2ab24f10e7c53` |
| typed source | `c33c8b0741ca5138a3237181f152daab3109fdf352265b2f6112b82c12429448` |
| typed raw distribution | `52f0114a214b170fe240aeca3ebae77ab6be8cf89825ad3f8633c70a99176b4c` |
| typed summary | `5a033cc8dde63013cbe084c5ad076578e506ca85db0c0b808039b0607e9e791c` |
