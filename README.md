# rsomics-vcf

`rsomics-vcf` is the rsomics product for VCF and BCF workflows.

The current implementation provides:

- `annotate`: edit headers and records or transfer typed fields from sorted
  VCF, BCF, BED, and tab-delimited annotation streams.
- `head`: stream ordered VCF headers and an optional record prefix from plain
  VCF, BGZF-compressed VCF, BCF, or standard input.
- `filter`: evaluate typed site and sample expressions, annotate failures, and
  apply genomic masks and gap filters.
- `index`: create and inspect CSI indexes for BGZF VCF or BCF and TBI indexes
  for BGZF VCF.
- `norm`: normalize against an indexed reference, split, join, atomize, and
  deduplicate variants while remapping typed INFO, FORMAT, and genotype fields.
- `query`: project typed site, INFO, FORMAT, and genotype fields from the same
  input formats with sample selection and transactional named output.
- `reheader`: replace headers, synchronize contigs from a FASTA index, or
  rename samples without changing the variant payload.
- `setgt`: select and replace typed genotypes across VCF and BCF encodings.
- `validate`: validate VCF 4.1–4.5 or BCF 2.2 structure, schema, typed values,
  cardinality, genotypes, and cross-field invariants.
- `view`: convert among VCF, BGZF VCF, raw BCF, and BGZF BCF while selecting
  records, samples, indexed regions, or streaming targets.

`annotate` applies checked header removals and renames, site-format ID edits,
typed fixed/INFO/FORMAT/GT transfer, explicit sample mapping, allele-aware
`Number=A`, `R`, and `G` remapping, target expressions, site marking, and
indexed regions. The annotation join advances over two coordinate-sorted
streams and retains only records that can still overlap the current target.
Named output is transactional; all four VCF/BCF encodings and bounded BGZF
compression workers use the same writer layer as the other commands.

```console
rsomics-vcf annotate calls.vcf.gz -a db.vcf.gz -c ID,INFO/AF
rsomics-vcf annotate calls.bcf -a depths.vcf.gz -c FORMAT/DP -s tumor,normal
rsomics-vcf annotate calls.vcf.gz -x INFO/OLD --rename-annotations names.tsv
```

Whole FORMAT transfer excludes GT unless `FORMAT/GT` is requested explicitly.
Missing definitions, incompatible schemas, invalid cardinalities, unavailable
samples, coordinate regressions, and impossible allele maps fail nonzero.
Experimental merge logic, dynamic source-column expressions, forced recovery,
automatic indexing, and provenance stamping are not accepted as placeholder
options.

The annotation release benchmark verifies complete normalized output against
bcftools 1.24 before timing. On an Apple M2 with three warmups and ten
alternating measured pairs, the 2,000,000-record interval join used 4.8 MB
median peak RSS versus 124.4 MB for bcftools while taking 3.895 versus 3.735
seconds. The 300,000-record, eight-sample typed transfer used 5.1 versus 8.3 MB
RSS but took 8.495 versus 1.560 seconds. This is a bounded-memory advantage,
not a general throughput claim; the complete distributions and fingerprints
are in `PERFORMANCE.md`.

`reheader` composes complete header replacement, FAI contig synchronization,
and positional or old-to-new sample renaming in that order. It preserves the
input encoding across plain VCF, BGZF VCF, raw BCF, and BGZF BCF. Named output
is atomic, and every requested edit is validated before commit.

```console
rsomics-vcf reheader -H header.vcfh -o renamed.vcf calls.vcf
rsomics-vcf reheader -f reference.fa.fai -N samples.tsv -o renamed.bcf calls.bcf
rsomics-vcf reheader -n tumor,normal calls.vcf.gz > renamed.vcf.gz
```

Plain VCF record bytes are copied unchanged. BGZF VCF rewrites only header
frames, preserves the remaining compressed record frames byte for byte, and
requires one canonical EOF block. BCF records are streamed through stable
edited dictionaries so removed definitions still referenced by a record fail
nonzero. Sample count mismatches, unknown mapping sources, duplicate final
names, malformed FAI rows, ordinary gzip input, truncated data, and unsafe
output aliases also fail instead of producing a partial artifact. Nonzero
`--threads` is accepted only for BGZF BCF.

The 0.5 release benchmark compares complete output with bcftools 1.24 before
timing. On a 2,000,000-record Apple M2 workload, plain VCF reheadering used
0.515 versus 3.075 seconds median wall time and 3.92 versus 6.39 MB median peak
RSS. BGZF reheadering used 4.64 versus 7.01 MB median peak RSS but took 0.040
versus 0.020 seconds. The complete distributions, equality hashes, commands,
and artifact fingerprints are in `PERFORMANCE.md`.

`setgt` edits genotypes through the product's typed record and expression
layers. Targets are `.`, `./x`, `./.`, `a`, `q`, `b:TAG<CMP>VALUE`, and
`r:FLOAT`; one principal target can compose with one random fraction.
Replacements are `.`, `0`, `0p`, `m`, `mp`, `M`, `Mp`, `X`, `p`, `u`, `i`,
and `c:GT`. Query targets require exactly one include or exclude expression,
and random selection can use a signed deterministic seed.

```console
rsomics-vcf setgt -t . -n 0 calls.vcf.gz
rsomics-vcf setgt -t q -i 'FMT/DP < 10' -n . -O z -o edited.vcf.gz calls.bcf
rsomics-vcf setgt -t a -t r:0.1 --seed 7 -n c:0/1 calls.vcf
```

Plain VCF, BGZF VCF, raw BCF, and BGZF BCF share the same edit engine. Named
output is transactional, and compression workers are accepted only for BGZF
outputs. Existing valid `INFO/AC` and `INFO/AN` values are recomputed after a
genotype change; absent tags remain absent, while malformed definitions or
values fail instead of being preserved stale.

Compatibility is defined against bcftools 1.24 except where failing loud or
typed per-sample behavior is safer. Ambiguous replacement spellings and
missing required `FORMAT/AD` values are errors. Query inversion edits the
selected sample, and partial-missing, binomial, and inversion decisions use
each sample's actual ploidy rather than a record-wide encoded width. A late
record error cannot truncate an existing named destination.

The 0.6 release benchmark validates all three operations across VCF, BGZF VCF,
and BCF after every measured command. On the 2,000,000-record, eight-sample
Apple M2 workload, median peak RSS is 28.50% to 37.02% below bcftools 1.24.
Median wall time is 3.20 to 16.60 times slower, so this is a bounded-memory
claim rather than a throughput claim. Complete distributions and fingerprints
are in `PERFORMANCE.md`.

`head` writes VCF text, preserves header order, normalizes the standard PASS
definition, removes BCF-internal `IDX` fields, and renders typed records
compatibly with bcftools 1.24. Unlike bcftools, malformed records and invalid
numeric arguments fail with a nonzero exit code.

```console
rsomics-vcf head variants.bcf
rsomics-vcf head -n 10 variants.vcf.gz
rsomics-vcf head -s 3 < variants.bcf
```

`-H` limits metadata header lines while retaining `-h` for the unified rsomics
help experience. `-s` starts output at `#CHROM` and then emits the requested
number of records.

`index` creates CSI by default, supports custom CSI minimum shifts and BGZF
decompression workers, and writes the completed index transactionally. It
rejects ordinary gzip, missing BGZF terminators, unsorted coordinates,
noncontiguous contig blocks, malformed records, and implicit replacement.
`--stats` and `--nrecords` read count metadata from a variant path or directly
from an existing index.

```console
rsomics-vcf index variants.vcf.gz
rsomics-vcf index --tbi variants.vcf.gz
rsomics-vcf index --min-shift 18 variants.bcf
rsomics-vcf index --stats variants.vcf.gz
rsomics-vcf index --nrecords variants.bcf.csi
```

`query` supports fixed columns, `%POS0`, `%END`, `%END0`, `%FIRST_ALT`,
`%TYPE`, `%INFO`, `%INFO/TAG`, `%FORMAT`, `%LINE`, array subscripts, and
per-sample loops with `%SAMPLE`, `%GT`, `%TGT`, `%IUPACGT`, and declared FORMAT
tags. `-H` and `-HH` produce indexed and plain projection headers. If a format
contains no newline, one is appended per record unless `-N` is used.

```console
rsomics-vcf query variants.bcf -f '%CHROM\t%POS\t%INFO/DP\n'
rsomics-vcf query variants.vcf.gz -s S1,S2 -f '%POS[\t%SAMPLE=%GT]\n'
rsomics-vcf query variants.bcf -f '%CHROM\t%POS0\t%END\n' -o variants.tsv
```

The current query contract is a single-input field projection. Region and
target selection, expression filtering and functions, multi-input `%MASK`,
`%PBINOM`, `%N_PASS`, `%TBCSQ`, `%VKX`, and undefined-tag fallback are outside
this contract.

`view` preserves typed records across all four encodings. It supports ordered
sample inclusion and exclusion, optional genotype removal, AC/AN recalculation,
FILTER, ID, allele-count, and current bcftools variant-type selection. Targets
stream through any supported input; regions require a CSI or TBI and use
position, record, or variant overlap semantics.

```console
rsomics-vcf view variants.bcf -O z -o variants.vcf.gz
rsomics-vcf view variants.vcf.gz -s tumor,normal -r chr1:1-100000
rsomics-vcf view variants.bcf -v snps,indels -f PASS -G
```

Named output is transactional, and JSON summaries require it so variant data
never shares standard output with the envelope. Expression filtering, allele
trimming and remapping, frequency/genotype predicates, output indexing, and
compression workers are outside this first stable `view` contract; they are
not accepted as partially implemented flags.

`norm` accepts every VCF and BCF input encoding supported by the product. It
left-aligns and trims alleles against indexed plain or BGZF FASTA, validates or
repairs REF mismatches, and can use GFF3 transcript orientation for HGVS
3-prime right alignment. Multiallelic splitting and joining preserve declared
`Number=A`, `R`, and `G` INFO and FORMAT values, genotypes, phasing, and mixed
ploidy. Complex variants can be atomized with explicit overlap and provenance
policies. Duplicate removal, expression-controlled transformation, indexed
regions, streaming targets, position or lexical local sorting, transactional
output, JSON summaries, and bounded BGZF compression workers compose with the
same pipeline.

```console
rsomics-vcf norm -f reference.fa calls.vcf.gz -O z -o normalized.vcf.gz
rsomics-vcf norm -m --keep-sum AD --split-overlaps missing calls.bcf
rsomics-vcf norm --join-multiallelic any --strict-filter calls.vcf
rsomics-vcf norm -f reference.fa -g transcripts.gff3.gz calls.vcf
```

Malformed allele-indexed fields fail instead of being silently discarded.
Automatic output indexing is intentionally excluded until the output and index
can be committed as one transaction.

The normalization release benchmark compares exact record bodies with
bcftools 1.24 before timing. It ran on an Apple M2 Mac mini with 8 GB RAM and
macOS 26.6.1, using three warmups and ten measured runs:

| Workload | Records | rsomics-vcf | bcftools | Throughput result |
|---|---:|---:|---:|---:|
| Reference-guided indel alignment | 500,000 | 0.920 ± 0.008 s | 1.170 ± 0.017 s | 1.27× faster |
| Typed eight-sample multiallelic split | 200,000 input, 400,000 output | 5.457 ± 0.061 s | 7.456 ± 0.047 s | 1.37× faster |

A single `/usr/bin/time -lp` run reported 6,799,360 versus 7,012,352 bytes RSS
for reference alignment, and 10,518,528 versus 8,437,760 bytes for typed
splitting. The committed `benchmarks/norm-vs-bcftools.sh` harness generates the
fixtures, verifies body equality, records hashes and machine provenance, and
captures the complete timing distributions.

`filter` evaluates the product's typed expression language over fixed, INFO,
FORMAT, genotype, calculated, arithmetic, logical, regex, file-set, and
statistical function values. It supports hard include or exclude selection,
soft FILTER annotation modes, failed-sample genotype replacement, masks,
SnpGap and IndelGap, indexed regions, streaming targets, every VCF/BCF output
encoding, and bounded BGZF compression workers. Named outputs use the shared
`rsomics-common` transaction and JSON summaries remain separate from variant
data.

```console
rsomics-vcf filter variants.vcf.gz -i 'QUAL >= 20 && INFO/DP >= 10'
rsomics-vcf filter variants.bcf -e 'FMT/DP < 10' -s LowDepth -S .
rsomics-vcf filter variants.vcf.gz -g 3:indel,mnp -G 5 -O z --threads 4 -o filtered.vcf.gz
```

The release benchmark uses bcftools 1.24 as the oracle on an Apple M2 Mac mini
with 8 GB RAM and macOS 26.6.1. On a 250,000-record, 7.5 MiB plain VCF
(`f3352403c8a071f09e71da3b38901bca091a73fdd6f7b0b2a49281095964512a`),
filtering `INFO/DP >= 20 && QUAL >= 30` to `/dev/null` took 57.16–59.24 ms
with a 58.24 ms estimate, versus 111.94–118.35 ms and a 114.91 ms estimate
for bcftools. Three warm `/usr/bin/time -l` runs used at most 3,702,784 bytes
RSS for rsomics and 6,782,976 bytes for bcftools. The committed Criterion
benchmark generates the input and records both implementations under the same
conditions.

`validate` accepts plain or gzip/BGZF-compressed VCF, raw or BGZF-compressed
BCF, and standard input. Diagnostics identify the record line and field, and
invalid input exits with status 1. By default it retains at most 100
diagnostics while still reporting the complete error and warning counts.

```console
rsomics-vcf validate variants.vcf.gz
rsomics-vcf validate variants.bcf --max-diagnostics 20
rsomics-vcf validate variants.vcf --require-evidence --json
```

`--require-evidence` additionally requires every record to carry `GT`,
`INFO/AF`, or both `INFO/AC` and `INFO/AN`. This is an explicit policy rather
than part of the VCF format validity contract.

## License

Licensed under either the Apache License, Version 2.0 or the MIT License.
