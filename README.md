# rsomics-vcf

`rsomics-vcf` is the rsomics product for VCF and BCF workflows.

The current implementation provides:

- `head`: stream ordered VCF headers and an optional record prefix from plain
  VCF, BGZF-compressed VCF, BCF, or standard input.
- `filter`: evaluate typed site and sample expressions, annotate failures, and
  apply genomic masks and gap filters.
- `index`: create and inspect CSI indexes for BGZF VCF or BCF and TBI indexes
  for BGZF VCF.
- `query`: project typed site, INFO, FORMAT, and genotype fields from the same
  input formats with sample selection and transactional named output.
- `validate`: validate VCF 4.1–4.5 or BCF 2.2 structure, schema, typed values,
  cardinality, genotypes, and cross-field invariants.
- `view`: convert among VCF, BGZF VCF, raw BCF, and BGZF BCF while selecting
  records, samples, indexed regions, or streaming targets.

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
