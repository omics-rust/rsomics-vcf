# rsomics-vcf

`rsomics-vcf` is the rsomics product for VCF and BCF workflows.

The current implementation provides:

- `head`: stream ordered VCF headers and an optional record prefix from plain
  VCF, BGZF-compressed VCF, BCF, or standard input.
- `index`: create and inspect CSI indexes for BGZF VCF or BCF and TBI indexes
  for BGZF VCF.
- `query`: project typed site, INFO, FORMAT, and genotype fields from the same
  input formats with sample selection and transactional named output.
- `validate`: validate VCF 4.1–4.5 or BCF 2.2 structure, schema, typed values,
  cardinality, genotypes, and cross-field invariants.

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
