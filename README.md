# rsomics-vcf

`rsomics-vcf` is the rsomics product for VCF and BCF workflows.

The current implementation provides:

- `head`: stream ordered VCF headers and an optional record prefix from plain
  VCF, BGZF-compressed VCF, BCF, or standard input.

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

## License

Licensed under either the Apache License, Version 2.0 or the MIT License.
