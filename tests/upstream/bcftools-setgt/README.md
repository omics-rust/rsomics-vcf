# bcftools setGT fixtures

These compact team-owned fixtures drive the live compatibility gate against
the pinned bcftools and HTSlib 1.24 build. They are not copied from the
bcftools source tree.

`core.vcf` omits genotype-derived INFO tags so ordinary genotype behavior can
be compared independently of rsomics reconciliation policy. It covers
biallelic and multiallelic records, haploid, diploid, and polyploid calls,
phased and unphased calls, and partial and complete missingness.

`malformed-ad.vcf` declares FORMAT/AD but omits it from the record. It records
the oracle's silent no-op separately from the rsomics fail-loud contract.

The gate generates all four input encodings and derives expected output by
executing the live oracle. Comparisons remove only tool provenance header
lines and normalize each result through bcftools `view --no-version -Ov`.
