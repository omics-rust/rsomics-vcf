# bcftools reheader fixtures

These compact fixtures are team-owned inputs written for the `bcftools 1.24`
live compatibility gate. They are not copied from the bcftools source tree.
The gate derives expected outputs by executing the pinned oracle rather than
committing generated golden files.

| File | Purpose | SHA-256 |
|---|---|---|
| `input.vcf` | Two-contig, two-sample typed input | `33c2b982b6e5fb8a2c3383e24a637a8b8ffe7a7cb6d4c88a03efb54b28dbc37b` |
| `replacement.vcfh` | Complete replacement header | `9530e76da8f64ebaabc0dfc104d12c26e536e23cf1bc6219b01a7c5a6c951c0d` |
| `reference.fai` | Retained and appended contigs | `f83766fe8d27a52105d82be83bed7d73d18125d064c575b74a33080dba51db7f` |
| `samples.txt` | Old-to-new sample mapping | `33954f9ef8c6b87808be4612053e39165eb45a040dc4fd0500483a6de1b10056` |

VCF outputs are compared as complete uncompressed text. BCF outputs are
decoded by bcftools, parsed as semantic headers, and compared with their typed
record text. Safety divergences are checked only against the rsomics contract;
oracle warnings or crashes are never treated as expected success.
