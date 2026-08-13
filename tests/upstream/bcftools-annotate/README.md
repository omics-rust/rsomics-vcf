# bcftools annotate fixtures

The files in the table below are unchanged copies from the bcftools 1.24
release archive, `test/`. The archive is pinned by SHA-256
`8caddc22610ee2851666047c859bb91da0c1e32d0c2ec553db6f153ad130e46f`.
They cover typed missing-value policies, allele-indexed values, VCF sample
mapping, site marking, and annotation renames.

| File | SHA-256 |
|---|---|
| `annotate.missing.vcf` | `7b0cb9ff1bc9f30278e23b48c1e4334bbe78b03dd357dbc4b04a37b60ae255fa` |
| `annotate.missing.tab` | `952df42468c5b086e70e0b856a9a830b78c69644e0af3c60521d7ee296fe57fd` |
| `annotate.AR.vcf` | `ab824f6e354f899a81d3ea74bc39c02d379e266dedee26b87efcb0ea714ff4a4` |
| `annotate.AR.tab` | `4f249b178c3f12a451fbcbd7350125e90df321de0bd9087490c28fb4f9ed96ce` |
| `annotate2.vcf` | `234be8671513276c22195c8283d7c4c7d3151205aded18633b8a05bf5c52d7d8` |
| `annots2.vcf` | `d96498afc7c0277c237f42432c25d8dedaddc377362879589f1e717c08235081` |
| `annots-mark.vcf` | `4610ab52b0801b1ec3b4353001c6a619fdd5fa0c39575c84a98364bb1ea1ccd5` |
| `annotate21.vcf` | `4342153fe76808d96924eb2e7d10f198dd5657e822f153038a0d77533b763a15` |
| `annotate21.txt` | `62f763f7d3471006df454a86b498ade7d32b22626146ab556ef8c113468805f0` |

The behavior reference is `vcfannotate.c` from the same release, SHA-256
`6be47073e1d549f2bcded27f4cf8952ccd03f90ad537088f418dfe8a5d645730`.
bcftools is copyright Genome Research Ltd. and offered under the MIT/Expat
license or GPL; these fixtures and the behavior study use the MIT/Expat
option recorded in `THIRD_PARTY_LICENSES.md`.

`regions.bed` predates this import. It is a team-owned compact BED regression
asset retained from the historical `rsomics-vcf-annotate` source pool and is
not represented as an upstream bcftools file.
