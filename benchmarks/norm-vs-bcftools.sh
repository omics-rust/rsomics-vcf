#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 /Volumes/.../benchmark-directory" >&2
    exit 2
fi

scratch="$1"
case "$scratch" in
    /Volumes/*) ;;
    *)
        echo "benchmark output must be on an external volume" >&2
        exit 2
        ;;
esac

binary="${RSOMICS_VCF_BIN:-$(pwd)/target/release/rsomics-vcf}"
bcftools="${BCFTOOLS_BIN:-bcftools}"
binary="$(realpath "$binary")"
bcftools="$(command -v "$bcftools")"
reference_records="${REFERENCE_RECORDS:-500000}"
split_records="${SPLIT_RECORDS:-200000}"
runs="${BENCH_RUNS:-10}"
warmup="${BENCH_WARMUP:-3}"

mkdir -p "$scratch"
reference="$scratch/reference.fa"
reference_input="$scratch/reference-$reference_records.vcf"
split_input="$scratch/split-$split_records.vcf"

awk -v records="$reference_records" 'BEGIN {
    motif = "CAAAAAAG"
    sequence_length = records * 8 + 200
    print ">chr1"
    for (i = 0; i < sequence_length; i++) {
        printf "%s", substr(motif, i % 8 + 1, 1)
        if (i % 80 == 79) print ""
    }
    if (sequence_length % 80) print ""
}' > "$reference"

reference_length=$((reference_records * 8 + 200))
printf 'chr1\t%s\t6\t80\t81\n' "$reference_length" > "$reference.fai"

awk -v records="$reference_records" -v length="$reference_length" 'BEGIN {
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=chr1,length=" length ">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    for (i = 0; i < records; i++) {
        print "chr1\t" i * 8 + 4 "\t.\tAA\tA\t50\tPASS\t."
    }
}' > "$reference_input"

awk -v records="$split_records" 'BEGIN {
    print "##fileformat=VCFv4.3"
    print "##contig=<ID=chr1,length=" records + 1 ">"
    print "##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">"
    print "##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">"
    print "##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">"
    print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihood\">"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (sample = 1; sample <= 8; sample++) printf "\tS%d", sample
    print ""
    value = "1/2:10,4,6:60,40,20,50,0,30"
    for (position = 1; position <= records; position++) {
        printf "chr1\t%d\t.\tA\tC,G\t50\tPASS\tIA=10,20;IR=5,3,2;IG=0,10,20,30,40,50\tGT:AD:PL", position
        for (sample = 1; sample <= 8; sample++) printf "\t%s", value
        print ""
    }
}' > "$split_input"

"$binary" --version > "$scratch/rsomics-version.txt"
"$bcftools" --version > "$scratch/bcftools-version.txt"

"$binary" norm --fasta-ref "$reference" "$reference_input" > "$scratch/reference-rsomics.vcf"
"$bcftools" norm --no-version --fasta-ref "$reference" "$reference_input" > "$scratch/reference-bcftools.vcf"
cmp <(sed '/^#/d' "$scratch/reference-rsomics.vcf") <(sed '/^#/d' "$scratch/reference-bcftools.vcf")

"$binary" norm --split-multiallelic "$split_input" > "$scratch/split-rsomics.vcf"
"$bcftools" norm --no-version -m -any "$split_input" > "$scratch/split-bcftools.vcf"
cmp <(sed '/^#/d' "$scratch/split-rsomics.vcf") <(sed '/^#/d' "$scratch/split-bcftools.vcf")
sed '/^#/d' "$scratch/reference-rsomics.vcf" | shasum -a 256 > "$scratch/reference-body.sha256"
sed '/^#/d' "$scratch/split-rsomics.vcf" | shasum -a 256 > "$scratch/split-body.sha256"

export RSOMICS_VCF_BIN="$binary"
export BCFTOOLS_BIN="$bcftools"
export REFERENCE_FASTA="$reference"
export REFERENCE_INPUT="$reference_input"
export SPLIT_INPUT="$split_input"

hyperfine \
    --warmup "$warmup" \
    --runs "$runs" \
    --export-json "$scratch/reference-hyperfine.json" \
    '"$RSOMICS_VCF_BIN" norm --fasta-ref "$REFERENCE_FASTA" "$REFERENCE_INPUT" > /dev/null' \
    '"$BCFTOOLS_BIN" norm --no-version --fasta-ref "$REFERENCE_FASTA" "$REFERENCE_INPUT" > /dev/null'

hyperfine \
    --warmup "$warmup" \
    --runs "$runs" \
    --export-json "$scratch/split-hyperfine.json" \
    '"$RSOMICS_VCF_BIN" norm --split-multiallelic "$SPLIT_INPUT" > /dev/null' \
    '"$BCFTOOLS_BIN" norm --no-version -m -any "$SPLIT_INPUT" > /dev/null'

/usr/bin/time -lp "$binary" norm --fasta-ref "$reference" "$reference_input" \
    > /dev/null 2> "$scratch/reference-rsomics-resource.txt"
/usr/bin/time -lp "$bcftools" norm --no-version --fasta-ref "$reference" "$reference_input" \
    > /dev/null 2> "$scratch/reference-bcftools-resource.txt"
/usr/bin/time -lp "$binary" norm --split-multiallelic "$split_input" \
    > /dev/null 2> "$scratch/split-rsomics-resource.txt"
/usr/bin/time -lp "$bcftools" norm --no-version -m -any "$split_input" \
    > /dev/null 2> "$scratch/split-bcftools-resource.txt"

{
    date -u '+utc=%Y-%m-%dT%H:%M:%SZ'
    uname -a
    sysctl -n hw.model
    sw_vers
    "$binary" --version
    "$bcftools" --version | sed -n '1,2p'
    hyperfine --version
    printf 'git_head=%s\n' "$(git rev-parse HEAD)"
    if git diff --quiet && git diff --cached --quiet; then
        echo 'git_dirty=false'
    else
        echo 'git_dirty=true'
    fi
    printf 'reference_records=%s\nsplit_records=%s\nruns=%s\nwarmup=%s\n' \
        "$reference_records" "$split_records" "$runs" "$warmup"
    printf 'reference_input_bytes=%s\n' "$(wc -c < "$reference_input" | tr -d ' ')"
    printf 'split_input_bytes=%s\n' "$(wc -c < "$split_input" | tr -d ' ')"
    shasum -a 256 "$binary" "$bcftools"
} > "$scratch/provenance.txt"

shasum -a 256 \
    "$reference" \
    "$reference.fai" \
    "$reference_input" \
    "$split_input" \
    "$scratch/reference-rsomics.vcf" \
    "$scratch/reference-bcftools.vcf" \
    "$scratch/split-rsomics.vcf" \
    "$scratch/split-bcftools.vcf" \
    "$scratch/reference-body.sha256" \
    "$scratch/split-body.sha256" \
    "$scratch/reference-hyperfine.json" \
    "$scratch/split-hyperfine.json" \
    > "$scratch/sha256.txt"
