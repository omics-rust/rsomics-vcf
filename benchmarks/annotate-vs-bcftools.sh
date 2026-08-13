#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage:" >&2
    echo "  $0 generate RSOMICS_VCF BCFTOOLS RESULT_DIRECTORY" >&2
    echo "  $0 run RSOMICS_VCF BCFTOOLS TARGET ANNOTATIONS COLUMNS RESULT_DIRECTORY" >&2
    exit 2
}

external_directory() {
    case "$1" in
        /Volumes/*) ;;
        *)
            echo "benchmark data and results must be on an external volume" >&2
            exit 2
            ;;
    esac
}

resolve_binary() {
    if [[ "$1" == */* ]]; then
        realpath "$1"
    else
        command -v "$1"
    fi
}

companion() {
    local parent
    parent="$(dirname "$1")"
    for candidate in "$parent/$2" "$parent/htslib-1.24/$2" "$parent/htslib/$2"; do
        if [[ -x "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return
        fi
    done
    command -v "$2"
}

canonicalize() {
    local source="$1"
    local destination="$2"
    local metadata="$destination.metadata"
    local body="$destination.body"
    awk -v metadata="$metadata" -v body="$body" '
        /^##bcftools_/ { next }
        /^##/ { print > metadata; next }
        { print > body }
    ' "$source"
    LC_ALL=C sort "$metadata" > "$destination"
    awk '{ print }' "$body" >> "$destination"
}

measure() {
    local pair="$1"
    local order="$2"
    local tool="$3"
    shift 3
    local timing="$result_directory/$tool-$pair.time"
    /usr/bin/time -lp "$@" > /dev/null 2> "$timing"
    local wall user system rss
    wall="$(awk '$1 == "real" { print $2 }' "$timing")"
    user="$(awk '$1 == "user" { print $2 }' "$timing")"
    system="$(awk '$1 == "sys" { print $2 }' "$timing")"
    rss="$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$pair" "$order" "$tool" "$wall" "$user" "$system" "$rss" >> "$raw_results"
}

summarize() {
    local tool="$1"
    local column="$2"
    local metric="$3"
    local values="$result_directory/$tool-$metric.values"
    awk -F '\t' -v tool="$tool" -v column="$column" \
        '$3 == tool { print $column }' "$raw_results" | LC_ALL=C sort -n > "$values"
    awk -v tool="$tool" -v metric="$metric" '
        { value[NR] = $1; sum += $1; sumsq += $1 * $1 }
        END {
            median = NR % 2 ? value[(NR + 1) / 2] : (value[NR / 2] + value[NR / 2 + 1]) / 2
            p99 = value[int(0.99 * NR + 0.999999)]
            mean = sum / NR
            stdev = NR > 1 ? sqrt((sumsq - sum * sum / NR) / (NR - 1)) : 0
            printf "%s\t%s\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n", tool, metric, NR, median, p99, mean, stdev
        }
    ' "$values" >> "$summary"
}

run_benchmark() {
    [[ $# -eq 6 ]] || usage
    binary="$(resolve_binary "$1")"
    bcftools="$(resolve_binary "$2")"
    rustc="$(resolve_binary "${RSOMICS_RUSTC:-rustc}")"
    target="$(realpath "$3")"
    annotations="$(realpath "$4")"
    columns="$5"
    result_directory="$6"
    external_directory "$target"
    external_directory "$annotations"
    external_directory "$result_directory"
    [[ "$(uname -s)" == "Darwin" ]] || {
        echo "resource measurements are calibrated for macOS /usr/bin/time" >&2
        exit 2
    }
    mkdir -p "$result_directory"
    result_directory="$(realpath "$result_directory")"

    runs="${BENCH_RUNS:-10}"
    warmups="${BENCH_WARMUPS:-3}"
    [[ "$runs" =~ ^[0-9]+$ && "$runs" -ge 10 ]] || {
        echo "BENCH_RUNS must be an integer of at least 10" >&2
        exit 2
    }
    [[ "$warmups" =~ ^[0-9]+$ && "$warmups" -ge 3 ]] || {
        echo "BENCH_WARMUPS must be an integer of at least 3" >&2
        exit 2
    }
    [[ -n "$columns" ]] || usage
    [[ "$($bcftools --version | sed -n '1p')" == "bcftools 1.24" ]] || {
        echo "bcftools 1.24 is required" >&2
        exit 2
    }

    ours_output="$result_directory/rsomics.vcf"
    oracle_output="$result_directory/bcftools.vcf"
    "$binary" annotate -a "$annotations" -c "$columns" "$target" > "$ours_output"
    "$bcftools" annotate --no-version -a "$annotations" -c "$columns" "$target" \
        > "$oracle_output"
    canonicalize "$ours_output" "$result_directory/rsomics.canonical.vcf"
    canonicalize "$oracle_output" "$result_directory/bcftools.canonical.vcf"
    cmp "$result_directory/rsomics.canonical.vcf" "$result_directory/bcftools.canonical.vcf"

    ours_command=("$binary" annotate -a "$annotations" -c "$columns" "$target")
    oracle_command=("$bcftools" annotate --no-version -a "$annotations" -c "$columns" "$target")
    for ((pair = 1; pair <= warmups; pair++)); do
        if ((pair % 2)); then
            "${ours_command[@]}" > /dev/null
            "${oracle_command[@]}" > /dev/null
        else
            "${oracle_command[@]}" > /dev/null
            "${ours_command[@]}" > /dev/null
        fi
    done

    raw_results="$result_directory/raw.tsv"
    printf 'pair\torder\ttool\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\n' \
        > "$raw_results"
    for ((pair = 1; pair <= runs; pair++)); do
        if ((pair % 2)); then
            measure "$pair" 1 rsomics "${ours_command[@]}"
            measure "$pair" 2 bcftools "${oracle_command[@]}"
        else
            measure "$pair" 1 bcftools "${oracle_command[@]}"
            measure "$pair" 2 rsomics "${ours_command[@]}"
        fi
    done

    summary="$result_directory/summary.tsv"
    printf 'tool\tmetric\tn\tmedian\tp99\tmean\tstdev\n' > "$summary"
    for tool in rsomics bcftools; do
        summarize "$tool" 4 wall_seconds
        summarize "$tool" 5 user_seconds
        summarize "$tool" 6 system_seconds
        summarize "$tool" 7 max_rss_bytes
    done

    repository="$(cd "$(dirname "$0")/.." && pwd)"
    {
        printf 'utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        uname -a
        sysctl -n hw.model
        sysctl -n machdep.cpu.brand_string
        printf 'physical_memory_bytes=%s\n' "$(sysctl -n hw.memsize)"
        sw_vers
        "$rustc" --version --verbose
        "$binary" --version
        "$bcftools" --version | sed -n '1,2p'
        printf 'git_head=%s\n' "$(git -C "$repository" rev-parse HEAD)"
        if [[ -z "$(git -C "$repository" status --porcelain)" ]]; then
            echo 'git_dirty=false'
        else
            echo 'git_dirty=true'
        fi
        printf 'runs=%s\nwarmups=%s\ncolumns=%s\n' "$runs" "$warmups" "$columns"
        printf 'rsomics_command='; printf '%q ' "${ours_command[@]}"; printf '\n'
        printf 'bcftools_command='; printf '%q ' "${oracle_command[@]}"; printf '\n'
        wc -c "$target" "$annotations" "$ours_output" "$oracle_output"
        shasum -a 256 "$binary" "$bcftools" "$rustc" "$target" "$annotations" \
            "$ours_output" "$oracle_output" "$result_directory/rsomics.canonical.vcf" \
            "$result_directory/bcftools.canonical.vcf" "$raw_results" "$summary"
    } > "$result_directory/provenance.txt"
}

generate_workloads() {
    [[ $# -eq 3 ]] || usage
    local binary bcftools workspace interval_records typed_records samples
    binary="$(resolve_binary "$1")"
    bcftools="$(resolve_binary "$2")"
    workspace="$3"
    external_directory "$workspace"
    mkdir -p "$workspace/inputs" "$workspace/results"
    workspace="$(realpath "$workspace")"
    interval_records="${INTERVAL_RECORDS:-2000000}"
    typed_records="${TYPED_RECORDS:-300000}"
    samples="${TYPED_SAMPLES:-8}"
    for value in "$interval_records" "$typed_records" "$samples"; do
        [[ "$value" =~ ^[0-9]+$ && "$value" -gt 0 ]] || {
            echo "generated workload sizes must be positive integers" >&2
            exit 2
        }
    done

    local interval_target="$workspace/inputs/interval-target.vcf"
    local interval_source="$workspace/inputs/interval-source.tsv"
    local interval_source_gz="$interval_source.gz"
    awk -v records="$interval_records" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Database value\">"
        print "##contig=<ID=chr1,length=" records * 2 + 10 ">"
        print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        for (i = 1; i <= records; i++)
            print "chr1\t" i * 2 "\t.\tA\tC\t50\tPASS\t."
    }' > "$interval_target"
    awk -v records="$interval_records" 'BEGIN {
        for (i = 1; i <= records; i++)
            print "chr1\t" i * 2 - 1 "\t" i * 2 "\t" i % 1000
    }' > "$interval_source"
    "$(companion "$bcftools" bgzip)" -c "$interval_source" > "$interval_source_gz"
    "$(companion "$bcftools" tabix)" -f -s 1 -b 2 -e 3 "$interval_source_gz"

    local typed_target_source="$workspace/inputs/typed-target.vcf"
    local typed_target="$typed_target_source.gz"
    local typed_source="$workspace/inputs/typed-source.vcf"
    local typed_source_gz="$typed_source.gz"
    awk -v records="$typed_records" -v samples="$samples" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">"
        print "##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">"
        print "##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">"
        print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihood\">"
        print "##contig=<ID=chr1,length=" records + 10 ">"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (sample = 1; sample <= samples; sample++) printf "\tS%d", sample
        print ""
        for (position = 1; position <= records; position++) {
            printf "chr1\t%d\t.\tA\tC,G\t50\tPASS\tIA=1,2;IR=3,4,5;IG=0,1,2,3,4,5\tGT:AD:PL", position
            for (sample = 1; sample <= samples; sample++)
                printf "\t1/2:10,4,6:60,40,20,50,0,30"
            print ""
        }
    }' > "$typed_target_source"
    awk -v records="$typed_records" -v samples="$samples" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##INFO=<ID=IA,Number=A,Type=Integer,Description=\"A\">"
        print "##INFO=<ID=IR,Number=R,Type=Integer,Description=\"R\">"
        print "##INFO=<ID=IG,Number=G,Type=Integer,Description=\"G\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depth\">"
        print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihood\">"
        print "##contig=<ID=chr1,length=" records + 10 ">"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (sample = samples; sample >= 1; sample--) printf "\tS%d", sample
        print ""
        for (position = 1; position <= records; position++) {
            printf "chr1\t%d\ts%d\tA\tC,G\t60\tPASS\tIA=10,20;IR=30,40,50;IG=0,10,20,30,40,50\tGT:AD:PL", position, position
            for (sample = samples; sample >= 1; sample--)
                printf "\t1|2:30,40,50:0,10,20,30,40,50"
            print ""
        }
    }' > "$typed_source"
    "$bcftools" view --no-version -Oz -o "$typed_target" "$typed_target_source"
    "$bcftools" index -f "$typed_target"
    "$bcftools" view --no-version -Oz -o "$typed_source_gz" "$typed_source"
    "$bcftools" index -f "$typed_source_gz"

    {
        printf 'interval_records=%s\ntyped_records=%s\ntyped_samples=%s\n' \
            "$interval_records" "$typed_records" "$samples"
        printf 'generator='; printf '%q ' "$0" generate "$binary" "$bcftools" "$workspace"; printf '\n'
        shasum -a 256 "$(realpath "$0")" "$interval_target" "$interval_source" \
            "$interval_source_gz" "$interval_source_gz.tbi" "$typed_target_source" \
            "$typed_target" "$typed_target.csi" "$typed_source" "$typed_source_gz" \
            "$typed_source_gz.csi"
    } > "$workspace/generation.txt"

    BENCH_RUNS="${BENCH_RUNS:-10}" BENCH_WARMUPS="${BENCH_WARMUPS:-3}" \
        "$0" run "$binary" "$bcftools" "$interval_target" "$interval_source_gz" \
        'CHROM,FROM,TO,INFO/SCORE' "$workspace/results/interval"
    BENCH_RUNS="${BENCH_RUNS:-10}" BENCH_WARMUPS="${BENCH_WARMUPS:-3}" \
        "$0" run "$binary" "$bcftools" "$typed_target" "$typed_source_gz" \
        'INFO/IA,INFO/IR,INFO/IG,FORMAT/GT,FORMAT/AD,FORMAT/PL' \
        "$workspace/results/typed"
}

[[ $# -gt 0 ]] || usage
mode="$1"
shift
case "$mode" in
    generate) generate_workloads "$@" ;;
    run) run_benchmark "$@" ;;
    *) usage ;;
esac
