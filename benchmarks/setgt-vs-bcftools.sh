#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<EOF
usage: $0 --records INT --samples INT --warmups INT --runs INT --results DIRECTORY
          [--binary FILE] [--bcftools FILE]
EOF
    exit 2
}

resolve_binary() {
    if [[ "$1" == */* ]]; then
        realpath "$1"
    else
        command -v "$1"
    fi
}

require_external() {
    case "$1" in
        /Volumes/*) ;;
        *)
            echo "benchmark scratch and results must be on an external volume" >&2
            exit 2
            ;;
    esac
}

require_positive_integer() {
    [[ "$2" =~ ^[0-9]+$ && "$2" -gt 0 ]] || {
        echo "$1 must be a positive integer" >&2
        exit 2
    }
}

require_nonnegative_integer() {
    [[ "$2" =~ ^[0-9]+$ ]] || {
        echo "$1 must be a nonnegative integer" >&2
        exit 2
    }
}

canonical_hash() {
    "$bcftools" view --no-version -Ov "$1" \
        | ruby -e '
            metadata = []
            emitted = false
            ARGF.each_line do |line|
              next if line.start_with?("##bcftools_")
              if line.start_with?("##")
                metadata << line
              else
                unless emitted
                  metadata.sort.each { |entry| print entry }
                  emitted = true
                end
                print line
              end
            end
            metadata.sort.each { |entry| print entry } unless emitted
          ' \
        | shasum -a 256 \
        | awk '{ print $1 }'
}

output_flag() {
    case "$1" in
        vcf) printf 'v' ;;
        bgzf) printf 'z' ;;
        bcf) printf 'b' ;;
        *) return 1 ;;
    esac
}

oracle_output_flag() {
    case "$1" in
        vcf) printf '%s' '-Ov' ;;
        bgzf) printf '%s' '-Oz' ;;
        bcf) printf '%s' '-Ob' ;;
        *) return 1 ;;
    esac
}

input_path() {
    case "$1" in
        vcf) printf '%s' "$scratch/input.vcf" ;;
        bgzf) printf '%s' "$scratch/input.vcf.gz" ;;
        bcf) printf '%s' "$scratch/input.bcf" ;;
        *) return 1 ;;
    esac
}

set_operation_arguments() {
    case "$1" in
        all_to_missing) operation_arguments=(-t a -n .) ;;
        missing_to_reference) operation_arguments=(-t . -n 0) ;;
        query_to_reference) operation_arguments=(-t q -i 'FMT/DP < 10' -n 0) ;;
        *) return 1 ;;
    esac
}

build_command() {
    local tool="$1"
    local operation="$2"
    local format="$3"
    local input="$4"
    local output="$5"
    set_operation_arguments "$operation"
    if [[ "$tool" == rsomics ]]; then
        command=("$binary" setgt "${operation_arguments[@]}" -O "$(output_flag "$format")" -o "$output" "$input")
    else
        command=("$bcftools" +setGT "$input" --no-version "$(oracle_output_flag "$format")" -o "$output" -- "${operation_arguments[@]}")
    fi
}

generate_inputs() {
    awk -v records="$records" -v samples="$samples" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##source=setgt-benchmark"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
        print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">"
        print "##contig=<ID=chr1,length=" records + 10 ">"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (sample = 1; sample <= samples; sample++) printf "\tS%d", sample
        print ""
        for (position = 1; position <= records; position++) {
            printf "chr1\t%d\t.\tA\tC\t50\tPASS\t.\tGT:AD:DP", position
            for (sample = 1; sample <= samples; sample++) {
                kind = (position + sample) % 5
                if (kind == 0) printf "\t./.:.,.:5"
                else if (kind == 1) printf "\t0/1:6,6:12"
                else if (kind == 2) printf "\t1/1:0,6:6"
                else if (kind == 3) printf "\t0/0:20,0:20"
                else printf "\t0/1:2,2:4"
            }
            print ""
        }
    }' > "$scratch/input.vcf"
    "$bcftools" view --no-version -Oz -o "$scratch/input.vcf.gz" "$scratch/input.vcf"
    "$bcftools" view --no-version -Ob -o "$scratch/input.bcf" "$scratch/input.vcf"
    {
        printf 'records=%s\n' "$records"
        printf 'samples=%s\n' "$samples"
        wc -c "$scratch"/input.vcf "$scratch"/input.vcf.gz "$scratch"/input.bcf
        shasum -a 256 "$scratch"/input.vcf "$scratch"/input.vcf.gz "$scratch"/input.bcf
    } > "$result_directory/generation.txt"
}

run_and_hash() {
    local tool="$1"
    local operation="$2"
    local format="$3"
    local output="$4"
    local input
    input="$(input_path "$format")"
    rm -f "$output"
    build_command "$tool" "$operation" "$format" "$input" "$output"
    "${command[@]}"
    canonical_hash "$output"
}

measure() {
    local case_name="$1"
    local pair="$2"
    local order="$3"
    local tool="$4"
    local operation="$5"
    local format="$6"
    local output="$scratch/measured-$tool.$format"
    local input timing wall user system rss semantic
    input="$(input_path "$format")"
    timing="$scratch/$case_name-$tool-$pair.time"
    rm -f "$output" "$timing"
    build_command "$tool" "$operation" "$format" "$input" "$output"
    /usr/bin/time -lp "${command[@]}" > /dev/null 2> "$timing"
    wall="$(awk '$1 == "real" { print $2 }' "$timing")"
    user="$(awk '$1 == "user" { print $2 }' "$timing")"
    system="$(awk '$1 == "sys" { print $2 }' "$timing")"
    rss="$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")"
    [[ -n "$wall" && -n "$user" && -n "$system" && -n "$rss" ]]
    semantic="$(canonical_hash "$output")"
    [[ "$semantic" == "$(<"$scratch/expected-$case_name.sha256")" ]] || {
        echo "$case_name $tool pair $pair changed semantics" >&2
        exit 1
    }
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_name" "$pair" "$order" "$tool" "$wall" "$user" "$system" "$rss" "$semantic" \
        >> "$raw_results"
    rm -f "$output" "$timing"
}

summarize() {
    local case_name="$1"
    local tool="$2"
    local column="$3"
    local metric="$4"
    awk -F '\t' -v case_name="$case_name" -v tool="$tool" -v column="$column" \
        '$1 == case_name && $4 == tool { print $column }' "$raw_results" \
        | LC_ALL=C sort -n \
        | awk -v case_name="$case_name" -v tool="$tool" -v metric="$metric" '
            { value[NR] = $1; sum += $1; sumsq += $1 * $1 }
            END {
                if (NR == 0) exit 1
                median = NR % 2 ? value[(NR + 1) / 2] : (value[NR / 2] + value[NR / 2 + 1]) / 2
                p99 = value[int(0.99 * NR + 0.999999)]
                mean = sum / NR
                stdev = NR > 1 ? sqrt((sumsq - sum * sum / NR) / (NR - 1)) : 0
                printf "%s\t%s\t%s\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n", case_name, tool, metric, NR, median, p99, mean, stdev
            }
        ' >> "$summary"
}

median_value() {
    awk -F '\t' -v case_name="$1" -v tool="$2" -v metric="$3" \
        '$1 == case_name && $2 == tool && $3 == metric { print $5 }' "$summary"
}

repository="$(cd "$(dirname "$0")/.." && pwd)"
records=2000000
samples=8
warmups=3
runs=10
result_directory=
binary="${CARGO_TARGET_DIR:-$repository/target}/release/rsomics-vcf"
bcftools=bcftools

while [[ $# -gt 0 ]]; do
    case "$1" in
        --records) records="${2:-}"; shift 2 ;;
        --samples) samples="${2:-}"; shift 2 ;;
        --warmups) warmups="${2:-}"; shift 2 ;;
        --runs) runs="${2:-}"; shift 2 ;;
        --results) result_directory="${2:-}"; shift 2 ;;
        --binary) binary="${2:-}"; shift 2 ;;
        --bcftools) bcftools="${2:-}"; shift 2 ;;
        -h|--help) usage ;;
        *) usage ;;
    esac
done

[[ -n "$result_directory" ]] || usage
require_positive_integer records "$records"
require_positive_integer samples "$samples"
require_nonnegative_integer warmups "$warmups"
require_positive_integer runs "$runs"
require_external "$result_directory"
[[ "$(uname -s)" == Darwin ]] || {
    echo "resource measurements are calibrated for macOS /usr/bin/time" >&2
    exit 2
}
binary="$(resolve_binary "$binary")"
bcftools="$(resolve_binary "$bcftools")"
[[ "$($bcftools --version | sed -n '1p')" == "bcftools 1.24" ]] || {
    echo "bcftools 1.24 is required" >&2
    exit 2
}
"$bcftools" plugin -l | grep -Fx setGT >/dev/null || {
    echo "the bcftools 1.24 setGT plugin is required" >&2
    exit 2
}
"$binary" --version >/dev/null

formal=false
if ((records >= 2000000 && samples >= 8 && warmups >= 3 && runs >= 10)); then
    formal=true
fi
git_dirty="$(git -C "$repository" status --porcelain)"
if [[ "$formal" == true && -n "$git_dirty" ]]; then
    echo "a formal run requires a clean product worktree" >&2
    exit 2
fi
if [[ -e "$result_directory" ]] && find "$result_directory" -mindepth 1 -print -quit | grep -q .; then
    echo "result directory must be absent or empty: $result_directory" >&2
    exit 2
fi
mkdir -p "$result_directory"
result_directory="$(realpath "$result_directory")"
scratch_root="${TMPDIR:-/Volumes/KIOXIA/Developments/tmp}"
require_external "$scratch_root"
mkdir -p "$scratch_root"
scratch="$(mktemp -d "$scratch_root/rsomics-vcf-setgt.XXXXXX")"
trap 'rm -rf "$scratch"' EXIT

generate_inputs

operations=(all_to_missing missing_to_reference query_to_reference)
formats=(vcf bgzf bcf)
equality="$result_directory/equality.tsv"
printf 'case\trsomics_sha256\tbcftools_sha256\n' > "$equality"
for operation in "${operations[@]}"; do
    for format in "${formats[@]}"; do
        case_name="$operation-$format"
        oracle_hash="$(run_and_hash bcftools "$operation" "$format" "$scratch/oracle.$format")"
        ours_hash="$(run_and_hash rsomics "$operation" "$format" "$scratch/rsomics.$format")"
        [[ "$ours_hash" == "$oracle_hash" ]] || {
            echo "$case_name differs from bcftools 1.24" >&2
            exit 1
        }
        printf '%s\n' "$oracle_hash" > "$scratch/expected-$case_name.sha256"
        printf '%s\t%s\t%s\n' "$case_name" "$ours_hash" "$oracle_hash" >> "$equality"
        rm -f "$scratch/oracle.$format" "$scratch/rsomics.$format"
    done
done

for operation in "${operations[@]}"; do
    for format in "${formats[@]}"; do
        case_name="$operation-$format"
        for ((pair = 1; pair <= warmups; pair++)); do
            if ((pair % 2)); then
                run_and_hash rsomics "$operation" "$format" "$scratch/warm-rsomics.$format" >/dev/null
                run_and_hash bcftools "$operation" "$format" "$scratch/warm-bcftools.$format" >/dev/null
            else
                run_and_hash bcftools "$operation" "$format" "$scratch/warm-bcftools.$format" >/dev/null
                run_and_hash rsomics "$operation" "$format" "$scratch/warm-rsomics.$format" >/dev/null
            fi
        done
        rm -f "$scratch/warm-rsomics.$format" "$scratch/warm-bcftools.$format"
    done
done

raw_results="$result_directory/raw.tsv"
printf 'case\tpair\torder\ttool\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\tsemantic_sha256\n' \
    > "$raw_results"
for operation in "${operations[@]}"; do
    for format in "${formats[@]}"; do
        case_name="$operation-$format"
        for ((pair = 1; pair <= runs; pair++)); do
            if ((pair % 2)); then
                measure "$case_name" "$pair" 1 rsomics "$operation" "$format"
                measure "$case_name" "$pair" 2 bcftools "$operation" "$format"
            else
                measure "$case_name" "$pair" 1 bcftools "$operation" "$format"
                measure "$case_name" "$pair" 2 rsomics "$operation" "$format"
            fi
        done
    done
done

summary="$result_directory/summary.tsv"
printf 'case\ttool\tmetric\tn\tmedian\tp99\tmean\tstdev\n' > "$summary"
for operation in "${operations[@]}"; do
    for format in "${formats[@]}"; do
        case_name="$operation-$format"
        for tool in rsomics bcftools; do
            summarize "$case_name" "$tool" 5 wall_seconds
            summarize "$case_name" "$tool" 6 user_seconds
            summarize "$case_name" "$tool" 7 system_seconds
            summarize "$case_name" "$tool" 8 max_rss_bytes
        done
    done
done

decision="$result_directory/decision.tsv"
printf 'case\trsomics_wall\tbcftools_wall\twall_ratio\trsomics_rss\tbcftools_rss\trss_ratio\tdecision\n' > "$decision"
passing_cases=0
for operation in "${operations[@]}"; do
    for format in "${formats[@]}"; do
        case_name="$operation-$format"
        ours_wall="$(median_value "$case_name" rsomics wall_seconds)"
        oracle_wall="$(median_value "$case_name" bcftools wall_seconds)"
        ours_rss="$(median_value "$case_name" rsomics max_rss_bytes)"
        oracle_rss="$(median_value "$case_name" bcftools max_rss_bytes)"
        case_decision="$(awk -v ow="$ours_wall" -v bw="$oracle_wall" -v or="$ours_rss" -v br="$oracle_rss" \
            'BEGIN { print (ow < bw || or < br) ? "pass" : "fail" }')"
        [[ "$case_decision" == pass ]] && passing_cases=$((passing_cases + 1))
        awk -v case_name="$case_name" -v ow="$ours_wall" -v bw="$oracle_wall" \
            -v or="$ours_rss" -v br="$oracle_rss" -v result="$case_decision" 'BEGIN {
                wall_ratio = bw > 0 ? ow / bw : -1
                rss_ratio = br > 0 ? or / br : -1
                printf "%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%s\n",
                    case_name, ow, bw, wall_ratio, or, br, rss_ratio, result
            }' >> "$decision"
    done
done
if [[ "$formal" == true ]]; then
    overall_decision=fail
    if ((passing_cases > 0)); then
        overall_decision=pass
    fi
else
    overall_decision=smoke
fi
printf 'overall\t-\t-\t-\t-\t-\t-\t%s\n' "$overall_decision" >> "$decision"

rustc="$(resolve_binary "${RSOMICS_RUSTC:-rustc}")"
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
    printf 'git_status_begin\n%s\ngit_status_end\n' "$git_dirty"
    printf 'records=%s\nsamples=%s\nwarmups=%s\nruns=%s\nformal=%s\n' \
        "$records" "$samples" "$warmups" "$runs" "$formal"
    printf 'bcftools_plugins=%s\n' "${BCFTOOLS_PLUGINS:-auto}"
    printf 'rsomics_command=%q setgt OPERATION -O FORMAT -o OUTPUT INPUT\n' "$binary"
    printf 'bcftools_command=%q +setGT INPUT --no-version FORMAT -o OUTPUT -- OPERATION\n' "$bcftools"
    shasum -a 256 "$binary" "$bcftools" "$rustc" "$repository/benchmarks/setgt-vs-bcftools.sh" \
        "$result_directory/generation.txt" "$raw_results" "$summary" "$decision" "$equality"
} > "$result_directory/provenance.txt"

if [[ "$formal" == true && "$overall_decision" != pass ]]; then
    exit 1
fi
