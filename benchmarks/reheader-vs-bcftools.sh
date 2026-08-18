#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "usage:" >&2
    echo "  $0 generate RSOMICS_VCF BCFTOOLS WORKSPACE" >&2
    echo "  $0 run RSOMICS_VCF BCFTOOLS WORKSPACE RESULT_DIRECTORY" >&2
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

require_oracle() {
    [[ "$($1 --version | sed -n '1p')" == "bcftools 1.24" ]] || {
        echo "bcftools 1.24 is required" >&2
        exit 2
    }
}

plain_body_hash() {
    awk '!/^#/' "$1" | shasum -a 256 | awk '{ print $1 }'
}

bgzf_body_hash() {
    gzip -dc "$1" | awk '!/^#/' | shasum -a 256 | awk '{ print $1 }'
}

raw_tail_hash() {
    ruby -rzlib -rstringio -rdigest - "$1" <<'RUBY'
path = ARGV.fetch(0)
canonical_eof = [
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
  0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
].pack("C*")
digest = Digest::SHA256.new
body_seen = false
eof_seen = false

File.open(path, "rb") do |io|
  loop do
    header = io.read(18)
    raise "missing canonical BGZF EOF" if header.nil? || header.empty?
    raise "partial BGZF header" unless header.bytesize == 18
    raise "invalid BGZF header" unless header.byteslice(0, 4) == "\x1f\x8b\x08\x04".b
    size = header.byteslice(16, 2).unpack1("v") + 1
    payload = io.read(size - 18)
    raise "partial BGZF payload" unless payload&.bytesize == size - 18
    raw = header + payload
    digest << raw if body_seen
    if raw == canonical_eof
      raise "bytes follow canonical BGZF EOF" unless io.read(1).nil?
      eof_seen = true
      break
    end
    decoded = Zlib::GzipReader.new(StringIO.new(raw)).read
    body_seen = true if decoded.match?(/(?:\A|\n)[^#\n]/)
  end
end

raise "BGZF contains no records" unless body_seen
raise "missing canonical BGZF EOF" unless eof_seen
puts digest.hexdigest
RUBY
}

generate_workload() {
    [[ $# -eq 3 ]] || usage
    local binary bcftools workspace records samples marker
    binary="$(resolve_binary "$1")"
    bcftools="$(resolve_binary "$2")"
    workspace="$3"
    external_directory "$workspace"
    require_oracle "$bcftools"
    "$binary" --version >/dev/null
    records="${SMOKE_RECORDS:-2000000}"
    samples="${REHEADER_SAMPLES:-8}"
    [[ "$records" =~ ^[0-9]+$ && "$records" -gt 0 ]] || {
        echo "record count must be a positive integer" >&2
        exit 2
    }
    [[ "$samples" =~ ^[0-9]+$ && "$samples" -gt 0 ]] || {
        echo "sample count must be a positive integer" >&2
        exit 2
    }

    mkdir -p "$workspace"
    workspace="$(realpath "$workspace")"
    marker="$workspace/generation.sha256"
    if [[ -f "$marker" ]]; then
        (cd "$workspace" && shasum -a 256 --check generation.sha256)
        return
    fi
    if find "$workspace" -mindepth 1 -print -quit | grep -q .; then
        echo "workspace is nonempty and has no complete generation manifest: $workspace" >&2
        exit 2
    fi
    mkdir -p "$workspace/inputs"

    awk -v records="$records" -v samples="$samples" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##source=reheader-benchmark"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">"
        print "##contig=<ID=chr1,length=" records + 10 ",assembly=benchmark>"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (sample = 1; sample <= samples; sample++) printf "\tS%d", sample
        print ""
        for (position = 1; position <= records; position++) {
            printf "chr1\t%d\t.\tA\tC\t50\tPASS\tDP=%d\tGT:DP", position, samples * 10
            for (sample = 1; sample <= samples; sample++) printf "\t0/1:10"
            print ""
        }
    }' > "$workspace/inputs/input.vcf"

    awk -v records="$records" -v samples="$samples" 'BEGIN {
        print "##fileformat=VCFv4.3"
        print "##source=reheader-benchmark-replacement"
        print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
        print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">"
        print "##contig=<ID=chr1,length=" records + 20 ",assembly=replacement>"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for (sample = 1; sample <= samples; sample++) printf "\tS%d", sample
        print ""
    }' > "$workspace/inputs/replacement.vcfh"

    printf 'chr1\t%s\t0\t0\t0\nchr2\t1000\t0\t0\t0\n' "$((records + 30))" \
        > "$workspace/inputs/reference.fai"
    awk -v samples="$samples" 'BEGIN {
        for (sample = 1; sample <= samples; sample++)
            printf "S%d\tN%d\n", sample, sample
    }' > "$workspace/inputs/samples.tsv"

    "$bcftools" view --no-version -Oz -o "$workspace/inputs/input.vcf.gz" \
        "$workspace/inputs/input.vcf"
    [[ "$(plain_body_hash "$workspace/inputs/input.vcf")" == \
        "$(bgzf_body_hash "$workspace/inputs/input.vcf.gz")" ]]
    {
        printf 'records=%s\nsamples=%s\n' "$records" "$samples"
        wc -c "$workspace"/inputs/*
    } > "$workspace/generation.txt"
    (
        cd "$workspace"
        shasum -a 256 inputs/* generation.txt > generation.sha256
    )
}

measure() {
    local path="$1"
    local pair="$2"
    local order="$3"
    local tool="$4"
    shift 4
    local timing="$result_directory/$path-$tool-$pair.time"
    /usr/bin/time -lp "$@" > /dev/null 2> "$timing"
    local wall user system rss
    wall="$(awk '$1 == "real" { print $2 }' "$timing")"
    user="$(awk '$1 == "user" { print $2 }' "$timing")"
    system="$(awk '$1 == "sys" { print $2 }' "$timing")"
    rss="$(awk '$2 == "maximum" && $3 == "resident" { print $1 }' "$timing")"
    [[ -n "$wall" && -n "$user" && -n "$system" && -n "$rss" ]]
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$path" "$pair" "$order" "$tool" "$wall" "$user" "$system" "$rss" \
        >> "$raw_results"
}

summarize() {
    local path="$1"
    local tool="$2"
    local column="$3"
    local metric="$4"
    local values="$result_directory/$path-$tool-$metric.values"
    awk -F '\t' -v path="$path" -v tool="$tool" -v column="$column" \
        '$1 == path && $4 == tool { print $column }' "$raw_results" | LC_ALL=C sort -n > "$values"
    awk -v path="$path" -v tool="$tool" -v metric="$metric" '
        { value[NR] = $1; sum += $1; sumsq += $1 * $1 }
        END {
            if (NR == 0) exit 1
            median = NR % 2 ? value[(NR + 1) / 2] : (value[NR / 2] + value[NR / 2 + 1]) / 2
            p99 = value[int(0.99 * NR + 0.999999)]
            mean = sum / NR
            stdev = NR > 1 ? sqrt((sumsq - sum * sum / NR) / (NR - 1)) : 0
            printf "%s\t%s\t%s\t%d\t%.9g\t%.9g\t%.9g\t%.9g\n", path, tool, metric, NR, median, p99, mean, stdev
        }
    ' "$values" >> "$summary"
}

median_value() {
    awk -F '\t' -v path="$1" -v tool="$2" -v metric="$3" \
        '$1 == path && $2 == tool && $3 == metric { print $5 }' "$summary"
}

run_gate() {
    [[ $# -eq 4 ]] || usage
    local binary bcftools workspace repository rustc runs warmups minimum_runs minimum_warmups git_dirty
    binary="$(resolve_binary "$1")"
    bcftools="$(resolve_binary "$2")"
    workspace="$(realpath "$3")"
    result_directory="$4"
    external_directory "$workspace"
    external_directory "$result_directory"
    [[ "$(uname -s)" == "Darwin" ]] || {
        echo "resource measurements are calibrated for macOS /usr/bin/time" >&2
        exit 2
    }
    require_oracle "$bcftools"
    repository="$(cd "$(dirname "$0")/.." && pwd)"
    [[ -z "$(git -C "$repository" status --porcelain)" || "${RSOMICS_BENCH_SMOKE:-0}" == 1 ]] || {
        echo "the product worktree must be clean" >&2
        exit 2
    }
    (cd "$workspace" && shasum -a 256 --check generation.sha256)
    if [[ -e "$result_directory" ]] && find "$result_directory" -mindepth 1 -print -quit | grep -q .; then
        echo "result directory must be absent or empty: $result_directory" >&2
        exit 2
    fi
    mkdir -p "$result_directory"
    result_directory="$(realpath "$result_directory")"

    runs="${BENCH_RUNS:-10}"
    warmups="${BENCH_WARMUPS:-3}"
    minimum_runs=10
    minimum_warmups=3
    if [[ "${RSOMICS_BENCH_SMOKE:-0}" == 1 ]]; then
        minimum_runs=2
        minimum_warmups=1
    fi
    [[ "$runs" =~ ^[0-9]+$ && "$runs" -ge "$minimum_runs" ]] || {
        echo "BENCH_RUNS must be an integer of at least $minimum_runs" >&2
        exit 2
    }
    [[ "$warmups" =~ ^[0-9]+$ && "$warmups" -ge "$minimum_warmups" ]] || {
        echo "BENCH_WARMUPS must be an integer of at least $minimum_warmups" >&2
        exit 2
    }

    local replacement="$workspace/inputs/replacement.vcfh"
    local fai="$workspace/inputs/reference.fai"
    local samples="$workspace/inputs/samples.tsv"
    local plain="$workspace/inputs/input.vcf"
    local bgzf="$workspace/inputs/input.vcf.gz"
    local ours_plain="$result_directory/rsomics.plain.vcf"
    local oracle_plain="$result_directory/bcftools.plain.vcf"
    local ours_bgzf="$result_directory/rsomics.vcf.gz"
    local oracle_bgzf="$result_directory/bcftools.vcf.gz"

    "$binary" reheader -H "$replacement" -f "$fai" -N "$samples" -o "$ours_plain" "$plain"
    "$bcftools" reheader -h "$replacement" -f "$fai" -N "$samples" -o "$oracle_plain" "$plain"
    "$binary" reheader -H "$replacement" -f "$fai" -N "$samples" -o "$ours_bgzf" "$bgzf"
    "$bcftools" reheader -h "$replacement" -f "$fai" -N "$samples" -o "$oracle_bgzf" "$bgzf"
    cmp "$ours_plain" "$oracle_plain"
    cmp <(gzip -dc "$ours_bgzf") <(gzip -dc "$oracle_bgzf")

    local input_plain_body ours_plain_body oracle_plain_body
    local input_bgzf_body ours_bgzf_body oracle_bgzf_body input_tail ours_tail
    input_plain_body="$(plain_body_hash "$plain")"
    ours_plain_body="$(plain_body_hash "$ours_plain")"
    oracle_plain_body="$(plain_body_hash "$oracle_plain")"
    input_bgzf_body="$(bgzf_body_hash "$bgzf")"
    ours_bgzf_body="$(bgzf_body_hash "$ours_bgzf")"
    oracle_bgzf_body="$(bgzf_body_hash "$oracle_bgzf")"
    input_tail="$(raw_tail_hash "$bgzf")"
    ours_tail="$(raw_tail_hash "$ours_bgzf")"
    [[ "$input_plain_body" == "$ours_plain_body" && "$input_plain_body" == "$oracle_plain_body" ]]
    [[ "$input_bgzf_body" == "$ours_bgzf_body" && "$input_bgzf_body" == "$oracle_bgzf_body" ]]
    [[ "$input_tail" == "$ours_tail" ]]
    {
        printf 'check\tinput\trsomics\tbcftools\n'
        printf 'plain_body_sha256\t%s\t%s\t%s\n' "$input_plain_body" "$ours_plain_body" "$oracle_plain_body"
        printf 'bgzf_body_sha256\t%s\t%s\t%s\n' "$input_bgzf_body" "$ours_bgzf_body" "$oracle_bgzf_body"
        printf 'bgzf_raw_tail_sha256\t%s\t%s\tnot-required\n' "$input_tail" "$ours_tail"
    } > "$result_directory/equality.tsv"

    local timed_ours_plain="$result_directory/timed-rsomics.plain.vcf"
    local timed_oracle_plain="$result_directory/timed-bcftools.plain.vcf"
    local timed_ours_bgzf="$result_directory/timed-rsomics.vcf.gz"
    local timed_oracle_bgzf="$result_directory/timed-bcftools.vcf.gz"
    ours_plain_command=("$binary" reheader -H "$replacement" -f "$fai" -N "$samples" -o "$timed_ours_plain" "$plain")
    oracle_plain_command=("$bcftools" reheader -h "$replacement" -f "$fai" -N "$samples" -o "$timed_oracle_plain" "$plain")
    ours_bgzf_command=("$binary" reheader -H "$replacement" -f "$fai" -N "$samples" -o "$timed_ours_bgzf" "$bgzf")
    oracle_bgzf_command=("$bcftools" reheader -h "$replacement" -f "$fai" -N "$samples" -o "$timed_oracle_bgzf" "$bgzf")

    local path pair
    local -a ours_command oracle_command
    for path in plain bgzf; do
        if [[ "$path" == plain ]]; then
            ours_command=("${ours_plain_command[@]}")
            oracle_command=("${oracle_plain_command[@]}")
        else
            ours_command=("${ours_bgzf_command[@]}")
            oracle_command=("${oracle_bgzf_command[@]}")
        fi
        for ((pair = 1; pair <= warmups; pair++)); do
            if ((pair % 2)); then
                "${ours_command[@]}"
                "${oracle_command[@]}"
            else
                "${oracle_command[@]}"
                "${ours_command[@]}"
            fi
        done
    done

    raw_results="$result_directory/raw.tsv"
    printf 'path\tpair\torder\ttool\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_bytes\n' \
        > "$raw_results"
    for path in plain bgzf; do
        if [[ "$path" == plain ]]; then
            ours_command=("${ours_plain_command[@]}")
            oracle_command=("${oracle_plain_command[@]}")
        else
            ours_command=("${ours_bgzf_command[@]}")
            oracle_command=("${oracle_bgzf_command[@]}")
        fi
        for ((pair = 1; pair <= runs; pair++)); do
            if ((pair % 2)); then
                measure "$path" "$pair" 1 rsomics "${ours_command[@]}"
                measure "$path" "$pair" 2 bcftools "${oracle_command[@]}"
            else
                measure "$path" "$pair" 1 bcftools "${oracle_command[@]}"
                measure "$path" "$pair" 2 rsomics "${ours_command[@]}"
            fi
        done
    done

    summary="$result_directory/summary.tsv"
    printf 'path\ttool\tmetric\tn\tmedian\tp99\tmean\tstdev\n' > "$summary"
    local tool
    for path in plain bgzf; do
        for tool in rsomics bcftools; do
            summarize "$path" "$tool" 5 wall_seconds
            summarize "$path" "$tool" 6 user_seconds
            summarize "$path" "$tool" 7 system_seconds
            summarize "$path" "$tool" 8 max_rss_bytes
        done
    done

    local ours_wall oracle_wall ours_rss oracle_rss result
    printf 'path\trsomics_wall\tbcftools_wall\twall_ratio\trsomics_rss\tbcftools_rss\trss_ratio\tdecision\n' \
        > "$result_directory/decision.tsv"
    for path in plain bgzf; do
        ours_wall="$(median_value "$path" rsomics wall_seconds)"
        oracle_wall="$(median_value "$path" bcftools wall_seconds)"
        ours_rss="$(median_value "$path" rsomics max_rss_bytes)"
        oracle_rss="$(median_value "$path" bcftools max_rss_bytes)"
        if [[ "${RSOMICS_BENCH_SMOKE:-0}" == 1 ]]; then
            result=smoke
        else
            result="$(awk -v ow="$ours_wall" -v bw="$oracle_wall" -v or="$ours_rss" -v br="$oracle_rss" \
                'BEGIN { print (ow < bw || or < br) ? "pass" : "fail" }')"
        fi
        awk -v path="$path" -v ow="$ours_wall" -v bw="$oracle_wall" -v or="$ours_rss" -v br="$oracle_rss" -v result="$result" \
            'BEGIN {
                wall_ratio = bw > 0 ? ow / bw : 0
                rss_ratio = br > 0 ? or / br : 0
                printf "%s\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%s\n", path, ow, bw, wall_ratio, or, br, rss_ratio, result
            }' \
            >> "$result_directory/decision.tsv"
    done

    rustc="$(resolve_binary "${RSOMICS_RUSTC:-rustc}")"
    if [[ -n "$(git -C "$repository" status --porcelain)" ]]; then
        git_dirty=true
    else
        git_dirty=false
    fi
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
        printf 'git_head=%s\ngit_dirty=%s\nruns=%s\nwarmups=%s\n' \
            "$(git -C "$repository" rev-parse HEAD)" "$git_dirty" "$runs" "$warmups"
        printf 'rsomics_plain_command='; printf '%q ' "${ours_plain_command[@]}"; printf '\n'
        printf 'bcftools_plain_command='; printf '%q ' "${oracle_plain_command[@]}"; printf '\n'
        printf 'rsomics_bgzf_command='; printf '%q ' "${ours_bgzf_command[@]}"; printf '\n'
        printf 'bcftools_bgzf_command='; printf '%q ' "${oracle_bgzf_command[@]}"; printf '\n'
        wc -c "$workspace"/inputs/* "$ours_plain" "$oracle_plain" "$ours_bgzf" "$oracle_bgzf"
        shasum -a 256 "$binary" "$bcftools" "$rustc" "$repository/benchmarks/reheader-vs-bcftools.sh" \
            "$workspace"/inputs/* "$ours_plain" "$oracle_plain" "$ours_bgzf" "$oracle_bgzf" \
            "$raw_results" "$summary" "$result_directory/decision.tsv" "$result_directory/equality.tsv"
    } > "$result_directory/provenance.txt"

    if [[ "${RSOMICS_BENCH_SMOKE:-0}" != 1 ]]; then
        awk -F '\t' 'NR > 1 && $8 != "pass" { failed = 1 } END { exit failed }' \
            "$result_directory/decision.tsv"
    fi
}

case "${1:-}" in
    generate)
        shift
        generate_workload "$@"
        ;;
    run)
        shift
        run_gate "$@"
        ;;
    *) usage ;;
esac
