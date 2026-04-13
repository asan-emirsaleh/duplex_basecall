#!/usr/bin/env bash
# Compare read lengths and counts between untrimmed and trimmed FASTQ data.
# Handles both standard adapter trimming (same read ID, shorter length) and
# mid-adapter splitting (new read IDs with a name suffix like :1/:2 or _part1/_part2).
#
# Outputs in --output dir:
#   untrimmed_stats.tsv         read_id, length, source_file
#   trimmed_stats.tsv           read_id, length, source_file
#   trimming_comparison.tsv     per-read delta for reads present in both (end-trimming only)
#   split_reads.tsv             new reads in trimmed not matching any untrimmed ID (mid-adapter splits)
#   counts_per_file.tsv         read counts before and after per source file
#   trimming_summary.tsv        aggregate statistics
#
# Requires: bioawk, sort, join, awk

PARSED_ARGUMENTS=$(getopt -o u:t:o: --long untrimmed:,trimmed:,output: -- "$@")
[[ $? -ne 0 ]] && exit 1
eval set -- "$PARSED_ARGUMENTS"

untrimmed=""
trimmed_dir=""
output_dir=""

while true; do
  case "$1" in
    -u|--untrimmed) untrimmed="$2";   shift 2 ;;
    -t|--trimmed)   trimmed_dir="$2"; shift 2 ;;
    -o|--output)    output_dir="$2";  shift 2 ;;
    --) shift; break ;;
    *) echo "Invalid option: $1"; exit 1 ;;
  esac
done

if [[ -z "$untrimmed" || -z "$trimmed_dir" || -z "$output_dir" ]]; then
    echo "Usage: $0 --untrimmed <file_or_dir> --trimmed <dir> --output <dir>"
    exit 1
fi

mkdir -p "${output_dir}"

# Extract read_id, length, source_file from a file or directory.
# Sorted on read_id (field 1).
extract_stats() {
    local src="$1"
    if [[ -f "$src" ]]; then
        bioawk -c fastx -v file="$(basename "$src")" \
            'BEGIN{OFS="\t"} {print $name, length($seq), file}' "$src"
    else
        for fq in "${src}"/*.fastq "${src}"/*.fastq.gz "${src}"/*.fq "${src}"/*.fq.gz; do
            [[ -e "$fq" ]] || continue
            bioawk -c fastx -v file="$(basename "$fq")" \
                'BEGIN{OFS="\t"} {print $name, length($seq), file}' "$fq"
        done
    fi | sort -k1,1
}

untrimmed_stats="${output_dir}/untrimmed_stats.tsv"
trimmed_stats="${output_dir}/trimmed_stats.tsv"
comparison="${output_dir}/trimming_comparison.tsv"
split_reads="${output_dir}/split_reads.tsv"
counts="${output_dir}/counts_per_file.tsv"
summary="${output_dir}/trimming_summary.tsv"

echo "Extracting stats..."
extract_stats "${untrimmed}"   > "${untrimmed_stats}"
extract_stats "${trimmed_dir}" > "${trimmed_stats}"

# --- End-trimming comparison (reads present in both, matched by exact ID) ---
echo -e "read_id\tuntrimmed_len\ttrimmed_len\tdelta" > "${comparison}"
join -t$'\t' -j1 \
    <(awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' "${untrimmed_stats}") \
    <(awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' "${trimmed_stats}") \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' \
    >> "${comparison}"

# --- Mid-adapter split reads (IDs in trimmed that have no exact match in untrimmed) ---
# Strategy: strip known suffixes (:1, :2, _part1, _part2) from trimmed IDs,
# then find those whose base ID exists in untrimmed (confirming they are splits,
# not entirely new reads from some other source).
echo -e "trimmed_id\tbase_id\tlength\tsource_file" > "${split_reads}"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3}' "${trimmed_stats}" \
    | awk -F'\t' 'BEGIN{OFS="\t"} {
        id = $1
        base = id
        # Strip common porechop / guppy split suffixes
        gsub(/:1$|:2$|_part1$|_part2$|;1$|;2$/, "", base)
        if (base != id) print id, base, $2, $3
    }' \
    > "${output_dir}/_split_candidates.tsv"

# Keep only those whose base ID is actually in the untrimmed set
awk -F'\t' '{print $1}' "${untrimmed_stats}" | sort -u > "${output_dir}/_untrimmed_ids.txt"
awk -F'\t' 'NR==FNR {known[$1]=1; next} known[$2] {print}' \
    "${output_dir}/_untrimmed_ids.txt" \
    "${output_dir}/_split_candidates.tsv" \
    >> "${split_reads}"
rm "${output_dir}/_split_candidates.tsv" "${output_dir}/_untrimmed_ids.txt"

# --- Read counts per source file ---
echo -e "file\tuntrimmed_count\ttrimmed_count" > "${counts}"
join -t$'\t' -j1 \
    <(awk -F'\t' '{print $3}' "${untrimmed_stats}" | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1) \
    <(awk -F'\t' '{print $3}' "${trimmed_stats}"   | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1) \
    >> "${counts}"

# --- Aggregate summary ---
n_untrimmed=$(awk 'END{print NR}' "${untrimmed_stats}")
n_trimmed=$(awk 'END{print NR}' "${trimmed_stats}")
n_split=$(awk 'NR>1' "${split_reads}" | wc -l)
n_compared=$(awk 'NR>1' "${comparison}" | wc -l)

awk -F'\t' 'NR>1 {
    n++
    if ($4 < 0) end_trimmed++
    if ($4 == 0) unchanged++
    delta += $4
    if ($4 < 0) bases_removed += -$4
}
END {
    printf "untrimmed_read_count\t%d\n",     '"${n_untrimmed}"'
    printf "trimmed_read_count\t%d\n",       '"${n_trimmed}"'
    printf "reads_compared\t%d\n",           n
    printf "reads_end_trimmed\t%d\n",        end_trimmed
    printf "reads_unchanged\t%d\n",          unchanged
    printf "reads_split_by_mid_adapter\t%d\n", '"${n_split}"'
    printf "mean_length_delta\t%.2f\n",      delta/n
    printf "total_bases_removed\t%d\n",      bases_removed
}' "${comparison}" > "${summary}"

cat "${summary}"