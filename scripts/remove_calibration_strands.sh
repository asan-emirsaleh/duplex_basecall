#!/usr/bin/env bash
# Align reads to a calibration reference (Lambda phage / PhiX).
# Produces per-type alignment statistics, separate calibration FASTQ files,
# and optionally excludes calibration reads from downstream data.
#
# Inputs (FASTQ directories):
#   --simplex       trimmed simplex reads
#   --duplex-distant trimmed duplex reads from distant (different-molecule) pairs
#   --duplex-split  trimmed duplex reads from split (self-ligated) pairs
#
# Outputs:
#   --calib-out     directory for calibration read FASTQs and alignment stats
#   --simplex-clean / --duplex-distant-clean / --duplex-split-clean
#                   clean read directories (identical to input if --keep is set)
#   --report        TSV: read_id, length, type, mapq
#
# Flags:
#   --keep          retain calibration reads in the clean output (do not exclude them).
#                   Use this when the reference organism shares sequences with the
#                   calibration spike-in (e.g. rRNA comparison), so mapped reads
#                   legitimately belong to the sample.
#
# Requires: minimap2, samtools, bioawk

PARSED_ARGUMENTS=$(getopt -o "" \
    --long simplex:,duplex-distant:,duplex-split:,\
calib-out:,simplex-clean:,duplex-distant-clean:,duplex-split-clean:,\
report:,ref:,threads:,keep \
    -- "$@")
[[ $? -ne 0 ]] && exit 1
eval set -- "$PARSED_ARGUMENTS"

simplex_dir=""
duplex_distant_dir=""
duplex_split_dir=""
calib_out=""
simplex_clean=""
duplex_distant_clean=""
duplex_split_clean=""
report=""
ref=""
threads=8
keep=0

while true; do
  case "$1" in
    --simplex)               simplex_dir="$2";          shift 2 ;;
    --duplex-distant)        duplex_distant_dir="$2";   shift 2 ;;
    --duplex-split)          duplex_split_dir="$2";     shift 2 ;;
    --calib-out)             calib_out="$2";             shift 2 ;;
    --simplex-clean)         simplex_clean="$2";         shift 2 ;;
    --duplex-distant-clean)  duplex_distant_clean="$2";  shift 2 ;;
    --duplex-split-clean)    duplex_split_clean="$2";    shift 2 ;;
    --report)                report="$2";                shift 2 ;;
    --ref)                   ref="$2";                   shift 2 ;;
    --threads)               threads="$2";               shift 2 ;;
    --keep)                  keep=1;                     shift   ;;
    --) shift; break ;;
    *) echo "Invalid option: $1"; exit 1 ;;
  esac
done

for var in simplex_dir duplex_distant_dir duplex_split_dir \
           calib_out simplex_clean duplex_distant_clean duplex_split_clean \
           report ref; do
    [[ -z "${!var}" ]] && { echo "Missing required argument: --${var//_/-}"; exit 1; }
done

mkdir -p "${calib_out}" "${simplex_clean}" "${duplex_distant_clean}" "${duplex_split_clean}"

echo -e "read_id\tlength\tmapq\ttype" > "${report}"

# Align one FASTQ to the calibration reference.
# Args: fq, clean_dir, read_type, stats_prefix
process_fastq() {
    local fq="$1"
    local clean_dir="$2"
    local read_type="$3"
    local base
    base=$(basename "${fq}" .fastq)

    local bam="${calib_out}/${base}.${read_type}.bam"
    local calib_ids="${calib_out}/${base}.${read_type}.calib_ids.txt"
    local calib_fq="${calib_out}/${base}.${read_type}.calib.fastq"
    local stats_file="${calib_out}/${base}.${read_type}.flagstat.txt"
    local aln_stats="${calib_out}/${base}.${read_type}.samtools_stats.txt"

    # Align
    minimap2 -ax map-ont -t "${threads}" "${ref}" "${fq}" 2>/dev/null \
        | samtools sort -@ "${threads}" -o "${bam}"
    samtools index "${bam}"

    # Alignment statistics
    samtools flagstat "${bam}" > "${stats_file}"
    samtools stats "${bam}"   > "${aln_stats}"

    # Extract calibration read IDs (mapped reads)
    samtools view -F 4 "${bam}" \
        | awk '{print $1}' | sort -u > "${calib_ids}"

    # Write calibration reads to their own FASTQ
    bioawk -c fastx -v ids="${calib_ids}" '
        BEGIN { while ((getline id < ids) > 0) keep[id] = 1 }
        keep[$name] { printf "@%s\n%s\n+\n%s\n", $name, $seq, $qual }
    ' "${fq}" > "${calib_fq}"

    # Append to report: read_id, length, mapq, type
    samtools view -F 4 "${bam}" \
        | awk -v type="${read_type}" 'BEGIN{OFS="\t"} {print $1, length($10), $5, type}' \
        >> "${report}"

    # Write clean output
    if [[ "${keep}" -eq 1 ]]; then
        # Keep mode: copy all reads as-is
        cp "${fq}" "${clean_dir}/$(basename "${fq}")"
    else
        # Exclude calibration reads
        bioawk -c fastx -v ids="${calib_ids}" '
            BEGIN { while ((getline id < ids) > 0) skip[id] = 1 }
            !skip[$name] { printf "@%s\n%s\n+\n%s\n", $name, $seq, $qual }
        ' "${fq}" > "${clean_dir}/$(basename "${fq}")"
    fi
}

process_dir() {
    local src_dir="$1"
    local clean_dir="$2"
    local read_type="$3"

    local found=0
    for fq in "${src_dir}"/*.fastq "${src_dir}"/*.fq; do
        [[ -e "$fq" ]] || continue
        found=1
        process_fastq "${fq}" "${clean_dir}" "${read_type}"
    done
    [[ "${found}" -eq 0 ]] && echo "Warning: no FASTQ files found in ${src_dir}" >&2
}

process_dir "${simplex_dir}"        "${simplex_clean}"         "simplex"
process_dir "${duplex_distant_dir}" "${duplex_distant_clean}"  "duplex-distant"
process_dir "${duplex_split_dir}"   "${duplex_split_clean}"    "duplex-split"

n_calib=$(awk 'NR>1' "${report}" | wc -l)
mode=$([[ "${keep}" -eq 1 ]] && echo "kept in dataset (--keep)" || echo "excluded from dataset")
echo "Calibration reads found: ${n_calib} — ${mode}" >&2
echo "Alignment stats written to: ${calib_out}/" >&2
echo "Report written to: ${report}" >&2