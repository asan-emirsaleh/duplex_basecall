#!/usr/bin/env bash

PARSED_ARGUMENTS=$(getopt -o i:o:s:t: --long input:,output:,search:,threads: -- "$@")
[[ $? -ne 0 ]] && exit 1
eval set -- "$PARSED_ARGUMENTS"

input_dir=""
output_dir=""
adapters=""
threads=8

while true; do
  case "$1" in
    -i|--input)   input_dir="$2";  shift 2 ;;
    -o|--output)  output_dir="$2"; shift 2 ;;
    -s|--search)  adapters="$2";   shift 2 ;;
    -t|--threads) threads="$2";    shift 2 ;;
    --) shift; break ;;
    *)  echo "Invalid option: $1"; exit 1 ;;
  esac
done

if [[ -z "$input_dir" || -z "$output_dir" ]]; then
    echo "Usage: $0 --input <dir> --output <dir> [--search <adapters_out>] [--threads <N>]"
    exit 1
fi

mkdir -p "${output_dir}"

# A single BAM file is expected per directory
for bam in "${input_dir}"/*.bam; do
    base=$(basename "$bam" .bam)

    # Detect adapters
    porechop --input "$bam" --threads "${threads}" > "${adapters}"

    # Trim adapters
    porechop \
        --input "$bam" \
        --output "${output_dir}/${base}.trimmed.bam" \
        --threads "${threads}"
done # bam