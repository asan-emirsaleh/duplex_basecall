#!/usr/bin/env bash
# Wrapper script for Dorado basecalling

dorado="_apps/dorado/bin/dorado"
models="./_apps/dorado/models"

# Use getopt for long options
PARSED_ARGUMENTS=$(getopt -o i:o:p: --long input:,output:,pore: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1
fi

eval set -- "$PARSED_ARGUMENTS"

# Default values for optional arguments
input_file=""
output_dir=""
pore="9.4.1"

# Parse options
while true; do
  case "$1" in
    -i|--input)
      input_file="$2"
      shift 2
      ;;
    -o|--output)
      output_dir="$2"
      shift 2
      ;;
    -p|--pore)
      pore="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "Invalid option: $1"
      exit 1
      ;;
  esac
done

# Display values
echo "Input file: $input_file"
echo "Output directory: $output_dir"
echo "Pore: $pore"

# Check if required arguments are provided
if [[ -z "$input_file" || -z "$output_dir" ]]; then
    echo "Usage: $0 --input <input_file> --output <output_dir> [--pore <pore>]"
    exit 1
fi

case ${pore} in
"9.4.1" )
  model="dna_r9.4.1_e8_sup@v3.6"
  ;;
"9.4.0" )
  model="dna_r9.4.1_e8_sup@v3.6"
  ;;
* )
  model="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
esac

mkdir -p "${models}"

[[ -d ${models}/${model} ]] || $dorado download --model ${model} --models-directory ${models}

# Run dorado basecalling
# note that for R9 pore data, trimming is set off by default,
# so `--trim` parameter should be always set to 'all'

gpu_option="cuda:all"

"${dorado}" basecaller \
    --models-directory "${models}" \
    --emit-moves \
    --output-dir "${output_dir}" \
    --trim all \
    -x "${gpu_option}" \
    "${models}/${model}" \
    "${input_file}"

