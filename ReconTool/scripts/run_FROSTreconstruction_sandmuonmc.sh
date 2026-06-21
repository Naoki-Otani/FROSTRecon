#!/usr/bin/env bash

set -euo pipefail

RECON="/home/nu/notani/FROSTRecon/ReconTool/FROST_reconstruction"

QUEUE="s"

INPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/4-FROSTResponse_mppccorrection"
OUTPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/5-FROSTRecon"
OUT_DIR="$OUTPUT_DIR/out"

mkdir -p "${OUTPUT_DIR}" "${OUT_DIR}"

shopt -s nullglob

input_files=("${INPUT_DIR}"/sandmuonmc_frostresponse_aftermppccorrection_*.root)

if [ ${#input_files[@]} -eq 0 ]; then
  echo "No input files found in ${INPUT_DIR}"
  exit 1
fi

for input_file in "${input_files[@]}"; do
  input_base=$(basename "${input_file}")

  # sandmuonmc_frostresponse_XXXX.root -> sandmuonmc_frostresponse_aftermppccorrection_XXXX.root
  suffix="${input_base#sandmuonmc_frostresponse_aftermppccorrection_}"
  output_file="${OUTPUT_DIR}/sandmuonmc_frostrecon_${suffix}"
  out_file="$OUT_DIR/sandmuonmc_frostrecon_${suffix%.root}.out"

  echo "Input : ${input_file}"
  echo "Output: ${output_file}"

  bsub -q "$QUEUE" -o "${out_file}" -N \
    "${RECON}" \
    --mc \
    --in "$input_file" \
    --out "$output_file"
done
