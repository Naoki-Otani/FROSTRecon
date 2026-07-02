#!/usr/bin/env bash

set -euo pipefail

MPPC_CORRECTION="/home/nu/notani/FROSTRecon/ProcessMC/src/mppc_correction"

QUEUE="h"

INPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/3-FROSTResponse"
OUTPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/4-FROSTResponse_mppccorrection"
OUT_DIR="$OUTPUT_DIR/out"

mkdir -p "${OUTPUT_DIR}" "${OUT_DIR}"

shopt -s nullglob

input_files=("${INPUT_DIR}"/artificial_sandmuonmc_frostresponse_*.root)

if [ ${#input_files[@]} -eq 0 ]; then
  echo "No input files found in ${INPUT_DIR}"
  exit 1
fi

for input_file in "${input_files[@]}"; do
  input_base=$(basename "${input_file}")

  # artificial_sandmuonmc_frostresponse_XXXX.root -> artificial_sandmuonmc_frostresponse_aftermppccorrection_XXXX.root
  suffix="${input_base#artificial_sandmuonmc_frostresponse_}"
  output_file="${OUTPUT_DIR}/artificial_sandmuonmc_frostresponse_aftermppccorrection_${suffix}"
  out_file="$OUT_DIR/artificial_sandmuonmc_frostresponse_aftermppccorrection_${suffix%.root}.out"

  echo "Input : ${input_file}"
  echo "Output: ${output_file}"

  bsub -q "$QUEUE" -o "${out_file}" -N \
    "${MPPC_CORRECTION}" \
    -i "${input_file}" \
    -o "${output_file}"
done
