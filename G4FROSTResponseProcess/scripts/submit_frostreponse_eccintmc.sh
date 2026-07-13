#!/usr/bin/env bash

set -euo pipefail

FROST_RESPONSE="/home/nu/notani/geant4/FROSTResponse/G4FROSTResponse/bin/frost_response"
MAC_TEMPLATE="frost_response_sandmuonmc.mac"

QUEUE="l"

INPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/ECCintMC/2-recon/track"
OUTPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/ECCintMC/3-FROSTResponse"
OUT_DIR="$OUTPUT_DIR/out"
TMP_DIR="$OUTPUT_DIR/tmp"

mkdir -p "${OUTPUT_DIR}" "${OUT_DIR}" "${TMP_DIR}"

shopt -s nullglob

input_files=("${INPUT_DIR}"/eccintmc_track_*.root)

if [ ${#input_files[@]} -eq 0 ]; then
  echo "No input files found in ${INPUT_DIR}"
  exit 1
fi

for input_file in "${input_files[@]}"; do
  input_base=$(basename "${input_file}")

  # eccintmc_track_XXXX.root -> eccintmc_frostresponse_XXXX.root
  suffix="${input_base#eccintmc_track_}"
  output_file="${OUTPUT_DIR}/eccintmc_frostresponse_${suffix}"
  out_file="$OUT_DIR/eccintmc_frostresponse_${suffix%.root}.out"

  tmp_mac="$TMP_DIR/eccintmc_frostresponse_${suffix%.root}.mac"

  # Replace only the /analysis/setFileName line.
  awk -v outfile="${output_file}" '
    /^\/analysis\/setFileName[[:space:]]+/ {
      print "/analysis/setFileName " outfile
      next
    }
    { print }
  ' "${MAC_TEMPLATE}" > "${tmp_mac}"

  echo "Input : ${input_file}"
  echo "Output: ${output_file}"

  bsub -q "$QUEUE" -o "${out_file}" -N \
    "${FROST_RESPONSE}" \
    -mac "${tmp_mac}" \
    -wagascimc "${input_file}"

  # rm -f "${tmp_mac}"
done
