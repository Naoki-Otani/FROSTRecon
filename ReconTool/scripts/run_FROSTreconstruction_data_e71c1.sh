#!/usr/bin/env bash
set -euo pipefail

SINGLE_MAP="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mapfunc/mapfunc_singlehit_4.5.root"
TWO_MAP="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mapfunc/mapfunc_twohit_4.5.root"
CHI2_THRESHOLD="1.23"

shopt -s nullglob

for dataset in e71c; do
    INPUT_DIR="/group/nu/ninja/work/otani/FROST_beamdata/${dataset}/rootfile_aftercalib"
    OUTPUT_DIR="/group/nu/ninja/work/otani/FROST_beamdata/${dataset}/rootfile_afterrecon"
  for input_file in "${INPUT_DIR}"/*.root; do
      base_name="$(basename "${input_file}")"

      # 例:
      # run00017_120000_129999_lightyield.root
      # -> run00017_120000_129999_recon.root
      output_name="${base_name%_lightyield.root}_recon.root"
      output_file="${OUTPUT_DIR}/${output_name}"

      echo "Processing: ${input_file}"
      echo "Output    : ${output_file}"

      ./FROST_reconstruction \
          --data \
          --in "${input_file}" \
          --out "${output_file}" \
          --single-map "${SINGLE_MAP}" \
          --two-map "${TWO_MAP}" \
          --chi2-threshold "${CHI2_THRESHOLD}"
  done
done
