#!/usr/bin/env bash
set -euo pipefail

RECO="./FROST_reconstruction"

IN="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSToutput/aftermppccorrection/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection.root"

OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/rootfile_afterrecon"
MAP_DIR="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mapfunc"

# Alpha values
ALPHAS=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0)

# chi2-threshold values corresponding to each alpha above
# Replace these example values with the ones you want to use.
CHI2_THRESHOLDS=(4.28 2.61 1.76 1.39 1.30 1.28 1.26 1.26 1.27 1.28 1.27 1.28 1.29)
mkdir -p "$OUT_DIR"

# Sanity check: both arrays must have the same length
if [ "${#ALPHAS[@]}" -ne "${#CHI2_THRESHOLDS[@]}" ]; then
  echo "[ERROR] ALPHAS and CHI2_THRESHOLDS must have the same length." >&2
  exit 1
fi

for i in "${!ALPHAS[@]}"; do
  a="${ALPHAS[$i]}"
  thr="${CHI2_THRESHOLDS[$i]}"

  OUT="${OUT_DIR}/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${a}.root"
  SINGLE_MAP="${MAP_DIR}/mapfunc_singlehit_${a}.root"
  TWO_MAP="${MAP_DIR}/mapfunc_twohit_${a}.root"

  echo "[RUN] alpha=${a}, chi2-threshold=${thr}"
  echo "      out=${OUT}"
  echo "      single-map=${SINGLE_MAP}"
  echo "      two-map=${TWO_MAP}"

  "$RECO" \
    --mc \
    --in "$IN" \
    --out "$OUT" \
    --single-map "$SINGLE_MAP" \
    --two-map "$TWO_MAP" \
    --chi2-threshold "$thr"
done

echo "[DONE] All reconstructions completed."
