#!/usr/bin/env bash
set -euo pipefail

# Path to the executable (adjust if needed)
MAPTOOL="./mapfunction_tool"

# Input patterns (do not change {label}; it is used by mapfunction_tool)
IN_X_PATTERN="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/aftermppccorrection/x/twohitmuonx{label}mmy-5mmto5mm_aftermppccorrection.root"
IN_Y_PATTERN="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/aftermppccorrection/y/twohitmuony{label}mmx-5mmto5mm_aftermppccorrection.root"

# Output directory
OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mapfunc"

# Alphas to run
ALPHAS=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0)

mkdir -p "$OUT_DIR"

for alpha in "${ALPHAS[@]}"; do
  out="${OUT_DIR}/mapfunc_twohit_${alpha}.root"
  echo "[RUN] alpha=${alpha}"
  "$MAPTOOL" \
    --mode two \
    --alpha "${alpha}" \
    --in-x-pattern "${IN_X_PATTERN}" \
    --in-y-pattern "${IN_Y_PATTERN}" \
    --out "${out}"
done

echo "[DONE] All runs completed."
