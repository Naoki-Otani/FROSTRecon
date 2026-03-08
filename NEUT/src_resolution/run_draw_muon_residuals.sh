#!/usr/bin/env bash
set -euo pipefail

MACRO="draw_muon_residuals.C"

BASE_IN="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/rootfile_afterrecon"
BASE_OUT="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/resolution"

ALPHAS=(1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0)

mkdir -p "$BASE_OUT"

for v in "${ALPHAS[@]}"; do
  input="${BASE_IN}/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${v}.root"
  output_pdf="${BASE_OUT}/resolution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${v}.pdf"
  output_root="${BASE_OUT}/resolution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${v}.root"

  echo "[RUN] value=${v}"
  echo "      input=${input}"
  echo "      output_pdf=${output_pdf}"
  echo "      output_root=${output_root}"

  root -l -q "${MACRO}(\"${input}\",\"${output_pdf}\",\"${output_root}\")"
done

echo "[DONE] All jobs finished."
