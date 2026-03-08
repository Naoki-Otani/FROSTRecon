#!/usr/bin/env bash

# Run draw_chi2hist for multiple threshold values.
# This script calls the ROOT macro in batch mode for:
#   1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5

set -euo pipefail

MACRO="draw_chi2hist.C"

BASE_IN="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/rootfile_afterrecon"
BASE_OUT="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/FROSTrecon_result/chi2"

for alpha in 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0; do
    input="${BASE_IN}/frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${alpha}.root"
    output_pdf="${BASE_OUT}/chi2_distribution_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${alpha}.pdf"
    output_txt="${BASE_OUT}/chi2_threshold_results_frost_neut_320kA_1.077e21pot_ECC12_aftermppccorrection_afterrecon_${alpha}.txt"

    echo "Running alpha = ${alpha}"
    root -l -b -q "${MACRO}(\"${input}\",\"${output_pdf}\",\"${output_txt}\")"
done

echo "All jobs finished."
