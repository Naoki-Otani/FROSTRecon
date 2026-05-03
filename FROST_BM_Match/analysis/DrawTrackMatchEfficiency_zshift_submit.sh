#!/bin/bash

set -eu

MACRO="DrawTrackMatchEfficiency.C"

for zshift in 0 1 2 3 4 5 -1 -2 -3 -4 -5; do
    INPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_new/zshift/zshift_${zshift}"
    OUTPUT_PDF="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/zshift/efficiency_zshift_${zshift}.pdf"
    OUTPUT_LOG="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/zshift/efficiency_zshift_${zshift}.log"

    excludedFiles1="BMPM_track_2025-11-29_13-46-59_Run0_afterTrackMatch_zshift${zshift}.root"
    excludedFiles2="BMPM_track_2025-11-30_13-11-36_Run0_afterTrackMatch_zshift${zshift}.root"

    OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/zshift/out"
    QUEUE="s"

    if [ ! -d "$INPUT_DIR" ]; then
        echo "Error: input directory does not exist: $INPUT_DIR"
        exit 1
    fi

    mkdir -p "$OUT_DIR"

    echo "Processing: $INPUT_DIR"
    bsub -q "$QUEUE" -o "$OUT_DIR"/$(basename "$OUTPUT_PDF" .pdf).out -N root -l -b -q "${MACRO}(\"${INPUT_DIR}\",\"${OUTPUT_PDF}\",\"${OUTPUT_LOG}\", std::vector<std::string>{\"${excludedFiles1}\", \"${excludedFiles2}\"})"

    echo "Done for zshift = $zshift."
done
