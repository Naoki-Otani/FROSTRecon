#!/bin/bash

set -eu

MACRO="alignment.C"

for zshift in 0 10 20 30 40 50 60 70 80 90 100 -10 -20 -30 -40 -50 -60 -70 -80 -90 -100; do
    INPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/latest20260504/rootfile_after_TrackMatch/zshift/zshift_${zshift}"
    OUTPUT_PDF="/group/nu/ninja/work/otani/FROSTReconData//BM_FROST/analysis_plot/alignment/alignment_zshift_${zshift}.pdf"
    OUTPUT_LOG="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment/alignment_zshift_${zshift}.log"
    OUTPUT_CSV="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment/alignment_zshift_${zshift}_fit.csv"

    excludedFiles1="BMPM_track_2025-11-29_13-46-59_Run0_afterTrackMatch_zshift${zshift}.root"
    excludedFiles2="BMPM_track_2025-11-30_13-11-36_Run0_afterTrackMatch_zshift${zshift}.root"

    OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/analysis_plot/alignment/out"
    QUEUE="s"

    if [ ! -d "$INPUT_DIR" ]; then
        echo "Error: input directory does not exist: $INPUT_DIR"
        exit 1
    fi

    mkdir -p "$OUT_DIR"

    echo "Processing: $INPUT_DIR"
    bsub -q "$QUEUE" -o "$OUT_DIR"/$(basename "$OUTPUT_PDF" .pdf).out -N root -l -b -q "${MACRO}(\"${INPUT_DIR}\",\"${OUTPUT_PDF}\",\"${OUTPUT_LOG}\",\"${OUTPUT_CSV}\", std::vector<std::string>{\"${excludedFiles1}\", \"${excludedFiles2}\"})"

    echo "Done for zshift = $zshift."
done
