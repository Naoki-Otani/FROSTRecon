#!/bin/bash

set -eu

for zshift in 0 1 2 3 4 5 -1 -2 -3 -4 -5; do
    # ========= path settings =========
    ROOT_IN_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_HitConverter_added2"
    ROOT_OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST/rootfile_after_TrackMatch_added2/zshift/zshift_${zshift}"
    OUT_DIR="$ROOT_OUT_DIR/out"

    SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
    TrackMatch="/home/nu/notani/FROSTRecon/FROST_BM_Match/install/bin/TrackMatch/TrackMatch"

    QUEUE="s"

    # ========= prepare directories =========
    mkdir -p "$ROOT_OUT_DIR" "$OUT_DIR"

    # ========= loop over all input files =========
    for input_file in "$ROOT_IN_DIR"/*.root; do
        # Safeguard for the case where the glob matches nothing
        [ -e "$input_file" ] || continue

        input_base=$(basename "$input_file")

        root_out_file="$ROOT_OUT_DIR/${input_base%_afterHitConverter.root}_afterTrackMatch_zshift$zshift.root"
        out_file="$OUT_DIR/${input_base%_afterHitConverter.root}_afterTrackMatch_zshift$zshift.out"

        echo "Submitting job for:"
        echo "  input ROOT file : $input_file"
        echo "  OUT : $out_file"
        echo "  output ROOT file : $root_out_file"

        bsub -q "$QUEUE" -o "$out_file" -N \
            apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
            "$TrackMatch" \
            "$input_file" \
            "$root_out_file" \
            $zshift 1
    done
done
