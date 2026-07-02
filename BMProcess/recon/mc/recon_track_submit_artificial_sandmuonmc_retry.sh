#!/bin/bash

set -eu

# ========= path settings =========
VERTEX_DIR="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/2-recon/vertex"

TRACK_DIR="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC_retry/2-recon/track"
OUT_DIR="$TRACK_DIR/out"
LOG_DIR="$TRACK_DIR/log"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/home/nu/notani/wagasci-babymind-track-reconstruction/bin/WagasciRecon"

QUEUE="h"

# ========= prepare directories =========
mkdir -p "$TRACK_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for vertex_file in "$VERTEX_DIR"/artificial_sandmuonmc_vertex_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$vertex_file" ] || continue

    vertex_base=$(basename "$vertex_file")
    # Common part obtained by removing "BMPMWG_vertex_"
    suffix="${vertex_base#artificial_sandmuonmc_vertex_}"
    out_file="$OUT_DIR/artificial_sandmuonmc_track_${suffix%.root}.out"
    log_file="$LOG_DIR/artificial_andmuonmc_track_${suffix%.root}.log"
    track_file="$TRACK_DIR/artificial_sandmuonmc_track_$suffix"

    echo "Submitting job for:"
    echo "  vertex  : $vertex_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  track: $track_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-track-mode \
        --input-file-path "$vertex_file" \
        --output-file-path "$track_file" \
        --log-file-path "$log_file" \
        --data-type mc \
        --copy-ttrees-except-tree
done
