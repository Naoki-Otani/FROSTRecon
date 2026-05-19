#!/bin/bash

set -eu

# ========= path settings =========
VERTEX_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-BMPMrecon/vertex"

TRACK_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-BMPMrecon/track"
OUT_DIR="$TRACK_DIR/out"
LOG_DIR="$TRACK_DIR/log"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/opt/wagasci_mc/WagasciTR/bin/WagasciRecon"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$TRACK_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for vertex_file in "$VERTEX_DIR"/BMPM_vertex_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$vertex_file" ] || continue

    vertex_base=$(basename "$vertex_file")
    # Common part obtained by removing "BMPM_vertex_"
    suffix="${vertex_base#BMPM_vertex_}"
    out_file="$OUT_DIR/BMPM_vertex_${suffix%.root}.out"
    log_file="$LOG_DIR/BMPM_vertex_${suffix%.root}.log"
    track_file="$TRACK_DIR/BMPM_track_$suffix"

    echo "Submitting job for:"
    echo "  BMPM_vertex  : $vertex_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  BMPM_track: $track_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-track-mode \
        --input-file-path "$vertex_file" \
        --output-file-path "$track_file" \
        --log-file-path "$log_file" \
        --data-type data
done
