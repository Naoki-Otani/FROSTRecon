#!/bin/bash

set -eu

# ========= path settings =========
SEED_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/3-BMPMWGrecon/seed/with_deadplane_correction"

VERTEX_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/3-BMPMWGrecon/vertex/with_deadplane_correction"
OUT_DIR="$VERTEX_DIR/out"
LOG_DIR="$VERTEX_DIR/log"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/home/nu/notani/wagasci-babymind-track-reconstruction/bin/WagasciRecon"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$VERTEX_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for seed_file in "$SEED_DIR"/BMPMWG_seed_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$seed_file" ] || continue

    seed_base=$(basename "$seed_file")
    # Common part obtained by removing "BMPMWG_seed_"
    suffix="${seed_base#BMPMWG_seed_}"

    out_file="$OUT_DIR/BMPMWG_vertex_${suffix%.root}.out"
    log_file="$LOG_DIR/BMPMWG_vertex_${suffix%.root}.log"
    vertex_file="$VERTEX_DIR/BMPMWG_vertex_$suffix"

    echo "Submitting job for:"
    echo "  BMPMWG_seed  : $seed_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  BMPMWG_vertex: $vertex_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-vertex-mode \
        --input-file-path "$seed_file" \
        --output-file-path "$vertex_file" \
        --log-file-path "$log_file" \
        --data-type data \
        --deadpln_bm_run15
done
