#!/bin/bash

set -eu

# ========= path settings =========
CONVERTED_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/2-BMPMWGconverter_divided/with_deadplane_correction"

SEED_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/3-BMPMWGrecon/seed/with_deadplane_correction"
OUT_DIR="$SEED_DIR/out"
LOG_DIR="$SEED_DIR/log"


SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/home/nu/notani/wagasci-babymind-track-reconstruction/bin/WagasciRecon"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$SEED_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for converted_file in "$CONVERTED_DIR"/BMPMWG_converted_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$converted_file" ] || continue

    converted_base=$(basename "$converted_file")
    # Common part obtained by removing "BMPMWG_converted_"
    suffix="${converted_base#BMPMWG_converted_}"

    out_file="$OUT_DIR/BMPMWG_seed_${suffix%.root}.out"
    log_file="$LOG_DIR/BMPMWG_seed_${suffix%.root}.log"
    seed_file="$SEED_DIR/BMPMWG_seed_$suffix"

    echo "Submitting job for:"
    echo "  BMPMWG_converted  : $converted_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  BMPMWG_seed: $seed_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-seed-mode \
        --input-file-path "$converted_file" \
        --output-file-path "$seed_file" \
        --log-file-path "$log_file" \
        --data-type data \
        --deadpln_bm_run15
done
