#!/bin/bash

set -eu

# ========= path settings =========
CONVERTED_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter_divided"

SEED_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-BMPMrecon/seed"
OUT_DIR="$SEED_DIR/out"
LOG_DIR="$SEED_DIR/log"


SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/opt/wagasci_mc/WagasciTR/bin/WagasciRecon"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$SEED_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for converted_file in "$CONVERTED_DIR"/BMPM_converted_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$converted_file" ] || continue

    converted_base=$(basename "$converted_file")
    # Common part obtained by removing "BMPM_converted_"
    suffix="${converted_base#BMPM_converted_}"

    out_file="$OUT_DIR/BMPM_seed_${suffix%.root}.out"
    log_file="$LOG_DIR/BMPM_seed_${suffix%.root}.log"
    seed_file="$SEED_DIR/BMPM_seed_$suffix"

    echo "Submitting job for:"
    echo "  BMPM_converted  : $converted_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  BMPM_seed: $seed_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-seed-mode \
        --input-file-path "$converted_file" \
        --output-file-path "$seed_file" \
        --log-file-path "$log_file" \
        --data-type data
done
