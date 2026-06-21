#!/bin/bash

set -eu

# ========= path settings =========
WAGASCIMC_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/1-WagasciMC"

SEED_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/2-recon/seed"
OUT_DIR="$SEED_DIR/out"
LOG_DIR="$SEED_DIR/log"


SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
RECON="/home/nu/notani/wagasci-babymind-track-reconstruction/bin/WagasciRecon"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$SEED_DIR" "$OUT_DIR" "$LOG_DIR"

# ========= loop over all converted files =========
for wagascimc_file in "$WAGASCIMC_DIR"/b2mc_sandmuon_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$wagascimc_file" ] || continue

    converted_base=$(basename "$wagascimc_file")
    # Common part obtained by removing "BMPMWG_converted_"
    suffix="${converted_base#b2mc_sandmuon_}"

    out_file="$OUT_DIR/sandmuonmc_seed_${suffix%.root}.out"
    log_file="$LOG_DIR/sandmuonmc_seed_${suffix%.root}.log"
    seed_file="$SEED_DIR/sandmuonmc_seed_$suffix"

    echo "Submitting job for:"
    echo "  WAGASCIMC file  : $wagascimc_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  seed: $seed_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$RECON" \
        --reconstruction-seed-mode \
        --input-file-path "$wagascimc_file" \
        --output-file-path "$seed_file" \
        --log-file-path "$log_file" \
        --data-type mc \
        --copy-ttrees-except-tree
done
