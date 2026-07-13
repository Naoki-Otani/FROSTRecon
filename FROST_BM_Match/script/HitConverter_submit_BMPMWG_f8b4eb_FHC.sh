#!/bin/bash

set -eu

# ========= path settings =========
BM_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMPMWGdata_f8b4eb/FHC"
FROST_DIR="/group/nu/ninja/work/otani/FROSTReconData/FROSTdata"

ROOT_OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BM_FROST_BMWGPM_f8b4eb/1-rootfile_after_HitConverter/FHC"
OUT_DIR="$ROOT_OUT_DIR/out"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
HITCONVERTER="/home/nu/notani/FROSTRecon/FROST_BM_Match/bin/HitConverter/HitConverter"

QUEUE="h"

# ========= prepare directories =========
mkdir -p "$ROOT_OUT_DIR" "$OUT_DIR"

# ========= loop over all BM files =========
for bm_file in "$BM_DIR"/*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$bm_file" ] || continue

    bm_base=$(basename "$bm_file")

    root_out_file="$ROOT_OUT_DIR/${bm_base%.root}_afterHitConverter.root"
    out_file="$OUT_DIR/${bm_base%.root}_afterHitConverter.out"

    echo "Submitting job for:"
    echo "  BM  : $bm_file"
    echo "  OUT : $out_file"
    echo "  output ROOT: $root_out_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$HITCONVERTER" \
        "$bm_file" \
        "$FROST_DIR" \
        "$root_out_file"
done
