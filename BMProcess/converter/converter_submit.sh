#!/bin/bash

set -eu

# ========= path settings =========
# BM_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD"
BM_DIR="/hsm/nu/wagasci/dhirata/bm_recon/2-BMBSD"
PM_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/4-PMBSD_BMformat"
OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/out"
LOG_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/log"
ROOT_OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter"
TMP_DIR="/group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/tmp"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
CONVERTER="/opt/wagasci_data_handling/WagasciConverter/bin/Converter"
GEOMETRY_DIR="/opt/wagasci_mc/WagasciMC/etc/wagasci/b2/geometry"

QUEUE="s"

# ========= prepare directories =========
mkdir -p "$OUT_DIR" "$LOG_DIR" "$ROOT_OUT_DIR" "$TMP_DIR"

# ========= loop over all BMBSD files =========
for bm_file in "$BM_DIR"/BMBSD_*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$bm_file" ] || continue

    bm_base=$(basename "$bm_file")

    # Common part obtained by removing "BMBSD_"
    suffix="${bm_base#BMBSD_}"

    pm_file="$PM_DIR/PMBSD_$suffix"
    out_file="$OUT_DIR/Converter_BMPM_${suffix%.root}.out"
    log_file="$LOG_DIR/Converter_BMPM_${suffix%.root}.log"
    root_out_file="$ROOT_OUT_DIR/BMPM_converted_$suffix"

    if [ ! -f "$pm_file" ]; then
        echo "Skip: corresponding PM file not found for $bm_base"
        echo "      expected: $pm_file"
        continue
    fi

    echo "Submitting job for:"
    echo "  BM  : $bm_file"
    echo "  PM  : $pm_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  ROOT: $root_out_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$CONVERTER" \
        --output-file-path "$root_out_file" \
        --proton-module \
        --proton-module-path "$pm_file" \
        --baby-mind \
        --baby-mind-path "$bm_file" \
        --log-file-path "$log_file" \
        --temporary-dir-path "$TMP_DIR" \
        --geometry-dir-path "$GEOMETRY_DIR"
done

# bsub -q s -o /group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/out/Converter_BMPM_2025-12-21_15-16-32_Run0.out -N apptainer run --bind /hsm/nu/ --bind /group/nu/ /home/nu/notani/wagasci_ana_0.2.4.sif /opt/wagasci_data_handling/WagasciConverter/bin/Converter --output-file-path /group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/BMPM_converted_2025-12-21_15-16-32_Run0.root --baby-mind --baby-mind-path /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD/BMBSD_2025-12-21_15-16-32_Run0.root --log-file-path /group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/log/Converter_BMPM_2025-12-21_15-16-32_Run0.log --temporary-dir-path /group/nu/ninja/work/otani/FROSTReconData/BMdata/3-BMPMconverter/tmp/ --geometry-dir-path /opt/wagasci_mc/WagasciMC/etc/wagasci/b2/geometry
