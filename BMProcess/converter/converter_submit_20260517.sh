#!/bin/bash

set -eu

# ========= path settings =========
BM_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/1-data_before_converter/Daily/BM"
PM_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/1-data_before_converter/Daily/PM"
WG_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/1-data_before_converter/Daily/WG"
ROOT_OUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/WGBMdata/2-BMPMWGconverter"
OUT_DIR="$ROOT_OUT_DIR/out"
LOG_DIR="$ROOT_OUT_DIR/log"
TMP_DIR="$ROOT_OUT_DIR/tmp"

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
    suffix_no_ext="${suffix%.root}"
    wg_directory="$WG_DIR/physics_run_${suffix_no_ext}_0"
    out_file="$OUT_DIR/Converter_BMPMWG_${suffix%.root}.out"
    log_file="$LOG_DIR/Converter_BMPMWG_${suffix%.root}.log"
    root_out_file="$ROOT_OUT_DIR/BMPMWG_converted_$suffix"

    if [ ! -f "$pm_file" ]; then
        echo "Skip: corresponding PM file not found for $bm_base"
        echo "      expected: $pm_file"
        continue
    fi

    # ========= WG option control =========
    wg_args=()
    use_wg=1

    case "$bm_base" in
        BMBSD_2025-12-03_00-00-00.root|\
        BMBSD_2026-01-16_00-00-00.root|\
        BMBSD_2026-01-20_00-00-00.root|\
        BMBSD_2026-01-21_00-00-00.root)
            use_wg=0
            ;;
    esac

    if [ "$use_wg" -eq 1 ]; then
        if [ ! -d "$wg_directory" ]; then
            echo "Skip: corresponding WG directory not found for $bm_base"
            echo "      expected: $wg_directory"
            continue
        fi

        wg_args=(
            --wagasci
            --wagasci-path "$wg_directory"
        )
    fi

    echo "Submitting job for:"
    echo "  BM  : $bm_file"
    echo "  PM  : $pm_file"

    if [ "$use_wg" -eq 1 ]; then
        echo "  WG  : $wg_directory"
    else
        echo "  WG  : disabled for this BM file"
    fi

    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  ROOT: $root_out_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$CONVERTER" \
        --output-file-path "$root_out_file" \
        --proton-module \
        --proton-module-path "$pm_file" \
        "${wg_args[@]}" \
        --baby-mind \
        --baby-mind-path "$bm_file" \
        --log-file-path "$log_file" \
        --temporary-dir-path "$TMP_DIR" \
        --geometry-dir-path "$GEOMETRY_DIR"
done
