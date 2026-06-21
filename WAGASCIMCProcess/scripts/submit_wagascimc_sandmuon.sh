#!/bin/bash

set -eu

# ========= path settings =========
NEUT_OUTPUT_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/0-Neut"

OUTPUT_ROOT_DIR="/group/nu/ninja/work/otani/FROSTReconData/SandmuonMC/1-WagasciMC"
OUT_DIR="$OUTPUT_ROOT_DIR/out"
LOG_DIR="$OUTPUT_ROOT_DIR/log"
OUT_GDML_DIR="$OUTPUT_ROOT_DIR/gdml"

GEOMETRY_DIR="/home/nu/notani/wagasci-babymind-monte-carlo/data/geometry"


SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
B2MC="/home/nu/notani/wagasci-babymind-monte-carlo/bin/B2MC"

QUEUE="l"

# ========= prepare directories =========
mkdir -p "$OUTPUT_ROOT_DIR" "$OUT_DIR" "$LOG_DIR" "$OUT_GDML_DIR"

# ========= loop over all neut output files ==========
for neut_output_file in "$NEUT_OUTPUT_DIR"/genev.*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$neut_output_file" ] || continue

    neut_output_base=$(basename "$neut_output_file")

    out_file="$OUT_DIR/b2mc_sandmuon_${neut_output_base%.root}.out"
    log_file="$LOG_DIR/b2mc_sandmuon_${neut_output_base%.root}.log"
    output_root_file="$OUTPUT_ROOT_DIR/b2mc_sandmuon_${neut_output_base%.root}.root"
    output_gdml_file="$OUT_GDML_DIR/b2mc_sandmuon.gdml"

    echo "Submitting job for:"
    echo "  INPUT_ROOT_FILE : $neut_output_file"
    echo "  OUT : $out_file"
    echo "  LOG : $log_file"
    echo "  OUTPUT_ROOT_FILE : $output_root_file"
    echo "  OUTPUT_GDML_FILE : $output_gdml_file"
    bsub -q "$QUEUE" -o "$out_file" -N \
        apptainer run --bind /hsm/nu/ --bind /group/nu/ "$SIF" \
        "$B2MC" \
        --input-file-path "$neut_output_file" \
        --neutrino-interaction-generator NEUTGEOM \
        --geometry-dir-path "$GEOMETRY_DIR" \
        --output-file-path "$output_root_file" \
        --output-gdml-file-path "$output_gdml_file" \
        --log-file-path "$log_file" \
        --ninja-geometry-version FROST
done
