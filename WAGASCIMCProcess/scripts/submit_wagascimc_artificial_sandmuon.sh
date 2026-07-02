#!/bin/bash

set -eu

OUTPUT_ROOT_DIR="/group/nu/ninja/work/otani/FROSTReconData/Artificial_sandmuonMC/1-WagasciMC"
OUT_DIR="$OUTPUT_ROOT_DIR/out"
LOG_DIR="$OUTPUT_ROOT_DIR/log"
OUT_GDML_DIR="$OUTPUT_ROOT_DIR/gdml"

GEOMETRY_DIR="/home/nu/notani/wagasci-babymind-monte-carlo/data/geometry"

SIF="/home/nu/notani/wagasci_ana_0.2.4.sif"
B2MC="/home/nu/notani/wagasci-babymind-monte-carlo/bin/B2MC"

QUEUE="l"
NUMBER_OF_EVENTS="1000"

mkdir -p "$OUTPUT_ROOT_DIR" "$OUT_DIR" "$LOG_DIR" "$OUT_GDML_DIR"

saved_gdml_file="$OUT_GDML_DIR/b2mc_artificial_sandmuon.gdml"

for i in {1..1000}; do
    out_file="$OUT_DIR/b2mc_artificial_sandmuon_${i}.out"
    log_file="$LOG_DIR/b2mc_artificial_sandmuon_${i}.log"
    output_root_file="$OUTPUT_ROOT_DIR/b2mc_artificial_sandmuon_${i}.root"

    if [ "$i" -eq 1 ]; then
        output_gdml_file="$saved_gdml_file"
        remove_gdml_command=":"
    else
        output_gdml_file="$OUT_GDML_DIR/b2mc_artificial_sandmuon_${i}.tmp.gdml"
        remove_gdml_command="rm -f \"$output_gdml_file\""
    fi

    echo "Submitting job for:"
    echo "  OUT              : $out_file"
    echo "  LOG              : $log_file"
    echo "  OUTPUT_ROOT_FILE : $output_root_file"
    echo "  OUTPUT_GDML_FILE : $output_gdml_file"

    bsub -q "$QUEUE" -o "$out_file" -N \
        bash -lc "
            set -e
            apptainer run --bind /hsm/nu/ --bind /group/nu/ --bind /home/nu/ \"$SIF\" \
                \"$B2MC\" \
                --frost-sand-muons \
                --number-of-events \"$NUMBER_OF_EVENTS\" \
                --geometry-dir-path \"$GEOMETRY_DIR\" \
                --output-file-path \"$output_root_file\" \
                --output-gdml-file-path \"$output_gdml_file\" \
                --log-file-path \"$log_file\" \
                --ninja-geometry-version FROST
            $remove_gdml_command
        "
done
