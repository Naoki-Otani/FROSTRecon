#!/bin/bash

set -eu

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
MACRO="divide_tree_every_1000.C"

OUT_DIR="$OUTPUT_DIR/out"
QUEUE="s"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory does not exist: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$OUT_DIR"

found_file=0

for input_file in "$INPUT_DIR"/*.root; do
    # Safeguard for the case where the glob matches nothing
    [ -e "$input_file" ] || continue

    found_file=1

    echo "Processing: $input_file"
    bsub -q "$QUEUE" -o "$OUT_DIR"/$(basename "$input_file" .root).out -N root -l -b -q "${MACRO}(\"${input_file}\",\"${OUTPUT_DIR}\")"
done

if [ "$found_file" -eq 0 ]; then
    echo "No ROOT files found in: $INPUT_DIR"
    exit 1
fi

echo "Done."
