#!/bin/bash

# Merge FROST ROOT files with names like:
#   run00018_310000_319999_recon.root
# into:
#   run00018_recon.root
#
# Usage:
#   ./merge_FROSTdata.sh input_dir output_dir

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_dir output_dir"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

# Extract unique run names such as run00018
find "$INPUT_DIR" -maxdepth 1 -type f -name "run*_recon.root" | \
    sed -E 's#^.*/##' | \
    sed -E 's#^(run[0-9]+)_.*_recon\.root$#\1#' | \
    sort -u | \
while read -r runname; do
    mapfile -t files < <(
        find "$INPUT_DIR" -maxdepth 1 -type f \
            -name "${runname}_*_recon.root" | \
        sort -V
    )

    if [ "${#files[@]}" -eq 0 ]; then
        continue
    fi

    output_file="${OUTPUT_DIR}/${runname}_recon.root"

    echo "Merging ${runname} -> ${output_file}"
    printf '  %s\n' "${files[@]}"

    hadd -f "$output_file" "${files[@]}"
done

echo "Done."
