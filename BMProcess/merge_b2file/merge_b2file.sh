#!/bin/bash

# Merge split ROOT files such as:
#   BMPM_track_2025-12-15_13-58-44_Run0_0.root
#   BMPM_track_2025-12-15_13-58-44_Run0_1.root
#   ...
# into:
#   BMPM_track_2025-12-15_13-58-44_Run0.root
#
# Usage:
#   ./merge_b2file.sh input_dir output_dir

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_dir output_dir"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"

# Find files matching *_<number>.root, extract the common base name,
# and process each base only once.
find "$INPUT_DIR" -maxdepth 1 -type f -name "*.root" | \
    sed -E 's#^.*/##' | \
    sed -E 's/_[0-9]+\.root$//' | \
    sort -u | \
while read -r base; do
    # Collect files in numeric order: _0.root, _1.root, _2.root, ...
    mapfile -t files < <(
        find "$INPUT_DIR" -maxdepth 1 -type f \
            -regextype posix-extended \
            -regex ".*/${base}_[0-9]+\.root" | \
        sort -V
    )

    if [ "${#files[@]}" -eq 0 ]; then
        continue
    fi

    output_file="${OUTPUT_DIR}/${base}.root"

    echo "Merging -> ${output_file}"
    printf '  %s\n' "${files[@]}"

    hadd -f "$output_file" "${files[@]}"
done

echo "Done."
