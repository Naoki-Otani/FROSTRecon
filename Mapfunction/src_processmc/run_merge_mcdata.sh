#!/usr/bin/env bash
set -euo pipefail

MACRO="merge_mcdata.C"

inDirArg=(
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/rawdata/x"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/rawdata/y"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/rawdata/x"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/rawdata/y"
)

outDirArg=(
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/mergedata/x"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/singlehit/mergedata/y"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/mergedata/x"
  "/group/nu/ninja/work/otani/FROSTReconData/MapFunc/mcdata/twohit/mergedata/y"
)

if [[ ${#inDirArg[@]} -ne ${#outDirArg[@]} ]]; then
  echo "[ERROR] inDirArg and outDirArg must have the same length." >&2
  exit 1
fi

for i in "${!inDirArg[@]}"; do
  in="${inDirArg[$i]}"
  out="${outDirArg[$i]}"

  echo "=== [$i] in=$in ==="
  echo "===     out=$out ==="

  if [[ ! -d "$in" ]]; then
    echo "[SKIP] input directory not found: $in"
    continue
  fi

  mkdir -p "$out"

  root -l -b -q "${MACRO}+(
    \"${in}\",
    \"${out}\"
  )"
done

echo "All done."
