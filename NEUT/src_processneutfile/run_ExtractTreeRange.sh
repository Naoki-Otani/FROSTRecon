#!/usr/bin/env bash
set -euo pipefail

INPUT="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/NEUToutput/rawroot/neut_320kA_1.077e21pot_ECC12.root"
OUTDIR="/group/nu/ninja/work/otani/FROSTReconData/FROST_NEUT/NEUToutput/dividedroot"

# Range settings: (0,99), (100,199), ... , (92200,92299)
START=0
END_MAX=92299
STEP=100

mkdir -p "$OUTDIR"

for ((first=START; first<=END_MAX; first+=STEP)); do
  last=$((first + STEP - 1))

  echo "Extracting entries ${first}-${last} ..."
  root -l -b -q "ExtractTreeRange.C(\"${INPUT}\",\"${OUTDIR}\",${first},${last})"
done

echo "Done."
