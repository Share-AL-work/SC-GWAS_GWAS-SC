#!/usr/bin/env bash
# Make .gs from every TSV in scdrs_tsv/
set -euo pipefail
cfg=${1:-config.yml}

TSV_DIR=$(yq '.paths.results_root' $cfg)/scdrs_tsv
GS_DIR=$(yq '.paths.results_root' $cfg)/scdrs_gs
mkdir -p "$GS_DIR"

for tsv in "$TSV_DIR"/*.tsv; do
  [ -e "$tsv" ] || { echo "No TSV in $TSV_DIR"; exit 1; }
  trait=$(basename "${tsv%.tsv}")
  scdrs munge-gs \
        --out_file   "$GS_DIR/${trait}.gs" \
        --zscore-file "$tsv" \
        --weight zscore \
        --n-max 1000
  echo "✓ GS for $trait"
done
