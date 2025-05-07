#!/usr/bin/env bash
# Compute scDRS scores for every (dataset × trait) pair
set -euo pipefail
cfg=${1:-config.yml}

SC_DATA=$(yq '.paths.data_root' $cfg)/singlecell
GS_DIR=$(yq '.paths.results_root' $cfg)/scdrs_gs
SCORE_DIR=$(yq '.paths.results_root' $cfg)/scdrs_scores
NCTRL=$(yq '.n_ctrl_scdrs' $cfg)

mkdir -p "$GS_DIR" "$SCORE_DIR"

for ds in $(yq '.singlecell_datasets[]' $cfg); do
  h5="$SC_DATA/${ds}.h5ad"
  for trait in $(yq '.gwas[]' $cfg); do
    gs_file="${GS_DIR}/${trait}.gs"
    echo "→ scDRS   $ds  ×  $trait"
    scdrs compute-score \
        --h5ad-file "$h5" \
        --h5ad-species human \
        --gs-file "$gs_file" \
        --gs-species human \
        --out-folder "$SCORE_DIR/${ds}/${trait}" \
        --n-ctrl "$NCTRL" \
        --flag-filter-data True \
        --flag-raw-count True
  done
done
