#!/usr/bin/env bash
# Chain: mBAT geneassoc → TSV → GS
set -euo pipefail
cfg=${1:-config.yml}

# 1. TSV
Rscript scripts/scdrs/mbat_to_tsv.R "$cfg"

# 2. GS
bash     scripts/scdrs/tsv_to_gs.sh  "$cfg"
