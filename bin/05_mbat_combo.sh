#!/usr/bin/env bash
# Run mBAT‑combo for every GWAS listed in config.yml
set -euo pipefail
cfg=${1:-config.yml}

GWAS_DIR=$(yq '.paths.data_root' $cfg)/gwas
OUTDIR=$(yq '.paths.results_root' $cfg)/mbat
GENE_COORD=$(yq '.gene_coord' $cfg)
WINDOW=$(yq '.window_kb' $cfg)

mkdir -p "$OUTDIR"

for trait in $(yq '.gwas[]' $cfg); do
  echo "→ mBAT‑combo  $trait"
  gcta64 --mbat \
         --bfile   "$(yq '.paths.ref' $cfg)/1000G_EUR_Phase3_plink/1000G.EUR.QC" \
         --sumstats "${GWAS_DIR}/${trait}.sumstats.gz" \
         --gene-coord "$GENE_COORD" \
         --out "${OUTDIR}/${trait}" \
         --window-kb "$WINDOW"
done
