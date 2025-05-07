#!/usr/bin/env bash
set -euo pipefail
cfg=${1:-config.yml}
metric=${2:-cepo}          # cepo | cellex | tdep | sclinker

DATA_ROOT=$(yq '.paths.data_root' "$cfg")
OUTDIR=$(yq '.paths.results_root' "$cfg")/metrics_csv

mkdir -p "$OUTDIR"
case $metric in
  cepo)       Rscript scripts/metrics/run_cepo.R      "$cfg" ;;
  cellex)     python  scripts/metrics/run_cellex.py   "$cfg" ;;
  tdep)       Rscript scripts/metrics/run_tdep_all.R  "$cfg" ;;
  sclinker)   python  scripts/metrics/run_sclinker.py "$cfg" ;;
  myMetric)   Rscript scripts/metrics/run_myMetric.R  "$cfg" ;;   # <── NEW
  *) echo "unknown metric"; exit 1 ;;
esac
