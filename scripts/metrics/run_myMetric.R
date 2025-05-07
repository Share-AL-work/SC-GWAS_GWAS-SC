#!/usr/bin/env Rscript
# minimal template — replace with real code
library(data.table); library(yaml)

cfg <- yaml::read_yaml(commandArgs(trailingOnly = TRUE)[1])
sc_root   <- file.path(cfg$paths$data_root, "singlecell")
out_root  <- file.path(cfg$paths$results_root, "metrics_csv")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (ds in cfg$singlecell_datasets) {
  sce_file <- file.path(sc_root, paste0(ds, ".h5ad"))
  # … compute your metric -> data.frame(gene, celltype1, celltype2, …)
  out_csv  <- file.path(out_root, paste0(ds, ".myMetric.csv.gz"))
  fwrite(metric_df, out_csv, sep = ",")
  message("✓ wrote ", out_csv)
}
