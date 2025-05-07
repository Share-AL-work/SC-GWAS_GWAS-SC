#!/usr/bin/env Rscript
# -------------------------------------------------------------
# Run EWCE‑TDEP metrics on every .h5ad file in DATA_DIR
# and save results to OUT_DIR/<groupName>/...
#
#  * No absolute paths – fully repo‑portable
#  * Works for both “all genes” and “*_pc.h5ad” files
#  * Adapt cell‑type column if needed (default: cell_type)
# -------------------------------------------------------------
suppressPackageStartupMessages({
  library(zellkonverter)   # readH5AD
  library(SingleCellExperiment)
  library(EWCE)
  library(DelayedArray)
  library(Matrix)
})

# ─────────────────────────────────────────────────────────────
# 1. CONFIG – edit if your repo uses different folders/column
# ─────────────────────────────────────────────────────────────
DATA_DIR <- "data"          # where your .h5ad files live
OUT_DIR  <- "tdep_out"      # where EWCE objects will be written
CELLTYPE_COL <- "cell_type" # column in colData with cell‑type labels
# ─────────────────────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Searching for .h5ad files in ", normalizePath(DATA_DIR))

h5_files <- list.files(DATA_DIR, pattern = "\\.h5ad$", full.names = TRUE)
if (length(h5_files) == 0)
  stop("No .h5ad files found in ", DATA_DIR)

process_one <- function(h5) {
  grp <- tools::file_path_sans_ext(basename(h5))  # groupName
  message("• Processing: ", grp)
  
  sce <- zellkonverter::readH5AD(h5)
  if (!CELLTYPE_COL %in% colnames(colData(sce)))
    stop("Column '", CELLTYPE_COL, "' not found in ", h5)
  
  # Ensure a counts layer exists for EWCE (dense matrix is fine)
  assay(sce, "counts") <- assay(sce, "X")
  
  # Convert DelayedArray to ordinary matrix for EWCE
  counts_mat <- as.matrix(counts(sce))
  
  # Clean cell‑type labels (spaces → _, hyphens → _)
  ct <- gsub("[\\s,-]", "_", colData(sce)[[CELLTYPE_COL]])
  colData(sce)[[CELLTYPE_COL]] <- ct
  
  # Annotation list (single level l1)
  annot_levels <- list(l1 = ct)
  
  # Run EWCE helper
  EWCE::generate_celltype_data(
    exp        = counts_mat,
    annotLevels= annot_levels,
    groupName  = grp,
    savePath   = OUT_DIR
  )
  message("  ↳ saved to ", file.path(OUT_DIR, grp))
}

# Main loop ----------------------------------------------------
for (f in h5_files) {
  tryCatch(process_one(f),
           error = function(e) message("  ⚠️  ", e$message))
}

message("Completed all datasets → ", normalizePath(OUT_DIR))
