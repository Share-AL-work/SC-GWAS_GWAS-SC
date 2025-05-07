#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# Convert a "gene × annotations" CSV into BED4 files (one BED per column)
# usage:
#   csv_matrix_to_bed.R  \
#       <matrix.csv.gz> <gene_coords.tsv> <window_kb> <out-dir>
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Need 4 args: matrix.csv gene_coords.tsv window_kb out_dir\n")
}

matrix_csv <- args[1]
coord_tsv  <- args[2]
window_kb  <- as.integer(args[3])
outdir     <- args[4]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── 1  read files ──────────────────────────────────────────────────────────
mat    <- fread(matrix_csv)             # gene, anno1, anno2, …
coords <- fread(coord_tsv)              # gene, chr, start, end, …

setnames(coords, 1:4, c("gene","chr","start","end"))

# ± window
if (window_kb > 0) {
  coords[, `:=`(start = pmax(start - window_kb*1000, 0),
                end   = end   + window_kb*1000)]
}

# ── 2  merge ───────────────────────────────────────────────────────────────
merged <- merge(coords, mat, by = "gene", all.x = FALSE)
annot_cols <- setdiff(names(merged), c("gene","chr","start","end"))

# ── 3  write BED per annotation ────────────────────────────────────────────
for (ann in annot_cols) {
  bed <- merged[, .(chr,
                    start = start - 1,   # 0‑based
                    end,
                    value = get(ann))]
  fwrite(bed[value != 0], file = file.path(outdir, paste0(ann, ".bed")),
         sep = "\t", col.names = FALSE)
  cat("✓ wrote", file.path(outdir, paste0(ann, ".bed")), "\n")
}
