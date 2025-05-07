#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Convert mBAT‑combo *.gene.assoc.mbat files to scDRS‑ready TSV
#   • filters for protein‑coding genes (optional)
#   • keeps top ≤ 1000 genes
#   • replaces P = 0 with smallest non‑zero
# usage:  mbat_to_tsv.R  <config.yml>
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(yaml); library(stats)
})

acat <- function(p){ qnorm(p/2, lower.tail = FALSE)*sign(p-0.5) } # p → Z

cfg   <- yaml::read_yaml(commandArgs(trailingOnly = TRUE)[1])
mbat_dir <- file.path(cfg$paths$results_root, "mbat")
tsv_dir  <- file.path(cfg$paths$results_root, "scdrs_tsv")
dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)

# optional protein‑coding list
pc_genes <- fread(cfg$gene_coord)$Gene

for(trait in cfg$gwas){
  f_in <- file.path(mbat_dir, paste0(trait, ".gene.assoc.mbat"))
  if(!file.exists(f_in)){ warning("missing ", f_in); next }
  dt <- fread(f_in)
  dt <- dt[Gene %in% pc_genes]
  dt[, P_mBATcombo := ifelse(P_mBATcombo == 0, min(P_mBATcombo[P_mBATcombo>0]), P_mBATcombo)]
  out <- dt[order(P_mBATcombo)][1:min(.N,1000),
         .(GENE = Gene, Z = acat(P_mBATcombo))]
  fwrite(out, file = file.path(tsv_dir, paste0(trait, ".tsv")), sep="\t")
  message("✓ TSV for ", trait)
}
