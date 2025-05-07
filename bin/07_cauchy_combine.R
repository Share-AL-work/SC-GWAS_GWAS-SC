#!/usr/bin/env Rscript
# Combine P‑values from LDSC, MAGMA‑GSEA, scDRS via ACAT
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(yaml)
})

cfg_file <- commandArgs(trailingOnly = TRUE)[1] %||% "config.yml"
cfg <- yaml::read_yaml(cfg_file)

res_root <- cfg$paths$results_root
traits   <- cfg$gwas
metrics  <- c("ldsc_continuous","magma_gsea","scdrs_mbat")

ACATO <- function(p){
  p <- p[!is.na(p) & p>0 & p<1]
  if(!length(p)) return(NA_real_)
  t_stat <- sum(tan((0.5 - p) * pi))/length(p)
  1 - pcauchy(t_stat)
}

combined <- list()

for(trait in traits){
  for(ds in cfg$singlecell_datasets){
    pvec <- c(
      fread(file.path(res_root,"ldsc",  paste0(ds,"_",trait,".tsv")))$P,
      fread(file.path(res_root,"magma", paste0(ds,"_",trait,".tsv")))$P,
      fread(file.path(res_root,"scdrs_scores",ds,trait,"summary.tsv"))$P
    )
    combined[[length(combined)+1]] <-
      data.table(dataset=ds, trait=trait, P_combo = ACATO(pvec))
  }
}

fwrite(rbindlist(combined),
       file = file.path(res_root,"pvalue_cauchy_combined.tsv"),
       sep="\t")
cat("✓ saved combined table\n")
