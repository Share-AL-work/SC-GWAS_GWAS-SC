#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# Convert a UCSC‑style BED4 file of *continuous* scores into an LDSC .annot.gz
#   usage:  make_ldsc_continuous_annot.R  \
#              <annotation.bed>           \
#              <plink.bim>                \
#              <out.prefix>               \
#              <mode>                     #  "full-annot" (default) | "thin-annot"
# The BED must be 0‑based half‑open with four columns:
#   chr  start  end  value
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Need at least three arguments: BED, BIM, OUT_PREFIX [mode]")

bed_file  <- args[1]
bim_file  <- args[2]
out_pref  <- args[3]              # WITHOUT .annot.gz & WITHOUT chromosome suffix
mode      <- ifelse(length(args) == 4, args[4], "full-annot")   # thin‑annot | full-annot

# 1  read PLINK .bim → SNP locations
bim <- fread(bim_file,
             col.names = c("CHR","SNP","CM","BP","A1","A2"))
bim.gr <- GRanges(seqnames = bim$CHR,
                  ranges   = IRanges(bim$BP, bim$BP),
                  SNP      = bim$SNP,
                  CM       = bim$CM)

# 2  read BED4
bed <- fread(bed_file, col.names = c("chr","start","end","value"))
bed$chr   <- sub("^chr","", bed$chr)
bed$start <- bed$start + 1               # convert 0‑based → 1‑based
bed.gr <- GRanges(seqnames = bed$chr,
                  ranges   = IRanges(bed$start, bed$end),
                  value    = bed$value)

# 3  overlap → attach score
hits <- findOverlaps(bim.gr, bed.gr)
scores <- rep(0, length(bim.gr))
scores[queryHits(hits)] <- mcols(bed.gr)$value[subjectHits(hits)]

# 4  compose LDSC table
if (mode == "thin-annot") {
  out <- data.table(ANNOT = scores)
} else {
  out <- data.table(
          CHR  = bim$CHR,
          BP   = bim$BP,
          SNP  = bim$SNP,
          CM   = bim$CM,
          ANNOT = scores)
}

dir.create(dirname(out_pref), showWarnings = FALSE, recursive = TRUE)
fwrite(out, file = gzfile(paste0(out_pref, ".annot.gz")),
       sep="\t", quote=FALSE)
cat("✓ saved ", paste0(out_pref, ".annot.gz"), "\n")
