#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# Build an LDSC *.M* file for a continuous annotation
#  * works chromosome‑wise or joint (all chr in one .annot.gz each)
#  * CLI mirrors the Perl script but is shorter
# usage:
#   make_quantile_M.R  \
#       --annot-prefix results/annot/dataset/metric.          # <prefix>.<chr>.annot.gz
#       --baseline-prefix data/ldsc/baseline_v1.1_thin_annot/baseline.
#       --frq-prefix data/ldsc/1000G_Phase3_frq/1000G.EUR.QC.   # needed only for MAF filter
#       --nb-quantile 5          # (default)
#       --exclude-zero TRUE      # drop SNPs with metric==0 before binning
#       --out  results/metric.q5_exclude_zero.q_M
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(R.utils)
})

opt <- list(
  make_option("--annot-prefix",      type="character"),
  make_option("--baseline-prefix",   type="character"),
  make_option("--frq-prefix",        type="character", default=NULL),
  make_option("--nb-quantile",       type="integer",   default=5),
  make_option("--exclude-zero",      type="logical",   default=TRUE),
  make_option("--maf",               type="double",    default=0.05),
  make_option("--out",               type="character"))
parser <- OptionParser(option_list = opt)
args   <- parse_args(parser)

if (is.null(args$`annot-prefix`) | is.null(args$out))
  stop("annot-prefix and out are mandatory")

chr_list <- sprintf("%d", 1:22)
# ── 1  read annotation for every chr ─────────────────────────────────────────
read_annot <- function(pref, chr) {
  f <- paste0(pref, chr, ".annot.gz")
  dt <- fread(cmd=paste("zcat", shQuote(f)))
  dt[, chr := chr]        # keep for later join
  dt
}

annot_dt <- rbindlist(lapply(chr_list, read_annot, pref=args$`annot-prefix`))

focal_name <- tail(names(annot_dt), 1)    # last column is your metric

# ── 2  optionally merge baseline annotations ────────────────────────────────
if (!is.null(args$`baseline-prefix`)) {
  base_dt <- rbindlist(lapply(chr_list, read_annot,
                              pref=args$`baseline-prefix`))
  base_dt <- base_dt[, !(c("CHR","BP","SNP","CM")), with=FALSE]
  annot_dt <- cbind(annot_dt, base_dt)
}

# ── 3  MAF filter (optional) ───────────────────────────────────────────────
if (!is.null(args$`frq-prefix`)) {
  frq <- rbindlist(lapply(chr_list, function(chr) {
    fread(paste0(args$`frq-prefix`, chr, ".frq"),
          col.names=c("CHR","A1","A2","MAF","NCHROBS","SNP"))
  }))
  annot_dt <- merge(annot_dt, frq[,.(SNP, MAF)], by="SNP")
  annot_dt <- annot_dt[MAF >= args$maf]
}

# ── 4  choose quantile boundaries ───────────────────────────────────────────
x <- annot_dt[[focal_name]]
if (args$`exclude-zero`) {
  bins <- quantile(x[x!=0], probs = seq(0,1,length.out=args$`nb-quantile`+1),
                   na.rm=TRUE)
  annot_dt <- annot_dt[x!=0]           # remove zeros before binning
} else {
  bins <- quantile(x, probs = seq(0,1,length.out=args$`nb-quantile`+1),
                   na.rm=TRUE)
}

annot_dt[, bin := cut(get(focal_name), breaks=bins,
                      include.lowest=TRUE, labels=FALSE)]

# ── 5  build the M matrix ───────────────────────────────────────────────────
cols <- setdiff(names(annot_dt), c("CHR","BP","SNP","CM","bin"))
M <- annot_dt[, lapply(.SD, sum), by=.(bin), .SDcols=cols]
setcolorder(M, c("bin", cols))

# prepend header lines exactly like Perl script
header <- c("## Intervals:", paste(bins, collapse=","))
writeLines(header, con = gzfile(args$out))
fwrite(M, file=gzfile(args$out), sep="\t", append=TRUE)
cat("✓ wrote", args$out, "\n")
