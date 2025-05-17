# SC-GWAS:GWAS-SC pipeline
This repo ties together:

* **a Cauchy combination analysis layer for below:**
* **mBAT‑combo based scDRS**
* **Binary‑sLDSC** (unmodified fork of <https://github.com/jbryois/scRNA_disease>)  
* **Continuous‑sLDSC** (Bash + R around pascaltimshel/ldsc @ d869cfd)  
* **MAGMA‑GSEA**(binary annotation,unmodified fork of <https://github.com/jbryois/scRNA_disease>),

## 0.1.Preprocess the GWAS summary data:
* COJO format: follow the section in https://github.com/zlintian/SBayesRC_pipeline
* LDSC format: follow https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial#step-0-create-conda-environment-for-munging

## 0.2.Preprocess the sc(n)RNAseq data:
You can find full details of the preprocessing step in the [preprocess documentation](https://github.com/Share-AL-work/SC-GWAS_GWAS-SC/tree/main/preprocess).


All absolute paths are replaced by a single editable `config.yml`.

#### Quick start (all stages) with Installation

```bash
git clone --recurse-submodules https://github.com:Share-AL-work/SC-GWAS_GWAS-SC.git
cd SC-GWAS_GWAS-SC
conda env create -f envs/ldsctools.yml
conda activate ldsctools
################################################################################
# 0.  (only once) create Conda environments
################################################################################
conda env create -f envs/metrics.yml     # R + python + cellex + Cepo
conda env create -f envs/ldsctools.yml   # R + python2 + bedtools + gcta
conda env create -f envs/scdrs.yml       # scDRS / scanpy stack

################################################################################
# 1.  Compute metrics  (repeat with any metric you want, even with new metrics)
################################################################################
conda activate metrics_env

bin/01_metrics.sh  config.yml  cepo      # Cepo
bin/01_metrics.sh  config.yml  cellex    # DET/NSI/GES/EP/ESmu
bin/01_metrics.sh  config.yml  tdep      # Top‑decile EP
bin/01_metrics.sh  config.yml  sclinker  # sc‑linker gene programs

################################################################################
# 2.  Convert each metrics‑CSV to BED4  (one‑off per dataset×metric)
################################################################################
Rscript bin/02_csv_to_bed.R  config.yml

################################################################################
# 3.  LDSC
################################################################################
conda activate ldsctools

# 3a. binary s‑LDSC (unmodified scRNA_disease pipeline)
bin/03_ldsc_binary.sh      config.yml

# 3b. continuous s‑LDSC (fast fork, Kevin‑Luo annot)
bin/03_ldsc_continuous.sh  config.yml

################################################################################
# 4.  MAGMA‑GSEA for binary metrics (reuse Bryois wrapper)
################################################################################
bin/04_magma_gsea.sh  config.yml      # (creates results/magma/...)

################################################################################
# 5.  mBAT‑combo 
# Manual of mBAT-combo:
# https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo
################################################################################
bin/05_mbat_combo.sh   config.yml     # *.gene.assoc.mbat into results/mbat

# 5b. convert mBAT output → TSV → weighted .gs   (for scDRS)
bin/05b_mbat_to_gs.sh  config.yml     # results/scdrs_gs/<trait>.gs

################################################################################
# 6.  scDRS (uses .gs from 5b)
# Manual of scDRS:
# https://martinjzhang.github.io/scDRS/
################################################################################
conda activate scdrs_env
bin/06_scdrs.sh        config.yml     # results/scdrs_scores/...

################################################################################
# 7.  Cauchy combination (LDSC + MAGMA + scDRS)
################################################################################
conda activate metrics_env
Rscript bin/07_cauchy_combine.R  config.yml
# --> results/pvalue_cauchy_combined.tsv
################################################################################
echo "✓ Pipeline finished — meta P‑values in results/pvalue_cauchy_combined.tsv"
```

####  Code to reproduce results of the paper is in Replicates_paper directory



























