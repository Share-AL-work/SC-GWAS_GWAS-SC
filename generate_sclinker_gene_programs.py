#!/usr/bin/env python
"""
Generate cell‑type‑enriched gene‑program matrices (p‑value, log‑fold‑change, score) for every .h5ad file in DATA_DIR.

Re‑implementation of 1.generateCelltypePrograms.ipynb (scgenetics,https://github.com/karthikj89/scgenetics/blob/master/src/1.generateCelltypePrograms.ipynb):
• No absolute paths
• Only uses the 'cell_type' column (edit CELLTYPE_COL if needed)
• Outputs CSVs to OUT_DIR, one set per dataset
"""

from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np

# ─────────────────────────────────────────────────────────────
# 1. CONFIG – edit if your repo uses different folders/column
# ─────────────────────────────────────────────────────────────
DATA_DIR     = Path("data")        # where the .h5ad files live
OUT_DIR      = Path("gene_programs")  # where CSVs will be written
CELLTYPE_COL = "cell_type"         # column with cell‑type labels
MIN_CELLS    = 10                  # min cells per cell‑type to keep
# ─────────────────────────────────────────────────────────────

OUT_DIR.mkdir(parents=True, exist_ok=True)


def write_matrices(adata: sc.AnnData, prefix: str) -> None:
    """
    Build 3 DataFrames (pval, logfold, score) with shape genes × cell‑types
    and write them as CSV files under OUT_DIR.
    """
    genes = list(adata.var_names)
    gene2idx = {g: i for i, g in enumerate(genes)}

    cell_subsets = adata.obs[CELLTYPE_COL].cat.categories
    cell2idx = {ct: i for i, ct in enumerate(cell_subsets)}

    # init matrices
    pval   = np.zeros((len(genes), len(cell_subsets)))
    logfc  = np.zeros_like(pval)
    score  = np.zeros_like(pval)

    de_key = f"{CELLTYPE_COL}_DE"

    for ct in cell_subsets:
        p_vals  = adata.uns[de_key]["pvals_adj"][ct]
        lfc     = adata.uns[de_key]["logfoldchanges"][ct]
        scores  = adata.uns[de_key]["scores"][ct]
        genes_ct = adata.uns[de_key]["names"][ct]

        for g, p, lf, sc_ in zip(genes_ct, p_vals, lfc, scores):
            idx_row = gene2idx.get(g)
            if idx_row is not None:
                idx_col = cell2idx[ct]
                pval[idx_row, idx_col]  = p
                logfc[idx_row, idx_col] = lf
                score[idx_row, idx_col] = sc_

    # DataFrames with index=genes, cols=celltype_L1
    colnames = [f"{ct}_L1" for ct in cell_subsets]
    (pd.DataFrame(pval,  index=genes, columns=colnames)
       .to_csv(OUT_DIR / f"{prefix}_pval.csv"))
    (pd.DataFrame(logfc, index=genes, columns=colnames)
       .to_csv(OUT_DIR / f"{prefix}_logfold.csv"))
    (pd.DataFrame(score, index=genes, columns=colnames)
       .to_csv(OUT_DIR / f"{prefix}_score.csv"))

    print(f"✓  {prefix}: matrices written ({len(genes)} genes × {len(colnames)} cts)")


def process_file(h5: Path):
    prefix = h5.stem
    print(f"\n• Processing: {prefix}")

    adata = sc.read_h5ad(h5)

    if CELLTYPE_COL not in adata.obs:
        raise KeyError(f"Column '{CELLTYPE_COL}' not found in {h5.name}")

    # ensure categorical
    adata.obs[CELLTYPE_COL] = adata.obs[CELLTYPE_COL].astype("category")

    # filter cell‑types with < MIN_CELLS
    ct_counts = adata.obs[CELLTYPE_COL].value_counts()
    keep_ct   = ct_counts[ct_counts >= MIN_CELLS].index
    adata     = adata[adata.obs[CELLTYPE_COL].isin(keep_ct)].copy()

    # rank genes (all genes)
    sc.tl.rank_genes_groups(
        adata,
        groupby=CELLTYPE_COL,
        key_added=f"{CELLTYPE_COL}_DE",
        use_raw=False,
        method="wilcoxon",
        n_genes=adata.n_vars,
    )

    write_matrices(adata, prefix)


def main():
    h5_files = sorted(DATA_DIR.glob("*.h5ad"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5ad files found in {DATA_DIR.resolve()}")

    for h5 in h5_files:
        try:
            process_file(h5)
        except Exception as e:
            print(f"⚠️  {h5.name}: {e}")


if __name__ == "__main__":
    main()
