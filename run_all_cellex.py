#!/usr/bin/env python
"""
Run CELLEX metrics (DET, NSI, GES, EP, ESμ) on all *.h5ad datasets
in DATA_DIR and write results to CELLEX_OUT.

•  No absolute paths — everything is relative to the repo.
•  Works for both 'all‑gene' and '*_pc.h5ad' protein‑coding files.
•  Ready for GitHub: you only commit the script; users configure paths.

Requires: cellex, scanpy
"""

from pathlib import Path
import pprint
import scanpy as sc
import cellex


# ─────────────────────────────────────────────────────────────
# 1 | Configure relative directories
#    (edit these two lines only)
# ─────────────────────────────────────────────────────────────
DATA_DIR   = Path("data")          # folder containing your .h5ad files
CELLEX_OUT = Path("cellex_out")    # where CSV results will be stored
# ─────────────────────────────────────────────────────────────

CELLEX_OUT.mkdir(parents=True, exist_ok=True)


def run_cellex(h5_path: Path) -> None:
    """Compute CELLEX metrics and save CSVs."""
    adata = sc.read_h5ad(h5_path)
    es = cellex.ESObject(
        data=adata.to_df().T,                     # genes × cells
        annotation=adata.obs["cell_type"],
        normalize=False,                          # counts already normalized
        verbose=True,
    )
    es.compute(verbose=True)

    prefix = h5_path.stem                        # filename w/o .h5ad
    es.save_as_csv(keys=["all"],
                   path=CELLEX_OUT,
                   file_prefix=prefix,
                   verbose=True)
    print(f"✔  Saved CELLEX results → {prefix}__*.csv")


def main() -> None:
    h5_files = sorted(DATA_DIR.glob("*.h5ad"))
    if not h5_files:
        raise FileNotFoundError(f"No .h5ad files found in {DATA_DIR.resolve()}")

    print("Running CELLEX on:")
    pprint.pp([p.name for p in h5_files], compact=True)
    print()

    for fp in h5_files:
        try:
            run_cellex(fp)
        except Exception as err:
            print(f"⚠️  {fp.name} failed: {err}")

    print("\nAll done — outputs are in:", CELLEX_OUT.resolve())


if __name__ == "__main__":
    main()
