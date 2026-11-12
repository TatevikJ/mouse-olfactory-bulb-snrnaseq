"""
Merge per-sample Parse Bio outputs into one h5ad with raw integer UMI counts.

  1. Read each per-sample triplet from INPUT_DIR/<raw_id>/:
     count_matrix.mtx (cells x genes, int UMI), all_genes.csv, cell_metadata.csv.
  2. Remap raw mouse IDs to canonical names (ASO11-20 -> ASO1-6, WT11-17 -> WT1-6)
     and set obs columns: sample_raw, sample, orig.ident, condition (WT/ASO).
  3. Use `bc_wells` as obs_names (already globally unique via `__s<N>` suffix)
     and `gene_name` as var_names (deduped via var_names_make_unique).
  4. Concat all 12 samples on the cell axis with outer join on genes
     (fill_value=0) and merge="same" for var columns.
  5. Add hSNCA_count and hSNCA_detected obs columns from `hSNCA_transgene`.
  6. Keep genes with count >0 in >= 3 cells.
  7. Compute `nCount_RNA` (X.sum per cell) and print a per-sample summary.
  8. Write gzip-compressed (level 7) to OUTPUT_FILE.
"""

# %% Imports
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.utils import make_index_unique
from scipy.io import mmread

# %% Config
ANALYSIS_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = ANALYSIS_ROOT / "GRCm39_v113_hSNCA_filtered_matrices"
OUTPUT_DIR = ANALYSIS_ROOT / "outputs" / "h5ad"
OUTPUT_FILE = OUTPUT_DIR / "aso_wt_data.h5ad"

# Raw mouse ID -> canonical name
SAMPLE_MAP = {
    "ASO11": "ASO1", "ASO13": "ASO2", "ASO15": "ASO3",
    "ASO16": "ASO4", "ASO17": "ASO5", "ASO20": "ASO6",
    "WT11":  "WT1",  "WT12":  "WT2",  "WT13":  "WT3",
    "WT14":  "WT4",  "WT15":  "WT5",  "WT17":  "WT6",
}

# %% Read each per-sample triplet
adatas = {}
for raw_id, canonical in SAMPLE_MAP.items():
    sdir = INPUT_DIR / raw_id
    print(f"[{canonical:<5}] reading {sdir.name}")

    # Parse Bio writes count_matrix as cells x genes; keep as int32 sparse CSR
    X = mmread(str(sdir / "count_matrix.mtx")).tocsr().astype(np.int32)

    var = pd.read_csv(sdir / "all_genes.csv")
    obs = pd.read_csv(sdir / "cell_metadata.csv")

    # Shape sanity: cells from cell_metadata, genes from all_genes
    assert X.shape == (obs.shape[0], var.shape[0]), (
        f"shape mismatch in {raw_id}: "
        f"X={X.shape} obs={obs.shape[0]} var={var.shape[0]}"
    )

    # Use bc_wells as obs index, gene_name (deduped) as var index
    obs.index = pd.Index(obs["bc_wells"].astype(str).values)
    var["gene_name"] = var["gene_name"].astype(str).str.strip()
    var.index = make_index_unique(pd.Index(var["gene_name"].values))

    a = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    # Sample-level metadata (overwrite incoming `sample` so it matches canonical)
    a.obs["sample_raw"] = raw_id
    a.obs["sample"] = canonical
    a.obs["orig.ident"] = canonical
    a.obs["condition"] = "WT" if canonical.startswith("WT") else "ASO"

    print(
        f"          shape: {a.shape}  X.max: {int(a.X.max())}  dtype: {a.X.dtype}"
    )
    adatas[canonical] = a

# %% Concatenate cells; align var by gene name
combined = ad.concat(
    adatas,
    axis=0,
    join="outer",       # keep all genes even if a sample is missing one
    merge="same",       # keep var columns only when identical across samples
    fill_value=0,
)
combined.obs_names_make_unique()

print(
    f"\nCombined: {combined.shape[0]:,} cells x {combined.shape[1]:,} genes"
)
print(f"X dtype: {combined.X.dtype}  X.max: {int(combined.X.max())}")
print("\nCells per sample:")
print(combined.obs["sample"].value_counts().sort_index())
print("\nCells per condition:")
print(combined.obs["condition"].value_counts())

# %% hSNCA transgene bookkeeping
gene_name = "hSNCA_transgene"
if gene_name in combined.var_names:
    col = combined.X[:, combined.var_names.get_loc(gene_name)]
    arr = col.toarray().flatten() if sp.issparse(col) else np.asarray(col).flatten()
    combined.obs["hSNCA_count"] = arr
    combined.obs["hSNCA_detected"] = (arr > 0).astype(int)

    print(f"\nhSNCA_transgene total UMI: {int(arr.sum()):,}")
    print("Detection rate (>=1 UMI) by condition:")
    print(combined.obs.groupby("condition", observed=True)["hSNCA_detected"].mean())
    print("Mean UMI per cell by condition:")
    print(combined.obs.groupby("condition", observed=True)["hSNCA_count"].mean())
else:
    print(f"\nWARNING: '{gene_name}' not found in var_names.")
    print(f"First 5 var_names: {list(combined.var_names[:5])}")

# %% Filter genes: keep those with count >0 in >= 3 cells
keep_genes = combined.X.getnnz(axis=0) >= 3
print(f"keeping {int(keep_genes.sum()):,} / {keep_genes.size:,} genes (count>0 in >=3 cells)")
combined = combined[:, keep_genes].copy()

# %% Per-sample nCount_RNA summary (X.sum per cell)
combined.obs["nCount_RNA"] = np.asarray(combined.X.sum(axis=1)).flatten()
print("\nnCount_RNA summary per sample:")
print(
    combined.obs.groupby("sample", observed=True)["nCount_RNA"]
    .describe()[["count", "mean", "std", "min", "25%", "50%", "75%", "max"]]
    .round(0)
    .astype(int)
)

# %% Write merged h5ad
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
combined.write_h5ad(OUTPUT_FILE, compression="gzip", compression_opts=7)
print(f"\nSaved: {OUTPUT_FILE}")
print(f"File size: {OUTPUT_FILE.stat().st_size / 1e9:.2f} GB")
