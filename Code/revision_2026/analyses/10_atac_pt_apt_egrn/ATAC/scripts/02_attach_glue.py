#!/usr/bin/env python3
"""Attach each GLUE co-embedding to the ArchR gene score matrix, separately.

Base    : gene_score_matrix.h5ad (ArchR GeneScoreMatrix, human hg38, HGNC symbols)
Attach  : glue/ and glue_v3/ glue_for_pseudotime.csv -> obsm['X_glue'] (30 dims)
Output  : gsm_with_glue.h5ad, gsm_with_glue_v3.h5ad
Seed    : 42 (unused; all operations deterministic)
"""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import anndata as ad, pandas as pd, numpy as np, scipy, h5py, sys, os

SEED = 42
GSM  = f"{KB_DATA}/New_ATAC_data/gene_score_matrix.h5ad"
A    = f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"
OUT  = f"{A}/objects"; os.makedirs(OUT, exist_ok=True)
LAT  = [str(i) for i in range(30)]

print("versions:", "anndata", ad.__version__, "| pandas", pd.__version__,
      "| numpy", np.__version__, "| scipy", scipy.__version__, flush=True)

print(f"\nreading {GSM} ...", flush=True)
gsm = ad.read_h5ad(GSM)
print(f"  gene score matrix: {gsm.shape[0]:,} nuclei x {gsm.shape[1]:,} genes", flush=True)
assert gsm.obs_names.is_unique, "gene score matrix obs_names not unique"

# The stored obsm DataFrames carry column names with "/" (level1 "VSM/P", level2 "VSMC/P").
# h5py treats "/" as a path separator, so current anndata refuses to write them. Rename
# rather than drop -- these are decoupler cell-type scores we want to keep.
renamed = {}
for k, v in gsm.obsm.items():
    if isinstance(v, pd.DataFrame):
        bad = [c for c in v.columns if "/" in str(c)]
        if bad:
            m = {c: str(c).replace("/", "_") for c in bad}
            gsm.obsm[k] = v.rename(columns=m)
            renamed[k] = m
if renamed:
    print("  sanitised obsm column names:", renamed, flush=True)
gsm.uns["obsm_columns_renamed"] = {k: list(m.items()) for k, m in renamed.items()}
print("  obsm carried:", list(gsm.obsm.keys()), flush=True)

for ver in ["glue", "glue_v3"]:
    print(f"\n══════ {ver} ══════", flush=True)
    df = pd.read_csv(f"{A}/{ver}/glue_for_pseudotime.csv", low_memory=False)
    df = df[~df.cellid.astype(str).str.split(r"[_#]").str[0].isin(["K95", "K96", "K98"])]  # AKI exclusion (aki_c5_exclusion_20261022)
    atac = df[df.modality == "ATAC"].copy()
    atac = atac.set_index("cellid")
    assert atac.index.is_unique, f"{ver}: duplicate ATAC cellids"

    missing = atac.index.difference(gsm.obs_names)
    assert len(missing) == 0, f"{ver}: {len(missing)} GLUE ATAC ids absent from gene score matrix"

    sub = gsm[atac.index].copy()          # subset AND order by the GLUE row order
    assert (sub.obs_names == atac.index).all(), "ordering mismatch after subset"

    sub.obsm["X_glue"] = atac[LAT].to_numpy(dtype=np.float32)
    # ATAC rows carry no level2 (that annotation exists only on the RNA side) and
    # RNA rows carry no cluster -- report missingness explicitly, never fill silently.
    for c in ["coembedding_id", "level1", "level2", "cluster"]:
        col = atac[c]
        n_miss = int(col.isna().sum())
        print(f"    obs glue_{c:<15} missing {n_miss:>7,} / {len(col):,} "
              f"({100*n_miss/len(col):.2f}%)", flush=True)
        sub.obs[f"glue_{c}"] = pd.Categorical(
            col.astype(object).where(col.notna(), "not stated").astype(str).values)
        sub.uns[f"glue_{c}_pct_missing"] = round(100 * n_miss / len(col), 4)
    sub.uns["glue_source"] = f"{A}/{ver}/glue_for_pseudotime.csv"
    sub.uns["glue_version"] = ver
    sub.uns["glue_n_latent"] = len(LAT)
    sub.uns["glue_rna_cells_in_coembedding"] = int((df.modality == "RNA").sum())
    sub.uns["gene_score_source"] = GSM
    sub.uns["genome_build"] = "hg38"
    sub.uns["gene_id_format"] = "HGNC symbol"
    sub.uns["seed"] = SEED

    nan = np.isnan(sub.obsm["X_glue"]).sum()
    print(f"  attached {sub.shape[0]:,} nuclei x 30 latent dims | NaN in X_glue: {nan}", flush=True)
    print(f"  dropped from gene score matrix (no embedding): {gsm.shape[0]-sub.shape[0]:,}", flush=True)
    print("  cluster composition:\n" +
          sub.obs["glue_cluster"].value_counts(dropna=False).to_string(), flush=True)

    p = f"{OUT}/gsm_with_{ver}.h5ad"
    sub.write_h5ad(p, compression="gzip")
    print(f"  wrote {p}  ({os.path.getsize(p)/1e9:.2f} GB)", flush=True)

print("\ndone", flush=True)
