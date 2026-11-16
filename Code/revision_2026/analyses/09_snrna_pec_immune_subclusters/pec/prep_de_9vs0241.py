#!/usr/bin/env python
"""Build the MAST DE input for contrast: leiden_0.8 cluster 9 vs {0,2,4,11} (new PEC data).

Mirrors de_analysis_pec.R's expectation: X = RNA counts (as.Seurat counts='X' then
LogNormalize). Source counts = pec_new_20260717_raw.h5ad (X rowsum == nCount_RNA exactly).
leiden_0.8 is joined from the sweep object by obs_names (100% shared). Subset to the
C_0to9_5cl population {0,2,4,9,11} so the R "Other" group is exactly {0,2,4,11}.
Donor latent variable = patient_ID (60 donors; the new data has no 'donor_id' column).
seed 42. Output: de_9vs0241/pec_de_input_leiden08_C5.h5ad
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os
import numpy as np, scanpy as sc

DIR = f"{KB_ANALYSIS}/pec_resubcluster_20260717"
OUT = f"{DIR}/de_9vs0241"; os.makedirs(OUT, exist_ok=True)
KEEP = ["0", "2", "4", "9", "11"]

raw = sc.read_h5ad(f"{DIR}/pec_new_20260717_raw.h5ad")
sweep = sc.read_h5ad(f"{DIR}/new_pec_annotation/pec_leiden_sweep.h5ad", backed="r")

# sanity: rowsum(X) == nCount_RNA -> confirms X is RNA counts
rs = np.asarray(raw.X.sum(1)).ravel()
assert np.corrcoef(rs, raw.obs["nCount_RNA"].to_numpy())[0, 1] > 0.999, "X is not RNA counts"

# join leiden_0.8 from the sweep (align by obs_names)
lab = sweep.obs["leiden_0.8"].astype(str).reindex(raw.obs_names)
assert lab.notna().all(), "unmatched cells when joining leiden_0.8"
raw.obs["leiden_0.8"] = lab.values

assert "patient_ID" in raw.obs.columns, "patient_ID missing"
raw.obs["patient_ID"] = raw.obs["patient_ID"].astype(str)

sub = raw[raw.obs["leiden_0.8"].isin(KEEP)].copy()
sub.obs["leiden_0.8"] = sub.obs["leiden_0.8"].astype("category")

# report composition + donor coverage (no imputation; nothing dropped silently)
print("subset:", sub.n_obs, "cells,", sub.n_vars, "genes")
print("per-cluster n:\n", sub.obs["leiden_0.8"].value_counts().sort_index().to_string())
tgt = sub.obs["leiden_0.8"] == "9"
print(f"Target (cl9) = {int(tgt.sum())}  |  Other (0,2,4,11) = {int((~tgt).sum())}")
print("donors (patient_ID): total", sub.obs["patient_ID"].nunique(),
      "| in Target", sub.obs.loc[tgt, "patient_ID"].nunique(),
      "| in Other", sub.obs.loc[~tgt, "patient_ID"].nunique())
# how many genes are all-zero in the subset (MAST will drop these anyway)
nz = np.asarray((sub.X > 0).sum(0)).ravel()
print(f"all-zero genes in subset: {(nz == 0).sum()}/{sub.n_vars} "
      f"(kept in file; FindMarkers min.pct=0 handles them)")

# keep obs lean + string-typed for zellkonverter/Seurat
sub.obs = sub.obs[["leiden_0.8", "patient_ID", "sample", "old_pec_subcluster",
                   "nCount_RNA", "nFeature_RNA"]].copy()
for c in ["sample", "old_pec_subcluster"]:
    sub.obs[c] = sub.obs[c].astype(str)
sub.var = sub.var[[]].copy()
del sub.obsm, sub.uns
outf = f"{OUT}/pec_de_input_leiden08_C5.h5ad"
sub.write_h5ad(outf)
print("wrote", outf)
