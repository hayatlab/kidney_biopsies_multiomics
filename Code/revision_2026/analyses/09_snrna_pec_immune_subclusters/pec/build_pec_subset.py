#!/usr/bin/env python
"""Phase 1: subset the newly annotated single-cell object (adata_20260717) to PEC
(level2=='PEC', 7,004 cells) and attach old-annotation provenance so the trajectory
notebook and the colleague-facing comparison figure can be built on a small object.

- New object is 7.1 GB -> read backed='r' and subset (light on RAM).
- Provenance join old<->new by sample+barcode (old index 'K89_BC-1-0' -> norm 'K89_BC-1';
  new cellid 'K1_BC-1'). KB## and K## are DIFFERENT samples (verified), no rename applied.
- old_main_class from the OLD scVI object (all 327,921 cells) -> provenance of gained cells.
- old_pec_subcluster (pec_subcluster_0.80) from pec_annotated.h5ad (the 4,707 old PECs)
  -> confusion map on the retained overlap.

Outputs: pec_resubcluster_20260717/pec_new_20260717_raw.h5ad
Seed 42. Species human (HGNC symbols). No imputation.
"""
from kbpaths import KB_ANALYSIS, KB_DATA, KB_PROJECT  # repository path config (config/kbpaths.py)
import re, numpy as np, pandas as pd, scanpy as sc, anndata as ad
import h5py

NEW   = f"{KB_DATA}/charlotte_new_data/adata_20260717_with_clinical.h5ad"
OLD   = f"{KB_DATA}/scVI_kidney_biopsy_with_clinical.h5ad"
PECAN = f"{KB_PROJECT}/trajectory_pec/pec_annotated.h5ad"
OUT   = f"{KB_ANALYSIS}/pec_resubcluster_20260717/pec_new_20260717_raw.h5ad"
np.random.seed(42)

def norm_old(s):  # 'K89_AAAC...-1-0' -> 'K89_AAAC...-1'
    return re.sub(r"-\d+$", "", s)

def readcol(f, key):
    d = f[f"obs/{key}"]
    if isinstance(d, h5py.Group):  # categorical
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in d["categories"][:]])
        codes = d["codes"][:]
        return np.where(codes >= 0, cats[codes.clip(0)], "NaN").astype(object)
    v = d[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v], dtype=object)

# --- old-annotation lookups (barcode key -> label) ------------------------------
print("[1/4] reading OLD scVI obs (main_class) ...", flush=True)
with h5py.File(OLD, "r") as f:
    old_key = np.array([norm_old(s) for s in readcol(f, "index")], dtype=object)
    old_mc  = readcol(f, "main_class")
old_main = dict(zip(old_key, old_mc))

print("[2/4] reading pec_annotated obs (pec_subcluster_0.80) ...", flush=True)
with h5py.File(PECAN, "r") as f:
    pk = np.array([norm_old(s) for s in readcol(f, f['obs'].attrs.get('_index','index'))], dtype=object)
    psub = readcol(f, "pec_subcluster_0.80")
old_sub = dict(zip(pk, psub))

# --- subset NEW to PEC (backed) -------------------------------------------------
print("[3/4] backed-reading NEW object and subsetting level2==PEC ...", flush=True)
adata = sc.read_h5ad(NEW, backed="r")
mask  = (adata.obs["level2"].astype(str).values == "PEC")
print(f"    level2==PEC cells: {int(mask.sum())}", flush=True)
pec = adata[mask].to_memory()
adata.file.close()

# --- attach provenance ----------------------------------------------------------
print("[4/4] attaching provenance ...", flush=True)
keys = np.asarray(pec.obs_names).astype(str)  # cellid is the obs index
pec.obs["old_main_class"]    = [old_main.get(k, "not_in_old") for k in keys]
pec.obs["old_pec_subcluster"] = [old_sub.get(k, "not_old_PEC") for k in keys]

def provenance(mc):
    if mc == "13_PEC":       return "retained_PEC"
    if mc == "not_in_old":   return "new_cell"
    if mc == "11_ATL_DTL":   return "gained_ATL_DTL"
    if mc == "1_TAL":        return "gained_TAL"
    if mc in ("4_aPT_aTAL_aDCT", "0_PT", "3_dPT_dTAL_dDCT"): return "gained_PT_axis"
    return "gained_other"
pec.obs["provenance"] = [provenance(m) for m in pec.obs["old_main_class"]]

print("\nprovenance breakdown:")
print(pec.obs["provenance"].value_counts().to_string())
print("\nold_main_class breakdown:")
print(pec.obs["old_main_class"].value_counts().head(15).to_string())

# X sanity (expected: library-normalized, not log, no counts layer)
xmax = float(pec.X.max()); xint = bool(np.allclose(pec.X.data[:5000], np.round(pec.X.data[:5000])))
print(f"\nX: shape={pec.shape} max={xmax:.2f} looks_integer={xint} layers={list(pec.layers.keys())} obsm={list(pec.obsm.keys())}")

pec.write_h5ad(OUT)
print(f"\n[done] wrote {OUT}  ({pec.n_obs} cells x {pec.n_vars} genes)")
