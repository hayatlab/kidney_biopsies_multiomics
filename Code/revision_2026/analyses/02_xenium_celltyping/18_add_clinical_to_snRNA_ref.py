"""Add donor-level clinical metadata to the snRNA reference (adata_20260717.h5ad).

Source of clinical fields : scVI_kidney_biopsy_with_clinical.h5ad (one consistent schema).
Join                      : reference obs['sample'] (K-numbers) <-> clinical obs['sample_ID'].
All 41 clinical fields are constant within a donor (verified) -> broadcast onto cells.

Coverage: 60/63 reference samples matched (96.3% of cells). Three samples have NO record in
this clinical object and are left NaN (not dropped, not imputed):
  * K68 - no clinical record in any available source.
  * K98 - no clinical record in any available source.
  * K95 - absent HERE, but present in 'Supplemental Data 1_metadata_update with slope.csv'
          (control_lowGFR, age 47, male, creatinine 257 umol/L, eGFR 26, donor of K95 and K98).
          NOT filled, to keep a single consistent source (h5ad-only, per user).

Writes a NEW file alongside the source (the original is NOT modified):
  adata_20260717_with_clinical.h5ad + adata_20260717_per_donor_clinical.csv + a .txt note.
"""
from __future__ import annotations
from kbpaths import KB_DATA  # repository path config (config/kbpaths.py)
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

DATA = Path(f"{KB_DATA}/charlotte_new_data")
REF = DATA / "adata_20260717.h5ad"
CLIN = Path(f"{KB_DATA}/scVI_kidney_biopsy_with_clinical.h5ad")
OUT = DATA / "adata_20260717_with_clinical.h5ad"
SIDE = DATA / "adata_20260717_per_donor_clinical.csv"
NOTE = DATA / "adata_20260717_with_clinical_README.txt"

REF_KEY = "sample"
CLIN_KEY = "sample_ID"

CLIN_FIELDS = [
    "disease", "disease_id", "Group", "Study",
    "age_at_biopsy", "age_group", "sex_clinical", "BMI",
    "systolic_blood_pressure", "diastolic_blood_pressure ", "hypertension_clinical",
    "diabetes", "HbA1c (%)", "heart_liver_disease",
    "egfr_at_time_of_biopsy", "egfr_group", "creatinine_at_time_of_biopsy",
    "IFTA", "ifta_group", "total_glomeruli", "global_gs", "prop_gs",
    "proteinuria_around_biopsy", "proteinuria_group",
    "total_slope", "total_intercept", "immediate_slope", "chronic_slope",
    "jm_slope", "jm_intercept", "classification",
    "treatment_before", "treatment_after", "ace_inhibition",
    "angiotensin_ii_inhibition", "statin", "corticosteroids_clinical",
    "dialysis", "transplant", "patient_ID", "ck_number_rna", "ck_number_atac",
]


def san(c: str) -> str:
    return (c.strip().replace(" ", "_").replace("(", "").replace(")", "")
            .replace("%", "pct").replace("/", "_"))


def main() -> int:
    # 1. per-donor clinical table (constant within donor -> first())
    co = ad.read_h5ad(str(CLIN), backed="r").obs
    fields = [f for f in CLIN_FIELDS if f in co.columns]
    print(f"[clinical] using {len(fields)} donor-level fields", flush=True)
    tab = co.groupby(co[CLIN_KEY].astype(str), observed=True)[fields].first()
    tab.columns = [san(c) for c in tab.columns]
    sfields = list(tab.columns)
    tab.to_csv(SIDE)
    print(f"[clinical] per-donor table {tab.shape} -> {SIDE.name}", flush=True)

    # 2. load reference (full, to write the new file)
    print(f"[ref] loading {REF} ...", flush=True)
    a = sc.read_h5ad(str(REF))
    smp = a.obs[REF_KEY].astype(str)
    matched = smp.isin(tab.index)
    print(f"[merge] cells matched: {matched.sum()}/{a.n_obs} ({matched.mean()*100:.1f}%)", flush=True)
    unmatched_samples = sorted(set(smp[~matched].unique()))
    print(f"[merge] reference samples with NO clinical record: {unmatched_samples}", flush=True)

    # 3. broadcast; keep numeric fields numeric (NaN preserved), else str ('nan' for missing)
    for f in sfields:
        vals = smp.map(tab[f])
        num = pd.to_numeric(vals, errors="coerce")
        if num.notna().sum() >= 0.5 * vals.notna().sum() and vals.notna().sum() > 0:
            a.obs[f] = num.to_numpy(dtype=float)          # numeric field
        else:
            a.obs[f] = vals.astype(object).where(vals.notna(), "nan").astype(str).to_numpy()
    # per-field NaN report
    print("\n[merge] % NaN per clinical field (after merge):", flush=True)
    for f in sfields:
        col = a.obs[f]
        n_nan = int(col.isna().sum()) if col.dtype.kind == "f" else int((col.astype(str) == "nan").sum())
        print(f"   {f:34s} {100*n_nan/a.n_obs:5.1f}% NaN", flush=True)

    # 4. write new file alongside (original untouched)
    print(f"\n[write] {OUT} (gzip) ...", flush=True)
    a.write_h5ad(str(OUT), compression="gzip")
    NOTE.write_text(
        "adata_20260717_with_clinical.h5ad\n=================================\n"
        f"= adata_20260717.h5ad + {len(sfields)} donor-level clinical fields from\n"
        "  scVI_kidney_biopsy_with_clinical.h5ad, joined on sample<->sample_ID.\n\n"
        f"Coverage: {matched.sum()}/{a.n_obs} cells ({matched.mean()*100:.1f}%); "
        f"{len(unmatched_samples)} samples with NO clinical record -> NaN: {unmatched_samples}.\n"
        "  K68, K98: no clinical record in any available source.\n"
        "  K95: absent here but present in 'Supplemental Data 1_metadata_update with slope.csv'\n"
        "       (control_lowGFR, age 47, male, creatinine 257 umol/L, eGFR 26, donor of K95 and K98);\n"
        "       left NaN to keep a single consistent source (h5ad-only).\n\n"
        "Per-donor table: adata_20260717_per_donor_clinical.csv\n"
        "Original adata_20260717.h5ad is unchanged.\n"
    )
    print(f"[done] wrote {OUT.name}, {SIDE.name}, {NOTE.name}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
