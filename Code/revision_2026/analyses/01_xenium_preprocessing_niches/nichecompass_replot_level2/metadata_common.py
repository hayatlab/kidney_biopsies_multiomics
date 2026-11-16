"""Shared: apply the coalesced per-core metadata (from
celltype_concordance/15_update_metadata_slope.py -> per_sample_metadata_updated.csv)
onto an AnnData's obs, writing to the EXISTING canonical column names (no doubling).

Maps each canonical per-core field to the obs column it should land in: existing
`supp_*` columns are refreshed in place; only the two genuinely-new fields
(ck_number_rna/atac) create new columns.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
from pathlib import Path
import numpy as np
import pandas as pd

PERCORE_CSV = Path(f"{KB_ANALYSIS}/celltype_concordance/"
                   "outputs/per_sample_metadata_updated.csv")

# canonical field (per-core CSV) -> target obs column (existing supp_* refreshed; ck_* new)
CANON2OBS = {
    "disease": "supp_disease", "age": "supp_age_at_biopsy", "sex": "supp_sex",
    "ifta_pct": "supp_IFTA", "total_glomeruli": "supp_total_glomeruli",
    "global_gs": "supp_global_gs", "prop_gs": "supp_prop_gs",
    "date_of_biopsy": "supp_date_of_biopsy",
    "creatinine_umol_l": "supp_creatinine_at_time_of_biopsy",
    "treatment_before": "supp_treatment_before", "ace_inhibition": "supp_ace_inhibition",
    "angiotensin_ii_inhibition": "supp_angiotensin_ii_inhibition", "statin": "supp_statin",
    "corticosteroids": "supp_corticosteroids", "treatment_after": "supp_treatment_after",
    "dialysis": "supp_dialysis", "transplant": "supp_transplant",
    "date_dialysis_transplant": "supp_date_of_dialysis_transplant",
    "egfr_at_biopsy": "supp_egfr_at_time_of_biopsy",
    "proteinuria": "supp_proteinuria_around_biopsy",
    "total_slope": "supp_total_slope", "total_intercept": "supp_total_intercept",
    "immediate_slope": "supp_immediate_slope", "chronic_slope": "supp_chronic_slope",
    "classification": "supp_classification", "jm_slope": "supp_jm_slope",
    "jm_intercept": "supp_jm_intercept", "xenium_available": "supp_xenium_available",
    "cold_ischemia_min": "supp_Cold ischemia time (min)", "dbd_dcd": "supp_DBD/DCD",
    "ck_number_rna": "ck_number_rna", "ck_number_atac": "ck_number_atac",
    "ethnicity": "Etnicity",
}


def load_percore(path: Path = PERCORE_CSV) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def apply_core_metadata(adata, percore: pd.DataFrame | None = None,
                        tma_key: str = "TMA_ID", verbose: bool = True) -> list[str]:
    """Broadcast per-core metadata onto obs by TMA_ID, writing canonical obs columns.

    Existing values are kept only where the per-core table has no value for that core
    (so we never blank out a field). Returns the list of obs columns written."""
    if percore is None:
        percore = load_percore()
    tma = adata.obs[tma_key].astype(str)
    written = []
    for canon, obscol in CANON2OBS.items():
        if canon not in percore.columns:
            continue
        vals = tma.map(percore[canon])                     # core -> value, NaN if absent
        vals.index = adata.obs.index
        # Overwrite unconditionally (changed 2026-07-22). This used to keep the
        # existing obs value wherever the per-core table had none -- but those obs
        # values are the output of the earlier broken fuzzy join, so the fallback
        # quietly resurrected mis-assigned clinical data for exactly the cores whose
        # record is genuinely missing. Missing must stay missing.
        adata.obs[obscol] = vals.values
        written.append(obscol)
    if verbose:
        nz = tma.isin(percore.index).sum()
        print(f"[apply_core_metadata] {len(written)} obs columns refreshed/added; "
              f"{nz}/{adata.n_obs} cells have a matched core", flush=True)
    return written
