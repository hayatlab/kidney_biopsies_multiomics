"""STEP 1: write the updated per-core metadata into xenium_labeled.h5ad obs.

Refreshes the existing `supp_*` columns in place (no doubling), adds the two new
fields (ck_number_rna/atac), and re-derives egfr_group / ifta_group / EGFR / IFTA_pct
from the refreshed values. A backup MUST be made before running (see the .sh / cp step).
"""
import sys
from pathlib import Path
import numpy as np
import scanpy as sc
import config as C

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "01_xenium_preprocessing_niches" / "nichecompass_replot_level2"))
sys.path.insert(0, str(HERE.parent / "01_xenium_preprocessing_niches"))  # for the downstream module
from metadata_common import apply_core_metadata, load_percore, CANON2OBS   # noqa: E402
import nichecompass_5regions_downstream as ds               # noqa: E402


def main():
    print("[load]", C.XENIUM_LABELED, flush=True)
    a = sc.read_h5ad(C.XENIUM_LABELED)
    percore = load_percore()

    before = {c: int(a.obs[c].notna().sum()) if c in a.obs else 0
              for c in set(CANON2OBS.values())}
    written = apply_core_metadata(a, percore)
    # re-derive binned groups + numeric EGFR / IFTA_pct from refreshed supp_ values
    ds.derive_ifta_egfr_groups(a)
    # numeric eGFR-slope convenience column
    a.obs["egfr_slope"] = (a.obs["supp_total_slope"].astype(str)
                           .str.replace(",", ".", regex=False)
                           .apply(lambda s: float(s) if s.replace("-", "").replace(".", "").isdigit() else np.nan))

    print("\ncoverage (cells non-NaN) BEFORE -> AFTER for a few key fields:")
    for c in ["supp_disease", "supp_age_at_biopsy", "supp_egfr_at_time_of_biopsy",
              "supp_total_slope", "supp_classification", "supp_Cold ischemia time (min)",
              "supp_DBD/DCD", "ck_number_rna", "ck_number_atac"]:
        aft = int(a.obs[c].notna().sum()) if c in a.obs else 0
        print(f"  {c:36s} {before.get(c,0):8d} -> {aft:8d}")

    # keep obs writeable: object columns -> str (NaN preserved as 'nan' -> restore)
    for col in a.obs.columns:
        if a.obs[col].dtype == object:
            a.obs[col] = a.obs[col].astype(str)
    print("\n[write]", C.XENIUM_LABELED, flush=True)
    a.write_h5ad(C.XENIUM_LABELED, compression="gzip")
    print("[done] wrote %d cells; refreshed/added %d obs columns" % (a.n_obs, len(written)), flush=True)


if __name__ == "__main__":
    main()
