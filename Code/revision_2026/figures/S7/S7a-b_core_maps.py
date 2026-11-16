#!/usr/bin/env python
"""Supp. Fig. 7a-b — cell-boundary maps of two Xenium cores: cell type (level2_refined), Factor 1 score (ULM) per cell
and NicheCompass niche (latent_leiden_0.4): a, the TIN core with low eGFR of ED 4d (K18, eGFR 23); b, a control core
with preserved kidney function (eGFR 60).
Panel        : All Supp Figures.pdf, Supplementary Fig. 7, panels a and b
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/geodataframes/cell_boundaries.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/02_build_geodataframes.py)
               $KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet (01_score_factors_xenium.py)
               $KB_ANALYSIS/vivien_takeover/outputs/authoritative_core_egfr.csv (10_build_authoritative_egfr.py)
Writes       : $KB_FIG_OUT/S7/S7a_core_maps_XEN116_K18.pdf, $KB_FIG_OUT/S7/S7b_core_maps_XEN008_control.pdf (+ .png)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::suppl7 (imported; the `names` and
               `prefix` keywords were added by patch rule build/patch_rules/fig_xenium.py, defaults unchanged); one
               call per core, cores selected by public id (XEN116 = K18 TIN; XEN008 = control) instead of the
               SUPPL7_CORES list (original run: 05_xenium_figures.sh, SLURM job 329443, 2026-09-07, which drew only
               these two panels)
Illustrator  : maps mirrored vertically; titles "<core> | <condition> | eGFR <x>" replaced by "TIN (low eGFR)" /
               "control (high eGFR)" (core names cropped); panel titles turned into legend titles ("cell type",
               "factor 1_ulm", "Niche"); the niche map of the TIN core omitted from panel a; dPT renamed aPT2
Environment  : conda env spatialdata
Verified     : 2026-11-12, DIFFERENT vs $KB_ANALYSIS/vivien_takeover/figures/xenium/suppl7_<K18 TIN core label>.pdf (0.16% of
               pixels) and suppl7_Appollo_d0_B5_CTRL.pdf (0.15%): only the suptitle band differs (core shown by its public id);
               rendered in memory with the original titles, both pages are PIXEL-IDENTICAL (PDF and PNG, text equal)

"high eGFR" for the control core is 60 ml/min/1.73 m2 (corrected per-core eGFR).
"""
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import core_display, load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402
from cores import core_label  # noqa: E402

PANELS = [("S7a", "XEN116"), ("S7b", "XEN008")]   # K18 TIN; control core
xf = load_script("03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py")


def main():
    np.random.seed(xf.SEED)                       # global state as in 05_xenium_figures.py::main
    warnings.filterwarnings("ignore", message=".*geographic CRS.*")
    sns.set_style("white")
    # ---- as in 05_xenium_figures.py::main (panel suppl7) ----
    obs = pd.read_parquet(xf.OUT / "xenium_factor_scores.parquet")
    gdf = gpd.read_parquet(xf.OUT / "geodataframes" / "cell_boundaries.parquet")
    gdf.index = gdf["sample_cell_id"]
    gdf[xf.FACTOR] = obs.set_index("sample_cell_id")[xf.FACTOR].reindex(gdf.index)
    eg, _, _, _ = xf.corrected_core_egfr(obs)        # NOT the corrupted EGFR column
    outdir = Path(fig_out("S7", "S7a_core_maps.pdf")).parent
    present = set(gdf["TMA_ID"])
    for panel, xen in PANELS:
        core = core_label(xen)
        if core not in present:
            raise SystemExit(f"Supp 7 core {xen} is not in the object")
        xf.suppl7(gdf, outdir, [core], eg, names={core: core_display(xen)}, prefix=f"{panel}_core_maps_")
    print(f"[S7a-b] wrote {outdir}/S7a_core_maps_*.pdf, S7b_core_maps_*.pdf", flush=True)


if __name__ == "__main__":
    main()
