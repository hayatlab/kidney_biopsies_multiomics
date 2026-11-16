#!/usr/bin/env python
"""Fig. 3f — distance-bounded Voronoi graph of one TIN core (K16, eGFR 5) coloured by cell type, and one 2-hop
neighbourhood around an aPT root cell drawn by cell type and by Factor 1 score (ULM).
Panel        : Main Figures.pdf, Figure 3, panel f
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/geodataframes/cell_boundaries.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/02_build_geodataframes.py)
               $KB_ANALYSIS/vivien_takeover/outputs/neighborhoods/voronoi_graphs_10.pkl
                 (analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py)
               $KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet (01_score_factors_xenium.py)
Writes       : $KB_FIG_OUT/Fig3/Fig3f_voronoi_2hop_neighbourhood.pdf (+ .png)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::fig3g (imported; the `stem`,
               `core_name` and `factor_cbar_kwds` keywords were added by patch rule build/patch_rules/fig_xenium.py,
               defaults unchanged; the core is selected by its internal label and shown by its public id XEN085;
               factor_cbar_kwds={} draws the full-height Factor 1 colour bar of the published panel - the shrunk colour
               bar of the shipped default was added after the published panel was drawn)
               (original run: 05_xenium_figures.sh, SLURM job 329429, 2026-09-07 00:37, which drew only this panel;
               core = FIG3G_CORE)
Illustrator  : core name removed from the title; the two neighbourhood views moved next to the core map and
               shrunk, their legend replaced by cell-type names written next to the cells; "2-hop neighbourhood
               (19 cells)" shortened; "Factor1_ulm" written "factor 1_ulm"; a red box and a scale bar added to the
               core map; dPT renamed aPT2
Environment  : conda env spatialdata
Verified     : 2026-11-12, DIFFERENT (0.10% of pixels) vs $KB_ANALYSIS/vivien_takeover/figures/xenium/fig3g_voronoi_2hop.pdf:
               93% of the differing pixels are the title (core shown as XEN085 K16 instead of its internal label); the rest
               (0.008%) are anti-aliasing specks where neighbourhood edges cross - fig3g draws the subgraph edges in set order,
               which follows Python's string-hash seed, so two runs of this script differ by ~0.005% (NEAR). Rendered in memory
               with the original title, the page has the original size and text and differs by 85-121 pixels (0.004-0.005%)
"""
import pickle
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

CORE = "XEN085"   # K16, TIN
xf = load_script("03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py")


def main():
    np.random.seed(xf.SEED)                       # global state as in 05_xenium_figures.py::main
    warnings.filterwarnings("ignore", message=".*geographic CRS.*")
    sns.set_style("white")
    # ---- as in 05_xenium_figures.py::main ----
    obs = pd.read_parquet(xf.OUT / "xenium_factor_scores.parquet", columns=["sample_cell_id", xf.FACTOR])
    gdf = gpd.read_parquet(xf.OUT / "geodataframes" / "cell_boundaries.parquet")
    gdf.index = gdf["sample_cell_id"]
    gdf[xf.FACTOR] = obs.set_index("sample_cell_id")[xf.FACTOR].reindex(gdf.index)
    with open(xf.OUT / "neighborhoods" / "voronoi_graphs_10.pkl", "rb") as f:
        graphs = pickle.load(f)
    out = Path(fig_out("Fig3", "Fig3f_voronoi_2hop_neighbourhood.pdf"))
    res = xf.fig3g(gdf, graphs, out.parent, core=core_label(CORE), stem=out.stem, core_name=core_display(CORE),
                   factor_cbar_kwds={})          # colour bar as published (full height)
    print({k: res[k] for k in ("n_cells", "n_edges", "inset_size", "inset_root_type", "inset_composition")})
    print(f"[Fig3f] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
