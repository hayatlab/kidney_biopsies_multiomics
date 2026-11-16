#!/usr/bin/env python
"""Fig. 3g — cell-type composition of the 2-hop Voronoi neighbourhood of every Xenium cell against the Factor 1 score
(ULM) of its root cell: per cell type a cubic least-squares spline (8 interior knots) of the proportions, clipped to
[0, 1], renormalised and stacked; x cropped to [-7.25, 5.25] (contiguous 0.25-unit bins with >= 500 neighbourhoods,
99.07% of the 258,375 neighbourhoods; analyses/03_xenium_factor1_neighbourhoods_cci/13_crop_justification.py).
Panel        : Main Figures.pdf, Figure 3, panel g
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/neighborhoods/nh_2hop_compo_10.csv
                 (analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py)
               $KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet (01_score_factors_xenium.py)
Writes       : $KB_FIG_OUT/Fig3/Fig3g_celltype_composition_vs_factor1.pdf (+ .png, _fitted.csv)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::ribbon (imported unchanged)
               (original run: 05_xenium_figures.sh, SLURM job 329425, 2026-09-07 -> fig3h_celltype_vs_factor1)
Illustrator  : title removed; "Factor1_ulm" written "factor 1_ulm"; y label "neighbourhood" -> "neighbourhoods";
               PT / Immune / FIB / aPT written on the ribbons; legend title removed; dPT renamed aPT2
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/fig3h_celltype_vs_factor1.pdf
               (the _fitted.csv is byte-identical to fig3h_celltype_vs_factor1_fitted.csv)
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402

xf = load_script("03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py")


def main():
    np.random.seed(xf.SEED)                       # global state as in 05_xenium_figures.py::main
    warnings.filterwarnings("ignore", message=".*geographic CRS.*")
    sns.set_style("white")
    # ---- as in 05_xenium_figures.py::main (panel fig3h) ----
    obs = pd.read_parquet(xf.OUT / "xenium_factor_scores.parquet", columns=["sample_cell_id", xf.FACTOR])
    compo = pd.read_csv(xf.OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0)
    levels = [c for c in compo.columns if c not in ("TMA_ID", "root_cell_type", "disease")]
    compo[xf.FACTOR] = obs.set_index("sample_cell_id")[xf.FACTOR].reindex(compo.index)
    out = Path(fig_out("Fig3", "Fig3g_celltype_composition_vs_factor1.pdf"))
    res = xf.ribbon(compo, levels, "Cell type", xf.FACTOR,
                    "Cell-type composition of 2-hop neighbourhoods vs Factor 1",
                    "Proportions in 2-hop neighbourhood", out.parent, out.stem)
    print(res)
    print(f"[Fig3g] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
