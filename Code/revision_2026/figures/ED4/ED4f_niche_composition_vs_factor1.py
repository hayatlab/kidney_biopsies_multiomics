#!/usr/bin/env python
"""ED 4f — NicheCompass niche composition of the 2-hop Voronoi neighbourhood of every Xenium cell against the Factor 1
score (ULM) of its root cell (stacked spline ribbon as in Fig. 3g; niches = latent_leiden_0.4, 10 niches; x cropped
to [-7.25, 5.25]).
Panel        : Extended Data Figures.pdf, Extended Data Fig. 4, panel f
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/neighborhoods/nh_2hop_niche_compo_10.csv
                 (analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py)
               $KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet (01_score_factors_xenium.py)
Writes       : $KB_FIG_OUT/ED4/ED4f_niche_composition_vs_factor1.pdf (+ .png, _fitted.csv)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::ribbon (imported unchanged)
               (original run: 05_xenium_figures.sh, SLURM job 329425, 2026-09-07 -> ed3f_niche_vs_factor1)
Illustrator  : title "vs Factor 1" written "vs factor 1"; x label "Factor1_ulm (cropped to [-7.2, 5.2])" replaced by
               "Factor 1 score (ULM)"; y label "neighbourhood" -> "neighbourhoods"; niches 1 and 4 labelled on the
               ribbons
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/ed3f_niche_vs_factor1.pdf
               (the _fitted.csv is byte-identical to ed3f_niche_vs_factor1_fitted.csv)
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
    # ---- as in 05_xenium_figures.py::main (panel ed3f) ----
    obs = pd.read_parquet(xf.OUT / "xenium_factor_scores.parquet", columns=["sample_cell_id", xf.FACTOR])
    nz = pd.read_csv(xf.OUT / "neighborhoods" / "nh_2hop_niche_compo_10.csv", index_col=0)
    nlv = [c for c in nz.columns if c not in ("TMA_ID", "root_niche", "disease")]
    nz[xf.FACTOR] = obs.set_index("sample_cell_id")[xf.FACTOR].reindex(nz.index)
    out = Path(fig_out("ED4", "ED4f_niche_composition_vs_factor1.pdf"))
    res = xf.ribbon(nz, nlv, "Niche", xf.FACTOR, "Niche composition of 2-hop neighbourhoods vs Factor 1",
                    "Proportions in 2-hop neighbourhood", out.parent, out.stem, cmap_name="tab10")
    print(res)
    print(f"[ED4f] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
