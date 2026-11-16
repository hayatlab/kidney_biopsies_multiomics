#!/usr/bin/env python
"""Fig. 3c — mean Factor 1 score (ULM) of all cells per Xenium core vs eGFR; each core is drawn as its spatial
thumbnail (cells coloured by Factor 1, frame coloured by condition); OLS fit with its 95% confidence band.
Panel        : Main Figures.pdf, Figure 3, panel c
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/xenium_factor_scores.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/01_score_factors_xenium.py)
               $KB_ANALYSIS/vivien_takeover/outputs/authoritative_core_egfr.csv
                 (analyses/03_xenium_factor1_neighbourhoods_cci/10_build_authoritative_egfr.py)
Writes       : $KB_FIG_OUT/Fig3/Fig3c_factor1_per_core_vs_egfr.pdf (+ .png)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::fig3d (imported; the `stem` and
               `write_table` keywords were added by patch rule build/patch_rules/fig_xenium.py, defaults unchanged;
               write_table=False skips fig3d_per_core.csv, whose rows are internal core labels with donor numbers)
               (original run: 05_xenium_figures.sh, SLURM job 329425, 2026-09-07)
Illustrator  : "Factor1_ulm" written "factor 1_ulm" and "per core" relabelled "per sample" in the y label; legend
               title "core outline" removed; statistics text and legend moved; canvas rescaled
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/fig3d_factor1_vs_egfr.pdf

Each point is a core (TMA punch), not a donor: several cores can come from one donor. fig3d also prints a donor-level
sensitivity analysis (original log: 111 donors, r = -0.400, p = 1.4e-05) that is not drawn in the panel.
"""
import json
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
    obs = pd.read_parquet(xf.OUT / "xenium_factor_scores.parquet")
    out = Path(fig_out("Fig3", "Fig3c_factor1_per_core_vs_egfr.pdf"))
    res = xf.fig3d(obs, out.parent, stem=out.stem, write_table=False)
    print(json.dumps(res, indent=1, default=float))
    print(f"[Fig3c] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
