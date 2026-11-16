#!/usr/bin/env python
"""Fig. 3h — local cell-cell communication from aPT (sender) to immune cells (target) along Factor 1: LOESS
(frac 0.4) of the per-root interaction score of each ligand-receptor pair against the Factor 1 score (ULM) of the
root cell, x cropped to [-7.25, 5.25] (the crop of Fig. 3g); the 12 strongest pairs are labelled.
Panel        : Main Figures.pdf, Figure 3, panel h
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/local_cci_results/cci_score_df.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/09_cci_postprocess.py, full run without --reuse,
                  from the chunks of 08_run_local_cci.py / 08_run_local_cci_array.sh; Factor 1 per root from
                  01_score_factors_xenium.py)
Writes       : $KB_FIG_OUT/Fig3/Fig3h_cci_aPT_to_Immune.pdf (+ .png, _loess.csv, Fig3h_cci_meta.json)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/09_cci_postprocess.py::plot (imported; the `outdir`,
               `pairs` and `stem` keywords were added by patch rule build/patch_rules/fig_xenium.py, defaults
               unchanged; only the aPT -> Immune pair of the three is drawn)
               (original run: 09_cci_postprocess.sh --reuse, SLURM job 329451, 2026-09-07 -> fig3i_aPT_to_Immune)
Illustrator  : "Immune" written "immune" in the title; "Factor1_ulm" written "factor 1_ulm"; pair labels recoloured
               and moved, leader lines removed
Environment  : conda env cci_env
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/cci/fig3i_aPT_to_Immune.pdf
               (the _loess.csv is byte-identical to fig3i_aPT_to_Immune_loess.csv)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

cp = load_script("03_xenium_factor1_neighbourhoods_cci/09_cci_postprocess.py")


def main():
    np.random.seed(cp.SEED)                       # global state as in 09_cci_postprocess.py::main
    sns.set_style("white")
    df = pd.read_parquet(cp.CCI / "cci_score_df.parquet")
    print(f"[reuse] cci_score_df.parquet {df.shape}", flush=True)
    out = Path(fig_out("Fig3", "Fig3h_cci_aPT_to_Immune.pdf"))
    cp.plot(df, argparse.Namespace(top_n=12), outdir=out.parent, pairs=[("aPT", "Immune")],
            stem="Fig3h_cci_{a}_to_{b}")
    rename(out.parent, {"09_cci_meta.json": "Fig3h_cci_meta.json"})
    print(f"[Fig3h] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
