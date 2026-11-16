#!/usr/bin/env python
"""ED 4c — cell-type composition of 2-hop Voronoi neighbourhoods: mean proportion of each cell type (columns) in the
2-hop neighbourhoods of root cells of each type (rows); diagonal (self-proportion) masked.
Panel        : Extended Data Figures.pdf, Extended Data Fig. 4, panel c
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/neighborhoods/nh_2hop_compo_10.csv
                 (analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py)
Writes       : $KB_FIG_OUT/ED4/ED4c_2hop_celltype_heatmap.pdf (+ .png, ED4c_2hop_celltype_matrix.csv)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::ed3b (imported; the `stem` and
               `matrix_stem` keywords were added by patch rule build/patch_rules/fig_xenium.py, defaults unchanged)
               (original run: 05_xenium_figures.sh, SLURM jobs 329333 / 329334, 2026-09-02)
Illustrator  : title "Cell-type composition ..." written "Cell type composition ..."; colour-bar label
               "Proportions" -> "proportions"; dPT renamed aPT2
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/ed3b_2hop_celltype_heatmap.pdf
               (the matrix CSV is byte-identical to ed3b_2hop_celltype_matrix.csv)
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
    # ---- as in 05_xenium_figures.py::main (panel ed3b) ----
    compo = pd.read_csv(xf.OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0)
    levels = [c for c in compo.columns if c not in ("TMA_ID", "root_cell_type", "disease")]
    out = Path(fig_out("ED4", "ED4c_2hop_celltype_heatmap.pdf"))
    res = xf.ed3b(compo, levels, out.parent, stem=out.stem, matrix_stem="ED4c_2hop_celltype_matrix")
    print(res)
    print(f"[ED4c] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
