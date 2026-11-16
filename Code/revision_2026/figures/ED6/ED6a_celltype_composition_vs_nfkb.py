#!/usr/bin/env python
"""ED 6a — cell-type composition of the 2-hop Voronoi neighbourhood of every Xenium cell against the NF-kB activity
(PROGENy NFkB, top-500 weights, ULM) of its root cell: stacked spline ribbon as in Fig. 3g, x cropped to [-1.52, 1.43]
(the Factor-1 crop rule applied to this score: contiguous 0.0686-unit bins with >= 500 neighbourhoods, 98.05%).
Panel        : Extended Data Figures.pdf, Extended Data Fig. 6, panel a
Reads        : $KB_ANALYSIS/vivien_takeover/outputs/neighborhoods/nh_2hop_compo_10.csv and nh_2hop_niche_compo_10.csv
                 (analyses/03_xenium_factor1_neighbourhoods_cci/03_voronoi_neighborhoods.py)
               $KB_ANALYSIS/vivien_takeover/outputs/xenium_nfkb_scores.parquet
                 (analyses/03_xenium_factor1_neighbourhoods_cci/26_score_nfkb_xenium.py --fetch-nets, via
                  26_score_nfkb_xenium.sh; PROGENy weights fetched from OmniPath at run time)
               $KB_ANALYSIS/vivien_takeover/figures/xenium/13_crop_justification_meta_nfkb_progeny.json
                 (13_crop_justification.py --score NFkB_progeny --stem nfkb_progeny)
Writes       : $KB_FIG_OUT/ED6/ED6a_celltype_composition_vs_nfkb.pdf (+ .png, _fitted.csv)
Plotting code: analyses/03_xenium_factor1_neighbourhoods_cci/05_xenium_figures.py::ribbon, called as in
               27_ed3f_nfkb.py::main (both imported unchanged; only the cell-type / PROGENy panel is drawn)
               (original run: python 27_ed3f_nfkb.py, SLURM job 329585, 2026-09-10 -> fig3h_celltype_vs_nfkb)
Illustrator  : title removed; "(cropped to [-1.5, 1.4])" removed from the x label; y label "neighbourhood" ->
               "neighborhood"; legend title removed and legend cut at the right edge; PT / Immune / FIB / aPT written
               on the ribbons; dPT renamed aPT2
Environment  : conda env spatialdata
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/vivien_takeover/figures/xenium/fig3h_celltype_vs_nfkb.pdf
               (the _fitted.csv is byte-identical to fig3h_celltype_vs_nfkb_fitted.csv)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script  # noqa: E402
from kbpaths import fig_out  # noqa: E402

COL = "NFkB_progeny"
nf = load_script("03_xenium_factor1_neighbourhoods_cci/27_ed3f_nfkb.py")   # imports 05_xenium_figures.py as nf.xf


def main():
    np.random.seed(nf.SEED)                       # global state as in 27_ed3f_nfkb.py::main
    sns.set_style("white")
    # ---- extracted from 27_ed3f_nfkb.py::main (score NFkB_progeny, cell-type table) ----
    nz = pd.read_csv(nf.OUT / "neighborhoods" / "nh_2hop_niche_compo_10.csv", index_col=0)
    ct = pd.read_csv(nf.OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0)
    clv = [c for c in ct.columns if c not in ("TMA_ID", "root_cell_type", "disease")]
    if not nz.index.equals(ct.index):
        raise SystemExit("niche and cell-type composition tables have different roots")
    sc_ = pd.read_parquet(nf.OUT / "xenium_nfkb_scores.parquet").set_index("sample_cell_id")
    stem13, _, _, xlab, suffix = nf.SCORES[COL]
    score = sc_[COL].reindex(nz.index)
    print(f"[{COL}] {int(score.isna().sum())} neighbourhoods have no root score -> left out of the fit", flush=True)
    crop, note, _ = nf.crop_for(stem13)
    ct[COL] = score
    out = Path(fig_out("ED6", "ED6a_celltype_composition_vs_nfkb.pdf"))
    res = nf.xf.ribbon(ct, clv, "Cell type", COL,
                       f"Cell-type composition of 2-hop neighbourhoods vs NF-κB activity{suffix}",
                       nf.YLAB, out.parent, out.stem, cmap_name="tab20", crop=crop, xlabel=xlab, crop_note=note)
    print(res)
    print(f"[ED6a] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
