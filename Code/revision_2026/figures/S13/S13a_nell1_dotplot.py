#!/usr/bin/env python
"""Supp. Fig. 13a — NELL1 expression per snRNA-seq cell type (level2, 363,656 nuclei): dot colour = mean log1p(CP10k)
over the nuclei of the type, dot size and label = % of nuclei with a nonzero NELL1 count.
Panel        : All Supp Figures.pdf, Supplementary Fig. 13, panel a (= the left column, "snRNA NELL1", of the
               three-column dot plot drawn here; the snRNA APP and Xenium APP columns are not in the figure)
Reads        : $KB_ANALYSIS/nell1_app_egfr/tables/expression_by_celltype.csv and 01_meta.json
                 (analyses/05_nell1_snrna/01_nell1_app_egfr.py via 01_nell1_app_egfr.sh)
               only when those are missing, or with --from-objects: $KB_DATA/charlotte_new_data/
                 adata_20260717_with_clinical.h5ad (analyses/02_xenium_celltyping/18_add_clinical_to_snRNA_ref.py) and
                 $KB_ANALYSIS/celltype_concordance/outputs/xenium_labeled.h5ad (analyses/02_xenium_celltyping/), from
                 which only the per-cell-type mean and % expressing are computed (01_nell1_app_egfr.py::extract,
                 ::celltype_summary; no donor aggregation is needed for this panel)
Writes       : $KB_FIG_OUT/S13/S13a_nell1_app_expression_by_celltype.pdf (+ .png)
Plotting code: analyses/05_nell1_snrna/01_nell1_app_egfr.py::dotplot (imported unchanged; written under its fixed
               name, then renamed) (original run: 01_nell1_app_egfr.sh, SLURM job 329592, 2026-09-10)
Illustrator  : only the left column kept; column title "snRNA / NELL1" -> "NELL1 snRNA-seq"; suptitle and the panel
               note removed; colour bar relabelled "mean log1p(CP10K)"; dPT renamed aPT2
Environment  : conda env nichecompass
Verified     : 2026-11-12, PIXEL-IDENTICAL vs $KB_ANALYSIS/nell1_app_egfr/figures/nell1_app_expression_by_celltype.pdf,
               both from the saved tables and with --from-objects
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_lib"))
from xenium_figs import load_script, rename  # noqa: E402
from kbpaths import fig_out  # noqa: E402

na = load_script("05_nell1_snrna/01_nell1_app_egfr.py")


def from_objects():
    """Per-cell-type summaries as in 01_nell1_app_egfr.py::main, from the objects (no donor aggregation)."""
    so, sG, slib, sgenes, _ = na.load_snrna(False)
    s = na.celltype_summary(so, sG, slib, sgenes, "snRNA")
    a = na.ad.read_h5ad(na.XENIUM)
    xG, xlib, xgenes, _ = na.extract(a.layers["counts"], a.var_names)
    xo = pd.DataFrame({"cell_type": a.obs["level2_refined"].astype(str).to_numpy()}, index=a.obs_names)
    x = na.celltype_summary(xo, xG, xlib, xgenes, "Xenium")
    return pd.concat([s, x]), {"not_on_panel": [g for g in na.GENES if g not in xgenes]}


def main():
    ap = argparse.ArgumentParser(description="Supp. Fig. 13a")
    ap.add_argument("--from-objects", action="store_true", help="recompute the summaries even if the tables exist")
    args = ap.parse_args()
    tab = na.HERE / "tables"
    if not args.from_objects and (tab / "expression_by_celltype.csv").exists() and (tab / "01_meta.json").exists():
        summ = pd.read_csv(tab / "expression_by_celltype.csv")
        xinfo = json.loads((tab / "01_meta.json").read_text())["Xenium"]
        print(f"[S13a] summaries read from {tab}", flush=True)
    else:
        summ, xinfo = from_objects()
        print("[S13a] summaries computed from the objects", flush=True)
    np.random.seed(na.SEED)                       # global state as in 01_nell1_app_egfr.py::main
    sns.set_style("white")
    outdir = Path(fig_out("S13", "S13a_nell1_app_expression_by_celltype.pdf")).parent
    na.dotplot(summ, outdir, "", xinfo)
    for f in rename(outdir, {f"nell1_app_expression_by_celltype.{e}": f"S13a_nell1_app_expression_by_celltype.{e}"
                             for e in ("pdf", "png")}):
        print(f"[S13a] wrote {f}", flush=True)


if __name__ == "__main__":
    main()
