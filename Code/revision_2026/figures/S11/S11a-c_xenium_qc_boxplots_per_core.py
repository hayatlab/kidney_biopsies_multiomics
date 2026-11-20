#!/usr/bin/env python
"""Supplementary Fig. 11a-c — Xenium quality metrics per cell, one box per core, for the five runs:
a total transcripts per cell, b unique transcripts (genes detected) per cell, c cell area.
Panel        : All Supp Figures.pdf, Supplementary Fig. 11, panels a-c
Reads        : $KB_RESULTS/xenium/leon/new_data/overall_with_glom_04_level2/xenium_with_niches_overall_with_glom.h5ad
               (analyses/01_xenium_preprocessing_niches; 258,375 cells x 468 genes, raw counts, obs `sample` = run,
               `cell_area`, `condition`, core column `TMA_ID`, or `xenium_core_id` in the de-identified object)
Writes       : $KB_FIG_OUT/S11/S11a_total_transcripts_per_core.pdf, S11b_unique_transcripts_per_core.pdf,
               S11c_cell_area_per_core.pdf, and S11a-c_box_statistics.csv (median, quartiles and cell number per box)
Plotting code: analyses/14_xenium_qc_per_core/qc_per_core.py (method of the notebook `Supp_fig_4a-c.ipynb` of Sidrah
               Maryam: scanpy QC metrics, seaborn box plots; see the module for the differences from the notebook)
Illustrator  : the 15 plots arranged in three columns (a, b, c) and five rows (runs)
Environment  : conda env nichecompass (scanpy, seaborn)
Verified     : 2026-11-20, the boxes equal the printed panels (quartiles, median and maximum of all 357 boxes measured
               in the vector graphics of the printed page; 15, 33, 16, 28 and 27 cores per run). The tick labels differ
               from the print: this script names cores by public ids, orders the boxes by public core id and takes
               label and box from the same table row. Compare boxes by their statistics, not by position.
"""
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import KB_RESULTS, fig_out  # noqa: E402
import cores  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402

spec = importlib.util.spec_from_file_location("qc_per_core", REPO / "analyses" / "14_xenium_qc_per_core" / "qc_per_core.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

OBJ = os.path.join(KB_RESULTS, "xenium", "leon", "new_data", "overall_with_glom_04_level2", "xenium_with_niches_overall_with_glom.h5ad")
RUNS = ["XETG00229__0041734__Region_1__20250416__104148", "XETG00229__0041734__Region_3__20250416__104148",
        "XETG00229__0042348__Region_1__20250416__104148", "XETG00229__0042348__Region_2__20250416__104148",
        "XETG00229__0042348__Region_3__20250416__104148"]                      # rows of the figure, top to bottom

adata = sc.read_h5ad(OBJ)
if "counts" in adata.layers:
    adata.X = adata.layers["counts"].copy()
if "xenium_core_id" in adata.obs:                       # de-identified object
    core_ids = adata.obs["xenium_core_id"].astype(str).tolist()
    sample_ids = adata.obs["sample_id"].astype(str).tolist() if "sample_id" in adata.obs else core_ids
else:                                                   # internal object: translate with the private key, never print the label
    pub = cores.public_core_ids(adata.obs["TMA_ID"].astype(str))
    core_ids, sample_ids = [p[0] for p in pub], [p[1] for p in pub]

qc = m.qc_table(adata, core_ids, sample_ids)
assert set(qc["run"]) == set(RUNS)
n = qc.groupby("run")["xenium_core_id"].nunique().reindex(RUNS).tolist()
assert n == [15, 33, 16, 28, 27] and len(qc) == 258375, (n, len(qc))

names = {"a": "S11a_total_transcripts_per_core.pdf", "b": "S11b_unique_transcripts_per_core.pdf", "c": "S11c_cell_area_per_core.pdf"}
stats = []
for panel, metric, title, ylabel in m.METRICS:
    fig = m.plot_metric(qc, metric, title, ylabel, RUNS)
    out = fig_out("S11", names[panel])
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("written:", out)
    s = qc.groupby(["run", "xenium_core_id", "sample_id", "condition"])[metric].describe()[["count", "25%", "50%", "75%", "max"]].reset_index()
    s.insert(0, "panel", panel); s.insert(1, "metric", metric)
    stats.append(s.rename(columns={"count": "n_cells", "25%": "q1", "50%": "median", "75%": "q3"}))
import pandas as pd  # noqa: E402
pd.concat(stats).to_csv(fig_out("S11", "S11a-c_box_statistics.csv"), index=False)
