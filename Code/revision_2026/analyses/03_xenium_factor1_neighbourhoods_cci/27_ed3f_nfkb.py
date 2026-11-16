#!/usr/bin/env python
"""ED 3f / Fig 3h variants: 2-hop neighbourhood composition vs NF-kB activity (Leon, 2026-09-10).

Same construction as ed3f_niche_vs_factor1 and fig3h_celltype_vs_factor1
(05_xenium_figures.py::ribbon, imported, not copied): per-level LSQUnivariateSpline(k=3,
8 interior knots) of the 2-hop neighbourhood proportions against the ROOT cell's score,
clipped to [0, 1], renormalised, stacked. Only the x-axis changes: the root cell's NF-kB
activity from 26_score_nfkb_xenium.py instead of Factor 1. Leon asked for the niche version
first, then "i wanted to compare it to cell type as well" -> the Fig 3h analogue too, with
the same palette as fig3h_celltype_vs_factor1 so the two can be read side by side.

Crop: the Factor-1 rule (Leon, 2026-09-07) -- the largest contiguous run of bins each holding
>= 500 neighbourhoods -- applied to each NF-kB score by 13_crop_justification.py, with a bin
width giving the same resolution as Factor 1 (the same number of bins across the 1st-99th
percentile). Niche and cell-type panels share the same 258,375 roots, so they share the crop.

Outputs (figures/xenium/):
  ed3f_niche_vs_nfkb.{pdf,png}                 niches, PROGENy NFkB (primary)
  ed3f_niche_vs_nfkb_collectri.{pdf,png}       niches, CollecTRI NF-kB (sensitivity)
  fig3h_celltype_vs_nfkb.{pdf,png}             cell types, PROGENy NFkB (primary)
  fig3h_celltype_vs_nfkb_collectri.{pdf,png}   cell types, CollecTRI NF-kB (sensitivity)
  nfkb_score_context.{pdf,png}                 NF-kB per cell type, and NF-kB vs Factor 1
  27_ed3f_nfkb_meta.json
Human GRCh38, HGNC symbols. Seed 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import importlib.util, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

HERE = Path(f"{KB_ANALYSIS}/vivien_takeover")
OUT, FIGX = HERE / "outputs", HERE / "figures" / "xenium"
SEED, THR = 42, 500
# score column -> (13 crop stem, niche stem, cell-type stem, axis label, title suffix)
SCORES = {
    "NFkB_progeny": ("nfkb_progeny", "ed3f_niche_vs_nfkb", "fig3h_celltype_vs_nfkb",
                     "NF-κB activity (PROGENy, ULM)", ""),
    "NFkB_collectri": ("nfkb_collectri", "ed3f_niche_vs_nfkb_collectri",
                       "fig3h_celltype_vs_nfkb_collectri",
                       "NF-κB TF activity (CollecTRI, ULM)", "  [sensitivity: CollecTRI]"),
}
YLAB = "Proportions in 2-hop neighbourhood"

_spec = importlib.util.spec_from_file_location("xenium_figures",
                                               Path(__file__).resolve().parent / "05_xenium_figures.py")
xf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(xf)


def crop_for(stem13):
    m = json.loads((FIGX / f"13_crop_justification_meta_{stem13}.json").read_text())
    c = m["candidates"][str(THR)]
    note = (f"largest contiguous run of {m['bin_width']:.3g}-unit bins with >= {THR} "
            f"neighbourhoods each (13_crop_justification.py --stem {stem13})")
    return (c["lo"], c["hi"]), note, m


def context_figure(sc_, f1, xlab):
    d = sc_.join(f1).dropna(subset=["NFkB_progeny"])
    order = (d.groupby("level2_refined")["NFkB_progeny"].median()
              .sort_values(ascending=False).index)
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5), gridspec_kw={"width_ratios": [1.35, 1]})
    sns.boxplot(data=d, x="level2_refined", y="NFkB_progeny", order=order, ax=axes[0],
                showfliers=False, color="#9ecae1", linewidth=0.8)
    axes[0].tick_params(axis="x", rotation=90)
    axes[0].set_xlabel(""); axes[0].set_ylabel(xlab, fontsize=12)
    axes[0].set_title("NF-κB activity per cell type  (box = IQR, whiskers 1.5 IQR, "
                      "outliers hidden)", fontsize=12)
    dd = d.dropna(subset=["Factor1_ulm"])
    hb = axes[1].hexbin(dd["Factor1_ulm"], dd["NFkB_progeny"], gridsize=90, bins="log",
                        cmap="viridis", mincnt=1, linewidths=0)
    rho = float(spearmanr(dd["Factor1_ulm"], dd["NFkB_progeny"]).statistic)
    axes[1].set_xlabel("Factor1_ulm", fontsize=12); axes[1].set_ylabel(xlab, fontsize=12)
    axes[1].set_title(f"NF-κB vs Factor 1, per cell  (Spearman ρ = {rho:.2f}, "
                      f"n = {len(dd):,})", fontsize=12)
    fig.colorbar(hb, ax=axes[1], label="cells (log)")
    for ax in axes:
        sns.despine(ax=ax)
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(FIGX / f"nfkb_score_context.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"spearman_nfkb_vs_factor1": rho, "n_cells": int(len(dd)),
            "celltype_order_by_median": list(order)}


def main():
    np.random.seed(SEED)
    sns.set_style("white")
    nz = pd.read_csv(OUT / "neighborhoods" / "nh_2hop_niche_compo_10.csv", index_col=0)
    nlv = [c for c in nz.columns if c not in ("TMA_ID", "root_niche", "disease")]
    ct = pd.read_csv(OUT / "neighborhoods" / "nh_2hop_compo_10.csv", index_col=0)
    clv = [c for c in ct.columns if c not in ("TMA_ID", "root_cell_type", "disease")]
    if not nz.index.equals(ct.index):
        raise SystemExit("niche and cell-type composition tables have different roots")
    sc_ = pd.read_parquet(OUT / "xenium_nfkb_scores.parquet").set_index("sample_cell_id")
    f1 = pd.read_parquet(OUT / "xenium_factor_scores.parquet",
                         columns=["sample_cell_id", "Factor1_ulm"]).set_index("sample_cell_id")
    print(f"[load] {len(nz)} 2-hop neighbourhoods, {len(nlv)} niches, {len(clv)} cell types; "
          f"{len(sc_)} scored cells", flush=True)

    meta = {}
    for col, (stem13, stem_n, stem_c, xlab, suffix) in SCORES.items():
        score = sc_[col].reindex(nz.index)
        n_nan = int(score.isna().sum())
        print(f"[{col}] {n_nan} neighbourhoods ({n_nan/len(nz)*100:.3f}%) have no root score "
              f"-> left out of the fit", flush=True)
        crop, note, m13 = crop_for(stem13)
        nz[col], ct[col] = score, score
        for tab, lv, legend, what, stem, cmap in [
                (nz, nlv, "Niche", "Niche", stem_n, "tab10"),
                (ct, clv, "Cell type", "Cell-type", stem_c, "tab20")]:
            r = xf.ribbon(tab, lv, legend, col,
                          f"{what} composition of 2-hop neighbourhoods vs NF-κB activity{suffix}",
                          YLAB, FIGX, stem, cmap_name=cmap, crop=crop, xlabel=xlab,
                          crop_note=note)
            r.update({"n_nhoods_without_score": n_nan, "bin_width": m13["bin_width"]})
            meta[stem] = r

    meta["context"] = context_figure(sc_, f1, SCORES["NFkB_progeny"][3])
    (FIGX / "27_ed3f_nfkb_meta.json").write_text(json.dumps(
        {"panels": meta, "crop_threshold": THR, "seed": SEED,
         "versions": {"matplotlib": matplotlib.__version__, "seaborn": sns.__version__,
                      "numpy": np.__version__, "pandas": pd.__version__}},
        indent=2, default=str))
    print(f"[done] {FIGX}", flush=True)


if __name__ == "__main__":
    main()
