#!/usr/bin/env python
"""TF-centred eGRN panels: one TF at the centre, the genes it is inferred to regulate on a ring.

Co-author comment (2026-10-20) on ED Fig. 3 f/g: "in a proposed GRN I would prefer to highlight
one (or few) TF(s) and the genes it is inferred to regulate, as opposed to one gene and the TFs
inferred to regulate that gene." Companion to 47 (TF ranking) -- the inverse of the old panel g.

Per TF: all its r > 0.7 targets on a ring, ordered by correlation; Factor-1 top-positive targets
drawn in orange and labelled, other targets small grey (labelled only when the TF has <= 40
targets, otherwise unlabelled so the panel stays legible). Edge width = TF-gene correlation.
Title carries total targets, Factor-1 targets, the enrichment FDR from 47 and the activity peak bin.

Default TF set = main-text TFs present in the network (REL, FOSL2, FOS, JUN, NFKB1, NFKB2, RELB,
HNF4A, RXRA) + the three strongest Factor-1-enriched TFs from 47 (TEAD4, EHF, ELF3) + STAT3.
Override with --tfs. One PDF page per TF + a 2 x 3 composite of the first six.

Inputs: scmega/figures_paper/egrn_paper_views_edges.csv, egrn_factor1_tf_side.csv (script 47),
cis_regulatory_pt_apt/outputs/task3/C4_factor1_top200_positive_genes.csv.
Outputs: scmega/figures_paper/egrn_tf_centric/<TF>.{pdf,png}, egrn_tf_centric_all.pdf,
         egrn_tf_centric_composite.{pdf,png}, 48_meta.json
Human hg38, HGNC symbols. Seed 42 (layout deterministic). Env agent_env, login node.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

A = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC")
FP = A / "scmega" / "figures_paper"
OUT = FP  # output folder (release patch: the figure entry scripts in figures/ redirect it; default = original folder); inputs are read from FP
F1 = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3"
          "/C4_factor1_top200_positive_genes.csv")
C_F1, C_TARGET, C_TF = "#e8820c", "#b9b9b9", "#4c72b0"
DEFAULT = ["TEAD4", "EHF", "ELF3", "STAT3", "REL", "FOSL2", "FOS", "JUN", "NFKB1", "NFKB2", "RELB", "HNF4A", "RXRA"]


def panel(ax, tf, sub, f1set, row, label_all):
    sub = sub.sort_values("correlation", ascending=False).reset_index(drop=True)
    n = len(sub)
    ang = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    isf1 = sub.gene.isin(f1set).to_numpy()
    stag = np.ones(n)
    if n > 40:
        stag = np.where(isf1, 1.0, 0.86)                 # Factor-1 targets on the outer ring
    AX, AY = 1.25, 1.0
    x, y = AX * stag * np.cos(ang), AY * stag * np.sin(ang)
    lw = 0.35 + 2.4 * (sub.correlation - 0.7) / 0.3
    for xi, yi, w, f in zip(x, y, lw, isf1):
        ax.plot([0, xi], [0, yi], color=C_F1 if f else "#c9d6e3", lw=w if f else 0.6 * w,
                alpha=0.9 if f else 0.7, zorder=2 if f else 1, solid_capstyle="round")
    ax.scatter(x[~isf1], y[~isf1], s=22 if n <= 40 else 11, c=C_TARGET, edgecolors="none", zorder=3)
    ax.scatter(x[isf1], y[isf1], s=70, c=C_F1, edgecolors="black", linewidths=0.6, zorder=4)
    for xi, yi, a, g, f, r in zip(x, y, ang, sub.gene, isf1, sub.correlation):
        if not (f or label_all):
            continue
        ca, sa = np.cos(a), np.sin(a)
        ha = "left" if ca > 0.1 else ("right" if ca < -0.1 else "center")
        va = "bottom" if sa > 0.55 else ("top" if sa < -0.55 else "center")
        ax.annotate(g, (xi, yi), xytext=(7 * ca, 7 * sa), textcoords="offset points", ha=ha, va=va,
                    fontsize=6.4 if f else 5.4, fontweight="bold" if f else "normal",
                    color="#7a4300" if f else "0.35", zorder=5)
    ax.scatter([0], [0], s=820, c=C_TF, edgecolors="black", linewidths=1.0, zorder=6)
    ax.annotate(tf, (0, 0), fontsize=10, fontweight="bold", ha="center", va="center", color="white", zorder=7)
    nf = int(isf1.sum())
    fdr = row.hypergeom_fdr_bh
    ax.set_title(f"{tf} — {n} targets in the r > 0.7 eGRN, {nf} of them Factor-1 top-positive genes\n"
                 f"Factor-1 enrichment BH-FDR {fdr:.2g}{'  *' if fdr < 0.05 else ''}; motif activity peaks at bin "
                 f"{row.activity_peak_bin:.0f}; module {'I (injury)' if row.module == 1 else 'II (differentiated PT)'}",
                 fontsize=8.6, fontweight="bold")
    ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.32, 1.32); ax.set_axis_off()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tfs", default=",".join(DEFAULT))
    args = ap.parse_args()
    np.random.seed(42)
    E = pd.read_csv(FP / "egrn_paper_views_edges.csv")
    T = pd.read_csv(FP / "egrn_factor1_tf_side.csv").set_index("tf")
    f1set = set(pd.read_csv(F1)["variable"])
    tfs = [t for t in args.tfs.split(",") if t]
    missing = [t for t in tfs if t not in set(E.tf)]
    tfs = [t for t in tfs if t in set(E.tf)]
    if missing:
        print(f"[warn] not in the network (not selected by SelectTFs at r > 0.7): {', '.join(missing)}")
    out = OUT / "egrn_tf_centric"; out.mkdir(exist_ok=True)
    leg = [Line2D([], [], marker="o", ls="", mfc=C_TF, mec="black", ms=9, label="transcription factor"),
           Line2D([], [], marker="o", ls="", mfc=C_F1, mec="black", ms=7, label="target that is a Factor-1 top-positive gene"),
           Line2D([], [], marker="o", ls="", mfc=C_TARGET, mec="none", ms=5, label="other target"),
           Line2D([], [], color="#c9d6e3", lw=2, label="edge: motif in a linked peak, r > 0.7 (width = r)")]
    with PdfPages(OUT / "egrn_tf_centric_all.pdf") as pdf:
        for tf in tfs:
            sub = E[E.tf == tf]
            fig, ax = plt.subplots(figsize=(9.0, 7.4))
            panel(ax, tf, sub, f1set, T.loc[tf], label_all=len(sub) <= 40)
            fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=7.5)
            fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.07)
            pdf.savefig(fig, facecolor="white")
            for ext in ("pdf", "png"):
                fig.savefig(out / f"{tf}.{ext}", dpi=200, facecolor="white")
            plt.close(fig)
    six = tfs[:6]
    fig, axes = plt.subplots(2, 3, figsize=(19.5, 12.4))
    for ax, tf in zip(axes.ravel(), six):
        panel(ax, tf, E[E.tf == tf], f1set, T.loc[tf], label_all=False)
    for ax in axes.ravel()[len(six):]:
        ax.set_axis_off()
    fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle("TF-centred views of the r > 0.7 scMEGA eGRN (PT → aPT, ArchR diffusion-map pseudotime): "
                 "each TF with the genes it is inferred to regulate; Factor-1 top-positive targets labelled",
                 fontsize=12, y=0.995)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.05, wspace=0.02, hspace=0.12)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"egrn_tf_centric_composite.{ext}", dpi=170, facecolor="white")
    plt.close(fig)
    (OUT / "48_meta.json").write_text(json.dumps({
        "script": "48_egrn_tf_centric_panels.py", "tfs": tfs, "not_in_network": missing, "composite": six,
        "per_tf": {t: {"n_targets": int((E.tf == t).sum()), "n_factor1": int(E[(E.tf == t)].gene.isin(f1set).sum())} for t in tfs},
    }, indent=2))
    print(f"[done] {len(tfs)} panels -> {out}; composite of {six}")


if __name__ == "__main__":
    main()
