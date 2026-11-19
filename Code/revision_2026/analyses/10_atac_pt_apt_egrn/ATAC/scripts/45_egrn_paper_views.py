#!/usr/bin/env python
"""Readable eGRN views of the r > 0.7 scMEGA network, for the paper's supplement.

Leon (2026-09-25): "in the atac archr trajectory we see two distinct networks in the 0.7
version -- give me an eGRN view in figures_paper for the supplementary, and for some of our
top positive Factor 1 genes a pruned version of the network to those individually, showing
one at a time."

Source network: scmega/egrn_tutorial/archr_dm/ -- script 43 with Mina's thresholds
(SelectTFs cor 0.7 / p 0.01, SelectGenes r 0.7 / FDR 1e-50, GRN edge r 0.7), on the ArchR
diffusion-map trajectory (17,637 PT->aPT pairs). 3,686 positive edges, 86 TFs, 388 targets.

The network has exactly TWO connected components, and they are not arbitrary:

  module I   310 nodes, 72 TFs -- AP-1 (FOS/FOSL1/FOSL2/JUN/JUND), NF-kB (NFKB1/NFKB2/REL/
             RELB), STAT1/3, IRF1, KLF/ELF/TEAD. The injury programme.
  module II  154 nodes, 14 TFs -- HNF4A, RXRA, VDR, ESRRA, NR2F1/6, RARA, RORC, NR1I3, PPAR-
             type nuclear receptors. The differentiated proximal-tubule metabolic programme.

They share no edge, which is the point of the figure: the two programmes are regulated by
disjoint TF sets over this trajectory.

Panel A of the full view draws each component separately (Kamada-Kawai per component) with
TFs coloured by the pseudotime bin where their motif activity peaks (`activity_peak_bin`,
1-100, early = PT-like, late = aPT-like) and sized by out-degree. Targets are small and grey,
except the Factor-1 top-positive genes, which are drawn in orange and labelled -- the link to
the per-gene panels.

Per-gene views: for every Factor-1 top-positive gene that is a target in the network, one
panel with the gene at the centre and its regulating TFs on a ring ordered by activity peak
bin. Edge width is the TF-gene correlation over the smoothed pseudotime bins; the number
beside each TF is how many linked peaks carry that TF's motif (>= 1 by construction).
Written one gene per file AND as a single multipage PDF.

Factor-1 genes: the positive half of the top 200 |loading| genes of the PT view of MCFA
Factor 1 (cis_regulatory_pt_apt/outputs/task3/C4_factor1_top200_positive_genes.csv, 80 genes;
46 of them are targets in this network).

Outputs (scmega/figures_paper/):
  egrn_archr_dm_two_modules.{pdf,png}
  egrn_factor1_targets_overview.{pdf,png}
  egrn_factor1_genes/<GENE>.{pdf,png}   and  egrn_factor1_genes_all.pdf
  egrn_paper_views_nodes.csv, egrn_paper_views_edges.csv, egrn_factor1_regulators.csv,
  45_egrn_views_meta.json

Human hg38, HGNC symbols. Seed 42.  Env: agent_env (light; login node).
Usage: python scripts/45_egrn_paper_views.py [--run archr_dm] [--genes GENE1,GENE2]
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import argparse, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

A = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC")
F1 = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3"
          "/C4_factor1_top200_positive_genes.csv")
OUT = A / "scmega" / "figures_paper"  # output folder (release patch: the figure entry scripts in figures/ redirect it; default = original folder)
SEED = 42
C_TARGET, C_F1, CMAP = "#c9c9c9", "#e8820c", "yellow_green_blue"
MODULE_NAMES = {0: "module I — AP-1 / NF-κB / STAT injury programme",
                1: "module II — HNF4A / RXRA / VDR differentiated-PT programme"}


def load(run):
    d = A / "scmega" / "egrn_tutorial" / run
    E = pd.read_csv(d / "grn_edges.csv").rename(columns={"weights": "correlation"})
    Eall = pd.read_csv(d / "grn_edges_all.csv")
    E = E.merge(Eall[["tf", "gene", "n_peaks"]], on=["tf", "gene"], how="left")
    C = pd.read_csv(d / "tf_centrality.csv").set_index("tf")
    meta = json.loads((d / "meta.json").read_text())
    return d, E, C, meta


def components(E):
    G = nx.from_pandas_edgelist(E, "tf", "gene", create_using=nx.Graph)
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    which = {n: i for i, c in enumerate(comps) for n in c}
    return G, comps, which


def tf_colour(C, tfs):
    """Colour TFs by the pseudotime bin where their motif activity peaks (early -> late)."""
    bins = C["activity_peak_bin"].reindex(tfs).astype(float)
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=1, vmax=100)
    return [cmap(norm(b)) if np.isfinite(b) else "#999999" for b in bins], cmap, norm


def draw_component(ax, E, tfset, f1set, C, title, ls=1.0):
    G = nx.from_pandas_edgelist(E, "tf", "gene", create_using=nx.Graph)
    pos = nx.kamada_kawai_layout(G)
    xy = pd.DataFrame(pos, index=["x", "y"]).T
    segs = [[xy.loc[a].to_numpy(), xy.loc[b].to_numpy()] for a, b in zip(E.tf, E.gene)]
    ax.add_collection(LineCollection(segs, colors="#b8b8b8", linewidths=0.30, alpha=0.45,
                                     zorder=1))
    tg = [n for n in G if n not in tfset]
    plain = [n for n in tg if n not in f1set]
    hit = [n for n in tg if n in f1set]
    ax.scatter(xy.loc[plain, "x"], xy.loc[plain, "y"], s=9, c=C_TARGET, linewidths=0, zorder=2)
    if hit:
        ax.scatter(xy.loc[hit, "x"], xy.loc[hit, "y"], s=46, c=C_F1, zorder=4,
                   edgecolors="white", linewidths=0.7)
    tfs = [n for n in G if n in tfset]
    cols, cmap, norm = tf_colour(C, tfs)
    deg = E.groupby("tf").size().reindex(tfs).fillna(0)
    size = 45 + 330 * (deg / max(deg.max(), 1))
    ax.scatter(xy.loc[tfs, "x"], xy.loc[tfs, "y"], s=size, c=cols, zorder=5,
               edgecolors="black", linewidths=0.55)
    for n in tfs:
        ax.annotate(n, (xy.loc[n, "x"], xy.loc[n, "y"]), fontsize=7.2 * ls, fontweight="bold",
                    xytext=(0, 7 * ls), textcoords="offset points", ha="center", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.72))
    for n in hit:
        ax.annotate(n, (xy.loc[n, "x"], xy.loc[n, "y"]), fontsize=6.4 * ls, color="#8a4b00",
                    xytext=(0, -9 * ls), textcoords="offset points", ha="center", zorder=6)
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left")
    ax.set_axis_off(); ax.autoscale_view()
    return cmap, norm


def gene_panel(ax, gene, sub, C, loading):
    """One Factor-1 gene at the centre, its regulating TFs on a ring ordered by peak bin."""
    sub = sub.copy()
    sub["bin"] = C["activity_peak_bin"].reindex(sub.tf).to_numpy()
    sub = sub.sort_values(["bin", "correlation"], na_position="last").reset_index(drop=True)
    n = len(sub)
    # elliptical ring, wider than tall (labels run horizontally, so horizontal room is what
    # is scarce), with every other node pushed out by 12 % so neighbouring labels never abut
    ang = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    if n > 24:      # three radial levels keep even 39 regulators' labels apart
        stag = np.array([1.0, 1.13, 1.26])[np.arange(n) % 3]
    elif n > 12:
        stag = np.where(np.arange(n) % 2 == 0, 1.0, 1.13)
    else:
        stag = np.ones(n)
    AX, AY = 1.30, 1.00
    x, y = AX * stag * np.cos(ang), AY * stag * np.sin(ang)
    lw = 0.6 + 3.2 * (sub.correlation - 0.7) / max(sub.correlation.max() - 0.7, 1e-9)
    for xi, yi, w in zip(x, y, lw):
        ax.plot([0, xi], [0, yi], color="#8fb8d8", lw=w, zorder=1, alpha=0.85,
                solid_capstyle="round")
    cols, cmap, norm = tf_colour(C, sub.tf.tolist())
    ax.scatter(x, y, s=190, c=cols, edgecolors="black", linewidths=0.7, zorder=3)
    for xi, yi, a, tf, npk, r in zip(x, y, ang, sub.tf, sub.n_peaks, sub.correlation):
        ca, sa = np.cos(a), np.sin(a)
        ha = "left" if ca > 0.10 else ("right" if ca < -0.10 else "center")
        va = "bottom" if sa > 0.55 else ("top" if sa < -0.55 else "center")
        ax.annotate(f"{tf}\nr {r:.2f} · {int(npk)} peak{'s' if npk != 1 else ''}",
                    (xi, yi), fontsize=6.6, ha=ha, va=va,
                    xytext=(13 * ca, 13 * sa), textcoords="offset points", zorder=5)
    ax.scatter([0], [0], s=760, c=C_F1, edgecolors="black", linewidths=1.0, zorder=4)
    ax.annotate(gene, (0, 0), fontsize=10, fontweight="bold", ha="center", va="center",
                color="white", zorder=6)
    lim = 1.42 + 0.72 * stag.max()
    ax.set_xlim(-lim * 1.38, lim * 1.38); ax.set_ylim(-lim, lim)
    ax.set_title(f"{gene} — {n} regulator{'s' if n != 1 else ''} in the r > 0.7 eGRN"
                 + (f"\nFactor 1 PT loading {loading:.2f}" if np.isfinite(loading) else ""),
                 fontsize=10.5, fontweight="bold")
    ax.set_axis_off()
    return cmap, norm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="archr_dm")
    ap.add_argument("--genes", default="", help="comma-separated override of the gene panels")
    ap.add_argument("--label-scale", type=float, default=1.0,
                    help="scale TF / target label font size in the two-module network (13 Nov 2026: 1.4 for the "
                         "revised ED3m); != 1 writes egrn_archr_dm_two_modules_labels<scale>x.*")
    args = ap.parse_args()
    np.random.seed(SEED)

    d, E, C, meta = load(args.run)
    G, comps, which = components(E)
    tfset = set(E.tf)
    f1 = pd.read_csv(F1).set_index("variable")["Loading"]
    f1set = set(f1.index) & set(E.gene)
    out = Path(OUT); out.mkdir(parents=True, exist_ok=True)
    print(f"[net] {args.run}: {len(E)} edges, {len(tfset)} TFs, {E.gene.nunique()} targets, "
          f"{len(comps)} components {[len(c) for c in comps]}", flush=True)
    print(f"[f1 ] {len(f1)} Factor-1 top-positive genes, {len(f1set)} of them targets here",
          flush=True)

    # ------------------------------------------------------------------ full view
    fig, axes = plt.subplots(1, len(comps), figsize=(13.0 + 8.6 * (len(comps) - 1), 11.8),
                             gridspec_kw={"width_ratios": [1.0] + [0.66] * (len(comps) - 1)})
    axes = np.atleast_1d(axes)
    for i, comp in enumerate(comps):
        sub = E[E.tf.isin(comp) & E.gene.isin(comp)]
        n_tf = sub.tf.nunique(); n_tg = sub.gene.nunique()
        n_f1 = len(set(sub.gene) & f1set)
        cmap, norm = draw_component(
            axes[i], sub, tfset, f1set, C,
            f"{MODULE_NAMES.get(i, f'module {i+1}')}\n{len(comp)} nodes — {n_tf} TFs, "
            f"{n_tg} target genes, {len(sub):,} edges — {n_f1} Factor-1 top-positive genes",
            ls=args.label_scale)
    sm = plt.cm.ScalarMappable(cmap=plt.get_cmap("viridis"), norm=plt.Normalize(1, 100))
    cax = fig.add_axes([0.34, 0.075, 0.32, 0.013])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("pseudotime bin of peak TF motif activity   (1 = PT-like  →  100 = aPT-like)",
                 fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=10, mfc="#5ec962", mec="black", mew=0.6,
               label="transcription factor (size = number of targets, fill = activity peak bin)"),
        Line2D([], [], marker="o", ls="", ms=6, mfc=C_TARGET, mec="none", label="target gene"),
        Line2D([], [], marker="o", ls="", ms=9, mfc=C_F1, mec="white",
               label="target gene that is a Factor-1 top-positive gene"),
        Line2D([], [], color="#b8b8b8", lw=1.4,
               label="TF–gene edge: motif in a linked peak, correlation r > 0.7 over the "
                     "smoothed pseudotime")],
        loc="lower center", ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(
        f"scMEGA eGRN on the ArchR PT→aPT diffusion-map trajectory, r > 0.7 — two disjoint modules\n"
        f"{len(E):,} positive TF–gene edges, {len(tfset)} TFs, {E.gene.nunique()} targets over "
        f"{meta['pairs']:,} RNA–ATAC pairs; SelectTFs cor {meta['thresholds']['SelectTFs_cor']} / "
        f"p {meta['thresholds']['SelectTFs_p']}, peak-to-gene r "
        f"{meta['thresholds']['SelectGenes_cor']} / FDR {meta['thresholds']['SelectGenes_fdr']:g}, "
        f"edge r {meta['thresholds']['GRN_edge_cor']}.  The two components share no edge, and every "
        f"Factor-1 top-positive gene in the network sits in module I.",
        fontsize=13, y=0.985)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.905, bottom=0.115, wspace=0.02)
    stem = ("egrn_archr_dm_two_modules" if args.label_scale == 1.0
            else f"egrn_archr_dm_two_modules_labels{args.label_scale:g}x")
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=200, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[fig] {out}/{stem}.pdf", flush=True)

    # ------------------------------------------------------------------ per-gene panels
    genes = ([g.strip() for g in args.genes.split(",")] if args.genes
             else [g for g in f1.sort_values(ascending=False).index if g in f1set])
    missing = [g for g in genes if g not in set(E.gene)]
    if missing:
        print(f"[warn] not targets in this network, skipped: {', '.join(missing)}", flush=True)
    genes = [g for g in genes if g in set(E.gene)]
    gdir = out / "egrn_factor1_genes"; gdir.mkdir(exist_ok=True)

    reg_rows = []
    with PdfPages(out / "egrn_factor1_genes_all.pdf") as pdf:
        for g in genes:
            sub = E[E.gene == g]
            fig, ax = plt.subplots(figsize=(9.4, 8.0))
            cmap, norm = gene_panel(ax, g, sub, C, float(f1.get(g, np.nan)))
            cax = fig.add_axes([0.37, 0.055, 0.26, 0.014])
            cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                              orientation="horizontal")
            cb.set_label("pseudotime bin of peak TF motif activity (1 = PT-like → 100 = aPT-like)",
                         fontsize=6.6)
            cb.ax.tick_params(labelsize=6)
            ax.text(0.5, -0.045, f"module {which[g] + 1} — {MODULE_NAMES.get(which[g], '')}\n"
                                f"edge width = TF–gene correlation; TF fill = pseudotime bin of "
                                f"peak motif activity",
                    transform=ax.transAxes, ha="center", va="top", fontsize=7.6, color="0.35")
            fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.14)
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            for ext in ("pdf", "png"):
                fig.savefig(gdir / f"{g}.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            r = sub.copy()
            r["gene_factor1_loading"] = float(f1.get(g, np.nan))
            r["module"] = which[g] + 1
            r["tf_activity_peak_bin"] = C["activity_peak_bin"].reindex(r.tf).to_numpy()
            r["tf_n_targets"] = C["n_targets"].reindex(r.tf).to_numpy()
            reg_rows.append(r)
    reg = pd.concat(reg_rows, ignore_index=True) if reg_rows else pd.DataFrame()
    reg.to_csv(out / "egrn_factor1_regulators.csv", index=False)
    print(f"[fig] {len(genes)} per-gene panels in {gdir} (+ egrn_factor1_genes_all.pdf)",
          flush=True)

    # ------------------------------------------------------------------ overview of the F1 genes
    top = (reg.groupby("gene")
              .agg(n_regulators=("tf", "nunique"), max_r=("correlation", "max"),
                   mean_r=("correlation", "mean"), module=("module", "first"),
                   loading=("gene_factor1_loading", "first"))
              .sort_values("n_regulators", ascending=False))
    fig, ax = plt.subplots(figsize=(9.5, 0.26 * len(top) + 2.2))
    col = np.where(top["module"] == 1, "#c44e52", "#4c72b0")
    ax.barh(range(len(top)), top["n_regulators"], color=col, edgecolor="none")
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=7.5)
    ax.invert_yaxis(); ax.set_xlabel("number of TFs regulating the gene in the r > 0.7 eGRN")
    ax.set_title(f"Factor-1 top-positive genes in the eGRN ({len(top)} of {len(f1)} are targets)",
                 fontsize=11, fontweight="bold")
    ax.legend(handles=[Line2D([], [], color="#c44e52", lw=6, label="module I (injury)"),
                       Line2D([], [], color="#4c72b0", lw=6, label="module II (differentiated PT)")],
              frameon=False, fontsize=9, loc="lower right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"egrn_factor1_targets_overview.{ext}", dpi=200, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)

    # ------------------------------------------------------------------ tables + provenance
    nodes = pd.DataFrame({"id": list(G)})
    nodes["is_tf"] = nodes.id.isin(tfset)
    nodes["module"] = nodes.id.map(which) + 1
    nodes["out_degree"] = nodes.id.map(E.groupby("tf").size()).fillna(0).astype(int)
    nodes["in_degree"] = nodes.id.map(E.groupby("gene").size()).fillna(0).astype(int)
    nodes["factor1_top_positive"] = nodes.id.isin(set(f1.index))
    nodes["factor1_loading"] = nodes.id.map(f1)
    for c in ["activity_peak_bin", "n_targets", "motif_expression_r", "pct_PT", "pct_aPT",
              "importance", "importance_rank"]:
        nodes[c] = nodes.id.map(C[c]) if c in C.columns else np.nan
    nodes.to_csv(out / "egrn_paper_views_nodes.csv", index=False)
    E.assign(module=E.tf.map(which) + 1).to_csv(out / "egrn_paper_views_edges.csv", index=False)

    (out / "45_egrn_views_meta.json").write_text(json.dumps({
        "script": "45_egrn_paper_views.py", "run": args.run, "source": str(d), "seed": SEED,
        "thresholds": meta["thresholds"], "pairs": meta["pairs"],
        "n_edges": int(len(E)), "n_tfs": int(len(tfset)), "n_targets": int(E.gene.nunique()),
        "components": [{"module": i + 1, "n_nodes": len(c),
                        "n_tfs": len(set(c) & tfset),
                        "tfs": sorted(set(c) & tfset),
                        "n_factor1_top_positive": len(set(c) & f1set)}
                       for i, c in enumerate(comps)],
        "factor1_gene_list": str(F1), "n_factor1_genes": int(len(f1)),
        "n_factor1_genes_in_network": int(len(f1set)),
        "factor1_genes_not_in_network": sorted(set(f1.index) - f1set),
        "layout": "kamada_kawai per component; per-gene panels = ring ordered by TF activity peak bin",
        "networkx": nx.__version__, "pandas": pd.__version__,
    }, indent=2))
    print(f"[done] {out}", flush=True)


if __name__ == "__main__":
    main()
