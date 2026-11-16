#!/usr/bin/env python
"""TF-side view of the Factor-1 target genes in the r > 0.7 scMEGA eGRN (Ivan's comment on ED Fig. 3g).

Ivan (2026-10-16, on ED Fig. 3 panel g = ITGB3 with 39 regulators): "a reviewer might indicate that
our analysis is too loose. Biologically it is unlikely that so many TFs can regulate a single gene.
I do like F a lot. What I would do instead is to get now all the genes in F and their TFs, plot
the final network and evaluate page-rank/degree statistics. ... A simpler alternative is just to
have a similar plot as f) with the TFs upstream of these 46 genes."

Inputs (all produced by 43_scmega_egrn_tutorial.R / 45_egrn_paper_views.py, run archr_dm):
  scmega/figures_paper/egrn_paper_views_edges.csv   3,686 TF->gene edges, r > 0.7, with module
  scmega/figures_paper/egrn_factor1_regulators.csv  the 678 edges into the 46 Factor-1 targets
  scmega/egrn_tutorial/archr_dm/tf_centrality.csv   scMEGA centrality (PageRank, betweenness, peak bin)
  scmega/egrn_tutorial/archr_dm/grn_edges_all.csv   motif-in-linked-peak candidates before the r filter
  scmega/egrn_tutorial/archr_dm/peak_to_gene_links.csv
  cis_regulatory_pt_apt/outputs/task3/C4_factor1_top200_positive_genes.csv  (80 genes, PT view loadings)

What is computed
  1. Sub-network S = every r > 0.7 edge whose target is one of the 46 Factor-1 genes (69 TFs, 678 edges).
  2. Per TF: number of the 46 genes it regulates (= out-degree in S), total out-degree in the full
     network, hypergeometric enrichment of Factor-1 genes among its targets (population = 388 targets,
     46 of them Factor-1; BH over the 86 TFs of the network), weighted out-degree, HITS hub score in S,
     PageRank on the REVERSED sub-network (gene -> TF), and scMEGA's own PageRank/importance on the
     full directed network for comparison.
  3. Per gene: in-degree in S (= the bars of panel f), HITS authority score.
  4. Edge-inflation diagnostics behind Ivan's point: fraction of motif-in-linked-peak candidates that
     pass r > 0.7, linked peaks per gene, candidates and passes for ITGB3.

Why not plain PageRank on the TF -> gene graph: in a bipartite graph where TFs only have out-edges,
PageRank of every TF equals the teleport mass (1 - d)/N -- identical for all of them. scMEGA's
importance is non-zero only for the 10 TFs that are themselves targets of another TF (123 TF -> TF
edges); 76 of 86 TFs get importance 0 (RUN_LOG 2026-09-18). PageRank on the reversed graph and the
HITS hub score are the versions of "PageRank" that rank regulators here; both are functions of
(weighted) out-degree and of how many *other* TFs share each target, so they must be read as degree
statistics, not as evidence of hierarchy.

Outputs (scmega/figures_paper/):
  egrn_factor1_tf_side.csv, egrn_factor1_gene_side.csv, egrn_factor1_subnetwork_edges.csv
  egrn_factor1_tf_ranking.{pdf,png}       the "panel f for TFs" (Ivan's simpler alternative)
  egrn_factor1_subnetwork.{pdf,png}       bipartite network, 69 TFs -> 46 genes
  egrn_factor1_tf_specificity.{pdf,png}   Factor-1 targets vs all targets per TF, with enrichment
  egrn_factor1_adjacency.{pdf,png}        clustered TF x gene adjacency matrix (the readable form of the hairball)
  47_egrn_factor1_tf_side_meta.json

Human hg38, HGNC symbols. Seed 42 (layouts deterministic). Env agent_env, login node (seconds).
Usage: python scripts/47_egrn_factor1_tf_side.py [--highlight TF1,TF2,...]
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
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy import stats  # noqa: E402

A = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC")
FP = A / "scmega" / "figures_paper"
OUT = FP  # output folder (release patch: the figure entry scripts in figures/ redirect it; default = original folder); inputs are read from FP
RUN = A / "scmega" / "egrn_tutorial" / "archr_dm"
F1 = Path(f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/cis_regulatory_pt_apt/outputs/task3"
          "/C4_factor1_top200_positive_genes.csv")
SEED = 42
C_F1, C_TARGET = "#e8820c", "#c9c9c9"
# TFs the manuscript / scMEGA step already singles out (scMEGA importance >= 0.1 in run archr_dm
# plus the NF-kB members in the network); overridden by --highlight
DEFAULT_HIGHLIGHT = ["ELF3", "STAT1", "FOSL2", "KLF6", "CREB5", "HNF4A", "VDR", "NFAT5", "EGR1",
                     "RXRA", "RELB", "REL", "NFKB1", "NFKB2"]


def bh(p):
    return stats.false_discovery_control(np.asarray(p, float), method="bh")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--highlight", default=",".join(DEFAULT_HIGHLIGHT))
    args = ap.parse_args()
    hl = [t for t in args.highlight.split(",") if t]
    np.random.seed(SEED)

    E = pd.read_csv(FP / "egrn_paper_views_edges.csv")           # full r > 0.7 network
    C = pd.read_csv(RUN / "tf_centrality.csv").set_index("tf")
    Eall = pd.read_csv(RUN / "grn_edges_all.csv")                # candidates before r filter
    P2G = pd.read_csv(RUN / "peak_to_gene_links.csv")
    f1 = pd.read_csv(F1).set_index("variable")["Loading"]
    tfs_all = sorted(set(E.tf)); targets_all = sorted(set(E.gene))
    f1set = sorted(set(f1.index) & set(targets_all))
    S = E[E.gene.isin(f1set)].copy()
    assert len(S) > 0 and len(f1set) > 0, (len(S), S.tf.nunique(), len(f1set))  # Sept sizes were 678 edges / 69 TFs / 46 genes; rerun sizes printed below
    print(f"[sub] {len(S)} edges, {S.tf.nunique()} TFs -> {len(f1set)} Factor-1 genes "
          f"(full network {len(E)} edges, {len(tfs_all)} TFs, {len(targets_all)} targets)")

    # ---------------------------------------------------------------- per-TF table
    N, K = len(targets_all), len(f1set)
    out_full = E.groupby("tf").size()
    t = pd.DataFrame(index=pd.Index(tfs_all, name="tf"))
    t["n_factor1_targets"] = S.groupby("tf").size().reindex(t.index).fillna(0).astype(int)
    t["n_targets_total"] = out_full.reindex(t.index).astype(int)
    t["frac_targets_factor1"] = t.n_factor1_targets / t.n_targets_total
    t["expected_factor1"] = t.n_targets_total * K / N
    t["fold_enrichment"] = t.n_factor1_targets / t.expected_factor1
    t["hypergeom_p"] = [stats.hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
                        for k, n in zip(t.n_factor1_targets, t.n_targets_total)]
    t["hypergeom_fdr_bh"] = bh(t.hypergeom_p)
    t["weighted_outdeg_factor1"] = S.groupby("tf").correlation.sum().reindex(t.index).fillna(0)
    t["mean_r_factor1"] = S.groupby("tf").correlation.mean().reindex(t.index)
    t["module"] = E.groupby("tf").module.first().reindex(t.index)
    for c in ["activity_peak_bin", "pagerank_raw", "importance", "importance_rank", "n_regulating_tfs",
              "motif_expression_r", "pct_PT", "pct_aPT"]:
        t[f"scmega_{c}" if c in ("pagerank_raw", "importance", "importance_rank") else c] = C[c].reindex(t.index)

    # graph statistics on S: HITS hubs (TFs) / authorities (genes); PageRank on the reversed graph
    G = nx.DiGraph()
    G.add_weighted_edges_from(zip(S.tf, S.gene, S.correlation))
    hubs, auth = nx.hits(G, max_iter=1000, tol=1e-10)
    prr = nx.pagerank(G.reverse(copy=True), alpha=0.85, weight="weight")
    prf = nx.pagerank(G, alpha=0.85, weight="weight")
    t["hits_hub_S"] = pd.Series(hubs).reindex(t.index).fillna(0)
    t["pagerank_reversed_S"] = pd.Series(prr).reindex(t.index).fillna(0)
    t["pagerank_forward_S"] = pd.Series(prf).reindex(t.index).fillna(0)
    # rank columns (1 = top) among the 69 TFs in S
    inS = t.n_factor1_targets > 0
    for c in ["n_factor1_targets", "weighted_outdeg_factor1", "hits_hub_S", "pagerank_reversed_S"]:
        t.loc[inS, f"rank_{c}"] = t.loc[inS, c].rank(ascending=False, method="min").astype(int)
    t["highlight"] = t.index.isin(hl)
    t = t.sort_values(["n_factor1_targets", "weighted_outdeg_factor1"], ascending=False)
    t.to_csv(OUT / "egrn_factor1_tf_side.csv")

    g = pd.DataFrame(index=pd.Index(f1set, name="gene"))
    g["factor1_loading"] = f1.reindex(g.index)
    g["n_regulators"] = S.groupby("gene").size().reindex(g.index)
    g["n_candidate_tfs_motif"] = Eall[Eall.gene.isin(f1set)].groupby("gene").size().reindex(g.index)
    g["n_linked_peaks"] = P2G[P2G.gene.isin(f1set)].groupby("gene").size().reindex(g.index)
    g["hits_authority_S"] = pd.Series(auth).reindex(g.index)
    g["module"] = S.groupby("gene").module.first().reindex(g.index)
    g = g.sort_values("n_regulators", ascending=False)
    g.to_csv(OUT / "egrn_factor1_gene_side.csv")
    S.assign(tf_activity_peak_bin=S.tf.map(C.activity_peak_bin)).to_csv(
        OUT / "egrn_factor1_subnetwork_edges.csv", index=False)

    # ---------------------------------------------------------------- inflation diagnostics
    pass_rate = float((Eall.correlation > 0.7).mean())
    cand = Eall.groupby("gene").size()
    diag = {
        "candidate_edges_motif_in_linked_peak": int(len(Eall)),
        "candidate_edges_passing_r_gt_0.7": int((Eall.correlation > 0.7).sum()),
        "pass_rate": round(pass_rate, 4),
        "pass_rate_positive_candidates": round(float((Eall.correlation[Eall.correlation > 0] > 0.7).mean()), 4),
        "linked_peaks_per_gene_median": float(P2G.groupby("gene").size().median()),
        "linked_peaks_per_gene_max": int(P2G.groupby("gene").size().max()),
        "candidate_tfs_per_gene_median": float(cand.median()),
        "regulators_per_target_median_all": float(E.groupby("gene").size().median()),
        "regulators_per_target_median_factor1": float(S.groupby("gene").size().median()),
        "ITGB3": {"linked_peaks": int((P2G.gene == "ITGB3").sum()),
                  "candidate_tfs": int((Eall.gene == "ITGB3").sum()),
                  "regulators_r_gt_0.7": int((E.gene == "ITGB3").sum()),
                  "candidate_tfs_per_peak": "not recoverable from the exported tables (motif matrix lives in the Seurat object)"},
    }
    print("[diag]", json.dumps(diag, indent=1))

    cmap = plt.get_cmap("viridis"); norm = plt.Normalize(1, 100)
    tt = t[inS].copy()
    n = len(tt)
    cols = [cmap(norm(b)) if np.isfinite(b) else "#999999" for b in tt.activity_peak_bin]

    # ---------------------------------------------------------------- figure 1: TF ranking (panel f for TFs)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.2, 0.215 * n + 2.0), sharey=True,
                                  gridspec_kw={"width_ratios": [1.0, 0.62], "wspace": 0.04})
    y = np.arange(n)
    ax.barh(y, tt.n_factor1_targets, color=cols, edgecolor="black", linewidth=0.3)
    ax.set_yticks(y); ax.set_yticklabels(tt.index, fontsize=7.2)
    for lab, h in zip(ax.get_yticklabels(), tt.highlight):
        if h:
            lab.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlabel(f"number of the {K} Factor-1 target genes the TF regulates (r > 0.7 eGRN)", fontsize=8.5)
    ax.set_title(f"TFs upstream of the {K} Factor-1 genes ({n} TFs)",
                 fontsize=10, fontweight="bold", loc="left")
    ax.set_ylim(n - 0.4, -0.6)
    # right: total targets (grey) with the Factor-1 share, and the enrichment marker
    ax2.barh(y, tt.n_targets_total, color="#dcdcdc", edgecolor="none")
    ax2.barh(y, tt.n_factor1_targets, color="#8a8a8a", edgecolor="none")
    exp = tt.expected_factor1
    ax2.plot(exp, y, color="black", lw=0.8, ls=":", zorder=3)
    sig = tt.hypergeom_fdr_bh < 0.05
    ax2.scatter(tt.n_targets_total[sig] + 3, y[sig], marker="*", s=34, color="#D55E00", zorder=4)
    ax2.set_xlabel("all targets of the TF (grey) and the Factor-1 share (dark)", fontsize=8.5)
    ax2.set_title("all targets vs Factor-1 share", fontsize=10, fontweight="bold", loc="left")
    for a in (ax, ax2):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        a.tick_params(axis="x", labelsize=7.5)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = fig.add_axes([0.13, 0.035, 0.30, 0.012])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("pseudotime bin of peak TF motif activity (1 = PT-like, 100 = aPT-like)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    fig.legend(handles=[
        Line2D([], [], color="black", lw=0.8, ls=":",
               label=f"expected Factor-1 targets if the TF's targets were drawn at random ({K}/{N} = {K / N:.1%})"),
        Line2D([], [], marker="*", ls="", ms=8, color="#D55E00",
               label="Factor-1 genes over-represented among the TF's targets (hypergeometric, BH-FDR < 0.05 over 86 TFs)"),
        Line2D([], [], ls="", label="bold = TF named in the manuscript main text")],
        loc="lower right", bbox_to_anchor=(0.995, 0.005), frameon=False, fontsize=7, handlelength=1.6)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.955, bottom=0.085)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"egrn_factor1_tf_ranking.{ext}", dpi=200, facecolor="white")
    plt.close(fig)

    # ---------------------------------------------------------------- figure 2: bipartite network
    tf_order = tt.sort_values(["activity_peak_bin", "n_factor1_targets"], ascending=[True, False]).index.tolist()
    gene_order = g.sort_values("factor1_loading", ascending=False).index.tolist()
    xt = np.linspace(0, 1, len(tf_order)); xg = np.linspace(0.04, 0.96, len(gene_order))
    pos = {tf: (x, 1.0) for tf, x in zip(tf_order, xt)}
    pos.update({ge: (x, 0.0) for ge, x in zip(gene_order, xg)})
    fig, ax = plt.subplots(figsize=(15.5, 7.2))
    segs = [[pos[a], pos[b]] for a, b in zip(S.tf, S.gene)]
    alpha = 0.10 + 0.55 * (S.correlation - 0.7) / 0.3
    ax.add_collection(LineCollection(segs, colors=[(0.45, 0.55, 0.65, a) for a in alpha],
                                     linewidths=0.45, zorder=1))
    stf = 40 + 420 * tt.n_factor1_targets.reindex(tf_order) / tt.n_factor1_targets.max()
    ctf = [cmap(norm(b)) if np.isfinite(b) else "#999999" for b in tt.activity_peak_bin.reindex(tf_order)]
    ax.scatter(xt, np.ones(len(xt)), s=stf, c=ctf, edgecolors="black", linewidths=0.6, zorder=3)
    sg = 30 + 300 * g.n_regulators.reindex(gene_order) / g.n_regulators.max()
    ax.scatter(xg, np.zeros(len(xg)), s=sg, c=C_F1, edgecolors="black", linewidths=0.6, zorder=3)
    for tf, x in zip(tf_order, xt):
        ax.annotate(tf, (x, 1.0), xytext=(0, 12), textcoords="offset points", rotation=90,
                    ha="center", va="bottom", fontsize=6.8,
                    fontweight="bold" if tf in hl else "normal")
    for ge, x in zip(gene_order, xg):
        ax.annotate(ge, (x, 0.0), xytext=(0, -12), textcoords="offset points", rotation=90,
                    ha="center", va="top", fontsize=6.8, color="#7a4300")
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.42, 1.42); ax.set_axis_off()
    ax.text(-0.025, 1.0, "TFs\n(ordered by\npeak bin)", ha="right", va="center", fontsize=8)
    ax.text(-0.025, 0.0, "Factor-1\ngenes\n(by loading)", ha="right", va="center", fontsize=8)
    ax.set_title(f"eGRN sub-network of the {K} Factor-1 top-positive target genes: {n} TFs, {len(S)} edges "
                 f"(r > 0.7, motif in a linked peak), all in module I\n"
                 f"TF size = number of Factor-1 targets; gene size = number of regulators; "
                 f"edge opacity = TF-gene correlation", fontsize=10.5)
    cax = fig.add_axes([0.40, 0.045, 0.22, 0.014])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("pseudotime bin of peak TF motif activity (1 = PT-like, 100 = aPT-like)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.08)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"egrn_factor1_subnetwork.{ext}", dpi=200, facecolor="white")
    plt.close(fig)

    # ---------------------------------------------------------------- figure 3: specificity scatter
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    xx = np.linspace(0, t.n_targets_total.max() * 1.05, 50)
    ax.plot(xx, xx * K / N, color="black", lw=0.9, ls=":", label=f"expected at random ({K / N:.1%})")
    allt = t.copy()
    ax.scatter(allt.n_targets_total, allt.n_factor1_targets,
               s=28, c=[cmap(norm(b)) if np.isfinite(b) else "#999999" for b in allt.activity_peak_bin],
               edgecolors=np.where(allt.hypergeom_fdr_bh < 0.05, "#D55E00", "black"),
               linewidths=np.where(allt.hypergeom_fdr_bh < 0.05, 1.6, 0.5), zorder=3)
    lab = allt[(allt.hypergeom_fdr_bh < 0.05) | allt.highlight | (allt.n_factor1_targets >= 17)]
    placed = []
    for tf, r in lab.sort_values(["n_factor1_targets", "n_targets_total"]).iterrows():
        k = sum(1 for (px, py) in placed if abs(px - r.n_targets_total) < 9 and abs(py - r.n_factor1_targets) < 0.8)
        dx, dy = (3, 3) if k == 0 else ((3, -9) if k == 1 else (-3, 9))
        ha = "left" if k != 2 else "right"
        ax.annotate(tf, (r.n_targets_total, r.n_factor1_targets), xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.8, ha=ha, fontweight="bold" if r.highlight else "normal")
        placed.append((r.n_targets_total, r.n_factor1_targets))
    ax.set_xlabel("all targets of the TF in the r > 0.7 eGRN"); ax.set_ylabel(f"of which Factor-1 top-positive genes ({K} in the network)")
    ax.set_title("Which TFs are specific to the Factor-1 genes, and which are simply hubs?", fontsize=10.5, fontweight="bold")
    ax.legend(handles=[Line2D([], [], color="black", ls=":", label=f"expected at random ({K / N:.1%} of targets)"),
                       Line2D([], [], marker="o", ls="", mfc="white", mec="#D55E00", mew=1.6, ms=7,
                              label="over-represented, BH-FDR < 0.05 (86 TFs)"),
                       Line2D([], [], ls="", label="bold = highlighted TF; colour = peak bin")],
              frameon=False, fontsize=8, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"egrn_factor1_tf_specificity.{ext}", dpi=200, facecolor="white")
    plt.close(fig)

    # ---------------------------------------------------------------- figure 4: clustered adjacency heatmap
    from scipy.cluster import hierarchy as hc
    M = S.pivot_table(index="tf", columns="gene", values="correlation", aggfunc="max").reindex(index=tt.index, columns=f1set)
    B = M.notna().astype(float)
    ro = hc.leaves_list(hc.linkage(B.to_numpy(), "average", metric="jaccard"))
    co = hc.leaves_list(hc.linkage(B.to_numpy().T, "average", metric="jaccard"))
    M = M.iloc[ro, co]
    fig, ax = plt.subplots(figsize=(0.26 * K + 2.6, 0.19 * n + 1.8))
    im = ax.imshow(M.to_numpy(), cmap="Blues", vmin=0.7, vmax=1.0, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(K)); ax.set_xticklabels(M.columns, rotation=90, fontsize=6.8, color="#7a4300")
    ax.set_yticks(range(n)); ax.set_yticklabels(M.index, fontsize=6.8)
    for lab_, h in zip(ax.get_yticklabels(), t.highlight.reindex(M.index)):
        if h:
            lab_.set_fontweight("bold")
    nf = t.n_factor1_targets.reindex(M.index)
    for i, v in enumerate(nf):
        ax.text(K - 0.3, i, str(int(v)), fontsize=6, va="center", ha="left", color="0.3")
    ax.set_xlim(-0.5, K + 1.2)
    ax.set_title(f"TF x Factor-1-gene adjacency of the r > 0.7 eGRN ({n} TFs x {K} genes, {len(S)} edges)\n"
                 f"rows/columns clustered on edge presence (Jaccard, average linkage); cell = TF-gene correlation; "
                 f"right margin = number of Factor-1 targets", fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.06)
    cb.set_label("TF-gene correlation r (edges only; white = no edge)", fontsize=7); cb.ax.tick_params(labelsize=6.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"egrn_factor1_adjacency.{ext}", dpi=200, facecolor="white")
    plt.close(fig)

    # ---------------------------------------------------------------- provenance
    top = tt.head(15)
    (OUT / "47_egrn_factor1_tf_side_meta.json").write_text(json.dumps({
        "script": "47_egrn_factor1_tf_side.py", "seed": SEED, "run": "archr_dm",
        "sub_network": {"edges": int(len(S)), "tfs": int(n), "genes": int(K)},
        "full_network": {"edges": int(len(E)), "tfs": len(tfs_all), "targets": N},
        "enrichment": {"test": "hypergeometric, population = 388 targets, 46 Factor-1 targets; BH over 86 TFs",
                       "n_fdr_lt_0.05": int((t.hypergeom_fdr_bh < 0.05).sum()),
                       "tfs_fdr_lt_0.05": t.index[t.hypergeom_fdr_bh < 0.05].tolist()},
        "top15_by_n_factor1_targets": top.n_factor1_targets.to_dict(),
        "top10_hits_hub": tt.sort_values("hits_hub_S", ascending=False).head(10).index.tolist(),
        "top10_pagerank_reversed": tt.sort_values("pagerank_reversed_S", ascending=False).head(10).index.tolist(),
        "spearman_rank_nF1_vs_hits_hub": float(stats.spearmanr(tt.n_factor1_targets, tt.hits_hub_S).statistic),
        "spearman_rank_nF1_vs_pagerank_reversed": float(stats.spearmanr(tt.n_factor1_targets, tt.pagerank_reversed_S).statistic),
        "forward_pagerank_tf_sd": float(tt.pagerank_forward_S.std()),
        "highlight": hl, "diagnostics": diag,
        "networkx": nx.__version__, "pandas": pd.__version__,
    }, indent=2, default=float))
    print(t[["n_factor1_targets", "n_targets_total", "fold_enrichment", "hypergeom_fdr_bh", "hits_hub_S",
             "pagerank_reversed_S", "activity_peak_bin", "scmega_importance"]].head(25).round(3).to_string())
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
