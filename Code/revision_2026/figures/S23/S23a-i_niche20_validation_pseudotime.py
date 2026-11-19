#!/usr/bin/env python
"""Supp. Fig. 23a-i — niche-20 net score along the crescent-formation trajectory of an external glomerular dataset (GSE294965, Sultana et al.).
Panel        : All Supp Figures.pdf, Supplementary Figure 23, panels a-i
               a-d PCA of glomeruli coloured by disease / DPT pseudotime / k-means cluster / niche-20 net score
               e-f % disease per k-means cluster / per pseudotime quadrant; g pseudotime vs net score; h-i net score per
               k-means cluster / per quadrant
Reads        : $KB_ANALYSIS/validation_data/niche20_score_validation/option1/adata_polygon_reduced_with_score.h5ad
               (analyses/08_glomerular_niche20_crescent/niche20_validation_pseudotime.py: polygon-level AnnData with PCA,
               DPT pseudotime, quadrants, k-means clusters and the niche-20 up/down/net scores)
               $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv (de_crescent_pec_mast.R;
               only to rebuild the signature gene lists, checked against the saved signature_genes_used.tsv)
               $KB_ANALYSIS/validation_data/niche20_score_validation/option1/summary.tsv (statistics checked against it)
Writes       : $KB_FIG_OUT/S23/S23a-d_pca_overview.pdf, S23e_heatmap_disease_kmeans.pdf, S23f_heatmap_disease_quadrant.pdf,
               S23g_scatter_pseudotime_vs_niche20.pdf, S23h-i_score_by_paper_grouping.pdf
Plotting code: extracted verbatim from analyses/08_glomerular_niche20_crescent/niche20_validation_pseudotime.py, section
               "5. Plots" (the script is not importable: it runs the whole analysis at import). Statistics in the titles
               (Spearman, Kruskal-Wallis) are recomputed from the saved table with the original code and asserted equal
               to the saved summary.tsv. Nothing upstream (polygon reduction, PCA, DPT, k-means) is recomputed.
               The strip-plot jitter of h-i draws from NumPy's global RNG, whose state in the original run was set by
               the two sc.tl.score_genes(random_state=42) calls; they are re-run on the saved object only to restore
               that state (their scores are asserted equal to the saved ones; the saved scores are plotted).
               The 95% CI whiskers of h-i are seaborn bootstraps with an unseeded generator and differ run to run.
Illustrator  : a: citation "Sultana Z, Khatri R, Yousefi B, et al. (PMID: 41028563)" added; b: title "dpt_pseudotime
               (crescent formation trajectory)" with arrow; c: title shortened to "KMeans"; d: title "niche 20 net score
               from this study"; g: "polygon-level: " removed from the title; i: "Kruskal" added to the title
Environment  : conda env spatial_env
Verified     : 2026-11-12, a-d, e, f, g PIXEL-IDENTICAL; h-i NEAR (0.004% of pixels, at three CI whisker tips)
               vs $KB_ANALYSIS/validation_data/niche20_score_validation/option1/{pca_overview,heatmap_disease_kmeans,
               heatmap_disease_quadrant,scatter_pseudotime_vs_niche20,score_by_paper_grouping}.pdf
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "config"))
from kbpaths import KB_ANALYSIS, fig_out  # noqa: E402

import warnings  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.stats import kruskal, spearmanr  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)

# constants as in niche20_validation_pseudotime.py
TAG = "niche20"
NET = f"{TAG}_net"
SEED = 42
DISEASE_ORDER = ["Cntrl", "SLE", "ANCA", "GBM"]
DISEASE_PALETTE = {"Cntrl": "#0173B2", "SLE": "#E3A5D6", "ANCA": "#029E73", "GBM": "#D55E00"}
RUN = Path(f"{KB_ANALYSIS}/validation_data/niche20_score_validation/option1")
DE_PATH = f"{KB_ANALYSIS}/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv"


def restore_rng_state(adata):
    """Re-run the two score_genes calls of the original (same gene lists, same seed) on the saved object, so that NumPy's
    global RNG is in the state the original strip plots started from; check the scores against the saved ones."""
    de = pd.read_csv(DE_PATH, sep=",")
    panel = set(adata.var_names)
    sig_up = [g for g in de[(de["avg_log2FC"] > 0) & (de["p_val_adj"] < 0.05)].sort_values("avg_log2FC", ascending=False)["gene"] if g in panel]
    sig_dn = [g for g in de[(de["avg_log2FC"] < 0) & (de["p_val_adj"] < 0.05)].sort_values("avg_log2FC")["gene"] if g in panel]
    saved = pd.read_csv(RUN / "signature_genes_used.tsv", sep="\t", header=None, index_col=0)[1]
    assert sig_up == saved["niche20_up"].split(";") and sig_dn == saved["niche20_dn"].split(";"), "signature lists differ"
    tmp = adata.copy()
    sc.tl.score_genes(tmp, gene_list=sig_up, score_name="chk_up", random_state=SEED, ctrl_size=min(50, len(sig_up) * 5))
    sc.tl.score_genes(tmp, gene_list=sig_dn, score_name="chk_dn", random_state=SEED, ctrl_size=min(50, len(sig_dn) * 5))
    for new, old in (("chk_up", f"{TAG}_up"), ("chk_dn", f"{TAG}_dn")):
        assert np.allclose(tmp.obs[new].values, adata.obs[old].values, rtol=0, atol=1e-12), f"{new} != saved {old}"


def main() -> None:
    sc.settings.verbosity = 0
    adata = sc.read_h5ad(RUN / "adata_polygon_reduced_with_score.h5ad")
    restore_rng_state(adata)

    # ---------- statistics shown in the titles (code of section 4 of the original) ----------
    rs_pt, p_pt = spearmanr(adata.obs[NET], adata.obs["Pseudotime"])
    groups_km = [adata.obs[NET][adata.obs["KMeans_cluster"] == c].dropna() for c in ["C1", "C2", "C3", "C4"]]
    H_km, p_km = kruskal(*groups_km)
    groups_q = [adata.obs[NET][adata.obs["Quadrant"] == q].dropna() for q in ["Q1", "Q2", "Q3", "Q4"]]
    H_q, p_q = kruskal(*groups_q)
    summ = pd.read_csv(RUN / "summary.tsv", sep="\t", header=None, index_col=0)[1]
    for k, val in (("spearman_score_vs_pseudotime", rs_pt), ("spearman_score_vs_pseudotime_p", p_pt),
                   ("kruskal_kmeans_p", p_km), ("kruskal_quadrant_p", p_q)):
        assert np.isclose(float(summ[k]), val, rtol=1e-9, atol=0), f"{k}: {val} vs saved {summ[k]}"
    print(f"[S23] Spearman rho={rs_pt:.3f} p={p_pt:.2e}; Kruskal k-means p={p_km:.2e}; quadrants p={p_q:.2e}")
    ct_km = pd.crosstab(adata.obs["Disease"], adata.obs["KMeans_cluster"])
    ct_q = pd.crosstab(adata.obs["Disease"], adata.obs["Quadrant"])

    # ---------- 5. Plots (verbatim; output paths changed) ----------
    plt.rcParams["pdf.fonttype"] = 42

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2))
    sc.pl.pca(adata, color="Disease", ax=axes[0], show=False, palette=DISEASE_PALETTE,
              title="Disease", legend_loc="right margin")
    sc.pl.pca(adata, color="Pseudotime", ax=axes[1], show=False, cmap="viridis",
              title="dpt_pseudotime")
    sc.pl.pca(adata, color="KMeans_cluster", ax=axes[2], show=False, title="KMeans (PC1-ordered C1->C4)")
    sc.pl.pca(adata, color=NET, ax=axes[3], show=False, cmap="RdBu_r",
              title="niche-20 net score")
    fig.tight_layout()
    fig.savefig(fig_out("S23", "S23a-d_pca_overview.pdf"), bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.6))
    for ax, key, order, title in [
        (axes[0], "KMeans_cluster", ["C1", "C2", "C3", "C4"], f"by kmeans cluster (Kruskal p={p_km:.1e})"),
        (axes[1], "Quadrant",       ["Q1", "Q2", "Q3", "Q4"], f"by pseudotime quadrant (p={p_q:.1e})"),
    ]:
        df_p = adata.obs[[key, NET]].dropna()
        df_p[key] = pd.Categorical(df_p[key], categories=order, ordered=True)
        sns.violinplot(data=df_p, x=key, y=NET, order=order, palette="rocket_r",
                       cut=0, inner=None, ax=ax, linewidth=0.5)
        sns.stripplot(data=df_p, x=key, y=NET, order=order, size=2, color="black", alpha=0.5, ax=ax)
        sns.pointplot(data=df_p, x=key, y=NET, order=order, color="black",
                      errorbar=("ci", 95), markers="_", linestyles="", ax=ax)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("niche-20 net score")
        ax.axhline(0, color="grey", lw=0.5, ls=":")
    fig.tight_layout()
    fig.savefig(fig_out("S23", "S23h-i_score_by_paper_grouping.pdf"))
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    df_p = adata.obs[["Pseudotime", NET, "Disease"]].dropna()
    disease_cats = [c for c in DISEASE_ORDER if c in df_p["Disease"].unique()]
    df_p["Disease"] = pd.Categorical(df_p["Disease"], categories=disease_cats, ordered=True)
    sns.scatterplot(data=df_p, x="Pseudotime", y=NET,
                    hue="Disease", hue_order=disease_cats, palette=DISEASE_PALETTE,
                    s=10, alpha=0.8, ax=ax)
    xs = df_p["Pseudotime"].values
    ys = df_p[NET].values
    b, a = np.polyfit(xs, ys, 1)
    xx = np.linspace(xs.min(), xs.max(), 50)
    ax.plot(xx, a + b * xx, color="black", lw=1, ls="--")
    ax.set_title(f"polygon-level: ρ={rs_pt:.2f} (Spearman, p={p_pt:.1e})", fontsize=9)
    ax.legend(fontsize=7, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(fig_out("S23", "S23g_scatter_pseudotime_vs_niche20.pdf"))
    plt.close(fig)

    def heatmap_pct(crosstab, title, out_pdf):
        pct = crosstab.div(crosstab.sum(0), axis=1) * 100
        fig, ax = plt.subplots(figsize=(4, 3))
        sns.heatmap(pct, annot=True, fmt=".0f", cmap="Reds", ax=ax, cbar_kws={"label": "% of column"})
        ax.set_title(title, fontsize=9)
        fig.tight_layout()
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
    heatmap_pct(ct_q, "% disease per pseudotime quadrant", fig_out("S23", "S23f_heatmap_disease_quadrant.pdf"))
    heatmap_pct(ct_km, "% disease per kmeans cluster", fig_out("S23", "S23e_heatmap_disease_kmeans.pdf"))


if __name__ == "__main__":
    main()
