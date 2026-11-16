"""
Reproduce Sultana et al. 2025 (notebook 08_pca_pseudotime/) polygon-level pipeline
and test whether the per-glom **niche-20** PEC signature recapitulates their
pseudotime / k-means crescent-like grouping.

This is `niche9_validation_pseudotime.py` with ONE substantive change: the signature source.
  niche9  : niche9_pec/DE/niche9_vs_others_pecOnly/markers_MAST_niche9_vs_others_pecOnly.tsv
            (in_glom_niche run, latent_leiden_0.6, niche 9; DE universe = 313 genes)
  niche20 : celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv
            (manual_glom run at res 1.0 -> crescent = niche 20; DE universe = 468 genes)
Everything else -- polygon reduction, PCA, DPT, quadrants, k-means, stats, plots -- is byte-for-byte
the same logic, so any difference in the result is attributable to the signature alone.

*** READ THIS BEFORE INTERPRETING THE OUTPUT ***
Restricted to the GSE294965 480-gene panel, the two signatures are ~92% IDENTICAL:
    UP  : 23 genes each, Jaccard 0.92 -- differ only by COL16A1 (niche9) vs IRF8 (niche20)
    DOWN: 22 genes each, Jaccard 0.91 -- the same two genes, swapped the other way
i.e. COL16A1 and IRF8 trade directions and nothing else moves. So this rerun is expected to
reproduce the niche9 numbers almost exactly (niche9: rho=0.41 vs pseudotime). That is a
CONSISTENCY/ROBUSTNESS check -- it is NOT independent validation of the niche20 signature, and it must
not be reported as a second, corroborating result. The reason the panels collapse together is that the
niche9 DE universe (313 genes) is a strict SUBSET of the niche20 DE universe (468), and the GSE294965
panel only intersects 47 / 60 of those respectively.
The overlap is recomputed and written to `signature_overlap_vs_niche9.tsv` at runtime rather than
trusted from this comment.

Faithful to Sultana et al. (compute_pseudotime.ipynb, cluster.ipynb):
  1. data = data[data.obs.is_in_polygon]
  2. drop NaN Biopsy_ID
  3. df = data.to_df(layer="counts"); df["polygon_flags"]; groupby("polygon_flags").mean()
  4. AnnData(df) -> normalize_total -> log1p -> sc.tl.pca -> sc.pp.neighbors
  5. iroot = first Cntrl polygon -> sc.tl.dpt
  6. quadrants = pd.cut(pseudotime, 4)
  7. KMeans(n_clusters=4, random_state=42, n_init=10) on X_pca

Outputs: review_analyses/validation_data/niche20_score_validation/option1/
(the niche9 outputs under niche9_score_validation/option1/ are left untouched for comparison)
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.stats import spearmanr, pearsonr, kruskal, mannwhitneyu

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 2

TAG = "niche20"
H5AD = f"{KB_ANALYSIS}/validation_data/GSE294965_processed_data.h5ad"
# niche-20 = manual_glom NicheCompass run at res 1.0, crescent niche 20 (241 cells, 79% PEC+Imm);
# see celltype_concordance/PLAN.md "Crescent-PEC signature concordance (2026-07-10)".
DE_PATH = f"{KB_ANALYSIS}/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv"
DE_SEP = ","                       # niche20 file is CSV; the niche9 file was TSV
DE_NICHE9 = f"{KB_ANALYSIS}/niche9_pec/DE/niche9_vs_others_pecOnly/markers_MAST_niche9_vs_others_pecOnly.tsv"
OUT = Path(f"{KB_ANALYSIS}/validation_data/niche20_score_validation/option1")
OUT.mkdir(exist_ok=True, parents=True)

SEED = 42
np.random.seed(SEED)
DISEASE_ORDER = ["Cntrl", "SLE", "ANCA", "GBM"]
DISEASE_PALETTE = {"Cntrl": "#0173B2", "SLE": "#E3A5D6", "ANCA": "#029E73", "GBM": "#D55E00"}

# ---------- 1. Reproduce polygon-reduction ----------
print(f"[load] {H5AD}")
data = sc.read_h5ad(H5AD)
print(f"[load] full: {data.shape}")
data = data[data.obs["is_in_polygon"]]
print(f"[load] is_in_polygon: {data.shape}")
data = data[data.obs["Biopsy_ID"].astype(str) != "nan"]
print(f"[load] biopsy not NaN: {data.shape}")

print("[agg] groupby('polygon_flags') of counts layer (mean)")
df = data.to_df(layer="counts")
df["polygon_flags"] = data.obs["polygon_flags"].astype(str).tolist()
df = df.groupby("polygon_flags").agg("mean")
print(f"[agg] polygon-level matrix: {df.shape}")

adata = AnnData(df)
adata.obs["Disease"] = [e.split("_")[-2] for e in adata.obs.index]
print(f"[agg] disease counts: {adata.obs['Disease'].value_counts().to_dict()}")

# ---------- 2. PCA / neighbors / pseudotime / kmeans ----------
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata)

adata.obs["Disease"] = adata.obs["Disease"].astype("category").cat.reorder_categories(
    [c for c in DISEASE_ORDER if c in adata.obs["Disease"].unique()]
)
adata.uns["iroot"] = int(np.flatnonzero(adata.obs["Disease"] == "Cntrl")[0])
sc.tl.dpt(adata)
adata.obs["Pseudotime"] = adata.obs["dpt_pseudotime"]
adata.obs["Quadrant"] = pd.cut(adata.obs["Pseudotime"], 4, labels=["Q1", "Q2", "Q3", "Q4"])

km = KMeans(n_clusters=4, random_state=SEED, n_init=10)
adata.obs["KMeans_cluster_raw"] = km.fit_predict(adata.obsm["X_pca"]).astype(str)

pc1 = adata.obsm["X_pca"][:, 0]
cluster_pc1_means = (
    pd.DataFrame({"k": adata.obs["KMeans_cluster_raw"].values, "pc1": pc1})
    .groupby("k")["pc1"].mean().sort_values()
)
remap = {raw: f"C{i+1}" for i, raw in enumerate(cluster_pc1_means.index)}
adata.obs["KMeans_cluster"] = adata.obs["KMeans_cluster_raw"].map(remap).astype(
    pd.CategoricalDtype(categories=["C1", "C2", "C3", "C4"], ordered=True)
)
print("[kmeans] cluster sizes:", adata.obs["KMeans_cluster"].value_counts().to_dict())

# ---------- 3. niche-20 signature scoring on polygon AnnData ----------
de = pd.read_csv(DE_PATH, sep=DE_SEP)
panel = set(adata.var_names)
sig_up = [g for g in de[(de["avg_log2FC"] > 0) & (de["p_val_adj"] < 0.05)].sort_values("avg_log2FC", ascending=False)["gene"] if g in panel]
sig_dn = [g for g in de[(de["avg_log2FC"] < 0) & (de["p_val_adj"] < 0.05)].sort_values("avg_log2FC")["gene"] if g in panel]
print(f"[sig] {TAG}: UP on panel: {len(sig_up)}, DOWN on panel: {len(sig_dn)}")
print(f"[sig] DE universe {len(set(de.gene))} genes -> {len(set(de.gene) & panel)} on the 480-gene panel")

# --- how much of this is actually new relative to niche9? computed, not assumed ---
de9 = pd.read_csv(DE_NICHE9, sep="\t")
up9 = set(g for g in de9[(de9.avg_log2FC > 0) & (de9.p_val_adj < 0.05)].gene if g in panel)
dn9 = set(g for g in de9[(de9.avg_log2FC < 0) & (de9.p_val_adj < 0.05)].gene if g in panel)
u20, d20 = set(sig_up), set(sig_dn)
def jac(a, b):
    return len(a & b) / len(a | b) if (a | b) else float("nan")
ov = pd.DataFrame([
    {"direction": "UP", "n_niche20": len(u20), "n_niche9": len(up9), "n_shared": len(u20 & up9),
     "jaccard": round(jac(u20, up9), 3), "only_niche20": ";".join(sorted(u20 - up9)) or "-",
     "only_niche9": ";".join(sorted(up9 - u20)) or "-"},
    {"direction": "DOWN", "n_niche20": len(d20), "n_niche9": len(dn9), "n_shared": len(d20 & dn9),
     "jaccard": round(jac(d20, dn9), 3), "only_niche20": ";".join(sorted(d20 - dn9)) or "-",
     "only_niche9": ";".join(sorted(dn9 - d20)) or "-"},
])
ov.to_csv(OUT / "signature_overlap_vs_niche9.tsv", sep="\t", index=False)
print("\n[overlap] niche20 vs niche9 ON THE GSE294965 PANEL -- if jaccard is ~0.9 this rerun is a\n"
      "          consistency check, NOT independent validation:")
print(ov.to_string(index=False))
pd.Series({"niche20_up": ";".join(sig_up), "niche20_dn": ";".join(sig_dn)}).to_csv(
    OUT / "signature_genes_used.tsv", sep="\t", header=False)

sc.tl.score_genes(adata, gene_list=sig_up, score_name=f"{TAG}_up",
                  random_state=SEED, ctrl_size=min(50, len(sig_up) * 5))
sc.tl.score_genes(adata, gene_list=sig_dn, score_name=f"{TAG}_dn",
                  random_state=SEED, ctrl_size=min(50, len(sig_dn) * 5))
adata.obs[f"{TAG}_net"] = adata.obs[f"{TAG}_up"] - adata.obs[f"{TAG}_dn"]
NET = f"{TAG}_net"

poly_table = adata.obs.copy()
poly_table["PC1"] = adata.obsm["X_pca"][:, 0]
poly_table["PC2"] = adata.obsm["X_pca"][:, 1]
poly_table["PC3"] = adata.obsm["X_pca"][:, 2]
poly_table.to_csv(OUT / "polygon_table.tsv", sep="\t")
print(f"[save] polygon_table: {OUT/'polygon_table.tsv'}")

adata.write_h5ad(OUT / "adata_polygon_reduced_with_score.h5ad")

# ---------- 4. Statistical tests ----------
print(f"\n[stats] correlations of {NET} with paper trajectory metrics")
rs_pt, p_pt = spearmanr(adata.obs[NET], adata.obs["Pseudotime"])
rp_pt, pp_pt = pearsonr(adata.obs[NET], adata.obs["Pseudotime"])
rs_pc1, p_pc1 = spearmanr(adata.obs[NET], poly_table["PC1"])
print(f"  vs Pseudotime: Spearman r={rs_pt:.3f} p={p_pt:.2e}; Pearson r={rp_pt:.3f} p={pp_pt:.2e}")
print(f"  vs PC1:        Spearman r={rs_pc1:.3f} p={p_pc1:.2e}")

groups_km = [adata.obs[NET][adata.obs["KMeans_cluster"] == c].dropna() for c in ["C1", "C2", "C3", "C4"]]
H_km, p_km = kruskal(*groups_km)
print(f"  Kruskal across KMeans (C1-C4): H={H_km:.2f} p={p_km:.2e}")

groups_q = [adata.obs[NET][adata.obs["Quadrant"] == q].dropna() for q in ["Q1", "Q2", "Q3", "Q4"]]
H_q, p_q = kruskal(*groups_q)
print(f"  Kruskal across Pseudotime quadrants: H={H_q:.2f} p={p_q:.2e}")

def pairwise(groups, labels, ref="C1"):
    rows = []
    iref = labels.index(ref)
    for i, lab in enumerate(labels):
        if lab == ref:
            continue
        x = groups[i]
        c = groups[iref]
        if len(x) and len(c):
            U, p = mannwhitneyu(x, c, alternative="greater")
            rows.append({"comparison": f"{lab}_vs_{ref}", "n": len(x), "n_ref": len(c),
                         "median": float(np.median(x)), "median_ref": float(np.median(c)),
                         "U": float(U), "p_one_sided_greater": float(p)})
    return pd.DataFrame(rows)
pw_km = pairwise(groups_km, ["C1", "C2", "C3", "C4"], "C1")
pw_q = pairwise(groups_q, ["Q1", "Q2", "Q3", "Q4"], "Q1")
pw_km.to_csv(OUT / "pairwise_kmeans_vs_C1.tsv", sep="\t", index=False)
pw_q.to_csv(OUT / "pairwise_quadrant_vs_Q1.tsv", sep="\t", index=False)
print("[stats] pairwise vs C1:\n", pw_km)
print("[stats] pairwise vs Q1:\n", pw_q)

ct_km = pd.crosstab(adata.obs["Disease"], adata.obs["KMeans_cluster"])
ct_q = pd.crosstab(adata.obs["Disease"], adata.obs["Quadrant"])
ct_km.to_csv(OUT / "crosstab_disease_kmeans.tsv", sep="\t")
ct_q.to_csv(OUT / "crosstab_disease_quadrant.tsv", sep="\t")
print("[stats] disease x kmeans:\n", ct_km)
print("[stats] disease x quadrant:\n", ct_q)

summary = {
    "signature": TAG,
    "de_source": DE_PATH,
    "n_polygons": int(adata.n_obs),
    "spearman_score_vs_pseudotime": float(rs_pt),
    "spearman_score_vs_pseudotime_p": float(p_pt),
    "pearson_score_vs_pseudotime": float(rp_pt),
    "spearman_score_vs_PC1": float(rs_pc1),
    "kruskal_kmeans_H": float(H_km), "kruskal_kmeans_p": float(p_km),
    "kruskal_quadrant_H": float(H_q), "kruskal_quadrant_p": float(p_q),
    "n_up_panel": len(sig_up), "n_dn_panel": len(sig_dn),
    "jaccard_up_vs_niche9": float(jac(u20, up9)),
    "jaccard_dn_vs_niche9": float(jac(d20, dn9)),
}
pd.Series(summary).to_csv(OUT / "summary.tsv", sep="\t", header=False)

# ---------- 5. Plots ----------
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
fig.savefig(OUT / "pca_overview.pdf", bbox_inches="tight")
fig.savefig(OUT / "pca_overview.png", dpi=200, bbox_inches="tight")
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
fig.savefig(OUT / "score_by_paper_grouping.pdf")
fig.savefig(OUT / "score_by_paper_grouping.png", dpi=200)
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
fig.savefig(OUT / "scatter_pseudotime_vs_niche20.pdf")
fig.savefig(OUT / "scatter_pseudotime_vs_niche20.png", dpi=200)
plt.close(fig)

def heatmap_pct(crosstab, title, fname):
    pct = crosstab.div(crosstab.sum(0), axis=1) * 100
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(pct, annot=True, fmt=".0f", cmap="Reds", ax=ax, cbar_kws={"label": "% of column"})
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=200, bbox_inches="tight")
    fig.savefig(OUT / fname.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
heatmap_pct(ct_q, "% disease per pseudotime quadrant", "heatmap_disease_quadrant.png")
heatmap_pct(ct_km, "% disease per kmeans cluster", "heatmap_disease_kmeans.png")

print(f"\n[done] outputs in {OUT}")
