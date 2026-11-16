#!/usr/bin/env python
"""Leiden-resolution sweep on the newly-annotated PEC set, for MANUAL annotation/exclusion.

Input : pec_new_20260717_raw.h5ad  (7,004 cells, level2=="PEC", X=library-normalized, obsm['X_harmony'])

EMBEDDING (the important change): the incoming `X_harmony` is the GLOBAL 30-dim embedding computed
on the full 363,656-cell object with GLOBAL HVGs, merely sliced to PEC. Those components were
optimized to separate whole kidney cell types, so PEC-internal structure is under-resolved -- on that
embedding the PEC UMAP is one continuous blob at every resolution (kept for reference in
../sweep_globalharmony_reference/). Here the embedding is rebuilt ON THE SUBSET:

    log1p -> HVG(2000, batch-aware over patient_ID) -> scale(max 10) -> PCA(30, seed 42)
          -> harmony(batch=patient_ID) -> neighbors(k=15) -> UMAP           [computed ONCE]
          -> leiden at each resolution (igraph, n_iterations=2, seed 42)    [only this varies]

30 harmony dims are used deliberately, matching the global embedding's 30, so the two sweeps differ
only in WHERE the variance came from, not in dimensionality.

Score : decoupler ULM, net = manual-glom niche positive markers (padj<0.05 & log2FC>0)

Because neighbors/UMAP are shared across resolutions, the embedding is IDENTICAL in every subfolder
-- only the partition changes. That is what makes them directly comparable.

Per resolution -> new_pec_annotation/res_<r>/:
  umap_res<r>.{pdf,png}                     UMAP by cluster (+ niche score + provenance for context)
  markers_dotplot_res<r>.{pdf,png}          top-5 Wilcoxon markers per cluster
  markers_res<r>.csv                        full Wilcoxon table (top 100/cluster, logFC/pval/padj/pct)
  ulm_niche_score_res<r>.{pdf,png}          niche ULM score: per-cluster violin
  ulm_niche_score_per_cluster_res<r>.csv    mean/median/sd niche score per cluster
  niche_fixed_top10_bottom10_dotplot_res<r>.* READING A: the niche's own top-10/bottom-10 genes (by
                                             spatial log2FC) shown across clusters
  niche_per_cluster_top10_bottom10_res<r>.csv READING B: per cluster, the 10 highest / 10 lowest
                                             expressed niche genes. NOTE the "lowest expressed" block
                                             is restricted to niche genes detected somewhere
                                             (mean log1p >= 0.10 in >=1 cluster) -- without that filter
                                             it just returns dropout genes (mean ~0.000), identical
                                             across clusters and uninformative.
  niche_gene_mean_expr_res<r>.csv           full niche-gene x cluster mean log1p matrix
  custom_genes_dotplot_res<r>.{pdf,png}     requested panel (SLC4A11/PTPRQ/ADAMTS19/DDN/ITGB3/CFH)
  custom_genes_umap_res<r>.{pdf,png}        same panel as UMAP feature plots, % detected in titles
  custom_genes_mean_expr_res<r>.csv         requested panel: mean log1p per cluster
  custom_genes_pct_detected_res<r>.csv      requested panel: % cells detected per cluster
  cluster_summary_res<r>.csv                n, niche score, provenance mix, top markers

Human, HGNC symbols. Seed 42 throughout. Nothing is excluded here -- exclusion is the user's call.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import os, warnings
import numpy as np, pandas as pd, scanpy as sc, decoupler as dc
import scanpy.external as sce
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 1
sc.settings.autoshow = False
np.random.seed(42)
SEED = 42

DIR  = f"{KB_ANALYSIS}/pec_resubcluster_20260717"
RAW  = f"{DIR}/pec_new_20260717_raw.h5ad"
MARK = f"{KB_ANALYSIS}/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv"
ROOT = f"{DIR}/new_pec_annotation"
RESOLUTIONS = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
BATCH   = "patient_ID"
N_HVG   = 2000
N_PCS   = 30
N_TOP   = 10             # top-N / bottom-N niche genes reported per cluster and in the fixed dotplot
EXPR_MIN = 0.10          # mean log1p floor for "expressed" (same rule as de_trajectory_screen)

# User-requested gene panel (SLC4A11 was listed twice; deduplicated). Detection in this object is
# VERY uneven -- DDN 0.6% of cells, ADAMTS19 1.4%, PTPRQ 7.6% vs ITGB3 54%, CFH 49%, SLC4A11 29%.
# The sparse ones are plotted as asked but will be near-blank rows: that is the data, not a bug.
# Detection rate is stamped onto the figure so it cannot be misread downstream.
CUSTOM_GENES = ["SLC4A11", "PTPRQ", "ADAMTS19", "DDN", "ITGB3", "CFH"]
os.makedirs(ROOT, exist_ok=True)

def save(fig, stem):
    fig.patch.set_facecolor("white")
    for ax in fig.get_axes():
        ax.set_facecolor("white")
    fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(stem + ".png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ---------------------------------------------------------------- load + normalize
adata = sc.read_h5ad(RAW)
print(f"loaded {adata.shape}  obsm={list(adata.obsm)}  X.max={adata.X.max():.1f}")
adata.obsm["X_harmony_global"] = adata.obsm.pop("X_harmony")   # keep, but do NOT cluster on it
adata.var["highly_variable_global"] = adata.var.pop("highly_variable")
adata.layers["libnorm"] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()        # full-gene log1p, used for markers/ULM/dotplots
print("log1p done; X.max=%.2f" % adata.X.max())
print(f"{BATCH}: {adata.obs[BATCH].nunique()} levels; "
      f"cells/patient min={adata.obs[BATCH].value_counts().min()} "
      f"max={adata.obs[BATCH].value_counts().max()}")

# ---------------------------------------------------------------- PEC-specific embedding
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat", batch_key=BATCH)
n_hvg = int(adata.var.highly_variable.sum())
overlap = int((adata.var.highly_variable & adata.var.highly_variable_global).sum())
print(f"\nPEC-subset HVGs: {n_hvg} (overlap with global HVGs: {overlap} = {100*overlap/n_hvg:.0f}%)")

emb = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(emb, max_value=10)
sc.tl.pca(emb, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
vr = emb.uns["pca"]["variance_ratio"]
print(f"PCA {N_PCS} comps; variance explained: PC1={vr[0]:.3f} cum={vr.sum():.3f}")
sce.pp.harmony_integrate(emb, key=BATCH, basis="X_pca", adjusted_basis="X_pca_harmony",
                         random_state=SEED)
adata.obsm["X_pca"] = emb.obsm["X_pca"]
adata.obsm["X_pca_harmony"] = emb.obsm["X_pca_harmony"]
del emb

fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.plot(np.arange(1, N_PCS + 1), vr, "o-", ms=4)
ax.set_xlabel("PC"); ax.set_ylabel("variance ratio")
ax.set_title(f"PEC-subset PCA ({n_hvg} HVGs, {N_PCS} PCs, cum={vr.sum():.2f})")
fig.tight_layout(); save(fig, f"{ROOT}/pca_variance_ratio")

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
print("subset PCA -> harmony -> neighbors -> UMAP done (shared across all resolutions)")

# ---------------------------------------------------------------- niche signature + ULM
m = pd.read_csv(MARK)
sig = m[m.p_val_adj < 0.05].copy()
pos = sig[sig.avg_log2FC > 0]
net = pd.DataFrame({"source": "crescent", "target": pos.gene.values, "weight": 1.0})
present = net.target.isin(adata.var_names)
missing = sorted(net.target[~present].tolist())
print(f"\nniche signature: {len(net)} positive markers, {int(present.sum())} present "
      f"({len(missing)} missing: {missing})")
net = net[present]
dc.run_ulm(mat=adata, net=net, source="source", target="target", weight="weight",
           min_n=3, use_raw=False, verbose=True)
adata.obs["crescent_lesion_score"] = adata.obsm["ulm_estimate"]["crescent"].values
print("crescent_lesion_score: min=%.2f med=%.2f max=%.2f"
      % (adata.obs.crescent_lesion_score.min(), adata.obs.crescent_lesion_score.median(),
         adata.obs.crescent_lesion_score.max()))

sig_present = sig[sig.gene.isin(adata.var_names)].copy()
niche_genes = sig_present.gene.tolist()
fixed_top = sig_present.nlargest(N_TOP, "avg_log2FC").gene.tolist()
fixed_bot = sig_present.nsmallest(N_TOP, "avg_log2FC").gene.tolist()
print(f"\nniche genes present (padj<0.05): {len(niche_genes)}")
print(f"  fixed TOP{N_TOP} (highest log2FC):    {fixed_top}")
print(f"  fixed BOTTOM{N_TOP} (lowest log2FC):  {fixed_bot}")

niche_expr = sc.get.obs_df(adata, keys=niche_genes)
lfc = sig_present.set_index("gene")["avg_log2FC"]
PROV = "provenance" if "provenance" in adata.obs.columns else None

# ---- user-requested custom panel: report detection honestly, never drop silently
custom = [g for g in CUSTOM_GENES if g in adata.var_names]
absent = [g for g in CUSTOM_GENES if g not in adata.var_names]
det = {}
for g in custom:
    x = adata[:, g].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    det[g] = {"mean_log1p": round(float(x.mean()), 4),
              "pct_detected": round(float(100 * (x > 0).mean()), 1),
              "n_cells_detected": int((x > 0).sum())}
det_df = pd.DataFrame(det).T
det_df["niche_log2FC"] = [round(float(lfc[g]), 2) if g in lfc.index else np.nan for g in custom]
det_df["in_niche_DE"] = [g in lfc.index for g in custom]
det_df.to_csv(f"{ROOT}/custom_genes_detection.csv")
print(f"\ncustom gene panel: {len(custom)}/{len(CUSTOM_GENES)} present"
      + (f"; ABSENT from object: {absent}" if absent else "; none absent"))
print(det_df.to_string())
sparse = [g for g in custom if det[g]["pct_detected"] < 10]
if sparse:
    print(f"  !! sparse (<10% of cells detected), will render near-blank: "
          + ", ".join(f"{g} {det[g]['pct_detected']}%" for g in sparse))
det_note = " | ".join(f"{g} {det[g]['pct_detected']}%" for g in custom)

# ---------------------------------------------------------------- sweep
overview = []
for r in RESOLUTIONS:
    key = f"leiden_{r}"
    tag = f"res{r}"
    out = f"{ROOT}/res_{r}"
    os.makedirs(out, exist_ok=True)
    sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    cats = list(adata.obs[key].cat.categories)
    counts = adata.obs[key].value_counts().sort_index()
    print(f"\n{'='*70}\nleiden {r} -> {len(cats)} clusters: {counts.to_dict()}")
    overview.append({"resolution": r, "n_clusters": len(cats),
                     "min_cluster": int(counts.min()), "max_cluster": int(counts.max()),
                     "sizes": ";".join(f"{c}:{n}" for c, n in counts.items())})

    # ---- UMAP
    colors = [key, "crescent_lesion_score"] + ([PROV] if PROV else [])
    fig, axs = plt.subplots(1, len(colors), figsize=(6.2 * len(colors), 5.4))
    axs = np.atleast_1d(axs)
    for ax, c in zip(axs, colors):
        sc.pl.umap(adata, color=c, ax=ax, show=False, frameon=False,
                   legend_loc="on data" if c == key else "right margin",
                   legend_fontsize=9, legend_fontoutline=2,
                   cmap="viridis" if c == "crescent_lesion_score" else None,
                   title={key: f"leiden {r} ({len(cats)} clusters)",
                          "crescent_lesion_score": "niche ULM score"}.get(c, c))
    fig.suptitle(f"PEC (level2==PEC, n={adata.n_obs}) - subset HVG/PCA/harmony - leiden {r}",
                 fontsize=13)
    fig.tight_layout()
    save(fig, f"{out}/umap_{tag}")

    # ---- Wilcoxon markers
    sc.tl.rank_genes_groups(adata, key, method="wilcoxon", n_genes=100, key_added=f"rgg_{key}")
    mk = sc.get.rank_genes_groups_df(adata, group=None, key=f"rgg_{key}")
    mk.to_csv(f"{out}/markers_{tag}.csv", index=False)

    top5 = {c: mk[mk.group == c].sort_values("scores", ascending=False).names.head(5).tolist()
            for c in cats}
    seen, dot_genes = set(), []
    for c in cats:
        for g in top5[c]:
            if g not in seen:
                seen.add(g); dot_genes.append(g)
    dp = sc.pl.dotplot(adata, var_names=top5, groupby=key, standard_scale="var",
                       return_fig=True,
                       figsize=(max(10, 0.30 * len(dot_genes)), 0.55 * len(cats) + 2.2))
    dp.add_totals().style(cmap="Reds")
    dp.make_figure()                       # dp.fig is None until this is called
    fig = dp.fig
    fig.suptitle(f"top-5 Wilcoxon markers per cluster - leiden {r}", fontsize=12, y=1.02)
    save(fig, f"{out}/markers_dotplot_{tag}")

    # ---- ULM niche score
    fig, ax = plt.subplots(figsize=(max(6, 0.85 * len(cats) + 2), 4.6))
    sc.pl.violin(adata, "crescent_lesion_score", groupby=key, ax=ax, show=False,
                 rotation=0, stripplot=False, inner="box")
    ax.axhline(0, ls="--", lw=1, c="grey")
    ax.set_title(f"niche ULM crescent score per cluster - leiden {r}")
    ax.set_xlabel("cluster"); ax.set_ylabel("ULM estimate")
    fig.tight_layout()
    save(fig, f"{out}/ulm_niche_score_{tag}")

    sdf = adata.obs.groupby(key, observed=True)["crescent_lesion_score"]
    pd.DataFrame({"n": sdf.size(), "mean": sdf.mean().round(3),
                  "median": sdf.median().round(3), "sd": sdf.std().round(3)}) \
      .sort_values("mean", ascending=False).to_csv(f"{out}/ulm_niche_score_per_cluster_{tag}.csv")

    # ---- niche expression per cluster
    ne = niche_expr.copy()
    ne[key] = adata.obs[key].values
    mean_niche = ne.groupby(key, observed=True).mean()          # cluster x niche gene
    mean_niche.round(3).T.to_csv(f"{out}/niche_gene_mean_expr_{tag}.csv")

    # READING A: niche's own fixed top-N / bottom-N (by spatial log2FC) across clusters
    dp = sc.pl.dotplot(adata, var_names={f"niche TOP{N_TOP} (up in crescent)": fixed_top,
                                         f"niche BOTTOM{N_TOP} (down in crescent)": fixed_bot},
                       groupby=key, standard_scale="var", return_fig=True,
                       figsize=(0.42 * 2 * N_TOP + 3, 0.55 * len(cats) + 2.2))
    dp.add_totals().style(cmap="RdBu_r")
    dp.make_figure()
    fig = dp.fig
    fig.suptitle(f"niche fixed top-{N_TOP} / bottom-{N_TOP} genes (by spatial log2FC) - leiden {r}",
                 fontsize=11, y=1.03)
    save(fig, f"{out}/niche_fixed_top{N_TOP}_bottom{N_TOP}_dotplot_{tag}")

    # READING B: per cluster, N highest / N lowest EXPRESSED niche genes.
    # The "lowest" block is restricted to genes detected somewhere -- otherwise it returns pure
    # dropout (mean ~0.000), identical in every cluster and uninformative.
    detected = mean_niche.columns[(mean_niche.max(axis=0) >= EXPR_MIN)]
    rows = []
    for c in cats:
        v = mean_niche.loc[c].sort_values(ascending=False)
        for rank, (g, val) in enumerate(v.head(N_TOP).items(), 1):
            rows.append({"cluster": c, "direction": f"top{N_TOP}_expressed", "rank": rank, "gene": g,
                         "mean_log1p": round(float(val), 3),
                         "spatial_log2FC": round(float(lfc[g]), 2)})
        vd = mean_niche.loc[c, detected].sort_values(ascending=False)
        for rank, (g, val) in enumerate(vd.tail(N_TOP)[::-1].items(), 1):
            rows.append({"cluster": c, "direction": f"bottom{N_TOP}_expressed_detectedonly",
                         "rank": rank, "gene": g, "mean_log1p": round(float(val), 3),
                         "spatial_log2FC": round(float(lfc[g]), 2)})
    pd.DataFrame(rows).to_csv(f"{out}/niche_per_cluster_top{N_TOP}_bottom{N_TOP}_{tag}.csv",
                              index=False)

    # ---- user-requested custom gene panel: dotplot + UMAP feature plots + per-cluster means
    dp = sc.pl.dotplot(adata, var_names=custom, groupby=key, return_fig=True,
                       figsize=(0.75 * len(custom) + 3, 0.55 * len(cats) + 2.2))
    dp.add_totals().style(cmap="Reds")
    dp.make_figure()
    fig = dp.fig
    fig.suptitle(f"requested gene panel - leiden {r}\n% of cells detected: {det_note}",
                 fontsize=10, y=1.05)
    save(fig, f"{out}/custom_genes_dotplot_{tag}")

    ncol = 3
    nrow = int(np.ceil(len(custom) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 4.4 * nrow))
    for ax, g in zip(np.atleast_1d(axs).ravel(), custom):
        sc.pl.umap(adata, color=g, ax=ax, show=False, frameon=False, cmap="viridis",
                   title=f"{g}  ({det[g]['pct_detected']}% of cells, n={det[g]['n_cells_detected']})")
    for ax in np.atleast_1d(axs).ravel()[len(custom):]:
        ax.axis("off")
    fig.suptitle(f"requested gene panel on UMAP (log1p) - leiden {r} embedding", fontsize=13)
    fig.tight_layout()
    save(fig, f"{out}/custom_genes_umap_{tag}")

    ce = sc.get.obs_df(adata, keys=custom)
    ce[key] = adata.obs[key].values
    cmean = ce.groupby(key, observed=True).mean().round(3)
    cpct = ce.groupby(key, observed=True).agg(lambda s: round(100 * (s > 0).mean(), 1))
    cmean.T.to_csv(f"{out}/custom_genes_mean_expr_{tag}.csv")
    cpct.T.to_csv(f"{out}/custom_genes_pct_detected_{tag}.csv")

    # ---- cluster summary
    rows = []
    for c in cats:
        sub = adata.obs[adata.obs[key] == c]
        row = {"cluster": c, "n": len(sub),
               "niche_score_mean": round(sub.crescent_lesion_score.mean(), 3),
               "top_markers": ",".join(top5[c])}
        if PROV:
            prov = sub[PROV].value_counts()
            row["top_prov"] = f"{prov.index[0]}({prov.iloc[0]})"
            row["prov_mix"] = ";".join(f"{k}:{v}" for k, v in prov.items() if v > 0)
        rows.append(row)
    summ = pd.DataFrame(rows).sort_values("niche_score_mean", ascending=False)
    summ.to_csv(f"{out}/cluster_summary_{tag}.csv", index=False)
    pd.set_option("display.width", 240, "display.max_colwidth", 70)
    print(summ.to_string(index=False))
    print(f"[saved] {out}")

# ---------------------------------------------------------------- overview
ov = pd.DataFrame(overview)
ov.to_csv(f"{ROOT}/resolution_overview.csv", index=False)
print(f"\n{'='*70}\nRESOLUTION OVERVIEW\n{ov.to_string(index=False)}")

ncol = 5
nrow = int(np.ceil((len(RESOLUTIONS) + 1) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 5.2 * nrow))
axs = axs.ravel()
for ax, r in zip(axs, RESOLUTIONS):
    sc.pl.umap(adata, color=f"leiden_{r}", ax=ax, show=False, frameon=False,
               legend_loc="on data", legend_fontsize=8, legend_fontoutline=2,
               title=f"leiden {r} ({adata.obs[f'leiden_{r}'].nunique()} clusters)")
sc.pl.umap(adata, color="crescent_lesion_score", ax=axs[len(RESOLUTIONS)], show=False, frameon=False,
           cmap="viridis", title="niche ULM score")
for ax in axs[len(RESOLUTIONS) + 1:]:
    ax.axis("off")
fig.suptitle(f"PEC leiden sweep - subset HVG/PCA/harmony({BATCH}) "
             f"(n={adata.n_obs}, shared UMAP, seed {SEED})", fontsize=15)
fig.tight_layout()
save(fig, f"{ROOT}/resolution_sweep_overview")

# side-by-side: what the global embedding hid
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
tmp = adata.copy(); tmp.obsm["X_umap"] = tmp.obsm["X_umap"]
sc.pl.umap(adata, color="leiden_0.8", ax=axs[1], show=False, frameon=False, legend_loc="on data",
           legend_fontsize=9, legend_fontoutline=2, title="subset HVG/PCA/harmony (this run)")
sc.pl.embedding(adata, basis="X_harmony_global", color="leiden_0.8", ax=axs[0], show=False,
                frameon=False, legend_loc="on data", legend_fontsize=9, legend_fontoutline=2,
                title="global harmony dims 1-2 (previous)")
fig.suptitle("same leiden-0.8 labels on the two embeddings", fontsize=13)
fig.tight_layout(); save(fig, f"{ROOT}/embedding_comparison_res0.8")

adata.write_h5ad(f"{ROOT}/pec_leiden_sweep.h5ad")
print(f"\n[done] all outputs under {ROOT}")
print(f"[done] object with all partitions: {ROOT}/pec_leiden_sweep.h5ad")
