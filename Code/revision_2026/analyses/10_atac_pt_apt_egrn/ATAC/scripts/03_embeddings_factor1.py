#!/usr/bin/env python3
"""UMAP / PHATE / DiffMap on the GLUE latent space, coloured by cluster and by Factor 1.

Both GLUE versions, side by side, so the cluster separation and a plausible trajectory
direction can be read off by eye.

Factor 1 score = plain X . W  (log-normalised gene scores against the PT-view Factor 1
loadings), matching score_factor_wsum() used for the published Factor-2 comparison.

Species/build: human hg38. Gene IDs: HGNC symbols. Seed: 42.
"""
from kbpaths import KB_ANALYSIS  # repository path config (config/kbpaths.py)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, phate, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

SEED = 42
np.random.seed(SEED)
sc.settings.n_jobs = 8
A    = f"{KB_ANALYSIS}/aki_c5_exclusion_20261022/atac_rerun/ATAC"
LOAD = f"{KB_ANALYSIS}/loadings~hvg~5.csv"
OUT  = f"{A}/figures"; os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{A}/tables", exist_ok=True)

print("versions: scanpy", sc.__version__, "| anndata", ad.__version__,
      "| phate", phate.__version__, "| numpy", np.__version__, flush=True)

L = pd.read_csv(LOAD)
W_pt = L[(L.Factor == "Factor1") & (L.view == "PT")].set_index("variable")["Loading"]
W_pt = W_pt[~W_pt.index.duplicated()]
print(f"Factor1 PT-view loadings: {len(W_pt):,} genes "
      f"({(W_pt>0).sum():,} positive, {(W_pt<0).sum():,} negative)", flush=True)

CLUST_COLORS = {"pt_C1":"#4C72B0","pt_C2":"#DD8452","pt_C3":"#55A868",
                "pt_C4":"#C44E52","pt_C5":"#8172B3","pt_C6":"#937860"}

def embed(a, tag):
    sc.pp.neighbors(a, use_rep="X_glue", n_neighbors=15, random_state=SEED)
    print(f"  [{tag}] neighbors done", flush=True)
    sc.tl.umap(a, random_state=SEED)
    print(f"  [{tag}] umap done", flush=True)
    sc.tl.diffmap(a, n_comps=15)
    print(f"  [{tag}] diffmap done", flush=True)
    op = phate.PHATE(n_components=2, random_state=SEED, n_jobs=8, verbose=0)
    a.obsm["X_phate"] = op.fit_transform(a.obsm["X_glue"])
    print(f"  [{tag}] phate done", flush=True)
    return a

def score_factor1(a):
    b = a.copy()
    b.X = b.X.astype(np.float32)
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    genes = [g for g in W_pt.index if g in b.var_names]
    miss  = len(W_pt) - len(genes)
    print(f"  Factor1 genes present in gene score matrix: {len(genes):,}/{len(W_pt):,} "
          f"(missing {miss:,}, {100*miss/len(W_pt):.1f}%)", flush=True)
    M = b[:, genes].X
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
    s = M @ W_pt.loc[genes].to_numpy(dtype=np.float32)
    del b, M
    return s, len(genes), miss

results = {}
for ver in ["glue", "glue_v3"]:
    print(f"\n══════ {ver} ══════", flush=True)
    a = ad.read_h5ad(f"{A}/objects/gsm_with_{ver}.h5ad")
    print(f"  {a.shape[0]:,} nuclei", flush=True)
    a.obs["Factor1"], n_g, n_miss = score_factor1(a)
    a = embed(a, ver)

    a.obs["cluster"] = a.obs["glue_cluster"].astype(str)
    df = pd.DataFrame({
        "cellid": a.obs_names, "cluster": a.obs["cluster"].values,
        "Factor1": a.obs["Factor1"].values,
        "umap1": a.obsm["X_umap"][:,0],   "umap2": a.obsm["X_umap"][:,1],
        "phate1": a.obsm["X_phate"][:,0], "phate2": a.obsm["X_phate"][:,1],
        "dc1": a.obsm["X_diffmap"][:,1],  "dc2": a.obsm["X_diffmap"][:,2],
        "dc3": a.obsm["X_diffmap"][:,3],
    })
    df.to_csv(f"{A}/tables/embeddings_{ver}.csv", index=False)

    # Factor1 by cluster: mean, 95% CI, n  (effect size + CI, not a bare p-value)
    g = df.groupby("cluster")["Factor1"]
    summ = pd.DataFrame({"n": g.size(), "mean": g.mean(), "sd": g.std()})
    summ["sem"] = summ["sd"]/np.sqrt(summ["n"])
    summ["ci95_lo"] = summ["mean"] - 1.96*summ["sem"]
    summ["ci95_hi"] = summ["mean"] + 1.96*summ["sem"]
    summ = summ.sort_values("mean")
    summ.to_csv(f"{A}/tables/factor1_by_cluster_{ver}.csv")
    print(f"  Factor1 by cluster (sorted):\n{summ.round(3).to_string()}", flush=True)
    results[ver] = dict(df=df, summ=summ, n_genes=n_g, n_missing=n_miss, n_cells=a.shape[0])
    del a

# ── figure: rows = version x colouring, cols = embedding ──────────────────────
EMB = [("umap1","umap2","UMAP"), ("phate1","phate2","PHATE"), ("dc1","dc2","DiffMap DC1/DC2")]
fig, axes = plt.subplots(4, 3, figsize=(16.5, 21))
for r_ver, ver in enumerate(["glue","glue_v3"]):
    df = results[ver]["df"]
    v = np.percentile(np.abs(df.Factor1 - df.Factor1.median()), 99)
    norm = TwoSlopeNorm(vcenter=float(df.Factor1.median()),
                        vmin=float(df.Factor1.median()-v), vmax=float(df.Factor1.median()+v))
    order = np.random.default_rng(SEED).permutation(len(df))
    for c,(x,y,name) in enumerate(EMB):
        # cluster
        ax = axes[2*r_ver, c]
        for cl in sorted(df.cluster.unique()):
            m = df.cluster == cl
            ax.scatter(df.loc[m,x], df.loc[m,y], s=1.6, lw=0, alpha=.55,
                       c=CLUST_COLORS.get(cl,"0.5"), label=f"{cl} ({m.sum():,})", rasterized=True)
        ax.set_title(f"{ver} — {name} — cluster", fontsize=11, fontweight="bold")
        if c == 0:
            leg = ax.legend(markerscale=7, fontsize=8, loc="best", frameon=True)
            for h in leg.legend_handles: h.set_alpha(1)
        # Factor1
        ax2 = axes[2*r_ver+1, c]
        sc_ = ax2.scatter(df[x].values[order], df[y].values[order], s=1.6, lw=0,
                          c=df.Factor1.values[order], cmap="RdBu_r", norm=norm, rasterized=True)
        ax2.set_title(f"{ver} — {name} — Factor 1", fontsize=11, fontweight="bold")
        plt.colorbar(sc_, ax=ax2, fraction=.046, pad=.02)
        for a_ in (ax, ax2):
            a_.set_xticks([]); a_.set_yticks([])
            a_.set_xlabel(x); a_.set_ylabel(y)
fig.suptitle("GLUE latent space — cluster separation and Factor 1, both versions\n"
             f"Factor 1 = log-normalised ATAC gene scores · PT-view Factor 1 loadings "
             f"({results['glue']['n_genes']:,} genes)", fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0,0,1,.975])
fig.savefig(f"{OUT}/glue_embeddings_factor1.pdf", dpi=200)
fig.savefig(f"{OUT}/glue_embeddings_factor1.png", dpi=140)
print(f"\nwrote {OUT}/glue_embeddings_factor1.pdf/.png", flush=True)
print("done", flush=True)
