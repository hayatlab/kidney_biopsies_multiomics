"""Per-lineage level2 refinement. For each level1 lineage that splits into >1 level2
type, re-embed the lineage's cells, Leiden-subcluster, and assign each subcluster its
CONSENSUS (majority) transferred pred_level2. Non-splitting lineages keep their level1
label. Produces:
  - obs['level2_refined']      : final level2 annotation for every cell (no NaN)
  - obs['subcluster_level2']   : per-cell subcluster id (lineage-prefixed)
  - outputs/metrics/level2_subcluster_assignments.csv
  - outputs/metrics/level2_refined_distribution.csv
  - figures/umap_level2_subclusters.pdf  (per-lineage: subclusters | pred_level2)
Writes the two obs columns back into outputs/xenium_labeled.h5ad (in place)."""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config as C

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 7})

# level1 lineage -> level2 subtypes it is expected to contain
SPLIT = {"PT": ["PT", "aPT", "dPT"], "Stromal": ["FIB", "MES", "VSM-P"],
         "IC": ["IC-A", "IC-B"], "CNT_PC": ["CNT", "PC"]}
# per-lineage Leiden resolution: CNT_PC needs finer resolution to split CNT from PC
# (at res 1.0 the best subcluster is only 49% CNT; at 2.0 two subclusters reach a CNT majority).
RES = {"PT": 1.0, "Stromal": 1.0, "IC": 1.0, "CNT_PC": 2.0}


def subcluster(sub, res):
    sub.X = sub.layers["counts"].copy().astype(np.float32)
    sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
    sc.pp.pca(sub, n_comps=30, random_state=C.SEED)
    sc.pp.neighbors(sub, use_rep="X_pca", random_state=C.SEED)
    sc.tl.leiden(sub, resolution=res, key_added="subcl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=C.SEED)
    sc.tl.umap(sub, random_state=C.SEED)
    return sub


def main():
    C.ensure_dirs()
    a = sc.read_h5ad(C.XENIUM_LABELED)
    obs = a.obs

    # every cell gets a value (no NaN): default = its level1 label (== level2 for non-split lineages)
    level2_refined = obs["pred_level1"].astype(str).copy()
    subcl_id = obs["pred_level1"].astype(str).copy()

    rows, subs = [], {}
    for L, expected in SPLIT.items():
        sub = subcluster(a[(obs["pred_level1"] == L).to_numpy()].copy(), RES[L])
        cons = sub.obs.groupby("subcl", observed=True)["pred_level2"].agg(
            lambda s: s.value_counts().index[0])
        sub.obs["level2_consensus"] = sub.obs["subcl"].map(cons).astype(str)
        idx = sub.obs.index
        level2_refined.loc[idx] = sub.obs["level2_consensus"].values
        subcl_id.loc[idx] = (L + "_" + sub.obs["subcl"].astype(str)).values
        for s_id, g in sub.obs.groupby("subcl", observed=True):
            vc = g["pred_level2"].value_counts()
            rows.append({"lineage": L, "subcluster": f"{L}_{s_id}", "n_cells": int(len(g)),
                         "consensus_level2": vc.index[0], "purity": round(vc.iloc[0] / len(g), 3),
                         "composition": "; ".join(f"{k}:{v}" for k, v in vc.head(4).items() if v > 0)})
        subs[L] = sub
        print(f"{L}: {sub.obs['subcl'].nunique()} subclusters -> consensus level2 = "
              f"{sorted(sub.obs['level2_consensus'].unique())}")

    # --- dPT recovery: hierarchical split of the injured-PT (aPT/dPT) compartment ---
    # dPT is diluted in the full PT lineage (1.6%); focusing on the aPT+dPT compartment
    # (dPT ~10%) lets a confident dPT core form pure subclusters. Same consensus rule.
    inj = a[a.obs["pred_level2"].isin(["aPT", "dPT"]).to_numpy()].copy()
    inj.X = inj.layers["counts"].copy().astype(np.float32)
    sc.pp.normalize_total(inj, target_sum=1e4); sc.pp.log1p(inj)
    sc.pp.pca(inj, n_comps=30, random_state=C.SEED)
    sc.pp.neighbors(inj, use_rep="X_pca", random_state=C.SEED)
    sc.tl.leiden(inj, resolution=1.0, key_added="icl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=C.SEED)
    icons = inj.obs.groupby("icl", observed=True)["pred_level2"].agg(lambda s: s.value_counts().index[0])
    dpt_icl = icons[icons == "dPT"].index
    dpt_cells = inj.obs.index[inj.obs["icl"].isin(dpt_icl)]
    level2_refined.loc[dpt_cells] = "dPT"
    subcl_id.loc[dpt_cells] = ["aPTinj_" + str(c) for c in inj.obs.loc[dpt_cells, "icl"]]
    for s_id, g in inj.obs[inj.obs["icl"].isin(dpt_icl)].groupby("icl", observed=True):
        vc = g["pred_level2"].value_counts()
        rows.append({"lineage": "aPT_injury", "subcluster": f"aPTinj_{s_id}", "n_cells": int(len(g)),
                     "consensus_level2": "dPT", "purity": round(vc.get("dPT", 0) / len(g), 3),
                     "composition": "; ".join(f"{k}:{v}" for k, v in vc.head(4).items() if v > 0)})
    print(f"dPT recovery: {len(dpt_cells)} cells -> dPT from {len(dpt_icl)} dPT-majority injury subclusters")

    a.obs["level2_refined"] = pd.Categorical(level2_refined)
    a.obs["subcluster_level2"] = pd.Categorical(subcl_id)

    pd.DataFrame(rows).to_csv(C.MET / "level2_subcluster_assignments.csv", index=False)
    dist = a.obs["level2_refined"].value_counts().rename_axis("level2_refined").reset_index(name="n_cells")
    dist.to_csv(C.MET / "level2_refined_distribution.csv", index=False)
    print("\nlevel2_refined distribution:\n", dist.to_string(index=False))

    # per-lineage UMAPs: subclusters | pred_level2
    n = len(subs)
    fig, axes = plt.subplots(n, 2, figsize=(9, 3.4 * n))
    for i, (L, sub) in enumerate(subs.items()):
        sc.pl.umap(sub, color="subcl", ax=axes[i, 0], show=False, legend_loc="on data",
                   legend_fontsize=5, size=3, frameon=False, title=f"{L}: subclusters")
        sc.pl.umap(sub, color="pred_level2", ax=axes[i, 1], show=False, legend_fontsize=5,
                   size=3, frameon=False, title=f"{L}: pred_level2")
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(C.FIG / "umap_level2_subclusters.pdf", dpi=200)
    plt.close(fig)

    a.write_h5ad(C.XENIUM_LABELED)
    print("\nwrote level2_refined + subcluster_level2 into", C.XENIUM_LABELED)


if __name__ == "__main__":
    main()
