"""Phase B (Option B): fast label transfer via joint PCA + Harmony + KNN.

Approach: concatenate reference + query, normalize + log1p + PCA jointly,
Harmony-correct by `sample`, then transfer `predic` from reference to
query cells with a weighted KNN classifier on the Harmony-corrected PCA.
Finally recompute neighbours + UMAP on that joint corrected space and
write the 5-region merged h5ad.

Why this over scVI: the 468-gene panel is shared exactly between the
reference and the two new Region 1 regions, and the reference is already
Harmony-corrected in its published form — so the heavy generative
modelling scVI adds to handle feature/batch mismatch is not needed.
Harmony on 261k x 50 PCs + KNN finishes in ~15 min on CPU, versus scVI
which clocked ~30 min/epoch on the full reference (would not fit in a
reasonable SLURM window).
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier

REVIEW_DIR = Path(f"{KB_ANALYSIS}")
REF_H5AD = Path(
    f"{KB_XENIUM_RAW}/h5ad_files_only_transcript/"
    "marker_genes/clustered_annotated_final_file.h5ad"
)
QUERY_H5ADS = [
    REVIEW_DIR / "annotated_0041734_Region_1.h5ad",
    REVIEW_DIR / "annotated_0042348_Region_1.h5ad",
]
OUT_MERGED = Path(
    f"{KB_XENIUM_RAW}/h5ad_files_only_transcript/"
    "marker_genes/clustered_annotated_final_file_5regions.h5ad"
)
QC_DIR = REVIEW_DIR / "phaseB_qc"

LABELS_KEY = "predic"
BATCH_KEY = "sample"


def load_reference() -> ad.AnnData:
    print(f"loading reference {REF_H5AD}", flush=True)
    a = ad.read_h5ad(REF_H5AD)
    n0 = a.n_obs
    sc.pp.filter_cells(a, min_counts=1)
    print(f"  {a.shape} (dropped {n0 - a.n_obs} empty cells)", flush=True)
    print(f"  labels: {sorted(a.obs[LABELS_KEY].astype(str).unique())}", flush=True)
    print(f"  batches: {sorted(a.obs[BATCH_KEY].astype(str).unique())}", flush=True)
    return a


def load_queries() -> ad.AnnData:
    parts = []
    for p in QUERY_H5ADS:
        q = ad.read_h5ad(p)
        print(f"loaded query {p.name}: {q.shape}", flush=True)
        parts.append(q)
    merged = ad.concat(parts, join="outer", axis=0, merge="first")
    n0 = merged.n_obs
    sc.pp.filter_cells(merged, min_counts=1)
    print(f"merged query: {merged.shape} (dropped {n0 - merged.n_obs} empty cells)", flush=True)
    return merged


def prep_joint(ref: ad.AnnData, query: ad.AnnData) -> ad.AnnData:
    common = ref.var_names.intersection(query.var_names)
    assert len(common) == ref.n_vars == query.n_vars, \
        f"gene-set mismatch: ref={ref.n_vars}, query={query.n_vars}, common={len(common)}"
    ref = ref[:, common].copy()
    query = query[:, common].copy()

    # stash raw counts so we can put them back on the merged .X before saving
    ref.layers["counts"] = ref.X.copy()
    query.layers["counts"] = query.X.copy()

    ref.obs["_origin"] = "ref"
    query.obs["_origin"] = "query"

    ref_spatial = np.asarray(ref.obsm["spatial"])
    q_spatial = np.asarray(query.obsm["spatial"])

    joint = ad.concat([ref, query], join="outer", axis=0, merge="first")
    joint.obsm["spatial"] = np.vstack([ref_spatial, q_spatial])
    print(f"joint: {joint.shape}", flush=True)
    return joint


def compute_embedding(joint: ad.AnnData, n_pcs: int) -> str:
    print("\n=== normalize + log1p + PCA ===", flush=True)
    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)
    sc.pp.scale(joint, max_value=10)
    sc.pp.pca(joint, n_comps=n_pcs, zero_center=True)
    print(f"  X_pca: {joint.obsm['X_pca'].shape}", flush=True)

    use_rep = "X_pca"
    try:
        import harmonypy
        print(f"\n=== harmony on '{BATCH_KEY}' ({joint.obs[BATCH_KEY].nunique()} batches) ===", flush=True)
        ho = harmonypy.run_harmony(
            joint.obsm["X_pca"],
            joint.obs,
            vars_use=[BATCH_KEY],
            random_state=0,
        )
        # harmonypy.Z_corr is (n_pcs, n_cells); transpose to (n_cells, n_pcs)
        z = np.asarray(ho.Z_corr)
        if z.shape[0] == joint.n_obs:
            z_corr = z
        elif z.shape[1] == joint.n_obs:
            z_corr = z.T
        else:
            raise RuntimeError(f"harmonypy output shape {z.shape} incompatible with n_obs={joint.n_obs}")
        joint.obsm["X_pca_harmony"] = z_corr
        use_rep = "X_pca_harmony"
        print(f"  X_pca_harmony: {joint.obsm['X_pca_harmony'].shape}", flush=True)
    except Exception as exc:
        print(f"  harmony failed ({exc!r}); falling back to raw X_pca", flush=True)

    return use_rep


def knn_label_transfer(joint: ad.AnnData, use_rep: str, k: int) -> None:
    print(f"\n=== KNN label transfer (k={k}, rep={use_rep}) ===", flush=True)
    ref_mask = (joint.obs["_origin"].astype(str).values == "ref")
    X = joint.obsm[use_rep]

    joint.obs[LABELS_KEY] = joint.obs[LABELS_KEY].astype(str).replace({"nan": "Unassigned"})
    y_ref = joint.obs.loc[ref_mask, LABELS_KEY].astype(str).values

    clf = KNeighborsClassifier(n_neighbors=k, weights="distance", n_jobs=-1)
    clf.fit(X[ref_mask], y_ref)
    preds = clf.predict(X[~ref_mask])
    probs = clf.predict_proba(X[~ref_mask]).max(axis=1)

    # assign back into joint.obs
    joint.obs[LABELS_KEY] = joint.obs[LABELS_KEY].astype(str)
    joint.obs.loc[~ref_mask, LABELS_KEY] = preds
    conf = pd.Series(np.nan, index=joint.obs_names, dtype=float)
    conf.loc[joint.obs_names[~ref_mask]] = probs
    joint.obs["predic_confidence"] = conf

    print(f"  query predictions: {pd.Series(preds).value_counts().to_dict()}", flush=True)
    print(f"  confidence: mean={probs.mean():.2f}, median={np.median(probs):.2f}, "
          f"<0.5: {(probs<0.5).sum()}, <0.3: {(probs<0.3).sum()}", flush=True)


def finalize(joint: ad.AnnData, use_rep: str) -> ad.AnnData:
    print("\n=== neighbours + UMAP on joint corrected space ===", flush=True)
    sc.pp.neighbors(joint, use_rep=use_rep, n_neighbors=15)
    sc.tl.umap(joint)

    # restore raw counts on .X so the merged h5ad matches the reference's convention
    if "counts" in joint.layers:
        joint.X = joint.layers["counts"]
        del joint.layers["counts"]
    # log1p uns flag pollutes downstream log1p calls; drop it
    joint.uns.pop("log1p", None)
    return joint


def plot_qc(merged: ad.AnnData) -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)

    # 1) UMAP by sample and by predic
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    sc.pl.umap(merged, color="sample", ax=axes[0], show=False,
               legend_loc="right margin", size=2, frameon=False)
    sc.pl.umap(merged, color=LABELS_KEY, ax=axes[1], show=False,
               legend_loc="on data", size=2, frameon=False)
    fig.suptitle("5-region merged: UMAP on Harmony-corrected PCA")
    fig.savefig(QC_DIR / "umap_5regions_sample_predic.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # 2) Spatial per sample, coloured by predic
    samples = sorted(merged.obs["sample"].astype(str).unique())
    fig, axes = plt.subplots(1, len(samples), figsize=(6 * len(samples), 6))
    if len(samples) == 1:
        axes = [axes]
    all_cats = sorted(merged.obs[LABELS_KEY].astype(str).unique())
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(all_cats), 1)))
    colour_map = {c: cmap[i % len(cmap)] for i, c in enumerate(all_cats)}
    for ax, s in zip(axes, samples):
        sub = merged[merged.obs["sample"].astype(str) == s]
        xy = np.asarray(sub.obsm["spatial"])
        colours = np.array([colour_map[c] for c in sub.obs[LABELS_KEY].astype(str)])
        ax.scatter(xy[:, 0], xy[:, 1], c=colours, s=0.3, linewidths=0)
        ax.set_aspect("equal")
        ax.set_title(s, fontsize=8)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=colour_map[c], markersize=6, label=c)
        for c in all_cats
    ]
    fig.legend(handles=handles, loc="center right", fontsize=7, frameon=False,
               bbox_to_anchor=(1.02, 0.5))
    fig.suptitle("Spatial per sample, coloured by predicted cell type")
    fig.savefig(QC_DIR / "spatial_per_sample_predic.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # 3) Crosstab Disease x predic for new cells only
    new_mask = merged.obs["sample"].astype(str).isin([
        "XETG00229__0041734__Region_1__20250416__104148",
        "XETG00229__0042348__Region_1__20250416__104148",
    ])
    if "Disease" in merged.obs.columns:
        ct = pd.crosstab(
            merged.obs.loc[new_mask, "Disease"].astype(str),
            merged.obs.loc[new_mask, LABELS_KEY].astype(str),
        )
        ct.to_csv(QC_DIR / "crosstab_new_Disease_x_predic.csv")
        print("\n=== new cells: Disease x predic ===", flush=True)
        print(ct.to_string(), flush=True)

    # 4) predic counts by sample
    cts_by_sample = (
        merged.obs.groupby("sample", observed=True)[LABELS_KEY]
        .value_counts()
        .unstack(fill_value=0)
    )
    cts_by_sample.to_csv(QC_DIR / "predic_counts_by_sample.csv")

    # 5) KNN confidence histogram for new cells
    if "predic_confidence" in merged.obs.columns:
        conf = merged.obs.loc[new_mask, "predic_confidence"].dropna().to_numpy()
        if len(conf):
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.hist(conf, bins=40, edgecolor="black")
            ax.axvline(0.5, color="red", ls="--", lw=1)
            ax.set_xlabel("KNN top-class probability")
            ax.set_ylabel("cells")
            ax.set_title("Label transfer confidence (new cells only)")
            fig.savefig(QC_DIR / "knn_confidence_new_cells.png", dpi=160, bbox_inches="tight")
            plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--knn-k", type=int, default=15)
    ap.add_argument("--smoke", action="store_true",
                    help="subsample ref to 20k + query to 5k for a quick end-to-end check")
    args = ap.parse_args()

    sc.settings.verbosity = 2
    ref = load_reference()
    query = load_queries()

    if args.smoke:
        rng = np.random.default_rng(0)
        ref = ref[rng.choice(ref.n_obs, 20_000, replace=False)].copy()
        query = query[rng.choice(query.n_obs, 5_000, replace=False)].copy()
        print(f"smoke: ref->{ref.shape}, query->{query.shape}", flush=True)

    joint = prep_joint(ref, query)
    use_rep = compute_embedding(joint, n_pcs=args.n_pcs)
    knn_label_transfer(joint, use_rep=use_rep, k=args.knn_k)
    merged = finalize(joint, use_rep=use_rep)

    OUT_MERGED.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nwriting merged: {merged.shape} -> {OUT_MERGED}", flush=True)
    merged.write_h5ad(OUT_MERGED, compression="gzip")

    plot_qc(merged)

    print("\n=== DONE ===", flush=True)
    print("merged shape:", merged.shape, flush=True)
    print("samples:", merged.obs["sample"].value_counts().to_dict(), flush=True)
    print("predic dist (merged):", merged.obs[LABELS_KEY].value_counts().to_dict(), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
