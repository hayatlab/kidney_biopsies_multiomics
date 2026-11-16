"""Post-process a NicheCompass run with two figures requested by Leon:

1. Niche x cell-type row-normalised heatmap *without* numeric annotations
   (the original `nichecompass_5regions_downstream.py:heatmap_niche_vs_celltype`
   uses `annot=True, fmt=".2f"`; here we drop the numbers).
2. Hallmark pathway dotplot identical in design to the one used for the
   in_glom_niche 9 analysis (`hallmark_dedup_comparison.py`):
   color = MinMax-scaled mean ULM hallmark score per niche, size encodes
   whether the pathway is a top marker for that niche, '*' overlay on the
   accompanying clustermap.

This script is run-dir-agnostic: point it at any NicheCompass run dir
that contains an `adata.h5ad` and pass the niche column key.

Usage
-----
python -u review_analyses/glom_postprocess.py \\
    --run-dir <run_dir> \\
    [--niche-key latent_leiden_0.6] \\
    [--out-subdir glom_postprocess]

Output goes into `<run_dir>/<out-subdir>/`:
  niche_vs_celltype_heatmap_noannot.{pdf,png,csv}
  hallmark_marker_pathways_dedup_heatmap.{pdf,png}
  hallmark_marker_pathways_dotplot.{pdf,png}
  hallmark_score_means.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import anndata as ad
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


CT_KEY_DEFAULT = "predic"
N_HALLMARK_MARKERS = 3


def heatmap_niche_vs_celltype(adata: ad.AnnData,
                              niche_key: str,
                              ct_key: str,
                              out_path_base: Path) -> None:
    """Row-normalised crosstab (niche x cell type), no in-cell numbers."""
    ct = pd.crosstab(
        adata.obs[niche_key].astype(str),
        adata.obs[ct_key].astype(str),
        normalize="index",
    )
    # niche rows ordered by total cell count desc
    niche_totals = adata.obs[niche_key].astype(str).value_counts()
    ct = ct.loc[[n for n in niche_totals.index if n in ct.index]]
    # cell-type columns ordered by global frequency desc
    ct_totals = adata.obs[ct_key].astype(str).value_counts()
    ct = ct[[c for c in ct_totals.index if c in ct.columns]]

    fig, ax = plt.subplots(
        figsize=(max(10, 0.55 * ct.shape[1] + 2), 0.55 * ct.shape[0] + 2.5)
    )
    sns.heatmap(
        ct, cmap="hot_r", annot=False,
        cbar_kws={"label": "row prop"},
        linewidths=0.3, linecolor="white", ax=ax,
    )
    ax.set_xlabel(ct_key)
    ax.set_ylabel(f"Niche ({niche_key})")
    ax.set_title(f"Niche x {ct_key} (row-normalised)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = f"{out_path_base}.{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=200 if ext == "png" else None)
        print(f"  wrote {out}", flush=True)
    plt.close(fig)
    ct.round(4).to_csv(f"{out_path_base}.csv")
    print(f"  wrote {out_path_base}.csv", flush=True)


def hallmark_pathway_dotplot(adata: ad.AnnData,
                             niche_key: str,
                             out_dir: Path,
                             n_markers: int = N_HALLMARK_MARKERS, prefix: str = "") -> None:
    """Mirror of `hallmark_dedup_comparison.py` (dedup heatmap + dotplot). `prefix` is prepended to every
    output file name (figure entry scripts); default "" = the original names."""
    print("[hallmark] fetching gene sets", flush=True)
    hallmark = dc.op.hallmark(organism="human")
    print(f"[hallmark] {hallmark.shape[0]} (source,target) pairs, "
          f"{hallmark['source'].nunique()} pathways", flush=True)

    # decoupler needs .X on a log-normalised scale. adata.X on the saved
    # NicheCompass objects is RAW COUNTS, so score it the same way
    # nichecompass_5regions_downstream.pathway_activity() does — otherwise the
    # per-niche ULM ranking differs and the dedup heatmap silently loses
    # pathways the "normal" matrixplot shows (e.g. TNFA_SIGNALING_VIA_NFKB).
    data = adata.copy()
    if "counts" in data.layers:
        data.X = data.layers["counts"]
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

    # ULM scoring writes to data.obsm['score_ulm']
    print("[hallmark] running ULM (log-normalised)", flush=True)
    dc.mt.ulm(data=data, net=hallmark)
    score = dc.pp.get_obsm(adata=data, key="score_ulm")
    print(f"[hallmark] score adata: {score.shape}", flush=True)

    print("[hallmark] ranking markers per niche", flush=True)
    df = dc.tl.rankby_group(
        adata=score,
        groupby=niche_key,
        reference="rest",
        method="t-test_overestim_var",
    )
    df = df[df["stat"] > 0]
    # sort explicitly rather than trusting rankby_group's row order
    top_markers_per_niche = (
        df.sort_values(["group", "stat"], ascending=[True, False])
          .groupby("group", observed=True).head(n_markers)
          .groupby("group", observed=True)["name"].apply(list).to_dict()
    )
    all_pathways = list(dict.fromkeys(
        p for ps in top_markers_per_niche.values() for p in ps
    ))
    print(f"[hallmark] unique top-{n_markers} pathways: {len(all_pathways)}",
          flush=True)
    if not all_pathways:
        print("[hallmark] WARN: no positive markers; aborting dotplot", flush=True)
        return

    # per-niche mean ULM score (pathway x niche)
    X = score.X.toarray() if hasattr(score.X, "toarray") else np.asarray(score.X)
    score_df = pd.DataFrame(X, index=score.obs_names, columns=score.var_names)
    score_df["_niche"] = score.obs[niche_key].astype(str).values
    mat = score_df.groupby("_niche")[all_pathways].mean().T

    # order niche columns by integer if possible
    try:
        niche_order = sorted(mat.columns, key=lambda v: int(v))
    except ValueError:
        niche_order = sorted(mat.columns)
    mat = mat[niche_order]

    # min-max scale per pathway
    denom = (mat.max(axis=1) - mat.min(axis=1)).replace(0, 1)
    mat_z = mat.sub(mat.min(axis=1), axis=0).div(denom, axis=0)

    # top-marker boolean mask
    is_top = pd.DataFrame(False, index=all_pathways, columns=mat.columns)
    for niche, paths in top_markers_per_niche.items():
        n = str(niche)
        if n in is_top.columns:
            for p in paths:
                if p in is_top.index:
                    is_top.loc[p, n] = True

    mat.round(4).to_csv(out_dir / f"{prefix}hallmark_score_means.csv")
    print(f"  wrote {out_dir / (prefix + 'hallmark_score_means.csv')}", flush=True)

    # --- dedup clustermap, saved twice: without and with the '*' overlay ---
    # Both variants come from the SAME clustermap object, so the clustering,
    # row/col order, colours and figure size are byte-for-byte comparable and the
    # no-stars version is the starred one with the annotation layer removed.
    print("[hallmark] plotting dedup clustermap", flush=True)
    g = sns.clustermap(
        mat_z, cmap="RdBu_r",
        figsize=(0.55 * mat_z.shape[1] + 4, 0.32 * mat_z.shape[0] + 2),
        cbar_kws={"label": "Min-max scaled score"},
        linewidths=0.4, linecolor="white",
        dendrogram_ratio=(0.12, 0.12),
    )
    row_order = [mat_z.index[i] for i in g.dendrogram_row.reordered_ind]
    col_order = [mat_z.columns[i] for i in g.dendrogram_col.reordered_ind]
    ax = g.ax_heatmap
    ax.set_xlabel(f"Niche ({niche_key})")
    ax.set_ylabel("")

    g.fig.suptitle("Hallmark pathways (deduplicated)", y=1.02)
    g.fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dedup_heatmap_nostars.pdf",
                  bbox_inches="tight")
    g.fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dedup_heatmap_nostars.png",
                  bbox_inches="tight", dpi=200)
    print(f"  wrote {out_dir / (prefix + 'hallmark_marker_pathways_dedup_heatmap_nostars.pdf')}",
          flush=True)

    for r, p in enumerate(row_order):
        for c, n in enumerate(col_order):
            if is_top.loc[p, n]:
                ax.text(c + 0.5, r + 0.5, "*", ha="center", va="center",
                        color="black", fontsize=11, fontweight="bold")
    g.fig.suptitle(
        "Hallmark pathways (deduplicated)   *  = top marker for that niche",
        y=1.02,
    )
    g.fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dedup_heatmap.pdf",
                  bbox_inches="tight")
    g.fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dedup_heatmap.png",
                  bbox_inches="tight", dpi=200)
    plt.close(g.fig)
    print(f"  wrote {out_dir / (prefix + 'hallmark_marker_pathways_dedup_heatmap.pdf')}",
          flush=True)

    # --- dotplot, same row/col order as the clustermap ---
    print("[hallmark] plotting dotplot", flush=True)
    fig, ax = plt.subplots(
        figsize=(0.55 * len(col_order) + 4, 0.32 * len(row_order) + 2)
    )
    xs, ys, cs, ss = [], [], [], []
    for r, p in enumerate(row_order):
        for c, n in enumerate(col_order):
            xs.append(c); ys.append(r)
            cs.append(mat_z.loc[p, n])
            ss.append(280 if is_top.loc[p, n] else 35)
    sc_h = ax.scatter(xs, ys, c=cs, s=ss, cmap="RdBu_r",
                      edgecolors="black", linewidths=0.5, vmin=0, vmax=1)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)
    ax.set_xlim(-0.5, len(col_order) - 0.5)
    ax.set_ylim(len(row_order) - 0.5, -0.5)
    ax.set_xlabel(f"Niche ({niche_key})")
    ax.grid(True, linestyle=":", alpha=0.3)
    plt.colorbar(sc_h, ax=ax, label="Min-max scaled score", shrink=0.5)
    for s, lbl in [(35, "background"), (280, "top marker")]:
        ax.scatter([], [], s=s, c="lightgrey", edgecolors="black", label=lbl)
    ax.legend(title="", loc="upper left", bbox_to_anchor=(1.22, 0.45),
              frameon=False)
    ax.set_title("Hallmark pathways - dot plot "
                 "(color = score, size = top marker)")
    fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dotplot.pdf",
                bbox_inches="tight")
    fig.savefig(out_dir / f"{prefix}hallmark_marker_pathways_dotplot.png",
                bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  wrote {out_dir / (prefix + 'hallmark_marker_pathways_dotplot.pdf')}",
          flush=True)

    # traceability: which pathways were the top-N per niche
    pd.DataFrame(
        [{"niche": n, "rank": i + 1, "pathway": p}
         for n, ps in top_markers_per_niche.items() for i, p in enumerate(ps)]
    ).to_csv(out_dir / f"{prefix}hallmark_top{n_markers}_per_niche.csv", index=False)
    return top_markers_per_niche


def niche_vs_in_glom(adata: ad.AnnData,
                     niche_key: str,
                     out_csv: Path) -> None:
    """For runs where in_glom is present in obs, print + save which niche
    is most enriched for in_glom — the analog of the original niche 9."""
    if "in_glom" not in adata.obs.columns:
        print("[niche_vs_in_glom] in_glom not in adata.obs — skipping", flush=True)
        return
    ct = (
        adata.obs.groupby(niche_key, observed=True)["in_glom"]
        .agg(n_total="count", n_in_glom="sum")
    )
    ct["pct_in_glom"] = 100.0 * ct["n_in_glom"] / ct["n_total"]
    total_in_glom = int(ct["n_in_glom"].sum())
    ct["share_of_in_glom"] = 100.0 * ct["n_in_glom"] / max(total_in_glom, 1)
    ct = ct.sort_values("pct_in_glom", ascending=False)
    ct.to_csv(out_csv)
    print(f"[niche_vs_in_glom] saved {out_csv}", flush=True)
    print(ct.to_string(), flush=True)
    if not ct.empty:
        top = ct.index[0]
        print(f"[niche_vs_in_glom] top niche by pct_in_glom: {top!r} "
              f"({ct.loc[top, 'pct_in_glom']:.1f}% of its cells, "
              f"{ct.loc[top, 'share_of_in_glom']:.1f}% of all in_glom)",
              flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True,
                    help="NicheCompass run directory containing adata.h5ad")
    ap.add_argument("--niche-key", default=None,
                    help="obs column with niche labels. Auto-detects "
                         "latent_leiden_{0.4,0.6} if not given.")
    ap.add_argument("--ct-key", default=CT_KEY_DEFAULT)
    ap.add_argument("--out-subdir", default="glom_postprocess",
                    help="subdirectory under run-dir for outputs")
    ap.add_argument("--n-markers", type=int, default=N_HALLMARK_MARKERS,
                    help="top-N hallmark markers per niche for dotplot")
    args = ap.parse_args()

    adata_path = args.run_dir / "adata.h5ad"
    if not adata_path.exists():
        print(f"[error] adata not found: {adata_path}", flush=True)
        return 1

    out_dir = args.run_dir / args.out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run] run_dir={args.run_dir}", flush=True)
    print(f"[run] out_dir={out_dir}", flush=True)

    print(f"[load] reading {adata_path}", flush=True)
    adata = ad.read_h5ad(adata_path)
    print(f"[load] adata: {adata.shape}", flush=True)

    # auto-detect niche key if not provided
    niche_key = args.niche_key
    if niche_key is None:
        for candidate in ("latent_leiden_0.6", "latent_leiden_0.4"):
            if candidate in adata.obs.columns:
                niche_key = candidate
                break
        if niche_key is None:
            avail = [c for c in adata.obs.columns if c.startswith("latent_leiden")]
            raise KeyError(f"no niche key found. Available: {avail}")
    print(f"[run] niche_key={niche_key}", flush=True)
    print(f"[run] ct_key={args.ct_key}", flush=True)

    if args.ct_key not in adata.obs.columns:
        raise KeyError(f"{args.ct_key} not in adata.obs")

    # 1. heatmap (no numbers)
    print("\n[fig] niche vs cell-type heatmap (no annotations)", flush=True)
    heatmap_niche_vs_celltype(
        adata, niche_key=niche_key, ct_key=args.ct_key,
        out_path_base=out_dir / "niche_vs_celltype_heatmap_noannot",
    )

    # 2. niche vs in_glom diagnostic
    print("\n[fig] niche x in_glom enrichment", flush=True)
    niche_vs_in_glom(
        adata, niche_key=niche_key,
        out_csv=out_dir / "niche_vs_in_glom.csv",
    )

    # 3. hallmark pathway dotplot
    print("\n[fig] hallmark pathway dotplot", flush=True)
    hallmark_pathway_dotplot(
        adata, niche_key=niche_key, out_dir=out_dir,
        n_markers=args.n_markers,
    )

    print(f"\n=== DONE outputs in {out_dir} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
