"""Downstream NicheCompass analyses on the trained 5-region model.

Reproduces every analysis from
`review_analyses/nichecompass_filtered_tma.ipynb` and adds a few:

 1. Re-computes niche Leiden (res 0.4) on the latent graph.
 2. Niche composition plots (cell type / condition / sample) with
    absolute counts + proportions side-by-side (already in the main
    script but rerun here so the enriched adata lines up).
 3. Niche x cell-type row-normalized crosstab heatmap.
 4. Niche x enriched-GP mean-activity heatmap (MinMax scaled, viridis).
 5. MSigDB Hallmark pathway activity per niche (decoupler ULM +
    rankby_group, top-3 per niche, matrixplot).
 6. MSigDB Reactome pathway activity per niche (same machinery).
 7. Per-TMA spatial niche grid (3 cols).
 8. Condition x niche binomial + Fisher enrichment; signed -log10(FDR)
    heatmap + clustermap; enriched-pair CSVs.
 9. Supplemental patient metadata merge (TMA_ID -> sample_ID).
10. IFTA_group and EGFR_group Fisher enrichment + bubble / effect-size /
    ECDF / ridgeline / lollipop plots (the cell-88 framework).
11. ROI per-cell flag ingestion from the three `*_roi.csv` files ->
    per-niche %True heatmap + per-flag barplots.
12. Per-niche within-niche DE (cell 98): `rank_genes_groups` on `predic`
    inside each niche-subset; CSVs + ranked / heatmap / dotplot.
13. Writes the fully enriched adata and a concise niche-annotation
    dict (numeric niche -> dominant cell type + top enriched conditions
    + top hallmark pathway).

This is run AFTER `nichecompass_5regions.py` finishes training.
"""
from __future__ import annotations
from kbpaths import KB_ANALYSIS, KB_PROJECT  # repository path config (config/kbpaths.py)

import argparse
import sys
import traceback
import warnings
from pathlib import Path
from typing import Iterable

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import binomtest, fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# -----------------------------------------------------------------------
PROJECT_DIR = Path(f"{KB_PROJECT}")
ARTIFACTS_DIR = PROJECT_DIR / "niche_compass_results/artifacts/single_sample"
ROI_DIR = PROJECT_DIR / "review_analyses/regions_of_interest"
SUPPL_CSV = PROJECT_DIR / "niche_compass_results/Supplemental.csv"

CELL_TYPE_KEY = "predic"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
TMA_KEY = "TMA_ID"
SPATIAL_KEY = "spatial"
LATENT_KEY = "nichecompass_latent"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
DIFF_KEY = "nichecompass_differential_gp_test_results"
LEIDEN_RES = 0.4
NICHE_KEY = f"latent_leiden_{LEIDEN_RES}"
LOG_BF = 2.3

# ROI files — cell 90 pattern
ROI_FILES = [
    ROI_DIR / "1734__Region_3_roi.csv",
    ROI_DIR / "42348__Region_2_roi.csv",
    ROI_DIR / "42348__Region_3_roi.csv",
]


# ======================================================================
# 0. load model + recompute niches
# ======================================================================
def load_trained(run_dir: Path):
    from nichecompass.models import NicheCompass
    model_dir = run_dir / "model"
    print(f"[load] {model_dir}", flush=True)
    model = NicheCompass.load(
        dir_path=str(model_dir),
        adata=None,
        adata_file_name="adata.h5ad",
        gp_names_key=GP_NAMES_KEY,
    )
    adata = model.adata
    print(f"[load] adata={adata.shape}", flush=True)

    if NICHE_KEY not in adata.obs.columns:
        print(f"[load] recomputing neighbors+umap+leiden ({LEIDEN_RES}) on latent", flush=True)
        sc.pp.neighbors(adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)
        sc.tl.umap(adata, neighbors_key=LATENT_KEY)
        sc.tl.leiden(adata=adata, resolution=LEIDEN_RES,
                     key_added=NICHE_KEY, neighbors_key=LATENT_KEY)
    print(f"[load] niches: {adata.obs[NICHE_KEY].value_counts().to_dict()}", flush=True)
    return model, adata


# ======================================================================
# 1. composition plots with abs counts + proportions
# ======================================================================
def plot_niche_composition(adata, niche_key: str, cat_key: str,
                           out_path_base: Path, title_suffix: str = "",
                           colour_override: dict | None = None) -> None:
    df = adata.obs[[niche_key, cat_key]].astype(str).copy()
    counts = (df.groupby([niche_key, cat_key]).size().unstack(fill_value=0))
    niche_totals = counts.sum(axis=1).sort_values(ascending=False)
    counts = counts.loc[niche_totals.index]
    cat_order = counts.sum(axis=0).sort_values(ascending=False).index
    counts = counts[cat_order]
    props = counts.div(counts.sum(axis=1), axis=0)

    if colour_override:
        colour_map = {c: colour_override.get(c, "#999999") for c in cat_order}
    else:
        palette = plt.cm.tab20(np.linspace(0, 1, max(len(cat_order), 1)))
        colour_map = {c: palette[i % len(palette)] for i, c in enumerate(cat_order)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
    bottom = np.zeros(len(counts.index))
    xpos = np.arange(len(counts.index))
    for cat in cat_order:
        axes[0].bar(xpos, counts[cat].values, bottom=bottom, color=colour_map[cat],
                    label=cat, width=0.78, edgecolor="white", linewidth=0.4)
        bottom += counts[cat].values
    totals = counts.sum(axis=1).values
    tmax = totals.max() if len(totals) else 1
    for x, t in zip(xpos, totals):
        axes[0].text(x, t + tmax * 0.015, f"n={int(t):,}",
                     ha="center", va="bottom", fontsize=8)
    axes[0].set_ylim(0, tmax * 1.12)
    axes[0].set_ylabel("Cell count")
    axes[0].set_title(f"Absolute counts per niche{title_suffix}")

    bottom = np.zeros(len(props.index))
    for cat in cat_order:
        axes[1].bar(xpos, props[cat].values, bottom=bottom, color=colour_map[cat],
                    width=0.78, edgecolor="white", linewidth=0.4)
        bottom += props[cat].values
    axes[1].set_ylim(0, 1.0)
    axes[1].set_ylabel("Proportion")
    axes[1].set_title(f"Proportions per niche{title_suffix}")
    for x, t in zip(xpos, totals):
        axes[1].text(x, -0.035, f"n={int(t):,}", ha="center", va="top",
                     fontsize=7, color="#555")

    for ax in axes:
        ax.set_xticks(xpos)
        ax.set_xticklabels([str(n) for n in counts.index])
        ax.set_xlabel(f"Niche ({niche_key})")
        ax.grid(axis="y", linestyle=":", alpha=0.4)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colour_map[c]) for c in cat_order]
    fig.legend(handles, list(cat_order), title=cat_key,
               loc="center left", bbox_to_anchor=(1.01, 0.5),
               frameon=False, fontsize=8)
    fig.suptitle(f"Niche composition by {cat_key}", fontsize=14, y=1.02)
    fig.tight_layout()
    base = str(out_path_base)
    for ext in ("svg", "png"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)
    counts.to_csv(f"{base}.csv")
    props.round(4).to_csv(f"{base}_proportions.csv")
    print(f"  wrote {base}.{{svg,png,csv}}", flush=True)


# ======================================================================
# 3. niche x cell-type row-normalized heatmap
# ======================================================================
def heatmap_niche_vs_celltype(adata, fig_dir: Path, out_name: str = "niche_vs_celltype_heatmap") -> None:
    ct = pd.crosstab(adata.obs[NICHE_KEY].astype(str),
                     adata.obs[CELL_TYPE_KEY].astype(str), normalize="index")
    fig, ax = plt.subplots(figsize=(max(10, 0.5 * ct.shape[1]), 0.6 * ct.shape[0] + 2))
    sns.heatmap(ct, cmap="hot_r", annot=True, fmt=".2f", cbar_kws={"label": "row prop"},
                linewidths=0.3, linecolor="white", ax=ax)
    ax.set_xlabel(CELL_TYPE_KEY)
    ax.set_ylabel("niche")
    ax.set_title("Niche x cell-type (row-normalized)")
    fig.tight_layout()
    out = fig_dir / out_name
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)
    ct.round(4).to_csv(f"{out}.csv")
    print(f"  wrote {out}.{{pdf,png,csv}}", flush=True)


# ======================================================================
# 4. GP activity heatmap (niche x enriched GPs)
# ======================================================================
def gp_activity_heatmap(model, adata, fig_dir: Path) -> None:
    try:
        enriched = model.run_differential_gp_tests(
            cat_key=NICHE_KEY,
            selected_cats=None,
            comparison_cats="rest",
            log_bayes_factor_thresh=LOG_BF,
            key_added=DIFF_KEY,
        )
    except Exception as exc:
        print(f"  [diff_gp] failed: {exc!r}", flush=True)
        traceback.print_exc()
        return
    print(f"  enriched GPs: {len(enriched)}", flush=True)

    # latent obsm columns align with ACTIVE_GP_NAMES_KEY (active_gp_thresh_ratio filtering)
    active_names = list(adata.uns[ACTIVE_GP_NAMES_KEY]) if ACTIVE_GP_NAMES_KEY in adata.uns else list(adata.uns[GP_NAMES_KEY])
    latent = np.asarray(adata.obsm[LATENT_KEY])
    if latent.shape[1] != len(active_names):
        # fall back: assume gp_names ordering and truncate to latent width
        print(f"  [gp_heatmap] latent width {latent.shape[1]} != active_names {len(active_names)}; "
              f"falling back to first {latent.shape[1]} gp_names", flush=True)
        all_names = list(adata.uns[GP_NAMES_KEY])
        active_names = all_names[:latent.shape[1]]
    name_to_idx = {n: i for i, n in enumerate(active_names)}
    cols_present = [g for g in enriched if g in name_to_idx]
    missing = [g for g in enriched if g not in name_to_idx]
    if missing:
        print(f"  [gp_heatmap] {len(missing)} enriched GPs not in active list (skipped): {missing[:5]}…", flush=True)
    if not cols_present:
        print("  [gp_heatmap] no enriched GP columns found in latent", flush=True)
        return
    idx = [name_to_idx[g] for g in cols_present]
    sub = latent[:, idx]
    niches = adata.obs[NICHE_KEY].astype(str).values
    mean_df = pd.DataFrame(sub, columns=cols_present).groupby(niches).mean()

    # MinMax scale per GP column
    mn = mean_df.min(axis=0)
    mx = mean_df.max(axis=0)
    rng = (mx - mn).replace(0, 1)
    scaled = (mean_df - mn) / rng

    fig, ax = plt.subplots(figsize=(min(24, 0.22 * scaled.shape[1] + 4),
                                    max(4, 0.4 * scaled.shape[0] + 2)))
    sns.heatmap(scaled, cmap="viridis", cbar_kws={"label": "MinMax scaled"},
                linewidths=0, ax=ax)
    ax.set_xlabel("Enriched GPs")
    ax.set_ylabel("niche")
    ax.set_title(f"Niche x enriched-GP activity (MinMax scaled; log_BF>{LOG_BF}, n={len(cols_present)})")
    fig.tight_layout()
    out = fig_dir / "enriched_gps_heatmap"
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}.png", bbox_inches="tight", dpi=140)
    plt.close(fig)
    mean_df.to_csv(f"{out}_means.csv")
    scaled.to_csv(f"{out}_minmax.csv")
    print(f"  wrote {out}.{{pdf,png}} + means / minmax csvs", flush=True)

    # enriched-GP summary
    try:
        gp_summary = model.get_gp_summary()
        ordered = gp_summary[gp_summary["gp_name"].isin(enriched)].copy()
        ordered["gp_name"] = pd.Categorical(ordered["gp_name"],
                                             categories=enriched, ordered=True)
        ordered = ordered.sort_values("gp_name")
        ordered.to_csv(fig_dir / f"log_bayes_factor_{LOG_BF}_niche_enriched_gps_summary.csv",
                       index=False)
    except Exception as exc:
        print(f"  [gp_summary] skipped: {exc!r}", flush=True)


# ======================================================================
# 5/6. Hallmark + Reactome pathway activity (decoupler ULM)
# ======================================================================
def pathway_activity(adata, fig_dir: Path, collection_fn, tag: str) -> None:
    try:
        import decoupler as dc
    except Exception as exc:
        print(f"  [{tag}] decoupler not installed: {exc!r}", flush=True)
        return

    print(f"\n[{tag}] loading collection", flush=True)
    net = collection_fn(dc)
    print(f"  net rows: {len(net)}  unique pathways: {net['source'].nunique()}", flush=True)

    # decoupler needs .X in a sensible scale; use counts layer if present,
    # else adata.X (already raw counts on our merged file).
    data = adata.copy()
    if "counts" in data.layers:
        data.X = data.layers["counts"]
    # log-normalise
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

    try:
        dc.mt.ulm(data=data, net=net, verbose=False)
    except Exception as exc:
        print(f"  [{tag}] dc.mt.ulm failed: {exc!r}", flush=True)
        traceback.print_exc()
        return

    score = dc.pp.get_obsm(data, key="score_ulm")
    # rank pathways per niche
    try:
        ranked = dc.tl.rankby_group(
            score, groupby=NICHE_KEY, reference="rest",
            method="t-test_overestim_var",
        )
    except Exception as exc:
        print(f"  [{tag}] rankby_group failed: {exc!r}", flush=True)
        return
    ranked = ranked[ranked["stat"] > 0].copy()
    ranked.to_csv(fig_dir / f"{tag}_rankby_niche_full.csv", index=False)

    # top-3 per niche
    source_markers = (
        ranked.sort_values(["group", "stat"], ascending=[True, False])
        .groupby("group", observed=True)
        .head(3)
    )
    top_dict = {str(g): df["name"].tolist()
                for g, df in source_markers.groupby("group", observed=True)}
    print(f"  [{tag}] top-3 per niche:")
    for g, lst in top_dict.items():
        print(f"    niche {g}: {lst}")

    # matrixplot
    try:
        score.obs[NICHE_KEY] = data.obs[NICHE_KEY].values
        mplot = sc.pl.matrixplot(
            score, var_names=top_dict, groupby=NICHE_KEY,
            dendrogram=True, standard_scale="var",
            colorbar_title="Z-scaled ULM score",
            cmap="RdBu_r", swap_axes=True, return_fig=True,
        )
        out = fig_dir / f"{tag}_marker_pathways"
        mplot.savefig(f"{out}.pdf", bbox_inches="tight")
        mplot.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
        plt.close("all")
        print(f"  wrote {out}.{{pdf,png}}", flush=True)
    except Exception as exc:
        print(f"  [{tag}] matrixplot failed: {exc!r}", flush=True)
        traceback.print_exc()


def hallmark_fn(dc):
    return dc.op.hallmark(organism="human")


def reactome_fn(dc):
    markers = dc.op.resource(name="MSigDB")
    markers = markers[markers["collection"] == "reactome_pathways"].copy()
    markers = markers.drop_duplicates(subset=["geneset", "genesymbol"])
    return markers[["geneset", "genesymbol"]].rename(
        columns={"geneset": "source", "genesymbol": "target"}
    )


# ======================================================================
# 7. per-TMA spatial grid
# ======================================================================
def spatial_niche_grid(adata, fig_dir: Path) -> None:
    samples = sorted(adata.obs[SAMPLE_KEY].astype(str).unique())
    n_cols = min(3, len(samples))
    n_rows = int(np.ceil(len(samples) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.5 * n_cols, 5.5 * n_rows),
                             squeeze=False)
    # stable palette from existing uns if available else tab20
    cats = sorted(adata.obs[NICHE_KEY].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else x)
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(cats), 1)))
    colour_map = {c: palette[i % len(palette)] for i, c in enumerate(cats)}
    for idx, s in enumerate(samples):
        r, c = divmod(idx, n_cols)
        ax = axes[r][c]
        sub = adata[adata.obs[SAMPLE_KEY].astype(str) == s]
        xy = np.asarray(sub.obsm[SPATIAL_KEY])
        colours = np.array([colour_map[p] for p in sub.obs[NICHE_KEY].astype(str)])
        ax.scatter(xy[:, 0], xy[:, 1], c=colours, s=0.3, linewidths=0)
        ax.set_aspect("equal")
        ax.set_title(s.replace("XETG00229__", ""), fontsize=9)
    # hide empty axes
    for idx in range(len(samples), n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r][c].axis("off")
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=colour_map[c], markersize=6, label=c)
               for c in cats]
    fig.legend(handles=handles, loc="center right", fontsize=8,
               bbox_to_anchor=(1.02, 0.5), title="niche", frameon=False)
    fig.suptitle("Niches per TMA (all samples)", fontsize=14)
    fig.tight_layout()
    out = fig_dir / "niches_per_TMA_ID"
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out}.{{pdf,png}}", flush=True)


# ======================================================================
# 8. condition enrichment (binomial + Fisher)
# ======================================================================
def run_categorical_enrichment(
    adata, cat_key: str, fig_dir: Path, prefix: str,
    min_niche_cells: int = 10,
) -> pd.DataFrame:
    """Niche x cat_key binomial + Fisher enrichment, with BH FDR.

    Returns the Fisher-test DataFrame (for downstream plots). Writes:
      - {prefix}_counts_by_niche.csv
      - {prefix}_proportions_by_niche.csv
      - {prefix}_enrichment_binomial.csv
      - {prefix}_enrichment_fisher.csv
      - {prefix}_proportions_by_niche_heatmap.{pdf,png}
      - {prefix}_enrichment_signedlog10fdr_heatmap.{pdf,png}
      - {prefix}_enrichment_signedlog10fdr_clustermap.{pdf,png}
    """
    obs = adata.obs[[NICHE_KEY, cat_key]].astype(str).dropna()
    obs = obs[obs[cat_key] != "nan"]
    counts = pd.crosstab(obs[NICHE_KEY], obs[cat_key])
    counts = counts.loc[counts.sum(axis=1) >= min_niche_cells]
    counts = counts.loc[niche_sorted(counts.index)]
    col_order = ordered_cats(cat_key, counts.columns)
    if col_order:
        counts = counts[col_order]
    props = counts.div(counts.sum(axis=1), axis=0)
    counts.to_csv(fig_dir / f"{prefix}_counts_by_niche.csv")
    props.round(4).to_csv(fig_dir / f"{prefix}_proportions_by_niche.csv")

    # proportions heatmap
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * counts.shape[1]),
                                    0.55 * counts.shape[0] + 2))
    sns.heatmap(props, cmap="mako", annot=True, fmt=".2f",
                cbar_kws={"label": "proportion"}, linewidths=0.3, ax=ax)
    ax.set_title(f"{cat_key} proportions per niche")
    fig.tight_layout()
    out = fig_dir / f"{prefix}_proportions_by_niche_heatmap"
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
    plt.close(fig)

    # binomial + fisher tests per (niche, cat)
    global_p = props.sum(axis=0) / props.shape[0]  # group-neutral global prop (equal niche weight)
    # Actually use total-cell-weighted global proportion:
    global_total = counts.sum(axis=0)
    global_p = global_total / global_total.sum()

    bin_rows, fish_rows = [], []
    total_per_niche = counts.sum(axis=1)
    total_per_cat = counts.sum(axis=0)
    grand = counts.values.sum()

    for niche in counts.index:
        n = int(total_per_niche.loc[niche])
        for cat in counts.columns:
            k = int(counts.loc[niche, cat])
            p0 = float(global_p[cat])
            if n == 0 or p0 == 0:
                continue
            # binomial greater
            try:
                pval_bin = binomtest(k, n, p=p0, alternative="greater").pvalue
            except Exception:
                pval_bin = np.nan
            bin_rows.append({"niche": niche, cat_key: cat, "k": k, "n": n,
                             "p0": p0, "obs_prop": k / n, "pval": pval_bin})
            # fisher 2x2: [[k, n-k], [cat_other, not_cat_other]]
            a = k
            b = int(total_per_cat[cat]) - k
            c_ = int(n) - k
            d = grand - a - b - c_
            try:
                oddsr, pval_f = fisher_exact([[a, b], [c_, d]], alternative="greater")
            except Exception:
                oddsr, pval_f = np.nan, np.nan
            prop_in = k / n if n else np.nan
            prop_out = b / (grand - n) if (grand - n) else np.nan
            fish_rows.append({
                "niche": niche, cat_key: cat, "k": k, "n": n,
                "p0": p0, "obs_prop": prop_in, "other_prop": prop_out,
                "log2_enrichment": np.log2((prop_in + 1e-9) / (prop_out + 1e-9)),
                "oddsr": oddsr, "pval": pval_f,
            })

    bin_df = pd.DataFrame(bin_rows)
    fish_df = pd.DataFrame(fish_rows)
    if len(bin_df):
        bin_df["pval_adj"] = multipletests(bin_df["pval"].fillna(1).values, method="fdr_bh")[1]
    if len(fish_df):
        fish_df["pval_adj"] = multipletests(fish_df["pval"].fillna(1).values, method="fdr_bh")[1]
    bin_df.to_csv(fig_dir / f"{prefix}_enrichment_binomial.csv", index=False)
    fish_df.to_csv(fig_dir / f"{prefix}_enrichment_fisher.csv", index=False)

    # signed-log10(FDR) heatmap
    if len(fish_df):
        fish_df["signed"] = np.sign(fish_df["log2_enrichment"]) * (-np.log10(fish_df["pval_adj"].clip(lower=1e-300)))
        pivot = fish_df.pivot(index="niche", columns=cat_key, values="signed")
        pivot = pivot.loc[counts.index, counts.columns]
        vmax = float(np.nanmax(np.abs(pivot.values))) if np.isfinite(pivot.values).any() else 1.0
        fig, ax = plt.subplots(figsize=(max(6, 0.6 * pivot.shape[1]),
                                        0.55 * pivot.shape[0] + 2))
        sns.heatmap(pivot, cmap="coolwarm", center=0, vmax=vmax, vmin=-vmax,
                    annot=True, fmt=".1f", cbar_kws={"label": "sign(log2FC) * -log10(FDR)"},
                    linewidths=0.3, ax=ax)
        ax.set_title(f"{cat_key} enrichment signed -log10(FDR) per niche")
        fig.tight_layout()
        out = fig_dir / f"{prefix}_enrichment_signedlog10fdr_heatmap"
        fig.savefig(f"{out}.pdf", bbox_inches="tight")
        fig.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
        plt.close(fig)

        # clustermap version
        try:
            cm = sns.clustermap(pivot.fillna(0), cmap="coolwarm", center=0,
                                vmax=vmax, vmin=-vmax, figsize=(max(8, 0.6 * pivot.shape[1] + 4),
                                                                0.55 * pivot.shape[0] + 4),
                                cbar_kws={"label": "sign(log2FC) * -log10(FDR)"})
            out = fig_dir / f"{prefix}_enrichment_signedlog10fdr_clustermap"
            cm.savefig(f"{out}.pdf", bbox_inches="tight")
            cm.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
            plt.close("all")
        except Exception as exc:
            print(f"  [{prefix}] clustermap failed: {exc!r}", flush=True)

        # significant pairs CSV
        sig = fish_df[(fish_df["pval_adj"] < 0.05) & (fish_df["obs_prop"] > fish_df["other_prop"])]
        sig.to_csv(fig_dir / f"{prefix}_enriched_pairs_fdr05.csv", index=False)
    return fish_df


# ======================================================================
# 9. Supplemental merge
# ======================================================================
def parse_supplemental() -> pd.DataFrame:
    if not SUPPL_CSV.exists():
        print(f"[suppl] {SUPPL_CSV} not found -> skipping", flush=True)
        return pd.DataFrame()
    df = pd.read_csv(SUPPL_CSV, sep=";", skiprows=1)
    print(f"[suppl] {df.shape} cols={list(df.columns)[:10]}…", flush=True)
    return df


def attach_supplemental(adata, suppl: pd.DataFrame) -> list[str]:
    """Attach supplemental metadata by TMA_ID keys; returns list of new obs columns."""
    if suppl.empty:
        return []

    def _norm(k: str) -> str:
        if k is None or (isinstance(k, float) and np.isnan(k)):
            return ""
        s = str(k).strip().replace("\xa0", " ").replace("  ", " ")
        return s

    # Build alias index from supplemental sample_ID / ck_number_rna / ck_number_atac / patient_ID
    alias_cols = [c for c in ["sample_ID", "ck_number_rna", "ck_number_atac", "patient_ID"]
                  if c in suppl.columns]
    alias_to_idx = {}
    for i, row in suppl.iterrows():
        for c in alias_cols:
            v = _norm(row[c])
            if v and v != "nan":
                alias_to_idx.setdefault(v.lower(), i)

    def _keys_for(tma: str) -> list[str]:
        """Build candidate keys from a TMA_ID string like 'Combi E2 <pathology no.> IgAN'."""
        if not tma:
            return []
        s = _norm(tma)
        keys = [s]
        # tokens
        toks = s.split()
        # plain tokens
        keys.extend(toks)
        # join of slide + position if 'Utrecht1 A5', 'Apollot0_A1' etc.
        for i in range(len(toks)):
            for j in range(i + 1, min(i + 3, len(toks)) + 1):
                keys.append("".join(toks[i:j]))
                keys.append("_".join(toks[i:j]))
        # K-number normalisation (K094 -> K94)
        for t in list(keys):
            if t.startswith("K"):
                digits = "".join(ch for ch in t if ch.isdigit())
                if digits:
                    keys.append("K" + str(int(digits)))
        # H-number like H<yy>-<n> -> keep
        return list({k.lower() for k in keys if k})

    new_cols = [c for c in suppl.columns if c not in alias_cols]
    obs_new = {f"supp_{c}": [np.nan] * adata.n_obs for c in new_cols}
    matched = 0
    unmatched = []
    tma_series = adata.obs[TMA_KEY].astype(str).values
    for row_idx, tma in enumerate(tma_series):
        cand_keys = _keys_for(tma)
        hit = None
        for k in cand_keys:
            if k in alias_to_idx:
                hit = alias_to_idx[k]
                break
        if hit is None:
            unmatched.append(tma)
            continue
        matched += 1
        for c in new_cols:
            obs_new[f"supp_{c}"][row_idx] = suppl.at[hit, c]

    for col, vals in obs_new.items():
        adata.obs[col] = vals
    n_unique_tma = adata.obs[TMA_KEY].nunique()
    n_unmatched_unique = len(set(unmatched))
    print(f"[suppl] matched {matched}/{len(tma_series)} cells "
          f"({n_unique_tma - n_unmatched_unique}/{n_unique_tma} unique TMAs)", flush=True)
    if unmatched:
        samples_ex = sorted(set(unmatched))[:10]
        print(f"[suppl] example unmatched TMAs: {samples_ex}", flush=True)
    return list(obs_new.keys())


# Canonical display order for ordinal clinical bins. Without this, every plot
# below does `.astype(str)` + crosstab and gets alphabetical order
# ('30-60', '60-90', '<30', '>90'), which is meaningless for an ordered scale.
ORDINAL_CAT_ORDERS: dict[str, list[str]] = {
    "egfr_group": ["<30", "30-60", "60-90", ">90"],
    "ifta_group": ["<5", "5-25", "26-50", ">50"],
}
_NA_TOKENS = {"nan", "NaN", "NA", "None", "<NA>", ""}


def niche_sorted(values) -> list[str]:
    """Niche labels in numeric order ('2' before '10'), not lexical order."""
    return sorted((str(v) for v in values),
                  key=lambda x: (0, int(x)) if x.isdigit() else (1, x))


def ordered_cats(cat_key: str, present) -> list[str] | None:
    """Canonical order for `cat_key` restricted to the labels actually present.

    Returns None for keys with no defined ordinal order, so callers keep their
    previous (frequency / alphabetical) behaviour.
    """
    order = ORDINAL_CAT_ORDERS.get(cat_key)
    if order is None:
        return None
    present = [str(c) for c in present if str(c) not in _NA_TOKENS]
    known = [c for c in order if c in present]
    extra = sorted(c for c in present if c not in order)
    return known + extra


def derive_ifta_egfr_groups(adata) -> None:
    """Bin IFTA and egfr into ordered categorical groups.

    IFTA: <5, 5-25, 26-50, >50.
    eGFR: <30, 30-60, 60-90, >90 (KDIGO-style; half-open bins, so '30-60' is
    30 <= x < 60 and '>90' is x >= 90). Both are stored as *ordered* pandas
    Categoricals so downstream plots keep the clinical order.
    """
    ifta_col = None
    for cand in ("supp_IFTA", "IFTA"):
        if cand in adata.obs.columns:
            ifta_col = cand
            break
    egfr_col = None
    for cand in ("supp_egfr_at_time_of_biopsy", "supp_eGFR_biopsy", "eGFR_biopsy"):
        if cand in adata.obs.columns:
            egfr_col = cand
            break

    if ifta_col is not None:
        def _bin_ifta(v):
            s = str(v).strip()
            if s in ("", "nan", "NaN", "NA", "None"):
                return np.nan
            if s.startswith("<"):
                return "<5"
            if s.startswith(">"):
                return ">50"
            try:
                x = float(s.replace(",", "."))
            except Exception:
                return np.nan
            if x < 5:
                return "<5"
            if x <= 25:
                return "5-25"
            if x <= 50:
                return "26-50"
            return ">50"

        adata.obs["ifta_group"] = pd.Categorical(
            adata.obs[ifta_col].apply(_bin_ifta),
            categories=ORDINAL_CAT_ORDERS["ifta_group"], ordered=True)
        # numeric ifta_pct if we can
        def _ifta_num(v):
            s = str(v).strip().replace(",", ".")
            if s in ("", "nan", "NA"):
                return np.nan
            s = s.replace("<", "").replace(">", "")
            try:
                return float(s)
            except Exception:
                return np.nan
        adata.obs["IFTA_pct"] = adata.obs[ifta_col].apply(_ifta_num)
    if egfr_col is not None:
        def _bin(v):
            # float("nan") SUCCEEDS, and every comparison against NaN is False, so a
            # missing eGFR used to fall through the whole ladder and land in the top
            # bin -- silently counting cores with no eGFR as the healthiest ones
            # (fixed 2026-07-22). Missing must stay missing.
            s = str(v).strip()
            if s.lower() in ("", "nan", "na", "none", "<na>"):
                return np.nan
            try:
                x = float(s.replace(",", "."))
            except Exception:
                return np.nan
            if not np.isfinite(x):
                return np.nan
            if x < 30:
                return "<30"
            if x < 60:
                return "30-60"
            if x < 90:
                return "60-90"
            return ">90"
        adata.obs["egfr_group"] = pd.Categorical(
            adata.obs[egfr_col].apply(_bin),
            categories=ORDINAL_CAT_ORDERS["egfr_group"], ordered=True)
        def _num(v):
            try:
                return float(str(v).replace(",", "."))
            except Exception:
                return np.nan
        adata.obs["EGFR"] = adata.obs[egfr_col].apply(_num)


# ======================================================================
# 10. Cell-88-style derived plots for any fisher_df
# ======================================================================
def plot_bubble_enrichment(fisher_df: pd.DataFrame, cat_key: str, out_path_base: Path,
                           title: str, min_n: int = 10) -> None:
    if fisher_df.empty:
        return
    df = fisher_df.copy()
    df["signed"] = np.sign(df["log2_enrichment"]) * (-np.log10(df["pval_adj"].clip(lower=1e-300)))
    # split alpha on n
    df["alpha"] = np.where(df["n"] < min_n, 0.25, 0.9)
    vmax = float(np.nanmax(np.abs(df["signed"].values))) if np.isfinite(df["signed"].values).any() else 1
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cats = ordered_cats(cat_key, df[cat_key].unique()) or sorted(df[cat_key].unique())
    df = df[df[cat_key].astype(str).isin(cats)]
    niches = sorted(df["niche"].unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    x_idx = {c: i for i, c in enumerate(cats)}
    y_idx = {n: i for i, n in enumerate(niches)}
    fig, ax = plt.subplots(figsize=(0.7 * len(cats) + 3, 0.45 * len(niches) + 3))
    sc = ax.scatter(
        [x_idx[c] for c in df[cat_key]], [y_idx[n] for n in df["niche"]],
        c=df["signed"], cmap="coolwarm", norm=norm,
        s=(df["k"].clip(lower=1) ** 0.5) * 6,
        alpha=df["alpha"],
        edgecolors="black", linewidths=0.3,
    )
    for _, row in df.iterrows():
        if row["pval_adj"] < 0.05 and row["signed"] > 0:
            ax.text(x_idx[row[cat_key]], y_idx[row["niche"]], "*",
                    ha="center", va="center", fontsize=10, color="black")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=45, ha="right")
    ax.set_yticks(range(len(niches)))
    ax.set_yticklabels(niches)
    ax.set_xlabel(cat_key)
    ax.set_ylabel("niche")
    ax.set_title(title)
    cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
    cbar.set_label("sign(log2FC) * -log10(FDR)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}",
                    bbox_inches="tight", dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_stacked_proportions(adata, cat_key: str, out_path_base: Path,
                             min_total: int = 20) -> None:
    counts = pd.crosstab(adata.obs[NICHE_KEY].astype(str),
                         adata.obs[cat_key].astype(str))
    counts = counts.loc[counts.sum(axis=1) >= min_total]
    if counts.empty:
        return
    counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
    col_order = ordered_cats(cat_key, counts.columns)
    if col_order:
        counts = counts[col_order]
    props = counts.div(counts.sum(axis=1), axis=0)
    # reuse main composition plotter style
    plot_niche_composition_from_props(counts, props, cat_key, out_path_base,
                                      cat_order=col_order)


def plot_niche_composition_from_props(counts, props, cat_key: str,
                                      out_path_base: Path,
                                      cat_order=None) -> None:
    if cat_order is None:
        cat_order = counts.sum(axis=0).sort_values(ascending=False).index
    counts = counts[cat_order]
    props = props[cat_order]
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(cat_order), 1)))
    colour_map = {c: palette[i % len(palette)] for i, c in enumerate(cat_order)}
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
    xpos = np.arange(len(counts.index))
    bottom = np.zeros(len(counts.index))
    for cat in cat_order:
        axes[0].bar(xpos, counts[cat].values, bottom=bottom, color=colour_map[cat],
                    width=0.78, label=cat, edgecolor="white", linewidth=0.4)
        bottom += counts[cat].values
    totals = counts.sum(axis=1).values
    tmax = totals.max()
    for x, t in zip(xpos, totals):
        axes[0].text(x, t + tmax * 0.015, f"n={int(t):,}",
                     ha="center", va="bottom", fontsize=8)
    axes[0].set_ylim(0, tmax * 1.12)
    axes[0].set_ylabel("Cell count")
    axes[0].set_title(f"Absolute counts per niche")
    bottom = np.zeros(len(props.index))
    for cat in cat_order:
        axes[1].bar(xpos, props[cat].values, bottom=bottom, color=colour_map[cat],
                    width=0.78, edgecolor="white", linewidth=0.4)
        bottom += props[cat].values
    axes[1].set_ylim(0, 1.0)
    axes[1].set_ylabel("Proportion")
    axes[1].set_title(f"Proportions per niche")
    for ax in axes:
        ax.set_xticks(xpos)
        ax.set_xticklabels([str(n) for n in counts.index])
        ax.set_xlabel(f"Niche ({NICHE_KEY})")
        ax.grid(axis="y", linestyle=":", alpha=0.4)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colour_map[c]) for c in cat_order]
    fig.legend(handles, list(cat_order), title=cat_key,
               loc="center left", bbox_to_anchor=(1.01, 0.5),
               frameon=False, fontsize=8)
    fig.suptitle(f"Niche composition by {cat_key}")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}",
                    bbox_inches="tight", dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_effectsize_heatmap(adata, cat_key: str, out_path_base: Path) -> None:
    counts = pd.crosstab(adata.obs[NICHE_KEY].astype(str),
                         adata.obs[cat_key].astype(str))
    if counts.empty:
        return
    counts = counts.loc[niche_sorted(counts.index)]
    col_order = ordered_cats(cat_key, counts.columns)
    if col_order:
        counts = counts[col_order]
    props = counts.div(counts.sum(axis=1), axis=0)
    global_p = counts.sum(axis=0) / counts.sum().sum()
    effect = props.subtract(global_p, axis=1)
    vmax = float(np.nanmax(np.abs(effect.values))) if np.isfinite(effect.values).any() else 0.5
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * effect.shape[1]),
                                    0.55 * effect.shape[0] + 2))
    sns.heatmap(effect, cmap="coolwarm", center=0, vmax=vmax, vmin=-vmax,
                annot=True, fmt=".2f", cbar_kws={"label": "obs prop - global prop"},
                linewidths=0.3, ax=ax)
    ax.set_title(f"{cat_key} effect size per niche (prop − baseline)")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}", bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_continuous_ecdf(adata, value_key: str, out_path_base: Path,
                         sample_level: bool = True) -> None:
    ser = adata.obs[value_key].astype(float).dropna()
    if ser.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    for niche, df in adata.obs.groupby(NICHE_KEY, observed=True):
        if sample_level:
            vals = df.groupby(TMA_KEY, observed=True)[value_key].median().dropna().values
        else:
            vals = df[value_key].dropna().values
        if len(vals) == 0:
            continue
        xs = np.sort(vals.astype(float))
        ys = np.arange(1, len(xs) + 1) / len(xs)
        ax.plot(xs, ys, label=f"{niche} (n={len(xs)})", lw=1.4)
    ax.set_xlabel(value_key + (" (per-TMA median)" if sample_level else ""))
    ax.set_ylabel("ECDF")
    ax.set_title(f"{value_key} ECDF by niche")
    ax.grid(ls=":", alpha=0.4)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}", bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_violin_by_niche(adata, value_key: str, out_path_base: Path) -> None:
    df = adata.obs[[NICHE_KEY, TMA_KEY, value_key]].copy()
    df[value_key] = pd.to_numeric(df[value_key], errors="coerce")
    agg = df.groupby([NICHE_KEY, TMA_KEY], observed=True)[value_key].median().reset_index().dropna()
    if agg.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 0.4 * agg[NICHE_KEY].nunique() + 3))
    # order by median desc
    order = agg.groupby(NICHE_KEY, observed=True)[value_key].median().sort_values(ascending=False).index
    sns.violinplot(data=agg, x=value_key, y=NICHE_KEY, order=order, color="lightgrey", inner=None, ax=ax)
    sns.stripplot(data=agg, x=value_key, y=NICHE_KEY, order=order, color="tab:red",
                  size=3, ax=ax)
    ax.set_xlabel(value_key)
    ax.set_ylabel("niche")
    ax.set_title(f"{value_key} per TMA, grouped by niche")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}", bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_lollipop(adata, cat_key: str, niche_of_interest: str,
                  out_path_base: Path, top_n: int = 20) -> None:
    counts = pd.crosstab(adata.obs[NICHE_KEY].astype(str),
                         adata.obs[cat_key].astype(str))
    if niche_of_interest not in counts.index:
        return
    col_order = ordered_cats(cat_key, counts.columns)
    if col_order:
        counts = counts[col_order]
    props = counts.div(counts.sum(axis=1), axis=0)
    global_p = counts.sum(axis=0) / counts.sum().sum()
    diff = props.loc[niche_of_interest] - global_p
    # ordinal categories keep their clinical order (worst at the bottom of the
    # y-axis); everything else stays ranked by effect size as before
    diff = diff.reindex(col_order[::-1]) if col_order else diff.sort_values()
    diff = diff.iloc[:top_n] if len(diff) > top_n else diff
    fig, ax = plt.subplots(figsize=(9, 0.35 * len(diff) + 2))
    ax.hlines(y=diff.index, xmin=0, xmax=diff.values, color="steelblue")
    ax.scatter(diff.values, diff.index, color="steelblue")
    ax.axvline(0, ls=":", color="grey", lw=1)
    ax.set_xlabel("prop_in − global")
    ax.set_title(f"Niche {niche_of_interest}: {cat_key} enrichment vs global")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_path_base}.{ext}", bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


# ======================================================================
# 11. ROI flags
# ======================================================================
def ingest_roi_flags(adata, fig_dir: Path) -> list[str]:
    print(f"\n[roi] reading {len(ROI_FILES)} ROI csvs", flush=True)
    combined = {}
    all_flag_cols: set[str] = set()
    for p in ROI_FILES:
        if not p.exists():
            print(f"[roi] missing {p}", flush=True)
            continue
        df = pd.read_csv(p)
        flag_cols = [c for c in df.columns if c.startswith("in_")]
        all_flag_cols.update(flag_cols)
        for _, row in df.iterrows():
            cid = row.get("cell_id")
            smp = row.get("sample")
            if pd.isna(cid):
                continue
            key = (str(smp), str(cid))
            if key not in combined:
                combined[key] = {c: False for c in flag_cols}
            for c in flag_cols:
                val = row.get(c, False)
                if isinstance(val, str):
                    val = val.strip().lower() in ("true", "1", "yes", "t")
                combined[key][c] = combined[key].get(c, False) or bool(val)
        print(f"  {p.name} -> {len(df)} rows, flags={flag_cols}", flush=True)

    print(f"[roi] unique cell keys collected: {len(combined)}", flush=True)
    all_flag_cols = sorted(all_flag_cols)
    # attach to adata.obs
    keys = list(zip(adata.obs[SAMPLE_KEY].astype(str), adata.obs["cell_id"].astype(str)))
    for c in all_flag_cols:
        vals = [combined.get(k, {}).get(c, False) for k in keys]
        adata.obs[c] = vals
    matched = sum(1 for k in keys if k in combined)
    print(f"[roi] attached flags: {all_flag_cols}", flush=True)
    print(f"[roi] matched {matched}/{adata.n_obs} cells", flush=True)

    # per-niche %True table
    roi_dir = fig_dir / "roi_flags"
    roi_dir.mkdir(exist_ok=True)
    summary = {}
    for c in all_flag_cols:
        s = adata.obs[[NICHE_KEY, c]].copy()
        s[c] = s[c].astype(bool)
        summary[c] = s.groupby(NICHE_KEY, observed=True)[c].mean() * 100.0
    summary_df = pd.DataFrame(summary)
    summary_df.round(2).to_csv(roi_dir / "percent_true_by_niche_all_flags.csv")
    # heatmap
    if not summary_df.empty:
        fig, ax = plt.subplots(figsize=(max(6, 0.6 * summary_df.shape[1]),
                                        0.55 * summary_df.shape[0] + 2))
        sns.heatmap(summary_df, cmap="viridis", annot=True, fmt=".1f",
                    cbar_kws={"label": "% True"}, linewidths=0.3, ax=ax)
        ax.set_title("ROI flag % True per niche")
        fig.tight_layout()
        out = roi_dir / "percent_true_heatmap"
        fig.savefig(f"{out}.pdf", bbox_inches="tight")
        fig.savefig(f"{out}.png", bbox_inches="tight", dpi=160)
        plt.close(fig)
        # per-flag barplots
        for c in all_flag_cols:
            fig, ax = plt.subplots(figsize=(8, 5))
            s = summary_df[c].sort_values(ascending=False)
            ax.bar(s.index.astype(str), s.values, color="steelblue")
            ax.set_xlabel("niche")
            ax.set_ylabel("% True")
            ax.set_title(f"ROI flag {c} per niche")
            for i, v in enumerate(s.values):
                ax.text(i, v + 0.5, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
            fig.tight_layout()
            for ext in ("pdf", "png"):
                fig.savefig(roi_dir / f"pct_by_niche__{c}.{ext}",
                            bbox_inches="tight", dpi=160 if ext == "png" else None)
            plt.close(fig)

    # cell_id x flag long table (for downstream reuse)
    long_rows = []
    for (smp, cid), flags in combined.items():
        r = {"sample": smp, "cell_id": cid}
        r.update(flags)
        long_rows.append(r)
    pd.DataFrame(long_rows).to_csv(roi_dir / "combined_in_flags_by_cell_id.csv", index=False)
    print(f"[roi] wrote summary + per-flag plots to {roi_dir}", flush=True)
    return all_flag_cols


# ======================================================================
# 12. per-niche within-niche DE
# ======================================================================
def per_niche_celltype_markers(adata, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n[de] within-niche DE on {CELL_TYPE_KEY}", flush=True)
    # normalise in a copy so we don't mutate the main adata
    ad2 = adata.copy()
    if "counts" not in ad2.layers:
        ad2.layers["counts"] = ad2.X.copy()
    sc.pp.normalize_total(ad2, target_sum=1e4)
    sc.pp.log1p(ad2)
    ad2.layers["lognorm"] = ad2.X.copy()
    ad2.X = ad2.layers["counts"]

    combined_tidy = []
    sc.settings.figdir = out_dir
    for niche in sorted(ad2.obs[NICHE_KEY].astype(str).unique(),
                        key=lambda x: int(x) if x.isdigit() else x):
        sub = ad2[ad2.obs[NICHE_KEY].astype(str) == niche].copy()
        # drop tiny cell-type groups
        ct_counts = sub.obs[CELL_TYPE_KEY].astype(str).value_counts()
        keep = ct_counts[ct_counts >= 5].index.tolist()
        sub = sub[sub.obs[CELL_TYPE_KEY].astype(str).isin(keep)].copy()
        if sub.n_obs < 20 or sub.obs[CELL_TYPE_KEY].nunique() < 2:
            print(f"  niche {niche}: skipping ({sub.n_obs} cells, {sub.obs[CELL_TYPE_KEY].nunique()} cts)", flush=True)
            continue
        try:
            sc.tl.rank_genes_groups(
                sub, groupby=CELL_TYPE_KEY, method="wilcoxon", n_genes=200,
                use_raw=False, layer="lognorm", corr_method="benjamini-hochberg",
                tie_correct=True, pts=True,
            )
        except Exception as exc:
            print(f"  niche {niche}: rank_genes failed: {exc!r}", flush=True)
            continue

        # tidy
        rgg = sub.uns["rank_genes_groups"]
        groups = rgg["names"].dtype.names
        rows = []
        for g in groups:
            for i in range(len(rgg["names"][g])):
                rows.append({
                    "niche": niche, "group": g,
                    "gene": rgg["names"][g][i],
                    "score": rgg["scores"][g][i],
                    "logfoldchange": rgg["logfoldchanges"][g][i],
                    "pval": rgg["pvals"][g][i],
                    "pval_adj": rgg["pvals_adj"][g][i],
                    "pct_in_group": rgg["pts"].loc[rgg["names"][g][i], g] if "pts" in rgg else np.nan,
                    "pct_out_group": (
                        rgg["pts_rest"].loc[rgg["names"][g][i], g]
                        if "pts_rest" in rgg else np.nan
                    ),
                })
        tidy = pd.DataFrame(rows)
        tidy.to_csv(out_dir / f"markers_celltypes_within_cluster_{niche}.csv", index=False)
        combined_tidy.append(tidy)

        try:
            sc.pl.rank_genes_groups(sub, n_genes=min(20, 200), sharey=False,
                                    save=f"_cluster_{niche}_ranked.png", show=False)
            sc.pl.rank_genes_groups_heatmap(sub, n_genes=10, swap_axes=True,
                                            dendrogram=False, use_raw=False,
                                            save=f"_cluster_{niche}_heatmap.png", show=False)
            # top 5 per group, dedup
            top5 = []
            for g in groups:
                top5.extend([rgg["names"][g][i] for i in range(min(5, len(rgg["names"][g])))])
            top5 = list(dict.fromkeys(top5))
            sc.pl.dotplot(sub, var_names=top5, groupby=CELL_TYPE_KEY,
                          standard_scale="var", use_raw=False,
                          save=f"_cluster_{niche}_dotplot.png", show=False)
        except Exception as exc:
            print(f"  niche {niche}: plot failed: {exc!r}", flush=True)
        plt.close("all")

    if combined_tidy:
        all_df = pd.concat(combined_tidy, ignore_index=True)
        all_df.to_csv(out_dir / "markers_celltypes_within_each_cluster_combined.csv", index=False)
        filt = all_df[(all_df["pval_adj"] < 0.05) & (all_df["logfoldchange"] > 0.25)]
        filt.to_csv(out_dir / "markers_celltypes_within_each_cluster_filtered.csv", index=False)
        print(f"  combined DE: {all_df.shape}; filtered: {filt.shape}", flush=True)


# ======================================================================
# 13. niche annotation summary table
# ======================================================================
def write_niche_annotation(adata, fig_dir: Path) -> None:
    niches = sorted(adata.obs[NICHE_KEY].astype(str).unique(),
                    key=lambda x: int(x) if x.isdigit() else x)
    rows = []
    for n in niches:
        mask = adata.obs[NICHE_KEY].astype(str) == n
        sub = adata.obs[mask]
        ctop = sub[CELL_TYPE_KEY].astype(str).value_counts(normalize=True).head(3)
        dtop = sub[CONDITION_KEY].astype(str).value_counts(normalize=True).head(3)
        rows.append({
            "niche": n,
            "n_cells": int(mask.sum()),
            "top_cell_types": "; ".join(f"{k} ({v*100:.1f}%)" for k, v in ctop.items()),
            "top_conditions": "; ".join(f"{k} ({v*100:.1f}%)" for k, v in dtop.items()),
        })
    pd.DataFrame(rows).to_csv(fig_dir / "niche_annotation_summary.csv", index=False)
    print(f"  wrote niche_annotation_summary.csv", flush=True)


# ======================================================================
# main
# ======================================================================
def main() -> int:
    global LEIDEN_RES, NICHE_KEY, CELL_TYPE_KEY
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=str, default=str(ARTIFACTS_DIR / "22042026_221229"))
    ap.add_argument("--leiden-res", type=float, default=LEIDEN_RES)
    ap.add_argument("--de-dir-name", type=str, default="niche_marker_results")
    ap.add_argument("--out-h5ad-name", type=str, default="xenium_with_niches.h5ad")
    ap.add_argument("--skip-de", action="store_true", help="Skip within-niche DE (slow)")
    ap.add_argument("--skip-pathways", action="store_true", help="Skip hallmark/reactome ULM")
    ap.add_argument("--cell-type-key", type=str, default=CELL_TYPE_KEY,
                    help="obs column for cell type (default: predic). If absent it is "
                         "merged from --labeled-h5ad on (sample, cell_id).")
    ap.add_argument("--out-dir", type=str, default="",
                    help="output base for figures/DE/h5ad (default: run_dir)")
    ap.add_argument("--labeled-h5ad", type=str,
                    default=f"{KB_ANALYSIS}/"
                            "celltype_concordance/outputs/xenium_labeled.h5ad")
    args = ap.parse_args()

    LEIDEN_RES = args.leiden_res
    NICHE_KEY = f"latent_leiden_{LEIDEN_RES}"
    CELL_TYPE_KEY = args.cell_type_key
    print(f"[main] leiden res={LEIDEN_RES}, niche_key={NICHE_KEY}, "
          f"cell_type_key={CELL_TYPE_KEY}", flush=True)

    run_dir = Path(args.run_dir)
    out_base = Path(args.out_dir) if args.out_dir else run_dir
    fig_dir = out_base / "figures"
    de_dir = out_base / args.de_dir_name
    fig_dir.mkdir(parents=True, exist_ok=True)

    model, adata = load_trained(run_dir)

    # attach refined cell typing if the requested key is absent (e.g. level2_refined)
    if CELL_TYPE_KEY not in adata.obs.columns:
        sys.path.insert(0, str(Path(__file__).resolve().parent / "nichecompass_replot_level2"))
        from level2_common import merge_level2
        adata, _ = merge_level2(adata, Path(args.labeled_h5ad), key=CELL_TYPE_KEY)

    def _safe(tag, fn, *a, **kw):
        print(f"\n--- {tag} ---", flush=True)
        try:
            return fn(*a, **kw)
        except Exception as exc:
            print(f"[{tag}] FAILED: {exc!r}", flush=True)
            traceback.print_exc()
            return None

    # 1. composition plots
    def _override(cat_key):
        ck = f"{cat_key}_colors"
        if ck in adata.uns and pd.api.types.is_categorical_dtype(adata.obs[cat_key]):
            return dict(zip(adata.obs[cat_key].cat.categories, adata.uns[ck]))
        return None
    for cat_key in (CELL_TYPE_KEY, CONDITION_KEY, SAMPLE_KEY):
        _safe(f"composition: {cat_key}", plot_niche_composition,
              adata, NICHE_KEY, cat_key,
              fig_dir / f"res_{LEIDEN_RES}_niche_composition_{cat_key}",
              title_suffix=f"  (n={adata.n_obs:,} cells)",
              colour_override=_override(cat_key))

    # 3. niche x cell-type heatmap
    _safe("niche x cell-type heatmap", heatmap_niche_vs_celltype, adata, fig_dir)

    # 4. GP activity heatmap
    _safe("GP activity heatmap", gp_activity_heatmap, model, adata, fig_dir)

    # 5+6. pathway activity
    if not args.skip_pathways:
        _safe("Hallmark pathway activity", pathway_activity, adata, fig_dir, hallmark_fn, "hallmark")
        _safe("Reactome pathway activity", pathway_activity, adata, fig_dir, reactome_fn, "reactome")

    # 7. per-TMA spatial grid
    _safe("per-TMA spatial grid", spatial_niche_grid, adata, fig_dir)

    # 8. condition enrichment
    fisher_cond = _safe("condition enrichment", run_categorical_enrichment,
                        adata, CONDITION_KEY, fig_dir, "condition")

    # 9. clinical metadata merge.
    # Use the identifier-resolved per-core table (celltype_concordance/
    # 15_update_metadata_slope.py). The legacy `attach_supplemental()` below matches
    # on plain name tokens -- TMA well 'b4', group words 'apollo'/'fsgs'/'dn' -- and
    # takes the first key that hits, which mis-assigned clinical records to 44% of
    # cores (see PLAN 2026-07-22). It is kept only as a fallback and warns loudly.
    def _merge_clinical(adata):
        sys.path.insert(0, str(Path(__file__).resolve().parent / "nichecompass_replot_level2"))
        from metadata_common import apply_core_metadata, load_percore
        apply_core_metadata(adata, load_percore())

    if _safe("clinical metadata (identifier join)", _merge_clinical, adata) is None:
        print("[suppl] WARNING falling back to the LEGACY token matcher -- its "
              "clinical assignments are known to be unreliable", flush=True)
        suppl = _safe("supplemental parse", parse_supplemental)
        if suppl is not None and not suppl.empty:
            _safe("supplemental attach", attach_supplemental, adata, suppl)
    _safe("derive IFTA/EGFR groups", derive_ifta_egfr_groups, adata)

    # 10a. IFTA / EGFR enrichment (cell-88 framework)
    for gkey in ("ifta_group", "egfr_group"):
        if gkey not in adata.obs.columns:
            continue
        fish = _safe(f"{gkey} enrichment", run_categorical_enrichment,
                     adata, gkey, fig_dir, gkey, min_niche_cells=20)
        if fish is None or fish.empty:
            continue
        _safe(f"{gkey} bubble", plot_bubble_enrichment,
              fish, gkey, fig_dir / f"{gkey}_bubble",
              f"{gkey} enrichment per niche (bubble)", min_n=20)
        _safe(f"{gkey} stacked", plot_stacked_proportions, adata, gkey,
              fig_dir / f"{gkey}_stacked_bars", min_total=50)
        _safe(f"{gkey} effect heatmap", plot_effectsize_heatmap,
              adata, gkey, fig_dir / f"{gkey}_effectsize_heatmap")
        pick = fish.sort_values("pval_adj").iloc[0]
        _safe(f"{gkey} lollipop", plot_lollipop, adata, gkey, pick["niche"],
              fig_dir / f"{gkey}_lollipop_{pick['niche']}", 20)
    # 10b. Continuous EGFR / IFTA_pct ECDF + violin
    for value_key, tag in (("EGFR", "EGFR"), ("IFTA_pct", "IFTA")):
        if value_key not in adata.obs.columns:
            continue
        _safe(f"{tag} ECDF", plot_continuous_ecdf, adata, value_key,
              fig_dir / f"{tag}_ecdf_samplelevel", True)
        _safe(f"{tag} violin", plot_violin_by_niche, adata, value_key,
              fig_dir / f"{tag}_ridgeline_violin")

    # 11. ROI flags
    _safe("ROI flag ingestion", ingest_roi_flags, adata, fig_dir)

    # 13. Niche annotation summary
    _safe("niche annotation summary", write_niche_annotation, adata, fig_dir)

    # 12. within-niche DE (slowest; last so a failure doesn't trash prior outputs)
    if not args.skip_de:
        _safe("within-niche DE", per_niche_celltype_markers, adata, de_dir)

    # persist enriched adata
    out_h5ad = out_base / args.out_h5ad_name
    print(f"\n--- writing enriched adata -> {out_h5ad} ---", flush=True)
    # clean up bool columns to avoid anndata write issues
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str)
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"[done] run_dir={run_dir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
