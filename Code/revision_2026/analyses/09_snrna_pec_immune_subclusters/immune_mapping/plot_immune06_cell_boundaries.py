#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Xenium cell-boundary plots coloured by the transferred immune leiden_0.6 subclusters.

Faithful adaptation of cell_boundary_plots.ipynb (PEC section): geopandas polygon plotting, obs joined
onto the cell-boundary GeoDataFrame by `cell_id`, one figure per TMA core, unmatched cells light grey.
Changes vs PEC: input AnnData = the immune mapping output; colour by `predictions` (18 subclusters);
context panel = Xenium `cell_type_manual` (Immune/PL/B/B). Env: nichecompass (geopandas)."""
from kbpaths import KB_ANALYSIS, KB_DATA  # repository path config (config/kbpaths.py)
import os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib.patches import Patch
import anndata as ad

# ---------------------- Config ----------------------
RES_DIR = f"{KB_DATA}/processed/xenium/"
CELL_BOUNDARIES_PATH = os.path.join(RES_DIR, "geodataframes/cell_boundaries.parquet")
NUCLEUS_BOUNDARIES_PATH = os.path.join(RES_DIR, "geodataframes/nucleus_boundaries.parquet")
ANNDATA_PATH = f"{KB_ANALYSIS}/immune_mapping_leiden06/prediction/xenium_immune_with_prediction.h5ad"
OUT_BASE = f"{KB_ANALYSIS}/immune_mapping_leiden06/cell_boundary_plots_immune_0_60"

COL_CTX = "cell_type_manual"   # context panel (coarse Xenium immune annotation)
COL_PRED = "predictions"       # transferred immune leiden_0.6 subcluster (main)
GDF_JOIN_COL = "cell_id"
FORCE_OBS_CELL_ID_COL = "cell_id"

NA_COLOR = "#D3D3D3"; EDGE_COLOR = "black"; EDGE_WIDTH = 0.15; FILL_ALPHA = 0.85
NUC_COLOR = "grey"; NUC_WIDTH = 0.15; DPI = 200
SUB_CTX = "cell_type_manual"; SUB_PRED = "predictions"; SUB_SIDE = "side_by_side__ctx__vs__predictions"


# ---------------------- Utils ----------------------
def ensure_dirs():
    for d in [OUT_BASE, os.path.join(OUT_BASE, SUB_CTX), os.path.join(OUT_BASE, SUB_PRED),
              os.path.join(OUT_BASE, SUB_SIDE), os.path.join(OUT_BASE, "_logs")]:
        os.makedirs(d, exist_ok=True)

def norm_id(series):
    s = pd.Series(series, copy=True).astype(str).str.strip()
    s = s.str.replace(r"\.0$", "", regex=True).str.replace(r"\s+", "", regex=True)
    return s

def safe_filename(s):
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in str(s))

def _catsort(cats):
    """numeric-aware sort so predictions 0..17 order correctly (not lexicographic)."""
    return sorted(cats, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))

def make_palette(categories):
    categories = list(categories)
    base = list(plt.get_cmap("tab20").colors)
    colors = base[:len(categories)] if len(categories) <= len(base) \
        else base + list(sns.color_palette("husl", len(categories) - len(base)))
    return {cat: colors[i] for i, cat in enumerate(categories)}

def unify_to_str(cat_series):
    s = pd.Series(cat_series, copy=True)
    if pd.api.types.is_categorical_dtype(s):
        s = s.astype("object")
    return s.where(s.isna(), s.astype(str))

def build_legend_handles(categories, palette_dict, present_set):
    handles, labels = [], []
    for cat in categories:
        if cat in present_set:
            handles.append(Patch(facecolor=palette_dict[cat], edgecolor="none")); labels.append(str(cat))
    return handles, labels

def plot_core_categorical(gdf_core, gdf_core_nuc, column, palette_dict, out_path, title=None, legend_title=None):
    vals = unify_to_str(gdf_core[column]); colors = vals.map(palette_dict).fillna(NA_COLOR)
    fig, ax = plt.subplots(figsize=(14, 12))
    gdf_core.plot(color=colors, edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH, alpha=FILL_ALPHA, ax=ax)
    if gdf_core_nuc is not None and len(gdf_core_nuc) > 0:
        gdf_core_nuc.boundary.plot(color=NUC_COLOR, linewidth=NUC_WIDTH, ax=ax, zorder=10)
    present = set(vals.dropna().unique().tolist())
    handles, labels = build_legend_handles(list(palette_dict.keys()), palette_dict, present)
    if handles:
        leg = ax.legend(handles=handles, labels=labels, title=legend_title or column,
                        bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, fontsize=10)
        if leg and leg.get_title(): leg.get_title().set_fontsize(11)
    ax.set_ylim(ax.get_ylim()[::-1]); ax.set_title(title or column, fontsize=16)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in ax.spines.values(): spine.set_visible(False)
    plt.tight_layout(); fig.savefig(out_path, bbox_inches="tight", dpi=DPI); plt.close(fig)

def plot_core_side_by_side(gdf_core, gdf_core_nuc, col_left, pal_left, col_right, pal_right, out_path, suptitle=None):
    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    for ax, col, pal in [(axes[0], col_left, pal_left), (axes[1], col_right, pal_right)]:
        vals = unify_to_str(gdf_core[col]); colors = vals.map(pal).fillna(NA_COLOR)
        gdf_core.plot(color=colors, edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH, alpha=FILL_ALPHA, ax=ax)
        if gdf_core_nuc is not None and len(gdf_core_nuc) > 0:
            gdf_core_nuc.boundary.plot(color=NUC_COLOR, linewidth=NUC_WIDTH, ax=ax, zorder=10)
        handles, labels = build_legend_handles(list(pal.keys()), pal, set(vals.dropna().unique().tolist()))
        if handles:
            leg = ax.legend(handles=handles, labels=labels, title=col, bbox_to_anchor=(1.0, 0.5),
                            loc="center left", frameon=False, fontsize=10)
            if leg and leg.get_title(): leg.get_title().set_fontsize(11)
        ax.set_ylim(ax.get_ylim()[::-1]); ax.set_title(col, fontsize=16)
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        for spine in ax.spines.values(): spine.set_visible(False)
    if suptitle: fig.suptitle(suptitle, fontsize=18)
    plt.tight_layout(); fig.savefig(out_path, bbox_inches="tight", dpi=DPI); plt.close(fig)


# ---------------------- Main ----------------------
def main():
    ensure_dirs()
    print("[Info] Loading GeoDataFrames ...", flush=True)
    gdf_cells = gpd.read_parquet(CELL_BOUNDARIES_PATH)
    gdf_nuc = gpd.read_parquet(NUCLEUS_BOUNDARIES_PATH)
    assert GDF_JOIN_COL in gdf_cells.columns, f"{GDF_JOIN_COL} not in cell boundaries: {list(gdf_cells.columns)}"

    print("[Info] Loading AnnData ...", flush=True)
    adata = ad.read_h5ad(ANNDATA_PATH)
    for col in [COL_PRED, FORCE_OBS_CELL_ID_COL]:
        assert col in adata.obs.columns, f"'{col}' missing in obs: {list(adata.obs.columns)}"
    # COL_CTX (cell_type_manual) comes from the boundary parquet itself (covers all 198k cells)
    assert COL_CTX in gdf_cells.columns, f"'{COL_CTX}' missing in cell boundaries: {list(gdf_cells.columns)}"

    # join ONLY the new prediction label onto ALL cell boundaries by cell_id (immune cells coloured,
    # everything else -> NA/grey). cell_type_manual already exists on gdf_cells -> don't re-merge it
    # (would collide into cell_type_manual_x/_y).
    if COL_PRED in gdf_cells.columns:
        gdf_cells = gdf_cells.drop(columns=[COL_PRED])
    obs_df = adata.obs[[COL_PRED]].copy()
    obs_df["join_id"] = norm_id(adata.obs[FORCE_OBS_CELL_ID_COL])
    obs_df = obs_df.dropna(subset=["join_id"])
    obs_df = obs_df[~obs_df["join_id"].duplicated(keep="first")]
    gdf_cells["join_id"] = norm_id(gdf_cells[GDF_JOIN_COL])

    merged = gdf_cells.merge(obs_df[["join_id", COL_PRED]], on="join_id", how="left")
    merged = gpd.GeoDataFrame(merged, geometry="geometry", crs=gdf_cells.crs).drop(columns=["join_id"])
    n_pred = int(merged[COL_PRED].notna().sum())
    print(f"[Info] matched {n_pred}/{len(merged)} boundary cells to a prediction "
          f"(should ~= {adata.n_obs} immune cells)", flush=True)
    pd.DataFrame({"n_boundary_cells": [len(merged)], "n_matched_predictions": [n_pred],
                  "n_anndata_cells": [adata.n_obs]}).to_csv(
        os.path.join(OUT_BASE, "_logs", "match_report.csv"), index=False)
    if n_pred == 0:
        sys.exit("[Error] no cell_id matches — check FORCE_OBS_CELL_ID_COL / boundary cell_id")

    # palettes (numeric-aware for the 18 predictions)
    ctx_cats = _catsort(unify_to_str(merged[COL_CTX]).dropna().unique().tolist())
    pred_cats = _catsort(unify_to_str(merged[COL_PRED]).dropna().unique().tolist())
    pal_ctx = make_palette(ctx_cats); pal_pred = make_palette(pred_cats)
    print(f"[Info] predictions categories ({len(pred_cats)}): {pred_cats}", flush=True)

    assert "TMA_ID" in merged.columns, "TMA_ID not in cell boundaries"
    cores = merged["TMA_ID"].dropna().unique().tolist()
    print(f"[Info] {len(cores)} TMA cores to plot.", flush=True)
    for i, core in enumerate(cores, 1):
        gdf_core = merged.loc[merged["TMA_ID"] == core].copy()
        gdf_core_nuc = gdf_nuc.loc[gdf_nuc["TMA_ID"] == core].copy() if "TMA_ID" in gdf_nuc.columns else None
        cs = safe_filename(core)
        plot_core_categorical(gdf_core, gdf_core_nuc, COL_PRED, pal_pred,
                              os.path.join(OUT_BASE, SUB_PRED, f"predictions_{cs}.png"),
                              title=f"{core} — immune leiden_0.6 predictions", legend_title="leiden_0.6")
        plot_core_categorical(gdf_core, gdf_core_nuc, COL_CTX, pal_ctx,
                              os.path.join(OUT_BASE, SUB_CTX, f"celltypemanual_{cs}.png"),
                              title=f"{core} — {COL_CTX}", legend_title=COL_CTX)
        plot_core_side_by_side(gdf_core, gdf_core_nuc, COL_CTX, pal_ctx, COL_PRED, pal_pred,
                               os.path.join(OUT_BASE, SUB_SIDE, f"ctx_vs_predictions_{cs}.png"), suptitle=str(core))
        print(f"[{i}/{len(cores)}] {core}", flush=True)
    print("[Done] plots ->", OUT_BASE, flush=True)


if __name__ == "__main__":
    main()
