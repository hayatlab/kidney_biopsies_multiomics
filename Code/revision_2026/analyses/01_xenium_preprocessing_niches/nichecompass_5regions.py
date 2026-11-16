"""NicheCompass on the 5-region merged Xenium aggregate.

Adapted from `$KB_PROJECT/niche_slurm.py` (the submitted
TMA-filtered run). Changes vs the original:

* Input: 5-region merged h5ad produced by Phase B.
* Cell-type key: `predic` (the label-transferred class; `cell_type_manual`
  is NaN on the new Region 1 cells).
* Adds "Combi D2 <pathology no.> IgAN" to the TMA exclusion list — that core
  failed QC (median 44 counts/cell vs 197 across the other 41734 R1 cores,
  EC-marker expression depleted 4-5x, 70% of cells mis-labelled as EC).
* Ensures `layers["counts"]` exists (merged h5ad only has counts on `.X`).
* Passes `lambda_l1_addon` into `.train()` (was defined but dropped in the
  original). Sets `load_timestamp` after save. Guards the CCK-network
  step so a missing GP no longer crashes the whole downstream block.
* Composition plots: side-by-side absolute counts + proportions with
  per-niche totals annotated above each bar, consistent colour palette
  across panels, sorted by niche size; saved as both SVG and PNG.
"""
from __future__ import annotations
from kbpaths import KB_PROJECT, KB_XENIUM_RAW  # repository path config (config/kbpaths.py)
from cores import excluded_core_labels  # private core lists (config/cores.py)

import argparse
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData

import nichecompass as nc
from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    compute_communication_gp_network,
    create_new_color_dict,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
    filter_and_combine_gp_dict_gps_v2,
    generate_enriched_gp_info_plots,
    visualize_communication_gp_network,
)


# ----- defaults (match niche_slurm.py / notebook) ------------------------
INPUT_H5AD = Path(
    f"{KB_XENIUM_RAW}/h5ad_files_only_transcript/"
    "marker_genes/clustered_annotated_final_file_5regions.h5ad"
)
PROJECT_DIR = Path(f"{KB_PROJECT}")
GP_DATA_DIR = PROJECT_DIR / "niche_compass_results/data/gene_programs"
ARTIFACTS_DIR = PROJECT_DIR / "niche_compass_results/artifacts/single_sample"

# obs / obsm / obsp / var / layer keys
CELL_TYPE_KEY = "predic"
SAMPLE_KEY = "sample"
CONDITION_KEY = "condition"
TMA_KEY = "TMA_ID"
SPATIAL_KEY = "spatial"
ADJ_KEY = "spatial_connectivities"
COUNTS_KEY = "counts"

# NicheCompass mask keys (match notebook / niche_slurm.py defaults)
GP_TARGETS_MASK_KEY = "nichecompass_gp_targets"
GP_TARGETS_CATEGORIES_MASK_KEY = "nichecompass_gp_targets_categories"
GP_SOURCES_MASK_KEY = "nichecompass_gp_sources"
GP_SOURCES_CATEGORIES_MASK_KEY = "nichecompass_gp_sources_categories"
GP_NAMES_KEY = "nichecompass_gp_names"
ACTIVE_GP_NAMES_KEY = "nichecompass_active_gp_names"
LATENT_KEY = "nichecompass_latent"

# graph + model + training hyperparameters
N_NEIGHBORS = 4
N_EPOCHS = 400
N_EPOCHS_ALL_GPS = 25
LR = 1e-3
LAMBDA_EDGE_RECON = 5e5
LAMBDA_GENE_EXPR_RECON = 300.0
LAMBDA_L1_MASKED = 0.0
LAMBDA_L1_ADDON = 30.0
EDGE_BATCH_SIZE = 1024
N_SAMPLED_NEIGHBORS = 4
LEIDEN_RES = 0.4

# existing TMA cores to drop (same list as niche_slurm.py) + the new QC-failed core
# Cores excluded before NicheCompass: unassigned cells, the Apollo d0 A2 control core and ten disease cores that
# failed QC (Utrecht cores of K4 DN, K24 TIN, K29 MN, K23 TIN, K25 TIN, K71 FSGS, K46 MN, K78 FSGS and K60 IgAN,
# and the Combi D2 IgAN core). Their internal labels contain pathology numbers and are read from a private list
# (config/cores.py); none of these cores is part of the de-identified Zenodo object.
TO_EXCLUDE = set(excluded_core_labels("nichecompass_blacklist"))


# ----- helpers -----------------------------------------------------------
def load_adata() -> AnnData:
    print(f"[load] reading {INPUT_H5AD}", flush=True)
    a = sc.read_h5ad(INPUT_H5AD)
    print(f"[load] full: {a.shape}", flush=True)

    # counts layer: Phase B left raw counts on .X only
    if COUNTS_KEY not in a.layers:
        print(f"[load] adding layers['{COUNTS_KEY}'] from .X", flush=True)
        a.layers[COUNTS_KEY] = a.X.copy()

    # filter TMA_ID blocklist + NaN
    before = a.n_obs
    mask = (~a.obs[TMA_KEY].astype(str).isin(TO_EXCLUDE)) & a.obs[TMA_KEY].notna() \
        & (a.obs[TMA_KEY].astype(str) != "nan")
    a = a[mask].copy()
    if pd.api.types.is_categorical_dtype(a.obs[TMA_KEY]):
        a.obs[TMA_KEY] = a.obs[TMA_KEY].cat.remove_unused_categories()
    print(f"[load] after TMA filter: {a.shape} (dropped {before - a.n_obs})", flush=True)

    # sanity: predic must be non-null
    nnull = a.obs[CELL_TYPE_KEY].isna().sum()
    if nnull:
        print(f"[load] WARNING: {nnull} cells with null {CELL_TYPE_KEY}, dropping", flush=True)
        a = a[a.obs[CELL_TYPE_KEY].notna()].copy()

    print(f"[load] samples: {sorted(a.obs[SAMPLE_KEY].astype(str).unique())}", flush=True)
    print(f"[load] conditions: {sorted(a.obs[CONDITION_KEY].astype(str).unique())}", flush=True)
    print(f"[load] {CELL_TYPE_KEY} cats: {a.obs[CELL_TYPE_KEY].nunique()}", flush=True)
    return a


def build_gp_dict() -> dict:
    GP_DATA_DIR.mkdir(parents=True, exist_ok=True)
    print("[gp] loading omnipath / nichenet / mebocost gene programs", flush=True)

    omnipath_csv = GP_DATA_DIR / "omnipath_lr_network.csv"
    # Reuse cached CSV if present (matches notebook behaviour).
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species="human",
        load_from_disk=omnipath_csv.exists(),
        save_to_disk=not omnipath_csv.exists(),
        lr_network_file_path=str(omnipath_csv),
        plot_gp_gene_count_distributions=False,
    )
    print(f"[gp] omnipath gps: {len(omnipath_gp_dict)}", flush=True)

    nichenet_files = {
        "lr_network_file_path": str(GP_DATA_DIR / "nichenet_lr_network_v2_human.csv"),
        "ligand_target_matrix_file_path": str(GP_DATA_DIR / "nichenet_ligand_target_matrix_v2_human.csv"),
    }
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species="human",
        version="v2",
        keep_target_genes_ratio=1.0,
        max_n_target_genes_per_gp=250,
        load_from_disk=all(Path(p).exists() for p in nichenet_files.values()),
        save_to_disk=not all(Path(p).exists() for p in nichenet_files.values()),
        plot_gp_gene_count_distributions=False,
        **nichenet_files,
    )
    print(f"[gp] nichenet gps: {len(nichenet_gp_dict)}", flush=True)

    mebocost_dir = GP_DATA_DIR / "metabolite_enzyme_sensor_gps"
    mebocost_dir.mkdir(parents=True, exist_ok=True)
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        species="human",
        dir_path=str(mebocost_dir),
        plot_gp_gene_count_distributions=False,
    )
    print(f"[gp] mebocost gps: {len(mebocost_gp_dict)}", flush=True)

    combined = filter_and_combine_gp_dict_gps_v2(
        gp_dicts=[omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict],
        verbose=False,
    )
    print(f"[gp] combined filtered gps: {len(combined)}", flush=True)
    return combined


def build_spatial_graph(adata: AnnData) -> None:
    print(f"[graph] spatial_neighbors(n_neighs={N_NEIGHBORS}, coord_type=generic)", flush=True)
    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        spatial_key=SPATIAL_KEY,
        n_neighs=N_NEIGHBORS,
    )
    # symmetrise
    adata.obsp[ADJ_KEY] = adata.obsp[ADJ_KEY].maximum(adata.obsp[ADJ_KEY].T)


def fit_nichecompass(adata: AnnData, gp_dict: dict) -> NicheCompass:
    print("[nc] add_gps_from_gp_dict_to_adata", flush=True)
    add_gps_from_gp_dict_to_adata(
        gp_dict=gp_dict,
        adata=adata,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        gp_names_key=GP_NAMES_KEY,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
    )

    print("[nc] init model", flush=True)
    model = NicheCompass(
        adata,
        counts_key=COUNTS_KEY,
        adj_key=ADJ_KEY,
        gp_names_key=GP_NAMES_KEY,
        active_gp_names_key=ACTIVE_GP_NAMES_KEY,
        gp_targets_mask_key=GP_TARGETS_MASK_KEY,
        gp_targets_categories_mask_key=GP_TARGETS_CATEGORIES_MASK_KEY,
        gp_sources_mask_key=GP_SOURCES_MASK_KEY,
        gp_sources_categories_mask_key=GP_SOURCES_CATEGORIES_MASK_KEY,
        latent_key=LATENT_KEY,
        conv_layer_encoder="gatv2conv",
        active_gp_thresh_ratio=0.01,
    )

    print(f"[nc] train  n_epochs={N_EPOCHS}  n_epochs_all_gps={N_EPOCHS_ALL_GPS}", flush=True)
    model.train(
        n_epochs=N_EPOCHS,
        n_epochs_all_gps=N_EPOCHS_ALL_GPS,
        lr=LR,
        lambda_edge_recon=LAMBDA_EDGE_RECON,
        lambda_gene_expr_recon=LAMBDA_GENE_EXPR_RECON,
        lambda_l1_masked=LAMBDA_L1_MASKED,
        lambda_l1_addon=LAMBDA_L1_ADDON,
        edge_batch_size=EDGE_BATCH_SIZE,
        n_sampled_neighbors=N_SAMPLED_NEIGHBORS,
        use_cuda_if_available=True,
        verbose=False,
    )
    return model


# ----- improved composition plot -----------------------------------------
def plot_niche_composition(
    adata: AnnData,
    niche_key: str,
    cat_key: str,
    out_path_base: Path,
    title_suffix: str = "",
    max_legend_cols: int = 1,
) -> None:
    """Side-by-side stacked bars: absolute counts + proportions, with
    per-niche totals annotated above each bar."""
    df = adata.obs[[niche_key, cat_key]].astype(str).copy()
    counts = (df.groupby([niche_key, cat_key]).size().unstack(fill_value=0))
    # order niches by total cell count desc, categories by global freq desc
    niche_totals = counts.sum(axis=1).sort_values(ascending=False)
    counts = counts.loc[niche_totals.index]
    cat_order = counts.sum(axis=0).sort_values(ascending=False).index
    counts = counts[cat_order]
    props = counts.div(counts.sum(axis=1), axis=0)

    palette = plt.cm.tab20(np.linspace(0, 1, max(len(cat_order), 1)))
    colour_map = {c: palette[i % len(palette)] for i, c in enumerate(cat_order)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)

    # --- absolute counts ---
    bottom = np.zeros(len(counts.index))
    xpos = np.arange(len(counts.index))
    for cat in cat_order:
        vals = counts[cat].values
        axes[0].bar(xpos, vals, bottom=bottom, color=colour_map[cat], label=cat,
                    width=0.78, edgecolor="white", linewidth=0.4)
        bottom += vals
    # annotate totals above each bar
    totals = counts.sum(axis=1).values
    tmax = totals.max() if len(totals) else 1
    for x, t in zip(xpos, totals):
        axes[0].text(x, t + tmax * 0.015, f"n={int(t):,}",
                     ha="center", va="bottom", fontsize=8)
    axes[0].set_ylim(0, tmax * 1.12)
    axes[0].set_ylabel("Cell count")
    axes[0].set_title(f"Absolute counts per niche{title_suffix}")

    # --- proportions ---
    bottom = np.zeros(len(props.index))
    for cat in cat_order:
        vals = props[cat].values
        axes[1].bar(xpos, vals, bottom=bottom, color=colour_map[cat],
                    width=0.78, edgecolor="white", linewidth=0.4)
        bottom += vals
    axes[1].set_ylim(0, 1.0)
    axes[1].set_ylabel("Proportion")
    axes[1].set_title(f"Proportions per niche{title_suffix}")
    # totals under proportion bars for cross-reference
    for x, t in zip(xpos, totals):
        axes[1].text(x, -0.035, f"n={int(t):,}", ha="center", va="top",
                     fontsize=7, color="#555")

    for ax in axes:
        ax.set_xticks(xpos)
        ax.set_xticklabels([str(n) for n in counts.index], rotation=0)
        ax.set_xlabel(f"Niche ({niche_key})")
        ax.grid(axis="y", linestyle=":", alpha=0.4)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colour_map[c]) for c in cat_order]
    fig.legend(handles, list(cat_order), title=cat_key,
               loc="center left", bbox_to_anchor=(1.01, 0.5),
               ncol=max_legend_cols, frameon=False, fontsize=8)

    fig.suptitle(f"Niche composition by {cat_key}", fontsize=14, y=1.02)
    fig.tight_layout()
    base = str(out_path_base)  # do NOT use Path.with_suffix: the resolution "0.4" in the name breaks it
    for ext in ("svg", "png"):
        p = f"{base}.{ext}"
        fig.savefig(p, bbox_inches="tight", dpi=160 if ext == "png" else None)
        print(f"  saved {p}", flush=True)
    plt.close(fig)

    out_csv = f"{base}.csv"
    out_csv_prop = f"{base}_proportions.csv"
    counts.to_csv(out_csv)
    props.round(4).to_csv(out_csv_prop)
    print(f"  saved {out_csv}", flush=True)
    print(f"  saved {out_csv_prop}", flush=True)


# ----- spatial / umap niche panels ---------------------------------------
def plot_niches_space(adata: AnnData, niche_key: str, fig_dir: Path, spot_size: float) -> None:
    latent_cluster_colors = create_new_color_dict(adata=adata, cat_key=niche_key)
    samples = sorted(adata.obs[SAMPLE_KEY].astype(str).unique())

    fig = plt.figure(figsize=(16, 14))
    title = fig.suptitle(t=f"NicheCompass niches in latent + physical space",
                         y=0.96, fontsize=18)
    spec1 = gridspec.GridSpec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[3, 2])
    spec2 = gridspec.GridSpec(ncols=len(samples), nrows=2,
                              width_ratios=[1] * len(samples), height_ratios=[3, 2])

    axs = [fig.add_subplot(spec1[0])]
    sc.pl.umap(adata=adata, color=[niche_key], palette=latent_cluster_colors,
               title="Niches (latent UMAP)", ax=axs[0], show=False, size=2)

    for idx, s in enumerate(samples):
        axs.append(fig.add_subplot(spec2[len(samples) + idx]))
        sub = adata[adata.obs[SAMPLE_KEY].astype(str) == s]
        sc.pl.spatial(adata=sub, color=[niche_key], palette=latent_cluster_colors,
                      spot_size=spot_size, title=s.replace("XETG00229__", ""),
                      legend_loc=None, ax=axs[idx + 1], show=False)
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.98, 0.5))
    axs[0].get_legend().remove()
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    out = fig_dir / f"res_{LEIDEN_RES}_niches_latent_physical_space"
    for ext in ("svg", "png"):
        fig.savefig(f"{out}.{ext}",
                    bbox_extra_artists=(lgd, title), bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


def plot_celltypes_space(adata: AnnData, cell_type_key: str, fig_dir: Path,
                         spot_size: float) -> None:
    cell_type_colors = create_new_color_dict(adata=adata, cat_key=cell_type_key)
    samples = sorted(adata.obs[SAMPLE_KEY].astype(str).unique())

    fig = plt.figure(figsize=(16, 14))
    title = fig.suptitle(t=f"Cell types ({cell_type_key}) in latent + physical space",
                         y=0.96, fontsize=18)
    spec1 = gridspec.GridSpec(ncols=1, nrows=2, width_ratios=[1], height_ratios=[3, 2])
    spec2 = gridspec.GridSpec(ncols=len(samples), nrows=2,
                              width_ratios=[1] * len(samples), height_ratios=[3, 2])
    axs = [fig.add_subplot(spec1[0])]
    sc.pl.umap(adata=adata, color=[cell_type_key], palette=cell_type_colors,
               title="Cell types (latent UMAP)", ax=axs[0], show=False, size=2)
    for idx, s in enumerate(samples):
        axs.append(fig.add_subplot(spec2[len(samples) + idx]))
        sub = adata[adata.obs[SAMPLE_KEY].astype(str) == s]
        sc.pl.spatial(adata=sub, color=[cell_type_key], palette=cell_type_colors,
                      spot_size=spot_size, title=s.replace("XETG00229__", ""),
                      legend_loc=None, ax=axs[idx + 1], show=False)
    handles, labels = axs[0].get_legend_handles_labels()
    lgd = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.98, 0.5))
    axs[0].get_legend().remove()
    plt.subplots_adjust(wspace=0.2, hspace=0.25)
    out = fig_dir / "cell_types_latent_physical_space"
    for ext in ("svg", "png"):
        fig.savefig(f"{out}.{ext}",
                    bbox_extra_artists=(lgd, title), bbox_inches="tight",
                    dpi=160 if ext == "png" else None)
    plt.close(fig)


# ----- main --------------------------------------------------------------
def main() -> int:
    global LEIDEN_RES, N_EPOCHS
    ap = argparse.ArgumentParser()
    ap.add_argument("--leiden-res", type=float, default=LEIDEN_RES)
    ap.add_argument("--n-epochs", type=int, default=N_EPOCHS)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    LEIDEN_RES = args.leiden_res
    N_EPOCHS = args.n_epochs

    ts = datetime.now().strftime("%d%m%Y_%H%M%S")
    run_dir = ARTIFACTS_DIR / ts
    model_dir = run_dir / "model"
    fig_dir = run_dir / "figures"
    model_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run] timestamp={ts}", flush=True)
    print(f"[run] model_dir={model_dir}", flush=True)
    print(f"[run] fig_dir={fig_dir}", flush=True)

    adata = load_adata()
    if args.smoke:
        rng = np.random.default_rng(0)
        keep = rng.choice(adata.n_obs, min(10_000, adata.n_obs), replace=False)
        adata = adata[keep].copy()
        N_EPOCHS = 3
        print(f"[smoke] adata -> {adata.shape}, n_epochs={N_EPOCHS}", flush=True)

    gp_dict = build_gp_dict()
    build_spatial_graph(adata)
    model = fit_nichecompass(adata, gp_dict)
    model.save(dir_path=str(model_dir), overwrite=True,
               save_adata=True, adata_file_name="adata.h5ad")
    load_timestamp = ts
    print(f"[nc] saved model + adata to {model_dir}", flush=True)

    # neighbours + umap on latent + niche leiden
    latent_cluster_key = f"latent_leiden_{LEIDEN_RES}"
    sc.pp.neighbors(model.adata, use_rep=LATENT_KEY, key_added=LATENT_KEY)
    sc.tl.umap(model.adata, neighbors_key=LATENT_KEY)
    sc.tl.leiden(adata=model.adata, resolution=LEIDEN_RES,
                 key_added=latent_cluster_key, neighbors_key=LATENT_KEY)
    print(f"[leiden] niches: {model.adata.obs[latent_cluster_key].value_counts().to_dict()}", flush=True)

    # compute a sensible spot_size from spatial extent
    xy = np.asarray(model.adata.obsm[SPATIAL_KEY])
    extent = max(xy[:, 0].ptp(), xy[:, 1].ptp())
    spot_size = max(8.0, float(extent) / 400.0)

    # latent / spatial panels
    plot_celltypes_space(model.adata, CELL_TYPE_KEY, fig_dir, spot_size)
    plot_niches_space(model.adata, latent_cluster_key, fig_dir, spot_size)

    # improved composition plots with absolute counts + proportions
    for cat_key, leg_cols in [
        (CELL_TYPE_KEY, 1),
        (CONDITION_KEY, 1),
        (SAMPLE_KEY, 1),
    ]:
        print(f"\n[composition] niche x {cat_key}", flush=True)
        plot_niche_composition(
            model.adata,
            niche_key=latent_cluster_key,
            cat_key=cat_key,
            out_path_base=fig_dir / f"res_{LEIDEN_RES}_niche_composition_{cat_key}",
            title_suffix=f"  (n={model.adata.n_obs:,} cells)",
            max_legend_cols=leg_cols,
        )

    # --- active GP report ---
    active_gps = model.get_active_gps()
    print(f"\n[gps] total {len(model.adata.uns[GP_NAMES_KEY])}, active {len(active_gps)}", flush=True)
    gp_summary_df = model.get_gp_summary()
    gp_summary_df.to_csv(fig_dir / "gp_summary.csv", index=False)

    # differential GP test: each niche vs rest
    log_bf = 2.3
    diff_key = "nichecompass_differential_gp_test_results"
    diff_done = False
    try:
        diff_gps = model.run_differential_gp_tests(
            cat_key=latent_cluster_key,
            selected_cats=None,
            comparison_cats="rest",
            log_bayes_factor_thresh=log_bf,
            key_added=diff_key,
        )
        pd.DataFrame(diff_gps).to_csv(
            fig_dir / f"log_bayes_factor_{log_bf}_niche_enriched_gps_summary.csv",
            index=False,
        )
        diff_done = True
    except Exception as exc:
        print(f"[diff_gps] skipped: {exc!r}", flush=True)
        traceback.print_exc()

    # enriched GP info plots (top 10) — depends on diff_gps results
    if diff_done:
        try:
            samples = sorted(model.adata.obs[SAMPLE_KEY].astype(str).unique())
            generate_enriched_gp_info_plots(
                model=model,
                plot_label=f"res_{LEIDEN_RES}",
                sample_key=SAMPLE_KEY,
                differential_gp_test_results_key=diff_key,
                cat_key=latent_cluster_key,
                cat_palette=create_new_color_dict(adata=model.adata, cat_key=latent_cluster_key),
                n_top_enriched_gp_start_idx=0,
                n_top_enriched_gp_end_idx=10,
                feature_spaces=samples,
                n_top_genes_per_gp=3,
                n_top_peaks_per_gp=0,
                save_figs=True,
                figure_folder_path=str(fig_dir),
                spot_size=spot_size,
            )
        except Exception as exc:
            print(f"[gp_info_plots] skipped: {exc!r}", flush=True)
            traceback.print_exc()

    # CCK network: only if the GP actually exists post-training
    try:
        gp_names = list(model.adata.uns[GP_NAMES_KEY])
        target_gp = "CCK_combined_GP"
        if target_gp in gp_names:
            network_df = compute_communication_gp_network(
                gp_list=[target_gp], model=model,
                group_key=latent_cluster_key,
            )
            visualize_communication_gp_network(
                adata=model.adata, network_df=network_df,
                figsize=(10, 10),
                save=True, file_path=str(fig_dir / f"gp_network_{target_gp}.svg"),
            )
        else:
            print(f"[cck] GP '{target_gp}' not in gp_names; skipping", flush=True)
    except Exception as exc:
        print(f"[cck] skipped: {exc!r}", flush=True)

    print(f"\n=== DONE ts={load_timestamp} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
