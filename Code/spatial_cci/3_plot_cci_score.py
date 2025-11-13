#
# Created on Sun Sep 21 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script computes and visualizes local cell–cell interaction (CCI) patterns from Xenium
# kidney biopsy spatial transcriptomics data. It integrates precomputed Voronoi neighborhoods,
# ligand–receptor (LR) interaction scores, and latent variable Factor2_ulm to identify spatial
# and functional dependencies across cell types.
#
# Datasets used/loaded:
# - Cell and nucleus boundaries: geodataframes/cell_boundaries.parquet, geodataframes/nucleus_boundaries.parquet
# - Processed Xenium object: xenium/xenium_processed.h5ad
# - Voronoi-based neighborhood graphs: neighborhoods/voronoi_graphs_10.pkl
# - Local ligand–receptor interaction scores: local_cci_results/cci_score_df.parquet
#
# Objects generated:
# - 2-hop adjacency dictionaries (adj_2hop_dict) for all tissue cores
# - Subsetted and harmonized AnnData object (adata)
# - Mean ligand–receptor score table (cci_avg_df)
# - LOESS smoothing cache objects (vs_factor2/loess_cache_<param>.pkl, vs_mean_factor2/loess_cache_<param>.pkl)
# - Correlation matrices of LR scores vs Factor2_ulm (factor2_corr/correlations.pkl)
#
# Datasets saved:
# - loess_cache_<param>.pkl (smoothed interaction–Factor2 trends)
# - correlations.pkl (Spearman correlation between interaction score and Factor2_ulm)
#
# Plots generated (by output folder):
#
# lr_heatmap/:
#   - Heatmaps of mean ligand–receptor interaction scores per cell-type pair
#   - Axes: ligands (rows), receptors (columns)
#   - Titles: "<source> → <target>"
#
# vs_factor2/:
#   - KDE plots of Factor2_ulm distribution across neighborhoods per cell-type pair
#   - Smoothed (LOESS) ligand–receptor interaction curves vs Factor2_ulm
#   - Figures labeled by cell-type pair (aPT → Immune, Immune → Fib/MyoFib, aPT → Fib/MyoFib)
#
# vs_mean_factor2/:
#   - LOESS smoothed ligand–receptor interaction curves vs mean Factor2_ulm
#
# factor2_corr/:
#   - Heatmap of correlation coefficients (score vs Factor2_ulm) across all LR and CT pairs
#   - Clustered correlation heatmap (Ward linkage on rows and columns)
#   - Individual heatmaps per LR pair showing correlation patterns across cell-type pairs
#
# The script performs:
# 1. Harmonization of cell annotations and data structures across spatial and transcriptomic data.
# 2. Construction of 2-hop cell neighborhood graphs and adjacency matrices.
# 3. Selection or exclusion of cell-type pairs of interest for targeted analyses.
# 4. Computation and visualization of mean LR interaction strengths.
# 5. Estimation of smoothed relationships between LR scores and Factor2_ulm using LOESS.
# 6. Quantification of correlations between Factor2_ulm and interaction scores, with clustering.
# --------------------------------------------------------------------------------------------


import re
import scanpy as sc
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import seaborn as sns
import os
import numpy as np
import geopandas as gpd
import networkx as nx
from matplotlib.colors import ListedColormap
import pickle
import warnings
import sys
import argparse
warnings.filterwarnings(
    "ignore", message="nopython is set for njit and is ignored")

# Load helper functions ----------------------------
sys.path.append(os.path.abspath("../utils/"))
from running_on_slurm import running_on_slurm
from spatial_cci_helpers import *

# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/spatial_cci/"


# -------------------------------------------------------------------------------
# Load and prepare data ---------------------------------------------------------
# -------------------------------------------------------------------------------

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(
    res_dir+"geodataframes/nucleus_boundaries.parquet")

# Replace aPT/PEC with aPT in gdf_cell_boundaries using cat.rename_categories
gdf_cell_boundaries["cell_type_manual"] = gdf_cell_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)
gdf_nucleus_boundaries["cell_type_manual"] = gdf_nucleus_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)

# Get f2 mapping and other info from gdf_cell_boundaries ----------------
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()
cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories.tolist()
factor2_map = gdf_cell_boundaries.set_index("sample_cell_id")["Factor2_ulm"]

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")


# Load adata ---------------------------------------------------
adata = sc.read_h5ad(res_dir + "xenium/xenium_processed.h5ad")

adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(
    str) + "_" + adata.obs["cell_id"].astype(str)

# Print cores in adata.obs["TMA_ID"] and not in cores -----
set(adata.obs["TMA_ID"].unique().tolist()) - set(cores)
# {'Utrecht 2 A5 T03-04110 K29 MN',
#  'Utrecht1 B5 T18-01374 K25 TIN',
#  'Utrecht2 C3 T17-11176 K24 TIN'}

# Remove them from adata ------
adata = adata[adata.obs["TMA_ID"].isin(cores), :].copy()

# Replace aPT/PEC with aPT in adata.obs["cell_type_manual"]
adata.obs["cell_type_manual"] = adata.obs["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)

# Load Voronoi graph ---------------------------------------------------
max_threshold = 10
with open(f"{res_dir}neighborhoods/voronoi_graphs_{max_threshold}.pkl", "rb") as f:
    voronoi_graphs = pickle.load(f)

# Create adj list ---------------------------------------------------
adj_2hop_list = {
    core: {
        n: list(
            set(nx.single_source_shortest_path_length(
                voronoi_graphs[core], n, cutoff=2).keys()) - {n}
        )
        for n in voronoi_graphs[core].nodes()
    }
    for core in cores
}
adj_2hop_dict = {
    root: neighbors
    for core in adj_2hop_list.keys()
    for root, neighbors in adj_2hop_list[core].items()
}

del adj_2hop_list, voronoi_graphs

# Load script parameters ---------------------------------------------------

if running_on_slurm():
    end_root = 198361
    # Parse parameter called "cts": either "ct_of_interest", "rest", or "all"
    parser = argparse.ArgumentParser(description="Plot local CCI results")
    parser.add_argument("--cts", type=str, default="all", help="Cell types parameter: ct_of_interest, rest, or all")
    # Load parameter plots_to_make: comma-separated list of plots to make
    parser.add_argument(
        "--plots_to_make", type=str,
        default="lr_heatmap,vs_factor2,factor2_corr",
        help="Plots to make: comma-separated list"
    )
    args = parser.parse_args()
    param_cts = args.cts
    plots_to_make = args.plots_to_make.split(",")

else:
    end_root = 25000
    param_cts = "ct_of_interest"
    plots_to_make = ["lr_heatmap", "vs_factor2", "factor2_corr"]

print("Running on SLURM:", running_on_slurm())

# Print all parameters
print("Parameters:")
print(" - Cell types:", param_cts)
print(" - Plots to make:", plots_to_make)


# Load CCI results ---------------------------------------------------

cci_score_df = pd.read_parquet(
    fig_dir + "local_cci_results/cci_score_df.parquet"
)
assert cci_score_df["root"].nunique() == 198361

# Subset to the number of unique roots specified by "end"
cci_score_df = cci_score_df[
    cci_score_df["root"].isin(
        cci_score_df["root"].unique().tolist()[:end_root]
    )
]

print("Number of roots:", cci_score_df["root"].nunique())


# Subsetting for cell types of interest if relevant -----------------
ct_of_interest = ["aPT", "Immune", "Fib/MyoFib"]
cell_type_pairs = [('aPT', 'Immune'),
                   ('Immune', 'Fib/MyoFib'),
                   ('aPT', 'Fib/MyoFib')]

if param_cts != "all":
    mask = (((cci_score_df["source"] == "aPT")
            & (cci_score_df["target"] == "Immune")) |
            ((cci_score_df["source"] == "Immune")
             & (cci_score_df["target"] == "Fib/MyoFib")) |
            ((cci_score_df["source"] == "aPT") &
             (cci_score_df["target"] == "Fib/MyoFib")))
    if param_cts == "ct_of_interest":
        cci_score_df = cci_score_df[mask].copy()
    elif param_cts == "rest":
        cci_score_df = cci_score_df[~mask].copy()

# Compute average
cci_avg_df = cci_score_df.groupby(['source', 'target', 'gene_a', 'gene_b']).agg({
    'score': 'mean'}).reset_index()


# 1. Plot heamaps results exclusively for cell types of interest ----------------------------------------------

if "lr_heatmap" in plots_to_make:

    # Filtering average values -----------------
    score_thresh = 0.1
    cci_avg_df.loc[cci_avg_df["score"] < score_thresh, "score"] = 0.0
    # Generalize heatmap computation for all ct_pairs in cci_score_df

    heatmaps = {}
    ct_pairs = cci_score_df["ct_pair"].unique().tolist()

    for ct_pair in ct_pairs:
        # Compute heatmaps to plots
        src, tgt = ct_pair.split("_to_")
        df = cci_avg_df[(cci_avg_df["source"] == src) & (cci_avg_df["target"] == tgt)].copy()
        if df.empty:
            continue
        ligands = df["gene_a"].unique().tolist()
        receptors = df["gene_b"].unique().tolist()
        heatmap = pivot_cci_score_df(df, ligands, receptors)
        # Remove rows and columns with all zeros
        heatmap = heatmap.loc[~(heatmap == 0).all(axis=1), :]
        heatmap = heatmap.loc[:, ~(heatmap == 0).all(axis=0)]
        # Set zero values to NaN for better visualization
        heatmap = heatmap.replace(0, np.nan)
        heatmaps[ct_pair] = heatmap
        figsize = (10, 8)
        # Plot heatmap
        title = f"{src} → {tgt}"
        fname = f"cci_heatmap_{src.replace('/', '-')}_to_{tgt.replace('/', '-')}.png"

        plt.figure(figsize=figsize)
        sns.heatmap(heatmap, cmap="YlGnBu", annot=False)
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0, ha="right")
        plt.colorbar(label="Local CCI Score")
        plt.title(title)
        plt.xlabel("Receptor")
        plt.ylabel("Ligand")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, "lr_heatmap", fname), dpi=300)
        plt.show()


# 2. Plot vs Factor2_ulm ---------

# 2.1 Visualize each neighborhood's factor2 density -------

out_dir = os.path.join(fig_dir, "vs_factor2", "nh_factor2_density")

for ct_pair in cci_score_df["ct_pair"].unique():

    sub = cci_score_df[cci_score_df["ct_pair"] == ct_pair]
    # keep one entry per unique root
    uniq = sub.drop_duplicates(subset="root")
    vals = uniq["Factor2_ulm"].dropna().astype(float)

    if vals.size < 2:
        print(f"Skipping {ct_pair}: not enough unique roots with Factor2_ulm")
        continue

    plt.figure(figsize=(6, 4))
    sns.kdeplot(vals, fill=True, bw_method="scott", color="C0", alpha=0.7)
    sns.rugplot(vals, height=0.05, color="k", alpha=0.3)

    plt.xlim(-4, 10)
    plt.xlabel("Factor2_ulm")
    plt.ylabel("Density")
    plt.title(f"Density of Factor2_ulm (unique roots)\n{ct_pair.replace('_to_', ' → ')}")
    plt.tight_layout()

    safe_name = ct_pair.replace("/", "-").replace(" ", "_")
    out_fn = os.path.join(out_dir, f"neighborhood_factor2_density_{safe_name}.png")
    plt.savefig(out_fn, dpi=300, bbox_inches="tight")
    plt.close()

# 2.2 vs Factor2 for all lr pairs -----------------

# 2.2.1 Precompute Smooth Estimates -----------------

# LOESS estimate -----------------
if running_on_slurm() & (not os.path.exists(fig_dir + f"vs_factor2/loess_cache_{param_cts}.pkl")):
    loess_fraction = 0.4
    loess_cache = {}
    for pair_label in cci_score_df["ct_pair"].unique():
        df = cci_score_df[cci_score_df["ct_pair"] == pair_label].copy()
        loess_cache[pair_label] = {}
        for ga in df["gene_a"].unique():
            for gb in df["gene_b"].unique():
                sub = df[(df["gene_a"] == ga) & (df["gene_b"] == gb)]
                sub = sub.dropna(subset=["Factor2_ulm", "score"])
                if sub.empty:
                    continue
                agg = sub.groupby("Factor2_ulm")["score"].mean().sort_index()
                x, y = agg.index.values, agg.values
                _, y_s = compute_loess(x, y, frac=loess_fraction)
                # Remove all values of x and y where x is not in -4, 8
                mask = (x >= -4) & (x <= 8)
                x, y_s = x[mask], y_s[mask]
                loess_cache[pair_label][(ga, gb)] = (x, y_s)
    # Save loess_cache to pickle
    with open(fig_dir + f"vs_factor2/loess_cache_{param_cts}.pkl", "wb") as f:
        pickle.dump(loess_cache, f)

# Load loess_cache from pickle
with open(fig_dir + f"vs_factor2/loess_cache_{param_cts}.pkl", "rb") as f:
    loess_cache = pickle.load(f)



# LOESS estimate for mean_Factor2_ulm -----------------
loess_mean_path = os.path.join(
    fig_dir, "vs_mean_factor2", f"loess_cache_{param_cts}.pkl"
)
if running_on_slurm() & (not os.path.exists(loess_mean_path)):
    loess_fraction = 0.4
    loess_cache_mean = {}
    for pair_label in cci_score_df["ct_pair"].unique():
        df = cci_score_df[cci_score_df["ct_pair"] == pair_label].copy()
        loess_cache_mean[pair_label] = {}
        for ga in df["gene_a"].unique():
            for gb in df["gene_b"].unique():
                sub = df[(df["gene_a"] == ga) & (df["gene_b"] == gb)]
                sub = sub.dropna(subset=["mean_Factor2_ulm", "score"])
                if sub.empty:
                    continue
                agg = sub.groupby("mean_Factor2_ulm")["score"].mean().sort_index()
                x = agg.index.values.astype(float)
                y = agg.values.astype(float)
                if len(x) == 0:
                    continue
                try:
                    _, y_s = compute_loess(x, y, frac=loess_fraction)
                except Exception:
                    continue
                mask = (x >= -4) & (x <= 8)
                x_masked, y_s_masked = x[mask], y_s[mask]
                if len(x_masked) == 0:
                    continue
                loess_cache_mean[pair_label][(ga, gb)] = (x_masked, y_s_masked)
    with open(loess_mean_path, "wb") as f:
        pickle.dump(loess_cache_mean, f)

with open(loess_mean_path, "rb") as f:
    loess_cache_mean = pickle.load(f)


# 2.2.2 Plot LR score vs factor2: one figure per cell type pair -----------------
y_min_dist = 0.2
ymax_thresh_label = 0.2

# Plot Loess estimate -------------
ct = cell_type_pairs[2]
for ct in cell_type_pairs:

    ct_pair = f"{ct[0]}_to_{ct[1]}"
    palette = sns.color_palette("tab20", n_colors=len(loess_cache[ct_pair]))
    color_map = {k: palette[i] for i, k in enumerate(loess_cache[ct_pair].keys())}

    fig, ax = plt.subplots(figsize=(10, 6))

    # --- 1. Collect candidate label positions ---
    labels = []  # each entry = [y_text, x_text, gene_a, gene_b]
    for (gene_a, gene_b), (x_s, y_s) in loess_cache[ct_pair].items():

        # If y_s has more than one third of zeros, dont plot
        if np.sum(y_s == 0) > len(y_s) / 3:
            print(f"Skipping {gene_a} → {gene_b} due to too many zeros in loess")
        else:
            ax.plot(x_s, y_s, color=color_map[(gene_a, gene_b)], linewidth=2, alpha=0.9)

        if np.nanmax(np.abs(y_s)) >= 0.2:  # only label if relevant
            x_txt = max(x_s) * 1.04
            y_txt = y_s[-1]
            labels.append([y_txt, x_txt, gene_a, gene_b])

    # --- 2. Sort by y coordinate and enforce vertical spacing ---
    labels = sorted(labels, key=lambda x: x[0])
    adjusted = []
    for i, (y, x, ga, gb) in enumerate(labels):
        if i == 0:
            new_y = y
        else:
            prev_y = adjusted[-1][0]
            new_y = max(y, prev_y + y_min_dist)
        adjusted.append([new_y, x, ga, gb])

    # --- 3. Draw adjusted labels ---
    for (new_y, x_txt, ga, gb) in adjusted:
        ax.text(x_txt, new_y, f"{ga} → {gb}",
                va='center', ha='left', fontsize=8)

    # --- 4. Formatting and save ---
    plt.xlim(-4, 13)
    # Remove frame of the plot
    for spine in ax.spines.values():
        spine.set_visible(False)
    ct_pair_to_plot = ct_pair.replace("_to_", " → ")
    plt.title(f"{ct_pair_to_plot}", fontsize=10)
    plt.xlabel("Factor2")
    plt.ylabel("Interaction score (loess smoothing)")
    plt.tight_layout(rect=[0, 0, 0.78, 1])

    # --- 4b. Remove last x-axis tick safely ---
    xticks = ax.get_xticks()
    if len(xticks) > 0:
        new_ticks = xticks[:-2]
        ax.set_xticks(new_ticks)               # set tick positions
        ax.set_xticklabels([str(t) for t in new_ticks])  # set tick labels

    out_fn = os.path.join(
        fig_dir, 'vs_factor2',
        f"loess_{ct_pair.replace('/', '-').replace(' ', '_')}.png"
    )
    plt.savefig(out_fn, dpi=300, bbox_inches="tight")
    plt.show()


# 3. Compute correlation with Factor2_ulm -----------------

if "factor2_corr" in plots_to_make:

    # Compute the correlations if not already done
    if running_on_slurm() and (not os.path.exists(os.path.join(fig_dir, "factor2_corr", "correlations.pkl"))):
        corr_df = (
            cci_score_df.groupby(["ct_pair", "lr_pair"])
            .apply(lambda x: x["score"].corr(x["Factor2_ulm"], method="spearman"))
            .reset_index(name="corr_score_Factor2")
        )
        corr_df.to_pickle(os.path.join(fig_dir, "factor2_corr", "correlations.pkl"))
    else:
        corr_df = pd.read_pickle(os.path.join(fig_dir, "factor2_corr", "correlations.pkl"))

    corr_mat = corr_df.pivot(index='ct_pair', columns='lr_pair', values='corr_score_Factor2')

    # 3.1.1 Plot heatmap of correlation matrix -----------------

    fig, ax = plt.subplots(figsize=(12, 20))
    im = ax.imshow(corr_mat.values, aspect='auto', interpolation='nearest',
                   vmin=-1, vmax=1, cmap='RdYlGn')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Correlation (score vs Factor2_ulm)")
    ax.set_xticks(np.arange(corr_mat.shape[1]))
    ax.set_yticks(np.arange(corr_mat.shape[0]))
    ax.set_xticklabels(corr_mat.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(corr_mat.index, fontsize=8)
    ax.set_xlabel("Ligand–Receptor Pair (lr_pair)", fontsize=10)
    ax.set_ylabel("Cell-Type Pair (ct_pair)", fontsize=10)
    ax.set_title("Correlation between score and Factor2_ulm", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "factor2_corr", "correlation_heatmap.png"))
    plt.show()

    # 3.1.2 Plot clustered correlation matrix -----------------

    corr_mat_filled = corr_mat.fillna(0)

    # a) Cluster rows (ct_pair)
    # pairwise distances
    row_dist = pdist(corr_mat_filled.values, metric='euclidean')
    # hierarchical clustering
    row_linkage = linkage(row_dist, method='ward')
    # optimal row ordering
    row_order = leaves_list(row_linkage)
    corr_mat_clustered = corr_mat_filled.iloc[row_order]

    # b) Cluster columns (lr_pair)
    col_dist = pdist(corr_mat_clustered.values.T, metric='euclidean')
    col_linkage = linkage(col_dist, method='ward')
    col_order = leaves_list(col_linkage)
    corr_mat_clustered = corr_mat_clustered.iloc[:, col_order]

    # c) Plot clustered heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(corr_mat_clustered.values,
                   aspect='auto', interpolation='nearest',
                   vmin=-1, vmax=1, cmap='RdYlGn')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Correlation (score vs Factor2_ulm)")
    ax.set_xticks(np.arange(corr_mat_clustered.shape[1]))
    ax.set_yticks(np.arange(corr_mat_clustered.shape[0]))
    ax.set_xticklabels(corr_mat_clustered.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(corr_mat_clustered.index, fontsize=8)
    ax.set_xlabel("Ligand–Receptor Pair (lr_pair)", fontsize=10)
    ax.set_ylabel("Cell-Type Pair (ct_pair)", fontsize=10)
    ax.set_title("Correlation (clustered by hierarchical clustering)", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "factor2_corr",
                "correlation_heatmap_clustered.png"))
    plt.show()

    # 3.1.3 Plot one correlation matrix per lr_pair -----------------

    out_dir = os.path.join(fig_dir, "factor2_corr", "per_lr")

    corr_df["source"] = corr_df["ct_pair"].apply(lambda x: x.split("_to_")[0])
    corr_df["target"] = corr_df["ct_pair"].apply(lambda x: x.split("_to_")[1])

    all_sources = sorted(corr_df["source"].unique())
    all_targets = sorted(corr_df["target"].unique())

    for lr in corr_df["lr_pair"].unique():
        sub = corr_df[corr_df["lr_pair"] == lr][["source", "target", "corr_score_Factor2"]].copy()
        if sub.empty:
            continue
        mat = sub.pivot(index="source", columns="target", values="corr_score_Factor2")
        mat.shape
        if mat.isna().all().all():
            continue

        fig, ax = plt.subplots(figsize=(max(6, len(all_targets) * 0.4), max(4, len(all_sources) * 0.4)))
        im = ax.imshow(mat.values, aspect="auto", interpolation="nearest", vmin=-1, vmax=1, cmap="RdYlGn")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Correlation (score vs Factor2_ulm)")

        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(mat.index, fontsize=8)

        ax.set_xlabel("Target", fontsize=9)
        ax.set_ylabel("Source", fontsize=9)
        ax.set_title(f"Correlation matrix for {lr}", fontsize=10)
        plt.tight_layout()

        safe_lr = re.sub(r'[^0-9A-Za-z_.-]', '_', lr)
        out_fn = os.path.join(out_dir, f"corr_matrix_{safe_lr}.png")
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close(fig)


# 3.2 Compute correlation with mean_Factor2_ulm -----------------

if "factor2_corr" in plots_to_make:

    # Compute the correlations if not already done
    if running_on_slurm() and (not os.path.exists(os.path.join(mean_out_dir, "correlations.pkl"))):
        corr_mean_df = (
            cci_score_df.groupby(["ct_pair", "lr_pair"])
            .apply(lambda x: x["score"].corr(x["mean_Factor2_ulm"], method="spearman"))
            .reset_index(name="corr_score_mean_Factor2")
        )
        corr_mean_df.to_pickle(os.path.join(mean_out_dir, "correlations.pkl"))
    else:
        corr_mean_df = pd.read_pickle(os.path.join(mean_out_dir, "correlations.pkl"))

    corr_mat_mean = corr_mean_df.pivot(index='ct_pair', columns='lr_pair', values='corr_score_mean_Factor2')

    # 3.2.1 Plot heatmap of correlation matrix -----------------
    fig, ax = plt.subplots(figsize=(12, 20))
    im = ax.imshow(corr_mat_mean.values, aspect='auto', interpolation='nearest',
                   vmin=-1, vmax=1, cmap='RdYlGn')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Correlation (score vs mean_Factor2_ulm)")
    ax.set_xticks(np.arange(corr_mat_mean.shape[1]))
    ax.set_yticks(np.arange(corr_mat_mean.shape[0]))
    ax.set_xticklabels(corr_mat_mean.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(corr_mat_mean.index, fontsize=8)
    ax.set_xlabel("Ligand–Receptor Pair (lr_pair)", fontsize=10)
    ax.set_ylabel("Cell-Type Pair (ct_pair)", fontsize=10)
    ax.set_title("Correlation between score and mean_Factor2_ulm", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(mean_out_dir, "correlation_heatmap.png"))
    plt.show()

    # 3.2.2 Plot clustered correlation matrix -----------------
    corr_mat_mean_filled = corr_mat_mean.fillna(0)

    # a) Cluster rows (ct_pair)
    row_dist = pdist(corr_mat_mean_filled.values, metric='euclidean')
    row_linkage = linkage(row_dist, method='ward')
    row_order = leaves_list(row_linkage)
    corr_mat_mean_clustered = corr_mat_mean_filled.iloc[row_order]

    # b) Cluster columns (lr_pair)
    col_dist = pdist(corr_mat_mean_clustered.values.T, metric='euclidean')
    col_linkage = linkage(col_dist, method='ward')
    col_order = leaves_list(col_linkage)
    corr_mat_mean_clustered = corr_mat_mean_clustered.iloc[:, col_order]

    # c) Plot clustered heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(corr_mat_mean_clustered.values, aspect='auto', interpolation='nearest',
                   vmin=-1, vmax=1, cmap='RdYlGn')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Correlation (score vs mean_Factor2_ulm)")
    ax.set_xticks(np.arange(corr_mat_mean_clustered.shape[1]))
    ax.set_yticks(np.arange(corr_mat_mean_clustered.shape[0]))
    ax.set_xticklabels(corr_mat_mean_clustered.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(corr_mat_mean_clustered.index, fontsize=8)
    ax.set_xlabel("Ligand–Receptor Pair (lr_pair)", fontsize=10)
    ax.set_ylabel("Cell-Type Pair (ct_pair)", fontsize=10)
    ax.set_title("Correlation (clustered by hierarchical clustering) - mean_Factor2_ulm", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(mean_out_dir, "correlation_heatmap_clustered.png"))
    plt.show()

    # 3.2.3 Plot one correlation matrix per lr_pair -----------------
    per_lr_out = os.path.join(mean_out_dir, "per_lr")
    os.makedirs(per_lr_out, exist_ok=True)

    corr_mean_df["source"] = corr_mean_df["ct_pair"].apply(lambda x: x.split("_to_")[0])
    corr_mean_df["target"] = corr_mean_df["ct_pair"].apply(lambda x: x.split("_to_")[1])

    all_sources_mean = sorted(corr_mean_df["source"].unique())
    all_targets_mean = sorted(corr_mean_df["target"].unique())

    for lr in corr_mean_df["lr_pair"].unique():
        sub = corr_mean_df[corr_mean_df["lr_pair"] == lr][["source", "target", "corr_score_mean_Factor2"]].copy()
        if sub.empty:
            continue
        mat = sub.pivot(index="source", columns="target", values="corr_score_mean_Factor2")
        if mat.isna().all().all():
            continue

        fig, ax = plt.subplots(figsize=(max(6, len(all_targets_mean) * 0.4), max(4, len(all_sources_mean) * 0.4)))
        im = ax.imshow(mat.values, aspect="auto", interpolation="nearest", vmin=-1, vmax=1, cmap="RdYlGn")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Correlation (score vs mean_Factor2_ulm)")

        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(mat.index, fontsize=8)

        ax.set_xlabel("Target", fontsize=9)
        ax.set_ylabel("Source", fontsize=9)
        ax.set_title(f"Correlation matrix for {lr} (mean_Factor2_ulm)", fontsize=10)
        plt.tight_layout()

        safe_lr = re.sub(r'[^0-9A-Za-z_.-]', '_', lr)
        out_fn = os.path.join(per_lr_out, f"corr_matrix_{safe_lr}.png")
        plt.savefig(out_fn, dpi=300, bbox_inches="tight")
        plt.close(fig)
