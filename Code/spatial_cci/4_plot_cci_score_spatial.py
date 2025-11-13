#
# Created on Sun Oct 19 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script generates spatial expression maps of ligand–receptor (LR) pairs across cell-type
# interactions in Xenium kidney biopsy spatial transcriptomics data. It visualizes ligand and
# receptor expression across cell-type pairs and tissue cores, highlighting spatial colocalization
# patterns between interacting cell populations.
#
# Datasets used/loaded:
# - Cell and nucleus boundaries: geodataframes/cell_boundaries.parquet, geodataframes/nucleus_boundaries.parquet
# - Circularized cell outlines: geodataframes/cell_circles.parquet
# - Processed Xenium AnnData object: xenium/xenium_processed.h5ad
# - (Optional, precomputed) LR pair list: spatial_cci/lr_pairs.txt
# - (If missing) LR interaction scores: local_cci_results/all_cores/interaction_scores_*.parquet
#
# Datasets/objects generated:
# - lr_pairs.txt: unique ligand–receptor pairs derived from CCI score files
# - adata: filtered and harmonized AnnData object aligned with spatial boundaries
# - gdf_cell_boundaries / gdf_nucleus_boundaries: harmonized with renamed categories (aPT/PEC → aPT)
# - Expression tables per LR pair and tissue core, merged with spatial geometries
#
# Figures generated (by output folder):
#
# spatial_plot/wo_pt1_pt2/:
#   - For each LR pair and each core:
#       - 2×2 subplot figure showing expression of ligand (orange) and receptor (blue)
#         across three main interacting pairs: (aPT→Fib/Myofib), (aPT→Immune), (Immune→Fib/Myofib)
#       - Fourth panel: spatial map of cell-type annotation
#       - Overlaid nucleus boundaries for spatial context
#       - Colorbars for each ligand and receptor, with reversed y-axis for spatial alignment
#
# spatial_plot/with_pt1_pt2/:
#   - For each LR pair and each core:
#       - Grid of spatial plots for all interactions:
#           (aPT→Fib/Myofib), (aPT→Immune), (Immune→Fib/Myofib),
#           (PT1→Immune), (PT2→Immune)
#       - Ligand (orange, "Oranges" colormap) and receptor (blue, "Blues" colormap)
#         plotted together for each cell-type pair
#       - Dual stacked colorbars per figure (ligand on top, receptor below)
#       - Titles colored according to corresponding colormaps
#
# The script performs:
# 1. Loading and harmonization of spatial geometries, annotations, and expression data.
# 2. Generation of a unified ligand–receptor pair list if absent.
# 3. Filtering of AnnData and boundary tables to include consistent tissue cores.
# 4. Computation of per-cell ligand and receptor expression for each LR pair.
# 5. Spatial visualization of LR signal across tissue sections for specific cell-type pairs,
#    with and without proximal tubule subtypes (PT1, PT2).
# 6. Automated figure generation for all LR pairs and tissue cores with structured output folders.
# --------------------------------------------------------------------------------------------


import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import geopandas as gpd
from matplotlib.colors import ListedColormap
import warnings
import sys
from math import ceil
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
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
gdf_cell_circles = gpd.read_parquet(
    res_dir+"geodataframes/cell_circles.parquet")


# Replace aPT/PEC with aPT in gdf_cell_boundaries using cat.rename_categories
gdf_cell_boundaries["cell_type_manual"] = gdf_cell_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)
gdf_nucleus_boundaries["cell_type_manual"] = gdf_nucleus_boundaries["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)

# Unname index
gdf_cell_boundaries.index.name = None
gdf_nucleus_boundaries.index.name = None

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


# Replace aPT/PEC with aPT in adata.obs["cell_type_manual"]
adata.obs["cell_type_manual"] = adata.obs["cell_type_manual"].cat.rename_categories(
    {"aPT/PEC": "aPT"}
)


adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(
    str) + "_" + adata.obs["cell_id"].astype(str)

# Print cores in adata.obs["TMA_ID"] and not in cores -----
set(adata.obs["TMA_ID"].unique().tolist()) - set(cores)
# {'Utrecht 2 A5 T03-04110 K29 MN',
#  'Utrecht1 B5 T18-01374 K25 TIN',
#  'Utrecht2 C3 T17-11176 K24 TIN'}

# Remove them from adata ------
adata = adata[adata.obs["TMA_ID"].isin(cores), :].copy()


# If it does not exist yet, create and save lr_pairs.txt --------------
if not os.path.exists(fig_dir + "lr_pairs.txt"):

    # Load LR list from results ------
    if running_on_slurm():
        end = 198361
    else:
        end = 25000

    print("Running on SLURM:", running_on_slurm())
    print("End:", end)

    cci_score_df = load_continuous_chunks(
        "interaction_scores",
        folder=fig_dir + "local_cci_results/all_cores/",
        start=0, end=end
    )

    # convert source and target as strings in cci_score_df
    cci_score_df["lr_pair"] = cci_score_df["gene_a"] + "_" + cci_score_df["gene_b"]

    # Save cci_score_df["lr_pair"].unique() as a list to text
    with open("lr_pairs.txt", "w") as f:
        for lr_pair in cci_score_df["lr_pair"].unique():
            f.write(lr_pair + "\n")
    lr_pairs = list(reversed(cci_score_df["lr_pair"].unique().tolist()))

else:
    with open(fig_dir + "lr_pairs.txt", "r") as f:
        lr_pairs = [line.strip() for line in f.readlines()]


# Update adata.obs_names to sample_cell_id ------------------------
adata.obs["sample_cell_id"] = adata.obs["TMA_ID"].astype(
    str) + "_" + adata.obs["cell_id"].astype(str)
adata.obs_names = adata.obs["sample_cell_id"].astype(str)


# 4. Spatial plot -------------

# 4.1 One figure per LR pair and per core (wo PT1 and PT2) ---------

out_dir_wo_pt1_pt2 = os.path.join(fig_dir, "spatial_plot", "wo_pt1_pt2")
os.makedirs(out_dir_wo_pt1_pt2, exist_ok=True)

ct_pairs = [("aPT", "Fib/Myofib"),
            ("aPT", "Immune"),
            ("Immune", "Fib/Myofib")]
unique_cts = {ct for pair in ct_pairs for ct in pair}

ct_pair = ct_pairs[0]
lr_pair = "APP_CD74"
core = cores[0]

for lr_pair in lr_pairs:
    gene_a, gene_b = lr_pair.split("_")

    # Compute expression once per lr_pair
    genes = [gene_a, gene_b]
    expr = adata[:, genes].to_df()
    expr[lr_pair] = expr[gene_a] * expr[gene_b]
    expr.reset_index(inplace=True)  # contains sample_cell_id

    for core in cores:

        gdf_sub = gdf_cell_boundaries.loc[gdf_cell_boundaries.TMA_ID == core]
        gdf_sub_nuc = gdf_nucleus_boundaries.loc[gdf_nucleus_boundaries.TMA_ID == core]
        gdf_sub = gdf_sub.merge(expr, left_on="sample_cell_id",
                                right_on="sample_cell_id", how="left")

        # Precompute cell-type masks
        masks = {ct: (gdf_sub["cell_type_manual"] == ct) for ct in unique_cts}

        # Create figure
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
        axes = axes.flatten()
        # Compute signal to plot: L for source cells and R for target cells
        ct_pair = ct_pairs[0]
        for ax_ind, ct_pair in enumerate(ct_pairs):
            mA = masks[ct_pair[0]]
            mB = masks[ct_pair[1]]
            gdf_sub[f"lr_A_{ct_pair}"] = gdf_sub[gene_a].where(mA)
            gdf_sub[f"lr_B_{ct_pair}"] = gdf_sub[gene_b].where(mB)

            lr_ct_pair_dict = {ct_pair[0]: gene_a, ct_pair[1]: gene_b}
            gdf_sub[f"lr_signal_{ct_pair}"] = np.where(
                gdf_sub["cell_type_manual"] == ct_pair[0], gdf_sub[gene_a],
                np.where(gdf_sub["cell_type_manual"] == ct_pair[1], gdf_sub[gene_b],
                         np.nan)
            )
            gdf_sub[f"lr_signal_{gene_a}"] = np.where(
                gdf_sub["cell_type_manual"] == ct_pair[0], gdf_sub[gene_a], np.nan
            )
            gdf_sub[f"lr_signal_{gene_b}"] = np.where(
                gdf_sub["cell_type_manual"] == ct_pair[1], gdf_sub[gene_b], np.nan
            )

            ax = axes[ax_ind]
            plot_geneb = gdf_sub.plot(
                column=f"lr_signal_{gene_b}",
                cmap="Blues",
                legend=True,
                linewidth=0.1,
                alpha=0.7,
                edgecolor="black",
                ax=ax,
                legend_kwds={
                    "shrink": 0.5,
                    "aspect": 20
                }
            )
            cbar2 = plot_geneb.get_figure().axes[-1]
            cbar2.set_title(gene_b, fontsize=8)
            plot_gene_a = gdf_sub.plot(
                column=f"lr_signal_{gene_a}",
                cmap="Oranges",
                legend=True,
                linewidth=0.1,
                alpha=0.7,
                edgecolor="black",
                ax=ax,
                legend_kwds={
                    "shrink": 0.5,
                    "aspect": 20
                }
            )
            cbar1 = plot_gene_a.get_figure().axes[-1]
            cbar1.set_title(gene_a, fontsize=8)
            gdf_sub.geometry.plot(  # Add borders of NaN cells
                color="none",
                linewidth=0.1,
                alpha=0.7,
                edgecolor="black",
                ax=ax,
            )
            gdf_sub_nuc.boundary.plot(
                color="grey",
                linewidth=0.1,
                ax=ax,
                zorder=10
            )
            ax.set_title(f"{ct_pair[0]} → {ct_pair[1]}", fontsize=12)
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.axis('off')

        # Subplot for cell types
        ax_celltype = axes[-1]
        gdf_sub.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            alpha=0.6,
            edgecolor="black",
            ax=ax_celltype,
            legend_kwds={
                'bbox_to_anchor': (0.98, 0.5),
                'loc': 'center left',
                'fontsize': 9
            }
        )
        gdf_sub_nuc.boundary.plot(
            color="grey",
            linewidth=0.15,
            ax=ax_celltype,
            zorder=10
        )
        ax_celltype.set_title("Cell types", fontsize=12)
        ax_celltype.axis('off')
        ax_celltype.set_ylim(ax_celltype.get_ylim()[::-1])

        plt.suptitle(f"{gene_a} to {gene_b} on {core}", fontsize=16, y=0.95)
        plt.tight_layout()

        core_filename = core.replace(" ", "_")
        lr_pair_filename = lr_pair.replace("/", "-")
        os.makedirs(out_dir_wo_pt1_pt2 + f"{lr_pair_filename}/", exist_ok=True)
        plt.savefig(
            out_dir_wo_pt1_pt2 + f"{lr_pair_filename}/{core_filename}.png",
            dpi=300,
            bbox_inches="tight"
        )
plt.close()


# 4.2 One figure per LR pair and per core (with PT1 and PT2) ---------

out_dir_with_pt1_pt2 = os.path.join(fig_dir, "spatial_plot", "with_pt1_pt2")
os.makedirs(out_dir_with_pt1_pt2, exist_ok=True)

# Update ct_pairs to include PT1 and PT2
ct_pairs = [
    ("aPT", "Fib/Myofib"),
    ("aPT", "Immune"),
    ("Immune", "Fib/Myofib"),
    ("PT1", "Immune"),
    ("PT2", "Immune"),
]
unique_cts = {ct for pair in ct_pairs for ct in pair}

lr_pair = "APP_CD74"
core = cores[0]


for lr_pair in lr_pairs:
    gene_a, gene_b = lr_pair.split("_")

    # prepare expression table for this LR pair
    genes = [gene_a, gene_b]
    expr = adata[:, genes].to_df()
    expr.reset_index(inplace=True)  # sample_cell_id column present

    a_min, a_max = expr[gene_a].min(), expr[gene_a].max()
    b_min, b_max = expr[gene_b].min(), expr[gene_b].max()

    for core in cores:

        # subset geometries for this core
        gdf_sub = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core].copy()
        gdf_sub_nuc = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core].copy()

        # merge expression into gdf_sub
        gdf_sub = gdf_sub.merge(expr, on="sample_cell_id", how="left")

        # precompute boolean masks for relevant cell types
        masks = {ct: (gdf_sub["cell_type_manual"] == ct) for ct in unique_cts}

        # plotting grid
        n = len(ct_pairs)
        ncols = 2
        nrows = ceil(n / ncols)
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4 * ncols, 3.5 * nrows))
        axes = np.array(axes).reshape(-1)

        # keep track of which axes are used
        used_axes = []

        for i, ct_pair in enumerate(ct_pairs):
            ax = axes[i]
            used_axes.append(ax)
            src, tgt = ct_pair

            # create per-ct signals: ligand only in source cells, receptor only in target cells
            gdf_sub[f"lr_a_{i}"] = gdf_sub[gene_a].where(masks.get(src, False))
            gdf_sub[f"lr_b_{i}"] = gdf_sub[gene_b].where(masks.get(tgt, False))

            # plot receptor (blue) first, then ligand (orange) on top
            gdf_sub.plot(
                column=f"lr_b_{i}",
                cmap="Blues",
                vmin=b_min,
                vmax=b_max,
                legend=False,
                linewidth=0.05,
                alpha=0.8,
                edgecolor="black",
                ax=ax,
            )
            gdf_sub.plot(
                column=f"lr_a_{i}",
                cmap="Oranges",
                vmin=a_min,
                vmax=a_max,
                legend=False,
                linewidth=0.05,
                alpha=0.8,
                edgecolor="black",
                ax=ax,
            )

            # add nucleus boundaries
            gdf_sub_nuc.boundary.plot(
                color="grey",
                linewidth=0.1,
                ax=ax,
                zorder=10
            )

            # create colored title with gene names colored according to their colormaps
            norm_a = Normalize(vmin=a_min, vmax=a_max)
            norm_b = Normalize(vmin=b_min, vmax=b_max)
            cmap_a = plt.get_cmap("Oranges")
            cmap_b = plt.get_cmap("Blues")
            # pick color at max value (top of colormap)
            try:
                color_a = cmap_a(norm_a(a_max*0.7))
                color_b = cmap_b(norm_b(b_max*0.7))
            except Exception:
                color_a = cmap_a(1.0)
                color_b = cmap_b(1.0)
            # place two colored text pieces slightly left and right of center
            title_y = 1.02
            ax.text(0.46, title_y, f"{gene_a} in {src}", color=color_a,
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=12)
            ax.text(0.50, title_y, " to ", color="black",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=12)
            ax.text(0.54, title_y, f"{gene_b} in {tgt}", color=color_b,
                    transform=ax.transAxes, ha="left", va="bottom", fontsize=12)
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.axis("off")

        # turn off any unused axes
        for j in range(n, axes.size):
            axes[j].axis("off")
        # create two shared colorbars stacked vertically (gene_a above, gene_b below)
        plt.subplots_adjust(right=0.82, hspace=0.25)

        # ScalarMappables for colorbars
        sm_b = ScalarMappable(norm=Normalize(vmin=b_min, vmax=b_max), cmap="Blues")
        sm_b.set_array([])  # required for colorbar
        sm_a = ScalarMappable(norm=Normalize(vmin=a_min, vmax=a_max), cmap="Oranges")
        sm_a.set_array([])

        # layout for stacked colorbars to the right of the figure
        x = 0.84
        width = 0.02
        height = 0.18
        gap = 0.03

        # center vertically relative to plotting area used above (approximation)
        center_y = 0.15 + 0.7 / 2.0
        total_h = 2 * height + gap
        bottom_y = center_y - total_h / 2.0
        top_y = bottom_y + height + gap

        # gene_a (ligand) on top, gene_b (receptor) below
        cbar_ax_a = fig.add_axes([x, top_y, width, height])
        cbar_ax_b = fig.add_axes([x, bottom_y, width, height])

        cbar_a = fig.colorbar(sm_a, cax=cbar_ax_a)
        cbar_a.set_label(gene_a, fontsize=9)
        cbar_b = fig.colorbar(sm_b, cax=cbar_ax_b)
        cbar_b.set_label(gene_b, fontsize=9)

        # move labels from y-axis to axis title above each colorbar
        cbar_a.set_label("")  # remove y-axis label
        cbar_a.ax.set_title(gene_a, fontsize=9, pad=4)

        cbar_b.set_label("")  # remove y-axis label
        cbar_b.ax.set_title(gene_b, fontsize=9, pad=4)

        plt.suptitle(core, fontsize=14, y=0.98)
        plt.tight_layout(rect=[0, 0, 0.82, 1.0])

        core_filename = core.replace(" ", "_")
        lr_pair_filename = lr_pair.replace("/", "-")
        save_dir = os.path.join(out_dir_with_pt1_pt2, lr_pair_filename)
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(
            os.path.join(save_dir, f"{core_filename}.png"),
            dpi=600,
            bbox_inches="tight"
        )
        plt.close(fig)
