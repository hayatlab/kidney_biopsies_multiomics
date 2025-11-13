#
# Created on Fri Aug 22 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script visualizes cell and nucleus boundaries from Xenium kidney biopsy spatial data.
#
# Datasets used/loaded:
# - Cell boundaries: geodataframes/cell_circles.parquet
# - Nucleus boundaries: nucleus_boundaries.parquet
# - Cell circles: cell_circles.parquet
#
# Plots generated:
#
# 1. Folder: cell_boundaries/
#    - For each sample, a whole-slide plot of cell boundaries colored by cell type (cell_type_manual).
#      Output: whole_tma_boundaries_{sample}.png
#
# 2. Folder: cell_boundaries/cell_type/
#    - For each TMA core (TMA_ID), a plot of cell boundaries colored by cell type (cell_type_manual),
#      with nucleus boundaries overlaid in grey.
#      Output: cell_type_{core}.png
#
# 3. Folder: cell_boundaries/factor2/
#    - For each TMA core, a plot of cell boundaries colored by Factor2_ulm values,
#      with nucleus boundaries overlaid in grey.
#      Output: Factor2_ulm_{core}.pdf
#
# 4. Folder: cell_boundaries/cell_type_factor2/
#    - For each TMA core, a two-panel plot:
#        Left: cell boundaries colored by cell type (cell_type_manual)
#        Right: cell boundaries colored by Factor2_ulm
#      Both panels overlay nucleus boundaries in grey.
#      Output: cell_type_factor2_{core}.pdf
#
# The script uses custom color palettes for cell types, inverts y-axes for correct spatial orientation,
# and saves all plots to the specified output directories.
# --------------------------------------------------------------------------------------------

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import geopandas as gpd
import seaborn as sns

# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/cell_boundaries/"

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(res_dir+"geodataframes/nucleus_boundaries.parquet")
gdf_cell_circles = gpd.read_parquet(res_dir+"geodataframes/cell_circles.parquet")

sample = "XETG00229__0042348__Region_3__20250416__104148"
samples = gdf_cell_boundaries["sample"].unique().tolist()
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()

cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")

# Get min and max of Factor2_ulm ---------------------------------------
min_factor2 = gdf_cell_boundaries["Factor2_ulm"].min()
max_factor2 = gdf_cell_boundaries["Factor2_ulm"].max()

# Analysis of Factor2 values ------------------------------------------------
metadata_df = gdf_cell_boundaries[["Factor2_ulm", "cell_type_manual"]].copy()
metadata_df = metadata_df[~metadata_df["Factor2_ulm"].isna()]

# Remove cell types with too few values
value_counts = metadata_df["cell_type_manual"].value_counts()
cts_to_rm = value_counts[value_counts < 10].index.tolist()
metadata_df = metadata_df[~metadata_df["cell_type_manual"].isin(cts_to_rm)]
metadata_df["cell_type_manual"] = metadata_df["cell_type_manual"].cat.remove_categories(cts_to_rm)

cell_types_density = metadata_df["cell_type_manual"].cat.categories
my_pal_dict = dict(zip(cell_types, my_pal))
palette_metadata = [my_pal_dict[ct] for ct in cell_types_density]

# Stacked density plot of Factor2_ulm by cell type
fig, ax = plt.subplots(figsize=(12, 8))
sns.kdeplot(
    data=metadata_df,
    x="Factor2_ulm",
    hue="cell_type_manual",
    multiple="stack",
    palette=palette_metadata,
    linewidth=0.1,
    color="grey",
    fill=True,
    alpha=0.7,
    ax=ax,
)
leg = ax.get_legend()
leg.set_title("Cell type")
plt.title("Distribution of Factor2 by cell type", fontsize=16)
plt.tight_layout()
plt.savefig(f"{fig_dir}Factor2_ulm_stacked_density_by_cell_type.png",
            bbox_inches='tight', dpi=150)
plt.show()


# Plot cell boundaries for whole slides -------------------------------------------
for sample in samples:
    fig, ax = plt.subplots(figsize=(12, 12))
    gdf_cell_boundaries[gdf_cell_boundaries["sample"] == sample].plot(
        column="cell_type_manual",
        categorical=True,
        cmap=my_cmap,
        legend=True,
        linewidth=0.1,
        edgecolor="none",
        ax=ax,
        legend_kwds={
            'bbox_to_anchor': (1.02, 0.5),
            'loc': 'center left',
            'fontsize': 10
        }
    )
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.title(f"Cell boundaries colored by cell type -- {sample}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}whole_tma_boundaries_{sample.replace(' ', '_')}.png",
                bbox_inches='tight')
    plt.show()


# 2. Plot cell boundaries per TMA_ID -------------------------------------
core = "Utrecht1 A6 T19-05689 K9 DN"

for core in cores:

    gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
    gdf_core_nuc = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core]

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    # Left: cell_type_manual
    gdf_core.plot(
        column="cell_type_manual",
        categorical=True,
        legend=True,
        cmap=my_cmap,
        linewidth=0.3,
        alpha=0.6,
        edgecolor="black",
        ax=axes[0],
        legend_kwds={
            'bbox_to_anchor': (1, 0.5),
            'loc': 'center left',
            'fontsize': 10
        }
    )
    gdf_core_nuc.boundary.plot(
        color="grey",
        linewidth=0.15,
        ax=axes[0],
        zorder=10
    )
    axes[0].set_ylim(axes[0].get_ylim()[::-1])
    axes[0].set_title("Cell type", fontsize=16)
    axes[0].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[0].spines.values():
        spine.set_visible(False)
    # Right: Factor2_ulm
    gdf_core.plot(
        column="Factor2_ulm",
        cmap="viridis",
        legend=True,
        vmax=max_factor2,
        vmin=min_factor2,
        linewidth=0.3,
        alpha=0.6,
        edgecolor="black",
        ax=axes[1],
        legend_kwds={
            "shrink": 0.5,
            "aspect": 20
        }
    )
    gdf_core_nuc.boundary.plot(
        color="grey",
        linewidth=0.15,
        ax=axes[1],
        zorder=10
    )
    axes[1].set_ylim(axes[1].get_ylim()[::-1])
    axes[1].set_title("Factor2", fontsize=16)
    axes[1].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[1].spines.values():
        spine.set_visible(False)
    plt.suptitle(f"{core}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}cell_type_factor2/same_f2_scale/cell_type_factor2_{core.replace(' ', '_')}.pdf")
    plt.show()

