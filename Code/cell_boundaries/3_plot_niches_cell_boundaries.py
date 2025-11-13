#
# Created on Wed Sep 03 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script visualizes spatial niches and cell-type boundaries from Xenium kidney biopsy data.
#
# Datasets used/loaded:
# - Cell boundaries: geodataframes/cell_boundaries.parquet
# - Nucleus boundaries: geodataframes/nucleus_boundaries.parquet
# - Niche annotations (latent_leiden clusters): xenium_with_niches.h5ad
# - Additional niche annotations:
#     * xenium_with_niches_in_glom_niche.h5ad
#     * xenium_with_niches_in_glom_manual.obs.csv
#
# Plots generated:
#
# 1. Folder: niche_spatial_plot/whole_core/
#    - For each TMA core with available niche annotations:
#        Left panel: niches (latent_leiden_0.4) from Xenium data
#        Right panel: cell types (cell_type_manual)
#      Both panels overlay nucleus boundaries in grey.
#      Output: {core}.pdf
#
# 2. Folder: niche_spatial_plot/glom_niche/
#    - For each TMA core with glomerular niche annotations (latent_leiden_0.6):
#        Left panel: niches (latent_leiden_0.6)
#        Right panel: cell types (cell_type_manual)
#      Both panels overlay nucleus boundaries in grey.
#      Output: {core}_glom.pdf
#
# 3. Folder: niche_spatial_plot/glom_niche_manual/
#    - For each TMA core with manually defined glomerular niches:
#        Left panel: niches (manual latent_leiden_0.6)
#        Right panel: cell types (cell_type_manual)
#      Both panels overlay nucleus boundaries in grey.
#      Output: {core}_manual.pdf
#
# The script merges spatial boundary GeoDataFrames with niche annotations from AnnData objects,
# detects missing or unmatched cores, and visualizes all valid intersections.
# Custom palettes are used for cell types and niche clusters. All spatial plots are saved
# with inverted y-axes for correct spatial orientation.
# --------------------------------------------------------------------------------------------


import pandas as pd
import scanpy as sc
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import geopandas as gpd
import seaborn as sns


# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir_whole_core = '/data/Kramann_lab/kidney_biopsy/results/xenium/leon/niche_spatial_plot/whole_core/'


# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(res_dir+"geodataframes/nucleus_boundaries.parquet")

gdf_cell_boundaries.index.names = [None]


# Load niches adata ----------------------------------------------------------------

# Remove uns object to make adata readable
# with h5py.File("/data/Kramann_lab/kidney_biopsy/results/xenium/leon/xenium_with_niches.h5ad", "r+") as f:
#     del f["uns"]["spatial_neighbors"]["params"]

adata_niches = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/results/xenium/leon/xenium_with_niches.h5ad")


# adata_niches_cores = set(adata_niches.obs["TMA_ID"].unique().tolist())
# gdf_cores = set(gdf_cell_boundaries["TMA_ID"].unique().tolist())
# gdf_cores - adata_niches_cores
# adata_niches_cores - gdf_cores

# cores_suspect = list(gdf_cores - adata_niches_cores)
# cores_suspect = [cores_suspect[i] for i in [0, 1, 4]]

# for core in cores_suspect:
#     gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
#     gdf_core.plot(
#         column="cell_type_manual",
#         categorical=True,
#         legend=True,
#         # cmap=my_cmap,
#         linewidth=0.3,
#         alpha=0.6,
#         edgecolor="black",
#         legend_kwds={
#             'bbox_to_anchor': (1, 0.5),
#             'loc': 'center left',
#             'fontsize': 10
#         }
#     )
#     plt.title(f"Core: {core}")
#     plt.show()

# cores_suspect

# ['Utrecht 2 A5 T03-04110 K29 MN',
#  'Utrecht1 B5 T18-01374 K25 TIN',
#  'Utrecht2 C3 T17-11176 K24 TIN']

# Define variable from gdf ----------------------------

cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories
samples = gdf_cell_boundaries["sample"].unique().tolist()
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()
adata_cores = adata_niches.obs["TMA_ID"].unique().tolist()

set(adata_cores) - set(cores)
cores_intersect = list(set(adata_cores) & set(cores))

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")

# Define color palette for niches ------------------------
colors_niches = adata_niches.uns["latent_leiden_0.4_colors"]
cmap_niches = ListedColormap(colors_niches)

# Merge adata_niches with geodataframe -------------------------------------------
# Prepare datasets
adata_niches.obs["sample_cell_id"] = (adata_niches.obs["TMA_ID"].astype(str) +
                                      "_" + adata_niches.obs["cell_id"].astype(str))
# Find cells missing in the geodataframe
missing_in_gdf = list(set(gdf_cell_boundaries["sample_cell_id"]) - set(adata_niches.obs["sample_cell_id"]))
gdf_cell_boundaries[gdf_cell_boundaries["sample_cell_id"].isin(missing_in_gdf)]["TMA_ID"].value_counts()

# Merge them
gdf_cell_boundaries = gdf_cell_boundaries.merge(adata_niches.obs[["sample_cell_id", "latent_leiden_0.4"]],
                                                left_on="sample_cell_id",
                                                right_on="sample_cell_id",
                                                how="left")


# 1. Get core which have no niche data -----------------
for core in cores:
    gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
    if gdf_core['latent_leiden_0.4'].isna().all():
        print(core)

        # Plot the core
        gdf_core.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            alpha=0.6,
            edgecolor="black",
            legend_kwds={
                'bbox_to_anchor': (1, 0.5),
                'loc': 'center left',
                'fontsize': 10
            }
        )


cores_wo_niche = ["Utecht2 A1 T07-18857 K12 TIN",
                  "Utrecht2 C3 T17-11176 K24 TIN",
                  "Utrecht 2 A5 T03-04110 K29 MN",
                  "Utrecht1 B5 T18-01374 K25 TIN",
                  "Utrecht1 D4 T13-08471 K51 IgAN"]

adata_niche_wo_data = adata_niches[adata_niches.obs["TMA_ID"].isin(cores_wo_niche)]
adata_niche_wo_data
adata_niche_wo_data.obs["latent_leiden_0.4"]


# 1.2. Plot Niches -----------------------------------------------------------------
core = "Utrecht 2 B2 T19-13338 K27 TIN"

for core in cores_intersect:

    gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
    gdf_core_nuc = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core]

    gdf_core["latent_leiden_0.4"]

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    # Left: latent_leiden_0.4
    gdf_core.plot(
        column="latent_leiden_0.4",
        categorical=True,
        legend=True,
        cmap=cmap_niches,
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
    axes[0].set_title("Niches", fontsize=16)
    axes[0].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[0].spines.values():
        spine.set_visible(False)
    # Right: cell_type_manual
    gdf_core.plot(
        column="cell_type_manual",
        categorical=True,
        legend=True,
        cmap=my_cmap,
        linewidth=0.3,
        alpha=0.6,
        edgecolor="black",
        ax=axes[1],
        legend_kwds={
            'bbox_to_anchor': (1, 0.5),
            'loc': 'center left',
            'fontsize': 10
        }
    )
    gdf_core_nuc.boundary.plot(
        color="grey",
        linewidth=0.15,
        ax=axes[1],
        zorder=10
    )
    axes[1].set_ylim(axes[1].get_ylim()[::-1])
    axes[1].set_title("Cell type", fontsize=16)
    axes[1].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[1].spines.values():
        spine.set_visible(False)
    plt.subplots_adjust(wspace=0.001)
    plt.suptitle(f"{core}", fontsize=16)
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir_whole_core}{core.replace(' ', '_')}.pdf")
    plt.show()


# gdf_core.plot(
#     column="cell_type_manual",
#     categorical=True,
#     legend=True,
#     cmap=my_cmap,
#     linewidth=0.3,
#     alpha=0.6,
#     edgecolor="black",
#     legend_kwds={
#         'bbox_to_anchor': (1, 0.5),
#         'loc': 'center left',
#         'fontsize': 10
#     }
# )


# 2. Load different adatas with niches ---------------------------------------------------

# Remove uns object to make adata readable
# with h5py.File("/data/Kramann_lab/kidney_biopsy/results/" +
# "xenium/leon/xenium_with_niches_in_glom_niche.h5ad", "r+") as f:
#     del f["uns"]["spatial_neighbors"]["params"]

adata_niches_glom = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/" +
    "results/xenium/leon/xenium_with_niches_in_glom_niche.h5ad")

print(adata_niches_glom.obs["TMA_ID"].value_counts().sort_values().to_string())
print(gdf_cell_boundaries["TMA_ID"].value_counts().sort_values().to_string())

for core in cores:
    adata_niches_glom_core = adata_niches_glom[adata_niches_glom.obs["TMA_ID"] == core]
    if adata_niches_glom_core.obs['latent_leiden_0.6'].isna().all():
        print(core)


adata_niches_glom.obs["sample_cell_id"] = (adata_niches_glom.obs["TMA_ID"].astype(str) +
                                           "_" + adata_niches_glom.obs["cell_id"].astype(str))
gdf_cell_boundaries_niches_glom = gdf_cell_boundaries.merge(
    adata_niches_glom.obs[["sample_cell_id", "latent_leiden_0.6"]],
    left_on="sample_cell_id",
    right_on="sample_cell_id",
    how="right"
)

cores_suspect_glom = ["Utecht2 A1  T07-18857 K12 TIN",
                      "Utrecht1 D4  T13-08471 K51 IgAN"]

# 2.1. Get core which have no niche data -----------------
for core in cores_suspect_glom:
    gdf_core = gdf_cell_boundaries_niches_glom[gdf_cell_boundaries_niches_glom["TMA_ID"] == core]
    if gdf_core['latent_leiden_0.6'].isna().all():
        print(core)
        # Plot the core
        gdf_core.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            alpha=0.6,
            edgecolor="black",
            legend_kwds={
                'bbox_to_anchor': (1, 0.5),
                'loc': 'center left',
                'fontsize': 10
            }
        )

# 2.2. Plot Niches -----------------------------------------------------------------
cores_intersect_glom = list(
    set(adata_niches_glom.obs["TMA_ID"].unique()) &
    set(gdf_cell_boundaries_niches_glom["TMA_ID"].unique())
)
colors_niches_glom = adata_niches_glom.uns["latent_leiden_0.6_colors"]
cmap_niches_glom = ListedColormap(colors_niches_glom)

fig_dir_glom_niche = '/data/Kramann_lab/kidney_biopsy/results/xenium/leon/niche_spatial_plot/glom_niche/'

for core in adata_niches_glom.obs["TMA_ID"].unique():
    gdf_core = gdf_cell_boundaries_niches_glom[gdf_cell_boundaries_niches_glom["TMA_ID"] == core]
    gdf_core_nuc = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core]

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    # Left: latent_leiden_0.6
    gdf_core.plot(
        column="latent_leiden_0.6",
        categorical=True,
        legend=True,
        cmap=cmap_niches_glom,
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
    axes[0].set_title("Niches", fontsize=16)
    axes[0].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[0].spines.values():
        spine.set_visible(False)
    # Right: cell_type_manual
    gdf_core.plot(
        column="cell_type_manual",
        categorical=True,
        legend=True,
        cmap=my_cmap,
        linewidth=0.3,
        alpha=0.6,
        edgecolor="black",
        ax=axes[1],
        legend_kwds={
            'bbox_to_anchor': (1, 0.5),
            'loc': 'center left',
            'fontsize': 10
        }
    )
    gdf_core_nuc.boundary.plot(
        color="grey",
        linewidth=0.15,
        ax=axes[1],
        zorder=10
    )
    axes[1].set_ylim(axes[1].get_ylim()[::-1])
    axes[1].set_title("Cell type", fontsize=16)
    axes[1].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[1].spines.values():
        spine.set_visible(False)
    plt.subplots_adjust(wspace=0.001)
    plt.suptitle(f"{core}", fontsize=16)
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir_glom_niche}{core.replace(' ', '_')}_glom.pdf")
    # plt.show()
    plt.close()


# 3. Load different adatas with niches manual ---------------------------------------------------
adata_metadata = pd.read_csv(
    "/data/Kramann_lab/kidney_biopsy/results/xenium/leon/" +
    "exports_x_and_metadata/xenium_with_niches_in_glom_manual.obs.csv",
    index_col=0
)
adata_metadata["sample_cell_id"] = (adata_metadata["TMA_ID"].astype(str) +
                                    "_" + adata_metadata["cell_id"].astype(str))
adata_metadata.index = adata_metadata["sample_cell_id"]
adata_metadata.index.names = [None]
adata_metadata["latent_leiden_0.6"]


# Merge latent_leiden_0.6 with gdf_cell_boundaries ------------------------
gdf_cell_boundaries_niches_glom_manual = gdf_cell_boundaries.merge(
    adata_metadata[["sample_cell_id", "latent_leiden_0.6"]],
    left_on="sample_cell_id",
    right_on="sample_cell_id",
    how="right"
)

fig_dir_glom_niche_manual = '/data/Kramann_lab/kidney_biopsy/results/xenium/leon/niche_spatial_plot/glom_niche_manual/'
cores_intersect_manual = list(
    set(adata_metadata["TMA_ID"].unique()) & 
    set(gdf_cell_boundaries_niches_glom_manual["TMA_ID"].unique())
)
colors_niches_manual = sns.color_palette('husl', adata_metadata["latent_leiden_0.6"].nunique())
cmap_niches_manual = ListedColormap(colors_niches_manual)

for core in cores_intersect_manual:
    gdf_core = gdf_cell_boundaries_niches_glom_manual[gdf_cell_boundaries_niches_glom_manual["TMA_ID"] == core]
    gdf_core_nuc = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core]

    fig, axes = plt.subplots(1, 2, figsize=(28, 12))
    # Left: latent_leiden_0.6
    gdf_core.plot(
        column="latent_leiden_0.6",
        categorical=True,
        legend=True,
        cmap=cmap_niches_manual,
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
    axes[0].set_title("Niches (manual)", fontsize=16)
    axes[0].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[0].spines.values():
        spine.set_visible(False)
    # Right: cell_type_manual
    gdf_core.plot(
        column="cell_type_manual",
        categorical=True,
        legend=True,
        cmap=my_cmap,
        linewidth=0.3,
        alpha=0.6,
        edgecolor="black",
        ax=axes[1],
        legend_kwds={
            'bbox_to_anchor': (1, 0.5),
            'loc': 'center left',
            'fontsize': 10
        }
    )
    gdf_core_nuc.boundary.plot(
        color="grey",
        linewidth=0.15,
        ax=axes[1],
        zorder=10
    )
    axes[1].set_ylim(axes[1].get_ylim()[::-1])
    axes[1].set_title("Cell type", fontsize=16)
    axes[1].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                        labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    for spine in axes[1].spines.values():
        spine.set_visible(False)
    plt.subplots_adjust(wspace=0.001)
    plt.suptitle(f"{core}", fontsize=16)
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir_glom_niche_manual}{core.replace(' ', '_')}_glom_manual.pdf")
    plt.close()
    