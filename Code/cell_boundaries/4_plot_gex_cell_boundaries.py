#
# Created on Fri Aug 22 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script visualizes gene expression spatially across cell boundaries in Xenium kidney biopsy data.
#
# Datasets used/loaded:
# - Processed Xenium AnnData object: xenium_processed.h5ad
# - Cell boundary GeoDataFrames: cell_circles.parquet, nucleus_boundaries.parquet, cell_circles.parquet
# - Podocyte DEG metadata: POD_DEGs_controlled_for_Factor2_in_Xenium_panel.csv
#
# Data processing:
# - Merges normalized gene expression (from AnnData) with cell boundary spatial data.
# - Associates each cell with its gene expression and spatial boundary.
#
# Plots generated:
# - For each disease (contrast) and tissue core, for each podocyte DEG gene:
#     - Spatial plot of gene expression (gene expression value per cell, colored by expression,
# overlaid on cell boundaries).
#     - Plots are stratified by disease and tissue core.
#
# Output:
# - PNG images of spatial gene expression for each gene, disease, and core.
# - Saved to: {fig_dir}gene_exp/pod_deg_genes/{disease}/{gene}_{core}.png
# --------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import geopandas as gpd
import scanpy as sc
import pandas as pd

# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/cell_boundaries/"

adata = sc.read_h5ad(res_dir + "/xenium_processed.h5ad")

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_circles.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(res_dir+"nucleus_boundaries.parquet")
gdf_cell_circles = gpd.read_parquet(res_dir+"cell_circles.parquet")

samples = gdf_cell_boundaries["sample"].unique().tolist()

gene_ls = adata.var_names

# Merge gene exp with gdf_cell_boundaries --------------------------------

# Assert that adata.X is the normalized counts
assert any(adata.X.data != adata.X.data.astype(int))

# Define dataframe of gexp
gexp_df = pd.DataFrame(adata.X.todense())
gexp_df.columns = adata.var_names
gexp_df.index = adata.obs["TMA_ID"].astype(str) + "_" + adata.obs["cell_id"].astype(str)
# replace the index with .replace(" ", "_")
gexp_df.index = gexp_df.index.str.replace(" ", "_")

# Merge gexp with cell boundaries ----
gdf_cell_boundaries_gexp = gdf_cell_boundaries.merge(
    gexp_df,
    left_on="sample_cell_id",
    right_index=True
)


# Estimated space required to store all gexp figures  --------------
n_plots = gdf_cell_boundaries_gexp["TMA_ID"].nunique() * gexp_df.shape[1]
fig_size = 700000
est_space_required = n_plots * fig_size / 10**9
print(est_space_required)  # It requires 29 GB of space


# Plot DEGs from podocytes -------------------------------------------

# Load genes dataframe
pod_deg_df = pd.read_csv(
    "/data/Kramann_lab/kidney_biopsy/data/metadata/POD_DEGs_controlled_for_Factor2_in_Xenium_panel.csv",
    index_col=0
)

pod_deg_genes = pod_deg_df.set_index("contrast")["gene"]
genes_of_interest = pod_deg_genes.groupby(
    pod_deg_genes.index).apply(list).to_dict()


core = "Utrecht1 A6 T19-05689 K9 DN"

# Plot gene exp for all core and all genes of interest --------------
for disease in pod_deg_df["contrast"].unique():
    gdf_disease = gdf_cell_boundaries_gexp[gdf_cell_boundaries_gexp["disease_id"] == disease]

    for core in gdf_disease["TMA_ID"].unique():
        gdf_core = gdf_disease[gdf_disease["TMA_ID"] == core]

        for gene in genes_of_interest[disease]:

            fig, ax = plt.subplots(figsize=(15, 12))
            gdf_core.plot(
                column=gene,
                cmap="viridis",
                legend=True,
                linewidth=0.3,
                alpha=0.6,
                edgecolor="black",
                ax=ax,
                # vmin=min_factor2,
                # vmax=max_factor2,
                legend_kwds={
                    "shrink": 0.5,
                    "aspect": 20,
                    "pad": 0.01   # smaller pad moves colorbar closer
                }
            )
            gdf_core.boundary.plot(
                color="grey",
                linewidth=0.15,
                # alpha=0.4,
                ax=ax,
                zorder=10
            )
            ax.set_ylim(ax.get_ylim()[::-1])
            plt.title(f"{gene} -- {core}", fontsize=16, pad=2)
            ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labeltop=False, labelleft=False, labelright=False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            plt.tight_layout()
            plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
            plt.savefig(f"{fig_dir}gene_exp/pod_deg_genes/{disease}/{gene}_{core.replace(' ', '_')}.png",
                        dpi=300)
            plt.show()


