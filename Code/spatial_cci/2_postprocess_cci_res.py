#
# Created on Sun Sep 21 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script performs **postprocessing and contextual integration of local CCI results**
# from Xenium kidney biopsy data, linking interaction metrics to tissue-level factors.
#
# Datasets used/loaded:
# - Spatial cell/nucleus boundaries:
#       geodataframes/{cell_boundaries, nucleus_boundaries}.parquet
# - Processed Xenium transcriptomic data:
#       xenium/xenium_processed.h5ad
# - Precomputed Voronoi neighborhood graphs:
#       neighborhoods/voronoi_graphs_10.pkl
# - Precomputed local CellPhoneDB-derived interaction scores:
#       local_cci_results/all_cores/interaction_scores_*.parquet
#
# Datasets generated/saved:
# - Aggregated and annotated interaction table including mean 2-hop Factor2_ulm values:
#       local_cci_results/cci_score_df.parquet
#
# The analysis pipeline:
# 1. Load cell boundaries and Xenium transcriptomic data, harmonizing cell-type labels
#    (merging aPT/PEC into aPT).
# 2. Construct per-core 2-hop adjacency graphs from precomputed Voronoi connectivity.
# 3. Load precomputed ligand–receptor interaction scores across all roots.
# 4. Integrate tissue-level “Factor2_ulm” scores by mapping them to each cell and
#    averaging within 2-hop microenvironments.
# 5. Append the resulting mean Factor2_ulm to each interaction entry, enabling correlation
#    analyses between local signaling intensity and tissue phenotypic gradients.
#
# Output:
# A single consolidated `cci_score_df.parquet` storing all per-neighborhood
# ligand–receptor interaction data annotated with local mean Factor2_ulm values,
# suitable for downstream visualization or statistical modeling.
# --------------------------------------------------------------------------------------------


from statsmodels.nonparametric.smoothers_lowess import lowess
import re
import scanpy as sc
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial.distance import pdist, squareform
import matplotlib.patches as mpatches
import seaborn as sns
import os
import numpy as np
import geopandas as gpd
import networkx as nx
from matplotlib.colors import ListedColormap
import pickle
import warnings
warnings.filterwarnings(
    "ignore", message="nopython is set for njit and is ignored")
import sys

# Load helper functions ----------------------------
# sys.path.append(os.path.abspath("/data/zd775156/kidney_biopsy/code/cci_utils/"))
# from simple_cpdb import *
sys.path.append(os.path.abspath("/data/zd775156/kidney_biopsy/code/utils/"))
from running_on_slurm import running_on_slurm
from spatial_cci_helpers import *
import argparse

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

# Load CCI results ---------------------------------------------------
end = 198361

cci_score_df = load_continuous_chunks(
    "interaction_scores",
    folder=fig_dir + "local_cci_results/all_cores/",
    start=0, end=end
)

print("Number of roots:", cci_score_df["root"].nunique())

# Investigate nhs with no scores -----------------
print(cci_score_df["root"].nunique())
print(gdf_cell_boundaries["sample_cell_id"].nunique())

# Rename aPT/PEC to aPT in cci_score_df
cci_score_df["source"] = cci_score_df["source"].replace("aPT/PEC", "aPT")
cci_score_df["target"] = cci_score_df["target"].replace("aPT/PEC", "aPT")


# Filter cell types of interest --------------------------------

# convert source and target as strings in cci_score_df
cci_score_df["source"] = cci_score_df["source"].astype(str)
cci_score_df["target"] = cci_score_df["target"].astype(str)

cci_score_df["Factor2_ulm"] = cci_score_df["root"].map(factor2_map)
cci_score_df["lr_pair"] = cci_score_df["gene_a"] + "_" + cci_score_df["gene_b"]
cci_score_df["ct_pair"] = cci_score_df["source"] + "_to_" + cci_score_df["target"]

print(cci_score_df["ct_pair"].unique())


# Compute Mean Factor2_ulm of 2-hop neighbors -----------------

results = []
for main_root, neighbor_list in adj_2hop_dict.items():
    # include the main root itself
    roots_to_take = [main_root] + neighbor_list

    # subset data
    sub = cci_score_df[cci_score_df["root"].isin(roots_to_take)]

    # compute mean Factor2_ulm
    mean_val = sub["Factor2_ulm"].mean()

    results.append({
        "neighborhood": main_root,
        "roots": roots_to_take,
        "mean_Factor2_ulm": mean_val
    })

results = pd.DataFrame(results)
cci_score_df["mean_Factor2_ulm"] = cci_score_df["root"].map(
    results.set_index("neighborhood")["mean_Factor2_ulm"]
)
del results

# Save cci_score_df with mean_Factor2_ulm
cci_score_df.to_parquet(
    fig_dir + "local_cci_results/cci_score_df.parquet"
)
