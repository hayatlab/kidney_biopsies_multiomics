#
# Created on Wed Aug 27 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script analyzes spatial relationships and neighborhood compositions of cells in Xenium
# kidney biopsy samples using cell boundary geodataframes.
#
# Datasets used/loaded:
# - Cell boundaries, nucleus boundaries, and cell circles (GeoDataFrames, .parquet)
# - Precomputed centroid distance matrices (distances.pkl)
#
# Main analyses and outputs:
# 1. Plots of cell centroids and polygons colored by Factor2_ulm for a selected core.
#    - Output: Displayed plot (not saved).
#
# 2. Construction of distance-bounded Voronoi adjacency graphs for each tissue core.
#    - For each core and distance threshold (e.g., 10 μm), builds a graph where nodes are cells
#      and edges connect Voronoi neighbors within the threshold distance.
#    - Output: NetworkX graphs saved as pickles (voronoi_graphs_{threshold}.pkl).
#
# 3. Visualization of Voronoi adjacency graphs overlaid on cell boundaries.
#    - For each core and threshold, plots cell polygons colored by cell_type_manual, overlays
#      Voronoi graph edges (white lines between centroids), and nucleus boundaries.
#    - Output folder: {fig_dir}/voronoi_adj_graphs/{threshold}_mum/
#      - Files: voronoi_{core}.png (Voronoi graph for each core and threshold)
#
# 4. Computation and visualization of 2-hop cell neighborhood composition.
#    - For each cell in each core, extracts the 2-hop neighborhood from the Voronoi graph,
#      plots the subgraph (cell type and Factor2_ulm), and computes the cell type composition.
#    - Output folder: {fig_dir}/2hop_nh/{core}/
#      - Files: around_{cell_id}.png (2-hop neighborhood plots for first 10 cells per core)
#    - Output table: {res_dir}/nh_2hop_compo_{threshold}.csv
#      - Contains, for each cell, the counts of each cell type in its 2-hop neighborhood.
#
# The script enables spatial analysis of cell neighborhoods, visualizes spatial graphs, and
# quantifies local cell type environments in kidney tissue sections.
# --------------------------------------------------------------------------------------------

from scipy.spatial import Voronoi, voronoi_plot_2d
import networkx as nx
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
import scanpy as sc
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import pickle
import os


# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/neighborhoods/"

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(res_dir+"geodataframes/nucleus_boundaries.parquet")
gdf_cell_circles = gpd.read_parquet(res_dir+"geodataframes/cell_circles.parquet")


# Map sample_cell_id to nh_compo_df
# nh_compo_df["sample_cell_id"] = nh_compo_df.index.map(
#     gdf_cell_boundaries["sample_cell_id"])
# nh_compo_df.index = nh_compo_df["sample_cell_id"]
# gdf_cell_boundaries.index = gdf_cell_boundaries["sample_cell_id"]

# # Set sample_cell_id as index
# gdf_cell_boundaries = gdf_cell_boundaries.set_index("sample_cell_id")
# gdf_nucleus_boundaries = gdf_nucleus_boundaries.set_index("sample_cell_id")
# gdf_cell_circles = gdf_cell_circles.set_index("sample_cell_id")

# # Remove NaN values of factor2 ------------------------------------------------------
# print(gdf_cell_boundaries['Factor2_ulm'].isna().value_counts(), "\n")
# print(gdf_nucleus_boundaries['Factor2_ulm'].isna().value_counts(), "\n")
# print(gdf_cell_circles['Factor2_ulm'].isna().value_counts(), "\n")

# which_na_cell_boundaries = gdf_cell_boundaries[gdf_cell_boundaries['Factor2_ulm'].isna()].index.tolist()
# which_na_nucleus_boundaries = gdf_nucleus_boundaries[gdf_nucleus_boundaries['Factor2_ulm'].isna()].index.tolist()
# which_na_cell_circles = gdf_cell_circles[gdf_cell_circles['Factor2_ulm'].isna()].index.tolist()

# # Assert that cell_boundaries NAs contains all the other ones
# assert len(set(which_na_nucleus_boundaries) - set(which_na_cell_boundaries)) == 0
# assert len(set(which_na_cell_circles) - set(which_na_cell_boundaries)) == 0

# print(gdf_cell_boundaries['Factor2_ulm'].describe())
# # The mean of Factor2_ulm is not close to zero, so we can't replace NaNs with 0

# gdf_core = gdf_core[gdf_core['Factor2_ulm'].notna()]

# Additional filtering: remove cells which are more than 400 microns away from any other cell in the same core
# compute distances between gdf_cell_circles
# distances_circles = {}
# max_distances_circles = {}
# for core in cores:
# distances_circles[core] = (gdf_cell_circles[gdf_cell_circles["TMA_ID"] == core]['geometry'].
# apply(lambda x: gdf_cell_circles[gdf_cell_circles["TMA_ID"] == core]['geometry'].distance(x)))
#     # Compute the maximum

#     # Get index of any index which has a value of distances_circles[core].max() above 400
#     index_max_distance

# distances_circles
# distances_circles = gdf_cell_circles['geometry'].apply(lambda x: gdf_cell_circles['geometry'].distance(x))

# Compute centroids
gdf_cell_boundaries['centroid'] = gdf_cell_boundaries.geometry.centroid

cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories.tolist()


# Create mapping dict -----------------------------------------------
factor2_mapping = gdf_cell_boundaries["Factor2_ulm"]
egfr_mapping = gdf_cell_boundaries["egfr_at_time_of_biopsy"]
disease_mapping = gdf_cell_boundaries.drop_duplicates(
    subset=["TMA_ID"]).set_index("TMA_ID")["disease_id"].to_dict()

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")

# Seaborn palette
colors = [my_cmap(i / max(1, len(cell_types)-1)) for i in range(len(cell_types))]
my_palette = dict(zip(cell_types, colors))

cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()


# 1. Compute distances between centroids --------------------------

# Compute distances ---------------------
# distances = {}
# for core in cores:
#     gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
#     distances[core] = gdf_core.geometry.apply(
#         lambda x: gdf_core.geometry.distance(x))
# # Save distances as a pickle
# with open(res_dir + "distances.pkl", "wb") as f:
#     pickle.dump(distances, f)


print("--------------------------------------------------")
print("----- Finished computing distances ---------------")
print("--------------------------------------------------")

# Load distances
with open(res_dir + "distances.pkl", "rb") as f:
    distances = pickle.load(f)


# 2. Distance-based graph --------------------------------------------

# Define your distance threshold
threshold = 5  # Example: 10 micrometers

for core in cores:
    gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]

    # Initialize an empty list to store adjacency information
    adjacency_list = []

    # Iterate through each polygon
    for i, row in gdf_core.iterrows():
        # Get the distances from the current polygon to all others
        # distances = distances.reset_index(drop=True)
        dist_to_others = distances.loc[i]
        # dist_to_others = dist_to_others.sort_values().values.tolist()

        # Find indices of polygons within the threshold distance
        adjacent_indices = dist_to_others[dist_to_others < threshold].index.tolist()

        # Remove the current polygon from its own list
        adjacent_indices.remove(str(i))
        # convert to float
        adjacent_indices = [float(x) for x in adjacent_indices]

        # Append the list of adjacent polygons
        adjacency_list.append(adjacent_indices)

    # Add the adjacency information as a new column in the GeoDataFrame
    gdf_core['adjacency'] = adjacency_list

# Define Graph from adjacency list ------------------------------------------------
distance_graph = nx.Graph()
for idx, neighbors in zip(gdf_core.index, gdf_core['adjacency']):
    for j in neighbors:
        distance_graph.add_edge(idx, j)

positions = {i: (geom.centroid.x, geom.centroid.y)
             for i, geom in gdf_core.geometry.items()}

# Plot distance graph --------------------
fig, ax = plt.subplots(figsize=(20, 20))
gdf_core.plot(
    column="cell_type_manual",
    categorical=True,
    legend=True,
    cmap=my_cmap,
    linewidth=0.1,
    edgecolor="none",
    ax=ax,
    legend_kwds={
        'bbox_to_anchor': (1.02, 0.5),
        'loc': 'center left',
        'fontsize': 10
    }
)
for u, v in distance_graph.edges():
    x0, y0 = positions[u]
    x1, y1 = positions[v]
    ax.plot([x0, x1], [y0, y1], color="white", linewidth=0.5)
ax.set_aspect("auto")
plt.show()
plt.savefig(f"{fig_dir}adj_graphs/{core.replace(' ', '_')}.png",
            bbox_inches="tight")


# 3. Voronoi-based graph --------------------------------------------

max_threshold = 10
max_threshold = 15
max_threshold_ls = [10, 15, 20, 25]

voronoi_graphs_wo_bounding = {}
voronoi_diagrams = {}

for max_threshold in [10]:
    voronoi_graphs = {}
    gdf_cell_boundaries[f"voronoi_adj_list_{max_threshold}"] = np.nan

    for core in cores:

        gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
        # Compute Voronoi diagram
        coords = np.array([p.coords[0] for p in gdf_core.geometry.centroid])
        vor = Voronoi(coords)

        voronoi_diagrams[core] = vor

        # Build adjacency graph
        voronoi_graph = nx.Graph()
        for i, region_index in enumerate(vor.point_region):
            region = vor.regions[region_index]
            if -1 in region or len(region) == 0:  # skip infinite regions
                continue
            # For each vertex in region, connect neighboring points
            for ridge_points in vor.ridge_points:
                if i in ridge_points:
                    neighbor = ridge_points[0] if ridge_points[1] == i else ridge_points[1]
                    voronoi_graph.add_edge(i, neighbor)
        # Set node names in voronoi_graph to the index of distances
        voronoi_graph = nx.relabel_nodes(voronoi_graph, pd.Series(distances[core].index.values).to_dict())

        voronoi_graphs_wo_bounding[core] = voronoi_graph.copy()

        # Remove connexions between cells which are too far apart
        for u, v in voronoi_graph.edges():
            if distances[core][u].loc[v] > max_threshold:
                voronoi_graph.remove_edge(u, v)

        voronoi_graphs[core] = voronoi_graph.copy()

        # Add adjacency list of voronoi graph to gdf
        # gdf_cell_boundaries[f"voronoi_adj_list_{max_threshold}"].loc[gdf_core.index] = gdf_core.index.map(
        #     lambda x: list(voronoi_graphs[core].neighbors(x))
        # )

    # Save networkx graph separately as pickle
    with open(f"{res_dir}neighborhoods/voronoi_graphs_{max_threshold}.pkl", "wb") as f:
        pickle.dump(voronoi_graphs, f)

    # Save Voronoi diagrams without bounding as pickle
    with open(f"{res_dir}neighborhoods/voronoi_graphs_wo_bounding_{max_threshold}.pkl", "wb") as f:
        pickle.dump(voronoi_graphs_wo_bounding, f)

    # Save Voronoi diagrams as pickle
    with open(f"{res_dir}neighborhoods/voronoi_diagrams_{max_threshold}.pkl", "wb") as f:
        pickle.dump(voronoi_diagrams, f)


for max_threshold in [10]:

    # Load networkx graph as pickle
    with open(f"{res_dir}neighborhoods/voronoi_graphs_{max_threshold}.pkl", "rb") as f:
        voronoi_graphs = pickle.load(f)

    # Plot voronoi diagram for all cores ------------------------------------------

    for core in cores:

        gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
        # Compute Voronoi diagram
        coords = np.array([p.coords[0] for p in gdf_core.geometry.centroid])
        vor = Voronoi(coords)

        # Get xmin and xmax from gdf_core
        xmin, xmax = gdf_core.total_bounds[[0, 2]]
        ymin, ymax = gdf_core.total_bounds[[1, 3]]

        fig, ax = plt.subplots(figsize=(18, 18))
        voronoi_plot_2d(voronoi_diagrams[core], ax=ax, show_vertices=False, line_colors='black',
                        line_width=0.5, line_alpha=0.6, point_size=4)
        ax.axis('off')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_aspect('equal')
        ax.set_title(f"Voronoi diagram for core {core}")
        plt.savefig(f"{fig_dir}voronoi_diagrams/whole_core/voronoi_{core.replace(' ', '_')}.png",
                    bbox_inches="tight")
        # Add gdf_core boundaries
        gdf_core.boundary.plot(ax=ax, color="black", linewidth=0.5)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_aspect('equal')
        ax.set_title(f"Voronoi diagram for core {core}")
        plt.savefig(f"{fig_dir}voronoi_diagrams/whole_core/voronoi_with_cells_{core.replace(' ', '_')}.png",
                    bbox_inches="tight",
                    dpi=300)
        plt.close()

    # Plot voronoi graph for all cores ------------------------------------------

    for core in cores:

        gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
        gdf_nuc_core = gdf_nucleus_boundaries[gdf_nucleus_boundaries["TMA_ID"] == core]

        # Standard plot, white edges after bounding  ------------------------------------------
        fig, ax = plt.subplots(figsize=(20, 20))
        gdf_core.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            edgecolor="black",
            ax=ax,
            legend_kwds={
                'bbox_to_anchor': (1.02, 0.5),
                'loc': 'center left',
                'fontsize': 10
            }
        )
        gdf_nuc_core.boundary.plot(
            color="grey",
            linewidth=0.15,
            ax=ax,
            zorder=10
        )
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_title(f"{core} with Distance-Bounded Voronoi Graph (bounded at {max_threshold} μm)",
                     fontsize=16, y=0.98)
        ax.set_aspect("equal")
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for u, v in voronoi_graphs[core].edges():
            x0, y0 = gdf_core.loc[u].centroid.x, gdf_core.loc[u].centroid.y
            x1, y1 = gdf_core.loc[v].centroid.x, gdf_core.loc[v].centroid.y
            ax.plot([x0, x1], [y0, y1], color="white", linewidth=0.5)
        # Add scatterplot of centroids
        for i in gdf_core.index:
            x, y = gdf_core.loc[i].centroid.x, gdf_core.loc[i].centroid.y
            ax.scatter(x, y, color="white", s=3)
        plt.tight_layout()
        plt.subplots_adjust(left=0.01, bottom=0.01)
        plt.savefig(
            f"{fig_dir}voronoi_adj_graphs/{max_threshold}_mum/" +
            f"white_edges/voronoi_{core.replace(' ', '_')}.png",
            bbox_inches="tight",
            dpi=300
        )
        plt.close()

        # With grey edges, before bounding ------------------------------------------

        fig, ax = plt.subplots(figsize=(20, 20))
        gdf_core.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            edgecolor="black",
            ax=ax,
            legend_kwds={
                'bbox_to_anchor': (1.02, 0.5),
                'loc': 'center left',
                'fontsize': 10
            }
        )
        gdf_nuc_core.boundary.plot(
            color="grey",
            linewidth=0.15,
            ax=ax,
            zorder=10
        )
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_title(f"{core} with Distance-Bounded Voronoi Graph",
                     fontsize=16, y=0.98)
        ax.set_aspect("equal")
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for u, v in voronoi_graphs_wo_bounding[core].edges():
            x0, y0 = gdf_core.loc[u].centroid.x, gdf_core.loc[u].centroid.y
            x1, y1 = gdf_core.loc[v].centroid.x, gdf_core.loc[v].centroid.y
            ax.plot([x0, x1], [y0, y1], color="grey", linewidth=0.5)
        # Add scatterplot of centroids
        for i in gdf_core.index:
            x, y = gdf_core.loc[i].centroid.x, gdf_core.loc[i].centroid.y
            ax.scatter(x, y, color="white", s=3)
        plt.tight_layout()
        plt.subplots_adjust(left=0.01, bottom=0.01)
        os.makedirs(f"{fig_dir}voronoi_adj_graphs/{max_threshold}_mum/grey_edges_before_bounding/", exist_ok=True)
        plt.savefig(
            f"{fig_dir}voronoi_adj_graphs/{max_threshold}_mum/" +
            f"grey_edges_before_bounding/voronoi_{core.replace(' ', '_')}.png",
            bbox_inches="tight",
            dpi=300
        )
        plt.close()

        # With grey edges, after bounding ------------------------------------------

        fig, ax = plt.subplots(figsize=(20, 20))
        gdf_core.plot(
            column="cell_type_manual",
            categorical=True,
            legend=True,
            cmap=my_cmap,
            linewidth=0.3,
            edgecolor="black",
            ax=ax,
            legend_kwds={
                'bbox_to_anchor': (1.02, 0.5),
                'loc': 'center left',
                'fontsize': 10
            }
        )
        gdf_nuc_core.boundary.plot(
            color="grey",
            linewidth=0.15,
            ax=ax,
            zorder=10
        )
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_title(f"{core} with Distance-Bounded Voronoi Graph (bounded at {max_threshold} μm)",
                     fontsize=16, y=0.98)
        ax.set_aspect("equal")
        ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for u, v in voronoi_graphs[core].edges():
            x0, y0 = gdf_core.loc[u].centroid.x, gdf_core.loc[u].centroid.y
            x1, y1 = gdf_core.loc[v].centroid.x, gdf_core.loc[v].centroid.y
            ax.plot([x0, x1], [y0, y1], color="grey", linewidth=0.5)
        # Add scatterplot of centroids
        for i in gdf_core.index:
            x, y = gdf_core.loc[i].centroid.x, gdf_core.loc[i].centroid.y
            ax.scatter(x, y, color="white", s=3)
        plt.tight_layout()
        plt.subplots_adjust(left=0.01, bottom=0.01)
        os.makedirs(f"{fig_dir}voronoi_adj_graphs/{max_threshold}_mum/grey_edges_after_bounding/", exist_ok=True)
        plt.savefig(
            f"{fig_dir}voronoi_adj_graphs/{max_threshold}_mum/" +
            f"grey_edges_after_bounding/voronoi_{core.replace(' ', '_')}.png",
            bbox_inches="tight",
            dpi=300
        )
        plt.close()

    # # Create adj_list
    # adj_list = {
    #     core: {n: list(voronoi_graphs[core].neighbors(n)) for n in voronoi_graphs[core].nodes()}
    #     for core in cores
    # }
    # # Unlist the adj_list dict and make it one dict only
    # adj_list = {k: v for d in adj_list.values() for k, v in d.items()}


# Save gdf_cell_boundaries as parquet
gdf_cell_boundaries.to_parquet(
    f"{res_dir}geodataframes/cell_boundaries_graph.parquet")


print("--------------------------------------------------")
print("----- Finished computing Voronoi graphs ----------")
print("--------------------------------------------------")


# 4. Compute Neighborhood Celltype Composition ----------------------------------------

max_threshold = 10

for max_threshold in [10, 15, 20, 25]:

    with open(f"{res_dir}neighborhoods/voronoi_graphs_{max_threshold}.pkl", "rb") as f:
        voronoi_graphs = pickle.load(f)

    with open(f"{res_dir}neighborhoods/voronoi_diagrams.pkl", "rb") as f:
        voronoi_diagrams = pickle.load(f)

    # gdf_core["sample_cell_id"] = gdf_core["TMA_ID"].astype(str) + "_" + gdf_core["cell_id"].astype(str)
    # gdf_core.index = gdf_core["sample_cell_id"]

    # Remove index name
    gdf_core.index.names = [None]
    voronoi_graphs[core].nodes()

    nh_compo_ls = []

    for core in cores:
        gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
        voronoi_diagram_core = voronoi_diagrams[core]

        iter_counter = 0
        idx, row = gdf_core.iloc[0].name, gdf_core.iloc[0]

        for iter_counter, (idx, row) in enumerate(gdf_core.iterrows()):

            # Get 2-hop neighbors using networkx
            neighbors_2hop = set(nx.single_source_shortest_path_length(voronoi_graphs[core], idx, cutoff=2).keys())
            neighbors_2hop.discard(idx)

            # Subgraph of 2-hop neighborhood
            subgraph = voronoi_graphs[core].subgraph([idx] + list(neighbors_2hop))
            sub_gdf = gdf_core.loc[list(subgraph.nodes)]

            if iter_counter < 10:

                # 1. Plot the 10 first subgraphs ------------------------------

                os.makedirs(
                    f"{fig_dir}2hop_nh/{max_threshold}/cell_type_factor2/{core.replace(' ', '_')}",
                    exist_ok=True
                )

                fig, axes = plt.subplots(1, 2, figsize=(16, 8))
                # Left: Plot polygons colored by cell_type_manual
                sub_gdf.plot(
                    column="cell_type_manual",
                    categorical=True,
                    cmap=my_cmap,
                    legend=True,
                    linewidth=0.5,
                    edgecolor="black",
                    ax=axes[0],
                    aspect="equal",
                    legend_kwds={
                        'bbox_to_anchor': (1.02, 0.5),
                        'loc': 'center left',
                        'fontsize': 10
                    }
                )
                # Draw edges between centroids (left)
                sub_positions = {i: (geom.centroid.x, geom.centroid.y) for i, geom in sub_gdf.geometry.items()}
                for u, v in subgraph.edges():
                    x0, y0 = sub_positions[u]
                    x1, y1 = sub_positions[v]
                    axes[0].plot([x0, x1], [y0, y1], color="white", linewidth=1)
                # Draw centroids (left)
                for i, (x, y) in sub_positions.items():
                    axes[0].scatter(x, y, color="white", s=10)
                axes[0].scatter(sub_positions[idx][0],
                                sub_positions[idx][1], color="white", s=80)
                axes[0].set_title("Cell type")
                axes[0].axis('off')
                # Right: Plot polygons colored by Factor2_ulm
                norm = mcolors.Normalize(vmin=sub_gdf['Factor2_ulm'].min(), vmax=sub_gdf['Factor2_ulm'].max())
                sub_gdf.plot(
                    column="Factor2_ulm",
                    cmap="viridis",
                    legend=True,
                    linewidth=0.5,
                    edgecolor="black",
                    ax=axes[1],
                    aspect="equal",
                    norm=norm,
                    legend_kwds={
                        "shrink": 0.5,
                        "aspect": 20,
                        "pad": 0.01
                    }
                )
                # Draw edges between centroids (right)
                for u, v in subgraph.edges():
                    x0, y0 = sub_positions[u]
                    x1, y1 = sub_positions[v]
                    axes[1].plot([x0, x1], [y0, y1], color="white", linewidth=1)
                # Draw centroids (right)
                for i, (x, y) in sub_positions.items():
                    axes[1].scatter(x, y, color="white", s=10)
                axes[1].scatter(sub_positions[idx][0], sub_positions[idx][1],
                                color="white", s=80)
                axes[1].set_title("Factor2")
                axes[1].axis('off')

                plt.suptitle(f"2-hop Neighborhood around the cell {idx.split('_')[-1]}")
                plt.tight_layout()
                plt.savefig(
                    f"{fig_dir}2hop_nh/{max_threshold}/cell_type_factor2/" +
                    f"{core.replace(' ', '_')}/around_{idx.split('_')[-1]}.png")
                plt.close()
                # plt.show()

                # 2. Plot the 10 first subgraphs with diagrams ---------------------

                os.makedirs(
                    f"{fig_dir}voronoi_diagrams/subgraphs/{core.replace(' ', '_')}",
                    exist_ok=True
                )
                os.makedirs(
                    f"{fig_dir}2hop_nh/{max_threshold}/cell_type_only/" +
                    f"{core.replace(' ', '_')}",
                    exist_ok=True
                )

                xmin, xmax = sub_gdf.total_bounds[[0, 2]]
                ymin, ymax = sub_gdf.total_bounds[[1, 3]]

                fig, ax = plt.subplots(figsize=(8, 8))
                sub_gdf.plot(
                    column="cell_type_manual",
                    categorical=True,
                    cmap=my_cmap,
                    legend=True,
                    linewidth=0.5,
                    edgecolor="black",
                    ax=ax,
                    aspect="equal",
                    legend_kwds={
                        'bbox_to_anchor': (1.02, 0.5),
                        'loc': 'center left',
                        'fontsize': 10
                    }
                )
                # Draw edges between centroids (left)
                sub_positions = {i: (geom.centroid.x, geom.centroid.y)
                                 for i, geom in sub_gdf.geometry.items()}
                for u, v in subgraph.edges():
                    x0, y0 = sub_positions[u]
                    x1, y1 = sub_positions[v]
                    ax.plot([x0, x1], [y0, y1],
                            color="white", linewidth=1)
                # Draw centroids (left)
                for i, (x, y) in sub_positions.items():
                    ax.scatter(x, y, color="white", s=20)
                ax.scatter(sub_positions[idx][0],
                           sub_positions[idx][1], color="white", s=60)
                ax.set_title("Cell type")
                ax.axis('off')
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                for spine in ax.spines.values():
                    spine.set_visible(False)
                ax.set_aspect('equal')
                ax.set_title(
                    f"2-hop Neighborhood around the cell {idx.split('_')[-1]}"
                )

                plt.tight_layout()
                plt.savefig(
                    f"{fig_dir}2hop_nh/{max_threshold}/cell_type_only/" +
                    f"{core.replace(' ', '_')}/around_{idx.split('_')[-1]}.png")
                voronoi_plot_2d(voronoi_diagram_core, ax=ax, show_vertices=False, line_colors='black',
                                line_width=0.5, line_alpha=0.6, point_size=0)
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
                ax.set_title(
                    f"Voronoi diagram and 2-hop Neighborhood around the cell {idx.split('_')[-1]}"
                )
                plt.savefig(
                    f"{fig_dir}voronoi_diagrams/subgraphs/" +
                    f"{core.replace(' ', '_')}/around_{idx.split('_')[-1]}.png",
                    bbox_inches="tight"
                )
                plt.show()

            # Get cell type composition
            cell_type_counts = sub_gdf['cell_type_manual'].value_counts()
            cell_type_counts = cell_type_counts.reindex(cell_types, fill_value=0)

            # Save composition result to dataframe
            # row_values = [core, idx, row["cell_type_manual"], *cell_type_counts.values]
            # row_df = pd.DataFrame([row_values], columns=nh_compo_df.columns)
            # nh_compo_df = pd.concat([nh_compo_df, row_df], ignore_index=True)

            nh_compo_ls.append({
                "TMA_ID": core,
                "sample_cell_id": idx,
                "root_cell_type": row["cell_type_manual"],
                **dict(zip(cell_types, cell_type_counts.values))
            })

    nh_compo_df = pd.DataFrame(nh_compo_ls)

    # Save nh_compo_df
    nh_compo_df["disease"] = nh_compo_df["TMA_ID"].map(disease_mapping)
    nh_compo_df.to_csv(res_dir + f"nh_2hop_compo_{max_threshold}.csv")


# 5. Compute Neighborhood Niche Composition ----------------------------------------

# 5.1 Load niche info from adata object ------------------------------
adata_niches = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/results/xenium/leon/xenium_with_niches.h5ad",
    backed="r",
)
adata_niches.obs["sample_cell_id"] = (adata_niches.obs["TMA_ID"].astype(str) + "_" +
                                      adata_niches.obs["cell_id"].astype(str))
niche_mapping = adata_niches.obs.set_index(
    "sample_cell_id")["latent_leiden_0.4"].to_dict()
colors_niches = adata_niches.uns["latent_leiden_0.4_colors"]
adata_niches.file.close()

# 5.2 Add niche info to gdf_cell_boundaries ------------------------------
niches = sorted(list(set(niche_mapping.values())),
                key=lambda x: int(x.split("_")[-1]))
cmap_niches = ListedColormap(colors_niches)

# Map niche to gdf_cell_boundaries
gdf_cell_boundaries["niche"] = gdf_cell_boundaries.index.map(niche_mapping)

gdf_cell_boundaries["niche"].value_counts(dropna=False)

# Compute the TMA_IDs with have NaN niches
cores_wo_niche = gdf_cell_boundaries[gdf_cell_boundaries["niche"].isna()]["TMA_ID"].unique().tolist()
print(cores_wo_niche)
# ['Utecht2 A1 T07-18857 K12 TIN', 'Utrecht1 D4 T13-08471 K51 IgAN']

# There are only 2 cores without niches: remove them from gdf_cell_boundaries
gdf_cell_boundaries = gdf_cell_boundaries[~gdf_cell_boundaries["TMA_ID"].isin(cores_wo_niche)]
gdf_nucleus_boundaries = gdf_nucleus_boundaries[~gdf_nucleus_boundaries["TMA_ID"].isin(cores_wo_niche)]

cores_niche = gdf_nucleus_boundaries["TMA_ID"].unique().tolist()

with open(f"{res_dir}neighborhoods/voronoi_graphs_10.pkl", "rb") as f:
    voronoi_graphs = pickle.load(f)

nh_niche_compo_ls = []

for core in cores_niche:

    gdf_core = gdf_cell_boundaries[gdf_cell_boundaries["TMA_ID"] == core]
    os.makedirs(
        f"{fig_dir}2hop_nh_niche/{core.replace(' ', '_')}", exist_ok=True
    )

    for iter_counter, (idx, row) in enumerate(gdf_core.iterrows()):

        # Get 2-hop neighbors using networkx
        neighbors_2hop = set(nx.single_source_shortest_path_length(
            voronoi_graphs[core], idx, cutoff=2).keys())
        neighbors_2hop.discard(idx)

        # Subgraph of 2-hop neighborhood
        subgraph = voronoi_graphs[core].subgraph(
            [idx] + list(neighbors_2hop))
        sub_gdf = gdf_core.loc[list(subgraph.nodes)]

        # Plot the 10 first subgraphs
        if iter_counter < 10:

            fig, axes = plt.subplots(1, 2, figsize=(16, 8))
            # Left: Plot polygons colored by niche
            sub_gdf.plot(
                column="niche",
                categorical=True,
                cmap=cmap_niches,
                legend=True,
                linewidth=0.5,
                edgecolor="black",
                ax=axes[0],
                aspect="equal",
                legend_kwds={
                    'bbox_to_anchor': (1.02, 0.5),
                    'loc': 'center left',
                    'fontsize': 10
                }
            )
            # Draw edges between centroids (left)
            sub_positions = {i: (geom.centroid.x, geom.centroid.y)
                             for i, geom in sub_gdf.geometry.items()}
            for u, v in subgraph.edges():
                x0, y0 = sub_positions[u]
                x1, y1 = sub_positions[v]
                axes[0].plot([x0, x1], [y0, y1], color="white", linewidth=1)
            # Draw centroids (left)
            for i, (x, y) in sub_positions.items():
                axes[0].scatter(x, y, color="white", s=10)
            axes[0].scatter(sub_positions[idx][0],
                            sub_positions[idx][1], color="white", s=30)
            axes[0].set_title("Cell type")
            axes[0].axis('off')

            # Right: Plot polygons colored by Factor2_ulm
            norm = mcolors.Normalize(
                vmin=sub_gdf['Factor2_ulm'].min(), vmax=sub_gdf['Factor2_ulm'].max())
            sub_gdf.plot(
                column="Factor2_ulm",
                cmap="viridis",
                legend=True,
                linewidth=0.5,
                edgecolor="black",
                ax=axes[1],
                aspect="equal",
                norm=norm,
                legend_kwds={
                    "shrink": 0.5,
                    "aspect": 20,
                    "pad": 0.01
                }
            )
            # Draw edges between centroids (right)
            for u, v in subgraph.edges():
                x0, y0 = sub_positions[u]
                x1, y1 = sub_positions[v]
                axes[1].plot([x0, x1], [y0, y1], color="white", linewidth=1)
            # Draw centroids (right)
            for i, (x, y) in sub_positions.items():
                axes[1].scatter(x, y, color="white", s=10)
            axes[1].scatter(sub_positions[idx][0],
                            sub_positions[idx][1], color="white", s=30)
            axes[1].set_title("Factor2")
            axes[1].axis('off')

            plt.suptitle(
                f"2-hop Neighborhood around the cell {idx.split('_')[-1]}")
            plt.tight_layout()
            plt.savefig(
                f"{fig_dir}2hop_nh_niche/{core.replace(' ', '_')}/around_{idx.split('_')[-1]}.png")
            plt.close()
            # plt.show()

        # Get niche composition
        niche_counts = sub_gdf['niche'].value_counts()
        niche_counts = niche_counts.reindex(
            niches, fill_value=0)

        nh_niche_compo_ls.append({
            "TMA_ID": core,
            "sample_cell_id": idx,
            "root_niche": row["niche"],
            **dict(zip(niches, niche_counts.values))
        })

nh_niche_compo_df = pd.DataFrame(nh_niche_compo_ls)

# Save nh_niche_compo_df
nh_niche_compo_df["disease"] = nh_niche_compo_df["TMA_ID"].map(disease_mapping)
nh_niche_compo_df.to_csv(res_dir + "neighborhoods/nh_2hop_niche_compo.csv")


print("--------------------------------------------------------")
print("----- Finished computing 2-hop nh niche composition ----")
print("--------------------------------------------------------")


