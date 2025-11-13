#
# Created on Sat Aug 30 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script analyzes and visualizes cell type composition in 2-hop neighborhoods from Xenium
# kidney biopsy spatial transcriptomics data.
#
# Datasets used/loaded:
# - Cell boundaries and metadata: geodataframes/cell_boundaries.parquet
# - 2-hop neighborhood composition: nh_2hop_compo_new.csv
#
# Datasets saved:
# - Neighborhood composition with clinical variables: nh_compo_df_w_vars.csv
#
# Plots generated (by output folder):
#
# heatmap_compo.png, heatmap_per_core/, heatmap_per_disease/:
#   - Heatmaps of cell type proportions in 2-hop neighborhoods:
#       - Overall (all cells)
#       - Stratified by tissue core (TMA_ID)
#       - Stratified by disease group
#
# self_excluding/:
#   - Heatmap of cell type composition in 2-hop neighborhoods, excluding the root cell type
#   - Smoothed ribbon plot of self-excluded cell type proportions vs Factor2_ulm
#
# self_excluding/by_root_cell_type/:
#   - Smoothed ribbon plots of self-excluded cell type proportions vs Factor2_ulm,
#     stratified by root cell type
#
# vs_factor2/:
#   - Barplot: stacked cell type proportions vs Factor2_ulm
#   - Segment plot: cell type proportions as horizontal segments vs Factor2_ulm
#   - Violin plot: distribution of Factor2_ulm by root cell type
#   - Smoothed spline plots: cell type proportions vs Factor2_ulm (scatter + spline, faceted by cell type)
#   - Smoothed ribbon plot: stacked cell type proportions vs Factor2_ulm
#
# vs_factor2/by_disease/:
#   - Smoothed spline and ribbon plots of cell type proportions vs Factor2_ulm,
#     stratified by disease group
#
# vs_factor2/by_root_cell_type/:
#   - Smoothed ribbon plots of cell type proportions vs Factor2_ulm,
#     stratified by root cell type
#
# vs_<clinical_variable>/ (for each clinical variable of interest):
#   - Smoothed spline and ribbon plots of cell type proportions vs the variable
#
# The script computes and visualizes cell type composition in 2-hop neighborhoods globally,
# by core, by disease, and as a function of clinical/latent variables (e.g., Factor2_ulm),
# including self-excluding analyses and smoothed trends using splines.
# --------------------------------------------------------------------------------------------

# import joypy
# from cpsplines import ConstrainedPRegression
import scanpy as sc
from matplotlib.colors import ListedColormap
from scipy.interpolate import LSQUnivariateSpline
import numpy as np
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import os


max_threshold = 10

# Set directories ---------------------------
data_dir = "/data/Kramann_lab/kidney_biopsy/data/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir = f"/data/Kramann_lab/From_Vivien/kidney_biopsy/neighborhoods/2hop_nh_compo/{max_threshold}_mum/"

# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_boundaries.parquet")

# Define parameters ------
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()
cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories.tolist()
diseases = gdf_cell_boundaries["disease_id"].cat.categories.tolist()

# Create custom palette for cell_types ---------------------------------
colors_tab20 = plt.get_cmap("tab20")(range(20))
colors_husl = sns.color_palette('husl', len(cell_types) - 20)
my_pal = list(colors_tab20) + colors_husl
my_cmap = ListedColormap(my_pal, name="my_custom")

# For seaborn
colors = [my_cmap(i / max(1, len(cell_types)-1))
          for i in range(len(cell_types))]
my_palette = dict(zip(cell_types, colors))

# Load nh_compo_df -------
nh_compo_df = pd.read_csv(
    res_dir + f"neighborhoods/nh_2hop_compo_{max_threshold}.csv",
    index_col=0
)
nh_compo_df.index = nh_compo_df["sample_cell_id"]

# To be removed later, after having cleaned nh_compo_df
disease_mapping = gdf_cell_boundaries.drop_duplicates(
    subset=["TMA_ID"]).set_index("TMA_ID")["disease_id"].to_dict()
nh_compo_df["disease"] = nh_compo_df["TMA_ID"].map(disease_mapping)
nh_compo_df["sample_cell_id"] = nh_compo_df.index
if "sample_cell_id.1" in nh_compo_df.columns:
    del nh_compo_df["sample_cell_id.1"]

# Warning: the index of nh_compo_df is not unique, which may cause issues later on.
if len(nh_compo_df.index.unique()) != len(nh_compo_df.index):
    print("Warning: the index of nh_compo_df is not unique, which may cause issues later on.")
# Simple fix: only take the first occurrence of each index
nh_compo_df = nh_compo_df[~nh_compo_df.index.duplicated(keep='first')]


# 1. Compute 2-hop nh proportions ---------------------------------------

nh_compo_prop_df = nh_compo_df.copy()
nh_compo_prop_df[cell_types] = nh_compo_prop_df[cell_types].div(
    nh_compo_prop_df[cell_types].sum(axis=1), axis=0
)
# Assert the sum by row is close to one
assert all(np.isclose(nh_compo_prop_df[cell_types].sum(axis=1), 1))

# Group the rows by cell_id and average the proportions
nh_compo_prop_df = nh_compo_prop_df.drop(
    columns=["sample_cell_id", "TMA_ID", "disease"]).groupby('root_cell_type').mean().fillna(0)
# arrange the rows by the order in the list "cell_types"
nh_compo_prop_df = nh_compo_prop_df.reindex(cell_types, axis=0)
# Set NaN to variables which are not present in the 2-hop neighborhood of any cell
nan_rows = nh_compo_prop_df.index[~nh_compo_prop_df.index.isin(nh_compo_df['root_cell_type'])].tolist()
for col in nan_rows:
    nh_compo_prop_df[col] = np.nan

# Plot a heatmap of the values
plt.figure(figsize=(10, 8))
sns.heatmap(nh_compo_prop_df[cell_types], cmap="viridis", cbar_kws={'label': 'Proportion'})
plt.title("Cell Type Composition (2-hop Neighborhood)")
plt.xlabel("Cell Types in the 2-hop Neighborhood")
plt.ylabel("Cell")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{fig_dir}heatmap_compo.png")
plt.show()


# 2. Compute 2-hop nh proportions per core ---------------------------------------

# Get the max value of prop across the cores:
vmax = nh_compo_prop_df[cell_types].max().max()

for core in cores:
    nh_compo_core = nh_compo_df[nh_compo_df["TMA_ID"] == core].copy()
    nh_compo_prop_core = nh_compo_core.copy()
    nh_compo_prop_core[cell_types] = nh_compo_prop_core[cell_types].div(
        nh_compo_prop_core[cell_types].sum(axis=1), axis=0
    )

    # Assert the sum by row is close to one (ignore empty rows)
    assert all(np.isclose(nh_compo_prop_core[cell_types].
                          sum(axis=1)[nh_compo_prop_core[cell_types].sum(axis=1) > 0], 1))
    # Group by cell type and average
    nh_compo_prop_core = nh_compo_prop_core.drop(
        columns=["sample_cell_id", "TMA_ID", "sample_cell_id", "disease"], errors="ignore"
    ).groupby('root_cell_type').mean().fillna(0)
    nh_compo_prop_core = nh_compo_prop_core.reindex(cell_types, axis=0)

    # Set NaN to variables not present in the 2-hop neighborhood of any cell
    nan_rows = nh_compo_prop_core.index[~nh_compo_prop_core.index.isin(nh_compo_core['root_cell_type'])].tolist()
    for col in nan_rows:
        nh_compo_prop_core[col] = np.nan

    # Plot heatmap for this core
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        nh_compo_prop_core[cell_types], cmap="viridis",
        vmax=vmax,
        cbar_kws={'label': 'Proportion'})
    plt.title(f"Cell Type Composition (2-hop Neighborhood)\n{core}")
    plt.xlabel("Cell Types in the 2-hop Neighborhood")
    plt.ylabel("Cell")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{fig_dir}heatmap_per_core/{core.replace(' ', '_')}.png")
    plt.show()


# 3. Compute 2-hop nh proportions per disease ---------------------------------------

for disease in diseases:

    nh_compo_disease = nh_compo_df[nh_compo_df["disease"] == disease].copy()
    nh_compo_prop_disease = nh_compo_disease.copy()
    nh_compo_prop_disease[cell_types] = nh_compo_prop_disease[cell_types].div(
        nh_compo_prop_disease[cell_types].sum(axis=1), axis=0
    )

    # Assert the sum by row is close to one (ignore empty rows)
    assert all(np.isclose(nh_compo_prop_disease[cell_types].sum(axis=1)
                          [nh_compo_prop_disease[cell_types].sum(axis=1) > 0], 1))
    # Group by cell type and average
    nh_compo_prop_disease = nh_compo_prop_disease.drop(
        columns=["sample_cell_id", "TMA_ID", "disease"], errors="ignore"
    ).groupby('root_cell_type').mean().fillna(0)
    nh_compo_prop_disease = nh_compo_prop_disease.reindex(cell_types, axis=0)

    # Set NaN to variables not present in the 2-hop neighborhood of any cell
    nan_rows = (nh_compo_prop_disease.index[~nh_compo_prop_disease.
                                            index.isin(nh_compo_disease['root_cell_type'])].tolist())
    for col in nan_rows:
        nh_compo_prop_disease[col] = np.nan

    # Plot heatmap for this disease
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        nh_compo_prop_disease[cell_types], cmap="viridis",
        vmax=vmax,
        cbar_kws={'label': 'Proportion'})
    plt.title(f"Cell Type Composition (2-hop Neighborhood)\nDisease: {disease}")
    plt.xlabel("Cell Types in the 2-hop Neighborhood")
    plt.ylabel("Cell")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{fig_dir}heatmap_per_disease/{str(disease).replace(' ', '_')}.png")
    plt.show()
    # plt.show()


# 4. Heatmap of Self-Excluding Cell Type Composition --------------------------

nh_compo_prop_se_df = nh_compo_df.copy()
# For each row, set the count of the cell which is in cell_id to zero
for idx, row in nh_compo_prop_se_df.iterrows():
    nh_compo_prop_se_df.at[idx, row['root_cell_type']] = 0

# Compute proportions for self-excluding neighborhood composition
denom = nh_compo_prop_se_df[cell_types].sum(axis=1)
non_empty_rows = denom != 0
nh_compo_prop_se_df = nh_compo_prop_se_df.loc[non_empty_rows]
nh_compo_prop_se_df[cell_types] = nh_compo_prop_se_df[cell_types].div(
    denom.replace(0, np.nan), axis=0)

assert all(np.isclose(nh_compo_prop_se_df[cell_types].sum(axis=1), 1))

# Group by cell type and average
nh_compo_prop_se_df = nh_compo_prop_se_df.drop(
    columns=["sample_cell_id", "TMA_ID", "disease"], errors="ignore"
).groupby('root_cell_type').mean()
nh_compo_prop_se_df = nh_compo_prop_se_df.reindex(cell_types, axis=0)
nh_compo_prop_se_df.sum(axis=1)
# Set NaN to variables which are not present in the 2-hop neighborhood of any cell
nan_rows_ct = nh_compo_prop_se_df.index[nh_compo_prop_se_df.isna().any(axis=1)].tolist()
for col in nan_rows_ct:
    nh_compo_prop_se_df[col] = np.nan
# Convert data type to float
nh_compo_prop_se_df = nh_compo_prop_se_df.astype(float)
# Set the diagonal values to NaN
for ind_row, ind_col in zip(nh_compo_prop_se_df.index, nh_compo_prop_se_df.columns):
    nh_compo_prop_se_df.loc[ind_row, ind_col] = np.nan

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(nh_compo_prop_se_df, cmap="viridis", cbar_kws={'label': 'Proportion'})
plt.title("Cell Type Composition (2-hop Neighborhood, Self-Excluded)")
plt.xlabel("Cell Types in the 2-hop Neighborhood")
plt.ylabel("Cell")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{fig_dir}self_excluding/ct_compo_self_excluded.png")
plt.show()


# 5. Compute nhood composition vs quantitative variables -------

# Define nh_prop_df -----------------
nh_prop_df = nh_compo_df.copy()
nh_prop_df[cell_types] = nh_prop_df[cell_types].div(nh_prop_df[cell_types].sum(axis=1), axis=0)

# Assert no NaN indices
assert nh_prop_df.index.isna().sum() == 0

# Define sample_cell_id
# nh_prop_df["cell_id"] = [f.split('_')[-1] for f in nh_prop_df.index.tolist()]
# nh_prop_df["sample_cell_id"] = nh_prop_df['TMA_ID'].astype(str) + '_' + nh_prop_df["cell_id"].astype(str)
# nh_prop_df.index = nh_prop_df["sample_cell_id"]
# nh_prop_df.index.names = [None]

# Load adata with clinical variables ------------------------------
adata = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/data/scVI_kidney_biopsy_with_clinical.h5ad",
    backed="r")
adata_obs = adata.obs.drop_duplicates(subset=["ID_number", "classification", "jm_slope"])
del adata

adata_mapping = adata_obs[
    ["ID_number", "jm_slope", "classification", "proteinuria_group", "diabetes",
     "hypertension_clinical", "Proteinuria over time", "egfr_group"]].set_index("ID_number")
# Remove rows where index is NaN
adata_mapping = adata_mapping[~adata_mapping.index.isna()]

# Add clinical variables from single cell data -------------------------
nh_prop_df.loc[:, "ID_number"] = [f.split(" ")[-2] for f in nh_prop_df["TMA_ID"]]
nh_prop_df = nh_prop_df.merge(
    adata_mapping, left_on="ID_number", right_on="ID_number", how="left"
).set_index(nh_prop_df.index)
del nh_prop_df["ID_number"]

# Workaround to the bug in adata mapping --------------------------------
gdf_cell_boundaries.index = gdf_cell_boundaries["sample_cell_id"]
nh_prop_df.index = nh_prop_df["sample_cell_id"]
# Remove doublet index by keeping the first occurrence
nh_prop_df = nh_prop_df[~nh_prop_df.index.duplicated(keep='first')]
gdf_cell_boundaries = gdf_cell_boundaries[~gdf_cell_boundaries.index.duplicated(keep='first')]


# Add clinical variables from gdf_boundaries ----------------------------------
vars_from_gdf = ["egfr_at_time_of_biopsy",
                 "proteinuria_around_biopsy protein to creatine ratio",
                 *["Factor"+str(i)+"_ulm" for i in range(1, 6)]]
for var in vars_from_gdf:
    nh_prop_df[var] = nh_prop_df.index.map(gdf_cell_boundaries[var])

vars_of_interest = vars_from_gdf + ["jm_slope"]

# Save nh_prop_df with variables ------------- ------
nh_prop_df.to_csv(f"{fig_dir}nh_compo_df_{max_threshold}_w_vars.csv")

# 5.1 Compute composition as a function of Factor2 -------------
nh_compo_prop_f2_df = nh_prop_df.copy()
nh_compo_prop_f2_df = nh_compo_prop_f2_df[~nh_compo_prop_f2_df["Factor2_ulm"].isna()]
nh_compo_prop_f2_df = nh_compo_prop_f2_df.sort_values(by="Factor2_ulm")


# 5.2 Plot barplot w.r.t. Factor2 ----------------

segments = nh_compo_prop_f2_df[cell_types].to_dict()
segments = {
    cell_type: [f for f in inner_dict.values()]
    for cell_type, inner_dict in segments.items()
}
y_positions = nh_compo_prop_f2_df["Factor2_ulm"].values
bottom = np.zeros(len(nh_compo_prop_f2_df.index.tolist()))
height = 0.01

fig, ax = plt.subplots(figsize=(6, 4))
for name, vals in segments.items():
    ax.barh(y_positions, vals, left=bottom, height=height,
            alpha=0.2, label=name, color=my_palette[name])
    bottom += np.array(vals)
ax.set_yticks(y_positions)
ax.set_yticklabels("")
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/ct_prop_vs_factor2_bars.png")
plt.show()


# 5.3 Plot 2D points w.r.t Factor2 -------------------------

# Compute the cumulative sum of each row and replace it
nh_compo_prop_f2_df_end = nh_compo_prop_f2_df.copy()
nh_compo_prop_f2_df_end["sample_cell_id"] = nh_compo_prop_f2_df_end.index
nh_compo_prop_f2_df_end.loc[:, cell_types] = nh_compo_prop_f2_df[cell_types].cumsum(axis=1)

nh_compo_prop_f2_df_start = nh_compo_prop_f2_df.copy()
nh_compo_prop_f2_df_start["sample_cell_id"] = nh_compo_prop_f2_df_start.index
nh_compo_prop_f2_df_start[cell_types] = nh_compo_prop_f2_df[cell_types].shift(
    1, axis=1).cumsum(axis=1)
nh_compo_prop_f2_df_start[cell_types] = nh_compo_prop_f2_df_start[cell_types].fillna(0)

# Pivot wider
nh_compo_prop_f2_df_end_longer = nh_compo_prop_f2_df_end.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types, var_name="cell_type", value_name="proportion_end"
)

# Assert that all values are close to 1
assert np.allclose(nh_compo_prop_f2_df_start[cell_types[0]], 0)
assert np.allclose(nh_compo_prop_f2_df_end[cell_types[len(cell_types)-1]], 1)


nh_compo_prop_f2_df_start_longer = nh_compo_prop_f2_df_start.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types,
    var_name="cell_type", value_name="proportion_start"
)
# merge by cell_type and Factor2_ulm
nh_compo_prop_f2_df_longer_seg = pd.merge(
    nh_compo_prop_f2_df_start_longer,
    nh_compo_prop_f2_df_end_longer,
    on=["Factor2_ulm", "sample_cell_id", "cell_type"]
)

# Line plot
plt.figure(figsize=(12, 10))
for cell_type in cell_types:
    df_ct = nh_compo_prop_f2_df_longer_seg[nh_compo_prop_f2_df_longer_seg["cell_type"] == cell_type]
    plt.hlines(
        y=df_ct["Factor2_ulm"],
        xmin=df_ct["proportion_start"],
        xmax=df_ct["proportion_end"],
        color=my_palette[cell_type],
        label=cell_type,
        linewidth=0.2,
        alpha=0.4,
    )
# Only show legend once per cell type
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Type")
plt.xlabel("Proportion in 2-hop Neighborhood")
plt.ylabel("Factor2_ulm")
plt.title("Cell Type Proportion vs Factor2_ulm (Stacked Segments)")
plt.savefig(
    f"{fig_dir}vs_factor2/ct_prop_vs_factor2_segments.png")
plt.tight_layout()
plt.show()


# 5. Plot the distribution of Factor2 per cell type (no nh_compo) -----------

plt.figure(figsize=(12, 8))
plt.grid(True, axis='y', linestyle='--', alpha=0.5)
sns.violinplot(
    data=nh_compo_prop_f2_df,
    x="root_cell_type",
    y="Factor2_ulm",
    palette=my_palette,
    order=cell_types,
    cut=0
)
plt.xticks(rotation=45, ha='right')
plt.xlabel("Root Cell Type")
plt.ylabel("Factor2_ulm")
plt.title("Distribution of Factor2_ulm by Root Cell Type")
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/factor2_violin_by_root_cell_type.png")
plt.show()


# 5.4 Smooth spline compo: Factor2 --------------------------------------------

nh_compo_prop_f2_df_longer = nh_compo_prop_f2_df.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types, var_name="cell_type", value_name="proportion"
)

# Define x_seq for plotting
factor2_min = nh_compo_prop_f2_df_longer["Factor2_ulm"].min()
factor2_max = nh_compo_prop_f2_df_longer["Factor2_ulm"].max()

# Different choices for restricted Factor2 space --------------------
# x_seq = np.linspace(factor2_min + (factor2_max - factor2_min) * 0.05,
#                     factor2_max - (factor2_max - factor2_min) * 0.05, 100)
# x_seq = np.linspace(factor2_min, factor2_max, 100)
# x_seq = np.linspace(nh_compo_prop_f2_df_longer["Factor2_ulm"].quantile(0.01),
#                     nh_compo_prop_f2_df_longer["Factor2_ulm"].quantile(0.99), 100)
x_seq = np.linspace(-4, 8, 100)  # Chosen method: manually set the x maximal values

fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
for ax, cell_type in zip(axes.flatten(), cell_types):
    aaa = (nh_compo_prop_f2_df_longer[nh_compo_prop_f2_df_longer["cell_type"] == cell_type].
           sort_values(by="Factor2_ulm"))

    x = aaa["Factor2_ulm"].values
    y = aaa["proportion"].values

    knots = np.linspace(x_seq[0], x_seq[-1], num=10)[1:-1]
    spline = LSQUnivariateSpline(x, y, t=knots, k=3)

    prop_end_fit = np.clip(spline(x_seq), 0, 1)

    # Save fitted spline values
    fitted_spline_df[cell_type] = prop_end_fit

    sns.scatterplot(
        data=aaa,
        y="proportion",
        x="Factor2_ulm",
        hue="cell_type",
        palette=my_palette,
        alpha=0.7,
        s=1.5,
        edgecolor=None,
        ax=ax
    )
    ax.plot(x_seq, prop_end_fit, label='Fitted Spline', color="black", linewidth=2)
    ax.set_title(cell_type)
    ax.set_xlabel("Factor2")
    if ax not in axes[-1]:
        ax.set_xlabel("")
        ax.xaxis.set_visible(False)
    if ax not in axes[:, 0]:
        ax.set_ylabel("")
    ax.legend_.remove()
for ax in axes.flatten()[len(cell_types):]:
    ax.remove()
plt.suptitle("Proportions of cells in neighborhood vs Factor2", fontsize=16, y=1.01)
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/ct_prop_vs_factor2_spline.png")
plt.show()


# Normalize each row to sum to 1
fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
    fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0)


# 4.4 Smooth Ribbon Compo: Factor2 --------------------------------------------

fitted_spline_end = fitted_spline_df.copy()
fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(axis=1)
fitted_spline_end_longer = fitted_spline_end.melt(
    id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_end")

fitted_spline_start = fitted_spline_df.copy()
fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(1, axis=1).cumsum(axis=1)
fitted_spline_start_longer = fitted_spline_start.melt(
    id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_start")
fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

fitted_spline_merged = pd.merge(
    fitted_spline_start_longer,
    fitted_spline_end_longer,
    on=["Factor2_ulm", "cell_type"],
    how="outer"
)

plt.figure(figsize=(12, 8))
for cell_type in cell_types:
    df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == cell_type]
    plt.fill_between(
        df_ct["Factor2_ulm"],
        df_ct["fitted_prop_start"],
        df_ct["fitted_prop_end"],
        color=my_palette[cell_type],
        label=cell_type,
        alpha=0.8,
        linewidth=0
    )
plt.xlabel("Factor2 cropped to [-4, 8]")
plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
plt.title("Cell Type composition of 2-hop Neighborhoods vs Factor2")
plt.legend(bbox_to_anchor=(1.05, 0.85), loc='upper left', title="Cell Type")
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/ct_prop_vs_factor2_ribbon.png")
plt.show()


# 4.5 Smooth spline compo: Factor2 by disease --------------------------------------------

# disease = "CTRL"
for disease in diseases:
    # Get all cores for this disease
    cores_for_disease = gdf_cell_boundaries[gdf_cell_boundaries["disease_id"] == disease]["TMA_ID"].unique().tolist()
    # Subset nh_compo_prop_f2_df to only those cores
    nh_compo_prop_f2_df_disease = nh_compo_prop_f2_df[nh_compo_prop_f2_df["TMA_ID"].isin(cores_for_disease)].copy()
    if nh_compo_prop_f2_df_disease.empty:
        continue

    nh_compo_prop_f2_df_longer = nh_compo_prop_f2_df_disease.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types, var_name="cell_type", value_name="proportion"
    )

    x_seq = np.linspace(-4, 8, 100)  # Manually set the x maximal values

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for ax, cell_type in zip(axes.flatten(), cell_types):
        aaa = (nh_compo_prop_f2_df_longer[nh_compo_prop_f2_df_longer["cell_type"] == cell_type].
               sort_values(by="Factor2_ulm"))
        x = aaa["Factor2_ulm"].values
        y = aaa["proportion"].values
        if len(x) < 4:
            ax.set_visible(False)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[cell_type] = prop_end_fit

        sns.scatterplot(
            data=aaa,
            y="proportion",
            x="Factor2_ulm",
            hue="cell_type",
            palette=my_palette,
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, prop_end_fit, label='Fitted Spline', color="black", linewidth=1.5)
        ax.set_title(cell_type)
        if ax not in axes[-1]:
            ax.set_xlabel("")
            ax.xaxis.set_visible(False)
        if ax not in axes[:, 0]:
            ax.set_ylabel("")
        if hasattr(ax, "legend_") and ax.legend_:
            ax.legend_.remove()
    for ax in axes.flatten()[len(cell_types):]:
        ax.remove()
    plt.suptitle(f"Proportions of cells in neighborhood vs Factor2_ulm\nDisease: {disease}", fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_factor2/by_disease/" +
                f"ct_prop_vs_factor2_{str(disease).replace(' ', '_')}_spline.png")
    plt.show()
    # plt.show()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
        fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0
    )

    # 4.6 Smooth Ribbon Compo: Factor2 by disease --------------------------------------------
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "cell_type"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for cell_type in cell_types:
        df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == cell_type]
        plt.fill_between(
            df_ct["Factor2_ulm"],
            df_ct["fitted_prop_start"],
            df_ct["fitted_prop_end"],
            color=my_palette[cell_type],
            label=cell_type,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel("Factor2")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(f"Cell Type composition of 2-hop Neighborhoods vs Factor2\nDisease: {disease}")
    plt.legend(bbox_to_anchor=(1.05, 0.85), loc='upper left', title="Cell Type")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_factor2/by_disease/" +
                f"ct_prop_vs_factor2_{str(disease).replace(' ', '_')}_ribbon.png")
    plt.show()
    # plt.show()


print("--------------------------------------------------")
print("----- Finished plotting Factor2 plot -------------")
print("--------------------------------------------------")


# 4.7 Smooth Ribbon Compo: all vars of interest -------------------------------------------

var = vars_of_interest[0]
# var = "jm_slope"

for var in vars_of_interest:

    if var == "jm_slope":  # Don't include jm_slope
        continue
    if var not in nh_prop_df.columns:
        continue
    nh_compo_prop_var_df = nh_prop_df.copy()
    nh_compo_prop_var_df = nh_compo_prop_var_df[~nh_compo_prop_var_df[var].isna()]
    nh_compo_prop_var_df = nh_compo_prop_var_df.sort_values(by=var)
    if nh_compo_prop_var_df.empty:
        continue
    os.makedirs(f"{fig_dir}vs_{var}/", exist_ok=True)

    nh_compo_prop_var_df_longer = nh_compo_prop_var_df.melt(
        id_vars=[var, "sample_cell_id"], value_vars=cell_types, var_name="cell_type", value_name="proportion"
    )

    var_min = nh_compo_prop_var_df_longer[var].min()
    var_max = nh_compo_prop_var_df_longer[var].max()
    # By cropping extremes
    # x_seq = np.linspace(var_min + (var_max - var_min) * 0.05,
    #                     var_max - (var_max - var_min) * 0.05, 100)
    # By quantile
    x_seq = np.linspace(nh_compo_prop_var_df_longer[var].quantile(0.05),
                        nh_compo_prop_var_df_longer[var].quantile(0.95), 100)

    fitted_spline_df = pd.DataFrame({var: x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for ax, cell_type in zip(axes.flatten(), cell_types):
        aaa = (nh_compo_prop_var_df_longer[nh_compo_prop_var_df_longer["cell_type"] == cell_type].
               sort_values(by=var))
        x = aaa[var].values
        y = aaa["proportion"].values
        if len(x) < 4:
            ax.set_visible(False)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[cell_type] = prop_end_fit

        sns.scatterplot(
            data=aaa,
            y="proportion",
            x=var,
            hue="cell_type",
            palette=my_palette,
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, prop_end_fit, label='Fitted Spline', color="black", linewidth=1.5)
        ax.set_title(cell_type)
        if ax not in axes[-1]:
            ax.set_xlabel("")
            ax.xaxis.set_visible(False)
        if ax not in axes[:, 0]:
            ax.set_ylabel("")
        if hasattr(ax, "legend_") and ax.legend_:
            ax.legend_.remove()
    for ax in axes.flatten()[len(cell_types):]:
        ax.remove()
    plt.suptitle(f"Proportions of cells in neighborhood vs {var}", fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_{var}/ct_prop_vs_{var}.png")
    plt.show()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
        fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=[var], var_name="cell_type", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=[var], var_name="cell_type", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=[var, "cell_type"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for cell_type in cell_types:
        df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == cell_type]
        plt.fill_between(
            df_ct[var],
            df_ct["fitted_prop_start"],
            df_ct["fitted_prop_end"],
            color=my_palette[cell_type],
            label=cell_type,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"{var} cropped to [{x_seq.min():.2f}, {x_seq.max():.2f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(f"Cell Type composition of 2-hop Neighborhoods vs {var}")
    plt.legend(bbox_to_anchor=(1.05, 0.85), loc='upper left', title="Cell Type")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_{var}/ct_prop_vs_{var}_ribbon.png")
    plt.show()


# 4.8 Smooth Ribbon vs Factor2 by root_cell_type -------------------

os.makedirs(f"{fig_dir}vs_factor2/by_root_cell_type/", exist_ok=True)


for root_cell_type in cell_types:
    nh_compo_prop_f2_df_cond = nh_compo_prop_f2_df[nh_compo_prop_f2_df["root_cell_type"] == root_cell_type].copy()
    if nh_compo_prop_f2_df_cond.empty:
        continue

    nh_compo_prop_f2_df_longer = nh_compo_prop_f2_df_cond.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types,
        var_name="cell_type", value_name="proportion"
    )

    # Set x_seq for the ribbon plot ---------------------
    # factor2_min = nh_compo_prop_f2_df_longer["Factor2_ulm"].min()
    # factor2_max = nh_compo_prop_f2_df_longer["Factor2_ulm"].max()
    # x_seq = np.linspace(factor2_min + (factor2_max - factor2_min) * 0.05,
    #                     factor2_max - (factor2_min) * 0.05, 100)
    x_seq = np.linspace(-4, 8, 100)  # Chosen method: manually set the x maximal values

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    for cell_type in cell_types:
        aaa = (nh_compo_prop_f2_df_longer[nh_compo_prop_f2_df_longer["cell_type"] == cell_type]
               .sort_values(by="Factor2_ulm"))
        x = aaa["Factor2_ulm"].values
        y = aaa["proportion"].values
        if len(x) < 4:
            fitted_spline_df[cell_type] = np.full_like(x_seq, np.nan)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[cell_type] = prop_end_fit

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
        fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "cell_type"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for cell_type in cell_types:
        df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == cell_type]
        plt.fill_between(
            df_ct["Factor2_ulm"],
            df_ct["fitted_prop_start"],
            df_ct["fitted_prop_end"],
            color=my_palette[cell_type],
            label=cell_type,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"Factor2 cropped to [{x_seq.min():.1f}, {x_seq.max():.1f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title("2-hop Neighborhood Composition vs " +
              f"Factor2\nRoot Cell Type: {root_cell_type} ({nh_compo_prop_f2_df_cond.shape[0]} cells)")
    plt.legend(bbox_to_anchor=(1.05, 0.85), loc='upper left', title="Cell Type")
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir}vs_factor2/by_root_cell_type/ct_prop_vs_factor2_{root_cell_type.replace('/', '_')}_ribbon.png")
    plt.show()


print("--------------------------------------------------")
print("----- Finished plotting all vars plot ------------")
print("--------------------------------------------------")


# 4.9 Self-Excluding Smooth Ribbon Compo: Factor2 (not stratified by root_cell_type) ------

# Copy and self-exclude: set the count of the root_cell_type to zero for each row
nh_compo_prop_f2_se_df = nh_compo_prop_f2_df.copy()
for idx, row in nh_compo_prop_f2_se_df.iterrows():
    nh_compo_prop_f2_se_df.at[idx, row['root_cell_type']] = 0

# Compute proportions for self-excluding neighborhood composition
denom = nh_compo_prop_f2_se_df[cell_types].sum(axis=1)
non_empty_rows = denom != 0
nh_compo_prop_f2_se_df = nh_compo_prop_f2_se_df.loc[non_empty_rows]
nh_compo_prop_f2_se_df[cell_types] = nh_compo_prop_f2_se_df[cell_types].div(
    denom.replace(0, np.nan), axis=0
)

# Melt to long format
nh_compo_prop_f2_se_long = nh_compo_prop_f2_se_df.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types,
    var_name="cell_type", value_name="proportion"
)

# Different choices for restricted Factor2 space --------------------
# factor2_min = nh_compo_prop_f2_se_long["Factor2_ulm"].min()
# factor2_max = nh_compo_prop_f2_se_long["Factor2_ulm"].max()
# x_seq = np.linspace(
#     factor2_min + (factor2_max - factor2_min) * 0.05,
#     factor2_max - (factor2_max - factor2_min) * 0.05, 100
# )
x_seq = np.linspace(-4, 8, 100)  # Chosen method: manually set the x maximal values

fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

for cell_type in cell_types:
    aaa = (nh_compo_prop_f2_se_long[nh_compo_prop_f2_se_long["cell_type"] == cell_type]
           .sort_values(by="Factor2_ulm"))
    x = aaa["Factor2_ulm"].values
    y = aaa["proportion"].values
    if len(x) < 4:
        fitted_spline_df[cell_type] = np.full_like(x_seq, np.nan)
        continue
    knots = np.linspace(x[0], x[-1], num=10)[1:-1]
    try:
        spline = LSQUnivariateSpline(x, y, t=knots, k=3)
        prop_end_fit = np.clip(spline(x_seq), 0, 1)
    except Exception:
        prop_end_fit = np.full_like(x_seq, np.nan)
    fitted_spline_df[cell_type] = prop_end_fit

# Normalize each row to sum to 1
fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
    fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0
)

# Ribbon plot
fitted_spline_end = fitted_spline_df.copy()
fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(axis=1)
fitted_spline_end_longer = fitted_spline_end.melt(
    id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_end")

fitted_spline_start = fitted_spline_df.copy()
fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(1, axis=1).cumsum(axis=1)
fitted_spline_start_longer = fitted_spline_start.melt(
    id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_start")
fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

fitted_spline_merged = pd.merge(
    fitted_spline_start_longer,
    fitted_spline_end_longer,
    on=["Factor2_ulm", "cell_type"],
    how="outer"
)

plt.figure(figsize=(12, 8))
for cell_type in cell_types:
    df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == cell_type]
    plt.fill_between(
        df_ct["Factor2_ulm"],
        df_ct["fitted_prop_start"],
        df_ct["fitted_prop_end"],
        color=my_palette[cell_type],
        label=cell_type,
        alpha=0.8,
        linewidth=0
    )
plt.xlabel(f"Factor2 cropped to [{x_seq.min():.2f}, {x_seq.max():.2f}]")
plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
plt.title("2-hop Neighborhood Composition vs Factor2 (Self-Excluded)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Cell Type")
plt.tight_layout()
plt.savefig(f"{fig_dir}self_excluding/ct_prop_vs_factor2_ribbon_self_excluded.png")
plt.show()


# 4.10 Self-Excluding Smooth Ribbon Compo: Factor2 by root_cell_type -------------------

os.makedirs(f"{fig_dir}self_excluding/by_root_cell_type/", exist_ok=True)

# root_cell_type = cell_types[0]

for root_cell_type in cell_types:

    # Subset to rows with this root_cell_type
    nh_compo_prop_f2_se_df = nh_compo_prop_f2_df[
        nh_compo_prop_f2_df["root_cell_type"] == root_cell_type
    ].copy()
    if nh_compo_prop_f2_se_df.empty:
        continue

    # Self-exclude: set the count of the root_cell_type to zero for each row
    for idx, row in nh_compo_prop_f2_se_df.iterrows():
        nh_compo_prop_f2_se_df.at[idx, root_cell_type] = 0

    # Compute proportions for self-excluding neighborhood composition
    denom = nh_compo_prop_f2_se_df[cell_types].sum(axis=1)
    non_empty_rows = denom != 0
    nh_compo_prop_f2_se_df = nh_compo_prop_f2_se_df.loc[non_empty_rows]
    nh_compo_prop_f2_se_df[cell_types] = nh_compo_prop_f2_se_df[cell_types].div(
        denom.replace(0, np.nan), axis=0
    )

    # Melt to long format
    nh_compo_prop_f2_se_long = nh_compo_prop_f2_se_df.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=cell_types,
        var_name="cell_type", value_name="proportion"
    )

    factor2_min = nh_compo_prop_f2_se_long["Factor2_ulm"].min()
    factor2_max = nh_compo_prop_f2_se_long["Factor2_ulm"].max()
    x_seq = np.linspace(
        factor2_min + (factor2_max - factor2_min) * 0.05,
        factor2_max - (factor2_max - factor2_min) * 0.05, 100
    )
    # x_seq = np.linspace(-4, 8, 100)

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for cell_type in cell_types:
        aaa = (nh_compo_prop_f2_se_long[nh_compo_prop_f2_se_long["cell_type"] == cell_type]
               .sort_values(by="Factor2_ulm"))
        x = aaa["Factor2_ulm"].values
        y = aaa["proportion"].values
        if len(x) < 4:
            fitted_spline_df[cell_type] = np.full_like(x_seq, np.nan)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[cell_type] = prop_end_fit

        # Plot scatter and spline for each cell type in its own subplot
        ax = axes.flatten()[cell_types.index(cell_type)]
        sns.scatterplot(
            data=aaa,
            y="proportion",
            x="Factor2_ulm",
            hue="cell_type",
            palette=my_palette,
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, fitted_spline_df[cell_type], label='Fitted Spline', color="black", linewidth=1.5)
        ax.set_title(cell_type)
    if ax not in axes[-1]:
        ax.set_xlabel("")
        ax.xaxis.set_visible(False)
    if ax not in axes[:, 0]:
        ax.set_ylabel("")
    if hasattr(ax, "legend_") and ax.legend_:
        ax.legend_.remove()
    for ax in axes.flatten()[len(cell_types):]:
        ax.remove()
    plt.suptitle(
        "Self-Excluded 2-hop Neighborhood Composition vs Factor2" +
        f"\nRoot Cell Type: {root_cell_type} ({nh_compo_prop_f2_se_df.shape[0]} cells)",
        fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir}self_excluding/by_root_cell_type/" +
        f"ct_prop_vs_factor2_{root_cell_type.replace('/', '_').replace(' ', '_')}_self_excluded_spline.png"
    )
    plt.show()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, cell_types] = fitted_spline_df.loc[:, cell_types].div(
        fitted_spline_df.loc[:, cell_types].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, cell_types] = fitted_spline_end[cell_types].cumsum(
        axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, cell_types] = fitted_spline_start[cell_types].shift(
        1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="cell_type", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "cell_type"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for ct in cell_types:
        df_ct = fitted_spline_merged[fitted_spline_merged["cell_type"] == ct]
        plt.fill_between(
            df_ct["Factor2_ulm"],
            df_ct["fitted_prop_start"],
            df_ct["fitted_prop_end"],
            color=my_palette[ct],
            label=ct,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"Factor2 cropped to [{x_seq.min():.1f}, {x_seq.max():.1f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(
        "2-hop Neighborhood Composition (Self-Excluded) vs " +
        f"Factor2\nRoot Cell Type: {root_cell_type} ({nh_compo_prop_f2_se_df.shape[0]} cells)")
    plt.legend(bbox_to_anchor=(1.05, 0.85),
               loc='upper left', title="Cell Type")
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir}self_excluding/by_root_cell_type/" +
        f"ct_prop_vs_factor2_ribbon_{root_cell_type.replace('/', '_').replace(' ', '_')}_self_excluded.png")
    plt.show()


print("--------------------------------------------------")
print("----- Finished plotting self-excluding plot ------")
print("--------------------------------------------------")


