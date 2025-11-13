#
# Created on Wed Sep 03 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script analyzes and visualizes **niche composition** in 2-hop spatial neighborhoods
# from Xenium kidney biopsy data, focusing on spatial structure and clinical correlations.
#
# Datasets used/loaded:
# - Cell boundaries: geodataframes/cell_boundaries.parquet
# - Nucleus boundaries: geodataframes/nucleus_boundaries.parquet
# - Xenium dataset with niche assignments: xenium_with_niches.h5ad
# - 2-hop neighborhood niche composition: neighborhoods/nh_2hop_niche_compo.csv
# - Clinical annotations: scVI_kidney_biopsy_with_clinical.h5ad
#
# Datasets generated/saved:
# - Neighborhood niche composition with clinical and latent variables:
#       nh_niche_compo_df_w_vars.csv
#
# Plots generated (by output folder):
#
# heatmap_compo.png:
#   - Global heatmap of 2-hop neighborhood niche composition (proportional abundances per root niche)
#
# heatmap_per_core_niche/:
#   - Heatmaps of 2-hop neighborhood niche composition per tissue core (TMA_ID)
#
# heatmap_per_disease_niche/:
#   - Heatmaps of 2-hop neighborhood niche composition stratified by disease group
#
# self_excluding/:
#   - Self-excluded heatmap (root niche excluded from its own neighborhood)
#
# vs_factor2/:
#   - Stacked barplot: niche proportions vs Factor2_ulm
#   - Segment plot: niche composition transitions vs Factor2_ulm
#   - Violin plot: distribution of Factor2_ulm by root niche
#   - Scatter + spline plots (faceted by niche): smoothed niche proportions vs Factor2_ulm
#
# The analysis pipeline:
# 1. Load and harmonize Xenium and spatial boundary data.
# 2. Integrate niche assignments from AnnData into GeoDataFrame.
# 3. Identify cores lacking niche annotations and visualize them.
# 4. Compute 2-hop neighborhood niche composition globally, per core, and per disease.
# 5. Produce self-excluding versions (excluding root niche from neighborhood counts).
# 6. Merge with clinical and latent variables from scVI-derived annotations.
# 7. Explore quantitative trends vs latent/clinical variables (Factor2_ulm), using
#    stacked bars, segment plots, violin plots, and spline-smoothed trends.
#
# The script provides a comprehensive spatial and clinical characterization of local
# niche composition and structure within Xenium kidney biopsy data.
# --------------------------------------------------------------------------------------------


import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import geopandas as gpd
import h5py
import seaborn as sns
from scipy.interpolate import LSQUnivariateSpline
import os

max_threshold = 10

# Set directories ---------------------------

data_dir = "/data/zd775156/kidney_biopsy/data/"
xenium_path = "/data/Kramann_lab/From_Sidrah/new_data/output-XETG00229__0042348__Region_3__20250416__104148/"
res_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
fig_dir_whole_core = '/data/Kramann_lab/kidney_biopsy/results/xenium/leon/niche_spatial_plot/whole_core/'
fig_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/neighborhoods/2hop_nh_niche_compo_clustering/"


# Load xenium with boundaries ---------------------------------------------------
gdf_cell_boundaries = gpd.read_parquet(res_dir+"geodataframes/cell_boundaries.parquet")
gdf_nucleus_boundaries = gpd.read_parquet(res_dir+"geodataframes/nucleus_boundaries.parquet")

gdf_cell_boundaries.index.names = [None]


# Load niches adata ----------------------------------------------------------------

# Remove uns object to make adata readable
with h5py.File("/data/Kramann_lab/kidney_biopsy/results/xenium/leon/xenium_with_niches.h5ad", "r+") as f:
    del f["uns"]["spatial_neighbors"]["params"]

adata_niches = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/results/xenium/leon/xenium_with_niches.h5ad"
)

# Define variable from gdf ----------------------------

cell_types = gdf_cell_boundaries["cell_type_manual"].cat.categories
samples = gdf_cell_boundaries["sample"].unique().tolist()
cores = gdf_cell_boundaries["TMA_ID"].unique().tolist()
adata_cores = adata_niches.obs["TMA_ID"].unique().tolist()

set(adata_cores) - set(cores)
cores_intersect = list(set(adata_cores) & set(cores))

# Create custom palette for cell_types ------------------
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
gdf_cell_boundaries = gdf_cell_boundaries.merge(
    adata_niches.obs[["sample_cell_id", "latent_leiden_0.4"]],
    left_on="sample_cell_id",
    right_on="sample_cell_id",
    how="left"
)

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
adata_niche_wo_data.obs["latent_leiden_0.4"]

# 2. Plot nh_compo w.r.t. Factor 2 ---------------------------------------------------

# Load clinical variables from adata ------------------------------
adata = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/data/scVI_kidney_biopsy_with_clinical.h5ad",
    backed="r"
)
adata_obs = adata.obs.drop_duplicates(
    subset=["ID_number", "classification", "jm_slope"])
adata_mapping = adata_obs[
    ["ID_number", "jm_slope", "classification", "proteinuria_group", "diabetes",
     "hypertension_clinical", "Proteinuria over time", "egfr_group"]].set_index("ID_number")
adata_mapping = adata_mapping[~adata_mapping.index.isna()]
adata.file.close()


# Load nh_compo_df -------
nh_niche_compo_df = pd.read_csv(
    res_dir + "neighborhoods/nh_2hop_niche_compo.csv",
    index_col=0
)
nh_niche_compo_df.index = nh_niche_compo_df["sample_cell_id"]

# Handle non-unique index of nh_niche_compo_df ---------------------------------

# Warning: the index of nh_niche_compo_df is not unique, which may cause issues later on.
if len(nh_niche_compo_df.index.unique()) != len(nh_niche_compo_df.index):
    print("Warning: the index of nh_niche_compo_df is not unique, which may cause issues later on.")
    # Simple fix: only take the first occurrence of each index
    nh_niche_compo_df = nh_niche_compo_df[~nh_niche_compo_df.index.duplicated(keep='first')]

# 1. Compute 2-hop nh niche proportions ---------------------------------------

niches = sorted(adata_niches.obs["latent_leiden_0.4"].unique().tolist(),
                key=lambda x: int(x) if x.isdigit() else float('inf'))

nh_niche_compo_prop_df = nh_niche_compo_df.copy()
nh_niche_compo_prop_df[niches] = nh_niche_compo_prop_df[niches].div(
    nh_niche_compo_prop_df[niches].sum(axis=1), axis=0
)

# Assert the sum by row is close to one
assert all(np.isclose(nh_niche_compo_prop_df[niches].sum(axis=1), 1))

# Group the rows by cell_id and average the proportions
nh_niche_compo_prop_df = nh_niche_compo_prop_df.drop(
    columns=["sample_cell_id", "TMA_ID", "disease"]).groupby('root_niche').mean().fillna(0)


# Plot a heatmap of the values
plt.figure(figsize=(10, 8))
sns.heatmap(nh_niche_compo_prop_df[niches],
            cmap="viridis", cbar_kws={'label': 'Proportion'})
plt.title("Cell Type Composition (2-hop Neighborhood)")
plt.xlabel("Cell Types in the 2-hop Neighborhood")
plt.ylabel("Cell")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{fig_dir}heatmap_compo.png")
plt.show()


# --- Niche composition analysis (analogous to cell_type composition) ---

# 2. Compute 2-hop nh niche proportions per core ---------------------------------------

vmax_niche = nh_niche_compo_df[niches].max().max()

for core in cores:

    nh_niche_core = nh_niche_compo_df[nh_niche_compo_df["TMA_ID"] == core].copy()

    # Compute proportions
    nh_niche_prop_core = nh_niche_core.copy()
    nh_niche_prop_core[niches] = nh_niche_prop_core[niches].div(
        nh_niche_prop_core[niches].sum(axis=1), axis=0
    )
    assert all(np.isclose(nh_niche_prop_core[niches].sum(axis=1)[nh_niche_prop_core[niches].sum(axis=1) > 0], 1))

    # Group by root_niche and average
    nh_niche_prop_core = nh_niche_prop_core.drop(
        columns=["sample_cell_id", "TMA_ID", "disease"], errors="ignore"
    ).groupby('root_niche').mean().fillna(0)

    # Plot heatmap for this core
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        nh_niche_prop_core[niches], cmap="viridis",
        vmax=1,
        cbar_kws={'label': 'Proportion'})
    plt.title(f"Niche Composition (2-hop Neighborhood)\n{core}")
    plt.xlabel("Niches in the 2-hop Neighborhood")
    plt.ylabel("Root Niche")
    plt.tight_layout()
    os.makedirs(f"{fig_dir}heatmap_per_core_niche/", exist_ok=True)
    plt.savefig(f"{fig_dir}heatmap_per_core_niche/{core.replace(' ', '_')}.png")
    plt.close()

# 3. Compute 2-hop nh niche proportions per disease ---------------------------------------

diseases = nh_niche_compo_df["disease"].dropna().unique().tolist()

for disease in diseases:

    nh_niche_disease = nh_niche_compo_df[nh_niche_compo_df["disease"] == disease].copy()
    nh_niche_prop_disease = nh_niche_disease.copy()
    nh_niche_prop_disease[niches] = nh_niche_prop_disease[niches].div(
        nh_niche_prop_disease[niches].sum(axis=1), axis=0
    )
    assert all(np.isclose(nh_niche_prop_disease[niches].sum(axis=1)[nh_niche_prop_disease[niches].sum(axis=1) > 0], 1))

    # Group by root_niche and average
    nh_niche_prop_disease = nh_niche_prop_disease.drop(
        columns=["sample_cell_id", "TMA_ID", "disease"], errors="ignore"
    ).groupby('root_niche').mean().fillna(0)

    # Plot heatmap for this disease
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        nh_niche_prop_disease[niches], cmap="viridis",
        vmax=1,
        cbar_kws={'label': 'Proportion'})
    plt.title(f"Niche Composition (2-hop Neighborhood)\nDisease: {disease}")
    plt.xlabel("Niches in the 2-hop Neighborhood")
    plt.ylabel("Root Niche")
    # plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    os.makedirs(f"{fig_dir}heatmap_per_disease_niche/", exist_ok=True)
    plt.savefig(
        f"{fig_dir}heatmap_per_disease_niche/{str(disease).replace(' ', '_')}.png")
    # plt.close()

# 4. Heatmap of Self-Excluding Niche Composition --------------------------

nh_niche_prop_se_df = nh_niche_compo_df.copy()
# For each row, set the count of the root_niche to zero
for idx, row in nh_niche_prop_se_df.iterrows():
    nh_niche_prop_se_df.at[idx, row['root_niche']] = 0

# Compute proportions for self-excluding neighborhood composition
denom = nh_niche_prop_se_df[niches].sum(axis=1)
non_empty_rows = denom != 0
nh_niche_prop_se_df = nh_niche_prop_se_df.loc[non_empty_rows]
nh_niche_prop_se_df[niches] = nh_niche_prop_se_df[niches].div(
    denom.replace(0, np.nan), axis=0)

assert all(np.isclose(nh_niche_prop_se_df[niches].sum(axis=1), 1))

# Remove all columns after 'disease'
nh_niche_prop_se_df = nh_niche_prop_se_df[["sample_cell_id", "TMA_ID", "disease", "root_niche"] + niches]

# Group by root_niche and average
nh_niche_prop_se_df = nh_niche_prop_se_df.drop(
    columns=["sample_cell_id", "TMA_ID", "disease"], errors="ignore"
).groupby('root_niche').mean()

# Convert data type to float
nh_niche_prop_se_df = nh_niche_prop_se_df.astype(float)
# Set the diagonal values to NaN
for ind_row, ind_col in zip(nh_niche_prop_se_df.index, nh_niche_prop_se_df.columns):
    nh_niche_prop_se_df.loc[ind_row, ind_col] = np.nan

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(nh_niche_prop_se_df, cmap="viridis",
            cbar_kws={'label': 'Proportion'})
plt.title("Niche Composition (2-hop Neighborhood, Self-Excluded)")
plt.xlabel("Niches in the 2-hop Neighborhood")
plt.ylabel("Root Niche")
plt.tight_layout()
os.makedirs(f"{fig_dir}self_excluding/", exist_ok=True)
plt.savefig(f"{fig_dir}self_excluding/niche_compo_self_excluded.png")
plt.show()

# 5. Compute nhood niche composition vs quantitative variables -------

# Define nh_niche_prop_df -----------------
nh_niche_prop_df = nh_niche_compo_df.copy()
nh_niche_prop_df[niches] = nh_niche_prop_df[niches].div(
    nh_niche_prop_df[niches].sum(axis=1), axis=0)

# Assert no NaN indices
assert nh_niche_prop_df.index.isna().sum() == 0

# Load adata with clinical variables ------------------------------
adata = sc.read_h5ad(
    "/data/Kramann_lab/kidney_biopsy/data/scVI_kidney_biopsy_with_clinical.h5ad",
    backed="r")
adata_obs = adata.obs.drop_duplicates(
    subset=["ID_number", "classification", "jm_slope"])
adata.file.close()

adata_mapping = adata_obs[
    ["ID_number", "jm_slope", "classification", "proteinuria_group", "diabetes",
     "hypertension_clinical", "Proteinuria over time", "egfr_group"]].set_index("ID_number")
adata_mapping = adata_mapping[~adata_mapping.index.isna()]
adata_mapping["ID_number"] = adata_mapping["ID_number"].astype(str)
adata_mapping.index.names = [None]

# Add clinical variables from single cell data -------------------------
nh_niche_prop_df.loc[:, "ID_number"] = [
    f.split(" ")[-2] for f in nh_niche_prop_df["TMA_ID"]]
nh_niche_prop_df = nh_niche_prop_df.merge(
    adata_mapping, left_on="ID_number", right_on="ID_number", how="left"
)
nh_niche_prop_df = nh_niche_prop_df.set_index(nh_niche_prop_df.index)

del nh_niche_prop_df["ID_number"]

# Workaround to the bug in adata mapping --------------------------------
gdf_cell_boundaries.index = gdf_cell_boundaries["sample_cell_id"]
nh_niche_prop_df.index = nh_niche_prop_df["sample_cell_id"]
nh_niche_prop_df = nh_niche_prop_df[~nh_niche_prop_df.index.duplicated(keep='first')]
gdf_cell_boundaries = gdf_cell_boundaries[~gdf_cell_boundaries.index.duplicated(keep='first')]

# Add clinical variables from gdf_boundaries ----------------------------------
vars_from_gdf = ["egfr_at_time_of_biopsy",
                 "proteinuria_around_biopsy protein to creatine ratio",
                 *["Factor"+str(i)+"_ulm" for i in range(1, 6)]]
for var in vars_from_gdf:
    nh_niche_prop_df[var] = nh_niche_prop_df.index.map(gdf_cell_boundaries[var])

vars_of_interest = vars_from_gdf + ["jm_slope"]

# Save nh_niche_prop_df with variables ------------- ------
nh_niche_prop_df.to_csv(f"{fig_dir}nh_niche_compo_df_w_vars.csv")

# 5.1 Compute composition as a function of Factor2 -------------
nh_niche_compo_prop_f2_df = nh_niche_prop_df.copy()
nh_niche_compo_prop_f2_df = nh_niche_compo_prop_f2_df[~nh_niche_compo_prop_f2_df["Factor2_ulm"].isna()]
nh_niche_compo_prop_f2_df = nh_niche_compo_prop_f2_df.sort_values(by="Factor2_ulm")

# 5.2 Plot barplot w.r.t. Factor2 ----------------

segments = nh_niche_compo_prop_f2_df[niches].to_dict()
segments = {
    niche: [f for f in inner_dict.values()]
    for niche, inner_dict in segments.items()
}
y_positions = nh_niche_compo_prop_f2_df["Factor2_ulm"].values
bottom = np.zeros(len(nh_niche_compo_prop_f2_df.index.tolist()))
height = 0.01

fig, ax = plt.subplots(figsize=(6, 4))
for name, vals in segments.items():
    ax.barh(y_positions, vals, left=bottom, height=height,
            alpha=0.2, label=name, color=cmap_niches.colors[int(name)])
    bottom += np.array(vals)
ax.set_yticks(y_positions)
ax.set_yticklabels("")
plt.tight_layout()
os.makedirs(f"{fig_dir}vs_factor2/", exist_ok=True)
plt.savefig(f"{fig_dir}vs_factor2/niche_prop_vs_factor2_bars.png")
plt.show()

# 5.3 Plot 2D points w.r.t Factor2 -------------------------

nh_niche_compo_prop_f2_df_end = nh_niche_compo_prop_f2_df.copy()
nh_niche_compo_prop_f2_df_end["sample_cell_id"] = nh_niche_compo_prop_f2_df_end.index
nh_niche_compo_prop_f2_df_end.loc[:, niches] = nh_niche_compo_prop_f2_df[niches].cumsum(axis=1)

nh_niche_compo_prop_f2_df_start = nh_niche_compo_prop_f2_df.copy()
nh_niche_compo_prop_f2_df_start["sample_cell_id"] = nh_niche_compo_prop_f2_df_start.index
nh_niche_compo_prop_f2_df_start[niches] = nh_niche_compo_prop_f2_df[niches].shift(1, axis=1).cumsum(axis=1)
nh_niche_compo_prop_f2_df_start[niches] = nh_niche_compo_prop_f2_df_start[niches].fillna(0)

nh_niche_compo_prop_f2_df_end_longer = nh_niche_compo_prop_f2_df_end.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches, var_name="niche", value_name="proportion_end"
)
assert np.allclose(nh_niche_compo_prop_f2_df_start[niches[0]], 0)
assert np.allclose(nh_niche_compo_prop_f2_df_end[niches[-1]], 1)

nh_niche_compo_prop_f2_df_start_longer = nh_niche_compo_prop_f2_df_start.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches,
    var_name="niche", value_name="proportion_start"
)
nh_niche_compo_prop_f2_df_longer_seg = pd.merge(
    nh_niche_compo_prop_f2_df_start_longer,
    nh_niche_compo_prop_f2_df_end_longer,
    on=["Factor2_ulm", "sample_cell_id", "niche"]
)

plt.figure(figsize=(12, 10))
for niche in niches:
    df_niche = nh_niche_compo_prop_f2_df_longer_seg[nh_niche_compo_prop_f2_df_longer_seg["niche"] == niche]
    plt.hlines(
        y=df_niche["Factor2_ulm"],
        xmin=df_niche["proportion_start"],
        xmax=df_niche["proportion_end"],
        color=cmap_niches.colors[int(niche)],
        label=niche,
        linewidth=0.2,
        alpha=0.4,
    )
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', title="Niche")
plt.xlabel("Proportion in 2-hop Neighborhood")
plt.ylabel("Factor2_ulm")
plt.title("Niche Proportion vs Factor2_ulm (Stacked Segments)")
plt.savefig(f"{fig_dir}vs_factor2/niche_prop_vs_factor2_segments.png")
plt.tight_layout()
plt.show()

# 5. Plot the distribution of Factor2 per root_niche (no nh_compo) -----------

plt.figure(figsize=(12, 8))
plt.grid(True, axis='y', linestyle='--', alpha=0.5)
sns.violinplot(
    data=nh_niche_compo_prop_f2_df,
    x="root_niche",
    y="Factor2_ulm",
    palette=[cmap_niches.colors[int(n)] for n in niches],
    order=niches,
    cut=0
)
plt.xlabel("Root Niche")
plt.ylabel("Factor2_ulm")
plt.title("Distribution of Factor2_ulm by Root Niche")
plt.tight_layout()
os.makedirs(f"{fig_dir}vs_factor2/", exist_ok=True)
plt.savefig(f"{fig_dir}vs_factor2/factor2_violin_by_root_niche.png")
plt.show()

# 5.4 Smooth spline compo: Factor2 --------------------------------------------

nh_niche_compo_prop_f2_df_longer = nh_niche_compo_prop_f2_df.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches, var_name="niche", value_name="proportion"
)

factor2_min = nh_niche_compo_prop_f2_df_longer["Factor2_ulm"].min()
factor2_max = nh_niche_compo_prop_f2_df_longer["Factor2_ulm"].max()
x_seq = np.linspace(-4, 8, 100)

fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
for ax, niche in zip(axes.flatten(), niches):
    aaa = (nh_niche_compo_prop_f2_df_longer[nh_niche_compo_prop_f2_df_longer["niche"] == niche].
           sort_values(by="Factor2_ulm"))
    x = aaa["Factor2_ulm"].values
    y = aaa["proportion"].values
    knots = np.linspace(x_seq[0], x_seq[-1], num=10)[1:-1]
    spline = LSQUnivariateSpline(x, y, t=knots, k=3)
    prop_end_fit = np.clip(spline(x_seq), 0, 1)
    fitted_spline_df[niche] = prop_end_fit
    sns.scatterplot(
        data=aaa,
        y="proportion",
        x="Factor2_ulm",
        hue="niche",
        palette=[cmap_niches.colors[int(niche)]],
        alpha=0.7,
        s=1.5,
        edgecolor=None,
        ax=ax
    )
    ax.plot(x_seq, prop_end_fit, label='Fitted Spline',
            color="black", linewidth=2)
    ax.set_title(f"Niche {niche}")
    ax.set_xlabel("Factor2")
    if ax not in axes[-1]:
        ax.set_xlabel("")
        ax.xaxis.set_visible(False)
    if ax not in axes[:, 0]:
        ax.set_ylabel("")
    if hasattr(ax, "legend_") and ax.legend_:
        ax.legend_.remove()
for ax in axes.flatten()[len(niches):]:
    ax.remove()
plt.suptitle("Proportions of niches in neighborhood vs Factor2", fontsize=16, y=1.01)
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/niche_prop_vs_factor2_spline.png")
plt.show()

# Normalize each row to sum to 1
fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
    fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
)

# 4.4 Smooth Ribbon Compo: Factor2 --------------------------------------------

fitted_spline_end = fitted_spline_df.copy()
fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(axis=1)
fitted_spline_end_longer = fitted_spline_end.melt(
    id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_end")

fitted_spline_start = fitted_spline_df.copy()
fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(1, axis=1).cumsum(axis=1)
fitted_spline_start_longer = fitted_spline_start.melt(
    id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_start")
fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

fitted_spline_merged = pd.merge(
    fitted_spline_start_longer,
    fitted_spline_end_longer,
    on=["Factor2_ulm", "niche"],
    how="outer"
)

plt.figure(figsize=(8, 6))
for niche in niches:
    df_niche = fitted_spline_merged[fitted_spline_merged["niche"] == niche]
    plt.fill_between(
        df_niche["Factor2_ulm"],
        df_niche["fitted_prop_start"],
        df_niche["fitted_prop_end"],
        color=cmap_niches.colors[int(niche)],
        label=niche,
        alpha=0.8,
        linewidth=0
    )
plt.xlabel("Factor2 cropped to [-4, 8]")
plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
plt.title("Niche composition of 2-hop Neighborhoods vs Factor2")
plt.legend(bbox_to_anchor=(1.01, 0.75), loc='upper left', title="Niche")
plt.tight_layout()
plt.savefig(f"{fig_dir}vs_factor2/niche_prop_vs_factor2_ribbon.png")
plt.show()

# 4.5 Smooth spline compo: Factor2 by disease (niche composition) --------------------------------------------

for disease in diseases:
    # Get all cores for this disease
    cores_for_disease = gdf_cell_boundaries[gdf_cell_boundaries["disease_id"]
                                            == disease]["TMA_ID"].unique().tolist()
    # Subset nh_niche_compo_prop_f2_df to only those cores
    nh_niche_compo_prop_f2_df_disease = nh_niche_compo_prop_f2_df[nh_niche_compo_prop_f2_df["TMA_ID"].isin(
        cores_for_disease)].copy()
    if nh_niche_compo_prop_f2_df_disease.empty:
        continue

    nh_niche_compo_prop_f2_df_longer = nh_niche_compo_prop_f2_df_disease.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches, var_name="niche", value_name="proportion"
    )

    x_seq = np.linspace(-4, 8, 100)  # Manually set the x maximal values

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for ax, niche in zip(axes.flatten(), niches):
        aaa = (nh_niche_compo_prop_f2_df_longer[nh_niche_compo_prop_f2_df_longer["niche"] == niche].
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
        fitted_spline_df[niche] = prop_end_fit

        sns.scatterplot(
            data=aaa,
            y="proportion",
            x="Factor2_ulm",
            hue="niche",
            palette=[cmap_niches.colors[int(niche)]],
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, prop_end_fit, label='Fitted Spline',
                color="black", linewidth=1.5)
        ax.set_title(f"Niche {niche}")
        if ax not in axes[-1]:
            ax.set_xlabel("")
            ax.xaxis.set_visible(False)
        if ax not in axes[:, 0]:
            ax.set_ylabel("")
        if hasattr(ax, "legend_") and ax.legend_:
            ax.legend_.remove()
    for ax in axes.flatten()[len(niches):]:
        ax.remove()
    plt.suptitle(
        f"Proportions of niches in neighborhood vs Factor2_ulm\nDisease: {disease}", fontsize=16, y=1.01)
    plt.tight_layout()
    os.makedirs(f"{fig_dir}vs_factor2/by_disease/", exist_ok=True)
    plt.savefig(f"{fig_dir}vs_factor2/by_disease/" +
                f"niche_prop_vs_factor2_{str(disease).replace(' ', '_')}_spline.png")
    plt.close()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
        fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
    )

    # 4.6 Smooth Ribbon Compo: Factor2 by disease (niche composition) --------------------------------------------
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "niche"],
        how="outer"
    )

    plt.figure(figsize=(8, 6))
    for niche in niches:
        df_niche = fitted_spline_merged[fitted_spline_merged["niche"] == niche]
        plt.fill_between(
            df_niche["Factor2_ulm"],
            df_niche["fitted_prop_start"],
            df_niche["fitted_prop_end"],
            color=cmap_niches.colors[int(niche)],
            label=niche,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel("Factor2")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(
        f"Niche composition of 2-hop Neighborhoods vs Factor2\nDisease: {disease}")
    plt.legend(bbox_to_anchor=(1.01, 0.75),
               loc='upper left', title="Niche")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_factor2/by_disease/" +
                f"niche_prop_vs_factor2_{str(disease).replace(' ', '_')}_ribbon.png")
    plt.close()

print("--------------------------------------------------")
print("----- Finished plotting Factor2 plot -------------")
print("--------------------------------------------------")


# 5. Smooth Ribbon Compo: all vars of interest -------------------------------------------

for var in vars_of_interest:
    if var == "jm_slope":
        continue
    if var not in nh_niche_prop_df.columns:
        continue
    nh_niche_compo_prop_var_df = nh_niche_prop_df.copy()
    nh_niche_compo_prop_var_df = nh_niche_compo_prop_var_df[~nh_niche_compo_prop_var_df[var].isna()]
    nh_niche_compo_prop_var_df = nh_niche_compo_prop_var_df.sort_values(by=var)
    if nh_niche_compo_prop_var_df.empty:
        continue
    os.makedirs(f"{fig_dir}vs_{var}/", exist_ok=True)

    nh_niche_compo_prop_var_df_longer = nh_niche_compo_prop_var_df.melt(
        id_vars=[var, "sample_cell_id"], value_vars=niches, var_name="niche", value_name="proportion"
    )

    x_seq = np.linspace(
        nh_niche_compo_prop_var_df_longer[var].quantile(0.05),
        nh_niche_compo_prop_var_df_longer[var].quantile(0.95), 100
    )

    fitted_spline_df = pd.DataFrame({var: x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for ax, niche in zip(axes.flatten(), niches):
        aaa = (nh_niche_compo_prop_var_df_longer[nh_niche_compo_prop_var_df_longer["niche"] == niche].
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
        fitted_spline_df[niche] = prop_end_fit

        sns.scatterplot(
            data=aaa,
            y="proportion",
            x=var,
            hue="niche",
            palette=[cmap_niches.colors[int(niche)]],
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, prop_end_fit, label='Fitted Spline',
                color="black", linewidth=1.5)
        ax.set_title(f"Niche {niche}")
        if ax not in axes[-1]:
            ax.set_xlabel("")
            ax.xaxis.set_visible(False)
        if ax not in axes[:, 0]:
            ax.set_ylabel("")
        if hasattr(ax, "legend_") and ax.legend_:
            ax.legend_.remove()
    for ax in axes.flatten()[len(niches):]:
        ax.remove()
    plt.suptitle(
        f"Proportions of niches in neighborhood vs {var}", fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_{var}/niche_prop_vs_{var}.png")
    plt.close()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
        fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=[var], var_name="niche", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=[var], var_name="niche", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=[var, "niche"],
        how="outer"
    )

    plt.figure(figsize=(8, 6))
    for niche in niches:
        df_niche = fitted_spline_merged[fitted_spline_merged["niche"] == niche]
        plt.fill_between(
            df_niche[var],
            df_niche["fitted_prop_start"],
            df_niche["fitted_prop_end"],
            color=cmap_niches.colors[int(niche)],
            label=niche,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"{var} cropped to [{x_seq.min():.2f}, {x_seq.max():.2f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(f"Niche composition of 2-hop Neighborhoods vs {var}")
    plt.legend(bbox_to_anchor=(1.01, 0.75),
               loc='upper left', title="Niche")
    plt.tight_layout()
    plt.savefig(f"{fig_dir}vs_{var}/niche_prop_vs_{var}_ribbon.png")
    plt.close()

print("--------------------------------------------------")
print("----- Finished plotting niche composition plots --")
print("--------------------------------------------------")


# 4.8 Smooth Ribbon vs Factor2 by root_niche -------------------

os.makedirs(f"{fig_dir}vs_factor2/by_root_niche/", exist_ok=True)

for root_niche in niches:
    nh_niche_compo_prop_f2_df_cond = (nh_niche_compo_prop_f2_df[
        nh_niche_compo_prop_f2_df["root_niche"] == root_niche].copy())
    if nh_niche_compo_prop_f2_df_cond.empty:
        continue

    nh_niche_compo_prop_f2_df_longer = nh_niche_compo_prop_f2_df_cond.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches,
        var_name="niche", value_name="proportion"
    )

    x_seq = np.linspace(-4, 8, 100)  # Manually set the x maximal values

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    for niche in niches:
        aaa = (nh_niche_compo_prop_f2_df_longer[nh_niche_compo_prop_f2_df_longer["niche"] == niche]
               .sort_values(by="Factor2_ulm"))
        x = aaa["Factor2_ulm"].values
        y = aaa["proportion"].values
        if len(x) < 4:
            fitted_spline_df[niche] = np.full_like(x_seq, np.nan)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[niche] = prop_end_fit

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
        fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "niche"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for niche in niches:
        df_niche = fitted_spline_merged[fitted_spline_merged["niche"] == niche]
        plt.fill_between(
            df_niche["Factor2_ulm"],
            df_niche["fitted_prop_start"],
            df_niche["fitted_prop_end"],
            color=cmap_niches.colors[int(niche)],
            label=niche,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"Factor2 cropped to [{x_seq.min():.1f}, {x_seq.max():.1f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title("2-hop Neighborhood Niche Composition vs " +
              f"Factor2\nRoot Niche: {root_niche} ({nh_niche_compo_prop_f2_df_cond.shape[0]} cells)")
    plt.legend(bbox_to_anchor=(1.05, 0.85), loc='upper left', title="Niche")
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir}vs_factor2/by_root_niche/niche_prop_vs_factor2_{root_niche.replace('/', '_')}_ribbon.png")
    # plt.close()

print("--------------------------------------------------")
print("----- Finished plotting all vars plot ------------")
print("--------------------------------------------------")

# 4.9 Self-Excluding Smooth Ribbon Compo: Factor2 (not stratified by root_niche) ------

# Copy and self-exclude: set the count of the root_niche to zero for each row
nh_niche_compo_prop_f2_se_df = nh_niche_compo_prop_f2_df.copy()
for idx, row in nh_niche_compo_prop_f2_se_df.iterrows():
    nh_niche_compo_prop_f2_se_df.at[idx, row['root_niche']] = 0

# Compute proportions for self-excluding neighborhood composition
denom = nh_niche_compo_prop_f2_se_df[niches].sum(axis=1)
non_empty_rows = denom != 0
nh_niche_compo_prop_f2_se_df = nh_niche_compo_prop_f2_se_df.loc[non_empty_rows]
nh_niche_compo_prop_f2_se_df[niches] = nh_niche_compo_prop_f2_se_df[niches].div(
    denom.replace(0, np.nan), axis=0
)

# Melt to long format
nh_niche_compo_prop_f2_se_long = nh_niche_compo_prop_f2_se_df.melt(
    id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches,
    var_name="niche", value_name="proportion"
)

x_seq = np.linspace(-4, 8, 100)  # Manually set the x maximal values

fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

for niche in niches:
    aaa = (nh_niche_compo_prop_f2_se_long[nh_niche_compo_prop_f2_se_long["niche"] == niche]
           .sort_values(by="Factor2_ulm"))
    x = aaa["Factor2_ulm"].values
    y = aaa["proportion"].values
    if len(x) < 4:
        fitted_spline_df[niche] = np.full_like(x_seq, np.nan)
        continue
    knots = np.linspace(x[0], x[-1], num=10)[1:-1]
    try:
        spline = LSQUnivariateSpline(x, y, t=knots, k=3)
        prop_end_fit = np.clip(spline(x_seq), 0, 1)
    except Exception:
        prop_end_fit = np.full_like(x_seq, np.nan)
    fitted_spline_df[niche] = prop_end_fit

# Normalize each row to sum to 1
fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
    fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
)

# Ribbon plot
fitted_spline_end = fitted_spline_df.copy()
fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(axis=1)
fitted_spline_end_longer = fitted_spline_end.melt(
    id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_end")

fitted_spline_start = fitted_spline_df.copy()
fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(1, axis=1).cumsum(axis=1)
fitted_spline_start_longer = fitted_spline_start.melt(
    id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_start")
fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

fitted_spline_merged = pd.merge(
    fitted_spline_start_longer,
    fitted_spline_end_longer,
    on=["Factor2_ulm", "niche"],
    how="outer"
)

plt.figure(figsize=(12, 8))
for niche in niches:
    df_niche = fitted_spline_merged[fitted_spline_merged["niche"] == niche]
    plt.fill_between(
        df_niche["Factor2_ulm"],
        df_niche["fitted_prop_start"],
        df_niche["fitted_prop_end"],
        color=cmap_niches.colors[int(niche)],
        label=niche,
        alpha=0.8,
        linewidth=0
    )
plt.xlabel(f"Factor2 cropped to [{x_seq.min():.2f}, {x_seq.max():.2f}]")
plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
plt.title("2-hop Neighborhood Niche Composition vs Factor2 (Self-Excluded)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Niche")
plt.tight_layout()
plt.savefig(f"{fig_dir}self_excluding/niche_prop_vs_factor2_ribbon_self_excluded.png")
# plt.close()

# 4.10 Self-Excluding Smooth Ribbon Compo: Factor2 by root_niche -------------------

os.makedirs(f"{fig_dir}self_excluding/by_root_niche/", exist_ok=True)

for root_niche in niches:

    # Subset to rows with this root_niche
    nh_niche_compo_prop_f2_se_df = nh_niche_compo_prop_f2_df[
        nh_niche_compo_prop_f2_df["root_niche"] == root_niche
    ].copy()
    if nh_niche_compo_prop_f2_se_df.empty:
        continue

    # Self-exclude: set the count of the root_niche to zero for each row
    for idx, row in nh_niche_compo_prop_f2_se_df.iterrows():
        nh_niche_compo_prop_f2_se_df.at[idx, root_niche] = 0

    # Compute proportions for self-excluding neighborhood composition
    denom = nh_niche_compo_prop_f2_se_df[niches].sum(axis=1)
    non_empty_rows = denom != 0
    nh_niche_compo_prop_f2_se_df = nh_niche_compo_prop_f2_se_df.loc[non_empty_rows]
    nh_niche_compo_prop_f2_se_df[niches] = nh_niche_compo_prop_f2_se_df[niches].div(
        denom.replace(0, np.nan), axis=0
    )

    # Melt to long format
    nh_niche_compo_prop_f2_se_long = nh_niche_compo_prop_f2_se_df.melt(
        id_vars=["Factor2_ulm", "sample_cell_id"], value_vars=niches,
        var_name="niche", value_name="proportion"
    )

    factor2_min = nh_niche_compo_prop_f2_se_long["Factor2_ulm"].min()
    factor2_max = nh_niche_compo_prop_f2_se_long["Factor2_ulm"].max()
    x_seq = np.linspace(
        factor2_min + (factor2_max - factor2_min) * 0.05,
        factor2_max - (factor2_max - factor2_min) * 0.05, 100
    )

    fitted_spline_df = pd.DataFrame({"Factor2_ulm": x_seq})

    fig, axes = plt.subplots(nrows=6, ncols=4, figsize=(20, 22), sharey=True)
    for niche in niches:
        aaa = (nh_niche_compo_prop_f2_se_long[nh_niche_compo_prop_f2_se_long["niche"] == niche]
               .sort_values(by="Factor2_ulm"))
        x = aaa["Factor2_ulm"].values
        y = aaa["proportion"].values
        if len(x) < 4:
            fitted_spline_df[niche] = np.full_like(x_seq, np.nan)
            continue
        knots = np.linspace(x[0], x[-1], num=10)[1:-1]
        try:
            spline = LSQUnivariateSpline(x, y, t=knots, k=3)
            prop_end_fit = np.clip(spline(x_seq), 0, 1)
        except Exception:
            prop_end_fit = np.full_like(x_seq, np.nan)
        fitted_spline_df[niche] = prop_end_fit

        # Plot scatter and spline for each niche in its own subplot
        ax = axes.flatten()[niches.index(niche)]
        sns.scatterplot(
            data=aaa,
            y="proportion",
            x="Factor2_ulm",
            hue="niche",
            palette=[cmap_niches.colors[int(niche)]],
            alpha=0.7,
            s=2,
            edgecolor=None,
            ax=ax
        )
        ax.plot(x_seq, fitted_spline_df[niche], label='Fitted Spline', color="black", linewidth=1.5)
        ax.set_title(niche)
    if ax not in axes[-1]:
        ax.set_xlabel("")
        ax.xaxis.set_visible(False)
    if ax not in axes[:, 0]:
        ax.set_ylabel("")
    if hasattr(ax, "legend_") and ax.legend_:
        ax.legend_.remove()
    for ax in axes.flatten()[len(niches):]:
        ax.remove()
    plt.suptitle(
        "Self-Excluded 2-hop Neighborhood Niche Composition vs Factor2" +
        f"\nRoot Niche: {root_niche} ({nh_niche_compo_prop_f2_se_df.shape[0]} cells)",
        fontsize=16, y=1.01)
    plt.tight_layout()
    plt.savefig(
        f"{fig_dir}self_excluding/by_root_niche/" +
        f"niche_prop_vs_factor2_{root_niche.replace('/', '_').replace(' ', '_')}_self_excluded_spline.png"
    )
    # plt.close()
    plt.show()

    # Normalize each row to sum to 1
    fitted_spline_df.loc[:, niches] = fitted_spline_df.loc[:, niches].div(
        fitted_spline_df.loc[:, niches].sum(axis=1), axis=0
    )

    # Ribbon plot
    fitted_spline_end = fitted_spline_df.copy()
    fitted_spline_end.loc[:, niches] = fitted_spline_end[niches].cumsum(
        axis=1)
    fitted_spline_end_longer = fitted_spline_end.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_end")

    fitted_spline_start = fitted_spline_df.copy()
    fitted_spline_start.loc[:, niches] = fitted_spline_start[niches].shift(
        1, axis=1).cumsum(axis=1)
    fitted_spline_start_longer = fitted_spline_start.melt(
        id_vars=["Factor2_ulm"], var_name="niche", value_name="fitted_prop_start")
    fitted_spline_start_longer = fitted_spline_start_longer.fillna(0)

    fitted_spline_merged = pd.merge(
        fitted_spline_start_longer,
        fitted_spline_end_longer,
        on=["Factor2_ulm", "niche"],
        how="outer"
    )

    plt.figure(figsize=(12, 8))
    for n in niches:
        df_n = fitted_spline_merged[fitted_spline_merged["niche"] == n]
        plt.fill_between(
            df_n["Factor2_ulm"],
            df_n["fitted_prop_start"],
            df_n["fitted_prop_end"],
            color=cmap_niches.colors[int(n)],
            label=n,
            alpha=0.8,
            linewidth=0
        )
    plt.xlabel(f"Factor2 cropped to [{x_seq.min():.1f}, {x_seq.max():.1f}]")
    plt.ylabel("Smoothed Proportion in 2-hop Neighborhood")
    plt.title(
        "2-hop Neighborhood Niche Composition (Self-Excluded) vs " +
        f"Factor2\nRoot Niche: {root_niche} ({nh_niche_compo_prop_f2_se_df.shape[0]} cells)")
    plt.legend(bbox_to_anchor=(1.05, 0.85),
               loc='upper left', title="Niche")
    plt.tight_layout()
    os.makedirs(f"{fig_dir}self_excluding/by_root_niche/", exist_ok=True)
    plt.savefig(
        f"{fig_dir}self_excluding/by_root_niche/" +
        f"niche_prop_vs_factor2_ribbon_{root_niche.replace('/', '_').replace(' ', '_')}_self_excluded.png")
    # plt.close()

print("--------------------------------------------------")
print("----- Finished plotting self-excluding plot ------")
print("--------------------------------------------------")
