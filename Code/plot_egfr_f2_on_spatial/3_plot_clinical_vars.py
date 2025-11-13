#
# Created on Mon Jun 30 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script generates scatterplots relating LIANA factors and clinical variables to cell type
# abundances in Xenium kidney biopsy spatial transcriptomics data.
#
# Datasets used/loaded:
# - Xenium processed annotated data: xenium_processed.h5ad
#
# Main analyses and plots:
#
# 1. For each cell type, plots are generated showing:
#    - "Factor2" (mean LIANA factor per sample) vs. cell type abundance (proportion of cells)
#      - Stratified by disease_id (color-coded)
#      - Linear regression line with confidence interval, annotated with correlation, slope, and p-value
#      - Output location: {res_dir}/abundance_vs_factor2/Factor2_vs_{cell_type}_abundance.pdf
#
#    - "eGFR at time of biopsy" (clinical variable) vs. cell type abundance (proportion of cells)
#      - Stratified by disease_id (color-coded)
#      - Linear regression line with confidence interval, annotated with correlation, slope, and p-value
#      - Output location: {res_dir}/abundance_vs_egfr/egfr_vs_{cell_type}_abundance.pdf
#
# Additional processing:
# - Computes mean LIANA factor values per sample and per cell type
# - Computes cell type proportions per sample
# - Orders samples by clinical variables for potential faceting (not directly plotted here)
#
# Output:
# - For each cell type, two PDF scatterplots are saved:
#   - abundance_vs_factor2/Factor2_vs_{cell_type}_abundance.pdf
#   - abundance_vs_egfr/egfr_vs_{cell_type}_abundance.pdf
# --------------------------------------------------------------------------------------------

import seaborn as sns
import statsmodels.api as sm
from matplotlib import patches as mpatches
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

xenium_data_dir = "/data/Kramann_lab/From_Sidrah/new_data/h5ad_files_only_transcript/marker_genes/"
# data_dir = "/data/zd775156/kidney_biopsy/data/"
data_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
res_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/sample_level_plots/cell_type_manual/"


# --------------------------------------------------------------------------------------------
# Load Xenium data ---------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

xenium = sc.read_h5ad(data_dir + "xenium_processed.h5ad")


cell_types = sorted(xenium.obs["cell_type_manual"].cat.categories.tolist())

# Set xenium.obs["Factor2_ulm_mean"] to float
xenium.obs["Factor2_ulm_mean"] = xenium.obs["Factor2_ulm_mean"].astype(float)

# -------------------------------------------------------------------------------------------------
# Set plotting parameters -------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Get all unique cell types and assign them consistent colors
palette = sc.pl.palettes.default_102
colors = dict(zip(cell_types, palette[:len(cell_types)]))

# Precompute the color mapping
xenium.uns['cell_type_manual_colors'] = [colors[ct] for ct in cell_types]

# Create legend handles
legend_handles = [
    mpatches.Patch(color=colors[ct], label=ct) for ct in cell_types
]

# Define plotting grid ----------------------------------------------------------------------
disease_counts = xenium.obs[['disease_id', "TMA_ID"]].drop_duplicates().groupby('disease_id')['TMA_ID'].nunique()
disease_counts = disease_counts.reset_index()
print("Number of TMA_IDs for each disease_id:")
print(disease_counts)

# Create a subplot grid with 5 columns and 6 rows
ncols = 11
nrows = sum(disease_counts["TMA_ID"] // ncols + 1)

disease_ids = disease_counts["disease_id"].tolist()
tma_ids_by_disease = [xenium.obs.loc[xenium.obs["disease_id"] == d, "TMA_ID"].unique() for d in disease_ids]

n_tmas_by_disease = [len(tmas) for tmas in tma_ids_by_disease]
max_tmas = max(n_tmas_by_disease)

# Create a 2D array (nrows x ncols) with TMA_IDs -------------
tma_id_grid = []
for tma_ids in tma_ids_by_disease:
    for i in range(0, len(tma_ids), ncols):
        row = [None] * ncols
        chunk = tma_ids[i:i + ncols]
        row[:len(chunk)] = chunk
        tma_id_grid.append(row)
tma_id_grid = np.array(tma_id_grid, dtype=object)


# --------------------------------------------------------------------------------------------
# Plot sample-level Factor2 vs eGFR ----------------------------------------------------------
# --------------------------------------------------------------------------------------------

sample_level_res_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/sample_level_plots/"
sample_level_df = xenium.obs[["sample_ID", "Factor2_ulm_mean", "egfr_at_time_of_biopsy"]].drop_duplicates()

# Plot sample-level Factor2 vs eGFR
plt.figure(figsize=(6, 5))
sns.scatterplot(data=sample_level_df, y="Factor2_ulm_mean",
                x="egfr_at_time_of_biopsy", alpha=0.7)
plt.ylabel("Factor 2")
plt.xlabel("eGFR at Time of Biopsy")
plt.title("Factor 2 vs eGFR")
plt.tight_layout()
plt.savefig(f"{sample_level_res_dir}sample_level_Factor2_vs_egfr.pdf")


# --------------------------------------------------------------------------------------------
# Plot clinical vars, quantitative ------------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Define clinical variables to plot ---------------------------------------------------
clinical_vars = ["egfr_at_time_of_biopsy", "age_at_biopsy",
                 "proteinuria_around_biopsy protein to creatine ratio"]


# Define tma_id_grid ordered by values of clinical variables ------
tma_id_grid_ordered_dict_quant = {}

for clinical_var in clinical_vars:
    xenium_obs_ordered = xenium.obs[["disease_id", "TMA_ID", clinical_var]].sort_values(
        by=["disease_id", clinical_var],
        ascending=[True, False]
    ).drop_duplicates()

    # Get a list of TMA_IDs ordered by clinical_var
    tma_ids_by_disease_ordered = []
    for disease_id in xenium_obs_ordered["disease_id"].unique():
        tma_ids = xenium_obs_ordered.loc[
            xenium_obs_ordered["disease_id"] == disease_id, "TMA_ID"
        ].tolist()
        tma_ids_by_disease_ordered.append(tma_ids)

    # Create a 2D array (nrows x ncols) of TMA_IDs
    tma_id_grid_ordered = []
    for tma_ids in tma_ids_by_disease_ordered:
        for i in range(0, len(tma_ids), ncols):
            row = [None] * ncols
            chunk = tma_ids[i:i + ncols]
            row[:len(chunk)] = chunk
            tma_id_grid_ordered.append(row)
    tma_id_grid_ordered = np.array(tma_id_grid_ordered, dtype=object)

    tma_id_grid_ordered_dict_quant[clinical_var] = tma_id_grid_ordered

clinical_var = clinical_vars[0]

# --------------------------------------------------------------------------------------------
# Plot clinical vars, quantitative, not ordered by disease -----------------------------------
# --------------------------------------------------------------------------------------------

xenium_obs_ordered = xenium.obs[["disease_id", "TMA_ID", "egfr_at_time_of_biopsy"]].sort_values(
    by=["egfr_at_time_of_biopsy"],
    ascending=[False]
    ).drop_duplicates()

# Get a list of TMA_IDs ordered by clinical_var only (no disease_id sorting)
tma_ids_ordered_egfr = [xenium_obs_ordered["TMA_ID"].tolist()]

# Create a 2D array (nrows x ncols) of TMA_IDs
# Determine number of rows needed
nrows_by_egfr = len(tma_ids_ordered_egfr[0]) // ncols + (1 if len(tma_ids_ordered_egfr[0]) % ncols != 0 else 0)

tma_id_grid_ordered_egfr = []
tma_ids = tma_ids_ordered_egfr[0]
for i in range(nrows_by_egfr):
    row = [None] * ncols
    chunk = tma_ids[i * ncols:(i + 1) * ncols]
    row[:len(chunk)] = chunk
    tma_id_grid_ordered_egfr.append(row)
tma_id_grid_ordered_egfr = np.array(tma_id_grid_ordered_egfr, dtype=object)

disease_id_colormap = {
    "FSGS": "#1f77b4",  # blue
    "IgAN": "#ff7f0e",  # orange
    "MN": "#2ca02c",  # green
    "TIN": "#d62728",  # red
    "DN": "#9467bd",  # purple
    "CTRL": "#8c564b"  # brown
}

clinical_var = "egfr_at_time_of_biopsy"

# --------------------------------------------------------------------------------------------
# Plot LIANA factors at sample level ---------------------------------------------------------
# --------------------------------------------------------------------------------------------

factors = ["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"]


# Scatterplot Factor2 vs egfr ---------------------------------------------------------
# For all factors, one facet per factor (columns), one row per disease
xenium_long = xenium.obs.melt(
    id_vars=["disease_id", "egfr_at_time_of_biopsy"],
    value_vars=["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"],
    var_name="Factor",
    value_name="Factor_value"
)
# Assert that there are only NaN values for disease_id == CTRL
assert xenium_long.loc[xenium_long["disease_id"] == "CTRL", "Factor_value"].isna().all(), \
    "There are non-NaN values for CTRL disease_id in the LIANA factors."
# Remove CTRL disease_id from the DataFrame
xenium_long = xenium_long[xenium_long["disease_id"] != "CTRL"]


# Use the ordering from the eGFR values -------------------
tma_id_grid_ordered_factors = tma_id_grid_ordered_dict_quant["egfr_at_time_of_biopsy"]


# --------------------------------------------------------------------------------------------
# Plot LIANA factors at cell level ---------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Compute mean Factor2_ulm per TMA_ID and add as a new column to xenium.obs
factor2_ulm_mean_per_tma = xenium.obs.groupby("TMA_ID")["Factor2_ulm"].mean()
xenium.obs["Factor2_ulm_mean"] = xenium.obs["TMA_ID"].map(factor2_ulm_mean_per_tma)

# For each TMA_ID and for each cell type, compute the mean Factor2_ulm for that cell type
cell_type_factors = xenium.obs.groupby(["TMA_ID", "cell_type_manual"])["Factor2_ulm"].mean().unstack()
cell_type_factors = cell_type_factors.fillna(0)

# Pivot longer with TMA_ID and Factor2_ulm_mean_{cell_type} columns
cell_type_factors_long = cell_type_factors.reset_index().melt(
    id_vars=["TMA_ID"],
    var_name="cell_type_manual",
    value_name="Factor2_ulm_mean"
)
# Pivot wider to have one column for TMA_ID and one column for each cell type's Factor2_ulm_mean
cell_type_factors_wide = cell_type_factors_long.pivot(
    index="TMA_ID",
    columns="cell_type_manual",
    values="Factor2_ulm_mean"
).reset_index()


# --------------------------------------------------------------------------------------------
# Scatterplot of Factor2 vs cell type abundance ----------------------------------------------
# --------------------------------------------------------------------------------------------


# For each TMA_ID, compute the proportion of each cell type
cell_type_proportions = xenium.obs.groupby(
    ["TMA_ID", "cell_type_manual"]).size().unstack(fill_value=0)
cell_type_proportions = cell_type_proportions.div(
    cell_type_proportions.sum(axis=1), axis=0)

# Pivot longer into a dataframe with TMA_ID, cell_type, proportion
ct_props_long = cell_type_proportions.reset_index().melt(
    id_vars=["TMA_ID"],
    var_name="cell_type_manual",
    value_name="proportion"
)

# Merge with eGFR values
egfr_values = xenium.obs[[
    "TMA_ID", "Factor2_ulm_mean", 'disease_id']].drop_duplicates()
ct_props_long = ct_props_long.merge(
    egfr_values,
    on="TMA_ID",
    how="left"
)

# Ensure the column exists and drop rows with missing eGFR values
ct_props_long = ct_props_long.dropna(
    subset=["Factor2_ulm_mean"])
ct_props_long["Factor2_ulm_mean"] = ct_props_long["Factor2_ulm_mean"].astype(float)

# Plot scatterplot of eGFR vs cell type abundance, one figure per cell type and save
for cell_type in cell_types:

    cell_type_filesafe = cell_type.replace('/', '_').replace(" ", "_")

    data = ct_props_long[ct_props_long["cell_type_manual"] == cell_type]
    # Double-check the column exists in the filtered data
    plt.figure(figsize=(6, 4.5))
    sns.scatterplot(data=data, x="Factor2_ulm_mean",
                    y="proportion", alpha=0.7, hue="disease_id")
    # Remove legend title and set legend position
    plt.legend(title=None, loc="upper right",
               bbox_to_anchor=(1, 1),
               borderaxespad=0, frameon=True)
    # Drop NA values for regression
    reg_data = data.dropna(subset=["Factor2_ulm_mean", "proportion"])
    X = reg_data["Factor2_ulm_mean"]
    y = reg_data["proportion"]
    X_const = sm.add_constant(X)
    model = sm.OLS(y, X_const).fit()
    slope = model.params["Factor2_ulm_mean"]
    pval = model.pvalues["Factor2_ulm_mean"]
    corr = X_const["Factor2_ulm_mean"].corr(y)
    # Plot regression line with confidence interval
    x_pred = np.linspace(X.min(), X.max(), 100)
    X_pred = sm.add_constant(x_pred)
    y_pred = model.predict(X_pred)
    pred_int = model.get_prediction(X_pred).conf_int()
    plt.plot(x_pred, y_pred, color="black", label="Linear fit")
    plt.fill_between(
        x_pred, pred_int[:, 0], pred_int[:, 1], color="gray", alpha=0.2, label="95% CI")
    # Annotate slope and p-value
    # Choose text position based on slope sign
    if abs(slope) < 0.001:
        xpos, ypos, align, vertalign = 0.05, 0.95, "top", "left"
    elif slope < 0:
        xpos, ypos, align, vertalign = 0.05, 0.015, "bottom", "left"
    else:
        xpos, ypos, align, vertalign = 0.95, 0.015, "bottom", "right"
    plt.text(
        0.01, 0.99,
        f"Correlation: {corr:.2f}\nSlope: {slope:.3f}\nP-value: {pval:.2e}",
        transform=plt.gca().transAxes,
        verticalalignment='top',
        fontsize=10,
        horizontalalignment="left",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0)
    )
    plt.xlabel("Factor2")
    plt.ylabel(f"Abundance of {cell_type} cells")
    plt.title(f"Factor2 vs Abundance of {cell_type} cells")
    plt.tight_layout()
    plt.savefig(
        f"{res_dir}abundance_vs_factor2/Factor2_vs_{cell_type_filesafe}_abundance.pdf")


# --------------------------------------------------------------------------------------------
# Scatterplot of eGFR vs cell type abundance -------------------------------------------------
# --------------------------------------------------------------------------------------------

# For each TMA_ID, compute the proportion of each cell type
cell_type_proportions = xenium.obs.groupby(["TMA_ID", "cell_type_manual"]).size().unstack(fill_value=0)
cell_type_proportions = cell_type_proportions.div(cell_type_proportions.sum(axis=1), axis=0)

# Pivot longer into a dataframe with TMA_ID, cell_type, proportion
ct_props_long = cell_type_proportions.reset_index().melt(
    id_vars=["TMA_ID"],
    var_name="cell_type_manual",
    value_name="proportion"
)

# Merge with eGFR values
egfr_values = xenium.obs[["TMA_ID", "egfr_at_time_of_biopsy", 'disease_id']].drop_duplicates()
ct_props_long = ct_props_long.merge(
    egfr_values,
    on="TMA_ID",
    how="left"
)

# Ensure the column exists and drop rows with missing eGFR values
ct_props_long = ct_props_long.dropna(
    subset=["egfr_at_time_of_biopsy"])
ct_props_long["egfr_at_time_of_biopsy"] = ct_props_long["egfr_at_time_of_biopsy"].astype(float)

# Plot scatterplot of eGFR vs cell type abundance, one figure per cell type and save
for cell_type in cell_types:

    cell_type_filesafe = cell_type.replace('/', '_').replace(" ", "_")

    data = ct_props_long[ct_props_long["cell_type_manual"] == cell_type]

    # Plot figure
    plt.figure(figsize=(6, 4))
    sns.scatterplot(data=data, x="egfr_at_time_of_biopsy",
                    y="proportion", alpha=0.7, hue="disease_id")
    # Remove legend title and set legend position
    plt.legend(title=None, loc="upper right",
               bbox_to_anchor=(1, 1),
               borderaxespad=0, frameon=True)
    # Drop NA values for regression
    reg_data = data.dropna(subset=["egfr_at_time_of_biopsy", "proportion"])
    X = reg_data["egfr_at_time_of_biopsy"]
    y = reg_data["proportion"]
    X_const = sm.add_constant(X)
    model = sm.OLS(y, X_const).fit()
    slope = model.params["egfr_at_time_of_biopsy"]
    pval = model.pvalues["egfr_at_time_of_biopsy"]
    corr = X_const["egfr_at_time_of_biopsy"].corr(y)
    # Plot regression line with confidence interval
    x_pred = np.linspace(X.min(), X.max(), 100)
    X_pred = sm.add_constant(x_pred)
    y_pred = model.predict(X_pred)
    pred_int = model.get_prediction(X_pred).conf_int()
    plt.plot(x_pred, y_pred, color="black", label="Linear fit")
    plt.fill_between(x_pred, pred_int[:, 0], pred_int[:, 1], color="gray", alpha=0.2, label="95% CI")
    # Annotate slope and p-value
    plt.text(
        0.01, 0.99,
        f"Correlation: {corr:.2f}\nSlope: {slope:.3f}\nP-value: {pval:.2e}",
        fontsize=10,
        transform=plt.gca().transAxes,
        horizontalalignment="left",
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0)
    )
    plt.xlabel("eGFR at time of biopsy")
    plt.ylabel(f"Abundance of {cell_type} cells")
    plt.title(f"Scatterplot of eGFR vs {cell_type} Cells Abundance")
    plt.tight_layout()
    plt.savefig(f"{res_dir}abundance_vs_egfr/egfr_vs_{cell_type_filesafe}_abundance.pdf")


