#
# Created on Mon Jun 30 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script generates scatterplots and spatial figure grids to visualize the relationship
# between clinical variables and LIANA factors in Xenium kidney biopsy spatial transcriptomics data.
#
# Datasets used/loaded:
# - Xenium processed data: xenium_processed.h5ad
#
# Main analyses and plots:
# 1. Computes mean LIANA Factor2_ulm scores for each TMA_ID and cell type (especially aPT/PEC).
# 2. Plots a scatterplot of eGFR at time of biopsy (x-axis) vs. mean Factor2_ulm in aPT/PEC cells (y-axis),
#    colored by disease_id, with a linear regression fit and confidence interval.
#    - Output: xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.pdf
#             xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.png
# 3. Overlays spatial mini-plots (Scanpy spatial plots of Factor2_ulm) at each sample's scatterplot position,
#    with subplot borders colored by disease_id.
#    - Output: xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.pdf
#             xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.png
# 4. Repeats the above with a LOWESS-smoothed fit instead of linear regression.
#    - Output: xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_smooth.pdf
#             xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_smooth.png
# 5. Generates a stacked histogram of Factor2_ulm values by cell type.
#
# Summary:
# - The script visualizes the association between clinical eGFR and LIANA Factor2_ulm (aPT/PEC),
#   stratified by disease, and spatially contextualizes each sample with its spatial transcriptomics map.
# - All output plots are saved in the xenium_scatterplot_of_figures/ directory.
# --------------------------------------------------------------------------------------------

from itertools import combinations
import seaborn as sns
from matplotlib import patches as mpatches
import scanpy as sc
import numpy as np
from matplotlib.lines import Line2D
import statsmodels.api as sm
import matplotlib.pyplot as plt

xenium_data_dir = "/data/Kramann_lab/From_Sidrah/new_data/h5ad_files_only_transcript/marker_genes/"
data_dir = "/data/zd775156/kidney_biopsy/data/"
data_dir = "/data/Kramann_lab/kidney_biopsy/data/processed/xenium/"
res_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/"


# --------------------------------------------------------------------------------------------
# Utility functions --------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------


def resolve_overlaps(layout_df, width=0.15, height=0.02, max_iter=100, padding=0.005):
    coords = layout_df[["egfr_norm", "factor2_norm"]].values
    moved = coords.copy()

    for _ in range(max_iter):
        overlap_found = False
        for i, j in combinations(range(len(moved)), 2):
            x1, y1 = moved[i]
            x2, y2 = moved[j]

            # Check overlap
            dx = abs(x1 - x2)
            dy = abs(y1 - y2)
            if dx < width + padding and dy < height + padding:
                # Overlap detected
                overlap_found = True

                # Vector from j to i (repel away from each other)
                direction = np.array([x1 - x2, y1 - y2])
                if np.linalg.norm(direction) == 0:
                    # Completely overlapping; apply arbitrary direction
                    direction = np.random.rand(2)
                direction = direction / np.linalg.norm(direction)

                # Calculate minimal movement needed
                move_dist_x = (width + padding - dx) / 2
                move_dist_y = (height + padding - dy) / 2
                move_vector = direction * max(move_dist_x, move_dist_y)

                # Move i and j in opposite directions
                moved[i] += move_vector
                moved[j] -= move_vector
        if not overlap_found:
            break

    layout_df["egfr_resolved"] = moved[:, 0]
    layout_df["factor2_resolved"] = moved[:, 1]
    return layout_df


# --------------------------------------------------------------------------------------------
# Load Xenium data ---------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

xenium = sc.read_h5ad(data_dir + "xenium_processed.h5ad")

xenium.obs[xenium.obs["egfr_at_time_of_biopsy"].isna()][["TMA_ID", "egfr_at_time_of_biopsy"]].drop_duplicates()
# 378	Unassigned	NaN
# 60298	Utrecht2 E5 T08-05972 K121 FSGS	NaN
# 60350	Utrecht2 E6 T06-13851 K110 FSGS	NaN
# 138530	Utrecht1 E1 T13-00212 K111 FSGS	NaN


cell_types = xenium.obs["cell_type_manual"].cat.categories.tolist()

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
# Get the number of TMA_IDs for each disease_id
disease_counts = xenium.obs[['disease_id', "TMA_ID"]].drop_duplicates().groupby('disease_id')['TMA_ID'].nunique()
disease_counts = disease_counts.reset_index()
print("Number of TMA_IDs for each disease_id:")
print(disease_counts)

# Create a subplot grid with 5 columns and 6 rows
disease_ids = disease_counts["disease_id"].tolist()

# --------------------------------------------------------------------------------------------
# Plot clinical vars, quantitative ------------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Define clinical variables to plot ---------------------------------------------------
clinical_vars = ["egfr_at_time_of_biopsy", "age_at_biopsy",
                 "proteinuria_around_biopsy protein to creatine ratio"]


# --------------------------------------------------------------------------------------------
# Plot clinical vars, quantitative, not ordered by disease -----------------------------------
# --------------------------------------------------------------------------------------------

xenium_obs_ordered = xenium.obs[["disease_id", "TMA_ID", "egfr_at_time_of_biopsy"]].sort_values(
    by=["egfr_at_time_of_biopsy"],
    ascending=[False]
    ).drop_duplicates()

# Get a list of TMA_IDs ordered by clinical_var only (no disease_id sorting)
tma_ids_ordered_egfr = [xenium_obs_ordered["TMA_ID"].tolist()]

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

# --------------------------------------------------------------------------------------------
# Plot LIANA factors at cell level ---------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Compute mean Factor2_ulm per TMA_ID and add as a new column to xenium.obs
factor2_ulm_mean_per_tma = xenium.obs.groupby("TMA_ID")["Factor2_ulm"].mean()
xenium.obs["Factor2_ulm_mean"] = xenium.obs["TMA_ID"].map(factor2_ulm_mean_per_tma)


# Compute mean Factor2_ulm per TMA_ID and add as a new column to xenium.obs
factor2_ulm_mean_per_tma = xenium.obs.groupby("TMA_ID")["Factor2_ulm"].mean()
xenium.obs["Factor2_ulm_mean"] = xenium.obs["TMA_ID"].map(
    factor2_ulm_mean_per_tma)

# For each TMA_ID and for each cell type, compute the mean Factor2_ulm for that cell type
cell_type_factors = xenium.obs.groupby(["TMA_ID", "cell_type_manual"])[
    "Factor2_ulm"].mean().unstack()
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
# Rename all columns except cell_type_manual and TMA_ID by appending "Factor2_ulm_mean_"
cell_type_factors_wide.columns = [
    col if col in ["TMA_ID", "cell_type_manual"] else f"Factor2_ulm_mean_{col}"
    for col in cell_type_factors_wide.columns
]

# Left merge with xenium.obs on TMA_ID
xenium.obs = xenium.obs.merge(
    cell_type_factors_wide,
    on="TMA_ID",
    how="left"
)

# Check that every TMA_ID has one unique value of Factor2_ulm_mean_aPT/PEC
assert xenium.obs.groupby("TMA_ID")["Factor2_ulm_mean_aPT/PEC"].nunique().max() == 1, \
    "There are multiple unique values of Factor2_ulm_mean_aPT/PEC"


# --------------------------------------------------------------------------------------------
# Plot LIANA factors ULM with scatterplot of figures -----------------------------------------
# --------------------------------------------------------------------------------------------

factor_var = "Factor2_ulm_mean_aPT/PEC"

# DataFrame with layout positions and identifiers
layout_df = xenium.obs[["TMA_ID", factor_var,
                        "egfr_at_time_of_biopsy", "disease_id"]].drop_duplicates()
# Remove rows where egfr_at_time_of_biopsy or factor_var is NaN
layout_df = layout_df.dropna(subset=["egfr_at_time_of_biopsy", factor_var])
# Drop "CTRL" from the categories of disease_id
layout_df["disease_id"] = layout_df["disease_id"].cat.remove_unused_categories()
# Rename egfr_at_time_of_biopsy to x and factor_var to y
layout_df.rename(columns={
    factor_var: "factor2",
    "egfr_at_time_of_biopsy": "egfr"
}, inplace=True)
layout_df["factor2"] = layout_df["factor2"].astype(float)


# Normalize for placement in figure space
layout_df["egfr_norm"] = (layout_df["egfr"] - layout_df["egfr"].min()) / \
    (layout_df["egfr"].max() - layout_df["egfr"].min())
layout_df["factor2_norm"] = (layout_df["factor2"] - layout_df["factor2"].min()) / \
    (layout_df["factor2"].max() - layout_df["factor2"].min())

# Stacked histogram of Factor2_ulm, colored by disease_id
xenium.obs["Factor2_ulm"].hist(
    bins=50,
    by=xenium.obs["cell_type_manual"],
    stacked=True,
    # color=[disease_id_colormap[d] for d in xenium.obs["disease_id"].cat.categories],
    alpha=0.7,
    edgecolor='black',
    figsize=(10, 6)
)

# Fit linear model
X = layout_df["egfr"]
y = layout_df["factor2"]
X_with_const = sm.add_constant(X)
model = sm.OLS(y, X_with_const).fit()
# Get predictions and confidence intervals
pred = model.get_prediction(X_with_const)
pred_mean = pred.predicted_mean
ci_lower, ci_upper = pred.conf_int().T
# Plot regression line
x_sorted = np.argsort(X)

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=layout_df,
    x="egfr",
    y="factor2",
    hue="disease_id",
    alpha=0.7,
    edgecolor=None
)
plt.xlabel("eGFR at time of biopsy")
plt.ylabel("Factor2 from aPT/PEC cells")
plt.plot(X.iloc[x_sorted], pred_mean[x_sorted], color='black', label='Fit')
plt.fill_between(X.iloc[x_sorted], ci_lower[x_sorted], ci_upper[x_sorted], color='gray', alpha=0.3, label='95% CI')
# Annotate slope and p-value
slope = model.params[1]
pval = model.pvalues[1]
cor = np.corrcoef(X, y)[0, 1]
plt.text(
    0.99, 0.01,
    f"Correlation: {cor:.2f}\nSlope: {slope:.2f}\nP-value: {pval:.2e}",
    transform=plt.gca().transAxes,
    fontsize=10,
    horizontalalignment='right',
    verticalalignment='bottom',
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)
plt.legend(
        title="",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=True
    )
plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.pdf",
    bbox_inches='tight', dpi=150)


# ----------------------------------------------------------------------------------
# Real subplot positionning within the x-y grid ------------------------------------
# ----------------------------------------------------------------------------------


vmin = np.nanmin(xenium.obs["Factor2_ulm"].dropna().values)
vmax = np.nanmax(xenium.obs["Factor2_ulm"].dropna().values)

xrange = layout_df["egfr"].min(), layout_df["egfr"].max()
yrange = layout_df["factor2"].min(), layout_df["factor2"].max()

fig = plt.figure(figsize=(20, 20))

# Step 1: Create the main axis with x and y axes and ticks
main_ax = fig.add_subplot(111)  # This uses data coordinates
buffer_x = 0.05 * (xrange[1] - xrange[0])
buffer_y = 0.05 * (yrange[1] - yrange[0])
main_ax.set_xlim(xrange[0] - buffer_x, xrange[1] + buffer_x)
main_ax.set_ylim(yrange[0] - buffer_y, yrange[1] + buffer_y)

main_ax.set_xlabel("eGFR", fontsize=25)
main_ax.set_ylabel("Factor2 in aPT/PEC cells", fontsize=25)
main_ax.tick_params(labelsize=20)

# Step 2: add the fitted regression
plt.plot(X.iloc[x_sorted], pred_mean[x_sorted], color='black', label='Fit')
plt.fill_between(X.iloc[x_sorted], ci_lower[x_sorted],
                 ci_upper[x_sorted], color='gray', alpha=0.3, label='95% CI')
# Annotate slope, p-value, and correlation
slope = model.params[1]
pval = model.pvalues[1]
cor = np.corrcoef(X, y)[0, 1]
plt.text(
    0.97, 0.03,
    f"Correlation: {cor:.2f}\nSlope: {slope:.2f}\nP-value: {pval:.2e}",
    transform=plt.gca().transAxes,
    fontsize=18,
    horizontalalignment='right',
    verticalalignment='bottom',
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)

# Optional grid
main_ax.grid(True, linestyle="--", alpha=0.5)

# Step 3: Set subplot size (in figure coords)
subplot_width = subplot_height = 0.05

# Step 4: Loop through layout_df and place axes based on data coords
for i, row in layout_df.iterrows():
    tma_id = row["TMA_ID"]
    x_data = row["egfr"]
    y_data = row["factor2"]
    disease_id = row["disease_id"]

    # Convert data coords to figure coords
    display_coord = main_ax.transData.transform((x_data, y_data))
    figure_coord = fig.transFigure.inverted().transform(display_coord)
    x_fig, y_fig = figure_coord

    # Place new axes using figure coordinates
    ax = fig.add_axes([x_fig - subplot_width/2, y_fig -
                      subplot_height/2, subplot_width, subplot_height])

    # Plot Scanpy spatial data
    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
    sc.pl.spatial(
        slide_data,
        color="Factor2_ulm",
        spot_size=12,
        img=None,
        library_id=None,
        # img_key="hires",
        title="",
        show=False,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        colorbar_loc=None
    )
    # Custom styling
    ax.set_title("")
    # remove x and y labels
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(disease_id_colormap.get(disease_id, "black"))
        spine.set_linewidth(2)
    # Create custom legend handles for disease_id
    disease_legend_handles = [
        Line2D([0], [0], marker='o', color='w', label=disease,
               markerfacecolor=color, markersize=15)
        for disease, color in disease_id_colormap.items()
    ]
    # Place the legend at the right of the figure
    fig.legend(
        handles=disease_legend_handles,
        title="Disease ID",
        loc="center left",
        bbox_to_anchor=(0.8, 0.8),
        fontsize=18,
        title_fontsize=20,
        frameon=False
    )

plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.pdf",
    bbox_inches='tight', dpi=150)
plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr.png",
    bbox_inches='tight', dpi=150)


# ----------------------------------------------------------------------------------
# Same, with average Factor2 in all cells ------------------------------------------
# ----------------------------------------------------------------------------------

# Plot using the sample-level mean Factor2 (across all cells) on the y axis

# Prepare layout with mean Factor2 per TMA
layout_all = xenium.obs[["TMA_ID", "Factor2_ulm_mean", "egfr_at_time_of_biopsy", "disease_id"]].drop_duplicates()
layout_all = layout_all.dropna(subset=["egfr_at_time_of_biopsy", "Factor2_ulm_mean"])
layout_all["disease_id"] = layout_all["disease_id"].cat.remove_unused_categories()
layout_all = layout_all.rename(columns={"Factor2_ulm_mean": "factor2_mean", "egfr_at_time_of_biopsy": "egfr"})
layout_all["factor2_mean"] = layout_all["factor2_mean"].astype(float)

# Fit linear model
X_all = layout_all["egfr"]
y_all = layout_all["factor2_mean"]
Xc_all = sm.add_constant(X_all)
model_all = sm.OLS(y_all, Xc_all).fit()
pred_all = model_all.get_prediction(Xc_all)
pred_mean_all = pred_all.predicted_mean
ci_all_lower, ci_all_upper = pred_all.conf_int().T
x_idx_all = np.argsort(X_all)


# Simple scatter + regression
plt.figure(figsize=(8, 6))
sns.scatterplot(data=layout_all, x="egfr", y="factor2_mean", hue="disease_id", alpha=0.7, edgecolor=None)
plt.plot(X_all.iloc[x_idx_all], pred_mean_all[x_idx_all], color="black", label="Fit")
plt.fill_between(X_all.iloc[x_idx_all], ci_all_lower[x_idx_all], ci_all_upper[x_idx_all],
                 color="gray", alpha=0.3, label="95% CI")
slope_all = model_all.params[1]
pval_all = model_all.pvalues[1]
cor_all = np.corrcoef(X_all, y_all)[0, 1]
plt.xlabel("eGFR at time of biopsy")
plt.ylabel("Mean Factor2 (all cells) per sample")
plt.text(0.99, 0.01, f"Correlation: {cor_all:.2f}\nSlope: {slope_all:.2f}\nP-value: {pval_all:.2e}",
         transform=plt.gca().transAxes, fontsize=10, horizontalalignment="right", verticalalignment="bottom",
         bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
plt.legend(title="", loc="center left", bbox_to_anchor=(1, 0.5), frameon=True)
plt.savefig(res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_mean_all_cells.pdf",
            bbox_inches="tight", dpi=150)

# Spatial mini-plots overlaid on the scatter grid (mean across all cells)
vmin_all = np.nanmin(xenium.obs["Factor2_ulm"].dropna().values)
vmax_all = np.nanmax(xenium.obs["Factor2_ulm"].dropna().values)

xrange_all = (layout_all["egfr"].min(), layout_all["egfr"].max())
yrange_all = (layout_all["factor2_mean"].min(), layout_all["factor2_mean"].max())


fig_all = plt.figure(figsize=(20, 18))
main_ax_all = fig_all.add_subplot(111)
buf_x = 0.05 * (xrange_all[1] - xrange_all[0])
buf_y = 0.05 * (yrange_all[1] - yrange_all[0])
main_ax_all.set_xlim(xrange_all[0] - buf_x, xrange_all[1] + buf_x)
main_ax_all.set_ylim(yrange_all[0] - buf_y, yrange_all[1] + buf_y)
main_ax_all.set_xlabel("eGFR", fontsize=25)
main_ax_all.set_ylabel("Mean Factor2 (all cells) per sample", fontsize=25)
main_ax_all.tick_params(labelsize=20)

# Plot regression and CI on main axis
plt.plot(X_all.iloc[x_idx_all], pred_mean_all[x_idx_all], color="black", label="Fit")
plt.fill_between(X_all.iloc[x_idx_all], ci_all_lower[x_idx_all], ci_all_upper[x_idx_all],
                 color="gray", alpha=0.3, label="95% CI")
plt.text(0.97, 0.03, f"Correlation: {cor_all:.2f}\nSlope: {slope_all:.2f}\nP-value: {pval_all:.2e}",
         transform=plt.gca().transAxes, fontsize=18, horizontalalignment="right", verticalalignment="bottom",
         bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
main_ax_all.grid(True, linestyle="--", alpha=0.5)

# Step 3: Set subplot size (in figure coords) for mean-all-cells figure
subplot_width = subplot_height = 0.05

# Step 4: Loop through layout_all and place axes based on data coords
for i, row in layout_all.iterrows():
    tma_id = row["TMA_ID"]
    x_data = row["egfr"]
    y_data = row["factor2_mean"]
    disease_id = row["disease_id"]

    # Convert data coords to figure coords
    display_coord = main_ax_all.transData.transform((x_data, y_data))
    figure_coord = fig_all.transFigure.inverted().transform(display_coord)
    x_fig, y_fig = figure_coord

    # Place new axes using figure coordinates
    ax = fig_all.add_axes([x_fig - subplot_width / 2, y_fig - subplot_height / 2, subplot_width, subplot_height])

    # Plot Scanpy spatial data for this TMA
    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
    sc.pl.spatial(
        slide_data,
        color="Factor2_ulm",
        spot_size=12,
        img=None,
        library_id=None,
        title="",
        show=False,
        ax=ax,
        vmin=vmin_all,
        vmax=vmax_all,
        colorbar_loc=None
    )

    # Styling for the mini-plot
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(disease_id_colormap.get(disease_id, "black"))
        spine.set_linewidth(2)

# Add a single legend for disease colors
disease_legend_handles = [
    Line2D([0], [0], marker='o', color='w', label=disease, markerfacecolor=color, markersize=12)
    for disease, color in disease_id_colormap.items()
]
fig_all.legend(
    handles=disease_legend_handles,
    title="Disease ID",
    loc="center left",
    bbox_to_anchor=(0.8, 0.8),
    fontsize=14,
    title_fontsize=16,
    frameon=False
)

# Save the combined figure (mean across all cells) as PDF and PNG
fig_all.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_mean_all_cells.pdf",
    bbox_inches="tight",
    dpi=300
)
fig_all.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_mean_all_cells.png",
    bbox_inches="tight",
    dpi=300
)


# ----------------------------------------------------------------------------------
# Real subplot positionning within the x-y grid, smooth curve ----------------------
# ----------------------------------------------------------------------------------

fig = plt.figure(figsize=(20, 20))

# Step 1: Create the main axis with x and y axes and ticks
main_ax = fig.add_subplot(111)  # This uses data coordinates
buffer_x = 0.05 * (xrange[1] - xrange[0])
buffer_y = 0.05 * (yrange[1] - yrange[0])
main_ax.set_xlim(xrange[0] - buffer_x, xrange[1] + buffer_x)
main_ax.set_ylim(yrange[0] - buffer_y, yrange[1] + buffer_y)

main_ax.set_xlabel("eGFR", fontsize=25)
main_ax.set_ylabel("Factor2 in aPT/PEC cells", fontsize=25)
main_ax.tick_params(labelsize=20)

# Step 2: add the fitted regression
# Use LOWESS smoothing for a smooth curve fit
lowess = sm.nonparametric.lowess
smoothed = lowess(y, X, frac=0.5)
plt.plot(smoothed[:, 0], smoothed[:, 1], color='black', label='LOWESS fit',
         linewidth=2)
main_ax.grid(True, linestyle="--", alpha=0.5)

# Step 3: Set subplot size (in figure coords)
subplot_width = subplot_height = 0.05

# Step 4: Loop through layout_df and place axes based on data coords
for i, row in layout_df.iterrows():
    tma_id = row["TMA_ID"]
    x_data = row["egfr"]
    y_data = row["factor2"]
    disease_id = row["disease_id"]

    # Convert data coords to figure coords
    display_coord = main_ax.transData.transform((x_data, y_data))
    figure_coord = fig.transFigure.inverted().transform(display_coord)
    x_fig, y_fig = figure_coord

    # Place new axes using figure coordinates
    ax = fig.add_axes([x_fig - subplot_width/2, y_fig -
                      subplot_height/2, subplot_width, subplot_height])

    # Plot Scanpy spatial data
    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
    sc.pl.spatial(
        slide_data,
        color="Factor2_ulm",
        spot_size=12,
        img_key="hires",
        title="",
        show=False,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        alpha=0.2,
        colorbar_loc=None
    )
    # Custom styling
    ax.set_title("")
    # remove x and y labels
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(disease_id_colormap.get(disease_id, "black"))
        spine.set_linewidth(2)
    # Create custom legend handles for disease_id
    disease_legend_handles = [
        Line2D([0], [0], marker='o', color='w', label=disease,
               markerfacecolor=color, markersize=15)
        for disease, color in disease_id_colormap.items()
    ]
    # Place the legend at the right of the figure
    fig.legend(
        handles=disease_legend_handles,
        title="Disease ID",
        loc="center left",
        bbox_to_anchor=(0.8, 0.8),
        fontsize=18,
        title_fontsize=20,
        frameon=False
    )
plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_smooth.pdf",
    bbox_inches='tight', dpi=150)
plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_smooth.png",
    bbox_inches='tight', dpi=150)
