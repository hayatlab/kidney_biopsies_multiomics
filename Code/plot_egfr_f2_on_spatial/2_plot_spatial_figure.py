#
# Created on Mon Jun 30 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script generates spatial and scatter plots for Xenium kidney biopsy data, visualizing
# cell types, gene expression, clinical variables, and LIANA factors at both sample and cell level.
#
# Datasets used/loaded:
# - Xenium processed data: xenium_processed.h5ad
# - Genes of interest: POD_DEGs_controlled_for_Factor2_in_Xenium_panel.csv
#
# Output plots (by folder):
#
# xenium_cell_types/
# - all_tma_ids_grouped_by_disease.png: Spatial plots of cell types for all TMA_IDs, grouped by disease.
# - cell_types_by_dec_egfr.pdf: Spatial plots of cell types, ordered by decreasing eGFR, with subplot borders 
# colored by disease.
#
# xenium_clinical_vars/
# - {clinical_var}.pdf: Spatial plots of quantitative (e.g., egfr_at_time_of_biopsy, age_at_biopsy, proteinuria)
# and qualitative (e.g., sex, classification) clinical variables, ordered by value within each disease.
#
# xenium_liana_factors/
# - scatter_all_factors_vs_egfr_facet.pdf: Faceted scatterplots of LIANA factors (Factor1-5) vs eGFR, colored
# by disease.
# - {factor}.pdf: Spatial plots of LIANA factors (Factor1-5), ordered by decreasing eGFR per disease.
# - {factor}_{factor_type}.pdf: Spatial plots of LIANA factors at cell level (ULM and WSUM), ordered by decreasing eGFR.
#
# xenium_scatterplot_of_figures/
# - Factor2_ulm_vs_egfr_layout.pdf: Composite figure with spatial plots of Factor2_ulm for each TMA_ID, arranged
# by normalized eGFR and mean Factor2_ulm.
#
# Main visualizations:
# - Spatial plots of cell types, genes, clinical variables, and LIANA factors across TMA_IDs, grouped or ordered
# by disease or clinical values.
# - Scatterplots of LIANA factors vs eGFR, and Factor2_ulm vs eGFR, stratified by disease and cell type.
# - Composite scatterplot layout of spatial plots, positioned by clinical and factor values.
# --------------------------------------------------------------------------------------------


import seaborn as sns
from matplotlib import patches as mpatches
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

xenium_data_dir = "/data/Kramann_lab/From_Sidrah/new_data/h5ad_files_only_transcript/marker_genes/"
metadata_dir = "/data/zd775156/kidney_biopsy/metadata/"
data_dir = "/data/zd775156/kidney_biopsy/data/"
metadata_dir = "/data/zd775156/kidney_biopsy/metadata/"
res_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/"

# Load xenium data --------------------------------------------------------------------------------
xenium = sc.read_h5ad(data_dir + "xenium_processed.h5ad")
cell_types = xenium.obs["cell_type_transfered"].cat.categories.tolist()


# -------------------------------------------------------------------------------------------------
# Set plotting parameters -------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Get all unique cell types and assign them consistent colors
palette = sc.pl.palettes.default_102
colors = dict(zip(cell_types, palette[:len(cell_types)]))

# Precompute the color mapping
xenium.uns['cell_type_transfered_colors'] = [colors[ct] for ct in cell_types]

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

# -------------------------------------------------------------------------------------------------
# Plot cell types ---------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# Fill subplots using tma_id_grid: each cell contains a TMA_ID or None
fig, axs = plt.subplots(tma_id_grid.shape[0], tma_id_grid.shape[1],
                        figsize=(3 * tma_id_grid.shape[1], 3 * tma_id_grid.shape[0]))
for row in range(tma_id_grid.shape[0]):
    for col in range(tma_id_grid.shape[1]):
        ax = axs[row, col]
        tma_id = tma_id_grid[row, col]
        if tma_id is not None:
            slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
            sc.pl.spatial(
                slide_data,
                color="cell_type_transfered",
                spot_size=12,
                img_key="hires",
                title=f"{tma_id}",
                show=False,
                ax=ax,
                legend_loc=None
            )
            ax.title.set_size(8)
        else:
            ax.axis('off')
        ax.set_xlabel("")
        ax.set_ylabel("")
        if col == 0:
            # Add disease_id label to the first column of each row
            disease_id = xenium.obs.loc[xenium.obs["TMA_ID"] == tma_id_grid[row, 0], "disease_id"].iloc[0]
            ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                          fontsize=25, rotation=90, labelpad=20, va='center')
        ax.spines['top'].set_linewidth(0.2)
        ax.spines['right'].set_linewidth(0.2)
        ax.spines['left'].set_linewidth(0.2)
        ax.spines['bottom'].set_linewidth(0.2)
# Add a legend
legend_handles = [
    mpatches.Patch(color=colors[ct], label=ct) for ct in cell_types
]
legend_ax = fig.add_axes([0.85, 0.3, 0.1, 0.4])
legend_ax.legend(handles=legend_handles, loc='center', fontsize=15)
legend_ax.axis('off')
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.suptitle("Spatial plots of all TMA_IDs grouped by disease", fontsize=20, y=1.02)
plt.savefig(
    res_dir+"xenium_cell_types/all_tma_ids_grouped_by_disease.png",
    bbox_inches='tight', dpi=150
)


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

if True:
    for clinical_var in clinical_vars:

        tma_id_grid_ordered = tma_id_grid_ordered_dict_quant[clinical_var]

        # Get vmin and vmax for this variable across all TMA_IDs
        vmin = np.nanmin(xenium.obs[clinical_var].dropna().values)
        vmax = np.nanmax(xenium.obs[clinical_var].dropna().values)

        # Fill subplots using tma_id_grid: each cell contains a TMA_ID or None
        fig, axs = plt.subplots(tma_id_grid_ordered.shape[0], tma_id_grid_ordered.shape[1], figsize=(
            3 * tma_id_grid_ordered.shape[1], 3 * tma_id_grid_ordered.shape[0]))
        # fig.set_facecolor('black')  # Sets the figure background to black
        for row in range(tma_id_grid_ordered.shape[0]):
            for col in range(tma_id_grid_ordered.shape[1]):
                ax = axs[row, col]
                # ax.set_facecolor('black')  # Sets the axes background to black
                tma_id = tma_id_grid_ordered[row, col]
                if tma_id is not None:
                    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
                    sc.pl.spatial(
                        slide_data,
                        color=clinical_var,
                        spot_size=12,
                        img_key="hires",
                        title=f"{tma_id}",
                        show=False,
                        ax=ax,
                        vmin=vmin,
                        vmax=vmax,
                        colorbar_loc=None
                    )
                    ax.title.set_size(8)
                    # ax.title.set_color('white')
                else:
                    ax.axis('off')
                ax.set_xlabel("")
                ax.set_ylabel("")
                if col == 0:
                    # Add disease_id label to the first column of each row
                    disease_id = xenium.obs.loc[xenium.obs["TMA_ID"]
                                                == tma_id_grid[row, 0], "disease_id"].iloc[0]
                    ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                                  fontsize=20, rotation=90, labelpad=20, va='center', color='black')
                ax.spines['top'].set_linewidth(0.2)
                ax.spines['right'].set_linewidth(0.2)
                ax.spines['left'].set_linewidth(0.2)
                ax.spines['bottom'].set_linewidth(0.2)

        # Add a single colorbar from the first subplot (top-left)
        # Find the first non-None tma_id in the grid
        first_ax = None
        for row in range(tma_id_grid.shape[0]):
            for col in range(tma_id_grid.shape[1]):
                if tma_id_grid[row, col] is not None:
                    first_ax = axs[row, col]
                    break
            if first_ax is not None:
                break
        if first_ax is not None:
            # Get the mappable from the first axis
            # sc.pl.spatial creates a collection in ax.collections[-1]
            mappable = first_ax.collections[-1] if first_ax.collections else None
            if mappable is not None:
                cbar = fig.colorbar(mappable, ax=axs, orientation='vertical', fraction=0.02, pad=0.01)
                cbar.ax.tick_params(labelsize=12)

        plt.suptitle(
            f"{clinical_var.replace('_', ' ').capitalize()}, ordered by decreasing value",
            fontsize=20, y=1.02, color='black'
        )
        plt.tight_layout(rect=[0, 0, 0.88, 1])
        fig.savefig(
            f"{res_dir}/xenium_clinical_vars/{clinical_var}.pdf",
            bbox_inches='tight', dpi=150)


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

if True:
    for clinical_var in clinical_vars:

        # Get vmin and vmax for this variable across all TMA_IDs
        vmin = np.nanmin(xenium.obs[clinical_var].dropna().values)
        vmax = np.nanmax(xenium.obs[clinical_var].dropna().values)

        # Fill subplots using tma_id_grid: each cell contains a TMA_ID or None
        fig, axs = plt.subplots(tma_id_grid_ordered_egfr.shape[0], tma_id_grid_ordered_egfr.shape[1], figsize=(
            3 * tma_id_grid_ordered_egfr.shape[1], 3 * tma_id_grid_ordered_egfr.shape[0]))
        # fig.set_facecolor('black')  # Sets the figure background to black
        for row in range(tma_id_grid_ordered_egfr.shape[0]):
            for col in range(tma_id_grid_ordered_egfr.shape[1]):
                ax = axs[row, col]
                # ax.set_facecolor('black')  # Sets the axes background to black
                tma_id = tma_id_grid_ordered_egfr[row, col]

                ax.set_xlabel("")
                ax.set_ylabel("")
                for pos in ['top', 'right', 'left', 'bottom']:
                    ax.spines[pos].set_linewidth(4)
                    ax.spines[pos].set_alpha(0.5)

                if tma_id is not None:
                    disease_id = xenium.obs.loc[xenium.obs["TMA_ID"] == tma_id, "disease_id"].iloc[0]
                    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
                    sc.pl.spatial(
                        slide_data,
                        color=clinical_var,
                        spot_size=12,
                        img_key="hires",
                        title=f"{tma_id}",
                        show=False,
                        ax=ax,
                        vmin=vmin,
                        vmax=vmax,
                        colorbar_loc=None
                    )
                    ax.title.set_size(8)
                    # ax.title.set_color('white')
                    # Set color of ax.spines based on disease_id
                    spine_color = disease_id_colormap.get(disease_id, "black")
                    for spine in ax.spines.values():
                        spine.set_color(spine_color)
                    ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                                  fontsize=20, rotation=90, labelpad=20, va='center', color='black')
                    ax.set_xlabel("")
                else:
                    ax.axis('off')

        # Add a single colorbar from the first subplot (top-left)
        # Find the first non-None tma_id in the grid
        first_ax = None
        for row in range(tma_id_grid.shape[0]):
            for col in range(tma_id_grid.shape[1]):
                if tma_id_grid[row, col] is not None:
                    first_ax = axs[row, col]
                    break
            if first_ax is not None:
                break
        if first_ax is not None:
            # Get the mappable from the first axis
            # sc.pl.spatial creates a collection in ax.collections[-1]
            mappable = first_ax.collections[-1] if first_ax.collections else None
            if mappable is not None:
                cbar = fig.colorbar(
                    mappable, ax=axs, orientation='vertical', fraction=0.02, pad=0.01)
                cbar.ax.tick_params(labelsize=12)

                # Manually create color legend for disease_id_colormap
                legend_handles = [
                    mpatches.Patch(color=color, label=disease)
                    for disease, color in disease_id_colormap.items()
                ]
                legend_ax = fig.add_axes([0.88, 0.3, 0.1, 0.4])
                legend_ax.legend(handles=legend_handles, loc='center', fontsize=15, title="Disease")
                legend_ax.axis('off')

        plt.suptitle(
            f"{clinical_var.replace('_', ' ').capitalize()}, ordered by decreasing value",
            fontsize=20, y=1.02, color='black'
        )
        plt.tight_layout(rect=[0, 0, 0.88, 1])


# --------------------------------------------------------------------------------------------
# Plot clinical vars, qualitative ------------------------------------------------------------
# --------------------------------------------------------------------------------------------

clinical_vars_quali = ["sex", "classification"]

# Define tma_id_grid ordered by values of clinical variables ------
tma_id_grid_ordered_dict = {}

for clinical_var in clinical_vars_quali:
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

    tma_id_grid_ordered_dict[clinical_var] = tma_id_grid_ordered

if True:
    for clinical_var in clinical_vars_quali:

        tma_id_grid_ordered = tma_id_grid_ordered_dict[clinical_var]

        fig, axs = plt.subplots(tma_id_grid_ordered.shape[0], tma_id_grid_ordered.shape[1], figsize=(
            3 * tma_id_grid_ordered.shape[1], 3 * tma_id_grid_ordered.shape[0]))
        for row in range(tma_id_grid_ordered.shape[0]):
            for col in range(tma_id_grid_ordered.shape[1]):
                ax = axs[row, col]
                # ax.set_facecolor('black')  # Sets the axes background to black
                tma_id = tma_id_grid_ordered[row, col]
                if tma_id is not None:

                    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]

                    vals = slide_data.obs[clinical_var]
                    if vals.notna().sum() > 0:
                        # Normal plot
                        sc.pl.spatial(
                            slide_data,
                            color=clinical_var,
                            spot_size=12,
                            img_key="hires",
                            title=f"{tma_id}",
                            show=False,
                            ax=ax
                        )
                    else:

                        # Temporary column with constant value
                        slide_data.obs["_all_grey"] = 1
                        sc.pl.spatial(
                            slide_data,
                            color="_all_grey",
                            spot_size=12,
                            img_key="hires",
                            na_color="grey",
                            title=f"{tma_id}",
                            show=False,
                            ax=ax,
                            palette={"no_data": "grey"},
                            colorbar_loc=None
                        )
                        slide_data.obs.drop(columns=["_all_grey"], inplace=True)

                else:
                    ax.axis('off')
                ax.set_xlabel("")
                ax.set_ylabel("")
                if col == 0:
                    # Add disease_id label to the first column of each row
                    disease_id = xenium.obs.loc[xenium.obs["TMA_ID"]
                                                == tma_id_grid[row, 0], "disease_id"].iloc[0]
                    ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                                  fontsize=20, rotation=90, labelpad=20, va='center', color='black')
                ax.spines['top'].set_linewidth(0.2)
                ax.spines['right'].set_linewidth(0.2)
                ax.spines['left'].set_linewidth(0.2)
                ax.spines['bottom'].set_linewidth(0.2)

        plt.suptitle(f"{clinical_var.replace('_', ' ').capitalize()}",
                     fontsize=20, y=1.02, color='black')
        # plt.tight_layout(rect=[0, 0, 0.88, 1])
        fig.savefig(
            f"{res_dir}/xenium_clinical_vars/{clinical_var}.pdf",
            bbox_inches='tight', dpi=150)


# --------------------------------------------------------------------------------------------
# Plot LIANA factors at sample level ---------------------------------------------------------
# --------------------------------------------------------------------------------------------

factors = ["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"]


# Scatterplot Factor2 vs egfr ---------------------------------------------------------
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

g = sns.FacetGrid(
    xenium_long,
    hue="disease_id",
    col="Factor",
    col_wrap=3,
    margin_titles=True,
    height=3,
    sharex=False,
    sharey=False
)
g.map_dataframe(
    sns.scatterplot,
    x="Factor_value",
    y="egfr_at_time_of_biopsy",
    alpha=0.5
)
g.set_axis_labels("Factor value", "eGFR at time of biopsy")
g.fig.subplots_adjust(top=0.9)
g.add_legend(title="Disease")
g.fig.suptitle("Scatterplot of LIANA Factors vs eGFR at time of biopsy by disease", fontsize=16)
g.savefig(res_dir + "xenium_liana_factors/scatter_all_factors_vs_egfr_facet.pdf", bbox_inches='tight', dpi=150)

# Count the number of unique TMA_IDs which have both Factor2 and egfr_at_time_of_biopsy values
tma_ids_with_factor2_and_egfr = xenium.obs.loc[
    ~xenium.obs["Factor2"].isna() & ~xenium.obs["egfr_at_time_of_biopsy"].isna(),
    "TMA_ID"
].unique()
print(f"Number of unique TMA_IDs with both Factor2 and eGFR values: {len(tma_ids_with_factor2_and_egfr)}")

tma_ids_with_factor2 = xenium.obs.loc[
    ~xenium.obs["Factor2"].isna(),
    "TMA_ID"
].unique()
print(f"Number of unique TMA_IDs with Factor2 values: {len(tma_ids_with_factor2)}")

tma_ids_with_egfr = xenium.obs.loc[
    ~xenium.obs["egfr_at_time_of_biopsy"].isna(),
    "TMA_ID"
].unique()
print(f"Number of unique TMA_IDs with eGFR values: {len(tma_ids_with_egfr)}")


# Use the ordering from the eGFR values -------------------
tma_id_grid_ordered_factors = tma_id_grid_ordered_dict_quant["egfr_at_time_of_biopsy"]

if True:
    for factor in factors:

        # Get vmin and vmax for this variable across all TMA_IDs
        vmin = np.nanmin(xenium.obs[factor].dropna().values)
        vmax = np.nanmax(xenium.obs[factor].dropna().values)

        # Fill subplots using tma_id_grid -----------
        fig, axs = plt.subplots(tma_id_grid_ordered_factors.shape[0], tma_id_grid_ordered_factors.shape[1], figsize=(
            3 * tma_id_grid_ordered_factors.shape[1], 3 * tma_id_grid_ordered_factors.shape[0]))
        # fig.set_facecolor('black')  # Sets the figure background to black
        for row in range(tma_id_grid_ordered_factors.shape[0]):
            for col in range(tma_id_grid_ordered_factors.shape[1]):
                ax = axs[row, col]
                # ax.set_facecolor('black')  # Sets the axes background to black
                tma_id = tma_id_grid_ordered_factors[row, col]
                if tma_id is not None:
                    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
                    sc.pl.spatial(
                        slide_data,
                        color=factor,
                        spot_size=12,
                        img_key="hires",
                        title=f"{tma_id}",
                        show=False,
                        ax=ax,
                        vmin=vmin,
                        vmax=vmax,
                        colorbar_loc=None
                    )
                    ax.title.set_size(8)
                else:
                    ax.axis('off')
                ax.set_xlabel("")
                ax.set_ylabel("")
                if col == 0:
                    # Add disease_id label to the first column of each row
                    disease_id = xenium.obs.loc[xenium.obs["TMA_ID"]
                                                == tma_id_grid[row, 0], "disease_id"].iloc[0]
                    ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                                  fontsize=25, rotation=90, labelpad=20, va='center', color='black')
                ax.spines['top'].set_linewidth(0.2)
                ax.spines['right'].set_linewidth(0.2)
                ax.spines['left'].set_linewidth(0.2)
                ax.spines['bottom'].set_linewidth(0.2)

        # Add a single colorbar from the first subplot (top-left) -------------------------------
        first_ax = None
        for row in range(tma_id_grid.shape[0]):
            for col in range(tma_id_grid.shape[1]):
                if tma_id_grid[row, col] is not None:
                    first_ax = axs[row, col]
                    break
            if first_ax is not None:
                break
        if first_ax is not None:
            mappable = first_ax.collections[-1] if first_ax.collections else None
            if mappable is not None:
                cbar = fig.colorbar(
                    mappable, ax=axs, orientation='vertical', fraction=0.02, pad=0.01)
                cbar.ax.tick_params(labelsize=14)

        plt.suptitle(f"LIANA {factor}, ordered by decreasing eGFR per disease",
                     fontsize=25, y=1.02, color='black')
        plt.tight_layout(rect=[0, 0, 0.88, 1])

        fig.savefig(
            res_dir+f"xenium_liana_factors/{factor}.pdf",
            bbox_inches='tight', dpi=150)


# --------------------------------------------------------------------------------------------
# Plot cell types ordered by decreasing eGFR -------------------------------------------------
# --------------------------------------------------------------------------------------------

if True:

    # Fill subplots using tma_id_grid: each cell contains a TMA_ID or None
    fig, axs = plt.subplots(tma_id_grid_ordered_egfr.shape[0], tma_id_grid_ordered_egfr.shape[1], figsize=(
        3 * tma_id_grid_ordered_egfr.shape[1], 3 * tma_id_grid_ordered_egfr.shape[0]))

    for row in range(tma_id_grid_ordered_egfr.shape[0]):
        for col in range(tma_id_grid_ordered_egfr.shape[1]):
            ax = axs[row, col]
            # ax.set_facecolor('black')  # Sets the axes background to black
            tma_id = tma_id_grid_ordered_egfr[row, col]

            ax.set_xlabel("")
            ax.set_ylabel("")
            for pos in ['top', 'right', 'left', 'bottom']:
                ax.spines[pos].set_linewidth(5)
                ax.spines[pos].set_alpha(0.5)

            if tma_id is not None:
                disease_id = xenium.obs.loc[xenium.obs["TMA_ID"]
                                            == tma_id, "disease_id"].iloc[0]
                slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
                sc.pl.spatial(
                    slide_data,
                    color="cell_type_transfered",
                    spot_size=10,
                    img_key="hires",
                    title=f"{tma_id}",
                    show=False,
                    ax=ax,
                    vmin=vmin,
                    vmax=vmax,
                    legend_loc=None
                )
                ax.title.set_size(8)
                # ax.title.set_color('white')
                # Set color of ax.spines based on disease_id
                spine_color = disease_id_colormap.get(disease_id, "black")
                for spine in ax.spines.values():
                    spine.set_color(spine_color)
                ax.set_xlabel("")
                ax.set_ylabel("")
            else:
                ax.axis('off')

    # Add a single colorbar from the first subplot (top-left)
    first_ax = None
    for row in range(tma_id_grid.shape[0]):
        for col in range(tma_id_grid.shape[1]):
            if tma_id_grid[row, col] is not None:
                first_ax = axs[row, col]
                break
        if first_ax is not None:
            break
    if first_ax is not None:
        # Get the mappable from the first axis
        mappable = first_ax.collections[-1] if first_ax.collections else None
        if mappable is not None:
            # Manually create color legend for disease_id_colormap
            legend_handles = [
                mpatches.Patch(color=color, label=disease)
                for disease, color in disease_id_colormap.items()
            ]

            # Disease legend (lower, right)
            legend_ax = fig.add_axes([0.88, 0.3, 0.1, 0.2])
            legend_ax.legend(handles=legend_handles,
                             loc='center', fontsize=15, title="Disease")

            # Cell type legend (above disease, right)
            legend_handles_cell_types = [
                mpatches.Patch(color=colors[ct], label=ct) for ct in cell_types
            ]
            legend_ax_cell_types = fig.add_axes([0.88, 0.52, 0.1, 0.35])
            legend_ax_cell_types.legend(handles=legend_handles_cell_types, loc='center', fontsize=15, title="Cell type")
            legend_ax_cell_types.axis('off')
            legend_ax.axis('off')

    plt.suptitle(
        "Cell types ordered by decreasing eGFR value",
        fontsize=20, y=1.02, color='black'
    )
    plt.tight_layout(rect=[0, 0, 0.88, 1])

    fig.savefig(
        res_dir+"xenium_cell_types/cell_types_by_dec_egfr.pdf",
        bbox_inches='tight', dpi=150)


# --------------------------------------------------------------------------------------------
# Plot LIANA factors at cell level ---------------------------------------------------------
# --------------------------------------------------------------------------------------------


# Scatterplot of Factor2_ulm vs Factor2_ulm_mean -----

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=xenium.obs,
    x="Factor2",
    y="Factor2_ulm_mean",
    hue="disease_id",
    alpha=0.5
)
plt.gca().set_title("")  # Remove default title
plt.title("Scatterplot of Factor2 vs Mean Factor2 (ULM) by disease")
correlation = xenium.obs["Factor2"].corr(xenium.obs["Factor2_ulm_mean"])
plt.suptitle(f"Correlation: {correlation:.2f}", fontsize=12, y=1.02)

# Scatterplot of Factor2_ulm vs egfr_at_time_of_biopsy ------
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=xenium.obs,
    x="Factor2_ulm",
    y="egfr_at_time_of_biopsy",
    hue="disease_id",
    alpha=0.5
)
plt.gca().set_title("")  # Remove default title

# Facet by cell_type_transfered ------
g = sns.FacetGrid(
    xenium.obs,
    col="cell_type_transfered",
    col_wrap=4,
    hue="disease_id",
    height=3,
    sharex=True,
    sharey=True
)
g.map_dataframe(
    sns.scatterplot,
    x="Factor2_ulm",
    y="egfr_at_time_of_biopsy",
    alpha=0.5
)
g.set_axis_labels("Factor2 (ULM)", "eGFR at time of biopsy")
g.add_legend(title="Disease")
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("Scatterplot of Factor2 (ULM) vs eGFR at time of biopsy by disease and cell type")
plt.xlabel("Factor2 (ULM)")
plt.ylabel("eGFR at time of biopsy")
plt.title("Scatterplot of Factor2 (ULM) vs eGFR at time of biopsy by disease")
plt.legend(title="Disease")
plt.tight_layout()
plt.savefig(res_dir + "xenium_liana_factors/scatter_Factor2_ulm_vs_egfr.pdf", bbox_inches='tight', dpi=150)
plt.close()


if True:

    for factor_type in ["ulm", "wsum"]:  # Two methods to compute LIANA factors (ULM is the reference method)
        for factor in factors:

            # Get vmin and vmax for this variable across all TMA_IDs
            vmin = np.nanmin(xenium.obs[factor+'_'+factor_type].dropna().values)
            vmax = np.nanmax(xenium.obs[factor+'_'+factor_type].dropna().values)

            # Fill subplots using tma_id_grid: each cell contains a TMA_ID or None
            fig, axs = plt.subplots(tma_id_grid_ordered_egfr.shape[0], tma_id_grid_ordered_egfr.shape[1], figsize=(
                3 * tma_id_grid_ordered_egfr.shape[1], 3 * tma_id_grid_ordered_egfr.shape[0]))
            # fig.set_facecolor('black')  # Sets the figure background to black
            for row in range(tma_id_grid_ordered_egfr.shape[0]):
                for col in range(tma_id_grid_ordered_egfr.shape[1]):
                    ax = axs[row, col]
                    # ax.set_facecolor('black')  # Sets the axes background to black
                    tma_id = tma_id_grid_ordered_egfr[row, col]

                    ax.set_xlabel("")
                    ax.set_ylabel("")
                    for pos in ['top', 'right', 'left', 'bottom']:
                        ax.spines[pos].set_linewidth(5)
                        ax.spines[pos].set_alpha(0.5)

                    if tma_id is not None:
                        disease_id = xenium.obs.loc[xenium.obs["TMA_ID"]
                                                    == tma_id, "disease_id"].iloc[0]
                        slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]
                        sc.pl.spatial(
                            slide_data,
                            color=factor+"_"+factor_type,
                            spot_size=12,
                            img_key="hires",
                            title=f"{tma_id}",
                            show=False,
                            ax=ax,
                            vmin=vmin,
                            vmax=vmax,
                            colorbar_loc=None
                        )
                        ax.title.set_size(8)
                        # ax.title.set_color('white')
                        # Set color of ax.spines based on disease_id
                        spine_color = disease_id_colormap.get(disease_id, "black")
                        for spine in ax.spines.values():
                            spine.set_color(spine_color)
                        # ax.set_ylabel(slide_data.obs["disease_id"].iloc[0],
                        #               fontsize=25, rotation=90, labelpad=20, va='center', color='black')
                        ax.set_xlabel("")
                    else:
                        ax.axis('off')

            # Add a single colorbar from the first subplot (top-left)
            # Find the first non-None tma_id in the grid
            first_ax = None
            for row in range(tma_id_grid.shape[0]):
                for col in range(tma_id_grid.shape[1]):
                    if tma_id_grid[row, col] is not None:
                        first_ax = axs[row, col]
                        break
                if first_ax is not None:
                    break
            if first_ax is not None:
                # Get the mappable from the first axis
                # sc.pl.spatial creates a collection in ax.collections[-1]
                mappable = first_ax.collections[-1] if first_ax.collections else None
                if mappable is not None:
                    cbar = fig.colorbar(
                        mappable, ax=axs, orientation='vertical', fraction=0.02, pad=0.01)
                    cbar.ax.tick_params(labelsize=12)

                    # Manually create color legend for disease_id_colormap
                    legend_handles = [
                        mpatches.Patch(color=color, label=disease)
                        for disease, color in disease_id_colormap.items()
                    ]
                    legend_ax = fig.add_axes([0.88, 0.3, 0.1, 0.4])
                    legend_ax.legend(handles=legend_handles,
                                     loc='center', fontsize=15, title="Disease")
                    legend_ax.axis('off')

            plt.suptitle(
                f"LIANA {factor}, ordered by decreasing eGFR value",
                fontsize=20, y=1.02, color='black'
            )
            plt.tight_layout(rect=[0, 0, 0.88, 1])

            fig.savefig(
                res_dir+f"xenium_liana_factors/{factor}_{factor_type}.pdf",
                bbox_inches='tight', dpi=150)


# --------------------------------------------------------------------------------------------
# Plot LIANA factors ULM with scatterplot of figures -----------------------------------------
# --------------------------------------------------------------------------------------------

# DataFrame with layout positions and identifiers
layout_df = xenium.obs[["TMA_ID", "Factor2_ulm_mean", "egfr_at_time_of_biopsy", "disease_id"]].drop_duplicates()
# Remove rows where egfr_at_time_of_biopsy or Factor2_ulm_mean is NaN
layout_df = layout_df.dropna(subset=["egfr_at_time_of_biopsy", "Factor2_ulm_mean"])
# Drop "CTRL" from the categories of disease_id
layout_df["disease_id"] = layout_df["disease_id"].cat.remove_unused_categories()
# Rename egfr_at_time_of_biopsy to x and Factor2_ulm_mean to y
layout_df.rename(columns={
    "Factor2_ulm_mean": "factor2",
    "egfr_at_time_of_biopsy": "egfr"
}, inplace=True)


# Normalize for placement in figure space
layout_df["egfr_norm"] = (layout_df["egfr"] - layout_df["egfr"].min()) / \
    (layout_df["egfr"].max() - layout_df["egfr"].min())
layout_df["factor2_norm"] = (layout_df["factor2"] - layout_df["factor2"].min()) / \
    (layout_df["factor2"].max() - layout_df["factor2"].min())

# Scatterplot layout x and y ------
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=layout_df,
    x="egfr_norm",
    y="factor2_norm",
    hue="disease_id",
    alpha=0.7,
    edgecolor=None
)

fig = plt.figure(figsize=(20, 20))
width = height = 0.05  # size of each subplot in figure coords
for i, row in layout_df.iterrows():
    tma_id = row["TMA_ID"]
    x = row["egfr_norm"]
    y = row["factor2_norm"]  # invert y for display top-to-bottom

    vmin = np.nanmin(xenium.obs["Factor2_ulm"].dropna().values)
    vmax = np.nanmax(xenium.obs["Factor2_ulm"].dropna().values)

    # Define a new axes at this location
    ax = fig.add_axes([x - width/2, y - height/2, width, height])

    # Filter data for this TMA_ID
    slide_data = xenium[xenium.obs["TMA_ID"] == tma_id]

    # Plot into this axis
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
        colorbar_loc=None  # turn off colorbar per subplot
    )
plt.savefig(
    res_dir + "xenium_scatterplot_of_figures/Factor2_ulm_vs_egfr_layout.pdf",
    bbox_inches='tight', dpi=150)

