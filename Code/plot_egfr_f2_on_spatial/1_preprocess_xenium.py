#
# Created on Mon Jun 30 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script preprocesses spatial transcriptomics data from Xenium kidney biopsy samples.
# 
# Datasets used/loaded:
# - Xenium annotated data: clustered_annotated_final_file.h5ad
# - Clinical metadata: shared_for_xenium_analysis.xlsx and shared_for_xenium_analysis_with_kpmp_egfr.xlsx
# - LIANA factors: liana_factors.csv
# - LIANA loadings: loadings_adata_LIANA_all_all_5.csv
#
# Dataset saved:
# - Preprocessed Xenium data: xenium_processed.h5ad
#
# The script loads annotated Xenium data, cleans metadata, merges clinical and LIANA factor information,
# computes cell type-specific factor scores, filters low-quality samples, and generates summary plots.
# The processed data is saved for downstream analysis.
# --------------------------------------------------------------------------------------------


import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import decoupler as dc
import re
import seaborn as sns

xenium_data_dir = "/data/Kramann_lab/From_Sidrah/new_data/h5ad_files_only_transcript/marker_genes/"
metadata_dir = "/data/Kramann_lab/kidney_biopsy/data/metadata/"
data_dir = "/data/zd775156/kidney_biopsy/data/"
res_dir = "/data/Kramann_lab/From_Vivien/kidney_biopsy/"

# --------------------------------------------------------------------------------------------
# Load Xenium data ---------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

xenium = sc.read_h5ad(xenium_data_dir+"clustered_annotated_final_file.h5ad")

# Print the number of observations of xenium
print(xenium.n_obs)  # 200297

# Columns
print(xenium.obs.columns.tolist())

slides = xenium.obs['sample'].unique()
xenium.obs["cell_type_transfered"] = xenium.obs["predic"]
xenium.obs["disease_id"] = xenium.obs["condition"]

xenium.obs.columns

# Remove samples with NA or "Unassigned" condition --------------------------------------------------
xenium[xenium.obs['condition'].isin(['Unassigned'])].obs["TMA_ID"]
# These likely are dots without cores
xenium = xenium[~xenium.obs['condition'].isin(['Unassigned'])]

xenium[xenium.obs['condition'].isna()].obs["TMA_ID"]
xenium[xenium.obs['condition'].isna()].obs

xenium = xenium[~xenium.obs['condition'].isna()]

# Change "TMA_ID" by replacing double spaces with single spaces
xenium.obs["TMA_ID"] = xenium.obs["TMA_ID"].str.replace("  ", " ", regex=False)

cell_types = xenium.obs["cell_type_transfered"].cat.categories.tolist()

xenium.obs["cell_type_manual"].value_counts()

# Normalize and log1p transform the data ---------------------------------------------
xenium.raw = xenium  # Save the raw data
xenium.layers["counts"] = xenium.X.copy()  # Save counts in a layer
sc.pp.normalize_total(xenium)
# Log1p transform the data
sc.pp.log1p(xenium)
xenium.layers["logcounts"] = xenium.X.copy()  # Save log1p in a layer

# Rename PT-aPT cell annotations ---------------------------------------------
xenium.obs['cell_type_manual'] = xenium.obs['cell_type_manual'].astype(str)


xenium.obs['cell_type_manual'] = xenium.obs['cell_type_manual'].replace('PT/aPT', 'PT1')
xenium.obs['cell_type_manual'] = xenium.obs['cell_type_manual'].replace('PT', 'PT2')
xenium.obs["cell_type_manual"].value_counts(dropna=False)

# --------------------------------------------------------------------------------------------
# Load clinical var --------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

clinical_data = pd.read_excel(metadata_dir + "shared_for_xenium_analysis.xlsx")

clinical_data.columns
clinical_data.loc[clinical_data["egfr_at_time_of_biopsy"].isna(), :]["Xenium ID"]

clinical_data_new = pd.read_excel(
    metadata_dir + "shared_for_xenium_analysis_with_kpmp_egfr.xlsx",
    sheet_name="Sheet1")

clinical_data_new.shape

clinical_data_new["Unique TMA sample ID"].unique()
clinical_data_new["TMA Block Name"].unique()
clinical_data_new["Xenium ID"].unique()
clinical_data_kpmp_tmas = clinical_data_new["Xenium ID"].unique().tolist()

xenium_tmas = xenium.obs["TMA_ID"].unique()
xenium_appolo_tmas = [tma for tma in xenium_tmas if "Appollo" in tma]

# Intersection of TMA IDs in xenium and clinical_data_new
intersection = set(xenium_tmas).intersection(set(clinical_data_kpmp_tmas))

# Print tmas in clinical_data_new but not in xenium
tmas_not_in_xenium = set(clinical_data_kpmp_tmas) - set(xenium_tmas)
# {'Appollo d0 A5 CTRL', 'Appollo d0 A6 CTRL'}

# Compare clinical_data to clinical_data_new
aaa = clinical_data[["Xenium ID", "egfr_at_time_of_biopsy"]].set_index("Xenium ID")
bbb = clinical_data_new[["Xenium ID", "eGFR "]].set_index("Xenium ID")
bbb = bbb.rename(columns={"eGFR ": "egfr_at_time_of_biopsy"})

ccc = aaa.merge(
    bbb,
    left_index=True,
    right_index=True,
    how="outer"
)
ccc["egfr_at_time_of_biopsy"] = ccc["egfr_at_time_of_biopsy_x"].combine_first(ccc["egfr_at_time_of_biopsy_y"])
del ccc["egfr_at_time_of_biopsy_x"], ccc["egfr_at_time_of_biopsy_y"]

# Map ccc onto clinical_data
clinical_data["egfr_at_time_of_biopsy"] = clinical_data["Xenium ID"].map(ccc["egfr_at_time_of_biopsy"])

# Merge clinical data with xenium.obs on TMA_ID ------------

# Ensure that clinical does not have variables with the same name as xenium.obs
assert len(set(clinical_data.columns) & set(xenium.obs.columns)) == 0

# Ensure that all values of TMA_ID in xenium.obs are present in clinical_data["TMA Block ID"]
assert set(xenium.obs["TMA_ID"]).issubset(set(clinical_data["Xenium ID"]))

xenium.obs = xenium.obs.merge(
    clinical_data,
    left_on="TMA_ID",
    right_on="Xenium ID",
    how="left"
)

xenium.obs["TMA_ID"] = xenium.obs["TMA_ID"].astype(str)
xenium.obs["Xenium ID"] = xenium.obs["Xenium ID"].astype(str)
xenium.obs["TMA_ID"] = xenium.obs["Xenium ID"]


# Clean "proteinuria_around_biopsy protein to creatine ratio" ----------------


def clean_and_convert(val):
    if pd.isna(val):
        return np.nan
    if isinstance(val, str):
        val = val.strip().replace('\xa0', '')  # remove non-breaking spaces
        if val in ['NA', 'na', '']:
            return np.nan
        try:
            return float(val)
        except ValueError:
            return np.nan
    return val  # already numeric or other type


xenium.obs["proteinuria_around_biopsy protein to creatine ratio"] = (
    xenium.obs["proteinuria_around_biopsy protein to creatine ratio"]
    .apply(clean_and_convert)
)


# --------------------------------------------------------------------------------------------
# Load LIANA Factors -------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

liana_factors = pd.read_csv(
    metadata_dir + "factors_obs_adata_LIANA_all_all_5.csv", index_col=0)
liana_factors = liana_factors[["Factor"+str(i) for i in range(1, 6)]]
# Set index to sample_ID
liana_factors["sample_ID"] = liana_factors.index

# Get set difference between sample_IDs in xenium and liana_factors ------
aaa = xenium.obs["sample_ID"].unique()
# remove nan
aaa = aaa[~pd.isna(aaa)]
bbb = liana_factors["sample_ID"].unique()
# remove nan
bbb = bbb[~pd.isna(bbb)]
set_diff = set(aaa) - set(bbb)
list(set_diff)

# Merge liana_factors with xenium.obs on sample_ID ----------------
xenium.obs = xenium.obs.merge(
    liana_factors,
    left_on="sample_ID",
    right_on="sample_ID",
    how="left"
)

# -------------------------------------------------------------------------------------------------
# Plot LIANA factor loadings ----------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

loadings = pd.read_csv(
    metadata_dir + "loadings_adata_LIANA_all_all_5.csv", index_col=0
)

loadings_factor2 = loadings[loadings["Factor"] == "Factor2"]

loadings_cts = loadings_factor2["view"].unique()

# Split the strings in loadings_cts at the character "&"
loadings_cts_split = [
    [ct.strip() for ct in cts.split("&")] for cts in loadings_cts
]
# flatten loadings_cts_split into a list of unique cell types
loadings_cts_flat = set([ct for sublist in loadings_cts_split for ct in sublist])

xenium_ct = xenium.obs["cell_type_transfered"]

# Create a dict with keys xenium_ct and with values obtained from removing "XX_" from the values of xenium_ct
# (XX being a strig of numeric characters)
xenium_ct_dict = {ct: re.sub(r'^\d+_', '', ct) for ct in xenium_ct}

xenium.obs["view"] = xenium.obs["cell_type_transfered"].map(xenium_ct_dict)

# Print cell types present in Xenium and not in the loadings
cts_not_in_loadings = set(xenium.obs["view"]) - set(loadings_cts_flat)
# 'ATL_DTL', 'Epi_DTL', 'IC-A', 'IC-B', 'M-TAL', 'PL_B', 'SPP1+_aPT_PEC', 'muscle_mac'

# Compute the weighted sum by hand -----------------------------------------------------------


def score_factor(celltype, factor="Factor2"):
    loadings_factor = loadings[loadings["Factor"] == factor]

    if celltype != 'IC':
        xenium_subset = xenium[xenium.obs['view'] == celltype].copy()
        net = loadings_factor[loadings_factor["view"] == celltype].copy()
    elif celltype == 'IC':
        xenium_subset = xenium[xenium.obs['view'].isin(
            ["IC-A", "IC-B"])].copy()
        xenium_subset.obs['view'] = xenium_subset.obs['view'].replace(
            {"IC-A": "IC", "IC-B": "IC"})
        net = loadings_factor[loadings_factor["view"] == "IC"].copy()

    net = net[net["variable"].isin(xenium_subset.var_names)]

    # If there are empty observations, remove them ------

    X = xenium_subset.raw.X if xenium_subset.raw is not None else xenium_subset.X
    empty_mask = np.asarray((X == 0).sum(axis=1)).flatten() == X.shape[1]
    # empty_obs_names = xenium_subset.obs_names[empty_mask]
    xenium_subset = xenium_subset[~empty_mask].copy()

    res = dc.run_ulm(
        mat=xenium_subset,
        net=net,
        source='Factor',
        target='variable',
        weight='Loading',
        min_n=1,
        verbose=True,
        use_raw=False
    )
    res = xenium_subset.obsm["ulm_estimate"]

    return res


def score_factor_wsum(celltype, factor="Factor2"):
    loadings_factor = loadings[loadings["Factor"] == factor]

    if celltype != 'IC':
        xenium_subset = xenium[xenium.obs['view'] == celltype].copy()
        net = loadings_factor[loadings_factor["view"] == celltype]
    elif celltype == 'IC':
        xenium_subset = xenium[xenium.obs['view'].isin(["IC-A", "IC-B"])].copy()
        xenium_subset.obs['view'] = xenium_subset.obs['view'].replace(
            {"IC-A": "IC", "IC-B": "IC"})
        net = loadings_factor[loadings_factor["view"] == "IC"]

    net = net[net["variable"].isin(xenium_subset.var_names)]
    weight = net.set_index("variable")["Loading"]

    W = np.zeros(len(xenium_subset.var_names))
    mask = xenium_subset.var_names.isin(weight.index)
    W[mask] = weight.reindex(xenium_subset.var_names[mask]).values

    # Extract expression
    X = xenium_subset.X if isinstance(
        xenium_subset.X, np.ndarray) else xenium_subset.X.toarray()

    # Compute score
    xenium_subset.obs["scores"] = X.dot(W)

    return pd.DataFrame(xenium_subset.obs["scores"])


# Check that all values of xenium.obs["view"] are in np.unique(loadings["view"])
# assert set(xenium.obs["view"]).issubset(set(loadings["view"]))
cts_not_in_loadings = set(xenium.obs["view"]) - set(loadings["view"])

loadings_cts = [ct for ct in loadings["view"].unique() if "&" not in ct]

# Add factors scores to xenium.obs["Factor2_wsum"] -----------------------------
for factor in ["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"]:
    results = [score_factor_wsum(celltype=celltype, factor=factor).assign(view=celltype) for celltype in loadings_cts]
    results_allviews = pd.concat(results, axis=0)

    # Map scores from results_allviews into xenium.obs["Factor2_wsum"] by index alignment
    results_allviews.index = results_allviews.index.astype(int)
    xenium.obs[factor+"_wsum"] = xenium.obs.index.map(results_allviews["scores"])


# Run ULM to compute Factor scores ------------------------------------------------------


for factor in ["Factor1", "Factor2", "Factor3", "Factor4", "Factor5"]:
    results = [score_factor(celltype=celltype, factor=factor).assign(
        view=celltype) for celltype in loadings_cts]
    results_allviews = pd.concat(results, axis=0)

    # Map scores from results_allviews into xenium.obs["Factor2_wsum"] by index alignment
    results_allviews.index = results_allviews.index.astype(int)
    xenium.obs[factor+"_ulm"] = xenium.obs.index.map(results_allviews[factor])


# Compute mean Factor2_ulm per TMA_ID and add as a new column to xenium.obs
factor2_ulm_mean_per_tma = xenium.obs.groupby("TMA_ID")["Factor2_ulm"].mean()
xenium.obs["Factor2_ulm_mean"] = xenium.obs["TMA_ID"].map(
    factor2_ulm_mean_per_tma)
xenium.obs["Factor2_ulm_mean"] = xenium.obs["Factor2_ulm_mean"].astype(float)

# --------------------------------------------------------------------------------------------
# Preprocess spatial data --------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Check for TMA_IDs that appear in more than one sample --------------------------------------------
tma_sample_df = xenium.obs[["TMA_ID", "sample"]].drop_duplicates()
tma_counts = tma_sample_df.groupby("TMA_ID")["sample"].nunique()
shared_tmas = tma_counts[tma_counts > 1]
if not shared_tmas.empty:
    print("Warning: The following TMA_IDs are shared between multiple samples:")
    print(shared_tmas)
else:
    print("All TMA_IDs are unique to a single sample.")

# For each value of TMA_ID, subset xenium data and get the dimensions of the spatial data
tma_id_spat_dims = {}
for tma_id in xenium.obs["TMA_ID"].unique():
    tma_data = xenium[xenium.obs["TMA_ID"] == tma_id]
    print(f"TMA_ID: {tma_id}, Dimensions: {tma_data.shape}")
    tma_id_spat_dims[tma_id] = tma_data.shape

xenium.obs["cell_type_manual"].value_counts(dropna=False)

# Filter TMA_IDs too few cells ---------------------------------------------------
n_cell_thresh = 100

# Count the number of cells in each TMA_ID
tma_id_cell_counts = {
    tma_id: xenium[xenium.obs["TMA_ID"] == tma_id].shape[0]
    for tma_id in xenium.obs["TMA_ID"].unique()
}
tma_id_cell_counts_df = pd.DataFrame.from_dict(tma_id_cell_counts, orient='index', columns=['cell_count'])
tma_id_cell_counts_df = tma_id_cell_counts_df.sort_values(by='cell_count', ascending=True)
tma_id_cell_counts_df.head(10)

# Correct value of eGFR for TMA_ID Appollo d0 B1 CTRL --------------------

if xenium.obs[xenium.obs["TMA_ID"] ==
              "Appollo d0 B1 CTRL"]["egfr_at_time_of_biopsy"].unique() == 135:
    print("Appollo d0 B1 CTRL has egfr_at_time_of_biopsy = 135")
    print("Setting to correct value: 35")
    xenium.obs.loc[
        xenium.obs["TMA_ID"] == "Appollo d0 B1 CTRL", "egfr_at_time_of_biopsy"
    ] = 35

# Print the TMA_IDs with less than 100 cells
tiny_tma_ids = tma_id_cell_counts_df[tma_id_cell_counts_df['cell_count'] < n_cell_thresh].index.tolist()
# ['Utrecht1 D1 T15-05977 K78 FSGS',
#  'Utrecht2 B6 T18-04295 K60 IgAN',
#  'Utrecht1 B1 T18-13048 K46 MN',
#  'Appollo d0 A2 CTRL',
#  'Utrecht2 C5 T17-03068 K23 TIN',
#  'Utrecht1 D5 T13-04884 K71 FSGS',
#  'Utrecht2 C2 T18-04025 K4 DN']

# Remove TMA_IDs with too few cells from xenium data
xenium = xenium[~xenium.obs["TMA_ID"].isin(tiny_tma_ids)]

# Remove the corresponding entries from tma_id_spat_dims
tma_id_spat_dims = {tma_id: dims for tma_id, dims in tma_id_spat_dims.items() if tma_id not in tiny_tma_ids}
max_spat_dims = (
    max([dims[0] for dims in tma_id_spat_dims.values()]),
    max([dims[1] for dims in tma_id_spat_dims.values()])
)

# Plot histogram of tma_id_spat_dims ----------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
axs[0].hist([dims[0] for dims in tma_id_spat_dims.values()], bins=20, alpha=0.7)
axs[0].set_title('Histogram of TMA_ID Spatial Dimensions (width)')
axs[0].set_xlabel('Number of pixels')
axs[0].set_ylabel('Frequency')
axs[1].hist([dims[1] for dims in tma_id_spat_dims.values()], bins=20, alpha=0.7)
axs[1].set_title('Histogram of TMA_ID Spatial Dimensions (height)')
axs[1].set_xlabel('Number of pixels')
axs[1].set_ylabel('Frequency')
plt.tight_layout()

# Save unique TMA_IDs ---------------------------------------------------
tma_ids = xenium.obs["TMA_ID"].unique()
print(f"Number of unique TMA_IDs: {len(tma_ids)}")


# --------------------------------------------------------------------------------------------
# Rename manual cell types --------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

xenium.obs.columns

# Plot heatmap of crosstable between variables "cell_type_manual" and "cell_type_transfered"
crosstab = xenium.obs.groupby(
    ["cell_type_manual", "cell_type_transfered"]).size().unstack(fill_value=0)

plt.figure(figsize=(14, 10))
sns.heatmap(
    crosstab,
    annot=True,
    fmt="d",
    cmap="Blues",
    cbar=True,
    linewidths=0.5
)
plt.title("Xenium data: table of cell_type_manual vs cell_type_transfered")
plt.xlabel("cell_type_transfered")
plt.ylabel("cell_type_manual")
plt.tight_layout()
plt.savefig(res_dir + "sample_level_plots/xenium_crosstab_cell_type_manual_vs_transfered.pdf",
            bbox_inches='tight', dpi=150)


# --------------------------------------------------------------------------------------------
# Set default spatial parameters -------------------------------------------------------------
# --------------------------------------------------------------------------------------------

xenium.uns['spatial'] = {
    "library": {
        "images": {},  # omit or fill if you have an image
        "scalefactors": {
            "tissue_hires_scalef": 1.0,
            "spot_diameter_fullres": 1.0
        }
    }
}


# --------------------------------------------------------------------------------------------
# Save preprocessed xenium data --------------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Print n_obs of xenium_processed
print(xenium.n_obs)  # 199027

# Column names
print(xenium.obs.columns.tolist())

# Save xenium as h5ad file
xenium.write_h5ad(data_dir + "xenium_processed.h5ad")
