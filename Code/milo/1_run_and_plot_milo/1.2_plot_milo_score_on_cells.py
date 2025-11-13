#
# Created on Wed May 07 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script analyzes MILO neighborhood enrichment in kidney biopsy scRNA-seq data,
# mapping neighborhood-level scores to cells, computing significant cells, and identifying
# DEGs per cell type and subcluster.
#
# Input datasets/objects:
# - AnnData objects with clinical metadata:
#     - scVI_kidney_biopsy_with_clinical_no_TIN.h5ad
#     - scVI_kidney_biopsy_with_clinical.h5ad
# - Milo MultiData object:
#     - 01_no_cov/milo_obj.h5mu
#
# Data processing steps:
# - Remove MALAT1 gene
# - Standardize 'classification' values
# - Compute cell-level MILO scores (logFC and -log10(adj_pval)) by mapping neighborhood scores
# - Define main_subcluster and mapping dictionaries for subcluster <-> main_class relationships
#
# Visualizations generated:
# - Histograms:
#     - MILO logFC and -log10(adj_pval) at neighborhood and cell level
# - Scatterplots:
#     - logFC vs -log10(adj_pval) for neighborhoods and cells
# - UMAPs:
#     - Colored by milo_logFC, milo_adj_log10pval, and significance
# - Violin plots:
#     - MILO scores per main_class and main_subcluster
#     - MILO scores per disease_id stratified by main_class or subcluster
# - Barplots:
#     - Proportion and counts of significant MILO+ cells per main_class and subcluster
#     - Stacked barplots per logFC sign (up/down) per main_class or subcluster
# - Dotplots:
#     - Top DEGs per main_class and main_subcluster
# - Gene specificity analyses:
#     - Compute subcluster-specific logFC per gene
#     - Scatterplots of gene specificity vs -log10(adj_pval)
#
# DEGs and statistical outputs:
# - significant_cell_types.csv / significant_cell_types_subcluster.csv
# - milo_signif_DEGs.csv
# - milo_signif_DEGs_subcluster.csv
# - milo_signif_DEGs_subcluster_with_gene_specif.csv
# - top_100_genes_df.csv, top_100_genes_specif_df.csv, top_100_genes_specif_higher_7_df.csv
# - Proportion of DEGs that are significant per subcluster
# - Importance scores for DEGs (abs(logFC) * -log10(adj_pval))
#
# Key computations:
# - Mapping MILO neighborhood logFC and adjusted p-values to cells
# - Identification of significant MILO+ cells
# - DEGs computation using Wilcoxon test per main_class and per subcluster
# - Computation of gene specificity for DEGs
# - Summary statistics and proportion plots of significant cells and DEGs
#
# Output files:
# - PNG figures for all histograms, scatterplots, UMAPs, violin plots, barplots, dotplots, gene specificity plots
# - CSV files containing DEGs, top genes, significance and gene specificity scores
#
# This script provides a comprehensive pipeline to analyze MILO neighborhood enrichment
# at both cell type and subcluster resolution, including visualization and identification
# of significant genes and their specificity across subclusters.
# --------------------------------------------------------------------------------------------


from scipy.stats import ranksums
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr
from mudata import read_h5mu

proj_dir = '/home/zd775156/kidney_biopsy/'
data_dir = proj_dir + "results/milo/"
res_dir = proj_dir + "results/milo/"

# ----------------------------------------------------------------------------
# Load data  -----------------------------------------------------------------
# ----------------------------------------------------------------------------

adata = sc.read_h5ad(proj_dir+'data/scVI_kidney_biopsy_with_clinical_no_TIN.h5ad')

bdata = sc.read_h5ad(proj_dir+'data/scVI_kidney_biopsy_with_clinical.h5ad')
subcluster_df = bdata.obs["subcluster"]
del bdata

# Remove MALAT1 gene from adata
if "MALAT1" in adata.var_names:
    adata = adata[:, adata.var_names != "MALAT1"]

# Change Fast Progressor to Fast_Progressor
adata.obs['classification'] = adata.obs['classification'].str.replace(
    "Fast Progressor", "Fast_Progressor").astype("category")

mdata = read_h5mu(data_dir+"01_no_cov/milo_obj.h5mu")

nhood_logfc = mdata["milo"].var["logFC"]
nhood_pvals = -np.log10(mdata["milo"].var["SpatialFDR"])

# Compute cell-level milo scores by mapping nhood scores to cells
nhoods_matrix = mdata["rna"].obsm["nhoods"].toarray().astype("int")
print(nhoods_matrix.shape)
nhoods_matrix.sum(axis=1).shape

# Normalize nhoods_matrix so that each row sums to 1
nhoods_matrix = nhoods_matrix / nhoods_matrix.sum(axis=1, keepdims=True)
nhoods_matrix = csr_matrix(nhoods_matrix)

# Get cell pvals
adata.obs["milo_logFC"] = nhoods_matrix @ nhood_logfc
adata.obs["milo_adj_log10pval"] = nhoods_matrix @ nhood_pvals

adata.obs["logFC_sign"] = ["up" if x > 0 else "down" for x in adata.obs["milo_logFC"]]
adata.obs["milo_signif"] = adata.obs["milo_adj_log10pval"] > -np.log10(0.05)

classes_of_interest = ["11_ATL_DTL", "2_CNT_PC", "3_dPT_dTAL_dDCT", "4_aPT_aTAL_aDCT", "1_TAL"]

adata_sub = adata[adata.obs["main_class"].isin(classes_of_interest)].copy()
adata_sub.obs["subcluster"] = subcluster_df[adata_sub.obs_names].values
adata_sub.obs["subcluster"] = adata_sub.obs["subcluster"].astype("str")

# Define main_subcluster
adata_sub.obs["main_subcluster"] = adata_sub.obs["main_class"].str.split("_").str[0] + "." + adata_sub.obs["subcluster"]
adata_sub.obs["main_subcluster"].value_counts()
# Order them by alphabetical order
adata_sub.obs["main_subcluster"] = pd.Categorical(adata_sub.obs["main_subcluster"],
                                                  categories=sorted(adata_sub.obs["main_subcluster"].unique()))

# Create a dict to map subcluster names to main classes
tmp = adata_sub.obs[["subcluster", "main_class"]].drop_duplicates()
subcluster_to_main_class = tmp.set_index("subcluster")["main_class"].to_dict()

# Create a dict to map main_subcluster names to main classes
tmp = adata_sub.obs[["main_subcluster", "main_class"]].drop_duplicates()
main_subcluster_to_main_class = tmp.set_index("main_subcluster")["main_class"].to_dict()

# Create a dict to map subcluster names to main_subcluster names
tmp = adata_sub.obs[["subcluster", "main_subcluster"]].drop_duplicates()
subcluster_to_main_subcluster = tmp.set_index("subcluster")["main_subcluster"].to_dict()

# --------------------------------------------------------------------------------------
# Plot MILO scores distributions: nhood level and cell level
# --------------------------------------------------------------------------------------

# Histograms of milo pvals at neighborhood and cell level
fig, axs = plt.subplots(2, 1, figsize=(6, 8))
axs.flatten()[0].axvline(x=-np.log10(0.05), color='grey', linestyle='--')
nhood_pvals.hist(ax=axs.flatten()[0], bins=50)
axs.flatten()[1].axvline(x=-np.log10(0.05), color='grey', linestyle='--')
adata.obs["milo_adj_log10pval"].hist(ax=axs.flatten()[1], bins=50)
axs.flatten()[0].grid(False)
axs.flatten()[1].grid(False)
axs.flatten()[0].set_title("Neighborhoods level -log10(adj_pvals)")
axs.flatten()[1].set_title("Cell level -log10(adj_pvals)")
plt.tight_layout()
plt.savefig(res_dir+"milo_pval_hist.png")

# Histograms of milo logFC at neighborhood and cell level
fig, axs = plt.subplots(2, 1, figsize=(6, 8))
axs.flatten()[0].axvline(x=0, color='grey', linestyle='--')
nhood_logfc.hist(ax=axs.flatten()[0], bins=50)
axs.flatten()[1].axvline(x=0, color='grey', linestyle='--')
adata.obs["milo_logFC"].hist(ax=axs.flatten()[1], bins=50)
axs.flatten()[0].grid(False)
axs.flatten()[1].grid(False)
axs.flatten()[0].set_title("Neighborhoods level logFC")
axs.flatten()[1].set_title("Cell level logFC")
plt.tight_layout()
plt.savefig(res_dir+"milo_logFC_hist.png")

# Scatterplot of logFC vs -log10(adj_pval) at neighborhood level
mdata["milo"].var["log_adj_pval"] = -np.log10(mdata["milo"].var["SpatialFDR"])
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.scatter(mdata["milo"], x="logFC", y="log_adj_pval", ax=ax, show=False,
              color="black", size=10, alpha=1)
ax.set_xlabel("logFC")
ax.set_ylabel("-log10(adj_pval)")
ax.set_title("Milo scores at neighborhood level")
ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
ax.axvline(x=0, color='grey', linestyle='--')
plt.savefig(res_dir+"milo_logFC_vs_milo_adj_log10pval_nhood.png")

# Scatterplot of logFC vs -log10(adj_pval)
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.scatter(adata, x="milo_logFC", y="milo_adj_log10pval", ax=ax, show=False,
              color="black", size=10, alpha=1)
ax.set_xlabel("logFC")
ax.set_ylabel("-log10(adj_pval)")
ax.set_title("Milo scores at cell level")
ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
ax.axvline(x=0, color='grey', linestyle='--')
plt.savefig(res_dir+"milo_logFC_vs_milo_adj_log10pval.png")

# --------------------------------------------------------------------------------------
# Plot MILO scores on UMAP
# --------------------------------------------------------------------------------------

# UMAP colored by milo_adj_log10pval and milo_logFC
fig, axs = plt.subplots(2, 1, figsize=(8, 10))
sc.pl.umap(adata, color=["milo_adj_log10pval"], frameon=False, ncols=1, legend_loc="on data",
           size=2, ax=axs.flatten()[0], show=False)
sc.pl.umap(adata, color=["milo_logFC"], frameon=False, ncols=1, legend_loc="on data",
           add_outline=True, outline_color=('black', 'white'),
           cmap="RdBu_r", size=2, ax=axs.flatten()[1], show=False)
axs.flatten()[0].set_title("-log10(adj_pval)")
axs.flatten()[1].set_title("logFC")
plt.suptitle("Milo scores at cell level", y=0.99)
plt.tight_layout()
fig.savefig(res_dir+"umap_milo_score.png")

# UMAP colored by milo_signif and milo_logFC
fig, axs = plt.subplots(2, 1, figsize=(8, 10))
sc.pl.umap(adata, color=["milo_signif"], frameon=False, ncols=1,
           size=2, ax=axs.flatten()[0], show=False)
sc.pl.umap(adata, color=["main_class"], frameon=False, ncols=1,
           legend_loc='right margin', legend_fontsize='x-small',
           size=2, ax=axs.flatten()[1], show=False)
axs.flatten()[0].set_title("adj_pval < 0.05")
axs.flatten()[1].set_title("")
plt.suptitle("Milo: significant cells", y=1.0)
plt.tight_layout()
fig.savefig(res_dir+"umap_signif_cells.png")


# --------------------------------------------------------------------------------------
# Plot MILO scores on violin plots
# --------------------------------------------------------------------------------------

# Plot violin plots of milo_adj_log10pval and milo_logFC for each main_class
fig, axs = plt.subplots(2, 1, figsize=(10, 10))
axs.flatten()[0].axhline(y=-np.log10(0.05), color='grey', linestyle='--')
sc.pl.violin(adata, "milo_adj_log10pval", groupby="main_class", ax=axs.flatten()[0], show=False)
axs.flatten()[1].axhline(y=0, color='grey', linestyle='--')
sc.pl.violin(adata, "milo_logFC", groupby="main_class", ax=axs.flatten()[1], show=False)
axs.flatten()[0].set_xticklabels(axs.flatten()[0].get_xticklabels(), rotation=45, ha="right")
axs.flatten()[1].set_xticklabels(axs.flatten()[1].get_xticklabels(), rotation=45, ha="right")
axs.flatten()[0].set_title("-log10(adj_pval)")
axs.flatten()[1].set_title("logFC")
axs.flatten()[0].set_xlabel("")
axs.flatten()[1].set_xlabel("")
axs.flatten()[0].set_ylabel("")
axs.flatten()[1].set_ylabel("")
plt.tight_layout()
plt.suptitle("Milo scores per cell type", y=1.05)
fig.savefig(res_dir+"milo_scores_violin.png")

# Plot violin plots of milo_adj_log10pval and milo_logFC for each main_subcluster
fig, axs = plt.subplots(2, 1, figsize=(10, 10))
axs.flatten()[0].axhline(y=-np.log10(0.05), color='grey', linestyle='--')
sc.pl.violin(adata_sub, "milo_adj_log10pval", groupby="main_subcluster", ax=axs.flatten()[0], show=False)
axs.flatten()[1].axhline(y=0, color='grey', linestyle='--')
sc.pl.violin(adata_sub, "milo_logFC", groupby="main_subcluster", ax=axs.flatten()[1], show=False)
axs.flatten()[0].set_xticklabels(axs.flatten()[0].get_xticklabels(), rotation=45, ha="right")
axs.flatten()[1].set_xticklabels(axs.flatten()[1].get_xticklabels(), rotation=45, ha="right")
axs.flatten()[0].set_title("-log10(adj_pval)")
axs.flatten()[1].set_title("logFC")
axs.flatten()[0].set_xlabel("")
axs.flatten()[1].set_xlabel("")
axs.flatten()[0].set_ylabel("")
axs.flatten()[1].set_ylabel("")
plt.tight_layout()
plt.suptitle("Milo scores per cell subcluster", y=1.05)
fig.savefig(res_dir+"milo_scores_violin_subcluster.png")


# --------------------------------------------------------------------------------------
# Plot props of MILO+ cells
# --------------------------------------------------------------------------------------

# Compute props per main_class ------------------------------
props_main = adata.obs.groupby("main_class")["milo_signif"].mean()
props_main = props_main.reset_index()

# Compute props per main_subcluster ------------------------------
props_sub = adata_sub.obs.groupby("main_subcluster")["milo_signif"].mean()
props_sub = props_sub.reset_index()
props_sub = props_sub[props_sub["milo_signif"] != 0]
props_sub['main_subcluster'] = props_sub['main_subcluster'].astype("str")
props_sub['main_class'] = props_sub['main_subcluster'].map(main_subcluster_to_main_class)

# Compute counts per main_class ------------------------------
counts_main = adata.obs.groupby("main_class")["milo_signif"].sum()
counts_main = counts_main.reset_index()

# Compute counts per main_subcluster ------------------------------
counts_sub = adata_sub.obs.groupby("main_subcluster")["milo_signif"].sum()
counts_sub = counts_sub.reset_index()
counts_sub = counts_sub[counts_sub["milo_signif"] != 0]
counts_sub['main_subcluster'] = counts_sub['main_subcluster'].astype("str")
counts_sub['main_class'] = counts_sub['main_subcluster'].map(main_subcluster_to_main_class)

# Plot props per main_class ---------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
props_main.plot(kind="bar", x="main_class", y="milo_signif", ax=ax)
for i, v in enumerate(props_main["milo_signif"]):
    ax.text(i, v + 0.005, str(int(counts_main["milo_signif"][i])), ha='center', va='bottom')
ax.set_ylabel("Proportion of milo_signif cells")
ylim = ax.get_ylim()
ax.set_ylim(0, ylim[1] * 1.05)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_xlabel("")
ax.get_legend().remove()
ax.set_title("Proportion and counts of milo_signif cells per main class")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_props_main.png")

# Plot props per subcluster ---------------------------------
fig, ax = plt.subplots(figsize=(10, 5))
props_sub.plot(kind="bar", x="main_subcluster", y="milo_signif", ax=ax)
for i, v in enumerate(props_sub["milo_signif"]):
    ax.text(i, v + 0.005, str(int(counts_sub["milo_signif"].iloc[i])), ha='center', va='bottom')
ax.set_ylabel("Proportion of milo_signif cells")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ylim = ax.get_ylim()
ax.set_ylim(0, ylim[1] * 1.05)
ax.set_xlabel("")
ax.get_legend().remove()
ax.set_title("Proportion and counts of milo_signif cells per main:subcluster class")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_props_subcluster.png")


# Compute props per main_class, colored per logFC_sign -----------
up_signif = (adata.obs["logFC_sign"] == "up") & adata.obs["milo_signif"]
down_signif = (adata.obs["logFC_sign"] == "down") & adata.obs["milo_signif"]
down_signif.sum()
up_signif.sum()
# Count total cells per class
total = adata.obs["main_class"].value_counts()
# Count up and down milo_signif per class
up_counts = adata.obs.loc[up_signif, "main_class"].value_counts()
down_counts = adata.obs.loc[down_signif, "main_class"].value_counts()
# Combine into DataFrame and normalize
summary_df = pd.DataFrame({
    "logFC > 0 of fast progressors": up_counts / total,
    "logFC < 0 of fast progressors": down_counts / total
}).fillna(0).reset_index().rename(columns={"index": "main_class"})
# Stacked barplot of summary_df
fig, ax = plt.subplots(figsize=(10, 5))
summary_df.plot(kind="bar", x="main_class", stacked=True, ax=ax)
ax.set_ylabel("Proportion of milo_signif cells (adj_pval < 0.05)")
ax.set_xlabel("")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
plt.suptitle("Milo score: Proportion of milo_signif cells per main class")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_props_logFC.png")

# Compute props of milo_signif cells per subcluster and logFC_sign
up_signif = (adata_sub.obs["logFC_sign"] == "up") & adata_sub.obs["milo_signif"]
down_signif = (adata_sub.obs["logFC_sign"] == "down") & adata_sub.obs["milo_signif"]
# Count total cells per class
total = adata_sub.obs["subcluster"].value_counts()
# Count up and down milo_signif per class
up_counts = adata_sub.obs.loc[up_signif, "subcluster"].value_counts()
down_counts = adata_sub.obs.loc[down_signif, "subcluster"].value_counts()
# Combine into DataFrame and normalize
summary_df = pd.DataFrame({
    "logFC > 0 of fast progressors": up_counts / total,
    "logFC < 0 of fast progressors": down_counts / total
}).fillna(0).reset_index().rename(columns={"index": "subcluster"})
# Stacked barplot of summary_df
fig, ax = plt.subplots(figsize=(10, 5))
summary_df.plot(kind="bar", x="subcluster", stacked=True, ax=ax)
ax.set_ylabel("Proportion of milo_signif cells (adj_pval < 0.05)")
ax.set_xlabel("")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
plt.suptitle("Milo score: Proportion of milo_signif cells per subcluster")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_props_logFC_subcluster.png")


# --------------------------------------------------------------------------------------
# Compute DEGs of milo_signif cells per main_class
# --------------------------------------------------------------------------------------

# Check if there's at least one significant cell in each 'main_class'
significant_cell_types = adata.obs.groupby('main_class')['milo_signif'].apply(lambda x: x.any())
significant_cell_types.to_csv(res_dir+"significant_cell_types.csv")
# Print the cell types with at least one significant cell
print(significant_cell_types[significant_cell_types].index.tolist())

# Normalize and log-transform the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Set the 'milo_signif' column as the group for DEGs
adata.obs['milo_signif_str'] = adata.obs['milo_signif'].astype(str)

# Split the data by main_class and compute
adata_ct = [adata[adata.obs['main_class'] == cell_type].copy() for cell_type in adata.obs['main_class'].unique()]

results = {}
for ind, cell_type in enumerate(adata.obs['main_class'].unique()):
    # Skip if fewer than 2 unique values in 'milo_signif'
    if adata_ct[ind].obs['milo_signif_str'].nunique() != 2:
        print(f"Only one group in milo_signif for cell type {cell_type}, skipping...")
    elif 1 in adata_ct[ind].obs['milo_signif_str'].value_counts().values:
        print(f"Only one cell in a milo_signif group for cell type {cell_type}, skipping...")
    else:
        # Compute DEGs ---------------------
        sc.tl.rank_genes_groups(
            adata_ct[ind],
            groupby='milo_signif_str',
            method='wilcoxon',
            reference='False')
        results[cell_type] = sc.get.rank_genes_groups_df(adata_ct[ind], group='True')


# Plot DEGs ------------------------------------------------------
for cell_type in results.keys():
    df = results[cell_type]

    # Filter significant DEGs (adj p-value < 0.05)
    signif_genes = df[df["is_signif"]]["names"].iloc[:20].tolist()

    # Skip if no significant DEGs
    if len(signif_genes) == 0:
        print(f"No significant DEGs for {cell_type}")
        continue

    # Get the AnnData subset for this cell type
    # (Assuming adata_ct is a list of AnnData objects corresponding to cell_type order)
    idx = list(adata.obs['main_class'].unique()).index(cell_type)
    adata_subset = adata_ct[idx]

    # Create dotplot
    sc.pl.dotplot(
        adata_subset,
        var_names=signif_genes,
        groupby='milo_signif_str',  # The 2-group comparison
        title=f'{cell_type}: DEGs (True vs False)',
        # figsize=(4, 8 + 0.3*len(signif_genes)),  # Adjust height based on # of genes
        dendrogram=False,
        standard_scale='var',  # Scale expression per gene
        # cmap='RdBu_r',  # Diverging colormap
    )

deg_df = pd.concat(results.values(), keys=results.keys())
deg_df = deg_df.reset_index(level=0).rename(columns={"level_0": "main_class"})
# Save main_class DEGs to csv
deg_df.to_csv(res_dir+"milo_signif_DEGs.csv", index=False)

adata.obs[["main_class", "milo_signif"]].value_counts().sort_index()

print(deg_df["main_class"])
cell_type = "2_CNT_PC"

for cell_type in classes_of_interest:
    print(f"Top DEGs for {cell_type}:")
    df_ct = deg_df[deg_df["main_class"] == cell_type]
    # Select the top 20 DEGs with the highest value of abs(log fold change)
    top_degs = df_ct.loc[df_ct["logfoldchanges"].abs().nlargest(20).index]
    # Print top_degs as markdown table
    print(top_degs[["main_class", "names", "logfoldchanges", "pvals_adj"]].to_markdown(index=False))
    # Save top_degs as csv
    top_degs.to_csv(res_dir+f"top_degs_{cell_type}.csv", index=False)


# --------------------------------------------------------------------------------------
# Compute subcluster DEGs (of milo_signif cells per subcluster)
# --------------------------------------------------------------------------------------

# Check if there's at least one significant cell in each 'subcluster'
significant_cell_types = adata_sub.obs.groupby('main_subcluster')['milo_signif'].apply(lambda x: x.any())
significant_cell_types.to_csv(res_dir+"significant_cell_types_subcluster.csv")

# Normalize and log-transform the data
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)

# Set the 'milo_signif' column as the group for DEGs
adata_sub.obs['milo_signif_str'] = adata_sub.obs['milo_signif'].astype(str)

# Split the data by main_class and compute
adata_sub_ct = [adata_sub[adata_sub.obs['main_subcluster'] == cell_type].copy()
                for cell_type in adata_sub.obs['main_subcluster'].unique()]

# Compute DEGs ---------------------
results = {}
for ind, cell_type in enumerate(adata_sub.obs['main_subcluster'].unique()):
    # Skip if fewer than 2 unique values in 'milo_signif'
    if adata_sub_ct[ind].obs['milo_signif_str'].nunique() != 2:
        print(f"Only one group in milo_signif for cell type {cell_type}, skipping...")
    elif 1 in adata_sub_ct[ind].obs['milo_signif_str'].value_counts().values:
        print(f"Only one cell in a milo_signif group for cell type {cell_type}, skipping...")
    else:
        sc.tl.rank_genes_groups(
            adata_sub_ct[ind],
            groupby='milo_signif_str',
            method='wilcoxon',
            reference='False')
        results[cell_type] = sc.get.rank_genes_groups_df(adata_sub_ct[ind], group='True')

# Plot DEGs ------------------------------------------------------
fig, axs = plt.subplots(19, figsize=(10, 50))
for cell_type in results.keys():
    df = results[cell_type]
    df["is_signif"] = df["pvals_adj"] < 0.05

    # Filter significant DEGs (adj p-value < 0.05)
    signif_genes = df[df["is_signif"]]["names"].iloc[:20].tolist()

    # Skip if no significant DEGs
    if len(signif_genes) == 0:
        print(f"No significant DEGs for {cell_type}")
        ax = axs.flatten()[list(results.keys()).index(cell_type)]
        ax.grid(False)
    else:

        # Get the AnnData subset for this cell type
        # (Assuming adata_ct is a list of AnnData objects corresponding to cell_type order)
        idx = list(adata_sub.obs['main_subcluster'].unique()).index(cell_type)
        adata_subset = adata_sub_ct[idx]

        # Create dotplot
        sc.pl.dotplot(
            adata_subset,
            var_names=signif_genes,
            groupby='milo_signif_str',  # The 2-group comparison
            # title=f'{cell_type}: DEGs (True vs False)',
            ax=axs.flatten()[list(results.keys()).index(cell_type)],
            show=False,
            # figsize=(4, 8 + 0.3*len(signif_genes)),  # Adjust height based on # of genes
            dendrogram=False,
            standard_scale='var',  # Scale expression per gene
            # cmap='RdBu_r',  # Diverging colormap
        )

# Save DEGs -------------------------------------------------------
deg_df_sub = pd.concat(results.values(), keys=results.keys())
deg_df_sub = deg_df_sub.reset_index(level=0).rename(columns={"level_0": "main_subcluster"})

deg_df_sub["is_signif"] = deg_df_sub["pvals_adj"] < 0.05
deg_df_sub["log10_adj_pval"] = -np.log10(deg_df_sub["pvals_adj"])
deg_df_sub["main_subcluster"] = deg_df_sub["main_subcluster"].astype("category")

adata_sub.obs["main_subcluster"].nunique()
deg_df_sub["main_subcluster"].nunique()
len(deg_df_sub["main_subcluster"].unique())

# Plot proportion of DEGs that are significant  -----------

signif_deg_counts_df = deg_df_sub.groupby("main_subcluster")["is_signif"].agg(["sum", "count"])
signif_deg_counts_df = signif_deg_counts_df.reset_index()
signif_deg_counts_df["prop"] = signif_deg_counts_df["sum"] / signif_deg_counts_df["count"]

fig, ax = plt.subplots(figsize=(10, 5))
signif_deg_counts_df.plot(kind="bar", x="main_subcluster", y="prop", ax=ax)
for i, v in enumerate(signif_deg_counts_df["prop"]):
    ax.text(i, v + 0.001, str(int(signif_deg_counts_df["sum"].iloc[i])), ha='center', va='bottom')
ax.set_ylabel("Proportion of DEGs that are significant")
ax.set_xlabel("Subcluster")
ylim = ax.get_ylim()
ax.set_ylim(0, ylim[1] * 1.05)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_title("Proportion of DEGs that are significant per subcluster")
plt.tight_layout()
fig.savefig(res_dir+"milo_signif_DEGs_props_subcluster.png")

# Save significant DEGs ------------------------
deg_df_signif = deg_df_sub[deg_df_sub["is_signif"]]
deg_df_signif["main_class"] = deg_df_signif["main_subcluster"].map(main_subcluster_to_main_class)

# Compute importance_score = abs(logFC) * -log10(pval)
deg_df_signif["importance_score"] = deg_df_signif["logfoldchanges"].abs() * deg_df_signif["log10_adj_pval"]
# Sort by decreasing importance_score
deg_df_signif = deg_df_signif.sort_values(by="importance_score", ascending=False)

deg_df_signif.to_csv(res_dir+"milo_signif_DEGs_subcluster.csv", index=False)

# Violin plot of importance_score by subcluster
fig, ax = plt.subplots(figsize=(10, 5))
sns.violinplot(data=deg_df_signif, x="main_subcluster", y="log10_adj_pval", ax=ax)
ax.set_xlabel("Subcluster", fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_ylabel("Importance score")
ax.set_title("Importance score of DEGs per subcluster")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_DEGs_importance_score_subcluster.png")

# Compute gene specificity score --------------------
degs = deg_df_sub[deg_df_sub["is_signif"]]["names"].unique()

# Initiate a subcluster x gene matrix to zero
subcluster_gene_specif = np.zeros((len(deg_df_sub["main_subcluster"].unique()), len(degs)))

for ind, subcluster in enumerate(deg_df_sub["main_subcluster"].unique()):
    expr_subcluster = adata_sub[adata_sub.obs["main_subcluster"] == subcluster][:, degs].X.sum(axis=0)
    expr_rest = adata_sub[adata_sub.obs["main_subcluster"] != subcluster][:, degs].X.sum(axis=0)
    # compute the log2 fold change
    logFC = np.log2((expr_subcluster + 1e-6) / (expr_rest + 1e-6))
    subcluster_gene_specif[ind, :] = logFC

# Pivot subcluster_gene_specif to a DataFrame with pivot_longer
subcluster_gene_specif_df = pd.DataFrame(subcluster_gene_specif, columns=degs)
subcluster_gene_specif_df["main_subcluster"] = deg_df_sub["main_subcluster"].unique()
subcluster_gene_specif_df = subcluster_gene_specif_df.melt(id_vars=["main_subcluster"], var_name="genes",
                                                           value_name="gene_specif_signed_logFC")
subcluster_gene_specif_df["gene_specif_logFC"] = abs(subcluster_gene_specif_df["gene_specif_signed_logFC"])
subcluster_gene_specif_df["gene_specif_logFC"].hist(bins=50)

# Add gene specificity to deg_df_signif
deg_df_signif = deg_df_signif.merge(subcluster_gene_specif_df, left_on=["main_subcluster", "names"],
                                    right_on=["main_subcluster", "genes"], how="left")

fig, ax = plt.subplots(figsize=(15, 10))
sns.scatterplot(data=deg_df_signif, y="log10_adj_pval", x="gene_specif_logFC", hue="main_class",
                style="main_subcluster", ax=ax, size=0.5)
ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--')
ax.set_xlabel("Gene specificity (logFC of expression in subcluster vs all other cells)")
ax.set_ylabel("-log10(adj_pval) of MILO+ cells")
plt.title("Gene specificity vs -log10(adj_pval) of DEGs")
plt.tight_layout()
fig.savefig(res_dir+"milo_signif_DEGs_gene_specif.png")

# Get all genes with gene specificity
# Order deg_df_signif by main_subcluster and log10_adj_pval (descending)
deg_df_signif = deg_df_signif.sort_values(by=["main_subcluster", "log10_adj_pval"], ascending=[True, False])
deg_df_signif.to_csv(res_dir+"milo_signif_DEGs_subcluster_with_gene_specif.csv", index=False)

# Get the top 100 genes by log10_adj_pval
top_100_degs_pval = deg_df_signif.nlargest(100, "log10_adj_pval")[["main_subcluster", "names", "log10_adj_pval"]]
top_100_degs_pval["main_subcluster"].value_counts()
top_100_degs_pval.to_csv(res_dir+"top_100_genes_df.csv", index=False)

# Get the top 100 genes by gene specificity
top_100_degs_specif = deg_df_signif.nsmallest(
    100, "gene_specif_logFC")[["main_subcluster", "names", "gene_specif_logFC"]]
top_100_degs_specif["main_subcluster"].value_counts()
top_100_degs_specif.to_csv(res_dir+"top_100_genes_specif_df.csv", index=False)

# Get the top genes by log10_adj_pval with gene specificity > 7
top_100_degs_specif = deg_df_signif[deg_df_signif["gene_specif_logFC"] > 7].nlargest(100, "log10_adj_pval")[[
    "main_subcluster", "names", "log10_adj_pval", "gene_specif_logFC"]]
top_100_degs_specif["names"].iloc[0:10].values.tolist()
top_100_degs_specif.to_csv(res_dir+"top_100_genes_specif_higher_7_df.csv", index=False)

# Plot gene specificity vs log10_adj_pval for 2 subclusters of interest
adata_sub_11 = adata_sub[adata_sub.obs["main_class"] == "11_ATL_DTL"].copy()
adata_sub_11_1 = adata_sub[adata_sub.obs["main_subcluster"] == "11.1: ATL_DTL_2"].copy()
adata_sub_4_2 = adata_sub[adata_sub.obs["main_subcluster"] == "4.2: aTAL/aDCT"].copy()
adata_sub_4_3 = adata_sub[adata_sub.obs["main_subcluster"] == "4.3: aTAL/aDCT/aPT"].copy()

# --------------------------------------------------------------------------------------
# Plot MILO scores per disease_id
# --------------------------------------------------------------------------------------

# Plot milo_scores per disease_id for each class -----------------------
fig, axs = plt.subplots(3, 2, figsize=(20, 15))
classes_of_interest
for ind, myclass in enumerate(classes_of_interest):
    adata_tmp = adata[adata.obs["main_class"] == myclass].copy()
    sc.pl.violin(adata_tmp, "milo_adj_log10pval", groupby="disease_id", show=False,
                 ax=axs.flatten()[ind], size=2)
    axs.flatten()[ind].set_title(myclass)
axs.flatten()[-1].set_frame_on(False)
axs.flatten()[-1].set_xticks([])
axs.flatten()[-1].set_yticks([])
plt.suptitle("Milo scores per disease_id, main classes", y=0.95)
plt.savefig(res_dir+"milo_scores_per_disease_id_main.png")

# Plot milo_scores per disease_id for each main_subcluster -----------------------
subclusters_of_interest = deg_df_signif["main_subcluster"].unique().tolist()
fig, axs = plt.subplots(4, 4, figsize=(20, 15))
for ind, mymain_subcluster in enumerate(subclusters_of_interest):
    adata_tmp = adata_sub[adata_sub.obs["main_subcluster"] == mymain_subcluster].copy()
    sc.pl.violin(adata_tmp, "milo_adj_log10pval", groupby="disease_id", show=False,
                 ax=axs.flatten()[ind], size=2)
    axs.flatten()[ind].set_title(mymain_subcluster)
    axs.flatten()[ind].set_xlabel("")
plt.tight_layout()
plt.suptitle("Milo scores (log10(adj_pval) per disease_id, subclusters", y=1.03)
plt.savefig(res_dir+"milo_scores_per_disease_id_subcluster.png")


# --------------------------------------------------------------------------------------
# Interpret subcluster DEGs
# --------------------------------------------------------------------------------------

# Plot histograms of log10_adj_pval colored by subclusters -----------
fig, ax = plt.subplots(figsize=(10, 6))
sns.histplot(data=deg_df_signif, x="log10_adj_pval", hue="subcluster", bins=50, ax=ax)
ax.set_xlabel("-log10(adj_pval)")
ax.set_ylabel("Number of DEGs")
ax.set_title("Distribution of DEGs -log10(adj_pval) per subcluster")
plt.tight_layout()
fig.savefig(res_dir+"milo_signif_DEGs_hist.png")

# Compute props of DEGs above a certain log10_adj_pval threshold per subcluster -----------
all_subclusters = deg_df_signif['subcluster'].unique()  # List of all unique subclusters
log10_pval_vect = range(1, deg_df_signif["log10_adj_pval"].astype(float).max().astype(int), 1)

n_subclusters_per_pval_ls = []
for log10_pval_val in log10_pval_vect:
    # Get the DEGs with log10_adj_pval > log10_pval_val
    degs = deg_df_signif[deg_df_signif["log10_adj_pval"] > log10_pval_val]
    # log10_pval_val = 10
    # Get the number of DEGs
    n_degs = degs.shape[0]
    # Get the props of subclusters
    subcluster_counts = degs["subcluster"].value_counts(dropna=False)
    subcluster_counts = subcluster_counts.reindex(all_subclusters, fill_value=0)
    props = subcluster_counts / n_degs
    props = props.reset_index()
    props.rename(columns={"count": "prop"}, inplace=True)
    props["counts"] = subcluster_counts.values
    props["log10_pval_val"] = log10_pval_val
    n_subclusters_per_pval_ls.append(props)

n_subclusters_per_pval_df = pd.concat(n_subclusters_per_pval_ls, axis=0)
n_subclusters_per_pval_df

# Create the stacked plot of props of DEGs by subcluster -----------
df_pivoted = n_subclusters_per_pval_df.pivot_table(index='log10_pval_val', columns='subcluster',
                                                   values='prop', aggfunc='sum').fillna(0)

colors = plt.cm.tab20(np.linspace(0, 1, len(df_pivoted.columns)))
fig, ax = plt.subplots(figsize=(10, 6))
ax.stackplot(df_pivoted.index, df_pivoted.T.values, labels=df_pivoted.columns, colors=colors, alpha=0.7)
ax.set_xlabel('log10(p-value)')
ax.set_ylabel('Proportion of subclusters with DEGs > log10(p-value)')
ax.set_title('Stacked Ribbon of props of DEGs by Subcluster')
ax.legend(title='Subcluster', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
plt.savefig(res_dir+"milo_signif_DEGs_stacked_props.png")

# Create the stacked plot of counts of DEGs by subcluster -----------
df_pivoted_counts = n_subclusters_per_pval_df.pivot_table(index='log10_pval_val', columns='subcluster',
                                                          values='counts', aggfunc='sum').fillna(0)

colors_counts = plt.cm.tab20(np.linspace(0, 1, len(df_pivoted_counts.columns)))
fig, ax = plt.subplots(figsize=(10, 6))
ax.stackplot(df_pivoted_counts.index, df_pivoted_counts.T.values, labels=df_pivoted_counts.columns,
             colors=colors_counts, alpha=0.7)
ax.set_xlabel('log10(p-value)')
ax.set_ylabel('Counts')
ax.set_title('Stacked Ribbon of Counts of DEGs by Subcluster')
ax.legend(title='Subcluster', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
plt.savefig(res_dir+"milo_signif_DEGs_stacked_counts.png")

# Compute number of DEGs by pvalue
fig, ax = plt.subplots(figsize=(10, 6))
n_subclusters_per_pval_df.groupby("log10_pval_val").agg({"counts": "sum"}).reset_index()["counts"].hist(bins=100, ax=ax)
ax.set_xlabel("-log10(adj_pval)")
ax.set_ylabel("Number of DEGs")
ax.set_title("Number of DEGs by log10(p-value)")
plt.tight_layout()
plt.savefig(res_dir+"milo_signif_DEGs_by_pval.png")

# --------------------------------------------------------------------------------------
# Compute Correlation of genes with milo_log10pval
# --------------------------------------------------------------------------------------

# Get genes most correlated with milo_log10pval
# Compute correlation (Spearman or Pearson)
gene_names = adata.var_names

gene_names = gene_names[0:10]

X = adata.raw.X if adata.raw is not None else adata.X
gene_names = adata.raw.var_names if adata.raw is not None else adata.var_names
score = adata.obs["milo_adj_log10pval"].values

correlations = []

for i in range(10):
    expr = X[:, i].toarray().ravel()  # only one column at a time in memory
    corr, _ = pearsonr(expr, score)
    correlations.append(corr)

cor_df = pd.DataFrame({
    "gene": gene_names,
    "pearson_corr": correlations
}).sort_values("pearson_corr", ascending=False)

# Save df as csv
cor_df.to_csv(res_dir+"milo_log10pval_correlations.csv", index=False)


# Plot milo_scores per disease_id for each class with p-values
fig, axs = plt.subplots(3, 2, figsize=(20, 15))
classes_of_interest = adata.obs["main_class"].unique()  # Ensure this matches your classes

for ind, myclass in enumerate(classes_of_interest):
    # Subset data for the current class
    adata_tmp = adata[adata.obs["main_class"] == myclass].copy()
    ax = axs.flatten()[ind]

    # Generate violin plot
    sc.pl.violin(
        adata_tmp,
        "milo_adj_log10pval",
        groupby="disease_id",
        show=False,
        ax=ax,
        size=2
    )
    ax.set_title(myclass)

    # Compute pairwise Wilcoxon p-values
    disease_ids = adata_tmp.obs["disease_id"].unique()
    if len(disease_ids) < 2:
        continue  # Skip if only 1 group

    # Get all unique pairs of disease_ids
    pairs = list(combinations(disease_ids, 2))
    p_values = []

    for pair in pairs:
        group1 = adata_tmp[adata_tmp.obs["disease_id"] == pair[0]].obs["milo_adj_log10pval"]
        group2 = adata_tmp[adata_tmp.obs["disease_id"] == pair[1]].obs["milo_adj_log10pval"]
        _, pval = ranksums(group1, group2)
        p_values.append((pair, pval))

    # Annotate significant pairs (adjust threshold as needed)
    y_max = adata_tmp.obs["milo_adj_log10pval"].max() + 0.5
    y_step = 0.1 * y_max

    for i, (pair, pval) in enumerate(p_values):
        if pval < 0.05:  # Significance threshold
            x1, x2 = np.where(disease_ids == pair[0])[0][0], np.where(disease_ids == pair[1])[0][0]
            y_pos = y_max + (i * y_step)

            # Draw significance bar
            ax.plot([x1, x1, x2, x2], [y_pos-0.05*y_step, y_pos, y_pos, y_pos-0.05*y_step], lw=1, color='k')

            # Add p-value asterisks
            n_stars = sum([pval < cutoff for cutoff in [0.05, 0.01, 0.001]])
            ax.text((x1 + x2) * 0.5, y_pos + 0.05*y_step, '*' * n_stars,
                    ha='center', va='bottom', color='k')

# Adjust layout
axs.flatten()[-1].set_frame_on(False)
axs.flatten()[-1].set_xticks([])
axs.flatten()[-1].set_yticks([])
plt.suptitle("Milo scores per disease_id (Wilcoxon p-values: *<0.05, **<0.01, ***<0.001)", y=0.95)
plt.tight_layout()
plt.savefig(res_dir + "milo_scores_per_disease_id_main_with_pvals.png")
