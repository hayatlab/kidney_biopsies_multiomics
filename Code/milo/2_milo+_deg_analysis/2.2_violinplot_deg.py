#
# Created on Tue May 13 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script processes top DEGs from MILO analysis of kidney biopsy scRNA-seq data,
# computes correlation matrices, and retrieves gene annotation information from MyGene.info.
#
# Datasets used/loaded:
# - scVI-normalized AnnData: scVI_kidney_biopsy_with_clinical_no_TIN.h5ad
# - Top DEGs CSVs:
#     - top_10_DEGs_per_subcluster_by_logFC.csv
#     - top_10_DEGs_per_subcluster_by_pval.csv
#     - top_100_DEGs_subcluster_by_logFC.csv
#     - top_100_DEGs_subcluster_by_pval.csv
#
# Datasets/objects saved:
# - Subsets of AnnData for top DEGs (by logFC and p-value)
# - Correlation matrices (as PNG):
#     - correlation_matrix_top10_DEGs_by_pval.png
#     - correlation_matrix_top100_DEGs_by_pval.png
#     - correlation_matrix_top10_DEGs_by_logFC.png
#     - correlation_matrix_top100_DEGs_by_logFC.png
# - Gene annotation summaries (CSV) per DEG list:
#     - deg_info_top10_pval_summary.csv
#     - deg_info_top100_pval_summary.csv
#     - deg_info_top10_logFC_summary.csv
#     - deg_info_top100_logFC_summary.csv
#
# Analyses performed:
# - Load AnnData object and store raw counts in layer
# - Normalize total counts and log-transform data
# - Load top DEGs by subcluster (logFC and p-value) for top 10 and top 100 genes
# - Remove duplicates and subset AnnData to top DEGs
# - Compute and plot correlation matrices of gene expression for each DEG list
# - Retrieve detailed gene annotations from MyGene.info (symbol, summary, GO terms, UniProt IDs)
# - Save summarized gene annotation info to CSV per DEG list
#
# Visualization outputs:
# - Heatmaps of correlation matrices for each DEG list (top10/top100, logFC/p-value)
#
# The script links DEG identification with correlation structure analysis and
# gene annotation retrieval, producing both visualization and tabular summaries
# for downstream functional interpretation.
# --------------------------------------------------------------------------------------------


import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.sparse import csr_matrix

proj_dir = '/home/zd775156/kidney_biopsy/'
res_dir = proj_dir + "results/milo/"

# ----------------------------------------------------------------------
# Load data ------------------------------------------------------------
# ----------------------------------------------------------------------

adata = sc.read_h5ad(proj_dir+'data/scVI_kidney_biopsy_with_clinical_no_TIN.h5ad')

adata.layers["counts"] = adata.X.copy()  # Store raw counts in a layer
# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)  # Log-transform the data

# ----------------------------------------------------------------------
# Read DEG csv ---------------------------------------------------------
# ----------------------------------------------------------------------

top_degs_per_subcluster_logFC = pd.read_csv(res_dir+"top_10_DEGs_per_subcluster_by_logFC.csv")
top_degs_per_subcluster_pval = pd.read_csv(res_dir+"top_10_DEGs_per_subcluster_by_pval.csv")
top_degs_by_logFC = pd.read_csv(res_dir+"top_100_DEGs_subcluster_by_logFC.csv")
top_degs_pval = pd.read_csv(res_dir+"top_100_DEGs_subcluster_by_pval.csv")

deg_top10_pval = top_degs_per_subcluster_pval["names"].values.tolist()
deg_top10_pval = list(set(deg_top10_pval))  # Remove duplicates
print(len(deg_top10_pval), "DEGs found in top 10 DEGs by p-value")

deg_top100_pval = top_degs_pval["names"].values.tolist()
deg_top100_pval = list(set(deg_top100_pval))
print(len(deg_top100_pval), "DEGs found in top 100 DEGs by p-value")

deg_top10_logFC = top_degs_per_subcluster_logFC["names"].values.tolist()
deg_top10_logFC = list(set(deg_top10_logFC))  # Remove duplicates
print(len(deg_top10_logFC), "DEGs found in top 10 DEGs by logFC")

deg_top100_logFC = top_degs_by_logFC["names"].values.tolist()
deg_top100_logFC = list(set(deg_top100_logFC))  # Remove duplicates
print(len(deg_top100_logFC), "DEGs found in top 100 DEGs by logFC")

# Subset anndata to only include the top DEGs by p-value
adata_deg_top10_pval = adata[:, deg_top10_pval].copy()
adata_deg_top100_pval = adata[:, deg_top100_pval].copy()
adata_deg_top10_logFC = adata[:, deg_top10_logFC].copy()
adata_deg_top100_logFC = adata[:, deg_top100_logFC].copy()


# Plot correlation matrices of top 10 DEGs by p-value
def plot_correlation_matrix(adata, genes, title, savefile=None):
    """Plot a correlation matrix for the given genes in the anndata object."""
    data = adata[:, genes].X
    if isinstance(data, csr_matrix):
        data = data.toarray()  # Convert sparse matrix to dense if necessary
    corr_matrix = np.corrcoef(data, rowvar=False)
    plt.figure(figsize=(20, 16))
    sns.heatmap(corr_matrix, annot=False, fmt=".1f",
                cmap='coolwarm',
                xticklabels=genes, yticklabels=genes)
    plt.title(title)
    if savefile:
        plt.savefig(savefile, dpi=300, bbox_inches='tight')
    plt.show()


# Plot correlation matrix for top 10 DEGs by p-value
plot_correlation_matrix(adata_deg_top10_pval,
                        deg_top10_pval,
                        "Correlation Matrix of Top 10 DEGs by p-value",
                        res_dir+"correlation_matrix_top10_DEGs_by_pval.png")

# Plot correlation matrix for top 100 DEGs by p-value
plot_correlation_matrix(adata_deg_top100_pval,
                        deg_top100_pval,
                        "Correlation Matrix of Top 100 DEGs by p-value",
                        res_dir+"correlation_matrix_top100_DEGs_by_pval.png")

# Plot correlation matrix for top 10 DEGs by logFC
plot_correlation_matrix(adata_deg_top10_logFC,
                        deg_top10_logFC,
                        "Correlation Matrix of Top 10 DEGs by logFC",
                        res_dir+"correlation_matrix_top10_DEGs_by_logFC.png")

# Plot correlation matrix for top 100 DEGs by logFC
plot_correlation_matrix(adata_deg_top100_logFC,
                        deg_top100_logFC,
                        "Correlation Matrix of Top 100 DEGs by logFC",
                        res_dir+"correlation_matrix_top100_DEGs_by_logFC.png")


# ----------------------------------------------------------------------
# Retrieve gene info
# ----------------------------------------------------------------------


# Retrieve detailed gene information including GO terms and UniProt IDs
def get_gene_info_detailed(symbol):
    url = f"https://mygene.info/v3/query?q=symbol:{symbol}&species=human" + \
        "&fields=symbol,name,summary,entrezgene,go,uniprot,interpro"
    response = requests.get(url)
    hits = response.json().get("hits", [])
    return hits[0] if hits else None


# Get gene info from the top DEGs for each list -------------------------
deg_lists = [
    ("top10_pval", deg_top10_pval),
    ("top10_pval", deg_top100_pval),
    ("top10_logFC", deg_top10_logFC),
    ("top100_logFC", deg_top100_logFC)
]

for label, deg_list in deg_lists:
    deg_info = [get_gene_info_detailed(gene) for gene in deg_list]

    all_keys = set()
    for d in deg_info:
        if isinstance(d, dict):
            all_keys.update(d.keys())

    deg_info_filled = [
        d if isinstance(d, dict) else dict.fromkeys(all_keys, None)
        for d in deg_info
    ]

    deg_info_df = pd.DataFrame(deg_info_filled)

    # Save summary
    summary_df = deg_info_df[["symbol", "summary"]].copy()
    summary_df.to_csv(res_dir+f"deg_info_{label}_summary.csv", index=False)
