#
# Created on Tue May 13 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script processes MILO-derived subcluster DEGs from kidney biopsy scRNA-seq data,
# identifies top DEGs per subcluster, performs pathway enrichment, and generates plots.
#
# Input datasets/objects:
# - milo_signif_DEGs_subcluster.csv
# - milo_signif_DEGs_subcluster_with_gene_specif.csv
#
# Output datasets/files:
# - top_100_DEGs_subcluster_by_logFC.csv
# - top_100_DEGs_subcluster_by_pval.csv
# - top_10_DEGs_per_subcluster_by_logFC.csv
# - top_10_DEGs_per_subcluster_by_pval.csv
# - top_100_DEGs_per_subcluster_by_logFC.csv
# - top_100_DEGs_per_subcluster_by_pval.csv
# - enrichment_results_top_10_DEGs.csv
# - enrichment_results_top_100_DEGs.csv
# - enrichment_results_top_10_DEGs_all_pws.pdf
# - enrichment_results_top_100_DEGs_top10_pws.pdf
# - enrichment_results_top_100_DEGs_top20_pws.pdf
#
# Analyses performed:
# - Filter DEGs for significance (adjusted p-value < 0.05)
# - Identify top DEGs:
#     - Overall top 100 DEGs by logFC
#     - Top 100 DEGs by adjusted p-value (upregulated genes)
#     - Top 10 DEGs per subcluster by logFC and by p-value
#     - Top 100 DEGs per subcluster by logFC and by p-value
# - Merge with original DEG tables to retain all columns
# - Compute pathway enrichment using Scanpy and GProfiler:
#     - GO:BP and KEGG pathways for top DEGs
#     - Save enrichment results to CSV
#     - Barplots of -log10(p-value) for top enriched pathways
# - Generate plots:
#     - Scatterplots of logFC vs -log10(adjusted p-value)
#     - Histograms of log10(adjusted p-value) for top DEGs
#     - Volcano plots for each subcluster (logFC vs -log10 adjusted p-value)
# - Optional GSEA using gseapy prerank on ranked gene lists per subcluster
#
# Visualization outputs:
# - Barplots of pathway enrichment (-log10 p-value)
# - Scatterplots of logFC vs log10 adjusted p-value
# - Histograms of adjusted p-values
# - Volcano plots per subcluster
#
# The script systematically summarizes subcluster-level DEG significance,
# identifies top genes for downstream analysis, evaluates functional enrichment,
# and produces publication-ready plots for top DEGs and enriched pathways.
# --------------------------------------------------------------------------------------------


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import gseapy as gp
import textwrap

proj_dir = '/home/zd775156/kidney_biopsy/'
res_dir = proj_dir + "results/milo/"

deg_df_signif_sub = pd.read_csv(res_dir+"milo_signif_DEGs_subcluster.csv", index_col=0)

# Set variable subcluster as index
deg_df_signif_sub.reset_index(drop=True, inplace=True)
deg_df_signif_sub["subcluster"] = deg_df_signif_sub.index

deg_df_signif_sub = pd.read_csv(res_dir+"milo_signif_DEGs_subcluster_with_gene_specif.csv")
deg_df_signif_sub = deg_df_signif_sub[deg_df_signif_sub["pvals_adj"] < 0.05]

# Get the top 100 DEGs (by logfoldchanges, among signif pvals) ----------
top_n = 100
top_degs_by_logFC = deg_df_signif_sub[deg_df_signif_sub["pvals_adj"] < 0.05].sort_values(by='logfoldchanges',
                                                                                         ascending=False).head(top_n)

# Save DEGs to csv
top_degs_by_logFC.to_csv(res_dir+"top_100_DEGs_subcluster_by_logFC.csv")

# Get the top 100 DEGs by adjusted p-value (only upregulated genes) ----------
top_degs_pval = deg_df_signif_sub[deg_df_signif_sub["logfoldchanges"] > 0].sort_values(by='pvals_adj').head(top_n)

top_degs_pval["main_subcluster"].value_counts()
# Save to csv
top_degs_pval.to_csv(res_dir+"top_100_DEGs_subcluster_by_pval.csv")

# -----------------------------------------------------------------------------------------
# Get the top 10 DEGs for each subcluster (by logfoldchanges, among signif DEGs) ----------
# -----------------------------------------------------------------------------------------

top_n_per_subcluster = 10
# Get count of upregulated DEGs per subcluster
count = deg_df_signif_sub[deg_df_signif_sub["logfoldchanges"] > 0].groupby("main_subcluster").size()
# Get the top 10 main_subclusters by count of upregulated DEGs
top_subclusters = count.nlargest(top_n_per_subcluster).index.tolist()

# For each top_subcluster, get the top 10 DEGs by logfoldchanges, keeping all columns from deg_df_signif_sub
top_degs_per_subcluster_logFC = (
    deg_df_signif_sub[deg_df_signif_sub["main_subcluster"].isin(top_subclusters)]
    .sort_values(by='logfoldchanges', ascending=False)
    .groupby('main_subcluster', as_index=False)
    .apply(lambda x: x.head(top_n_per_subcluster))
    .reset_index(drop=True)
)
# Left join with the original deg_df_signif_sub to keep all columns
top_degs_per_subcluster_logFC = pd.merge(
    top_degs_per_subcluster_logFC,
    deg_df_signif_sub,
    on=['main_subcluster', 'names', 'logfoldchanges', 'pvals_adj'],
    how='left',
    suffixes=('', '_original')
)

# For each top_subcluster, get the top 10 DEGs by adjusted p-value, keeping all columns from deg_df_signif_sub
top_degs_per_subcluster_pval = (
    deg_df_signif_sub[deg_df_signif_sub["main_subcluster"].isin(top_subclusters)]
    .sort_values(by='pvals_adj')
    .groupby('main_subcluster', as_index=False)
    .apply(lambda x: x.head(top_n_per_subcluster))
    .reset_index(drop=True)
)
# Left join with the original deg_df_signif_sub to keep all columns
top_degs_per_subcluster_pval = pd.merge(
    top_degs_per_subcluster_pval,
    deg_df_signif_sub,
    on=['main_subcluster', 'names', 'logfoldchanges', 'pvals_adj'],
    how='left',
    suffixes=('', '_original')
)

# Save to csv -------------------
top_degs_per_subcluster_logFC.to_csv(res_dir+"top_10_DEGs_per_subcluster_by_logFC.csv", index=False)
top_degs_per_subcluster_pval.to_csv(res_dir+"top_10_DEGs_per_subcluster_by_pval.csv", index=False)


# ----------------------------------------------------------------
# Get pathways enriched in top DEGs ------------------------------
# ----------------------------------------------------------------

gene_list = top_degs_per_subcluster_logFC["names"].values.tolist()
enrich_df = sc.queries.enrich(gene_list, org="hsapiens", gprofiler_kwargs={"sources": ["GO:BP", "KEGG"]})
enrich_df["-log10(p_value)"] = -np.log10(enrich_df["p_value"])
enrich_df["source_name"] = enrich_df['source'] + ": " + enrich_df['name']
# Save the enrichment results to csv
enrich_df.to_csv(res_dir+"enrichment_results_top_10_DEGs.csv", index=False)

# Wrap y-labels for better readability
enrich_df['source_name_wrapped'] = enrich_df['source_name'].apply(lambda x: '\n'.join(textwrap.wrap(x, 30)))
sns.barplot(y='source_name_wrapped', x='-log10(p_value)', data=enrich_df)
plt.xlabel('-log10(p‑value)')
plt.ylabel('')
plt.tight_layout()
# Save the plot
plt.savefig(res_dir+"enrichment_results_top_10_DEGs_all_pws.pdf", dpi=300, bbox_inches='tight')


# -----------------------------------------------------------------------------------------
# Get the top 100 DEGs for each subcluster (by logfoldchanges, among signif DEGs) ----------
# -----------------------------------------------------------------------------------------

top_n_per_subcluster = 100
# Get count of upregulated DEGs per subcluster
count_100 = deg_df_signif_sub[deg_df_signif_sub["logfoldchanges"] > 0].groupby("main_subcluster").size()
# Get the top 100 main_subclusters by count of upregulated DEGs
top_subclusters_100 = count_100.nlargest(top_n_per_subcluster).index.tolist()

# For each top_subcluster, get the top 100 DEGs by logfoldchanges, keeping all columns from deg_df_signif_sub
top_degs_per_subcluster_logFC_100 = (
    deg_df_signif_sub[deg_df_signif_sub["main_subcluster"].isin(top_subclusters_100)]
    .sort_values(by='logfoldchanges', ascending=False)
    .groupby('main_subcluster', as_index=False)
    .apply(lambda x: x.head(top_n_per_subcluster))
    .reset_index(drop=True)
)
top_degs_per_subcluster_logFC_100 = pd.merge(
    top_degs_per_subcluster_logFC_100,
    deg_df_signif_sub,
    on=['main_subcluster', 'names', 'logfoldchanges', 'pvals_adj'],
    how='left',
    suffixes=('', '_original')
)

# For each top_subcluster, get the top 100 DEGs by adjusted p-value, keeping all columns from deg_df_signif_sub
top_degs_per_subcluster_pval_100 = (
    deg_df_signif_sub[deg_df_signif_sub["main_subcluster"].isin(top_subclusters_100)]
    .sort_values(by='pvals_adj')
    .groupby('main_subcluster', as_index=False)
    .apply(lambda x: x.head(top_n_per_subcluster))
    .reset_index(drop=True)
)
top_degs_per_subcluster_pval_100 = pd.merge(
    top_degs_per_subcluster_pval_100,
    deg_df_signif_sub,
    on=['main_subcluster', 'names', 'logfoldchanges', 'pvals_adj'],
    how='left',
    suffixes=('', '_original')
)

# Save to csv
top_degs_per_subcluster_logFC_100.to_csv(res_dir+"top_100_DEGs_per_subcluster_by_logFC.csv", index=False)
top_degs_per_subcluster_pval_100.to_csv(res_dir+"top_100_DEGs_per_subcluster_by_pval.csv", index=False)

top_degs_per_subcluster_logFC_100["main_subcluster"].value_counts()
top_degs_per_subcluster_pval_100["main_subcluster"].value_counts()

# PLot log10_adj_pval vs gene_specif_logFC_original in top_degs_per_subcluster_logFC_100 using seaborn scatterplot
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=top_degs_per_subcluster_pval_100,
    x='logfoldchanges',
    y='log10_adj_pval',
    hue='main_subcluster',
    alpha=0.7
)
plt.axhline(y=0.05, color='red', linestyle='--')
plt.axvline(x=0, color='black', linestyle='--')

# PLot histogram of log10_adj_pval for top_degs_per_subcluster_logFC_100
plt.figure(figsize=(10, 6))
sns.histplot(
    data=top_degs_per_subcluster_logFC_100,
    x='log10_adj_pval',
    bins=50,
    kde=True,
    stat='density',
    common_norm=False
)
plt.title('Histogram of log10 Adjusted P-values for Top 100 DEGs per Subcluster')
plt.savefig(res_dir+"histogram_log10_adj_pval_top_degs_per_subcluster_logFC_100.png", dpi=300, bbox_inches='tight')


# ----------------------------------------------------------------
# Get pathways enriched in top DEGs ------------------------------
# ----------------------------------------------------------------

gene_list_100 = top_degs_per_subcluster_logFC_100["names"].values.tolist()
enrich_df_100 = sc.queries.enrich(gene_list_100, org="hsapiens", gprofiler_kwargs={"sources": ["GO:BP", "KEGG"]})

enrich_df_100["-log10(p_value)"] = -np.log10(enrich_df_100["p_value"])
enrich_df_100["source_name"] = enrich_df_100['source'] + ": " + enrich_df_100['name']
# Save the enrichment results to csv
enrich_df_100.to_csv(res_dir+"enrichment_results_top_100_DEGs.csv", index=False)

enrich_df_100_top10 = enrich_df_100.nlargest(10, '-log10(p_value)')
# Wrap y-labels for better readability
enrich_df_100_top10['source_name_wrapped'] = (enrich_df_100_top10['source_name'].
                                              apply(lambda x: '\n'.join(textwrap.wrap(x, 30))))
sns.barplot(y='source_name_wrapped', x='-log10(p_value)', data=enrich_df_100_top10)
plt.xlabel('-log10(p‑value)')
plt.ylabel('')
plt.tight_layout()
# Save the plot
plt.savefig(res_dir+"enrichment_results_top_100_DEGs_top10_pws.pdf", dpi=300, bbox_inches='tight')

enrich_df_100_top20 = enrich_df_100.nlargest(20, '-log10(p_value)')
# Wrap y-labels for better readability
enrich_df_100_top20['source_name_wrapped'] = (enrich_df_100_top20['source_name'].
                                              apply(lambda x: '\n'.join(textwrap.wrap(x, 30))))
sns.barplot(y='source_name_wrapped', x='-log10(p_value)', data=enrich_df_100_top20)
plt.xlabel('-log10(p‑value)')
# set y axis text font size
plt.tick_params(axis='y', labelsize=8)
plt.ylabel('')
plt.tight_layout()
# Save the plot
plt.savefig(res_dir+"enrichment_results_top_100_DEGs_top20_pws.pdf", dpi=300, bbox_inches='tight')


# -----------------------------------------------------------------------------------------
# For each subcluster, plot the volcano plot of the DEGs -----------------------------------
# -----------------------------------------------------------------------------------------


# Get unique cell types
cell_types = deg_df_signif_sub["main_subcluster"].unique()
n_cell_types = len(cell_types)

# Set up the subplot grid
ncols = 3
nrows = int(np.ceil(n_cell_types / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False)

for idx, cell_type in enumerate(cell_types):
    row = idx // ncols
    col = idx % ncols
    ax = axes[row, col]

    # Filter for the current cell type
    cell_type_df = deg_df_signif_sub[deg_df_signif_sub["main_subcluster"] == cell_type].copy()
    cell_type_df["significant"] = cell_type_df["pvals_adj"] < 0.05

    # Volcano plot
    sns.scatterplot(
        data=cell_type_df,
        x='logfoldchanges',
        y=-np.log10(cell_type_df['pvals_adj']),
        hue='significant',
        alpha=0.7,
        ax=ax,
        legend=False
    )
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--')
    ax.axvline(x=0, color='black', linestyle='--')
    ax.set_title(f'Volcano Plot for {cell_type}')
    ax.set_xlabel('Log Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')

# Remove empty subplots
for idx in range(n_cell_types, nrows * ncols):
    fig.delaxes(axes[idx // ncols, idx % ncols])

plt.tight_layout()
plt.show()


# -----------------------------------------------------------------------------
# Scanpy's enrich analysis ----------------------------------------------------
# -----------------------------------------------------------------------------

for cell_type in deg_df_signif_sub["subcluster"].unique():
    # Filter for the current cell type
    cell_type_df = deg_df_signif_sub[deg_df_signif_sub["subcluster"] == cell_type]

    # Get the list of genes
    gene_list = cell_type_df["names"].tolist()

    # ORA -----------------------------------------------------
    sc.queries.enrich(gene_list, org="hsapiens")

    # GSEA -----------------------------------------------------
    # Example: Assume your dataframe contains 'gene', 'logFC' (log fold change), and 'p_val' columns
    # Rank the genes by logFC (or any other metric) for GSEA
    ranked_gene_list = (deg_df_signif_sub[['names', "pvals_adj", 'logfoldchanges']].
                        set_index('names').sort_values(by='pvals_adj'))

    # Save this dataframe in the format expected by gseapy
    ranked_gene_list = ranked_gene_list.reset_index()

    # Run the prerank method
    pre_res = gp.prerank(rnk=ranked_gene_list,
                         gene_sets='KEGG_2016',
                         outdir='gsea_results')
