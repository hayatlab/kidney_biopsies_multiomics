#
# Created on Sat May 17 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script performs Milo-based differential abundance (DA) analysis on scVI-integrated
# kidney biopsy scRNA-seq data, generating neighborhood graphs, DA statistics, and UMAP visualizations.
#
# Datasets used/loaded:
# - scVI-integrated AnnData: data/scVI_kidney_biopsy_with_clinical_no_TIN.h5ad
#
# Datasets saved:
# - Milo neighborhood object with DA test results: results/milo/milo_obj.h5mu
#
# Plots generated (by output folder):
#
# results/milo/:
#   - umap_classification_jm_slope.png:
#       UMAPs colored by clinical classification, main cell class, and jm_slope.
#   - umap_main_class.png:
#       UMAP colored by main cell class.
#   - nhood_size.png:
#       Histogram of the number of cells per Milo neighborhood.
#   - da_volcano_plot.png:
#       Volcano plot of logFC vs -log10(SpatialFDR) for DA testing (~classification).
#   - nhood_umap.png:
#       Milo neighborhood graph colored by logFC (DA results visualization).
#   - nhood_annotation_frac.png:
#       Histogram of dominant cell-type proportions per neighborhood.
#
# Analyses performed:
# 1. Data preparation:
#    - Loaded scVI-integrated data; standardized labels ("Fast_Progressor" naming).
#    - Removed MALAT1 gene.
#
# 2. Visualization:
#    - Generated UMAPs for classification, disease, main class, and jm_slope.
#
# 3. Milo setup and neighborhood construction:
#    - Initialized Milo object and loaded AnnData.
#    - Built neighbor graph using X_scVI embedding (n_neighbors=150).
#    - Defined 10% of cells as neighborhood seeds (prop=0.1).
#    - Plotted distribution of neighborhood sizes.
#
# 4. Differential abundance analysis:
#    - Performed DA test comparing Fast_Progressor vs stable:
#         design: "~classification"
#         contrast: "classificationFast_Progressor - classificationstable"
#    - Generated volcano and UMAP plots for DA results.
#
# 5. Annotation and output:
#    - Annotated neighborhoods by dominant cell type (main_class).
#    - Visualized neighborhood composition.
#    - Saved processed Milo object with DA results to disk.
#
# The script outputs Milo neighborhood statistics, DA results, annotated objects, and diagnostic plots
# summarizing neighborhood-level abundance shifts across clinical groups.
# --------------------------------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import pertpy as pt

proj_dir = '/home/zd775156/kidney_biopsy/'
res_dir = proj_dir + "results/milo/"

# ----------------------------------------------------------------------------
# Load data  -----------------------------------------------------------------
# ----------------------------------------------------------------------------

adata = sc.read_h5ad(proj_dir+'data/scVI_kidney_biopsy_with_clinical_no_TIN.h5ad')

# Change Fast Progressor to Fast_Progressor
adata.obs['classification'] = adata.obs['classification'].str.replace(
    "Fast Progressor", "Fast_Progressor").astype("category")

# Remove MALAT1 gene from adata
if "MALAT1" in adata.var_names:
    adata = adata[:, adata.var_names != "MALAT1"]

# -----------------------------------------------------------------------------
# Plot data on UMAPs ----------------------------------------------------------
# -----------------------------------------------------------------------------

sc.pl.umap(adata, color=["classification", "disease_id",
           "main_class"], show=False, frameon=False, ncols=1)

fig, axs = plt.subplots(3, 1, figsize=(10, 30))
sc.pl.umap(adata, color=["classification"], frameon=False, ncols=1, legend_loc="none",
           ax=axs.flatten()[0], show=False)
sc.pl.umap(adata, color=["main_class"], frameon=False, ncols=1, legend_loc="on data",
           ax=axs.flatten()[1], show=False)
sc.pl.umap(adata, color=["jm_slope"], frameon=False, ncols=1, legend_loc="on data",
           ax=axs.flatten()[2], show=False)
plt.tight_layout()
fig.savefig(res_dir+"umap_classification_jm_slope.png")

fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(adata, color=["main_class"], frameon=False, ncols=1, ax=ax, show=False)
plt.tight_layout()
plt.suptitle("UMAP with main class", y=1.05)
fig.savefig(res_dir+"umap_main_class.png")

# -----------------------------------------------------------------------------
# Initialize object for Milo analysis -----------------------------------------
# -----------------------------------------------------------------------------

milo = pt.tl.Milo()
mdata = milo.load(adata)

# Recommended value for milo: [135, 225]
sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150)

# -----------------------------------------------------------------------------
# Construct neighborhoods -----------------------------------------------------
# -----------------------------------------------------------------------------

# Randomly sample vertices on a KNN graph to define neighbourhoods of cells.

# milo.make_nhoods(mdata["rna"], prop=0.05, seed=0)
milo.make_nhoods(mdata["rna"], prop=0.1, seed=0)

# Print neighborhood information
(mdata["rna"][mdata["rna"].obs["nhood_ixs_refined"] != 0]
 .obs[["nhood_ixs_refined", "nhood_kth_distance"]])

# Plot number of cells in each neighbourhood
nhood_size = np.array(mdata["rna"].obsm["nhoods"].sum(0)).ravel()
plt.hist(nhood_size, bins=100)
plt.xlabel("# cells in nhood")
plt.ylabel("# nhoods")
plt.title("Distribution of Milo neighborhood size")
plt.savefig(res_dir+"nhood_size.png")
plt.show()

# -----------------------------------------------------------------------------
# Count cells in neighbourhoods
# -----------------------------------------------------------------------------

mdata = milo.count_nhoods(mdata, sample_col="sample_ID")

# change Fast Progressor to Fast_Progressor
mdata['rna'].obs['classification'] = mdata['rna'].obs['classification'].str.replace(
    "Fast Progressor", "Fast_Progressor").astype("category")

mdata["rna"].obs["Status"] = mdata["rna"].obs["classification"].cat.reorder_categories(["stable", "Fast_Progressor"])


# -----------------------------------------------------------------------------
# 1. Predict eGFR classification ---------------------------------------------------
# -----------------------------------------------------------------------------

# Differential abundance testing with GLM --------------------------------------
milo.da_nhoods(mdata, design="~classification",
               model_contrasts="classificationFast_Progressor-classificationstable")

old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10, 10]
plt.subplot(1, 2, 1)
plt.hist(mdata["milo"].var.PValue, bins=50)
plt.xlabel("P-Vals")
plt.subplot(1, 2, 2)
plt.plot(mdata["milo"].var.logFC, -np.log10(mdata["milo"].var.SpatialFDR), ".")
plt.xlabel("log-Fold Change")
plt.ylabel("- log10(Spatial FDR)")
plt.tight_layout()
plt.suptitle("Milo: DA testing (~sample)", y=1.05)
plt.rcParams["figure.figsize"] = old_figsize
plt.savefig(res_dir+"da_volcano_plot.png")


# Plot umap with logFC ----------------------------
milo.build_nhood_graph(mdata)
plt.rcParams["figure.figsize"] = [16, 16]
milo.plot_nhood_graph(
    mdata,
    alpha=0.1,  # SpatialFDR level (1%)
    min_size=0.1,  # Size of smallest dot
    title="Neighborhood Differential Abundance LogFC (~sample)",
    return_fig=True
)
plt.savefig(res_dir+"nhood_umap.png")
plt.show()

# Visualize result by celltype -----
milo.annotate_nhoods(mdata, anno_col="main_class")

plt.hist(mdata["milo"].var["nhood_annotation_frac"], bins=30)
plt.xlabel("celltype fraction")
plt.suptitle("Concentration of celltypes in neighbourhoods (~sample)", y=0.97)
plt.savefig(res_dir+"nhood_annotation_frac.png")

mdata["milo"].var["nhood_annotation"] = mdata["milo"].var["nhood_annotation"].astype("str")

# Save mdata object
mdata.write_h5mu(res_dir+"milo_obj.h5mu")
