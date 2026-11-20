#!/usr/bin/env Rscript
# Step 1 - ligand-receptor inference with LIANA (CellPhoneDB method) on the Xenium cells of one spatial niche.
#
# Analysis by Sidrah Maryam (`crosstalker.R` in the run folders, 22 Jul 2026). This file is that script with the two
# run folders as a parameter and the paths taken from config/; the computation is unchanged:
#   counts of all cells of the niche -> Seurat LogNormalize (scale factor 10,000) -> liana_wrap(method = "cellphonedb"),
#   cell types = `level2_refined`, all other LIANA settings at their defaults (in the published results cell types with
#   fewer than 5 cells are absent and the permutation p values have a resolution of 0.001).
# Differences from the script as found on disk:
#   - the copy in the niche 1 folder subsets by an undefined `cond`; the published niche 1 result was computed on all
#     cells of the niche (its ligand expression equals the mean over all 6,793 aPT cells), as done here;
#   - the result is written to <run>/CCI_CrossTalkeR/, where step 2 reads it;
#   - set.seed(42) is new: the published runs were unseeded, so permutation p values of a re-run differ slightly
#     (the LR score `lr.mean` does not depend on the permutations).
# Usage: Rscript 01_liana_cellphonedb.R <niche20|niche1>          Environment: R with Seurat, liana (see README)
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
source(file.path(KB_REPO, "analyses", "12_xenium_niche_cci_liana_crosstalker", "cci_runs.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(liana)
  library(tibble)
})

run <- cci_run(commandArgs(trailingOnly = TRUE)[1])
dir.create(run$out, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

S <- readRDS(file.path(run$dir, run$object))
Idents(S) <- S$level2_refined
data <- as.SingleCellExperiment(S)
sce_subset <- assays(data)[["counts"]]
metadata_sce <- tibble::tibble(
  cells = rownames(colData(data)),
  celltype = as.character(colData(data)$level2_refined)
)
message(run$name, ": ", ncol(sce_subset), " cells, ", nrow(sce_subset), " genes")
print(table(metadata_sce$celltype))

lia_seurat <- CreateSeuratObject(counts = sce_subset)
Idents(lia_seurat) <- metadata_sce$celltype
lia_seurat <- NormalizeData(lia_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

liana_test <- liana_wrap(lia_seurat, method = "cellphonedb")
saveRDS(liana_test, file = run$liana)
message("written: ", run$liana)
