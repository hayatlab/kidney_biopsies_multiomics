#!/usr/bin/env Rscript
# MAST DE of crescent-niche PEC vs other-niche PEC.
# Input: a small h5ad of PEC cells only, with X = raw integer counts and obs columns
#   'grp'  (values "crescent" / "other") and a donor column (default 'donor').
# Output: markers CSV (avg_log2FC = effect size, p_val_adj = BH-FDR, pct.1/pct.2), ALL genes
#   tested (logfc.threshold=0, min.pct=0) so the full effect-size vector is available for the
#   cross-analysis concordance comparison. MAST with donor as latent variable (per group standard).
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment)
  library(Seurat); library(SeuratObject); library(MAST)
})
args <- commandArgs(trailingOnly = TRUE)
in_h5ad <- args[1]; out_csv <- args[2]
donor   <- ifelse(length(args) >= 3, args[3], "donor")
cat(sprintf("[de] in=%s\n[de] out=%s\n[de] donor=%s\n", in_h5ad, out_csv, donor))

sce <- readH5AD(in_h5ad)                               # .X exported as raw counts -> assay "X"
counts <- assay(sce, 1)
md <- as.data.frame(colData(sce))
obj <- CreateSeuratObject(counts = counts, meta.data = md)
obj <- NormalizeData(obj, verbose = FALSE)             # log-normalise for MAST
stopifnot("grp" %in% colnames(obj@meta.data))
Idents(obj) <- "grp"
cat("[de] group sizes:\n"); print(table(obj$grp))
n_donor_cres <- length(unique(md[[donor]][md$grp == "crescent"]))
cat(sprintf("[de] donor levels in crescent: %d | other: %d\n",
            n_donor_cres, length(unique(md[[donor]][md$grp == "other"]))))

# donor is only identifiable as a covariate if the crescent group spans >1 donor;
# here the crescent lesion is a single biopsy, so drop the (collinear) donor term.
use_latent <- if (donor != "none" && n_donor_cres > 1) donor else NULL
cat(sprintf("[de] latent.vars = %s\n", ifelse(is.null(use_latent), "NONE (crescent is single-donor)", use_latent)))
m <- FindMarkers(obj, ident.1 = "crescent", ident.2 = "other",
                 test.use = "MAST", latent.vars = use_latent,
                 logfc.threshold = 0, min.pct = 0, verbose = FALSE)
m$gene <- rownames(m)
write.csv(m, out_csv, row.names = FALSE)
cat(sprintf("[de] wrote %d genes -> %s\n", nrow(m), out_csv))
cat("[de] top up in crescent:\n")
print(head(m[order(-m$avg_log2FC), c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")], 12))
