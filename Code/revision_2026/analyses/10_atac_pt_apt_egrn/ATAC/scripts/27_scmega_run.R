#!/usr/bin/env Rscript
# scMEGA on the GLUE co-embedding, PT -> aPT (dPT branch excluded).
# Inputs (scripts/25_scmega_export.py) in scmega/inputs{,_smoke}/: rna_counts.h5 (RNA PT + aPT; ambient-
#   corrected counts), atac_peaks.h5 (ATAC pt_C3 + pt_C6; peaks in >= 1% of cells; hg38), glue_latent.csv,
#   metadata.csv.
# Trajectory: pt_traj = 100 x rank(pseudotime_joint) within the PT -> aPT set (tables/trajectory_pt_traj.csv),
#   averaged over the two members of each pair and re-ranked across pairs onto 0-100.
# Normalisation: GetTrajectory() averages the 'data' slot per pseudotime bin and depth-normalises + log2s it
#   itself (ArchR convention), so the paired object keeps data = counts (v3 assays) -- no TF-IDF/NormalizeData.
# Motifs: JASPAR2020 CORE vertebrates (latest versions); names upper-cased (human-first) and de-duplicated so
#   they match HGNC symbols; chromVAR with BSgenome.Hsapiens.UCSC.hg38.
# Correlations are across 100 smoothed pseudotime bins -> autocorrelated; treat p-values as descriptive.
# Human hg38 | HGNC symbols | seed 42 | SMOKE=1 -> inputs_smoke, outputs tagged _smoke
# Steps 5-7 live in 27_downstream.R, shared with 32_scmega_traj.R (re-runs them on the ArchR trajectories).
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
options(Seurat.object.assay.version = "v3", future.globals.maxSize = 200 * 1024^3, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(scMEGA); library(ArchR)
  library(chromVAR); library(motifmatchr); library(JASPAR2020); library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38); library(EnsDb.Hsapiens.v86)
  library(Matrix); library(dplyr); library(ggplot2); library(ComplexHeatmap); library(BiocParallel)
})
set.seed(42)
A   <- paste0(KB_ANALYSIS, "/aki_c5_exclusion_20261022/atac_rerun/ATAC")
SM  <- Sys.getenv("SMOKE") == "1"; TAG <- if (SM) "_smoke" else ""
IN  <- file.path(A, "scmega", paste0("inputs", TAG))
OUT <- file.path(A, "scmega", "outputs")
FIG <- file.path(A, "scmega", "figures_dc1rank")   # DC1-rank set; main figure = ArchR diffusion-map set (Leon, 2026-09-14)
dir.create(OUT, FALSE, TRUE); dir.create(FIG, FALSE, TRUE)
o  <- function(x, ext) file.path(OUT, paste0(x, TAG, ext))
fg <- function(x, ext) file.path(FIG, paste0(x, TAG, ext))
vs <- function(p) as.character(packageVersion(p))
cat("Seurat", vs("Seurat"), "| Signac", vs("Signac"), "| scMEGA", vs("scMEGA"), "| ArchR", vs("ArchR"),
    "| chromVAR", vs("chromVAR"), "| JASPAR2020", vs("JASPAR2020"), "| smoke", SM, "\n")
register(MulticoreParam(workers = 8, progressbar = FALSE))

meta <- read.csv(file.path(IN, "metadata.csv")); rownames(meta) <- meta$cellid
ptt  <- read.csv(file.path(A, "tables", "trajectory_pt_traj.csv"))
meta$pt_traj <- 100 * ptt$pt_traj[match(meta$cellid, ptt$cellid)]
stopifnot(!anyNA(meta$pt_traj))
glue <- as.matrix(read.csv(file.path(IN, "glue_latent.csv"), row.names = 1, check.names = FALSE))
colnames(glue) <- paste0("GLUE_", seq_len(ncol(glue)))

# (1) RNA and ATAC objects ------------------------------------------------------------------------------
rna <- Read10X_h5(file.path(IN, "rna_counts.h5"))
obj.rna <- CreateSeuratObject(rna, assay = "RNA", meta.data = meta[colnames(rna), ])
atac <- Read10X_h5(file.path(IN, "atac_peaks.h5"))
# genome = "hg38" (a string) makes Signac/GenomeInfoDb download chromInfo.txt.gz from hgdownload.cse.ucsc.edu,
# which fails (job 329589: "cannot open URL"; curl from the login node: the SSL certificate does not match that
# host name). Pass the Seqinfo shipped with the local BSgenome package instead -- no download.
obj.atac <- CreateSeuratObject(CreateChromatinAssay(atac, sep = c("-", "-"),
                                                    genome = seqinfo(BSgenome.Hsapiens.UCSC.hg38),
                                                    min.cells = 0, min.features = 0),
                               assay = "ATAC", meta.data = meta[colnames(atac), ])
rm(rna, atac); invisible(gc())
cat(sprintf("RNA %d cells x %d genes | ATAC %d cells x %d peaks\n",
            ncol(obj.rna), nrow(obj.rna), ncol(obj.atac), nrow(obj.atac)))

# (2) co-embedding object: all cells, GLUE latent as the reduction, 'tech' for pairing -----------------
cells <- c(colnames(obj.rna), colnames(obj.atac))
dummy <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(1, length(cells)), dimnames = list("dummy", cells))
coembed <- CreateSeuratObject(dummy, assay = "RNA", min.cells = 0, min.features = 0, meta.data = meta[cells, ])
coembed$tech <- ifelse(coembed$modality == "RNA", "RNA", "ATAC")
coembed[["glue"]] <- CreateDimReducObject(embeddings = glue[cells, ], key = "GLUE_", assay = "RNA")

# (3) pair cells on the GLUE latent; paired object ------------------------------------------------------
df.pair <- PairCells(object = coembed, reduction = "glue", pair.by = "tech", ident1 = "ATAC", ident2 = "RNA",
                     assay = "RNA", pair.mode = "greedy", seed = 42)
write.csv(df.pair, o("pairs", ".csv"), row.names = FALSE)
pq <- data.frame(atac_label = meta[df.pair$ATAC, "label"], rna_label = meta[df.pair$RNA, "label"],
                 pt_atac = meta[df.pair$ATAC, "pt_traj"], pt_rna = meta[df.pair$RNA, "pt_traj"])
cat(sprintf("pairs %d | distinct ATAC %d of %d | distinct RNA %d of %d | Spearman(pt ATAC member, pt RNA member) %.3f\n",
            nrow(df.pair), length(unique(df.pair$ATAC)), ncol(obj.atac), length(unique(df.pair$RNA)), ncol(obj.rna),
            cor(pq$pt_atac, pq$pt_rna, method = "spearman")))
print(table(ATAC = pq$atac_label, RNA = pq$rna_label))
write.csv(pq, o("pairs_qc", ".csv"), row.names = FALSE)

obj.pair <- CreatePairedObject(df.pair = df.pair, obj.coembed = coembed, obj.rna = obj.rna, obj.atac = obj.atac,
                               rna.assay = "RNA", atac.assay = "ATAC")
rm(coembed); invisible(gc())
ix <- match(colnames(obj.pair), df.pair$cell_name)
obj.pair$rna_cellid  <- df.pair$RNA[ix];  obj.pair$atac_cellid <- df.pair$ATAC[ix]
obj.pair$atac_label  <- meta[obj.pair$atac_cellid, "label"]
m <- (meta[obj.pair$rna_cellid, "pt_traj"] + meta[obj.pair$atac_cellid, "pt_traj"]) / 2
# rank / n scaling -> values in (0, 100]. GetTrajectory() bins cells as (b - 1, b] from 0, so a value of exactly 0
# (the earlier (rank - 1) / (n - 1) scaling) would silently drop the earliest pair from every bin.
obj.pair$Trajectory <- 100 * rank(m, ties.method = "average") / length(m)
# UCSC-style names by hand: seqlevelsStyle<- "UCSC" needs the same chromInfo.txt.gz download and fails the same
# way (scripts/scmega_offline_check.R, job 329631). GetGRangesFromEnsDb keeps standard chromosomes only:
# 1-22, X, Y, MT -> chr1-22, chrX, chrY, chrM.
annot <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
sl <- seqlevels(annot)
annot <- renameSeqlevels(annot, setNames(ifelse(sl == "MT", "chrM", paste0("chr", sl)), sl))
stopifnot(all(seqlevels(annot) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38)))
genome(annot) <- "hg38"
Annotation(obj.pair[["ATAC"]]) <- annot
cat(sprintf("paired object: %d pairs x %d genes (RNA) / %d peaks (ATAC, >= 10 paired cells)\n",
            ncol(obj.pair), nrow(obj.pair[["RNA"]]), nrow(obj.pair[["ATAC"]])))

# (4) motifs + chromVAR ---------------------------------------------------------------------------------
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
nm  <- toupper(vapply(seq_along(pfm), function(i) name(pfm[[i]]), ""))
hum <- vapply(seq_along(pfm), function(i) any(grepl("Homo sapiens|9606", unlist(tags(pfm[[i]])$species))), TRUE)
ord <- order(!hum)                                   # human matrices first, so they win name ties
keep <- ord[!duplicated(nm[ord])]
# TFBSTools 1.40.0 exports the name() getter but no name<- setter (job 329652) -> set the slot directly
pl <- lapply(keep, function(i) { x <- pfm[[i]]; x@name <- nm[i]; x })
pfm2 <- do.call(PFMatrixList, pl); names(pfm2) <- names(pfm)[keep]
cat(sprintf("JASPAR2020 CORE vertebrates: %d motifs -> %d after upper-casing + de-duplicating names (%d human)\n",
            length(pfm), length(pfm2), sum(hum[keep])))
obj.pair <- AddMotifs(obj.pair, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm2, assay = "ATAC")
obj.pair <- RunChromVAR(obj.pair, genome = BSgenome.Hsapiens.UCSC.hg38, new.assay.name = "chromvar", assay = "ATAC")
# SelectTFs, GetTFGeneCorrelation and step 5 name the chromVAR rows BY POSITION from motif.names -> the order
# must match, or TF labels would be silently shuffled.
stopifnot(identical(rownames(obj.pair[["chromvar"]]), names(obj.pair@assays$ATAC@motifs@motif.names)))
# chromVAR z-scores are undefined (NA) in cells with (almost) no reads in a motif's peaks and its background peaks.
# One NA cell makes its pseudotime bin NA and smoothing spreads that over 7 bins; scMEGA's TrajectoryHeatmap then
# drops the motif while GetCorrelation (complete cases) still selects it -> GetTFGeneCorrelation fails (smoke job
# 329674: 5 pairs, ATAC median 113 counts). Exclude pairs whose ATAC member has any NA deviation -- reported, not imputed.
cv <- LayerData(obj.pair, assay = "chromvar", layer = "data")
na.cell <- Matrix::colSums(is.na(cv)) > 0
atac.n <- Matrix::colSums(LayerData(obj.pair, assay = "ATAC", layer = "counts"))
qq <- quantile(atac.n, c(0, .01, .05, .5))
cat(sprintf("ATAC counts per paired cell (exported peaks): min %.0f, 1%% %.0f, 5%% %.0f, median %.0f\n", qq[1], qq[2], qq[3], qq[4]))
cat(sprintf("chromVAR: %d NA z-scores (%.4f%%) in %d pairs x %d motifs -> excluding those %d pairs (%.2f%%; their ATAC counts median %s)\n",
            sum(is.na(cv)), 100 * mean(is.na(cv)), sum(na.cell), sum(Matrix::rowSums(is.na(cv)) > 0), sum(na.cell),
            100 * mean(na.cell), if (any(na.cell)) sprintf("%.0f", median(atac.n[na.cell])) else "n/a"))
if (any(na.cell)) {
  write.csv(data.frame(pair = colnames(obj.pair)[na.cell], rna_cellid = obj.pair$rna_cellid[na.cell],
                       atac_cellid = obj.pair$atac_cellid[na.cell], atac_label = obj.pair$atac_label[na.cell],
                       atac_counts = atac.n[na.cell], n_na_motifs = Matrix::colSums(is.na(cv))[na.cell]),
            o("excluded_pairs_chromvar_na", ".csv"), row.names = FALSE)
  obj.pair <- subset(obj.pair, cells = colnames(obj.pair)[!na.cell])
  obj.pair$Trajectory <- 100 * rank(obj.pair$Trajectory, ties.method = "average") / ncol(obj.pair)   # re-rank onto (0, 100]
}
stopifnot(!anyNA(LayerData(obj.pair, assay = "chromvar", layer = "data")))
rm(cv); invisible(gc())
saveRDS(obj.pair, o("obj_pair", ".rds"))

# Detection of each gene in the RNA input (all PT / aPT RNA cells, not only the paired ones): % of cells with >= 1 count.
rc <- LayerData(obj.rna, layer = "counts")
det <- data.frame(gene = rownames(rc),
                  pct_PT  = 100 * Matrix::rowMeans(rc[, obj.rna$label == "PT"] > 0),
                  pct_aPT = 100 * Matrix::rowMeans(rc[, obj.rna$label == "aPT"] > 0))
rm(rc); invisible(gc())
write.csv(det, o("rna_gene_detection", ".csv"), row.names = FALSE)

# (5)-(7) TF selection, peak-to-gene links, GRN -- shared with 32_scmega_traj.R
source(file.path(KB_REPO, "analyses", "10_atac_pt_apt_egrn", "ATAC", "scripts", "27_downstream.R"), local = TRUE)
