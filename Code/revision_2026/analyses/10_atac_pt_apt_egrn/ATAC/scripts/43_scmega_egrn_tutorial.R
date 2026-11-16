#!/usr/bin/env Rscript
# Mina Shaigan's scMEGA eGRN views (email of 2026-09-16; her notebook Charlotte_project/scMEGA_eGRN.ipynb, which follows the
# scMEGA myofibroblast-GRN tutorial) on our paired object, as a script:
#   (1) SelectTFs with a stricter TF activity-expression correlation, cor.cutoff = 0.7 (Ivan: use this correlation to reduce
#       the TFs). SelectTFs keeps r > cutoff only -- positive; strongly negative TFs are counted here but never selected
#   (2) SelectGenes with Mina's link thresholds (r > 0.7, FDR < 1e-50) -> GetTFGeneCorrelation (GRNHeatmap) -> GetGRN;
#       network edges = TF-gene r > 0.7 (Mina; the tutorial uses 0.4)
#   (3) TF centrality in that network exactly as scMEGA::GRNPlot computes it: PageRank and betweenness (igraph, directed
#       TF -> gene, edge weight = r), each z-scored over the TFs; importance = distance from the minimum of both ->
#       importance plot, PageRank vs betweenness, eGRN (GRNPlot). Out-degree (number of targets) added for reading
#   (4) AddTargetAssay (Seurat AddModuleScore of each TF's targets) + the PseudotimePlot curves (motif activity, TF
#       expression, target expression; each z-scored over the 100 smoothed bins) for the TFs with importance >= 0.1
#       (Mina's rule) + HNF4A, ELF3 (her examples) + the NF-kB TFs in the network. Target scores on log-normalised RNA
#       (AddModuleScore expects it; the paired object keeps counts in its data slot, which Mina's run used) -- the counts
#       version is computed as well and compared. Curves are rebuilt from the same scMEGA calls PseudotimePlot makes and
#       checked against PseudotimePlot for one TF
# TRAJ=archr_dm (main-figure pseudotime; pairs without an ArchR value removed, as in 32) | dc1rank (the object's default
#   Trajectory column -- the axis Mina's plots used). SMOKE=1 -> smoke object. Thresholds: TF_COR, P2G_COR, P2G_FDR, GRN_COR.
# Outputs: scmega/egrn_tutorial/<TRAJ>[_smoke]/ (tables, figures/, meta.json). Human hg38, HGNC symbols, JASPAR2020, seed 42.
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
options(Seurat.object.assay.version = "v3", stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(scMEGA); library(ArchR); library(SummarizedExperiment)
  library(dplyr); library(igraph); library(ggraph); library(ggplot2); library(ggrepel); library(patchwork)
  library(ComplexHeatmap)
})
set.seed(42)
A    <- paste0(KB_ANALYSIS, "/aki_c5_exclusion_20261022/atac_rerun/ATAC")
SM   <- Sys.getenv("SMOKE") == "1"; TAG0 <- if (SM) "_smoke" else ""
TRAJ <- Sys.getenv("TRAJ", "archr_dm"); stopifnot(TRAJ %in% c("archr_dm", "dc1rank"))
num  <- function(k, d) as.numeric(Sys.getenv(k, as.character(d)))
TF_COR <- num("TF_COR", 0.7); TF_P <- 0.01
P2G_COR <- num("P2G_COR", 0.7); P2G_FDR <- num("P2G_FDR", 1e-50); GRN_COR <- num("GRN_COR", 0.7)
IMP_MIN <- 0.1                                        # Mina: TFs with importance >= 0.1 get a pseudotime plot
NFKB <- c("RELA", "RELB", "REL", "NFKB1", "NFKB2")
OUT <- file.path(A, "scmega", "egrn_tutorial", paste0(TRAJ, TAG0, Sys.getenv("OUT_SUFFIX"))); FIG <- file.path(OUT, "figures")
dir.create(file.path(FIG, "per_tf"), FALSE, TRUE)
o <- function(x) file.path(OUT, x)
save_fig <- function(name, w, h, p = NULL, draw_fun = NULL, dir = FIG) {
  for (ext in c("pdf", "png")) {
    f <- file.path(dir, paste0(name, ".", ext))
    # cairo devices: pdf() cannot draw the arrow / >= glyphs in the labels (smoke job 329976)
    if (ext == "pdf") cairo_pdf(f, width = w, height = h) else png(f, width = w, height = h, units = "in", res = 200, type = "cairo")
    if (is.null(draw_fun)) print(p) else draw_fun()
    dev.off()
  }
}
PAL <- ArchR::paletteContinuous(set = "blueYellow")
meta <- list(script = "43_scmega_egrn_tutorial.R", run_date = as.character(Sys.Date()), trajectory = TRAJ, smoke = SM,
             versions = list(scMEGA = as.character(packageVersion("scMEGA")), igraph = as.character(packageVersion("igraph")),
                             Seurat = as.character(packageVersion("Seurat"))),
             thresholds = list(SelectTFs_p = TF_P, SelectTFs_cor = TF_COR, SelectGenes_cor = P2G_COR,
                               SelectGenes_fdr = P2G_FDR, GRN_edge_cor = GRN_COR, importance_min = IMP_MIN))
cat(sprintf("scMEGA %s | trajectory %s | smoke %s | TF r > %.2f, links r > %.2f FDR < %g, edges r > %.2f\n",
            meta$versions$scMEGA, TRAJ, SM, TF_COR, P2G_COR, P2G_FDR, GRN_COR))

obj <- readRDS(file.path(A, "scmega", "outputs", paste0("obj_pair", TAG0, ".rds")))
if (TRAJ == "archr_dm") {                              # as in 32_scmega_traj.R
  tr <- read.csv(file.path(A, "tables", paste0("archr_trajectories", TAG0, ".csv")))
  stopifnot(all(colnames(obj) %in% tr$pair))
  v <- tr$traj_archr_dm[match(colnames(obj), tr$pair)]
  cat(sprintf("archr_dm: %d of %d pairs have a value; %d without removed\n", sum(!is.na(v)), length(v), sum(is.na(v))))
  obj$Trajectory <- v
  obj <- subset(obj, cells = colnames(obj)[!is.na(v)])
}
stopifnot(!anyNA(obj$Trajectory), all(obj$Trajectory > 0 & obj$Trajectory <= 100))
meta$pairs <- ncol(obj)
det <- read.csv(file.path(A, "scmega", "outputs", paste0("rna_gene_detection", TAG0, ".csv")))
add_det <- function(d) {
  ix <- match(d$tfs, det$gene)
  d$pct_PT <- det$pct_PT[ix]; d$pct_aPT <- det$pct_aPT[ix]
  d$low_detection_lt5pct <- pmax(d$pct_PT, d$pct_aPT) < 5
  d
}

# (1) TFs: motif activity vs expression along the pseudotime ----------------------------------------------------------
trajMM <- GetTrajectory(obj, assay = "chromvar", trajectory.name = "Trajectory", slot = "data", smoothWindow = 7,
                        log2Norm = FALSE)
rownames(trajMM) <- obj@assays$ATAC@motifs@motif.names
trajRNA <- GetTrajectory(obj, assay = "RNA", trajectory.name = "Trajectory", slot = "data", smoothWindow = 7, log2Norm = TRUE)
df.all <- add_det(GetCorrelation(trajMM, trajRNA))       # every motif whose name is an RNA gene, both signs
write.csv(df.all[order(-df.all$correlation), ], o("tf_motif_expression_correlation_all.csv"), row.names = FALSE)
cut_tab <- do.call(rbind, lapply(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), function(k) {
  pos <- df.all$adj_p < TF_P & df.all$correlation > k
  data.frame(cor_cutoff = k, tfs_tested = nrow(df.all), selected = sum(pos),
             selected_low_detection = sum(pos & df.all$low_detection_lt5pct, na.rm = TRUE),
             negative_below_minus_cutoff = sum(df.all$adj_p < TF_P & df.all$correlation < -k))
}))
write.csv(cut_tab, o("tf_counts_by_cor_cutoff.csv"), row.names = FALSE)
print(cut_tab)
res <- SelectTFs(object = obj, tf.assay = "chromvar", rna.assay = "RNA", atac.assay = "ATAC",
                 trajectory.name = "Trajectory", p.cutoff = TF_P, cor.cutoff = TF_COR, return.heatmap = TRUE)
df.cor <- add_det(res$tfs)
sm <- assay(trajMM)
df.cor$activity_peak_bin <- unname(apply(sm[df.cor$tfs, , drop = FALSE], 1, which.max))
write.csv(df.cor, o("selected_tfs.csv"), row.names = FALSE)
neg <- df.all[df.all$adj_p < TF_P & df.all$correlation < -TF_COR, ]
write.csv(neg[order(neg$correlation), ], o("tfs_negative_not_selected.csv"), row.names = FALSE)
meta$tfs <- list(tested = nrow(df.all), selected = nrow(df.cor),
                 selected_low_detection = sum(df.cor$low_detection_lt5pct, na.rm = TRUE),
                 negative_not_selected = nrow(neg), negative_names = neg$tfs[order(neg$correlation)])
cat(sprintf("(1) TFs tested %d | selected (adj p < %g, r > %.2f) %d, %d detected in < 5%% of PT and aPT | r < -%.2f (not selected): %d\n",
            nrow(df.all), TF_P, TF_COR, nrow(df.cor), meta$tfs$selected_low_detection, TF_COR, nrow(neg)))
# SelectTFs' own heatmap (CorrelationHeatmap) labels only 50 rows; redrawn from the same smoothed trajectories and the same
# call with every TF labelled
ht_all <- suppressMessages(CorrelationHeatmap(trajectory1 = trajMM[df.cor$tfs, ], trajectory2 = trajRNA[df.cor$tfs, ],
                                              name1 = "TF activity", name2 = "Gene expression",
                                              labelTop1 = nrow(df.cor), labelTop2 = nrow(df.cor)))
hh <- max(6, 0.13 * nrow(df.cor) + 2)
save_fig(sprintf("tf_heatmap_cor%02d", round(100 * TF_COR)), 10, hh,
         draw_fun = function() draw(ht_all, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 12), "mm")))

# (2) peak-to-gene links, TF-gene correlation, GRN --------------------------------------------------------------------
sel_genes <- function(hm) SelectGenes(object = obj, atac.assay = "ATAC", rna.assay = "RNA", trajectory.name = "Trajectory",
                                      cor.cutoff = P2G_COR, fdr.cutoff = P2G_FDR, genome = "hg38", labelTop1 = 0,
                                      labelTop2 = 0, return.heatmap = hm)
# scMEGA's heatmap step fails on an empty link table (smoke job 329975) -> links without the heatmap, reported
res2 <- tryCatch(sel_genes(TRUE), error = function(e) {
  cat("SelectGenes heatmap failed:", conditionMessage(e), "-> links without heatmap\n"); sel_genes(FALSE) })
df.p2g <- res2$p2g
write.csv(df.p2g, o("peak_to_gene_links.csv"), row.names = FALSE)
meta$links <- list(links = nrow(df.p2g), genes = length(unique(df.p2g$gene)), peaks = length(unique(df.p2g$peak)))
cat(sprintf("(2) links (r > %.2f, FDR < %g): %d links, %d genes, %d peaks\n", P2G_COR, P2G_FDR, nrow(df.p2g),
            meta$links$genes, meta$links$peaks))
if (nrow(df.p2g) == 0) stop("no peak-to-gene links at these thresholds")
if (!is.null(res2$heatmap)) save_fig("p2g_heatmap", 10, 8, draw_fun = function() draw(res2$heatmap))
tfg <- GetTFGeneCorrelation(object = obj, tf.use = df.cor$tfs, gene.use = unique(df.p2g$gene),
                            tf.assay = "chromvar", gene.assay = "RNA", trajectory.name = "Trajectory")
write.csv(tfg, o("tf_gene_correlation.csv"), row.names = FALSE)
save_fig("tf_gene_correlation_heatmap", 22, 8,
         draw_fun = function() draw(GRNHeatmap(tfg, tf.timepoint = df.cor$time_point)))
mm <- obj@assays$ATAC@motifs@data
mn <- unlist(obj@assays$ATAC@motifs@motif.names)          # motif ID -> TF name, mapped by ID (as in 27_downstream.R)
stopifnot(all(colnames(mm) %in% names(mn)))
colnames(mm) <- unname(mn[colnames(mm)])
mm <- mm[unique(df.p2g$peak), intersect(unique(tfg$tf), colnames(mm)), drop = FALSE]
grn <- GetGRN(motif.matching = mm, df.cor = tfg, df.p2g = df.p2g)
ok <- !is.na(grn$gene) & !is.na(grn$correlation)          # motif peaks without a linked gene / genes without a correlation
grn <- grn[ok, ]
write.csv(grn, o("grn_edges_all.csv"), row.names = FALSE)
df.grn2 <- grn %>% subset(correlation > GRN_COR) %>% dplyr::select(tf, gene, correlation) %>%
  dplyr::rename(weights = correlation)
write.csv(df.grn2, o("grn_edges.csv"), row.names = FALSE)
meta$grn <- list(rows_dropped_na = sum(!ok), edges_all = nrow(grn), edges_positive = sum(grn$correlation > 0),
                 edges = nrow(df.grn2), tfs = length(unique(df.grn2$tf)), targets = length(unique(df.grn2$gene)))
cat(sprintf("(2) GRN: %d edges (%d NA rows dropped); r > %.2f: %d edges, %d TFs, %d targets\n", nrow(grn), sum(!ok),
            GRN_COR, nrow(df.grn2), meta$grn$tfs, meta$grn$targets))

# (3) centrality (GRNPlot's formula) ----------------------------------------------------------------------------------
g <- igraph::graph_from_data_frame(df.grn2, directed = TRUE)
pr <- page_rank(g, weights = E(g)$weights)$vector
bt <- betweenness(g, weights = E(g)$weights, normalized = TRUE)
cent <- data.frame(tf = V(g)$name, pagerank_raw = unname(pr[V(g)$name]), betweenness_raw = unname(bt[V(g)$name])) %>%
  subset(tf %in% df.grn2$tf) %>%
  mutate(pagerank = scale(pagerank_raw)[, 1], betweenness = scale(betweenness_raw)[, 1])
cent$importance <- sqrt((cent$pagerank - min(cent$pagerank))^2 + (cent$betweenness - min(cent$betweenness))^2)
cent$n_targets <- as.integer(table(df.grn2$tf)[cent$tf])
cent$n_regulating_tfs <- vapply(cent$tf, function(t) sum(df.grn2$gene == t), 1L)   # in-degree (edges from other TFs)
ix <- match(cent$tf, df.cor$tfs)
cent <- cbind(cent, df.cor[ix, c("correlation", "time_point", "activity_peak_bin", "pct_PT", "pct_aPT", "low_detection_lt5pct")])
names(cent)[names(cent) == "correlation"] <- "motif_expression_r"
cent <- cent[order(-cent$importance, -cent$n_targets), ]; cent$importance_rank <- seq_len(nrow(cent))
write.csv(cent, o("tf_centrality.csv"), row.names = FALSE)
imp_tfs <- cent$tf[cent$importance >= IMP_MIN]
plot_tfs <- unique(c(imp_tfs, intersect(c("HNF4A", "ELF3"), cent$tf), intersect(NFKB, cent$tf)))
meta$centrality <- list(tfs = nrow(cent), importance_ge_min = length(imp_tfs), importance_ge_min_tfs = imp_tfs,
                        importance_zero = sum(cent$importance == 0),
                        top10_by_out_degree = head(cent$tf[order(-cent$n_targets)], 10),
                        spearman_importance_vs_out_degree = unname(cor(cent$importance, cent$n_targets, method = "spearman")))
cat(sprintf("(3) %d TFs in the network; importance >= %.1f: %d (%s)\n", nrow(cent), IMP_MIN, length(imp_tfs),
            paste(imp_tfs, collapse = ", ")))

co <- cent; co$tf <- factor(co$tf, levels = co$tf)
p_imp <- ggplot(co, aes(tf, importance)) +
  geom_hline(yintercept = IMP_MIN, linetype = 2, colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = activity_peak_bin, shape = low_detection_lt5pct), size = 2.2) +
  scale_colour_gradientn(colours = PAL, limits = c(1, 100), name = "Motif activity\npeak (bin)") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("≥ 5%", "< 5%"),
                     name = "TF detected in\nPT or aPT cells") +
  labs(x = NULL, y = "scMEGA importance",
       title = sprintf("TF importance in the eGRN (%d TFs; %s pseudotime)", nrow(co),
                       c(archr_dm = "ArchR diffusion-map", dc1rank = "DC1-rank")[[TRAJ]])) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), plot.title = element_text(size = 11))
save_fig("tf_importance", max(8, 0.12 * nrow(co) + 2.5), 4.5, p_imp)

p_pb <- ggplot(cent, aes(pagerank, betweenness)) +
  geom_point(aes(colour = activity_peak_bin, size = n_targets), alpha = 0.85) +
  geom_text_repel(data = subset(cent, importance >= IMP_MIN), aes(label = tf), size = 3.3, max.overlaps = Inf,
                  min.segment.length = 0, box.padding = 0.35, seed = 42) +
  annotate("text", x = min(cent$pagerank), y = min(cent$betweenness), hjust = -0.08, vjust = 1.6, size = 3.1,
           colour = "grey30", label = sprintf("%d TFs here: no edge from another TF in the network (importance 0)",
                                              sum(cent$importance == 0))) +
  scale_colour_gradientn(colours = PAL, limits = c(1, 100), name = "Motif activity\npeak (bin)") +
  scale_size_area(max_size = 7, name = "Targets\n(out-degree)") +
  labs(x = "PageRank (z-score over TFs)", y = "Betweenness (z-score over TFs)",
       title = sprintf("eGRN centrality: labelled = importance ≥ %.1f (%d of %d TFs)", IMP_MIN, length(imp_tfs), nrow(cent))) +
  theme_classic(base_size = 11) + theme(plot.title = element_text(size = 11))
save_fig("pagerank_betweenness", 8.5, 6.5, p_pb)

tfs.timepoint <- setNames(df.cor$time_point, df.cor$tfs)[order(df.cor$time_point)]
p_net <- tryCatch(GRNPlot(df.grn2, tfs.timepoint = tfs.timepoint, show.tf.labels = TRUE, seed = 42,
                          plot.importance = FALSE, min.importance = 2, remove.isolated = FALSE),
                  error = function(e) { cat("GRNPlot failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(p_net)) save_fig("egrn_network", 15, 15, p_net)
# same layout (same seed and call), labels only for the TFs with pseudotime curves (importance >= 0.1, HNF4A, ELF3, NF-kB)
p_net2 <- tryCatch(GRNPlot(df.grn2, tfs.timepoint = tfs.timepoint, show.tf.labels = FALSE, seed = 42,
                           plot.importance = FALSE, min.importance = 2, remove.isolated = FALSE,
                           genes.highlight = plot_tfs, cols.highlight = "#d95f02"),
                   error = function(e) { cat("GRNPlot (key TFs) failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(p_net2)) save_fig("egrn_network_key_tfs", 11, 11, p_net2)

# (4) targets, pseudotime curves ---------------------------------------------------------------------------------------
drop_tt <- function(x) { x@meta.data <- x@meta.data[, !grepl("^tf_target_", colnames(x@meta.data)), drop = FALSE]; x }
obj <- AddTargetAssay(object = obj, rna.assay = "RNA", df.grn = df.grn2)          # Mina's run: RNA data slot = counts
obj[["target_counts"]] <- obj[["target"]]; obj <- drop_tt(obj)
obj[["RNAlog"]] <- CreateAssayObject(counts = LayerData(obj, assay = "RNA", layer = "counts"))
obj <- NormalizeData(obj, assay = "RNAlog", normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj <- AddTargetAssay(object = obj, rna.assay = "RNAlog", df.grn = df.grn2)       # used for the figures
obj <- drop_tt(obj)

meta$plotted_tfs <- plot_tfs
meta$nfkb_in_network <- intersect(NFKB, cent$tf); meta$nfkb_selected <- intersect(NFKB, df.cor$tfs)
zlong <- function(se, rows, label) {                      # PseudotimePlot's transformation (row z-score over the bins)
  m <- assay(se); rownames(m) <- toupper(rownames(m)); m <- m[rows, , drop = FALSE]
  m <- t(scale(t(m)))
  data.frame(tf = rep(rownames(m), times = ncol(m)), pseudotime = rep(seq_len(ncol(m)), each = nrow(m)),
             value = as.vector(m), data = label)
}
trajT  <- GetTrajectory(obj, assay = "target", trajectory.name = "Trajectory", slot = "data", smoothWindow = 7, log2Norm = FALSE)
trajTc <- GetTrajectory(obj, assay = "target_counts", trajectory.name = "Trajectory", slot = "data", smoothWindow = 7,
                        log2Norm = FALSE)
all_t <- cent$tf
curves_all <- rbind(zlong(trajMM, all_t, "TF activity"), zlong(trajRNA, all_t, "TF expression"),
                    zlong(trajT, all_t, "Targets expression"), zlong(trajTc, all_t, "Targets expression (counts)"))
write.csv(curves_all, gzfile(o("pseudotime_curves.csv.gz")), row.names = FALSE)
# check against scMEGA's PseudotimePlot for the first plotted TF
chk <- PseudotimePlot(object = obj, tf.use = plot_tfs[1])$data
mine <- subset(curves_all, tf == plot_tfs[1] & data != "Targets expression (counts)")
key <- function(d) paste(d$data, d$pseudotime)
meta$pseudotimeplot_check <- list(tf = plot_tfs[1], rows = nrow(chk),
                                  max_abs_diff = max(abs(chk$value - mine$value[match(key(chk), key(mine))])))
cat(sprintf("(4) PseudotimePlot check (%s): max |diff| = %.2e over %d values\n", plot_tfs[1],
            meta$pseudotimeplot_check$max_abs_diff, nrow(chk)))
wide <- tidyr::pivot_wider(curves_all, names_from = data, values_from = value)
per_tf <- wide %>% group_by(tf) %>% summarise(
  r_activity_expression = cor(`TF activity`, `TF expression`),
  r_activity_targets = cor(`TF activity`, `Targets expression`),
  r_expression_targets = cor(`TF expression`, `Targets expression`),
  r_targets_log_vs_counts = cor(`Targets expression`, `Targets expression (counts)`), .groups = "drop")
per_tf <- per_tf[match(cent$tf, per_tf$tf), ]
write.csv(per_tf, o("pseudotime_curve_correlations.csv"), row.names = FALSE)
meta$targets_log_vs_counts <- list(median_r = median(per_tf$r_targets_log_vs_counts),
                                   min_r = min(per_tf$r_targets_log_vs_counts),
                                   min_tf = per_tf$tf[which.min(per_tf$r_targets_log_vs_counts)],
                                   n_below_0.9 = sum(per_tf$r_targets_log_vs_counts < 0.9))
cat(sprintf("(4) targets curve, log-normalised vs counts: median r %.3f, min %.3f (%s), %d TFs < 0.9\n",
            meta$targets_log_vs_counts$median_r, meta$targets_log_vs_counts$min_r, meta$targets_log_vs_counts$min_tf,
            meta$targets_log_vs_counts$n_below_0.9))

COLS <- c(`TF activity` = "#D55E00", `TF expression` = "#0072B2", `Targets expression` = "#009E73")
lab <- setNames(sprintf("%s (%d targets%s)", cent$tf, cent$n_targets, ifelse(cent$low_detection_lt5pct, "; < 5% detected", "")),
                cent$tf)
cv <- subset(curves_all, tf %in% plot_tfs & data %in% names(COLS))
cv$data <- factor(cv$data, levels = names(COLS))
cv$panel <- factor(lab[cv$tf], levels = lab[plot_tfs])
curve_plot <- function(d, ncol = NULL) {
  ggplot(d, aes(pseudotime, value, colour = data)) +
    geom_point(size = 0.35, alpha = 0.25) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_colour_manual(values = COLS, name = NULL) +
    labs(x = "Pseudotime (PT → aPT)", y = "z-score over pseudotime bins") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9))
}
nc <- 4; nr <- ceiling(length(plot_tfs) / nc)
save_fig("tf_pseudotime_curves", 3.1 * nc, 2.5 * nr + 0.8,
         curve_plot(cv) + facet_wrap(~ panel, ncol = nc, scales = "free_y"))
for (t in plot_tfs) {
  save_fig(paste0("scMEGA_TF_", t), 5, 4, curve_plot(subset(cv, tf == t)) + ggtitle(lab[[t]]) +
             theme(plot.title = element_text(size = 10, face = "bold")), dir = file.path(FIG, "per_tf"))
}

jsonlite::write_json(meta, o("meta.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), o("sessionInfo.txt"))
cat("done\n")
