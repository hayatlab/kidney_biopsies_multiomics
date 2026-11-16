# Steps 5-7 of the scMEGA pipeline: TFs whose expression tracks their motif activity (panel c), peak-to-gene links,
# TF -> gene regulatory network. Sourced by 27_scmega_run.R (DC1-rank Trajectory) and 32_scmega_traj.R (ArchR
# trajectories). Expects in the calling environment: obj.pair (with $Trajectory in (0, 100], no NA), det (per gene:
# % of PT / aPT RNA cells with >= 1 count), and the output-path helpers o() and fg().

# (5) TFs whose expression tracks their motif activity (panel c) ----------------------------------------
trajMM  <- GetTrajectory(obj.pair, assay = "chromvar", trajectory.name = "Trajectory", slot = "data",
                         smoothWindow = 7, log2Norm = FALSE)
rownames(trajMM) <- obj.pair@assays$ATAC@motifs@motif.names
trajRNA <- GetTrajectory(obj.pair, assay = "RNA", trajectory.name = "Trajectory", slot = "data",
                         smoothWindow = 7, log2Norm = TRUE)
df.all <- GetCorrelation(trajMM, trajRNA)            # every motif whose name is an RNA gene, both signs
# Leon (2026-09-11): no TF expression filter, but flag TFs detected in < 5% of both PT and aPT cells (ambient level:
# the immune marker PTPRC is detected in 3-4% of PT / aPT cells).
add_det <- function(d) {
  ix <- match(d$tfs, det$gene)
  d$pct_PT <- det$pct_PT[ix]; d$pct_aPT <- det$pct_aPT[ix]
  d$low_detection_lt5pct <- pmax(d$pct_PT, d$pct_aPT) < 5
  d
}
df.all <- add_det(df.all)
write.csv(cbind(tf = rownames(df.all), df.all), o("tf_motif_expression_correlation_all", ".csv"), row.names = FALSE)
res <- SelectTFs(object = obj.pair, tf.assay = "chromvar", rna.assay = "RNA", atac.assay = "ATAC",
                 trajectory.name = "Trajectory", p.cutoff = 0.01, cor.cutoff = 0.3, return.heatmap = TRUE)
df.tf <- add_det(res$tfs)
# pseudotime bin (1-100) where the TF's smoothed motif activity peaks (scMEGA's time_point is only a heatmap row rank)
sm <- SummarizedExperiment::assay(trajMM)
df.tf$activity_peak_bin <- unname(apply(sm[df.tf$tfs, , drop = FALSE], 1, which.max))
write.csv(df.tf, o("selected_tfs", ".csv"), row.names = FALSE)
cat(sprintf("TFs with motif + expressed gene: %d | positively correlated (adj p < 0.01, r > 0.3): %d (%d detected in < 5%% of both PT and aPT cells)\n",
            nrow(df.all), nrow(df.tf), sum(df.tf$low_detection_lt5pct, na.rm = TRUE)))
for (ext in c(".pdf", ".png")) {
  if (ext == ".pdf") pdf(fg("panel_c_tf_heatmap", ext), width = 9, height = max(6, 0.16 * nrow(df.tf) + 2.5))
  else png(fg("panel_c_tf_heatmap", ext), width = 9, height = max(6, 0.16 * nrow(df.tf) + 2.5), units = "in", res = 200)
  draw(res$heatmap); dev.off()
}

# (6) peak-to-gene links along the trajectory -----------------------------------------------------------
res2 <- SelectGenes(object = obj.pair, atac.assay = "ATAC", rna.assay = "RNA", trajectory.name = "Trajectory",
                    genome = "hg38", labelTop1 = 0, labelTop2 = 0, return.heatmap = TRUE)
df.p2g <- res2$p2g; write.csv(df.p2g, o("peak_to_gene_links", ".csv"), row.names = FALSE)
cat(sprintf("peak-to-gene links (FDR < 1e-4, peak > 2 kb from TSS): %d links, %d genes, %d peaks\n",
            nrow(df.p2g), length(unique(df.p2g$gene)), length(unique(df.p2g$peak))))
if (nrow(df.p2g) > 0) for (ext in c(".pdf", ".png")) {
  if (ext == ".pdf") pdf(fg("p2g_heatmap", ext), width = 10, height = 8)
  else png(fg("p2g_heatmap", ext), width = 10, height = 8, units = "in", res = 200)
  draw(res2$heatmap); dev.off()
}

# (7) TF -> gene regulatory network ---------------------------------------------------------------------
if (nrow(df.p2g) > 0 && nrow(df.tf) > 0) {
  grn <- tryCatch({
    tfg <- GetTFGeneCorrelation(object = obj.pair, tf.use = df.tf$tfs, gene.use = unique(df.p2g$gene),
                                tf.assay = "chromvar", gene.assay = "RNA", trajectory.name = "Trajectory")
    mm <- obj.pair@assays$ATAC@motifs@data
    mn <- unlist(obj.pair@assays$ATAC@motifs@motif.names)          # motif ID -> TF name, mapped by ID
    stopifnot(all(colnames(mm) %in% names(mn)))
    colnames(mm) <- unname(mn[colnames(mm)])
    mm <- mm[unique(df.p2g$peak), intersect(unique(tfg$tf), colnames(mm)), drop = FALSE]
    GetGRN(motif.matching = mm, df.cor = tfg, df.p2g = df.p2g)
  }, error = function(e) { cat("GRN step failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(grn)) {
    # GetGRN left-joins motif hits to the links, so motif-carrying peaks without a linked gene return gene = NA,
    # and targets outside GetTFGeneCorrelation's top-variable genes return correlation = NA. Neither is an edge:
    # keep the raw table, report the count, drop them from the edge list.
    write.csv(grn, o("grn_raw", ".csv"), row.names = FALSE)
    ok <- !is.na(grn$gene) & !is.na(grn$correlation)
    cat(sprintf("GRN: %d rows from GetGRN; %d without a linked gene or a TF-gene correlation dropped\n", nrow(grn), sum(!ok)))
    grn <- grn[ok, ]
    write.csv(grn, o("grn_edges", ".csv"), row.names = FALSE)
    cat(sprintf("GRN: %d TF -> gene edges (%d with r > 0), %d TFs, %d targets\n", nrow(grn), sum(grn$correlation > 0),
                length(unique(grn$tf)), length(unique(grn$gene))))
  }
}
writeLines(capture.output(sessionInfo()), o("sessionInfo", ".txt"))
cat("done\n")
