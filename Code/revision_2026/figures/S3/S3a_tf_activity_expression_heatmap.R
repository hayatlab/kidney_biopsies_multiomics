#!/usr/bin/env Rscript
# Supplementary Fig. 3a — chromVAR TF motif activity (left) and TF expression (right) along the PT->aPT ArchR
# pseudotime for the 81 TFs selected by scMEGA::SelectTFs (motif activity vs expression r > 0.7, adj. p < 0.01):
# rows = TFs (ordered by scMEGA::CorrelationHeatmap), columns = 100 smoothed pseudotime bins, values row-scaled.
# Panel        : All Supp Figures.pdf, Supplementary Figure 3, panel a
# Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/selected_tfs.csv
#                  (the 81 selected TFs in SelectTFs order; analyses/10_atac_pt_apt_egrn/ATAC/scripts/
#                  43_scmega_egrn_tutorial.R, TRAJ=archr_dm, job 339468 of 26 Oct 2026)
#                $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm_paper/
#                  pseudotime_curves.csv.gz (rows "TF activity" / "TF expression": the smoothed trajectories, row
#                  z-scored; same script with OUT_SUFFIX=_paper P2G_FDR=1e-4, job 340411). The archr_dm run saves curves
#                  only for the 80 TFs that keep an edge in its network (not NR1H4); TF selection and these two
#                  trajectories are computed before the peak-to-gene step, so they do not depend on P2G_FDR: the two
#                  runs have identical selected_tfs.csv and bit-identical curves for the 80 shared TFs (checked
#                  2026-11-12), and the _paper table holds all 81
# Writes       : $KB_FIG_OUT/S3/S3a_tf_heatmap_cor70.{pdf,png}
# Plotting code: extracted verbatim from analyses/10_atac_pt_apt_egrn/ATAC/scripts/43_scmega_egrn_tutorial.R, section 1
#                (CorrelationHeatmap call, height rule and save_fig); the two trajectories are rebuilt as
#                SummarizedExperiments from the saved table instead of scMEGA::GetTrajectory on the 1.2 GB paired object.
#                Row scaling is idempotent and CorrelationHeatmap scales every row again, so the z-scored input gives the
#                same heatmap
# Illustrator  : the two colour keys ("TF activity", "Gene expression", -2 to 2) moved from below to above the
#                heatmaps; row-label font re-set; panel letter
# Environment  : conda env scmega2 (R 4.3.3, scMEGA 1.1.0, ArchR 1.0.3, ComplexHeatmap)
# Verified     : 2026-11-12, PIXEL-IDENTICAL (PNG byte-identical) vs
#                $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm/figures/tf_heatmap_cor70.pdf (= scmega/figures_paper/panel_c_tf_heatmap_archr_dm.pdf)
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
options(Seurat.object.assay.version = "v3", stringsAsFactors = FALSE)
suppressPackageStartupMessages({                       # the packages 43 attaches before drawing this heatmap
  library(Seurat); library(Signac); library(scMEGA); library(ArchR); library(SummarizedExperiment)
  library(ComplexHeatmap)
})
set.seed(42)
TF_COR <- 0.7                                           # 43's SelectTFs cut-off (file name tf_heatmap_cor70)
B <- file.path(KB_ANALYSIS, "aki_c5_exclusion_20261022", "atac_rerun", "ATAC", "scmega", "egrn_tutorial")
df.cor <- read.csv(file.path(B, "archr_dm", "selected_tfs.csv"))
cur <- read.csv(file.path(B, "archr_dm_paper", "pseudotime_curves.csv.gz"))
stopifnot(nrow(df.cor) == 81, !anyDuplicated(df.cor$tfs))
bins <- paste0("T.", 0:99, "_", 1:100)                  # GetTrajectory's column names (groupEvery = 1)
se_from <- function(label) {                            # one TF x 100-bin matrix per readout, rows in SelectTFs order
  d <- cur[cur$data == label & cur$tf %in% df.cor$tfs, ]
  m <- matrix(NA_real_, nrow(df.cor), 100, dimnames = list(df.cor$tfs, bins))
  m[cbind(match(d$tf, df.cor$tfs), d$pseudotime)] <- d$value
  stopifnot(!anyNA(m), nrow(d) == 100 * nrow(df.cor))
  SummarizedExperiment(assays = S4Vectors::SimpleList(smoothMat = m))
}
trajMM <- se_from("TF activity"); trajRNA <- se_from("TF expression")
FIG <- dirname(fig_out("S3", "S3a_tf_heatmap_cor70.pdf"))

# ---- verbatim from 43_scmega_egrn_tutorial.R ------------------------------------------------------------------------
save_fig <- function(name, w, h, p = NULL, draw_fun = NULL, dir = FIG) {
  for (ext in c("pdf", "png")) {
    f <- file.path(dir, paste0(name, ".", ext))
    # cairo devices: pdf() cannot draw the arrow / >= glyphs in the labels (smoke job 329976)
    if (ext == "pdf") cairo_pdf(f, width = w, height = h) else png(f, width = w, height = h, units = "in", res = 200, type = "cairo")
    if (is.null(draw_fun)) print(p) else draw_fun()
    dev.off()
  }
}
# SelectTFs' own heatmap (CorrelationHeatmap) labels only 50 rows; redrawn from the same smoothed trajectories and the same
# call with every TF labelled
ht_all <- suppressMessages(CorrelationHeatmap(trajectory1 = trajMM[df.cor$tfs, ], trajectory2 = trajRNA[df.cor$tfs, ],
                                              name1 = "TF activity", name2 = "Gene expression",
                                              labelTop1 = nrow(df.cor), labelTop2 = nrow(df.cor)))
hh <- max(6, 0.13 * nrow(df.cor) + 2)
save_fig(sprintf("S3a_tf_heatmap_cor%02d", round(100 * TF_COR)), 10, hh,
         draw_fun = function() draw(ht_all, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 12), "mm")))
# -----------------------------------------------------------------------------------------------------------------------
cat(sprintf("[S3a] %s\n", file.path(FIG, sprintf("S3a_tf_heatmap_cor%02d.pdf", round(100 * TF_COR)))))
