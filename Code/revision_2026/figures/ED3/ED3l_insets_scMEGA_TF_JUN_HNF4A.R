#!/usr/bin/env Rscript
# Extended Data Fig. 3l (insets) — scMEGA pseudotime curves of JUN and HNF4A along the PT->aPT ArchR pseudotime:
# z-scored TF motif activity, TF expression and expression of the TF's eGRN targets over the 100 smoothed pseudotime
# bins (points) with loess fits (lines).
# Panel        : Extended Data Figures.pdf, Extended Data Figure 3, panel l (the JUN and HNF4A insets above the network)
# Reads        : $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm_paper/
#                  pseudotime_curves.csv.gz and tf_centrality.csv (analyses/10_atac_pt_apt_egrn/ATAC/scripts/
#                  43_scmega_egrn_tutorial.R with TRAJ=archr_dm OUT_SUFFIX=_paper P2G_FDR=1e-4, job 340411 of
#                  9 Nov 2026; wrapper analyses/10_atac_pt_apt_egrn/job_43_egrn_paper_p2g1e-4.sh). NB the insets use the
#                  loose peak-to-gene FDR (1e-4; JUN 517 targets), the network beside them FDR 1e-50 (JUN 51 targets)
# Writes       : $KB_FIG_OUT/ED3/ED3l_inset_scMEGA_TF_JUN.{pdf,png} and ED3l_inset_scMEGA_TF_HNF4A.{pdf,png}
# Plotting code: extracted verbatim from analyses/10_atac_pt_apt_egrn/ATAC/scripts/43_scmega_egrn_tutorial.R, section 4
#                (save_fig, COLS, lab, cv, curve_plot and the per-TF loop that writes figures/per_tf/scMEGA_TF_<TF>);
#                the curves are read from the table that script saved instead of being recomputed from the 1.2 GB
#                paired object (43 would redo the whole eGRN)
# Illustrator  : shrunk to inset size; title "<TF> (<n> targets)" replaced by the TF name inside the plot; axis titles
#                "z-score" / "pseudotime"; the ggplot legend replaced by a separate legend whose keys use ggplot's default
#                colours (red "Targets expression", green "TF activity", blue "TF expression"). The curves keep the
#                colours of this script (orange = TF activity, green = Targets expression, blue = TF expression), so the
#                printed legend swaps TF activity and Targets expression (JUN: TF activity ends highest, z 2.29 vs
#                2.05 for targets in bins 96-100; drawn orange, labelled "Targets expression" by the printed legend)
# Environment  : conda env scmega2 (R 4.3.3, ggplot2)
# Verified     : 2026-11-12, PIXEL-IDENTICAL for both insets (PNG byte-identical) vs
#                $KB_ANALYSIS/aki_c5_exclusion_20261022/atac_rerun/ATAC/scmega/egrn_tutorial/archr_dm_paper/figures/per_tf/scMEGA_TF_{JUN,HNF4A}.pdf
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressPackageStartupMessages(library(ggplot2))
RUN <- file.path(KB_ANALYSIS, "aki_c5_exclusion_20261022", "atac_rerun", "ATAC", "scmega", "egrn_tutorial", "archr_dm_paper")
curves_all <- read.csv(file.path(RUN, "pseudotime_curves.csv.gz"), stringsAsFactors = FALSE)  # written by 43 (zlong)
cent <- read.csv(file.path(RUN, "tf_centrality.csv"), stringsAsFactors = FALSE)                # written by 43 (section 3)
plot_tfs <- c("JUN", "HNF4A")      # the two insets (43 draws every TF in meta.json$plotted_tfs, 36 in this run)
stopifnot(all(plot_tfs %in% cent$tf), is.logical(cent$low_detection_lt5pct), is.numeric(curves_all$value))
FIG <- dirname(fig_out("ED3", "ED3l_inset_scMEGA_TF_JUN.pdf"))

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
for (t in plot_tfs) {
  save_fig(paste0("ED3l_inset_scMEGA_TF_", t), 5, 4, curve_plot(subset(cv, tf == t)) + ggtitle(lab[[t]]) +
             theme(plot.title = element_text(size = 10, face = "bold")))
}
# -----------------------------------------------------------------------------------------------------------------------
cat(sprintf("[ED3l insets] %s\n", paste(file.path(FIG, paste0("ED3l_inset_scMEGA_TF_", plot_tfs, ".pdf")), collapse = ", ")))
