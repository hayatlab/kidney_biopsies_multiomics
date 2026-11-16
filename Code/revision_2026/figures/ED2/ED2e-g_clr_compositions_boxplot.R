#!/usr/bin/env Rscript
# Extended Data Fig. 2e-g — CLR-transformed abundance of aPT, PT and immune cells (snRNA-seq level-2 annotation, one
# point per sample) by condition (e; reference = control), eGFR at biopsy (f) and IF/TA (g); boxplots with Welch
# t-test brackets (ggpubr::stat_compare_means, unadjusted p).
# Panel        : Extended Data Figures.pdf, Extended Data Figure 2, panels e, f and g
# Reads        : $KB_ANALYSIS/Charlotte_help/outputs/clr~level2~new_annotation.csv and
#                  $KB_ANALYSIS/Charlotte_help/outputs/patient_metadata_new_annotation.csv (columns sample_ID, disease,
#                  classification, egfr_at_time_of_biopsy, IFTA, proteinuria_around_biopsy; no patient identifiers)
#                  (analyses/04_snrna_overview_composition/01_build_clr_inputs.R, which reads the per-sample tables of
#                  03_overview_figures_v2.py and the corrected sample metadata)
# Writes       : $KB_FIG_OUT/ED2/ED2e-g_clr_compositions_boxplot/ext2e_clr_compositions_boxplot_level2.{pdf,png} and
#                ext2e_stats_level2.csv (the numbers behind the brackets: n, mean difference with 95% CI, Welch t, p, BH)
# Plotting code: analyses/04_snrna_overview_composition/02_ext2e_clr_compositions_boxplot.R (Leon's adaptation of the
#                saezlab MCFA script ext2e_clr_compositions_boxplot.R), sourced unchanged with LEVEL = level2;
#                EXT2E_FIGDIR / EXT2E_TABLEDIR (release patch, default = the original folders) point to $KB_FIG_OUT
# Illustrator  : facet strip boxes removed, "1.73m2" set as "1.73m²", panel letters e / f / g placed on the three
#                patchwork panels
# Environment  : conda env scmega2 (R 4.3.3; ggplot2, ggpubr, patchwork, dplyr, tidyr)
# Verified     : 2026-11-12, DIFFERENT (1.2% of pixels, all of them jittered points; text layer equal, statistics
#                table byte-identical) vs $KB_ANALYSIS/Charlotte_help/figures/ext2e_clr_compositions_boxplot_level2.pdf
#                (2026-09-25). geom_jitter was unseeded in the original run, so the horizontal positions of the points
#                cannot be reproduced; this entry script seeds the jitter (42) so that its own output is reproducible
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
OUT <- fig_out("ED2", "ED2e-g_clr_compositions_boxplot")   # a folder: the script writes two figures and a table
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(EXT2E_FIGDIR = OUT, EXT2E_TABLEDIR = OUT)
set.seed(42)                       # the original script draws geom_jitter without a seed
# the script takes its level from commandArgs(trailingOnly = TRUE); this entry script is run without arguments ->
# "level2", the published version. Its plots are drawn only by ggsave (no auto-printing under source()).
source(file.path(KB_REPO, "analyses", "04_snrna_overview_composition", "02_ext2e_clr_compositions_boxplot.R"),
       local = new.env())
cat(sprintf("[ED2e-g] %s\n", file.path(OUT, "ext2e_clr_compositions_boxplot_level2.pdf")))
