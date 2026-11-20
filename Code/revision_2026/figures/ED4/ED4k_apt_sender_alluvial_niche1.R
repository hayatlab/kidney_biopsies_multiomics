#!/usr/bin/env Rscript
# Extended Data Fig. 4k — alluvial plot of ligand-receptor interactions with aPT as sender in niche 1 (Xenium):
# sender -> ligand -> receptor -> receiver, 16 flows coloured by LR score.
# Panel        : Extended Data Figures.pdf, Extended Data Fig. 4, panel k
# Reads        : $KB_RESULTS/xenium/sidrah/overall_glom_04_level2/CCI_CrossTalkeR/LR_data_final.Rds
#                  (analyses/12_xenium_niche_cci_liana_crosstalker/01_liana_cellphonedb.R + 02_crosstalker_report.R, run niche1;
#                   published run of 22 Jul 2026 by Sidrah Maryam)
# Writes       : $KB_FIG_OUT/ED4/ED4k_apt_sender_alluvial_niche1_flows.csv (always: the 16 plotted interactions)
#                $KB_FIG_OUT/ED4/ED4k_apt_sender_alluvial_niche1.pdf (only if the CrossTalkeR package is installed)
# Plotting code: CrossTalkeR::plot_sankey(ligand_cluster = "aPT", threshold = 16). With a threshold below the number of
#                rows, plot_sankey draws the threshold/2 highest- and the threshold/2 lowest-scoring interactions: here
#                the 8 highest (all APP) and the 8 lowest (EFNB2, BMP2, ADAM15) of the 1,753 significant aPT-sender
#                interactions. The call of the published plot (aPT_sender_niche1_5LR.pdf, 30 Jul 2026) was not saved;
#                the selection below was identified from the stratum sizes of that PDF and is asserted here.
# Illustrator  : title "Top 5 ligands in niche 1 with aPT as sender"; "|L" / "|R" removed from the labels; dPT relabelled aPT2
# Environment  : R (selection table: conda env de_env); the plot needs CrossTalkeR and ggalluvial
# Verified     : 2026-11-20, the 16 selected interactions reproduce the strata of the published PDF (ligands APP 8,
#                EFNB2 6, BMP2 1, ADAM15 1; receptors APLP2 6, EPHB2 3, CD74 2, EPHA4 2, EPHB4 1, ACVR2A*BMPR1B 1, ITGA5 1;
#                11 receivers). The plot itself was NOT re-drawn: CrossTalkeR is not installed in the environments of this
#                repository's author.
# Writes in this script: write.csv / ggsave below.

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
source(file.path(KB_REPO, "analyses", "12_xenium_niche_cci_liana_crosstalker", "cci_runs.R"))

run <- cci_run("niche1")
x <- readRDS(file.path(run$out, "LR_data_final.Rds"))
tb <- as.data.frame(attr(x, "tables")[[run$label]])
THRESHOLD <- 16

# the rows plot_sankey selects (slice_max / slice_min on LRScore, n = ceiling(threshold / 2), no ties in this table)
a <- tb[tb$type_gene_A == "Ligand" & tb$type_gene_B == "Receptor" & tb$source == "aPT", ]
a <- a[order(-a$LRScore), ]
stopifnot(nrow(a) == 1753, !any(duplicated(a$LRScore[c(1:9, (nrow(a) - 8):nrow(a))])))
sel <- rbind(cbind(selection = "8 highest LR scores", head(a, THRESHOLD / 2)),
             cbind(selection = "8 lowest LR scores", tail(a, THRESHOLD / 2)))
cnt <- function(v) { t <- table(v); t[order(-t, names(t))] }
stopifnot(identical(as.vector(cnt(sel$gene_A)), c(8L, 6L, 1L, 1L)), identical(names(cnt(sel$gene_A)), c("APP|L", "EFNB2|L", "ADAM15|L", "BMP2|L")),
          sum(sel$gene_B == "APLP2|R") == 6, sum(sel$gene_B == "EPHB2|R") == 3, sum(sel$gene_B == "CD74|R") == 2,
          sum(sel$gene_B == "EPHA4|R") == 2, length(unique(sel$gene_B)) == 7, length(unique(sel$target)) == 11)
write.csv(sel[, c("selection", "source", "gene_A", "gene_B", "target", "LRScore", "ligand.expr", "receptor.expr", "pvalue")],
          fig_out("ED4", "ED4k_apt_sender_alluvial_niche1_flows.csv"), row.names = FALSE)

if (requireNamespace("CrossTalkeR", quietly = TRUE) && requireNamespace("ggalluvial", quietly = TRUE)) {
  p <- CrossTalkeR::plot_sankey(lrobj_tbl = tb, target = NULL, ligand_cluster = "aPT", receptor_cluster = NULL,
                                plt_name = "aPT as sender, niche 1", threshold = THRESHOLD)
  ggplot2::ggsave(fig_out("ED4", "ED4k_apt_sender_alluvial_niche1.pdf"), p, width = 7, height = 6)
} else {
  message("CrossTalkeR / ggalluvial not installed: flow table written, plot skipped")
}
