#!/usr/bin/env Rscript
# Fig. 6d — cell-cell interaction network of glomerular niche 20 (Xenium): 6 cell types, 36 directed edges.
# Panel        : Main Figures.pdf, Figure 6, panel d
# Reads        : $KB_RESULTS/xenium/sidrah/whole_data_onlyglom_20/CCI_CrossTalkeR/LR_data_final.Rds
#                  (analyses/12_xenium_niche_cci_liana_crosstalker/01_liana_cellphonedb.R + 02_crosstalker_report.R, run niche20;
#                   published run of 22 Jul 2026 by Sidrah Maryam)
# Writes       : $KB_FIG_OUT/Fig6/Fig6d_cci_network_niche20_edges.csv (always: the plotted edge table)
#                $KB_FIG_OUT/Fig6/Fig6d_cci_network_niche20.pdf (only if the CrossTalkeR package is installed)
# Plotting code: CrossTalkeR::plot_cci on the cell-cell graph of the report object. Edge colour = weight (sum of the LR
#                scores of the cell-type pair, 0 to 484); edge width = `inter` (number of significant interactions of the
#                pair, rescaled to 0.1-1.1; legend "Percentage of the interactions").
#                The published panel is the cell-cell interaction plot of the CrossTalkeR report
#                (Single_overall_glom04_niche20.html); the arguments of plot_cci below are those needed to draw the same
#                graph and may differ in cosmetic settings from the report template.
# Illustrator  : node labels, legends ("% of the interactions", "Weight") and title "CCI: Glom niche 20" re-set
# Environment  : R with igraph (edge table: conda env de_env); the plot needs CrossTalkeR
# Verified     : 2026-11-20, edge table equal to the printed legends (weight maximum 484.29 = PEC -> PEC, `inter`
#                0.1 to 1.1, 36 edges) and recomputed from the interaction table (weight = sum of LR scores, inter =
#                rescaled count). The plot itself was NOT re-drawn: CrossTalkeR is not installed in the environments of
#                this repository's author.
# Writes in this script: write.csv / pdf below.

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
source(file.path(KB_REPO, "analyses", "12_xenium_niche_cci_liana_crosstalker", "cci_runs.R"))
suppressPackageStartupMessages(library(igraph))

run <- cci_run("niche20")
x <- readRDS(file.path(run$out, "LR_data_final.Rds"))
g <- attr(x, "graphs")[[run$label]]
tb <- as.data.frame(attr(x, "tables")[[run$label]])

# the plotted edge table, checked against the interaction table
e <- igraph::as_data_frame(g, "edges")
e$inter <- as.numeric(e$inter)                       # stored as a one-dimensional table
agg <- aggregate(LRScore ~ source + target, tb, sum)
n <- aggregate(LRScore ~ source + target, tb, length); names(n)[3] <- "n_interactions"
chk <- merge(merge(e, agg, by.x = c("from", "to"), by.y = c("source", "target")), n, by.x = c("from", "to"), by.y = c("source", "target"))
stopifnot(nrow(chk) == 36,
          isTRUE(all.equal(chk$weight, chk$LRScore.y)),
          isTRUE(all.equal(chk$inter, 0.1 + (chk$n_interactions - min(chk$n_interactions)) / diff(range(chk$n_interactions)))),
          abs(max(e$weight) - 484.29) < 0.01, isTRUE(all.equal(range(e$inter), c(0.1, 1.1))))
out <- data.frame(sender = chk$from, receiver = chk$to, n_interactions = chk$n_interactions,
                  pct_interactions_scaled = chk$inter, weight = chk$weight)
out <- out[order(-out$weight), ]
write.csv(out, fig_out("Fig6", "Fig6d_cci_network_niche20_edges.csv"), row.names = FALSE)

if (requireNamespace("CrossTalkeR", quietly = TRUE)) {
  pdf(fig_out("Fig6", "Fig6d_cci_network_niche20.pdf"), width = 7, height = 7)
  CrossTalkeR::plot_cci(graph = g, plt_name = "CCI: Glom niche 20",
                        colors = attr(x, "colors"), coords = attr(x, "coords")[igraph::V(g)$name, ],
                        emax = NULL, leg = TRUE, low = 0, high = 0, ignore_alpha = FALSE, log = FALSE,
                        efactor = 8, vfactor = 12)
  dev.off()
} else {
  message("CrossTalkeR is not installed: edge table written, plot skipped")
}
