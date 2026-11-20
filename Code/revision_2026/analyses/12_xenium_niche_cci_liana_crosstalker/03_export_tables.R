#!/usr/bin/env Rscript
# Step 3 - export the CrossTalkeR results of both niches as CSV tables for the figure scripts.
#
# Written 20 Nov 2026 (Claude Fable 5.2, for Leon Frommer) so that the figure scripts and the source data do not need
# the CrossTalkeR package: the LRObj written by step 2 is read through its attributes.
# Per run (niche20, niche1) it writes to $KB_CCI_EXPORT:
#   <run>_lr_table.csv        one row per significant interaction (source, target, gene_A, gene_B, LRScore, expression,
#                             p value, ...); for niche 20 this equals `LR_Report.csv`, the table read by the notebook of Fig. 6e
#   <run>_graph_edges.csv     edges of the cell-cell graph (from, to, LRScore, inter, weight)
#   <run>_cells_per_type.csv  number of cells per cell type in the input object
# Usage: Rscript 03_export_tables.R          Environment: R with igraph and SeuratObject (conda env de_env)
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
source(file.path(KB_REPO, "analyses", "12_xenium_niche_cci_liana_crosstalker", "cci_runs.R"))
suppressPackageStartupMessages({
  library(igraph)
  library(SeuratObject)
})
dir.create(CCI_EXPORT, showWarnings = FALSE, recursive = TRUE)

for (n in names(CCI_RUNS)) {
  run <- cci_run(n)
  x <- readRDS(file.path(run$out, "LR_data_final.Rds"))
  tb <- as.data.frame(attr(x, "tables")[[run$label]])
  write.csv(tb, file.path(CCI_EXPORT, paste0(n, "_lr_table.csv")), row.names = FALSE)
  g <- attr(x, "graphs")[[run$label]]
  write.csv(igraph::as_data_frame(g, "edges"), file.path(CCI_EXPORT, paste0(n, "_graph_edges.csv")), row.names = FALSE)
  S <- readRDS(file.path(run$dir, run$object))
  ct <- as.data.frame(table(cell_type = S[[]]$level2_refined))
  write.csv(ct, file.path(CCI_EXPORT, paste0(n, "_cells_per_type.csv")), row.names = FALSE)
  message(n, ": ", nrow(tb), " interactions, ", igraph::ecount(g), " edges, ", ncol(S), " cells")
}
