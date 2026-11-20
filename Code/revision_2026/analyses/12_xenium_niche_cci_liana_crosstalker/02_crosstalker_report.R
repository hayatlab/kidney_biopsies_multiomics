#!/usr/bin/env Rscript
# Step 2 - CrossTalkeR single-condition report from the LIANA result of one niche.
#
# Analysis by Sidrah Maryam (`cross_down.R` in the run folders, 22 Jul 2026), with the run as a parameter and the paths
# from config/; the computation is unchanged:
#   interactions with CellPhoneDB permutation p <= 0.05 are kept; LR score = `lr.mean` (mean of ligand and receptor
#   expression); `_` in complex names becomes `*`; CrossTalkeR::generate_report(threshold = 0) builds the cell-cell graph
#   (edge weight = sum of the LR scores of a cell-type pair; edge `inter` = number of interactions of the pair, rescaled
#   to 0.1-1.1) and the gene-cell graph, and writes LR_data.Rds, LR_data_final.Rds and Single_<report>.
# The network of Fig. 6d is the cell-cell interaction plot of this report.
# Usage: Rscript 02_crosstalker_report.R <niche20|niche1>          Environment: R with liana, CrossTalkeR (see README)
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
source(file.path(KB_REPO, "analyses", "12_xenium_niche_cci_liana_crosstalker", "cci_runs.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(CrossTalkeR)
})

run <- cci_run(commandArgs(trailingOnly = TRUE)[1])
liana_test <- readRDS(run$liana)

gud <- liana_test %>%
  filter(pvalue <= 0.05) %>%
  mutate(gene_A = gsub("_", "*", ligand.complex)) %>%
  mutate(gene_B = gsub("_", "*", receptor.complex)) %>%
  mutate(source = gsub("_", "", source)) %>%
  mutate(target = gsub("_", "", target)) %>%
  mutate(type_gene_A = rep("Ligand", length(.data$ligand))) %>%
  mutate(type_gene_B = rep("Receptor", length(.data$receptor)))

# list of conditions for CrossTalkeR (one condition = single-condition analysis)
l <- list()
l[[run$label]] <- gud

data_dis <- generate_report(lrpaths = l,
                            genes = NULL,
                            out_path = paste0(run$out, "/"),
                            threshold = 0,
                            out_file = run$report,
                            output_fmt = "html_document",
                            report = TRUE,
                            sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B", "lr.mean"))
