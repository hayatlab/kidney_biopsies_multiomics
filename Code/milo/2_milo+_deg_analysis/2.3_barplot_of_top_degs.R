#
# Created on Fri Jul 25 2025
#
# Copyright (c) 2025 Vivien Goepp
#

# --------------------------------------------------------------------------------------------
# This script analyzes subcluster-specific DEGs from kidney biopsy scRNA-seq data,
# performs gene set enrichment analysis (GSEA) using Reactome and Hallmark pathways,
# and visualizes top enriched pathways per subcluster.
#
# Datasets used/loaded:
# - Top DEGs per subcluster: top_10_DEGs_per_subcluster_by_pval.csv
# - All significant DEGs per subcluster: milo_signif_DEGs_subcluster.csv
# - Reactome gene sets: C2:CP:REACTOME (msigdbr)
# - Hallmark gene sets: H (msigdbr)
#
# Datasets/objects saved:
# - Barplots of top 10 DEGs per subcluster by adjusted p-value:
#     - barplot_top_10_DEGs_per_subcluster_by_pval.pdf
# - Barplots of top 10 Reactome pathways per subcluster (by NES):
#     - barplot_top_10_pathways_per_subcluster.pdf
# - Barplots of top 10 Hallmark pathways per subcluster (by NES):
#     - barplot_top_10_hallmark_pathways_per_subcluster.pdf
#
# Analyses performed:
# - Read and process top DEGs per subcluster; reorder factors for plotting
# - Generate horizontal barplots of top DEGs (log10 adjusted p-value) per subcluster
# - Prepare ranked gene lists per subcluster for GSEA
# - Run GSEA on Reactome pathways per subcluster
# - Extract top 10 Reactome pathways per subcluster by NES
# - Plot top Reactome pathways per subcluster with faceted barplots
# - Prepare Hallmark gene sets and run GSEA per subcluster
# - Extract top 10 Hallmark pathways per subcluster by NES
# - Plot top Hallmark pathways per subcluster with faceted barplots
#
# Visualization outputs:
# - Horizontal barplots of top DEGs per subcluster
# - Faceted barplots of top Reactome pathways per subcluster
# - Faceted barplots of top Hallmark pathways per subcluster
#
# The script integrates DEG identification with pathway-level functional interpretation,
# producing both DEG summary plots and pathway enrichment barplots per subcluster.
# --------------------------------------------------------------------------------------------


library(tidytext)
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(BiocManager)
library(dplyr)

res_dir = "/data/zd775156/kidney_biopsy/results/milo/"


# 1. Plot top 10 DEGs per subcluster by adjusted p-value ---------------------------

top_degs_per_subcluster_logFC = read.csv(
  paste0(res_dir, "top_10_DEGs_per_subcluster_by_pval.csv")
) %>% as_tibble()

top_degs_per_subcluster_logFC$main_subcluster %>% unique()

# Define order of main_subcluster
main_subcluster_order = 
  c("1.0: TAL", "1.1: spp1+ aTAL/dTL", "1.2: cTAL", 
    "2.0: CNT/CCD-PC", "3.0: TAL/PT", "3.1: PT/DCT", "4.0: DTL", 
    '4.2: aTAL/aDCT', "4.3: aTAL/aDCT/aPT", "11.1: ATL_DTL_2")

top_degs_per_subcluster_logFC = top_degs_per_subcluster_logFC %>%
  mutate(main_subcluster = factor(main_subcluster, levels = main_subcluster_order))

top_degs_per_subcluster_logFC %>% names
# Barplot of log10_adj_pval vs names

ggplot(
  top_degs_per_subcluster_logFC %>%
    mutate(names = reorder_within(names, log10_adj_pval, main_subcluster)),
  aes(x = names, y = log10_adj_pval)
) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  coord_flip() +
  labs(x = "", y = "log10 Adjusted P-value") +
  facet_wrap(~main_subcluster, scales = "free_y") +
  scale_x_reordered() +
  theme(strip.background = element_blank()) +
  ggtitle("Top 10 DEGs per Subcluster by Adjusted P-value")

ggsave(
  paste0(res_dir, "barplot_top_10_DEGs_per_subcluster_by_pval.pdf"),
  width = 14, height = 10
)

# 2. Compute GO enrichment and plot barplot of top enriched terms ---------------------------

# 2.0 Functions to prepare ranks and run GSEA
prepare_ranks <- function(df) {
  df %>%
    group_by(names) %>%
    slice_max(log10_adj_pval, n = 1) %>%
    ungroup() %>%
    arrange(desc(log10_adj_pval)) %>%
    {
      setNames(.$log10_adj_pval, .$names)
    }
}

# 4. Function to run GSEA
run_gsea <- function(ranks, pathways) {
  GSEA(
    ranks,
    TERM2GENE = pathways,
    scoreType = "pos", # since all values > 0
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1
  )
}

# 2.1 Load data ----------------

all_DEGs <- read.csv(paste0(res_dir, "milo_signif_DEGs_subcluster.csv")) %>% as_tibble()


# 1. Prepare Reactome gene sets for human
reactome_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
) %>%
  select(gs_name, gene_symbol)

# 2. Split your DEG table by subcluster
subclusters_list <- all_DEGs %>%
  split(.$main_subcluster)

subclusters_list <- subclusters_list[names(subclusters_list) %in% main_subcluster_order]

# 5. Loop over subclusters
gsea_results_reactome <- map(subclusters_list, function(df_sub) {
  ranks <- prepare_ranks(df_sub)
  gsea <- run_gsea(ranks, reactome_sets)
  # convert to data.frame and add subcluster info
  as.data.frame(gsea) %>% mutate(subcluster = df_sub$main_subcluster[1])
})
gsea_results_reactome_df <- bind_rows(gsea_results_reactome)

# 6. Take top 10 pathways per subcluster (by NES)
top_pathways_per_sub <- gsea_results_reactome_df %>%
  # remove rows with NA subcluster
  filter(!is.na(subcluster)) %>%
  group_by(subcluster) %>%
  slice_max(order_by = NES, n = 10) %>%
  ungroup()

# 7. Plot REACTOME pathways
main_subcluster_order_2 = c(
  "1.0: TAL", "1.1: spp1+ aTAL/dTL", "1.2: cTAL",
  "2.0: CNT/CCD-PC", "2.2: DCT/CNT/TAL/PT", "3.0: TAL/PT", "3.1: PT/DCT",
  "4.0: DTL", "4.1: aPT", "4.2: aTAL/aDCT", "4.3: aTAL/aDCT/aPT", "11.1: ATL_DTL_2"
)

top_pathways_per_sub <- top_pathways_per_sub %>%
  mutate(subcluster = factor(subcluster, levels = main_subcluster_order_2)) %>%
  group_by(subcluster) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = as.character(Description)) %>%
  # Remove starting "REACTOME_" and replace "_" with " "
  mutate(Description = gsub("^REACTOME_", "", Description)) %>%
  mutate(Description = gsub("_", " ", Description)) %>%
  # break long descriptions into multiple lines
  mutate(Description = str_wrap(Description, width = 20))

top_pathways_per_sub <- top_pathways_per_sub %>%
  mutate(Description = reorder_within(Description, NES, subcluster))

ggplot(top_pathways_per_sub, aes(y = Description, x = NES)) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  labs(y = "", x = "Normalized Enrichment Score (NES)") +
  # theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.background = element_blank()) +
  facet_wrap(~subcluster, scales = "free_y", ncol = 4) +
  scale_y_reordered() +   # <- removes the "___subcluster" suffix in labels
  # scale_y_reordered(labels = function(x) str_wrap(x, width = 30)) +
  # scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  ggtitle("Top Reactome Pathways per Subcluster (GSEA)")

ggsave(
  paste0(res_dir, "barplot_top_10_pathways_per_subcluster.pdf"),
  plot = last_plot(),
  width = 22, height = 16
)


# ---------------------------------------------------------------------------
# 3. Pathway activity of Hallmark gene sets ---------------------------
# ---------------------------------------------------------------------------

# 3. Pathway activity of Hallmark gene sets ---------------------------

# 3.1. Prepare Hallmark gene sets for human
hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  select(gs_name, gene_symbol)

# 3.2 Compute GSEA for Hallmark pathways per subcluster
gsea_results_hallmark <- map(subclusters_list, function(df_sub) {
  ranks <- prepare_ranks(df_sub)
  # Only run GSEA if there are valid gene symbols
  if (length(ranks) > 0 && all(!is.na(names(ranks)))) {
    gsea <- run_gsea(ranks, hallmark_sets)
    as.data.frame(gsea) %>% mutate(subcluster = df_sub$main_subcluster[1])
  } else {
    NULL
  }
}) %>% compact()
gsea_results_hallmark <- bind_rows(gsea_results_hallmark)

# 3.3 Plot gsea_results_hallmark
# Prepare for plotting: take top 10 Hallmark pathways per subcluster by NES
top_hallmark_per_sub <- gsea_results_hallmark %>%
  filter(!is.na(subcluster)) %>%
  group_by(subcluster) %>%
  slice_max(order_by = NES, n = 10) %>%
  ungroup() %>%
  mutate(Description = as.character(Description)) %>%
  # Remove starting "HALLMARK_" and replace "_" with " "
  mutate(Description = gsub("^HALLMARK_", "", Description)) %>%
  mutate(Description = gsub("_", " ", Description)) %>%
  mutate(Description = str_wrap(Description, width = 20)) %>%
  mutate(subcluster = factor(subcluster, levels = main_subcluster_order_2)) %>%
  mutate(Description = reorder_within(Description, NES, subcluster))

# Plot
ggplot(top_hallmark_per_sub, aes(y = Description, x = NES)) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  labs(y = "", x = "Normalized Enrichment Score (NES)") +
  theme(axis.text.y = element_text(size = 12)) +
  theme(strip.background = element_blank()) +
  facet_wrap(~subcluster, scales = "free_y", ncol = 4) +
  scale_y_reordered() +
  ggtitle("Top Hallmark Pathways per Subcluster (GSEA)")

ggsave(
  paste0(res_dir, "barplot_top_10_hallmark_pathways_per_subcluster.pdf"),
  plot = last_plot(),
  width = 22, height = 16
)
