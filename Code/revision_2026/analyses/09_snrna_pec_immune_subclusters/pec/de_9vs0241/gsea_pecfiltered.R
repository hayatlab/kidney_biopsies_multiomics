#!/usr/bin/env Rscript
# Re-run GSEA (Hallmark + Reactome) on the cl9-vs-{0,2,4,11} MAST DE AFTER removing canonical
# non-PEC lineage markers (the "_pecfiltered" table produced by filter_and_volcano_pecfiltered.py).
#
# Reuses the ranking + gene-set + GSEA logic of de_analysis_9vs0241.R VERBATIM (same seed 1234,
# minGSSize 10, maxGSSize 5000, BH, pvalueCutoff 0.05, msigdbr H + C2:CP:REACTOME). It reads the
# already-computed MAST table -- no MAST re-run. Outputs go to a NEW GSEA_pecfiltered/ dir so the
# original GSEA is preserved for side-by-side comparison.

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(clusterProfiler); library(msigdbr)
  library(enrichplot); library(tibble)
})

# ---- config (mirrors de_analysis_9vs0241.R) ----
DE_DIR   <- paste0(KB_ANALYSIS, "/pec_resubcluster_20260717/de_9vs0241/de_results/9_vs_0_2_4_11")
FILT_CSV <- file.path(DE_DIR, "DE", "markers_MAST_9_vs_0_2_4_11_pecfiltered.csv")
OUT_ROOT <- file.path(DE_DIR, "GSEA_pecfiltered")
contrast_title <- "leiden_0.8 cl9 vs cl0/2/4/11 (new PEC) - non-PEC markers removed"

seed_val <- 1234; show_categories <- 10
minGSSize <- 10; maxGSSize <- 5000; gsea_pcut <- 0.05
species <- "Homo sapiens"

safe_dir <- function(...) { p <- file.path(...); if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE); p }

initialize_gene_sets <- function(species = "Homo sapiens") {
  message(sprintf("Loading MSigDB Hallmark and Reactome for: %s", species))
  msig_h <- msigdbr(species = species, category = "H")
  hallmark_T2G <- msig_h %>% dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>% dplyr::distinct()
  hallmark_T2N <- msig_h %>% dplyr::select(gs_name, gs_description) %>%
    dplyr::rename(term = gs_name, name = gs_description) %>% dplyr::distinct()
  msig_r <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
  reactome_T2G <- msig_r %>% dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>% dplyr::distinct()
  reactome_T2N <- msig_r %>% dplyr::select(gs_name, gs_description) %>%
    dplyr::rename(term = gs_name, name = gs_description) %>% dplyr::distinct()
  list(hallmark = list(T2G = hallmark_T2G, T2N = hallmark_T2N),
       reactome = list(T2G = reactome_T2G, T2N = reactome_T2N))
}

shorten_pathway_names <- function(x, width = 35) {
  s <- x %>% str_replace("^REACTOME_", "") %>% str_replace("^HALLMARK_", "") %>% str_replace_all("_", " ")
  sapply(s, function(n) stringr::str_wrap(n, width = width))
}

rank_genes <- function(markers_df) {
  df <- markers_df %>% dplyr::select(gene, avg_log2FC) %>% dplyr::filter(!is.na(avg_log2FC))
  ranks <- df %>% dplyr::group_by(gene) %>% dplyr::summarise(stat = mean(avg_log2FC), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(stat))
  gl <- ranks$stat; names(gl) <- ranks$gene; gl
}

run_gsea_safe <- function(geneList, TERM2GENE, TERM2NAME, collection_name,
                          minGSSize = 10, maxGSSize = 5000, pvalueCutoff = 0.05, seed = 1234) {
  overlap <- intersect(names(geneList), TERM2GENE$gene)
  if (length(overlap) < minGSSize) { message(sprintf("Skip %s: overlap %d < %d", collection_name, length(overlap), minGSSize)); return(NULL) }
  gl <- sort(geneList[overlap], decreasing = TRUE)
  set.seed(seed)
  tryCatch(GSEA(geneList = gl, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff,
                pAdjustMethod = "BH", verbose = FALSE, seed = TRUE),
           error = function(e) { message(sprintf("GSEA failed %s: %s", collection_name, e$message)); NULL })
}

shorten_gsea_descriptions <- function(o) {
  if (is.null(o) || nrow(o@result) == 0) return(o)
  r <- o@result; r$Description <- shorten_pathway_names(r$ID); o@result <- r; o
}

dotplot_styled <- function(o, showCategory = 20, title = NULL, label_wrap_width = 40, base_size = 14) {
  if (is.null(o) || nrow(o@result) == 0) return(ggplot() + theme_void() + ggtitle(paste0(title, " (no enriched terms)")))
  o2 <- o
  o2@result$Description <- gsub("[\r\n]+", " ", o2@result$Description)  # flatten the width-35 pre-wrap
  # enrichplot::dotplot re-wraps y-labels at its own default (~30 chars) -> override the y scale
  # directly so long pathway names wrap at OUR width (fewer lines -> no vertical overlap).
  wrap_fun <- function(x) stringr::str_wrap(gsub("[\r\n]+", " ", x), width = label_wrap_width)
  dotplot(o2, showCategory = showCategory, title = title) +
    scale_y_discrete(labels = wrap_fun) + theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0, face = "bold"),
          axis.text.y = element_text(size = base_size * 0.9),
          axis.text.x = element_text(size = base_size * 0.9, angle = 0, vjust = 0.5),
          axis.title = element_text(size = base_size),
          legend.text = element_text(size = base_size * 0.9),
          legend.title = element_text(size = base_size))
}

save_gsea_outputs <- function(o, out_dir, dotplot_title, show_categories = 20,
                              fig_width = 10, fig_height = 8, base_size = 14, label_wrap_width = 40) {
  safe_dir(out_dir)
  readr::write_csv(if (!is.null(o) && nrow(o@result) > 0) o@result else tibble(), file.path(out_dir, "gsea_results.csv"))
  if (!is.null(o) && nrow(o@result) > 0) {
    p <- dotplot_styled(o, showCategory = show_categories, title = dotplot_title,
                        label_wrap_width = label_wrap_width, base_size = base_size)
    ggsave(file.path(out_dir, "dotplot.pdf"), p, width = fig_width, height = fig_height, device = cairo_pdf)
  } else message(sprintf("No enriched terms: %s", dotplot_title))
}

# ---- run ----
message("Reading filtered DE table: ", FILT_CSV)
markers <- readr::read_csv(FILT_CSV, show_col_types = FALSE)
message(sprintf("Filtered DE genes: %d", nrow(markers)))
gs <- initialize_gene_sets(species = species)
geneList <- rank_genes(markers)
message(sprintf("Ranked gene list: %d genes (range %.2f .. %.2f)", length(geneList), min(geneList), max(geneList)))

hallmark_dir <- safe_dir(OUT_ROOT, "Hallmark")
reactome_dir <- safe_dir(OUT_ROOT, "Reactome")

message("GSEA - Hallmark...")
gsea_h <- run_gsea_safe(geneList, gs$hallmark$T2G, gs$hallmark$T2N, "Hallmark",
                        minGSSize, maxGSSize, gsea_pcut, seed_val)
gsea_h <- shorten_gsea_descriptions(gsea_h)
save_gsea_outputs(gsea_h, hallmark_dir, paste0(contrast_title, " - Hallmark"), show_categories,
                  fig_width = 9, fig_height = 5.5, base_size = 13, label_wrap_width = 45)
message(sprintf("  Hallmark enriched terms: %d", if (is.null(gsea_h)) 0 else nrow(gsea_h@result)))

message("GSEA - Reactome...")
gsea_r <- run_gsea_safe(geneList, gs$reactome$T2G, gs$reactome$T2N, "Reactome",
                        minGSSize, maxGSSize, gsea_pcut, seed_val)
gsea_r <- shorten_gsea_descriptions(gsea_r)
# Reactome: top 10 terms, compact canvas so rows sit close together (little inter-row whitespace)
save_gsea_outputs(gsea_r, reactome_dir, paste0(contrast_title, " - Reactome"), show_categories,
                  fig_width = 10, fig_height = 6.0, base_size = 13, label_wrap_width = 50)
message(sprintf("  Reactome enriched terms: %d", if (is.null(gsea_r)) 0 else nrow(gsea_r@result)))

message("Done. Outputs -> ", OUT_ROOT)
