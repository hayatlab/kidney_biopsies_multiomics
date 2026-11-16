#!/usr/bin/env Rscript

# Volcano + GSEA for NICHE 20 (crescent PEC vs other PEC), the manual_glom
# NicheCompass res-1.0 crescent niche (241 cells, 79% PEC+Imm).
#
# This is the niche-20 analogue of niche9_de_pec_subset.R. Unlike the niche-9
# script it does NOT re-run MAST: the niche-20 DE table is already the canonical
# one produced by de_crescent_pec_mast.R (crescent-PEC grp=="crescent" vs
# other-PEC) and is used by every other niche-20 downstream (dotplots, pseudotime
# validation, comparison). We build the volcano + GSEA directly on that table so
# the niche-20 figures are perfectly consistent with those.
#
#   input  : celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv
#            (columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, gene; 468 genes)
#   outputs: niche20_pec/DE/niche20_vs_others_pecOnly/volcano_niche20_vs_others_pecOnly.pdf
#            niche20_pec/GSEA/niche20_vs_others_pecOnly/{Hallmark,Reactome,GO_BP}/{gsea_results.tsv,dotplot.pdf}
#
# Volcano + GSEA settings are copied verbatim from niche9_de_pec_subset.R
# (EnhancedVolcano avg_log2FC vs p_val_adj, pCut 0.05 FCcut 0.25; clusterProfiler
# GSEA preranked by avg_log2FC over Hallmark / Reactome (C2:CP:REACTOME) / GO_BP,
# seed 1234, minGSSize 10, maxGSSize 5000, BH, pvalueCutoff 0.05).
#
# NOTE (same caveat as niche 9): the DE universe is the ~468-gene Xenium panel,
# a very small background for GSEA, so Reactome / GO_BP frequently return no terms
# surviving FDR (see gsea_diagnostic.R). Hallmark is the informative collection here.

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(tibble)
})

# ---------------------------
# Config
# ---------------------------
de_csv   <- paste0(KB_ANALYSIS, "/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv")
base_dir <- paste0(KB_ANALYSIS, "/niche20_pec")
nm       <- "niche20_vs_others_pecOnly"

# Volcano thresholds
p_cutoff  <- 0.05
fc_cutoff <- 0.25

# GSEA settings
seed_val        <- 1234
show_categories <- 20
minGSSize <- 10
maxGSSize <- 5000
gsea_pcut <- 0.05

go_ontology <- "BP"
species <- "Homo sapiens"

# ---------------------------
# Helpers  (verbatim from niche9_de_pec_subset.R)
# ---------------------------
safe_dir <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  return(path)
}

initialize_gene_sets <- function(species = "Homo sapiens", go_ontology = "BP") {
  message(sprintf("Loading MSigDB gene sets for species: %s", species))
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

  msig_go <- msigdbr(species = species, category = "C5")
  patt <- switch(go_ontology,
                 "BP" = "(^GOBP_|^GO_BIOLOGICAL_PROCESS_)",
                 "CC" = "(^GOCC_|^GO_CELLULAR_COMPONENT_)",
                 "MF" = "(^GOMF_|^GO_MOLECULAR_FUNCTION_)",
                 "(^GOBP_|^GO_BIOLOGICAL_PROCESS_)")
  msig_go_sel <- msig_go %>% dplyr::filter(stringr::str_detect(gs_name, patt))
  go_T2G <- msig_go_sel %>% dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>% dplyr::distinct()
  go_T2N <- msig_go_sel %>% dplyr::select(gs_name, gs_description) %>%
    dplyr::rename(term = gs_name, name = gs_description) %>% dplyr::distinct()

  list(hallmark = list(T2G = hallmark_T2G, T2N = hallmark_T2N),
       reactome = list(T2G = reactome_T2G, T2N = reactome_T2N),
       go = list(T2G = go_T2G, T2N = go_T2N))
}

shorten_pathway_names <- function(pathway_names, width = 35) {
  shortened <- pathway_names %>%
    str_replace("^REACTOME_", "") %>% str_replace("^HALLMARK_", "") %>%
    str_replace("^GO_BIOLOGICAL_PROCESS_", "GO:BP ") %>%
    str_replace("^GO_CELLULAR_COMPONENT_", "GO:CC ") %>%
    str_replace("^GO_MOLECULAR_FUNCTION_", "GO:MF ") %>%
    str_replace("^GOBP_", "GO:BP ") %>% str_replace("^GOCC_", "GO:CC ") %>%
    str_replace("^GOMF_", "GO:MF ")
  shortened <- str_replace_all(shortened, "_", " ")
  sapply(shortened, function(name) stringr::str_wrap(name, width = width))
}

save_volcano <- function(markers_df, out_pdf, title, p_cutoff = 0.05, fc_cutoff = 0.25) {
  if (nrow(markers_df) == 0) return(invisible(NULL))
  pdf(out_pdf, width = 7, height = 6)
  print(EnhancedVolcano(markers_df, lab = markers_df$gene,
                        x = "avg_log2FC", y = "p_val_adj", title = title,
                        pCutoff = p_cutoff, FCcutoff = fc_cutoff,
                        xlab = "log2 fold change", ylab = "-log10 adjusted p",
                        drawConnectors = TRUE,
                        legendLabels = c("NS","|log2FC|","p-adj","p-adj & |log2FC|")))
  dev.off()
}

rank_genes <- function(markers_df, exclude_genes = NULL) {
  df <- markers_df %>% dplyr::select(gene, avg_log2FC) %>% dplyr::filter(!is.na(avg_log2FC))
  if (!is.null(exclude_genes) && length(exclude_genes) > 0) {
    df <- df %>% dplyr::filter(!tolower(gene) %in% tolower(exclude_genes))
  }
  ranks <- df %>% dplyr::group_by(gene) %>%
    dplyr::summarise(stat = mean(avg_log2FC), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(stat))
  geneList <- ranks$stat; names(geneList) <- ranks$gene; geneList
}

run_gsea_safe <- function(geneList, TERM2GENE, TERM2NAME, collection_name,
                          minGSSize = 10, maxGSSize = 5000, pvalueCutoff = 0.05, seed = 1234) {
  overlap <- intersect(names(geneList), TERM2GENE$gene)
  if (length(overlap) < minGSSize) {
    message(sprintf("Skipping GSEA for %s: insufficient overlap (%d < %d).",
                    collection_name, length(overlap), minGSSize))
    return(NULL)
  }
  geneList2 <- sort(geneList[overlap], decreasing = TRUE)
  set.seed(seed)
  tryCatch(GSEA(geneList = geneList2, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME,
                minGSSize = minGSSize, maxGSSize = maxGSSize,
                pvalueCutoff = pvalueCutoff, pAdjustMethod = "BH",
                verbose = FALSE, seed = TRUE),
           error = function(e) { message(sprintf("GSEA failed for %s: %s", collection_name, e$message)); NULL })
}

shorten_gsea_descriptions <- function(gsea_obj) {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) return(gsea_obj)
  res <- gsea_obj@result; res$Description <- shorten_pathway_names(res$ID)
  gsea_obj@result <- res; gsea_obj
}

dotplot_styled <- function(gsea_obj, showCategory = 20, title = NULL, label_wrap_width = 40, base_size = 14) {
  p <- dotplot(gsea_obj, showCategory = showCategory, title = title) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0, face = "bold"),
          axis.text.y = element_text(size = base_size * 0.9),
          axis.text.x = element_text(size = base_size * 0.9, angle = 0, vjust = 0.5),
          axis.title = element_text(size = base_size),
          legend.text = element_text(size = base_size * 0.9),
          legend.title = element_text(size = base_size))
  if (!is.null(gsea_obj) && nrow(gsea_obj@result) > 0 && "Description" %in% colnames(gsea_obj@result)) {
    gsea_obj2 <- gsea_obj
    gsea_obj2@result$Description <- stringr::str_wrap(gsea_obj2@result$Description, width = label_wrap_width)
    p <- dotplot(gsea_obj2, showCategory = showCategory, title = title) +
      theme_bw(base_size = base_size) +
      theme(plot.title = element_text(hjust = 0, face = "bold"),
            axis.text.y = element_text(size = base_size * 0.9),
            axis.text.x = element_text(size = base_size * 0.9, angle = 0, vjust = 0.5),
            axis.title = element_text(size = base_size),
            legend.text = element_text(size = base_size * 0.9),
            legend.title = element_text(size = base_size))
  }
  p
}

save_gsea_outputs <- function(gsea_obj, out_dir, dotplot_title, show_categories = 20,
                              fig_width = 10, fig_height = 8, base_size = 14, label_wrap_width = 40) {
  safe_dir(out_dir)
  res_tsv <- file.path(out_dir, "gsea_results.tsv")
  plot_pdf <- file.path(out_dir, "dotplot.pdf")
  res_tbl <- if (!is.null(gsea_obj) && nrow(gsea_obj@result) > 0) gsea_obj@result else tibble()
  readr::write_tsv(res_tbl, res_tsv)
  if (!is.null(gsea_obj) && nrow(gsea_obj@result) > 0) {
    pdf(plot_pdf, width = fig_width, height = fig_height)
    print(dotplot_styled(gsea_obj, showCategory = show_categories,
                         title = dotplot_title, label_wrap_width = label_wrap_width,
                         base_size = base_size))
    dev.off()
  } else {
    message(sprintf("No enriched terms for: %s", dotplot_title))
  }
}

# ---------------------------
# Load the niche-20 DE table (no MAST re-run)
# ---------------------------
message(sprintf("Reading niche-20 DE table: %s", de_csv))
markers <- readr::read_csv(de_csv, show_col_types = FALSE)
# read_csv keeps "pct.1"/"pct.2" as-is; ensure the columns EnhancedVolcano/GSEA need exist
stopifnot(all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(markers)))
markers <- markers %>% dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
message(sprintf("  %d genes; UP(lfc>%.2f & padj<%.2f)=%d, DOWN=%d",
                nrow(markers), fc_cutoff, p_cutoff,
                sum(markers$avg_log2FC > fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE),
                sum(markers$avg_log2FC < -fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE)))

de_root      <- safe_dir(base_dir, "DE", nm)
gsea_root    <- safe_dir(base_dir, "GSEA", nm)
hallmark_dir <- safe_dir(gsea_root, "Hallmark")
reactome_dir <- safe_dir(gsea_root, "Reactome")
go_dir       <- safe_dir(gsea_root, paste0("GO_", go_ontology))

# copy the DE table into the niche20_pec tree for self-containment
markers_out <- file.path(de_root, sprintf("markers_MAST_%s.csv", nm))
readr::write_csv(markers, markers_out)
message(sprintf("Wrote DE table copy: %s", markers_out))

# ---------------------------
# Volcano
# ---------------------------
volcano_pdf   <- file.path(de_root, sprintf("volcano_%s.pdf", nm))
volcano_title <- "niche 20 (crescent PEC) vs other PEC"
save_volcano(markers, volcano_pdf, volcano_title, p_cutoff = p_cutoff, fc_cutoff = fc_cutoff)
message(sprintf("Wrote volcano: %s", volcano_pdf))

# ---------------------------
# GSEA (Hallmark / Reactome / GO_BP), preranked by avg_log2FC
# ---------------------------
gs <- initialize_gene_sets(species = species, go_ontology = go_ontology)
geneList <- rank_genes(markers)

gsea_h <- shorten_gsea_descriptions(run_gsea_safe(geneList, gs$hallmark$T2G, gs$hallmark$T2N,
                                                  "Hallmark", minGSSize, maxGSSize, gsea_pcut, seed_val))
save_gsea_outputs(gsea_h, hallmark_dir, sprintf("%s - Hallmark", nm), show_categories)

gsea_r <- shorten_gsea_descriptions(run_gsea_safe(geneList, gs$reactome$T2G, gs$reactome$T2N,
                                                  "Reactome", minGSSize, maxGSSize, gsea_pcut, seed_val))
save_gsea_outputs(gsea_r, reactome_dir, sprintf("%s - Reactome", nm), show_categories)

gsea_go <- shorten_gsea_descriptions(run_gsea_safe(geneList, gs$go$T2G, gs$go$T2N,
                                                   paste0("GO_", go_ontology), minGSSize, maxGSSize, gsea_pcut, seed_val))
save_gsea_outputs(gsea_go, go_dir, sprintf("%s - GO %s", nm, go_ontology), show_categories)

message("Done. Outputs under:")
message(sprintf("  DE:   %s", de_root))
message(sprintf("  GSEA: %s", gsea_root))
