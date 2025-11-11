#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(zellkonverter)
  library(MAST)
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
file_path <- "/data/ep924610/kidney_biopsy/trajectory_pec/pec_annotated_for_DEGs.h5ad"
base_dir  <- "/data/ep924610/kidney_biopsy/trajectory_pec/de_results"

cluster_col <- "pec_subcluster_0.80"  
donor_var   <- "donor_id"              

# Volcano thresholds
p_cutoff  <- 0.05
fc_cutoff <- 0.25

# GSEA settings
seed_val        <- 1234
show_categories <- 20
minGSSize <- 10
maxGSSize <- 5000
gsea_pcut <- 0.05

# Species for msigdbr
species <- "Homo sapiens"

# ---------------------------
# Helpers
# ---------------------------
safe_dir <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  return(path)
}

choose_assay <- function(obj) {
  # Try to pick a reasonable assay to work with
  candidates <- c("RNA", "originalexp", "X", "counts")
  aa <- intersect(candidates, names(obj@assays))
  if (length(aa) > 0) return(aa[1])
  # Fallback to first assay if nothing matches
  return(names(obj@assays)[1])
}

initialize_gene_sets <- function(species = "Homo sapiens") {
  message(sprintf("Loading MSigDB Hallmark and Reactome gene sets for species: %s", species))

  # Hallmark
  msig_h <- msigdbr(species = species, category = "H")
  hallmark_T2G <- msig_h %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>%
    dplyr::distinct()
  hallmark_T2N <- msig_h %>%
    dplyr::select(gs_name, gs_description) %>%
    dplyr::rename(term = gs_name, name = gs_description) %>%
    dplyr::distinct()

  # Reactome
  msig_r <- msigdbr(species = species, category = "C2", subcategory = "CP:REACTOME")
  reactome_T2G <- msig_r %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename(term = gs_name, gene = gene_symbol) %>%
    dplyr::distinct()
  reactome_T2N <- msig_r %>%
    dplyr::select(gs_name, gs_description) %>%
    dplyr::rename(term = gs_name, name = gs_description) %>%
    dplyr::distinct()

  list(
    hallmark = list(T2G = hallmark_T2G, T2N = hallmark_T2N),
    reactome = list(T2G = reactome_T2G, T2N = reactome_T2N)
  )
}

shorten_pathway_names <- function(pathway_names, width = 35) {
  shortened <- pathway_names %>%
    str_replace("^REACTOME_", "") %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace_all("_", " ")
  sapply(shortened, function(name) stringr::str_wrap(name, width = width))
}

save_volcano <- function(markers_df, out_pdf, out_png, title,
                         p_cutoff = 0.05, fc_cutoff = 0.25) {
  if (nrow(markers_df) == 0) {
    message("No DE genes to plot.")
    return(invisible(NULL))
  }
  vp <- EnhancedVolcano(
    markers_df,
    lab = markers_df$gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    title = title,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    xlab = "log2 fold change",
    ylab = "-log10 adjusted p",
    drawConnectors = TRUE,
    legendLabels = c("NS","|log2FC|","p-adj","p-adj & |log2FC|")
  )
  ggsave(out_pdf, vp, width = 7, height = 6, device = cairo_pdf)
  ggsave(out_png, vp, width = 7, height = 6, dpi = 300)
}

rank_genes <- function(markers_df) {
  df <- markers_df %>% dplyr::select(gene, avg_log2FC) %>% dplyr::filter(!is.na(avg_log2FC))
  ranks <- df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(stat = mean(avg_log2FC), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(stat))
  geneList <- ranks$stat
  names(geneList) <- ranks$gene
  return(geneList)
}

run_gsea_safe <- function(geneList, TERM2GENE, TERM2NAME, collection_name,
                          minGSSize = 10, maxGSSize = 5000, pvalueCutoff = 0.05, seed = 1234) {
  overlap <- intersect(names(geneList), TERM2GENE$gene)
  if (length(overlap) < minGSSize) {
    message(sprintf("Skipping GSEA for %s: insufficient overlap (%d < %d).",
                    collection_name, length(overlap), minGSSize))
    return(NULL)
  }
  geneList2 <- geneList[overlap]
  geneList2 <- sort(geneList2, decreasing = TRUE)
  set.seed(seed)
  res <- tryCatch({
    GSEA(
      geneList      = geneList2,
      TERM2GENE     = TERM2GENE,
      TERM2NAME     = TERM2NAME,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      pvalueCutoff  = pvalueCutoff,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = TRUE
    )
  }, error = function(e) {
    message(sprintf("GSEA failed for %s: %s", collection_name, e$message))
    return(NULL)
  })
  res
}

shorten_gsea_descriptions <- function(gsea_obj) {
  if (is.null(gsea_obj)) return(gsea_obj)
  if (nrow(gsea_obj@result) == 0) return(gsea_obj)
  res <- gsea_obj@result
  res$Description <- shorten_pathway_names(res$ID)
  gsea_obj@result <- res
  gsea_obj
}

dotplot_styled <- function(gsea_obj, showCategory = 20, title = NULL, label_wrap_width = 40, base_size = 14) {
  if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0(title, " (no enriched terms)")))
  }
  gsea_obj2 <- gsea_obj
  gsea_obj2@result$Description <- stringr::str_wrap(gsea_obj2@result$Description, width = label_wrap_width)
  dotplot(gsea_obj2, showCategory = showCategory, title = title) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.y = element_text(size = base_size * 0.9),
      axis.text.x = element_text(size = base_size * 0.9, angle = 0, vjust = 0.5),
      axis.title  = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.9),
      legend.title= element_text(size = base_size)
    )
}

save_gsea_outputs <- function(gsea_obj, out_dir, dotplot_title, show_categories = 20,
                              fig_width = 10, fig_height = 8, base_size = 14, label_wrap_width = 40) {
  safe_dir(out_dir)
  res_csv <- file.path(out_dir, "gsea_results.csv")
  plot_pdf <- file.path(out_dir, "dotplot.pdf")
  res_tbl <- if (!is.null(gsea_obj) && nrow(gsea_obj@result) > 0) gsea_obj@result else tibble()
  readr::write_csv(res_tbl, res_csv)
  if (!is.null(gsea_obj) && nrow(gsea_obj@result) > 0) {
    p <- dotplot_styled(
      gsea_obj,
      showCategory = show_categories,
      title = dotplot_title,
      label_wrap_width = label_wrap_width,
      base_size = base_size
    )
    ggsave(plot_pdf, p, width = fig_width, height = fig_height, device = cairo_pdf)
  } else {
    message(sprintf("No enriched terms for: %s", dotplot_title))
  }
}

# ---------------------------
# Load data and preprocess
# ---------------------------
message("Reading AnnData and converting to Seurat...")
sce <- readH5AD(file_path)
seu <- as.Seurat(sce, counts = "X", data = NULL)

assay_to_use <- choose_assay(seu)
DefaultAssay(seu) <- assay_to_use
message(sprintf("Using assay: %s", assay_to_use))

message("Normalizing data (LogNormalize, scale.factor = 1e4)...")
seu <- NormalizeData(seu, assay = assay_to_use,
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)

# ---------------------------
# Basic checks
# ---------------------------
meta_cols <- colnames(seu@meta.data)
if (!(cluster_col %in% meta_cols)) {
  stop(sprintf("Column '%s' not found in metadata. Available columns: %s",
               cluster_col, paste(meta_cols, collapse = ", ")))
}
if (!(donor_var %in% meta_cols)) {
  stop(sprintf("Column '%s' not found in metadata. Available columns: %s",
               donor_var, paste(meta_cols, collapse = ", ")))
}

unique_levels <- sort(unique(seu@meta.data[[cluster_col]]))
message(sprintf("Found %d levels in %s: %s",
                length(unique_levels), cluster_col, paste(unique_levels, collapse = ", ")))

required_lvls <- c("7_aPEC", "2_aPECs")
missing_lvls <- setdiff(required_lvls, unique_levels)
if (length(missing_lvls) > 0) {
  stop(sprintf("Required levels missing in %s: %s", cluster_col, paste(missing_lvls, collapse = ", ")))
}

# ---------------------------
# Prepare outputs and gene sets
# ---------------------------
de_root   <- safe_dir(base_dir)
gsea_root <- de_root  # keep under the same base dir, organized per contrast

gs <- initialize_gene_sets(species = species)

# ---------------------------
# Contrast runner
# ---------------------------
run_contrast <- function(seu, targets, contrast_code, contrast_title) {
  message(sprintf("\n--- Running contrast: %s ---", contrast_title))

  # Prepare directories
  out_root      <- safe_dir(de_root, contrast_code)
  out_de_dir    <- safe_dir(out_root, "DE")
  out_gsea_root <- safe_dir(out_root, "GSEA")
  hallmark_dir  <- safe_dir(out_gsea_root, "Hallmark")
  reactome_dir  <- safe_dir(out_gsea_root, "Reactome")

  # Build grouping
  grp_name <- paste0("grp_", contrast_code)
  seu[[grp_name]] <- ifelse(seu@meta.data[[cluster_col]] %in% targets, "Target", "Other")
  Idents(seu) <- factor(seu[[grp_name]][,1])
  Idents(seu) <- droplevels(Idents(seu))

  present_lvls <- levels(Idents(seu))
  if (!all(c("Target", "Other") %in% present_lvls)) {
    stop(sprintf("Grouping failed for %s. Present levels: %s", contrast_title, paste(present_lvls, collapse=", ")))
  }

  # DE with MAST, donor_id as latent variable
  message("Running MAST DE (with donor_id as latent var)...")
  set.seed(seed_val)
  markers <- FindMarkers(
    object = seu,
    ident.1 = "Target",
    ident.2 = "Other",
    test.use = "MAST",
    min.pct = 0,
    logfc.threshold = 0,
    latent.vars = donor_var,
    verbose = FALSE
  ) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))

  # Save DE table (CSV)
  markers_csv <- file.path(out_de_dir, paste0("markers_MAST_", contrast_code, ".csv"))
  readr::write_csv(markers, markers_csv)
  message(sprintf("DE table saved: %s", markers_csv))

  # Volcano plots
  volcano_pdf <- file.path(out_de_dir, paste0("volcano_", contrast_code, ".pdf"))
  volcano_png <- file.path(out_de_dir, paste0("volcano_", contrast_code, ".png"))
  volcano_title <- paste0(contrast_title, " (MAST)")
  save_volcano(markers, volcano_pdf, volcano_png, volcano_title,
               p_cutoff = p_cutoff, fc_cutoff = fc_cutoff)
  message(sprintf("Volcano saved: %s ; %s", volcano_pdf, volcano_png))

  # GSEA (Hallmark and Reactome) using preranked list
  geneList <- rank_genes(markers)

  # Hallmark
  message("Running GSEA - Hallmark...")
  gsea_h <- run_gsea_safe(
    geneList = geneList,
    TERM2GENE = gs$hallmark$T2G,
    TERM2NAME = gs$hallmark$T2N,
    collection_name = "Hallmark",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = gsea_pcut,
    seed = seed_val
  )
  gsea_h <- shorten_gsea_descriptions(gsea_h)
  save_gsea_outputs(gsea_h, hallmark_dir,
                    dotplot_title = paste0(contrast_title, " - Hallmark"),
                    show_categories = show_categories)

  # Reactome
  message("Running GSEA - Reactome...")
  gsea_r <- run_gsea_safe(
    geneList = geneList,
    TERM2GENE = gs$reactome$T2G,
    TERM2NAME = gs$reactome$T2N,
    collection_name = "Reactome",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = gsea_pcut,
    seed = seed_val
  )
  gsea_r <- shorten_gsea_descriptions(gsea_r)
  save_gsea_outputs(gsea_r, reactome_dir,
                    dotplot_title = paste0(contrast_title, " - Reactome"),
                    show_categories = show_categories)

  invisible(list(markers = markers, gsea_h = gsea_h, gsea_r = gsea_r))
}

# ---------------------------
# Run requested contrasts
# ---------------------------
contrast_list <- list(
  list(
    code = "7_aPEC_vs_all",
    title = "7_aPEC vs all other (pec_subcluster_0.80)",
    targets = c("7_aPEC")
  ),
  list(
    code = "7_aPEC_2_aPECs_vs_all",
    title = "7_aPEC + 2_aPECs vs all other (pec_subcluster_0.80)",
    targets = c("7_aPEC", "2_aPECs")
  )
)

results <- list()
for (cc in contrast_list) {
  results[[cc$code]] <- run_contrast(
    seu = seu,
    targets = cc$targets,
    contrast_code = cc$code,
    contrast_title = cc$title
  )
}

# ---------------------------
# Save session info
# ---------------------------
sink(file.path(base_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("All analyses completed. Outputs saved under: ")
message(sprintf("  %s", base_dir))