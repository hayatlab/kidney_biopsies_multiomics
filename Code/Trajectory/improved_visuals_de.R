# Dependencies
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(EnhancedVolcano)
  library(scales)
})

# Utility: Safe directory creation
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# ---------------------------
# 1) Volcano plot with top 10 up/down labels
# ---------------------------

save_volcano_top10 <- function(
  de_csv,
  out_pdf,
  out_png,
  title,
  p_cutoff = 0.05,
  fc_cutoff = 0.25,
  top_n = 10,
  base_size = 14
) {
  markers_df <- readr::read_csv(de_csv, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(gene), !is.na(avg_log2FC), !is.na(p_val_adj))

  if (nrow(markers_df) == 0) {
    message("No DE genes to plot.")
    return(invisible(NULL))
  }

  # Determine top 10 up/down among significant genes
  sig_df <- markers_df %>% dplyr::filter(p_val_adj < p_cutoff)
  top_up <- sig_df %>% dplyr::arrange(dplyr::desc(avg_log2FC)) %>% dplyr::slice_head(n = top_n)
  top_down <- sig_df %>% dplyr::arrange(avg_log2FC) %>% dplyr::slice_head(n = top_n)
  select_genes <- unique(c(top_up$gene, top_down$gene))

  # Build volcano using EnhancedVolcano
  vp <- EnhancedVolcano(
    markers_df,
    lab = markers_df$gene,
    selectLab = select_genes,
    x = "avg_log2FC",
    y = "p_val_adj",
    title = title,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    xlab = "log2 fold change",
    ylab = "-log10 adjusted p",
    drawConnectors = TRUE,
    legendLabels = c("NS","|log2FC|","p-adj","p-adj & |log2FC|"),
    labSize = 4.0,
    xlim = c(-6, 6),
    pointSize = 1.2,
    colAlpha = 0.8,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'bottom',
  ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text  = element_text(size = base_size * 0.9),
      axis.title = element_text(size = base_size)
    )

  # Save outputs
  ggsave(out_pdf, vp, width = 10, height = 6, device = cairo_pdf)
  ggsave(out_png, vp, width = 10, height = 6, dpi = 300)
  invisible(list(plot = vp, selected = select_genes))
}

# ---------------------------
# 2) GSEA dotplot from CSV (no line breaks in pathway names)
#    Styled similarly to a clusterProfiler dotplot
# ---------------------------

# Creates a dotplot-like bubble chart with:
#   x = NES
#   y = Description (no line breaks)
#   size = setSize
#   color = p.adjust (smaller = red)
save_gsea_dotplot_from_csv <- function(
  gsea_csv,
  out_pdf,
  out_png,
  title = "Reactome GSEA",
  showCategory = 20,
  order_by = c("p.adjust", "NES", "absNES"),
  base_size = 14
) {
  order_by <- match.arg(order_by)

  gsea_df <- readr::read_csv(gsea_csv, show_col_types = FALSE)

  # Clean description: remove line breaks so it's one line per pathway
  gsea_df <- gsea_df %>%
    dplyr::mutate(
      Description = stringr::str_replace_all(Description, "[\\r\\n]+", " "),
      absNES = abs(NES)
    ) %>%
    dplyr::filter(!is.na(Description), !is.na(NES), !is.na(setSize), !is.na(p.adjust))

  if (nrow(gsea_df) == 0) {
    message("No GSEA rows found.")
    p <- ggplot() + theme_void() + ggtitle(paste0(title, " (no enriched terms)"))
    ggsave(out_pdf, p, width = 15, height = 8, device = cairo_pdf)
    ggsave(out_png, p, width = 15, height = 8, dpi = 300)
    return(invisible(NULL))
  }

  # Choose top categories
  df_top <- switch(
    order_by,
    "p.adjust" = gsea_df %>% arrange(p.adjust, desc(absNES)) %>% slice_head(n = showCategory),
    "NES"      = gsea_df %>% arrange(desc(NES), p.adjust) %>% slice_head(n = showCategory),
    "absNES"   = gsea_df %>% arrange(desc(absNES), p.adjust) %>% slice_head(n = showCategory)
  )

  # Ensure the y-axis ordering is stable (show strongest at top)
  df_top <- df_top %>%
    dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

  p <- ggplot(df_top, aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    scale_size_area(max_size = 10, name = "Set size") +
    scale_color_gradient(
      low = "#D73027", high = "#4575B4",
      name = "p.adjust", trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::label_number(accuracy = 0.0001, scale_cut = cut_short_scale())
    ) +
    labs(title = title, x = "NES", y = NULL) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.y = element_text(size = base_size * 0.9),
      axis.text.x = element_text(size = base_size * 0.9),
      axis.title  = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.9),
      legend.title= element_text(size = base_size)
    )

  ggsave(out_pdf, p, width = 15, height = 8, device = cairo_pdf)
  ggsave(out_png, p, width = 15, height = 8, dpi = 300)
  invisible(p)
}

# ---------------------------
# Example usage with your files
# ---------------------------

# Paths you provided
de_csv_path <- "/data/ep924610/kidney_biopsy/trajectory_pec/de_results/7_aPEC_2_aPECs_vs_all/DE/markers_MAST_7_aPEC_2_aPECs_vs_all.csv"
gsea_csv_path <- "/data/ep924610/kidney_biopsy/trajectory_pec/de_results/7_aPEC_2_aPECs_vs_all/GSEA/Reactome/gsea_results.csv"

# Output directories
out_dir_de <- "/data/ep924610/kidney_biopsy/trajectory_pec/de_results/7_aPEC_2_aPECs_vs_all/plots"
out_dir_gsea <- "/data/ep924610/kidney_biopsy/trajectory_pec/de_results/7_aPEC_2_aPECs_vs_all/GSEA/Reactome/plots"
safe_dir(out_dir_de)
safe_dir(out_dir_gsea)

# Volcano with top-10 up/down labels
save_volcano_top10(
  de_csv   = de_csv_path,
  out_pdf  = file.path(out_dir_de, "volcano_no_genes.pdf"),
  out_png  = file.path(out_dir_de, "volcano_no_genes.png"),
  title    = "7_aPEC_2_aPECs vs all (MAST)",
  p_cutoff = 0.05,
  fc_cutoff= 0.25,
  top_n    = 0
)

# GSEA dotplot from CSV (no line breaks in labels)
save_gsea_dotplot_from_csv(
  gsea_csv = gsea_csv_path,
  out_pdf  = file.path(out_dir_gsea, "reactome_gsea_dotplot_top20.pdf"),
  out_png  = file.path(out_dir_gsea, "reactome_gsea_dotplot_top20.png"),
  title    = "Reactome GSEA (one line per pathway)",
  showCategory = 20,
  order_by = "NES"  # or "NES" / "absNES"
)