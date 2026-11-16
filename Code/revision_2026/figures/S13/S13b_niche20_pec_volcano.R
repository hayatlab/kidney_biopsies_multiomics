#!/usr/bin/env Rscript
# Supp. Fig. 13b — Volcano plot of the MAST DE between Xenium PECs in glomerular niche 20 and PECs in all other glomerular niches.
# Panel        : All Supp Figures.pdf, Supplementary Figure 13, panel b
# Reads        : $KB_ANALYSIS/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv
#                (analyses/08_glomerular_niche20_crescent/de_crescent_pec_mast.R on the output of 00_build_pec_manualglom10.py)
# Writes       : $KB_FIG_OUT/S13/S13b_niche20_pec_volcano.pdf
# Plotting code: extracted verbatim from analyses/08_glomerular_niche20_crescent/niche20_volcano_gsea.R (save_volcano and
#                the table preparation before it; the script cannot be sourced without re-running its GSEA). Same
#                libraries loaded in the same order.
# Illustrator  : title replaced ("Xenium PEC glom niche 20 vs PEC all other glom niches"); "EnhancedVolcano" subtitle
#                removed; label positions adjusted
# Environment  : conda env de_env
# Verified     : 2026-11-12, DIFFERENT (label placement only) vs $KB_ANALYSIS/niche20_pec/DE/niche20_vs_others_pecOnly/volcano_niche20_vs_others_pecOnly.pdf
#                EnhancedVolcano draws the gene labels with ggrepel without a seed (the original run set none before the
#                volcano), so label and connector positions change from run to run; points, axes, thresholds, legend
#                and caption are the same (see build/verification/niches.csv).

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

de_csv    <- paste0(KB_ANALYSIS, "/celltype_concordance/outputs/crescent_check/de/markers_manualglom10.csv")
p_cutoff  <- 0.05
fc_cutoff <- 0.25

# ---- verbatim from niche20_volcano_gsea.R ----
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

markers <- readr::read_csv(de_csv, show_col_types = FALSE)
stopifnot(all(c("gene", "avg_log2FC", "p_val_adj") %in% colnames(markers)))
markers <- markers %>% dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))
message(sprintf("  %d genes; UP(lfc>%.2f & padj<%.2f)=%d, DOWN=%d",
                nrow(markers), fc_cutoff, p_cutoff,
                sum(markers$avg_log2FC > fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE),
                sum(markers$avg_log2FC < -fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE)))

volcano_pdf   <- fig_out("S13", "S13b_niche20_pec_volcano.pdf")
volcano_title <- "niche 20 (crescent PEC) vs other PEC"
save_volcano(markers, volcano_pdf, volcano_title, p_cutoff = p_cutoff, fc_cutoff = fc_cutoff)
message(sprintf("Wrote volcano: %s", volcano_pdf))
