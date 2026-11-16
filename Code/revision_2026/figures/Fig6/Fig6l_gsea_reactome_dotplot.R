#!/usr/bin/env Rscript
# Fig. 6l — GSEA (Reactome) of the snRNA-seq comparison aPEC (PEC subcluster 9) vs PEC subclusters 0, 2, 4 and 11
# (MAST log2FC ranking after removal of canonical non-PEC lineage markers): top 10 enriched terms.
# Panel        : Main Figures.pdf, Figure 6, panel l
# Reads        : $KB_ANALYSIS/pec_resubcluster_20260717/de_9vs0241/de_results/9_vs_0_2_4_11/GSEA_pecfiltered/Reactome/gsea_results.csv
#                  (analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/gsea_pecfiltered.R, wrapper
#                   run_gsea_pecfiltered.sh, job 328394 of 2026-07-30: clusterProfiler::GSEA (fgsea multilevel, seed
#                   1234) on markers_MAST_9_vs_0_2_4_11_pecfiltered.csv from filter_and_volcano_pecfiltered.py)
# Writes       : $KB_FIG_OUT/Fig6/Fig6l_gsea_reactome_dotplot.pdf
# Plotting code: analyses/09_snrna_pec_immune_subclusters/pec/de_9vs0241/gsea_pecfiltered.R::dotplot_styled and its
#                config values (contrast_title, show_categories, OUT_ROOT) — evaluated from the shipped script by name,
#                without running the GSEA; the saved result table is wrapped in a DOSE gseaResult object (dotplot()
#                uses only its result table) and plotted with the arguments of the original Reactome call. The GSEA
#                is NOT re-run (permutation based).
# Illustrator  : title replaced ("Reactome / aPEC (subcluster 9) vs PEC subclusters 0, 2, 4 and 11 (snRNA-seq)")
# Environment  : conda env de_env
# Verified     : 2026-11-12, PIXEL-IDENTICAL (text layer equal) vs $KB_ANALYSIS/pec_resubcluster_20260717/de_9vs0241/
#                de_results/9_vs_0_2_4_11/GSEA_pecfiltered/Reactome/dotplot.pdf
# Writes in this script: only ggsave(fig_out(...)) below.

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))
suppressPackageStartupMessages({                  # same packages, same order as gsea_pecfiltered.R
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(clusterProfiler); library(msigdbr)
  library(enrichplot); library(tibble)
})

# ---- evaluate only these top-level assignments of the shipped analysis script (no GSEA, no writes) ----
src <- file.path(KB_REPO, "analyses", "09_snrna_pec_immune_subclusters", "pec", "de_9vs0241", "gsea_pecfiltered.R")
wanted <- c("DE_DIR", "OUT_ROOT", "contrast_title", "show_categories", "dotplot_styled")
for (e in parse(src, keep.source = FALSE)) {
  if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% wanted) eval(e, envir = globalenv())
}
stopifnot(all(vapply(wanted, exists, logical(1))))

res <- readr::read_csv(file.path(OUT_ROOT, "Reactome", "gsea_results.csv"), show_col_types = FALSE)
res <- as.data.frame(res); rownames(res) <- res$ID
gsea_r <- new("gseaResult", result = res, organism = "UNKNOWN", setType = "UNKNOWN", geneSets = list(),
              geneList = numeric(), keytype = "UNKNOWN", permScores = matrix(), params = list(),
              gene2Symbol = character(), readable = FALSE, termsim = matrix(), method = character(), dr = list())

# arguments of the original call: save_gsea_outputs(gsea_r, reactome_dir, paste0(contrast_title, " - Reactome"),
#   show_categories, fig_width = 10, fig_height = 6.0, base_size = 13, label_wrap_width = 50)
p <- dotplot_styled(gsea_r, showCategory = show_categories, title = paste0(contrast_title, " - Reactome"),
                    label_wrap_width = 50, base_size = 13)
out <- fig_out("Fig6", "Fig6l_gsea_reactome_dotplot.pdf")
ggsave(out, p, width = 10, height = 6.0, device = cairo_pdf)
message("[Fig6l] wrote ", out)
