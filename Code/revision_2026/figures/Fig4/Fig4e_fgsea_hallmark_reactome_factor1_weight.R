#!/usr/bin/env Rscript
# Fig. 4e — gene-set enrichment (fgsea) of Hallmark and Reactome pathways on the 984 UK Biobank Olink proteins that are
# matched to the Factor 1 tubular weight, proteins ranked by that weight: the 10 terms with the smallest adjusted p in
# each direction (Hallmark and Reactome pooled), bars = normalised enrichment score (NES), dark = BH-adjusted p < 0.05.
# RE-DRAWING. The published panel was not made by code on this system: no script on disk writes it. This entry script
# re-draws it from Leon's fgsea table with the selection rule and style read off the published panel; terms, their
# order, the shading and the NES values were checked against the published page (see Verified).
# Panel        : Main Figures.pdf, Figure 4, panel e
# Reads        : $KB_ANALYSIS/factor1_reviewer_comment/outputs/fgsea_all_collections.csv, rows with ranking
#                  "Factor 1 weight, 984 matched proteins" and collection Hallmark or Reactome
#                  (analyses/07_factor1_ukb_protein_enrichment/04_composite_and_lowerleft.R, 8 Oct 2026, de_env:
#                  fgsea 1.32.2 with minSize 10, maxSize 500, eps 0; msigdbr 7.5.1 gene sets; seed 42; protein table from
#                  01_build_protein_table.R)
# Writes       : $KB_FIG_OUT/Fig4/Fig4e_fgsea_hallmark_reactome_factor1_weight.{pdf,png} and ..._terms.csv (the 20 rows
#                drawn: collection, pathway, label, NES, padj, pval, size)
# Plotting code: written for the release (re-drawing). Selection: per sign of NES, the 10 rows with the smallest padj,
#                Hallmark + Reactome pooled; bars ordered by NES. Labels: collection prefix removed, title case with
#                acronyms kept (ECM, UV, DN, IL6, JAK, STAT3, G2M) and "of" / "and" in lower case, "[collection]"
#                appended; one long label shortened exactly as printed ("... Remodeling and Clear."). Bar style as
#                mk() in analyses/07_factor1_ukb_protein_enrichment/05_figures_addendum.R (Okabe-Ito vermillion /
#                blue for NES > 0 / < 0, alpha 1 / 0.35 for padj < 0.05 / >= 0.05, theme_bw, legend at the bottom)
# Illustrator  : not applicable (no exported original). On the published page the bars render slightly different from
#                the two Okabe-Ito hues used here (dark bars about #D5581E / #3475BB instead of #D55E00 / #0072B2); label
#                font and panel proportions differ
# Environment  : conda env de_env (R 4.4.1, ggplot2 3.5.1, data.table 1.16.4)
# Verified     : 2026-11-12, RE-DRAWING (no original file exists, so no pixel comparison): the 20 labels (exact strings)
#                and their order equal the PDF text layer of Main Figures.pdf p. 4, the shading equals padj < 0.05 for all
#                20 bars, and the bar ends measured on a 600-dpi render give the table's NES within 0.004
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

RANKING <- "Factor 1 weight, 984 matched proteins"
fg <- fread(file.path(KB_ANALYSIS, "factor1_reviewer_comment", "outputs", "fgsea_all_collections.csv"))
w <- fg[ranking == RANKING & collection %in% c("Hallmark", "Reactome") & !is.na(NES)]
stopifnot(nrow(w) > 20)
pick <- function(d) {                                   # 10 smallest padj; refuse a tie across the cut
  d <- d[order(padj)]
  if (nrow(d) > 10 && d$padj[10] == d$padj[11]) stop("tie in padj at the 10th term")
  d[1:10]
}
sel <- rbind(pick(w[NES > 0]), pick(w[NES < 0]))

ACRONYMS <- c("ECM", "UV", "DN", "UP", "JAK", "TNFA", "NFKB", "KRAS", "TGF", "MYC", "UPR", "DNA", "RNA", "PI3K", "AKT",
              "MTOR", "MTORC1", "E2F", "WNT", "IFN", "HDL", "LDL", "TCA", "SLC", "ABC", "NMD", "GPCR", "MHC")
SMALL <- c("OF", "AND", "THE", "IN", "TO", "BY", "FOR", "ON", "AT", "WITH", "A", "AN", "OR", "VIA", "FROM")
SHORT <- c("Plasma Lipoprotein Assembly Remodeling and Clearance" = "Plasma Lipoprotein Assembly Remodeling and Clear.")
nice <- function(id) {
  x <- strsplit(sub("^(HALLMARK|REACTOME)_", "", id), "_")[[1]]
  x <- ifelse(x %in% ACRONYMS | grepl("[0-9]", x), x,
              ifelse(x %in% SMALL & seq_along(x) > 1, tolower(x), paste0(substr(x, 1, 1), tolower(substring(x, 2)))))
  x <- paste(x, collapse = " ")
  if (x %in% names(SHORT)) SHORT[[x]] else x
}
sel[, label := sprintf("%s [%s]", vapply(pathway, nice, ""), collection)]
sel <- sel[order(-NES)]
sel[, label := factor(label, levels = rev(label))]      # highest NES at the top
sel[, sig := padj < 0.05]
fwrite(sel[, .(collection, pathway, label, NES, padj, pval, size)],
       fig_out("Fig4", "Fig4e_fgsea_hallmark_reactome_factor1_weight_terms.csv"))

p <- ggplot(sel, aes(NES, label, fill = NES > 0, alpha = sig)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), breaks = c(FALSE, TRUE),
                     labels = c(`FALSE` = ">= 0.05", `TRUE` = "< 0.05"), name = "adjusted p-value") +
  labs(x = "Normalised enrichment score (fgsea)", y = NULL,
       title = "Hallmark & Reactome pathway, ranked by factor 1 tubular weight (984 matched proteins)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        plot.title.position = "plot", plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.y = element_text(size = 7.5, colour = "black"))
for (ext in c("pdf", "png"))
  ggsave(fig_out("Fig4", paste0("Fig4e_fgsea_hallmark_reactome_factor1_weight.", ext)), p, width = 6.6, height = 3.9,
         dpi = 300)
print(sel[, .(collection, label = as.character(label), NES = round(NES, 3), padj = signif(padj, 3), dark = sig)])
cat(sprintf("[Fig4e] %s\n", fig_out("Fig4", "Fig4e_fgsea_hallmark_reactome_factor1_weight.pdf")))
