#!/usr/bin/env Rscript
# Addendum figures: (5) top 40 by composite rank; (6) lower-left quadrant proteins + ORA; (7) fgsea tails on Cox z and on weight.
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressMessages({library(data.table); library(ggplot2); library(ggrepel); library(patchwork)})
root <- paste0(KB_ANALYSIS, "/factor1_reviewer_comment")
tab <- fread(file.path(root, "outputs/protein_table.csv")); ora <- fread(file.path(root, "outputs/ora_lowerleft.csv"))
fg  <- fread(file.path(root, "outputs/fgsea_all_collections.csv"))
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.key.size = unit(0.4, "cm"),
                                      plot.title = element_text(face = "bold", size = 10))
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999")
hpa_cols <- c(secreted = okabe[3], membrane = okabe[4], `secreted+membrane isoforms` = okabe[5], intracellular = okabe[7], unannotated = "white")

strip_fig <- function(top, title, ylab_col = "protein") {
  top <- copy(top); top[, lab := factor(get(ylab_col), rev(get(ylab_col)))]
  long <- melt(top[, .(lab, EMT = h_emt, `TNFA/NF-kB` = h_tnfa_nfkb, Apoptosis = h_apoptosis, TNFRSF = tnfrsf, p53 = h_p53)],
               id.vars = "lab", variable.name = "set", value.name = "member")
  pa <- ggplot(top, aes(z_stat, lab)) + geom_col(fill = okabe[1], width = 0.7) +
    geom_errorbarh(aes(xmin = log(LCI) / log(HR) * z_stat, xmax = log(HCI) / log(HR) * z_stat), height = 0, linewidth = 0.3, colour = "grey30") +
    labs(x = "Cox z (95 % CI of log HR rescaled to z)", y = NULL, title = title) + th
  pb <- ggplot(top, aes(weight, lab)) + geom_col(fill = okabe[2], width = 0.7) + labs(x = "Factor 1 weight", y = NULL) + th + theme(axis.text.y = element_blank())
  pc <- ggplot(long, aes(set, lab, fill = member)) + geom_tile(colour = "white") +
    scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey92"), guide = "none") + labs(x = NULL, y = NULL) + th +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  pd <- ggplot(top, aes("class", lab, fill = hpa_class)) + geom_tile(colour = "white") + scale_fill_manual(values = hpa_cols, name = "HPA class") +
    labs(x = NULL, y = NULL) + th + theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
  pa + pb + pc + pd + plot_layout(widths = c(3, 1.6, 2, 0.5), guides = "collect") & theme(legend.position = "bottom")
}

# ---- Fig 5: composite rank -------------------------------------------------------------------------------------------
top <- tab[!is.na(composite_rank)][order(composite_rank, -z_stat)][1:40]
top[, ylab := sprintf("%2d  %s  (w %d, z %d)", composite_order, protein, loading_rank, z_stat_rank)]
p5 <- strip_fig(top, "Top 40 by composite rank (loading rank + Cox-z rank; 504 proteins)\nLabel: rank  protein  (w = loading rank, z = Cox-z rank)", "ylab")
ggsave(file.path(root, "figures/fig5_top_composite.pdf"), p5, width = 7.6, height = 7.5)
ggsave(file.path(root, "figures/fig5_top_composite.png"), p5, width = 7.6, height = 7.5, dpi = 200)

# ---- Fig 6: lower-left --------------------------------------------------------------------------------------------------
tab[, x := fifelse(has_weight, weight, 0)]
tab[, ll_cat := fifelse(LL_primary, "weight < 0, z < -1.96 (n = 17)", fifelse(LL_file, "weight < 0, z < 0 (n = 84)",
               fifelse(protective & has_weight, "protective, weight > 0 (n = 12)", fifelse(protective, "protective, no weight (n = 46)", "other"))))]
tab[, ll_cat := factor(ll_cat, c("weight < 0, z < -1.96 (n = 17)", "weight < 0, z < 0 (n = 84)", "protective, weight > 0 (n = 12)", "protective, no weight (n = 46)", "other"))]
lab6 <- tab[LL_primary | (protective & has_weight & weight > 0 & z_stat < -2.3) | (protective & !has_weight & z_stat < -3.5), protein]
p6a <- ggplot(tab[z_stat < 3], aes(x, z_stat)) +
  geom_hline(yintercept = c(0, -1.96), linetype = c("solid", "dashed"), colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(data = tab[z_stat < 3 & ll_cat == "other"], colour = "grey80", size = 0.7, alpha = 0.6) +
  geom_point(data = tab[ll_cat != "other"], aes(colour = ll_cat, shape = has_weight), size = 1.7) +
  scale_colour_manual(values = okabe[c(2, 5, 1, 7)], name = NULL, limits = levels(tab$ll_cat)[1:4]) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), labels = c(`TRUE` = "has weight", `FALSE` = "no weight (drawn at x = 0)"), name = NULL) +
  geom_text_repel(data = tab[protein %in% lab6], aes(label = protein), size = 2.2, max.overlaps = 60, segment.size = 0.2, box.padding = 0.15, min.segment.length = 0.1) +
  coord_cartesian(xlim = c(-1.1, 1.0), ylim = c(-7, 3)) +
  labs(x = "Factor 1 tubular-mean weight (0 = no weight)", y = "Cox z, incident CKD (UK Biobank)", title = "Lower-left: plasma level falls with incident CKD") +
  th + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "cm")) + guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
# ORA dot plot: top terms per set vs panel
pick <- rbind(ora[set == "LL_primary" & universe == "U1_panel"][order(p.adjust)][1:6],
              ora[set == "LL_file" & universe == "U1_panel"][order(p.adjust)][1:6],
              ora[set == "LL_median" & universe == "U1_panel" & collection != "GO_BP"][order(p.adjust)][1:4],
              ora[set == "LL_median" & universe == "U1_panel" & collection == "GO_BP"][order(p.adjust)][1:6])
pick[, term := gsub("_", " ", sub("^(GOBP|GOCC|HALLMARK|REACTOME|KEGG)_", "", ID))]
pick[, term := paste0(tolower(substr(term, 1, 1)), tolower(substr(term, 2, 60)))]
pick[, term := sprintf("%s [%s] %d/%d", term, collection, k, M)]
pick[, set := factor(set, c("LL_primary", "LL_file", "LL_median"), c("weight < 0 & z < -1.96 (17)", "weight < 0 & z < 0 (84)", "below medians (252)"))]
pick[, term := factor(term, rev(unique(term)))]
p6b <- ggplot(pick, aes(OR, term, colour = set)) + geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = OR_lo, xmax = pmin(OR_hi, 60)), height = 0, linewidth = 0.3) + geom_point(aes(size = -log10(p.adjust))) +
  geom_point(data = pick[p.adjust >= 0.05], shape = 4, colour = "black", size = 1.5) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) + scale_colour_manual(values = okabe[c(2, 5, 6)], name = NULL) + scale_size_continuous(range = c(1, 4), name = "-log10 BH p") +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Odds ratio vs Olink panel (95 % CI; capped at 60); × = BH >= 0.05", y = NULL, title = "Over-representation vs Olink panel") +
  th + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, size = 7), axis.text.y = element_text(size = 6.5))
p6 <- free(p6a) / p6b + plot_layout(heights = c(1.1, 1.3))
ggsave(file.path(root, "figures/fig6_lowerleft.pdf"), p6, width = 7.2, height = 10)
ggsave(file.path(root, "figures/fig6_lowerleft.png"), p6, width = 7.2, height = 10, dpi = 200)

# ---- Fig 7: fgsea tails -------------------------------------------------------------------------------------------------
fz <- fg[grepl("Cox", ranking)]; fw <- fg[grepl("weight", ranking)]
sel <- function(d, n_neg = 12, n_pos = 10) rbind(d[NES < 0][order(padj)][1:n_neg], d[NES > 0][order(padj)][1:n_pos])
mk <- function(d, title) { d <- copy(d)[!is.na(pathway)]
  d[, term := gsub("_", " ", sub("^(GOBP|GOCC|HALLMARK|REACTOME|KEGG)_", "", pathway))]; d[, term := paste0(substr(term, 1, 1), tolower(substr(term, 2, 55)))]
  d[, term := sprintf("%s [%s, %d]", term, collection, size)]; d <- d[order(NES)]; d[, term := factor(term, unique(term))]; d[, sig := padj < 0.05]
  ggplot(d, aes(NES, term, fill = NES > 0, alpha = sig)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(`TRUE` = okabe[2], `FALSE` = okabe[1]), guide = "none") + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), labels = c(`TRUE` = "BH < 0.05", `FALSE` = "BH >= 0.05"), name = NULL) +
    labs(x = "Normalised enrichment score (fgsea)", y = NULL, title = title) + th + theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom", plot.title.position = "plot") }
# Full-collection version kept as fig7_fgsea_tails_allcollections; per-collection versions (user request 9 Oct 2026):
# fig7_fgsea_hallmark and fig7_fgsea_reactome, each with the Cox-z panel (A) and the weight panel (B).
for (cn in c("Hallmark", "Reactome")) {
  fzc <- fz[collection == cn & !is.na(NES)]; fwc <- fw[collection == cn & !is.na(NES)]
  k_z <- min(10, sum(fzc$NES < 0)); k_w <- min(10, sum(fwc$NES < 0))
  neg_z <- fzc[NES < 0][order(padj)][seq_len(k_z)]; pos_z <- fzc[NES > 0][order(padj)][1:10]
  neg_w <- fwc[NES < 0][order(padj)][seq_len(k_w)]; pos_w <- fwc[NES > 0][order(padj)][1:10]
  nz <- fzc[NES < 0 & padj < 0.05, .N]; nw <- fwc[NES < 0 & padj < 0.05, .N]
  pa <- mk(rbind(neg_z, pos_z), sprintf("A  %s, ranked by Cox z for incident CKD (all 2,923 proteins)\n    negative = plasma level falls with incident CKD; %d negative term(s) at BH < 0.05", cn, nz))
  pb <- mk(rbind(neg_w, pos_w), sprintf("B  %s, ranked by Factor 1 tubular weight (984 matched proteins)\n    negative = low weight in the tubular factor; %d negative term(s) at BH < 0.05", cn, nw))
  pp <- pa / pb
  ggsave(file.path(root, sprintf("figures/fig7_fgsea_%s.pdf", tolower(cn))), pp, width = 7.2, height = 9.5)
  ggsave(file.path(root, sprintf("figures/fig7_fgsea_%s.png", tolower(cn))), pp, width = 7.2, height = 9.5, dpi = 200)
}
cat("addendum figures written\n")
