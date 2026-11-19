#!/usr/bin/env Rscript
# ext2e_clr_compositions_boxplot.R  --  saezlab/kidney_biopsies_MCFA,
#   workflow/scripts/figures/ext2e_clr_compositions_boxplot.R
# redrawn on OUR new single-nucleus annotation (level2 of adata_20260717).  Leon, 2026-09-25.
#
# The plotting code below is the upstream script VERBATIM. Only the inputs and five things
# that the new annotation makes impossible to keep are changed; each is marked ## CHANGED.
# A sixth change (19 Nov 2026) brings the script in line with the panels printed in the final
# figures (see 6 below).
#
#   1. inputs -> outputs/clr~level2~new_annotation.csv (16 cell types, 63 samples, no
#      pseudocount: there is not one structural zero) and
#      outputs/patient_metadata_new_annotation.csv, built by 01_build_clr_inputs.R.
#   2. the `control_KPMP` block is removed. Upstream added the 11 external KPMP control
#      samples (the `KB*` ids) as their own disease group and used them as the reference of
#      every comparison. Those samples are not in the new object at all (63 samples, all
#      K-numbers), so the group does not exist.
#   3. therefore the reference of the disease comparisons is our own `control` (n = 7).
#   4. the `aPT -> aPT1` / `aPT_aTAL_aDCT -> aPT` rename is removed. Upstream needed it
#      because its `aPT` was the small old cluster 18 and the aPT-like compartment was the
#      combined cluster 4 `aPT_aTAL_aDCT`. In level2, `aPT` is aPT.
#   5. ylim. Upstream fixed it at (-2.5, 9.5) because its CLR carried a pseudocount floor
#      that pushed values far out. Ours spans [-3.0, 2.9], so that window would both squash
#      the panels and silently DROP the points below -2.5. The limit is therefore taken from
#      the data with headroom for the brackets, and is still shared across the three panels.
#   6. (19 Nov 2026) final panels, figure PDFs of 18 Nov 2026. The printed Extended Data
#      Fig. 2e-g were redrawn after 25 Sep with two differences, both reproduced here:
#      (a) the two control_AKI samples (K95, K98; two biopsies of one donor) are kept: their
#          own group "Control (AKI)" without a test in panel e, and in the eGFR group of
#          their eGFR in panel f (they have no IF/TA value, so panel g does not change);
#      (b) the eGFR group is the column `egfr_group` of Supplementary Data 1 of the paper,
#          which changes the group of 9 of the other 61 samples. The column is used as it
#          stands: for K5 and K53 the table gives an eGFR of 30.0 and the group "<30", so
#          thresholding `egfr_at_biopsy` would put them in 30-60 and not reproduce the panel.
#      With (a) and (b) all 30 printed p-values are reproduced (e 15, f 9, g 6). The samples
#      and statistics of the run of 25 Sep (control_AKI dropped, eGFR of
#      sample_metadata_corrected.csv) are obtained with EXT2E_VERSION=20260925. Groups are
#      ordered and named as printed (DN ... TIN, Control, Control (AKI)); the code that drew
#      the printed panels is not on this system, so styling may differ from it.
#
# The three panels, the facets (aPT / PT / Immune), theme_classic, the Set2 / Set1 palettes,
# geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5), the grey hline at 0, the 45-degree
# x labels, the t-test brackets via ggpubr::stat_compare_means and the patchwork assembly are
# untouched.
#
# Env: scmega2.  Usage: Rscript scripts/02_ext2e_clr_compositions_boxplot.R [level2|level1]

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
library(magrittr)
library(dplyr)
library(ggplot2)

HERE  <- paste0(KB_ANALYSIS, "/Charlotte_help")          ## CHANGED 1
LEVEL <- {a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "level2"}
data <- read.csv(file.path(HERE, "outputs", sprintf("clr~%s~new_annotation.csv", LEVEL)),
                 check.names = FALSE)
patient_metadata <- read.csv(file.path(HERE, "outputs", "patient_metadata_new_annotation.csv"))

## CHANGED 6: eGFR at biopsy of Supplementary Data 1; control_AKI samples kept
FINAL <- Sys.getenv("EXT2E_VERSION", "final") != "20260925"
if (FINAL) {
  SUPP <- Sys.getenv("KB_SUPP_DATA1", file.path(KB_ANALYSIS, "Supplementary Data 1.csv"))
  supp <- read.csv(SUPP, sep = ";", check.names = FALSE)
  i <- match(patient_metadata$sample_ID, supp$sample_id)
  stopifnot(!anyNA(i), all(supp$snRNA_seq[i] %in% c("True", "TRUE", TRUE, 1)),
            all(ifelse(patient_metadata$disease == "IgA", "IgAN", patient_metadata$disease) == supp$disease[i]),
            isTRUE(all.equal(patient_metadata$IFTA, supp$IFTA_pct[i] / 100)))   # condition and IF/TA agree in both tables
  cat(sprintf("[meta] eGFR at biopsy from Supplementary Data 1: %d of %d samples differ from the 25 Sep metadata; NA %d\n",
              sum(abs(patient_metadata$egfr_at_time_of_biopsy - supp$egfr_at_biopsy[i]) > 1e-9, na.rm = TRUE) +
                sum(is.na(patient_metadata$egfr_at_time_of_biopsy) != is.na(supp$egfr_at_biopsy[i])),
              nrow(patient_metadata), sum(is.na(supp$egfr_at_biopsy[i]))))
  patient_metadata$egfr_at_time_of_biopsy <- supp$egfr_at_biopsy[i]
  SUPP_GROUP <- c("\u226590" = ">= 90", "60-90" = "60-90", "30-60" = "30-60", "<30" = "< 30")
  stopifnot(all(supp$egfr_group[i] %in% names(SUPP_GROUP)))
  patient_metadata$egfr_group_supp <- unname(SUPP_GROUP[supp$egfr_group[i]])
} else patient_metadata$egfr_group_supp <- NA_character_
DROP <- if (FINAL) character(0) else "control_AKI"            # samples left out of all three panels
COND <- c("DN", "FSGS", "IgAN", "MN", "TIN", "Control", "Control (AKI)")   # order of panel e as printed
cond_label <- function(d) factor(dplyr::recode(d, control = "Control", control_AKI = "Control (AKI)"), levels = COND)

metadata <- patient_metadata %>%
  dplyr::select(sample_ID, disease, ) %>%
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  dplyr::distinct()
## CHANGED 2: no external KPMP controls in the new object -> no `control_KPMP` group
REFERENCE <- "Control"                                                          ## CHANGED 3 (named as printed, CHANGED 6)
DISEASES  <- c("DN", "FSGS", "IgAN", "MN", "TIN")

colnames(data) <- stringr::str_replace(colnames(data), "X[:digit:]+_", "")

to_plot1 <- data %>%
  tidyr::pivot_longer(!Sample, names_to = 'group', values_to = 'prop') %>%
  dplyr::left_join(metadata, by = c("Sample" = "sample_ID")) %>%
  dplyr::filter(!disease %in% DROP) %>%                                          ## CHANGED 6
  dplyr::mutate(disease = droplevels(cond_label(disease)))
stopifnot(!anyNA(to_plot1$disease))
## CHANGED 4: no aPT / aPT_aTAL_aDCT rename -- level2 `aPT` is aPT

## CHANGED 5: limits from the data, shared by all three panels (upstream: c(-2.5, 9.5))
YLIM <- c(floor(min(to_plot1$prop) * 2) / 2, ceiling(max(to_plot1$prop) * 2) / 2 + 3.0)
cat(sprintf("[data] %s: %d samples x %d cell types | CLR range [%.2f, %.2f] -> ylim [%.1f, %.1f]\n",
            LEVEL, dplyr::n_distinct(to_plot1$Sample), dplyr::n_distinct(to_plot1$group),
            min(to_plot1$prop), max(to_plot1$prop), YLIM[1], YLIM[2]))
cat(sprintf("[data] groups: %s\n", paste(sort(unique(to_plot1$disease)), collapse = ", ")))

(
  p2 <- to_plot1 %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    ggplot(aes(x = disease, y = prop, color = disease)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set2") +
    labs(x = element_blank(), y = 'CLR-transformed proportion', subtitle = "Cell type abundance by condition") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = lapply(DISEASES, function(d) c(REFERENCE, d)),
                               size = 3) +
    ylim(YLIM)
)

metadata <- patient_metadata %>%
  dplyr::select(sample_ID, disease, egfr_at_time_of_biopsy, egfr_group_supp, IFTA, proteinuria_around_biopsy, classification) %>%
  dplyr::mutate(disease = ifelse(disease == "IgA", "IgAN", disease)) %>%
  dplyr::distinct()

to_plot <- data %>%
  tidyr::pivot_longer(!Sample, names_to = 'group', values_to = 'prop') %>%
  dplyr::left_join(metadata, by = c("Sample" = "sample_ID")) %>%
  dplyr::filter(!disease %in% DROP) %>%                                          ## CHANGED 6
  dplyr::mutate(eGFR = ifelse(egfr_at_time_of_biopsy < 30,
                              "< 30",
                              ifelse(egfr_at_time_of_biopsy < 60,
                                     "30-60",
                                     ifelse(egfr_at_time_of_biopsy < 90,
                                            "60-90",
                                            ifelse(egfr_at_time_of_biopsy >= 90,
                                                   ">= 90",
                                                   NA)
                                     )
                              )
  )
  ) %>%
  dplyr::mutate(eGFR = if (FINAL) egfr_group_supp else eGFR) %>%                 ## CHANGED 6: group as given in Supplementary Data 1
  dplyr::mutate(eGFR = factor(eGFR, levels = c(">= 90", "60-90", "30-60", "< 30"))) %>%
  dplyr::mutate(IFTA_levels = ifelse(IFTA < 0.1,
                                     "< 10%",
                                     ifelse(IFTA < 0.5,
                                            "10-50%",
                                            ifelse(IFTA >= 0.5,
                                                   ">= 50%",
                                                   NA)))) %>%
  dplyr::mutate(IFTA_levels = factor(IFTA_levels, levels = c("< 10%", "10-50%", ">= 50%"))) %>%
  dplyr::mutate(PCR = ifelse(proteinuria_around_biopsy < 1,
                             "< 1 g/g",
                             ifelse(proteinuria_around_biopsy < 3,
                                    "1-3 g/g",
                                    ifelse(proteinuria_around_biopsy >= 3,
                                           ">= 3 g/g",
                                           NA)))) %>%
  dplyr::mutate(PCR = factor(PCR, levels = c("< 1 g/g", "1-3 g/g", ">= 3 g/g")))
## CHANGED 4 (again): no aPT rename

(
  p3 <- to_plot %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    tidyr::drop_na(eGFR) %>%
    ggplot(aes(x = eGFR,
                            y = prop, color = eGFR)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set1") +
    labs(x = element_blank(), y = 'CLR-transformed proportion', subtitle = "Cell type abundance by eGFR at biopsy (ml/min/1.73m2)") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = list(c(">= 90", "60-90"),
                                                  c(">= 90", "30-60"),
                                                  c(">= 90", "< 30")), size = 3)  +
    ylim(YLIM)
)

(
  p4 <- to_plot %>%
    dplyr::filter(group %in% c("aPT", "PT", "Immune")) %>%
    tidyr::drop_na(IFTA_levels) %>%
    ggplot(aes(x = IFTA_levels,
                            y = prop, color = IFTA_levels)) +
    geom_hline(yintercept = 0, colour = "lightgray") +
    facet_wrap(.~factor(group, levels = c("aPT", "PT", "Immune")), nrow=1) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) +
    theme_classic() + scale_colour_brewer(palette = "Set1") +
    labs(x = element_blank(), y = 'CLR-transformed proportion', subtitle = "Cell type abundance by IF/TA") +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggpubr::stat_compare_means(method = "t.test",
                               comparisons = list(c("< 10%", "10-50%"),
                                                  c("< 10%", ">= 50%")), size = 3)  +
    ylim(YLIM)
)

library(patchwork)
p <- p2+p3+p4

FIG <- Sys.getenv("EXT2E_FIGDIR", file.path(HERE, "figures"))  # release patch: the ED 2e-g entry script redirects it
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
for (ext in c("pdf", "png"))
  ggsave(file.path(FIG, sprintf("ext2e_clr_compositions_boxplot_%s.%s", LEVEL, ext)),
         p, width = 16, height = 5.2, dpi = 300)

# ---- the numbers behind the brackets, so nothing in the figure is untraceable ----------
suppressWarnings({
  stat_rows <- do.call(rbind, lapply(c("aPT", "PT", "Immune"), function(g) {
    d1 <- to_plot1 %>% dplyr::filter(group == g)
    do.call(rbind, lapply(DISEASES, function(dis) {
      a <- d1$prop[d1$disease == REFERENCE]; b <- d1$prop[d1$disease == dis]
      if (length(a) < 2 || length(b) < 2) return(NULL)
      tt <- t.test(b, a)
      data.frame(panel = "disease", cell_type = g, group1 = REFERENCE, group2 = dis,
                 n1 = length(a), n2 = length(b),
                 mean_diff = unname(diff(rev(tt$estimate))),
                 ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
                 t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
    }))
  }))
  for (spec in list(list("eGFR", c(">= 90", "60-90"), c(">= 90", "30-60"), c(">= 90", "< 30")),
                    list("IFTA_levels", c("< 10%", "10-50%"), c("< 10%", ">= 50%")))) {
    v <- spec[[1]]
    for (g in c("aPT", "PT", "Immune")) for (cmp in spec[-1]) {
      d1 <- to_plot %>% dplyr::filter(group == g)
      a <- d1$prop[!is.na(d1[[v]]) & d1[[v]] == cmp[1]]
      b <- d1$prop[!is.na(d1[[v]]) & d1[[v]] == cmp[2]]
      if (length(a) < 2 || length(b) < 2) next
      tt <- t.test(b, a)
      stat_rows <- rbind(stat_rows, data.frame(
        panel = v, cell_type = g, group1 = cmp[1], group2 = cmp[2],
        n1 = length(a), n2 = length(b), mean_diff = unname(diff(rev(tt$estimate))),
        ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2], t = unname(tt$statistic),
        df = unname(tt$parameter), p = tt$p.value))
    }
  }
})
stat_rows$p_BH_within_panel <- ave(stat_rows$p, stat_rows$panel, FUN = function(x) p.adjust(x, "BH"))
write.csv(stat_rows, file.path(Sys.getenv("EXT2E_TABLEDIR", file.path(HERE, "outputs")), sprintf("ext2e_stats_%s.csv", LEVEL)),
          row.names = FALSE)
print(stat_rows[order(stat_rows$panel, stat_rows$cell_type, stat_rows$p), ], row.names = FALSE, digits = 3)
cat(sprintf("\n[done] %s/ext2e_clr_compositions_boxplot_%s.pdf\n", FIG, LEVEL))
