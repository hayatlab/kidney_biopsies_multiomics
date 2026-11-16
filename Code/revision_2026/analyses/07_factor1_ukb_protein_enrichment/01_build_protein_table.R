#!/usr/bin/env Rscript
# Merge per-protein Cox results (UK Biobank, Olink Explore 3072) with the Factor 1 tubular-mean signature weights,
# define the quadrants, attach HPA protein classes and gene-family flags. Human, HGNC symbols. No randomness.
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressMessages({library(data.table)})
root <- paste0(KB_ANALYSIS, "/factor1_reviewer_comment")
raw  <- file.path(root, "raw/F1_Tubular_mean_results/Factor1")
cox  <- as.data.frame(readRDS(file.path(raw, "final_results_df.rds")))
cox$protein <- rownames(cox); setDT(cox)
bundle <- readRDS(file.path(raw, "signature_bundle.rds"))
sig <- as.data.table(bundle$signature); sig[, gene_up := toupper(gene)]
stopifnot(bundle$context$mode == "tubular_mean", !any(duplicated(sig$gene_up)), !any(duplicated(toupper(cox$protein))))

tab <- merge(cox, sig[, .(gene_up, weight)], by.x = "protein", by.y = "gene_up", all.x = TRUE)
tab[, has_weight := !is.na(weight)]
cat(sprintf("Olink proteins: %d; with Factor 1 weight: %d (%.1f%%); without: %d\n",
            nrow(tab), sum(tab$has_weight), 100*mean(tab$has_weight), sum(!tab$has_weight)))
cat(sprintf("Signature genes: %d; matched to a panel protein: %d\n", nrow(sig), sum(sig$gene_up %in% tab$protein)))

# Colleague's CSV must equal weight>0 & z>0
tb <- fread(file.path(raw, "top_biomarker_candidates_z_score.csv"))
stopifnot(setequal(toupper(tb$variable), tab[has_weight & weight > 0 & z_stat > 0, protein]))

med_w <- median(tab[has_weight == TRUE, weight]); med_z <- median(tab[has_weight == TRUE, z_stat])
cat(sprintf("medians among matched proteins: weight %.4f, z %.3f\n", med_w, med_z))
tab[, Q_primary := has_weight & weight > 0 & z_stat > 1.96]
tab[, Q_file    := has_weight & weight > 0 & z_stat > 0]
tab[, Q_median  := has_weight & weight > med_w & z_stat > med_z]
tab[, ckd_predictor := z_stat > 1.96]
tab[, quadrant := fifelse(!has_weight, "no weight",
                  fifelse(weight > 0 & z_stat > 1.96, "High loading, high risk",
                  fifelse(weight > 0, "High loading, low risk",
                  fifelse(z_stat > 1.96, "Low loading, high risk", "Low loading, low risk"))))]
print(tab[, .N, by = quadrant][order(-N)])

# HPA protein classes (downloaded 8 Oct 2026, see resources/HPA_download_date.txt)
hpa <- function(f) fread(file.path(root, "resources", f))$Gene
tab[, hpa_secreted      := protein %in% hpa("Predicted_secreted_proteins.tsv")]
tab[, hpa_membrane      := protein %in% hpa("Predicted_membrane_proteins.tsv")]
tab[, hpa_intracellular := protein %in% hpa("Predicted_intracellular_proteins.tsv")]
tab[, hpa_plasma        := protein %in% hpa("Plasma_proteins.tsv")]
tab[, hpa_any := hpa_secreted | hpa_membrane | hpa_intracellular]
cat(sprintf("Panel proteins with any HPA class annotation: %d / %d\n", sum(tab$hpa_any), nrow(tab)))
tab[, hpa_class := fifelse(hpa_secreted & !hpa_membrane, "secreted",
                   fifelse(hpa_membrane & !hpa_secreted, "membrane",
                   fifelse(hpa_secreted & hpa_membrane, "secreted+membrane isoforms",
                   fifelse(hpa_intracellular, "intracellular", "unannotated"))))]

# Gene families / curated sets (membership from MSigDB Hallmark, msigdbr 7.5.1)
suppressMessages(library(msigdbr))
H <- as.data.table(msigdbr(species = "Homo sapiens", category = "H"))
tab[, tnfrsf   := grepl("^TNFRSF", protein)]
tab[, tnfsf    := grepl("^TNFSF", protein)]
tab[, h_apoptosis := protein %in% H[gs_name == "HALLMARK_APOPTOSIS", gene_symbol]]
tab[, h_tnfa_nfkb := protein %in% H[gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB", gene_symbol]]
tab[, h_emt       := protein %in% H[gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol]]
tab[, h_p53       := protein %in% H[gs_name == "HALLMARK_P53_PATHWAY", gene_symbol]]

setorder(tab, -z_stat)
fwrite(tab, file.path(root, "outputs/protein_table.csv"))
cat("wrote outputs/protein_table.csv\n")
