#!/usr/bin/env Rscript
# Addendum: (a) composite ranking (loading rank + z rank within the colleague's 504 weight>0 & z>0 proteins);
# (b) lower-left quadrant membership and over-representation; (c) pathways that fall with incident CKD
#     (fgsea on Cox z, negative tail, five MSigDB collections; ORA of protective predictors z < -1.96).
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressMessages({library(data.table); library(msigdbr); library(clusterProfiler); library(fgsea)})
set.seed(42)
root <- paste0(KB_ANALYSIS, "/factor1_reviewer_comment")
tab <- fread(file.path(root, "outputs/protein_table.csv"))
cs  <- fread(file.path(root, "raw/F1_Tubular_mean_results/Factor1/top_biomarker_candidates_z_score.csv"))

# ---- (a) composite rank ---------------------------------------------------------------------------------------------
stopifnot(all(cs$composite_rank == cs$loading_rank + cs$z_stat_rank),
          all(cs$loading_rank == rank(-cs$mean_loading)), all(cs$z_stat_rank == rank(-cs$z_stat)))
cs[, protein := toupper(variable)]
stopifnot(all(cs$protein %in% tab$protein), nrow(cs) == sum(tab$Q_file))
tab <- merge(tab, cs[, .(protein, loading_rank, z_stat_rank, composite_rank)], by = "protein", all.x = TRUE)
tab[, composite_order := frank(composite_rank, ties.method = "min", na.last = "keep")]
cat("composite rank = loading rank + z rank, both descending within the", nrow(cs), "weight>0 & z>0 proteins\n")
cat("Spearman composite_rank vs z rank:", round(cor(cs$composite_rank, cs$z_stat_rank, method = "spearman"), 3),
    "; vs loading rank:", round(cor(cs$composite_rank, cs$loading_rank, method = "spearman"), 3), "\n")
top40c <- tab[!is.na(composite_rank)][order(composite_rank, -z_stat)][1:40]
cat("\nTop 40 by composite rank:\n")
print(top40c[, .(composite_order, protein, weight = round(weight, 3), loading_rank, HR = round(HR, 2), z = round(z_stat, 1), z_stat_rank,
                 Q_primary, hpa_class, EMT = h_emt, NFkB = h_tnfa_nfkb, apop = h_apoptosis, tnfrsf)])
cat("\nTop 40 by composite that are NOT in top 40 by z:", paste(setdiff(top40c$protein, tab[Q_primary == TRUE][order(-z_stat)][1:40, protein]), collapse = ", "), "\n")
cat("Top 40 by composite in Q_primary (z > 1.96):", sum(top40c$Q_primary), "of 40\n")

# ---- (b) lower-left definitions ------------------------------------------------------------------------------------
med_w <- tab[has_weight == TRUE, median(weight)]; med_z <- tab[has_weight == TRUE, median(z_stat)]
tab[, LL_primary := has_weight & weight < 0 & z_stat < -1.96]      # mirror of Q_primary
tab[, LL_file    := has_weight & weight < 0 & z_stat < 0]          # mirror of the colleague's CSV
tab[, LL_median  := has_weight & weight < med_w & z_stat < med_z]  # manuscript's median split
tab[, protective := z_stat < -1.96]                                # whole panel
cat("\nLower-left counts: LL_primary", sum(tab$LL_primary), "; LL_file", sum(tab$LL_file), "; LL_median", sum(tab$LL_median),
    "(medians of matched: weight", round(med_w, 4), ", z", round(med_z, 3), ")\n")
cat("Protective predictors (z < -1.96) on the whole panel:", sum(tab$protective), " of which with weight:", sum(tab$protective & tab$has_weight),
    " with negative weight:", sum(tab$LL_primary), "\n")
cat("Risk predictors (z > 1.96):", sum(tab$ckd_predictor), " -> ratio risk:protective", round(sum(tab$ckd_predictor) / sum(tab$protective), 1), "\n")
cat("\nLL_primary proteins:\n")
print(tab[LL_primary == TRUE][order(z_stat), .(protein, weight = round(weight, 3), HR = round(HR, 2), LCI = round(LCI, 2), HCI = round(HCI, 2), z = round(z_stat, 1), hpa_class)])
cat("\nProtective predictors WITHOUT weight (strongest 25):\n")
print(tab[protective == TRUE & has_weight == FALSE][order(z_stat)][1:25, .(protein, HR = round(HR, 2), LCI = round(LCI, 2), HCI = round(HCI, 2), z = round(z_stat, 1), hpa_class)])
cat("\nProtective predictors WITH POSITIVE weight (lower-right, matched):\n")
print(tab[protective == TRUE & has_weight == TRUE & weight > 0][order(z_stat), .(protein, weight = round(weight, 3), HR = round(HR, 2), z = round(z_stat, 1), hpa_class)])

# ---- gene sets ----------------------------------------------------------------------------------------------------------
coll <- list(
  Hallmark = msigdbr(species = "Homo sapiens", category = "H"),
  Reactome = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
  KEGG     = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
  GO_BP    = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
  GO_CC    = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC"))
t2g <- lapply(coll, function(d) unique(as.data.frame(d)[, c("gs_name", "gene_symbol")]))

universes <- list(U1_panel = tab$protein, U2_matched = tab[has_weight == TRUE, protein], U4_protective = tab[protective == TRUE, protein])
sets <- list(LL_primary = tab[LL_primary == TRUE, protein], LL_file = tab[LL_file == TRUE, protein],
             LL_median = tab[LL_median == TRUE, protein], protective_all = tab[protective == TRUE, protein])
f <- function(s) as.numeric(strsplit(s, "/")[[1]])
res <- list()
for (sn in names(sets)) for (un in names(universes)) for (cn in names(coll)) {
  if (sn == "protective_all" && un == "U4_protective") next
  genes <- intersect(sets[[sn]], universes[[un]])
  e <- enricher(genes, universe = universes[[un]], TERM2GENE = t2g[[cn]],
                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
  if (is.null(e) || nrow(as.data.frame(e)) == 0) next
  d <- as.data.table(as.data.frame(e))
  d[, `:=`(set = sn, universe = un, collection = cn, n_query = length(genes), n_universe = length(universes[[un]]))]
  res[[length(res) + 1]] <- d
}
ora <- rbindlist(res)
ora[, c("k", "n") := .(sapply(GeneRatio, function(s) f(s)[1]), sapply(GeneRatio, function(s) f(s)[2]))]
ora[, c("M", "N") := .(sapply(BgRatio, function(s) f(s)[1]), sapply(BgRatio, function(s) f(s)[2]))]
ora[, `:=`(a = k, b = n - k, c = M - k, d = N - M - (n - k))]
ora[, OR := (a * d) / (b * c)]
ora[, OR_lo := exp(log(OR) - 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))]
ora[, OR_hi := exp(log(OR) + 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))]
ora[, fold := (k / n) / (M / N)]
setorder(ora, set, universe, collection, p.adjust)
fwrite(ora[, .(set, universe, collection, ID, k, n, M, N, fold, OR, OR_lo, OR_hi, pvalue, p.adjust, geneID)],
       file.path(root, "outputs/ora_lowerleft.csv"))
cat("\nORA rows:", nrow(ora), "\nSignificant (BH < 0.05) per set x universe x collection:\n")
print(dcast(ora[p.adjust < 0.05, .N, by = .(set, universe, collection)], set + universe ~ collection, value.var = "N", fill = 0))
for (sn in names(sets)) { cat("\n--- ", sn, " vs U1 panel, BH < 0.1 (top 25):\n")
  print(ora[set == sn & universe == "U1_panel" & p.adjust < 0.1][order(p.adjust)][1:25, .(collection, ID, k, n, M, N, OR = round(OR, 2), p.adjust = signif(p.adjust, 2), geneID)]) }
cat("\n--- LL_file vs U2 matched, BH < 0.1:\n")
print(ora[set == "LL_file" & universe == "U2_matched" & p.adjust < 0.1][order(p.adjust)][1:20, .(collection, ID, k, n, M, N, OR = round(OR, 2), p.adjust = signif(p.adjust, 2))])
cat("\n--- LL_primary vs U4 protective (does negative weight add anything beyond being protective?), BH < 0.2:\n")
print(ora[set == "LL_primary" & universe == "U4_protective" & p.adjust < 0.2][order(p.adjust)][1:20, .(collection, ID, k, n, M, N, OR = round(OR, 2), p.adjust = signif(p.adjust, 2))])

# ---- (c) fgsea on Cox z, all collections (both tails) and on weight ----------------------------------------------------
rz <- setNames(tab$z_stat, tab$protein)
tw <- tab[has_weight == TRUE]; rw <- setNames(tw$weight, tw$protein)
fg <- rbindlist(lapply(names(coll), function(cn) {
  gs <- split(t2g[[cn]]$gene_symbol, t2g[[cn]]$gs_name)
  a <- as.data.table(fgsea(gs, rz, minSize = 10, maxSize = 500, eps = 0)); a[, `:=`(collection = cn, ranking = "Cox z, all 2923 proteins")]
  b <- as.data.table(fgsea(gs, rw, minSize = 10, maxSize = 500, eps = 0)); b[, `:=`(collection = cn, ranking = "Factor 1 weight, 984 matched proteins")]
  rbind(a, b) }))
fg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
setorder(fg, ranking, collection, padj)
fwrite(fg, file.path(root, "outputs/fgsea_all_collections.csv"))
cat("\nfgsea on Cox z: terms with NEGATIVE NES and padj < 0.1 (these fall with incident CKD):\n")
print(fg[grepl("Cox", ranking) & NES < 0 & padj < 0.1][order(padj), .(collection, pathway, NES = round(NES, 2), padj = signif(padj, 2), size, leadingEdge = substr(leadingEdge, 1, 60))])
cat("\nfgsea on Cox z: most negative NES per collection regardless of padj (top 5 each):\n")
print(fg[grepl("Cox", ranking) & NES < 0][order(collection, NES), head(.SD, 5), by = collection][, .(collection, pathway, NES = round(NES, 2), padj = signif(padj, 2), size)])
cat("\nfgsea on Cox z: positive NES padj < 0.05 count per collection:\n"); print(fg[grepl("Cox", ranking) & NES > 0 & padj < 0.05, .N, by = collection])
cat("\nfgsea on weight: NEGATIVE NES padj < 0.1 (programmes with low tubular-factor weight):\n")
print(fg[grepl("weight", ranking) & NES < 0 & padj < 0.1][order(padj), .(collection, pathway, NES = round(NES, 2), padj = signif(padj, 2), size, leadingEdge = substr(leadingEdge, 1, 60))])

# ---- HPA / family Fisher for lower-left sets -------------------------------------------------------------------------------
flags <- c("hpa_secreted", "hpa_membrane", "hpa_intracellular", "hpa_plasma", "tnfrsf", "h_apoptosis", "h_tnfa_nfkb", "h_emt")
ct <- list()
for (sn in names(sets)) for (un in c("U1_panel", "U2_matched")) for (fl in flags) {
  u <- tab[protein %in% universes[[un]]]; inq <- u$protein %in% sets[[sn]]
  m <- matrix(c(sum(inq & u[[fl]]), sum(inq & !u[[fl]]), sum(!inq & u[[fl]]), sum(!inq & !u[[fl]])), 2); ft <- fisher.test(m)
  ct[[length(ct) + 1]] <- data.table(set = sn, universe = un, flag = fl, in_set = m[1, 1], set_size = sum(inq), in_rest = m[1, 2], rest_size = sum(!inq),
      frac_set = m[1, 1] / sum(inq), frac_rest = m[1, 2] / sum(!inq), OR = unname(ft$estimate), OR_lo = ft$conf.int[1], OR_hi = ft$conf.int[2], p = ft$p.value) }
ct <- rbindlist(ct); ct[, p_BH := p.adjust(p, "BH"), by = .(set, universe)]
fwrite(ct, file.path(root, "outputs/class_tests_lowerleft.csv"))
cat("\nClass tests, LL_file and protective_all vs panel:\n")
print(ct[set %in% c("LL_file", "protective_all") & universe == "U1_panel", .(set, flag, in_set, set_size, frac_set = round(frac_set, 3), frac_rest = round(frac_rest, 3), OR = round(OR, 2), OR_lo = round(OR_lo, 2), OR_hi = round(OR_hi, 2), p = signif(p, 2))])

fwrite(tab, file.path(root, "outputs/protein_table.csv"))
fwrite(top40c[, .(composite_order, protein, weight, loading_rank, HR, LCI, HCI, z_stat, z_stat_rank, composite_rank, Q_primary, hpa_class, h_emt, h_tnfa_nfkb, h_apoptosis, tnfrsf)],
       file.path(root, "outputs/top40_composite.csv"))
fwrite(tab[LL_file == TRUE][order(z_stat), .(protein, weight, HR, LCI, HCI, z_stat, LL_primary, hpa_class, h_emt, h_tnfa_nfkb, h_apoptosis)], file.path(root, "outputs/lowerleft_proteins.csv"))
cat("\nprotein_table.csv updated with composite ranks and LL flags\n")
