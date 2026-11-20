#!/usr/bin/env Rscript
# Bulk RNA-seq of cultured human proximal tubule cells: DESeq2 differential expression, PCA, volcano plots, fgsea
# (Hallmark, GO BP, Reactome), g:Profiler GO enrichment and top-50 heat maps, per experiment and contrast.
# Script of Sidrah Maryam (received 20 Nov 2026), unchanged except for the two folder variables below, which can now be
# set through KB_BULK_DIR and KB_BULK_GMT_DIR (defaults = the locations of the published run on another cluster).
# Usage: Rscript Bulk.R <1|2|3>     See README.md for the experiments, inputs and what was (not) checked.

library(DESeq2)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(gprofiler2)
library(data.table)
library(readr)
library(matrixStats)

options(bitmapType = "cairo")  

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
experiment <- ifelse(length(args) >= 1, args[1], "1")


BASE_DIR    <- Sys.getenv("KB_BULK_DIR", unset = "/beegfs/data/KramannLab/RelA_Bulk_TS/Processed")
COUNTS_FILE <- file.path(BASE_DIR, "star_salmon", "salmon.merged.gene_counts.tsv")
GMT_DIR     <- Sys.getenv("KB_BULK_GMT_DIR", unset = "/beegfs/data/KramannLab/RelA_Bulk_TS/reference")
qval_cutoff <- 0.05
log2_cutoff <- 0.5
min_count_sum <- 10   

#   h.all.v2026.1.Hs.symbols.gmt   (Hallmark, used as-is)
#   c5.all.v2026.1.Hs.symbols.gmt  (filtered to GOBP_* only)
#   c2.all.v2026.1.Hs.symbols.gmt  (filtered to REACTOME_* only)

experiment_config <- list(
  "1" = list(
    outdir     = file.path(BASE_DIR, "experiment1/DEG_results"),
    conditions = list(
      "DMSO_Vehicle"      = c("TS390","TS391","TS392"),
      "DMSO_sAPP"         = c("TS393","TS394","TS395"),
      "Etoposide_Vehicle" = c("TS396","TS397","TS398"),
      "Etoposide_sAPP"    = c("TS399","TS400","TS401")
    ),
    contrasts = list(
      list(A = "DMSO_sAPP",         B = "DMSO_Vehicle"),
      list(A = "Etoposide_Vehicle", B = "DMSO_Vehicle"),
      list(A = "Etoposide_sAPP",    B = "DMSO_Vehicle"),
      list(A = "Etoposide_sAPP",    B = "Etoposide_Vehicle")
    )
  ),
  "2" = list(
    outdir     = file.path(BASE_DIR, "experiment2/DEG_results"),
    conditions = list(
      "Ca_rela_5_1" = c("TS402","TS403","TS404"),
      "Ca_rela_5_2" = c("TS405","TS406","TS407"),
      "Ca_rela_3_3" = c("TS408","TS409","TS410"),
      "CaNT"        = c("TS411","TS412","TS413")
    ),
    sample_overrides = list(
      "CaNT_ext" = c("TS411","TS412","TS413","TS414")
    ),
    contrasts = list(
      list(A = "Ca_rela_5_1", B = "CaNT"),
      list(A = "Ca_rela_5_2", B = "CaNT"),
      list(A = "Ca_rela_3_3", B = "CaNT_ext")
    )
  ),
  "3" = list(
    outdir     = file.path(BASE_DIR, "experiment3/DEG_results"),
    conditions = list(
      "CaNT1_rela"  = c("TS414","TS415","TS416"),
      "CaNT2_rela"  = c("TS417","TS418","TS419"),
      "Ca_rela_5_1" = c("TS420","TS421","TS422"),
      "Ca_rela_5_2" = c("TS423","TS424","TS425"),
      "Ca_rela_3_3" = c("TS426","TS427","TS428")
    ),
    contrasts = list(
      list(A = "Ca_rela_5_1", B = "CaNT1_rela"),
      list(A = "Ca_rela_5_2", B = "CaNT1_rela"),
      list(A = "Ca_rela_3_3", B = "CaNT1_rela"),
      list(A = "Ca_rela_5_1", B = "CaNT2_rela"),
      list(A = "Ca_rela_5_2", B = "CaNT2_rela"),
      list(A = "Ca_rela_3_3", B = "CaNT2_rela")
    )
  )
)

cfg <- experiment_config[[experiment]]
if (is.null(cfg)) stop("Unknown experiment number: ", experiment, " (expected 1, 2, or 3)")

outdir <- cfg$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


if (!file.exists(COUNTS_FILE)) {
  stop("Count matrix not found at: ", COUNTS_FILE,
       "\nCheck the star_salmon/ output directory for the exact filename.")
}

raw <- read.delim(COUNTS_FILE, header = TRUE, check.names = FALSE)
sample_cols <- setdiff(colnames(raw), c("gene_id", "gene_name"))

message("Sample columns found in count matrix (first 5): ",
        paste(head(sample_cols, 5), collapse = ", "), " ...")

raw_agg <- raw %>%
  filter(gene_name != "") %>%
  group_by(gene_name) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop")

dat_full <- as.data.frame(raw_agg[, sample_cols])
rownames(dat_full) <- raw_agg$gene_name

dat_full <- dat_full[rowSums(dat_full) > min_count_sum | rownames(dat_full) == "RELA", ]

sample_to_condition <- unlist(lapply(names(cfg$conditions), function(cond) {
  setNames(rep(cond, length(cfg$conditions[[cond]])), cfg$conditions[[cond]])
}))

# Match columns in the counts matrix to the config
missing_samples <- setdiff(names(sample_to_condition), colnames(dat_full))
if (length(missing_samples) > 0) {
  stop("These samples from the config are missing in the count matrix: ",
       paste(missing_samples, collapse = ", "),
       "\nActual columns available: ", paste(colnames(dat_full), collapse = ", "))
}

dat <- dat_full[, names(sample_to_condition)]
condition <- factor(sample_to_condition[colnames(dat)])
meta <- data.frame(sample = colnames(dat), condition = condition)

#DEseq

dds <- DESeqDataSetFromMatrix(countData = round(dat, 0), colData = meta, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

vsd <- vst(dds)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file.path(outdir, paste0("pca.experiment", experiment, ".pdf")), width = 10, height = 8)
print(
  ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = dds$sample), size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_light() +
    ggtitle(paste0("PCA -- Experiment ", experiment))
)
dev.off()


get_DEgene <- function(dds, A, B, folder) {
  contrast <- c("condition", A, B)
  res <- results(dds, contrast = contrast, alpha = qval_cutoff)
  pv <- res$pvalue[!is.na(res$pvalue) & res$baseMean > 10]
  pdf(file.path(folder, paste0("pvalue_hist.", A, ".vs.", B, ".pdf")), width = 5, height = 4)
  hist(pv, breaks = 50, col = "grey70", border = "white",
       main = paste0(A, " vs ", B), xlab = "Unadjusted p-value")
  dev.off()
  
  res_df <- as.data.frame(res) %>% rownames_to_column("gene_name")
  res_df <- res_df[!is.na(res_df$padj), ]
  # plot RELA
  if (!("RELA" %in% res_df$gene_name) && "RELA" %in% rownames(res)) {
    res_raw <- results(dds, contrast = contrast,
                       independentFiltering = FALSE, cooksCutoff = FALSE)
    rela_row <- as.data.frame(res_raw)["RELA", , drop = FALSE] %>%
      rownames_to_column("gene_name")
    if (!is.na(rela_row$padj[1])) {
      res_df <- bind_rows(res_df, rela_row)
    } else {
      message("  RELA has no plottable p-value for ", A, " vs ", B,
              " even with independent filtering / Cook's cutoff disabled -- point cannot be placed")
    }
  }
  
  top10_labels <- res_df %>%
    filter(padj < qval_cutoff & abs(log2FoldChange) > log2_cutoff) %>%
    arrange(padj) %>%
    slice_head(n = 10) %>%
    pull(gene_name)
  label_genes <- if ("RELA" %in% res_df$gene_name) {
    union(top10_labels, "RELA")
  } else {
    top10_labels
  }
  
  pdf(file.path(folder, paste0("volcano.", A, ".vs.", B, ".pdf")), width = 7, height = 7)
  p <- EnhancedVolcano(res_df, lab = res_df$gene_name,
                       selectLab = label_genes,
                       x = "log2FoldChange", y = "padj",
                       title = "", subtitle = "",
                       pointSize = 3.0, labSize = 4.0,
                       ylab = bquote(~-Log[10] ~ adjusted ~ italic(P)),
                       colAlpha = 1,
                       FCcutoff = log2_cutoff,
                       pCutoff = qval_cutoff,
                       legendPosition = "top",
                       legendLabSize = 8,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       max.overlaps = Inf)
  print(p)
  dev.off()
  
  write.table(res_df, file = file.path(folder, paste0("DESEQ2.", A, ".vs.", B, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return(res_df)
}

run_fgsea_collections <- function(ranks, folder, A, B) {
  gmt_files <- c(
    Hallmark = "h.all.v2026.1.Hs.symbols.gmt",
    GO_BP    = "c5.all.v2026.1.Hs.symbols.gmt",
    Reactome = "c2.all.v2026.1.Hs.symbols.gmt"
  )
  for (set_name in names(gmt_files)) {
    gmt_path <- file.path(GMT_DIR, gmt_files[set_name])
    if (!file.exists(gmt_path)) {
      warning("Missing gene set file, skipping ", set_name, ": ", gmt_path)
      next
    }
    pathways <- gmtPathways(gmt_path)
    if (set_name == "Reactome") pathways <- pathways[grepl("^REACTOME_", names(pathways))]
    if (set_name == "GO_BP")    pathways <- pathways[grepl("^GOBP_",     names(pathways))]
    message(set_name, ": ", length(pathways), " gene sets after filtering")

    fres <- fgseaMultilevel(pathways = pathways, stats = ranks)
    fres_tidy <- fres %>% as_tibble() %>% arrange(desc(NES))
    fres_out <- fres_tidy %>%
      mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))
    write.table(fres_out, file = file.path(folder, paste0("fgsea.", set_name, ".", A, ".vs.", B, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    plt_data <- if (set_name == "Hallmark") {
      fres_tidy %>% filter(padj < qval_cutoff)
    } else {
      fres_tidy %>% top_n(50, wt = abs(NES))
    }
    
    if (nrow(plt_data) == 0) {
      message("  No significant ", set_name, " pathways (padj < ", qval_cutoff, ") for ", A, " vs ", B, " -- skipping plot")
    } else {
      plt <- plt_data %>%
        ggplot(aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = padj < qval_cutoff)) +
        geom_text(aes(label = sprintf("NES %.2f, padj %.2g", NES, padj),
                      hjust = ifelse(NES >= 0, -0.05, 1.05)),
                  size = 2.8) +
        scale_y_continuous(expand = expansion(mult = 0.35)) +
        coord_flip() +
        labs(x = "Pathway", y = "Normalized Enrichment Score",
             title = paste0("GSEA ", set_name, " -- ", A, " vs ", B)) +
        theme_minimal()
      ggsave(file.path(folder, paste0("GSEA_", set_name, ".", A, ".vs.", B, ".pdf")),
             plot = plt, width = 15, height = 8, units = "in")
    }
  }
}

run_go_enrichment <- function(deg_df, folder, A, B) {
  deg_df_ord <- deg_df[order(deg_df$t, decreasing = TRUE), ]
  top500_up <- head(deg_df_ord, 500)$gene_name
  top500_dn <- tail(deg_df_ord, 500)$gene_name
  
  for (direction in c("up", "dn")) {
    genes <- if (direction == "up") top500_up else top500_dn
    go_res <- tryCatch(
      gost(query = genes, organism = "hsapiens", ordered_query = TRUE,
           significant = TRUE, correction_method = "g_SCS",
           domain_scope = "annotated", sources = NULL),
      error = function(e) NULL
    )
    if (is.null(go_res)) next
    res <- go_res$result
    res_bp <- res[res$term_size < 2000 & res$source == "GO:BP", ]
    res_bp <- res_bp[order(res_bp$p_value), ]
    top25 <- head(res_bp, 25)
    if (nrow(top25) == 0) next
    top25$log10_pvalue <- -log10(top25$p_value)
    
    plt <- ggplot(top25, aes(x = reorder(term_name, log10_pvalue), y = log10_pvalue)) +
      geom_bar(aes(fill = log10_pvalue), stat = "identity") +
      labs(title = paste0("Top 25 GO:BP (", direction, ") -- ", A, " vs ", B),
           x = "", y = "-log10(p-value)") +
      theme_minimal() +
      coord_flip()
    ggsave(file.path(folder, paste0("GO_BP.", direction, ".", A, ".vs.", B, ".pdf")),
           plot = plt, width = 10, height = 8, units = "in")
  }
}

# main
known_groups <- levels(meta$condition)

for (contrast in cfg$contrasts) {
  A <- contrast$A
  B <- contrast$B
  message("Running contrast: ", A, " vs ", B)
  
  folder <- file.path(outdir, paste0(A, ".vs.", B))
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  needs_custom <- !(A %in% known_groups) || !(B %in% known_groups)
  
  if (needs_custom) {
    get_group_samples <- function(grp) {
      if (grp %in% names(cfg$conditions)) return(cfg$conditions[[grp]])
      if (!is.null(cfg$sample_overrides) && grp %in% names(cfg$sample_overrides)) {
        return(cfg$sample_overrides[[grp]])
      }
      stop("Group '", grp, "' not found in conditions or sample_overrides for experiment ", experiment)
    }
    samplesA <- get_group_samples(A)
    samplesB <- get_group_samples(B)
    custom_samples <- c(samplesA, samplesB)
    
    missing_custom <- setdiff(custom_samples, colnames(dat_full))
    if (length(missing_custom) > 0) {
      stop("Samples for contrast ", A, " vs ", B, " missing from count matrix: ",
           paste(missing_custom, collapse = ", "))
    }
    
    custom_dat <- dat_full[, custom_samples]
    custom_condition <- factor(c(rep(A, length(samplesA)), rep(B, length(samplesB))),
                               levels = c(B, A))
    custom_meta <- data.frame(sample = custom_samples, condition = custom_condition)
    
    message("  -> custom contrast: fitting separate DESeq2 model on ",
            length(custom_samples), " samples (", paste(custom_samples, collapse = ", "), ")")
    
    contrast_dds <- DESeqDataSetFromMatrix(countData = round(custom_dat, 0),
                                           colData = custom_meta, design = ~condition)
    contrast_dds <- estimateSizeFactors(contrast_dds)
    contrast_dds <- estimateDispersions(contrast_dds)
    contrast_dds <- nbinomWaldTest(contrast_dds)
    contrast_meta <- custom_meta
  } else {
    contrast_dds <- dds
    contrast_meta <- meta
  }
  
  res_df <- get_DEgene(contrast_dds, A, B, folder)
  
  deg_df <- res_df %>%
    dplyr::rename(log2FC = log2FoldChange, adj.pvalue = padj, t = stat) %>%
    filter(!is.na(adj.pvalue)) %>%
    arrange(t)
  
  ranks <- deg_df %>%
    dplyr::select(gene_name, t) %>%
    filter(!is.na(t)) %>%
    distinct() %>%
    group_by(gene_name) %>%
    summarize(t = mean(t)) %>%
    deframe()
  
  run_fgsea_collections(ranks, folder, A, B)
  run_go_enrichment(deg_df, folder, A, B)
  
  express_mat <- t(counts(contrast_dds, normalized = TRUE))
  top50_df <- deg_df %>% arrange(adj.pvalue) %>% slice_head(n = 50) %>% arrange(desc(log2FC))
  top50_genes <- top50_df$gene_name
  samples_AB <- contrast_meta$sample[contrast_meta$condition %in% c(A, B)]
  condition_colors <- ifelse(contrast_meta$condition[contrast_meta$sample %in% samples_AB] == B, "black", "indianred3")

  hmcolors <- colorRampPalette(c("navyblue", "white", "firebrick3"))(256)
  
  pdf(file.path(folder, paste0("heatmap.top50.", A, ".vs.", B, ".pdf")), width = 12, height = 8)
  tryCatch({
    heatmap.2(t(express_mat[samples_AB, top50_genes]),
              scale = "row", density.info = "none", trace = "none",
              col = hmcolors, ColSideColors = condition_colors,
              Rowv = FALSE, dendrogram = "column",
              cexRow = 1.0, keysize = 1.25,
              main = paste0("DEG ", A, " vs ", B), margins = c(3, 10))
  }, error = function(e) message("Heatmap failed for ", A, " vs ", B, ": ", e$message))
  dev.off()
}

message("Done. Results in: ", outdir)