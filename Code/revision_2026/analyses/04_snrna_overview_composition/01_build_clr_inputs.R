#!/usr/bin/env Rscript
# Build the CLR composition table and the sample metadata that
# 02_ext2e_clr_compositions_boxplot.R consumes, from OUR new single-nucleus annotation.
# Leon, 2026-09-25.
#
# The saezlab workflow (github.com/saezlab/kidney_biopsies_MCFA) feeds its figure script a
# file `results/composition/clr~main_class~atlas.csv` built from the OLD 21 numbered clusters
# (0_PT, 1_TAL, ..., 18_aPT, ...). This rebuilds the same object from `level2` of
# adata_20260717 -- 16 cell types, 63 samples, 363,655 nuclei -- so the figure can be redrawn
# with nothing changed but the cell typing and the values underneath it.
#
# CLR, per sample, over all cell types of the chosen level:
#     clr(x_i) = log(x_i) - mean_j log(x_j)      on the compositional proportions
# No pseudocount is used or needed: every one of the 16 level2 types (and every level1 type)
# is present in every one of the 63 samples, so there is not a single structural zero. The
# old table carried a pseudocount floor that several samples shared; this one does not.
#
# Human GRCh38, HGNC symbols. Seed 42 (nothing stochastic here).  Env: scmega2.
# Usage: Rscript scripts/01_build_clr_inputs.R

source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
suppressPackageStartupMessages({ library(magrittr); library(dplyr) })
set.seed(42)

RA   <- KB_ANALYSIS
HERE <- file.path(RA, "Charlotte_help")
TAB  <- file.path(RA, "snrna_overview_20260831", "figures", "tables")
META <- file.path(RA, "vivien_takeover", "outputs", "metadata", "sample_metadata_corrected.csv")
OUT  <- file.path(HERE, "outputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

clr_from_proportions <- function(prop_df, level) {
  stopifnot("sample" %in% names(prop_df))
  m <- as.matrix(prop_df[, setdiff(names(prop_df), "sample")])
  rownames(m) <- prop_df$sample
  nz <- sum(m == 0)
  cat(sprintf("[%s] %d samples x %d cell types | structural zeros: %d\n",
              level, nrow(m), ncol(m), nz))
  if (nz > 0) stop("zeros present -- a pseudocount rule must be chosen before CLR")
  m <- m / rowSums(m)                                   # proportions, exactly
  lg <- log(m)
  clr <- lg - rowMeans(lg)
  out <- data.frame(Sample = rownames(clr), clr, check.names = FALSE)
  cat(sprintf("[%s] CLR range [%.2f, %.2f]; row sums %.2e (must be 0)\n",
              level, min(clr), max(clr), max(abs(rowSums(clr)))))
  out
}

nuc <- read.csv(file.path(TAB, "nuclei_per_sample.csv"))
for (level in c("level2", "level1")) {
  prop <- read.csv(file.path(TAB, sprintf("composition_%s_per_sample.csv", level)),
                   check.names = FALSE)
  clr <- clr_from_proportions(prop, level)
  f <- file.path(OUT, sprintf("clr~%s~new_annotation.csv", level))
  write.csv(clr, f, row.names = FALSE)
  cat(sprintf("[%s] wrote %s\n", level, f))
  # counts too, so anyone can redo the CLR under a different rule
  cnt <- prop
  n <- nuc$n_nuclei[match(cnt$sample, nuc$sample)]
  cnt[, setdiff(names(cnt), "sample")] <- round(cnt[, setdiff(names(cnt), "sample")] * n)
  write.csv(cnt, file.path(OUT, sprintf("counts~%s~new_annotation.csv", level)),
            row.names = FALSE)
}

# ---- metadata in the shape the figure script expects (its `patient_metadata_EDA.csv`) ----
m <- read.csv(META, check.names = FALSE)
pm <- m %>%
  transmute(sample_ID = sample,
            disease   = disease,
            classification = classification,
            egfr_at_time_of_biopsy = egfr_at_time_of_biopsy,
            # the saezlab script thresholds IFTA at 0.1 / 0.5, i.e. a FRACTION; ours is a
            # percentage (0-95), so it is divided by 100 and the script's cut-offs stand
            IFTA = IFTA / 100,
            proteinuria_around_biopsy = proteinuria_around_biopsy)
write.csv(pm, file.path(OUT, "patient_metadata_new_annotation.csv"), row.names = FALSE)
cat(sprintf("\n[meta] %d samples | disease: %s\n", nrow(pm),
            paste(sprintf("%s %d", names(table(pm$disease)), table(pm$disease)), collapse = ", ")))
cat(sprintf("[meta] eGFR NA %d | IFTA NA %d | proteinuria NA %d\n",
            sum(is.na(pm$egfr_at_time_of_biopsy)), sum(is.na(pm$IFTA)),
            sum(is.na(pm$proteinuria_around_biopsy))))
cat(sprintf("[done] %s\n", OUT))
