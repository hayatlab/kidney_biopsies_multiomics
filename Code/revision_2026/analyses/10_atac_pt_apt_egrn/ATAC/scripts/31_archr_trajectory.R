#!/usr/bin/env Rscript
# ArchR trajectory on the scMEGA paired object, side by side with the DC1-rank Trajectory (Leon, 2026-09-11).
# scMEGA::AddTrajectory is a line-by-line port of ArchR::addTrajectory (ArchR's own function needs an ArchRProject /
# Arrow files, which do not exist for these samples): cells of each backbone group within the 90% nearest their
# group centroid (pre.filter.quantile 0.9) get an initial time from their distance to the next group's centroid; a
# smoothing spline (df 250, spar 1) is fitted through them per dimension; every cell is placed at its nearest spline
# point; the 10% of cells farthest from the curve get no value (post.filter.quantile 0.9); the rest are scaled to
# 100 * quantile, i.e. (0, 100].
# Leon's choices: backbone = RNA label of the pair's RNA cell, PT -> aPT; ArchR defaults; two spaces
#   glue  30-dim GLUE latent
#   dm    joint diffusion map DC1-DC2 (script 12; the map in panel a; only DC1-DC2 were saved)
# Each pair sits at the mean of its two cells' coordinates (scMEGA's paired object keeps the RNA member's only),
# the same way the DC1-rank Trajectory averages the two members.
# Out: tables/archr_trajectories{_smoke}.csv (per pair), scmega/outputs/archr_vs_dc1rank_summary{_smoke}.csv. The
#      comparison figure is drawn from the table by 35_archr_vs_dc1rank_figure.R. Human hg38 | seed 42 | SMOKE=1 -> the
#      smoke paired object.
source(file.path(Sys.getenv("KB_REPO"), "config", "kbpaths.R"))  # repository path config
options(Seurat.object.assay.version = "v3", stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(Signac); library(scMEGA) })
set.seed(42)
A   <- paste0(KB_ANALYSIS, "/aki_c5_exclusion_20261022/atac_rerun/ATAC")
SM  <- Sys.getenv("SMOKE") == "1"; TAG <- if (SM) "_smoke" else ""
IN  <- file.path(A, "scmega", paste0("inputs", TAG))
o   <- function(x, ext) file.path(A, "scmega", "outputs", paste0(x, TAG, ext))
obj <- readRDS(o("obj_pair", ".rds"))
stopifnot(all(obj$label %in% c("PT", "aPT")))
cat(sprintf("paired object: %d pairs (%d PT, %d aPT by RNA label) | scMEGA %s | smoke %s\n", ncol(obj),
            sum(obj$label == "PT"), sum(obj$label == "aPT"), as.character(packageVersion("scMEGA")), SM))

# pair coordinates = mean of the two members
glue <- as.matrix(read.csv(file.path(IN, "glue_latent.csv"), row.names = 1, check.names = FALSE))
ho <- read.csv(file.path(A, "tables", "HANDOFF_joint_pseudotime.csv"))
dm <- as.matrix(ho[, c("joint_dc1", "joint_dc2")]); rownames(dm) <- ho$cellid
for (M in list(glue, dm)) stopifnot(all(obj$rna_cellid %in% rownames(M)), all(obj$atac_cellid %in% rownames(M)))
pm <- function(M, key) {
  x <- (M[obj$rna_cellid, , drop = FALSE] + M[obj$atac_cellid, , drop = FALSE]) / 2
  dimnames(x) <- list(colnames(obj), paste0(key, seq_len(ncol(x))))
  CreateDimReducObject(embeddings = x, key = key, assay = "RNA")
}
obj[["gluepair"]] <- pm(glue, "GLUEPAIR_"); obj[["dmpair"]] <- pm(dm, "DMPAIR_")

sp <- c(glue = "gluepair", dm = "dmpair")
for (s in names(sp))
  obj <- AddTrajectory(obj, name = paste0("archr_", s), trajectory = c("PT", "aPT"), group.by = "label",
                       reduction = sp[[s]], dims = NULL, pre.filter.quantile = 0.9, post.filter.quantile = 0.9,
                       dof = 250, spar = 1, seed = 42)

d <- data.frame(pair = colnames(obj), rna_cellid = obj$rna_cellid, atac_cellid = obj$atac_cellid,
                rna_label = obj$label, atac_label = obj$atac_label, traj_dc1rank = obj$Trajectory,
                traj_archr_glue = obj$archr_glue, traj_archr_dm = obj$archr_dm)
write.csv(d, file.path(A, "tables", paste0("archr_trajectories", TAG, ".csv")), row.names = FALSE)

# Summary: coverage, agreement with the DC1 rank and PT vs aPT position. Spearman rho with a Fisher-z 95% CI using
# the Fieller-Hartley-Pearson variance 1.06 / (n - 3); pairs share donors, so the CI is optimistic.
rho_ci <- function(x, y) {
  k <- !is.na(x) & !is.na(y); r <- cor(x[k], y[k], method = "spearman"); n <- sum(k)
  h <- 1.96 * sqrt(1.06 / (n - 3))
  c(rho = r, lo = tanh(atanh(r) - h), hi = tanh(atanh(r) + h), n = n)
}
ax <- c(dc1rank = "DC1 rank", archr_glue = "ArchR, GLUE latent", archr_dm = "ArchR, diffusion map")
S <- do.call(rbind, lapply(names(ax), function(a) {
  v <- d[[paste0("traj_", a)]]; na <- is.na(v); r <- rho_ci(d$traj_dc1rank, v)
  data.frame(axis = a, pairs_with_value = sum(!na), pairs_without = sum(na), pct_without = 100 * mean(na),
             without_PT = sum(na & d$rna_label == "PT"), without_aPT = sum(na & d$rna_label == "aPT"),
             rho_vs_dc1rank = r[["rho"]], ci95_lo = r[["lo"]], ci95_hi = r[["hi"]],
             median_PT = median(v[d$rna_label == "PT"], na.rm = TRUE),
             median_aPT = median(v[d$rna_label == "aPT"], na.rm = TRUE))
}))
write.csv(S, o("archr_vs_dc1rank_summary", ".csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
r_gd <- rho_ci(d$traj_archr_glue, d$traj_archr_dm)
cat(sprintf("ArchR GLUE vs ArchR diffusion map: Spearman rho %.3f (95%% CI %.3f-%.3f; %d pairs with both)\n",
            r_gd[["rho"]], r_gd[["lo"]], r_gd[["hi"]], as.integer(r_gd[["n"]])))
cat("done\n")
