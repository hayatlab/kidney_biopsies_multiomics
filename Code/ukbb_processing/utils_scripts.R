---
  title: "Factor 2 Proteomic Signature Projection in the UK Biobank"
author: "Sergio Alejandro Gómez Ochoa"
date: "2025-09-26"
---

calculate_weighted_signature_score <- function(
    proteomics_data,
    signature_df,
    gene_col,
    score_col,
    factor_name,
    id_col          = "Sample.ID",
    anchors         = NULL,      # character vector of anchor genes/proteins (case-insensitive)
    winsor_pct      = 0.99,      # winsorize |w| at this percentile; set NULL to skip
    top_q           = NULL,      # keep top q fraction by |w|, e.g. 0.20; set NULL to skip
    do_l2           = TRUE,      # L2-normalize weights
    do_cosine       = FALSE,     # additionally L2-normalize each sample row (cosine projection)
    log2_transform  = FALSE,     # set TRUE only if input is not already log/NPX
    zscale_columns  = TRUE,      # recommended
    median_impute   = FALSE,     # if TRUE, median-impute before scaling; otherwise z-then NA->0
    verbose         = TRUE
) {
  if (verbose) message(sprintf("--- Calculating scores for: %s ---", factor_name))
  stopifnot(is.data.frame(proteomics_data), is.data.frame(signature_df))
  stopifnot(id_col %in% colnames(proteomics_data))
  stopifnot(all(c(gene_col, score_col) %in% colnames(signature_df)))
  
  ## ---- Build named weight vector (lowercase keys) ----
  w_raw <- setNames(signature_df[[score_col]], tolower(as.character(signature_df[[gene_col]])))
  w_raw <- w_raw[is.finite(w_raw)]
  
  ## ---- Build numeric expression matrix X (samples x features) ----
  pd <- proteomics_data
  pd[[id_col]] <- as.character(pd[[id_col]])
  
  # Select numeric columns except id_col
  feature_cols <- setdiff(colnames(pd), id_col)
  X <- as.matrix(pd[, feature_cols, drop = FALSE])
  storage.mode(X) <- "double"
  
  if (log2_transform) {
    if (verbose) message("Log2-transforming intensities.")
    X <- log2(X)
  }
  
  # Optional: median impute BEFORE scaling (not recommended for heavy censoring)
  if (isTRUE(median_impute)) {
    if (verbose) message("Median-imputing missing values per feature (pre-scaling).")
    med <- apply(X, 2, function(v) median(v, na.rm = TRUE))
    for (j in seq_len(ncol(X))) {
      nas <- is.na(X[, j])
      if (any(nas)) X[nas, j] <- med[j]
    }
  }
  
  # Column standardization, then NA -> 0 to keep dot-product geometry stable
  col_mu  <- rep(0, ncol(X))
  col_sd  <- rep(1, ncol(X))
  if (isTRUE(zscale_columns)) {
    if (verbose) message("Z-scaling features (columns).")
    col_mu <- colMeans(X, na.rm = TRUE)
    col_sd <- apply(X, 2, sd, na.rm = TRUE)
    col_sd[!is.finite(col_sd) | col_sd == 0] <- 1
    X <- sweep(X, 2, col_mu, "-")
    X <- sweep(X, 2, col_sd, "/")
  }
  # Set remaining NAs to 0 (feature mean after z-scaling)
  X[is.na(X)] <- 0
  
  rownames(X) <- pd[[id_col]]
  
  ## ---- Align features to weights (case-insensitive) ----
  cn_lower <- tolower(colnames(X))
  common   <- intersect(names(w_raw), cn_lower)
  if (length(common) == 0) stop("No common proteins found between signature and proteomics.")
  
  Xc <- X[, match(common, cn_lower), drop = FALSE]
  colnames(Xc) <- common
  w <- w_raw[common]
  
  # Collapse potential duplicate columns
  if (any(duplicated(colnames(Xc)))) {
    grp <- split(seq_along(colnames(Xc)), colnames(Xc))
    Xc  <- do.call(cbind, lapply(grp, function(ii) {
      if (length(ii) == 1) Xc[, ii, drop = FALSE]
      else matrix(rowMeans(Xc[, ii, drop = FALSE]), ncol = 1)
    }))
    colnames(Xc) <- names(grp)
    # Realign weights to collapsed columns
    w <- w[colnames(Xc)]
  }
  stopifnot(identical(names(w), colnames(Xc)))
  
  ## ---- Weight regularization (winsorize + optional top-q threshold) ----
  winsor_cap <- NA_real_
  if (!is.null(winsor_pct)) {
    winsor_cap <- as.numeric(stats::quantile(abs(w), probs = winsor_pct, na.rm = TRUE))
    w <- pmax(pmin(w, winsor_cap), -winsor_cap)
  }
  if (!is.null(top_q)) {
    stopifnot(top_q > 0 && top_q <= 1)
    thr <- as.numeric(stats::quantile(abs(w), probs = 1 - top_q, na.rm = TRUE))
    w[abs(w) < thr] <- 0
    # Drop all-zero columns to avoid wasted compute
    keep <- which(w != 0)
    Xc   <- Xc[, keep, drop = FALSE]
    w    <- w[keep]
  }
  
  ## ---- Optional cosine normalization of rows (after alignment) ----
  if (isTRUE(do_cosine)) {
    row_norm <- sqrt(rowSums(Xc^2))
    row_norm[!is.finite(row_norm) | row_norm == 0] <- 1
    Xc <- Xc / row_norm
  }
  
  ## ---- L2-normalize weights ----
  l2_norm <- sqrt(sum(w^2))
  if (isTRUE(do_l2) && l2_norm > 0) w <- w / l2_norm
  
  ## ---- Optional anchor-based orientation (flip sign if needed) ----
  flipped <- FALSE
  anchor_stat <- NA_real_
  anchors_used <- character(0)
  if (!is.null(anchors) && length(anchors) > 0) {
    anc <- intersect(tolower(anchors), colnames(Xc))
    if (length(anc) > 0) {
      anchors_used <- anc
      s_raw <- as.numeric(Xc %*% w)
      A     <- rowMeans(Xc[, anc, drop = FALSE])
      # Prefer Pearson; if NA, fallback to Spearman
      r <- suppressWarnings(stats::cor(s_raw, A, use = "pairwise.complete.obs"))
      if (!is.finite(r)) {
        r <- suppressWarnings(stats::cor(s_raw, A, use = "pairwise.complete.obs", method = "spearman"))
      }
      if (!is.finite(r)) r <- 0
      anchor_stat <- r
      if (r < 0) { w <- -w; flipped <- TRUE }
    } else if (verbose) {
      message("Anchors provided, but none found in data; no orientation applied.")
    }
  }
  
  ## ---- Scores ----
  s <- as.numeric(Xc %*% w)
  names(s) <- rownames(Xc)
  
  out <- data.frame(
    Participant.ID = rownames(Xc),
    score          = s[rownames(Xc)],
    score_z        = as.numeric(scale(s[rownames(Xc)])),
    check.names    = FALSE
  )
  colnames(out)[2:3] <- c(paste0(factor_name, "_score"),
                          paste0(factor_name, "_score_z"))
  
  ## ---- Metadata for auditability ----
  attr(out, "weights_used")   <- w
  attr(out, "n_common")       <- length(w)
  attr(out, "winsor_cap")     <- winsor_cap
  attr(out, "l2_norm")        <- l2_norm
  attr(out, "flipped")        <- flipped
  attr(out, "anchor_stat")    <- anchor_stat
  attr(out, "anchors_used")   <- anchors_used
  attr(out, "col_scaling_mu") <- setNames(col_mu, feature_cols)
  attr(out, "col_scaling_sd") <- setNames(col_sd, feature_cols)
  attr(out, "top_q")          <- top_q
  attr(out, "do_cosine")      <- do_cosine
  
  if (verbose) {
    message(sprintf("Features used: %d; flipped: %s; anchor r: %s",
                    length(w), flipped, ifelse(is.na(anchor_stat), "NA", round(anchor_stat, 3))))
  }
  
  out
}




find_medication <- function(data, med_cols, drug_list) {
  all_meds <- apply(data[, med_cols], 1, paste, collapse = " ")
  pattern <- paste0("\\b", drug_list, "\\b", collapse = "|")
  as.integer(str_detect(tolower(all_meds), pattern))
}


# Helper function to find the minimum date, ignoring NA values
pmin_date <- function(...) {
  pmin(..., na.rm = TRUE)
}

# Helper function to parse dates, handling potential errors
to_Date <- function(x) {
  as.Date(x, format = "%Y-%m-%d")
}

