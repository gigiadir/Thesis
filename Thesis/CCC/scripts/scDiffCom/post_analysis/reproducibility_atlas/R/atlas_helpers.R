# Shared helpers for reproducibility atlas sections (not stage runners).

save.atlas.checkpoint <- function(atlas_env, stage_n,
                                  results_dir = get("results_dir", envir = .GlobalEnv)) {
  path <- file.path(results_dir, sprintf("stage%02d_atlas_env.rds", stage_n))
  saveRDS(atlas_env, path)
  message("Saved ", basename(path))
  invisible(path)
}

cci_union_per_cohort <- function(malignant_by_cohort) {
  lapply(malignant_by_cohort, function(gene_list) {
    sort(unique(unlist(lapply(gene_list, function(df) df$CCI))))
  })
}

gene_logfc_vector <- function(df, J) {
  vec <- setNames(rep(NA_real_, length(J)), J)
  if (nrow(df) == 0) return(vec)
  agg <- df %>%
    group_by(CCI) %>%
    summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop")
  vec[agg$CCI] <- agg$LOGFC
  vec
}

assert_tensor_alignment <- function(X, gene_universe, J, cohorts) {
  stopifnot(length(X) == length(cohorts))
  for (i in seq_along(cohorts)) {
    stopifnot(identical(rownames(X[[i]]), gene_universe))
    stopifnot(identical(colnames(X[[i]]), J))
    stopifnot(nrow(X[[i]]) == length(gene_universe))
    stopifnot(ncol(X[[i]]) == length(J))
  }
  if (length(cohorts) > 1) {
    for (i in seq(2, length(cohorts))) {
      stopifnot(identical(rownames(X[[1]]), rownames(X[[i]])))
      stopifnot(identical(colnames(X[[1]]), colnames(X[[i]])))
    }
  }
}

center_cohort_matrix <- function(mat, min_genes = 10, eps = 1e-8) {
  centered <- mat
  zero_sd <- character(0)
  for (j in seq_len(ncol(mat))) {
    vals <- mat[, j]
    finite <- is.finite(vals)
    if (sum(finite) < min_genes) next
    mu <- mean(vals[finite])
    sigma <- stats::sd(vals[finite])
    if (!is.finite(sigma) || sigma == 0) {
      zero_sd <- c(zero_sd, colnames(mat)[j])
      next
    }
    centered[finite, j] <- (vals[finite] - mu) / max(sigma, eps)
  }
  list(mat = centered, zero_sd = zero_sd)
}

cohort_pairs <- function(n_cohorts) {
  combn(n_cohorts, 2, simplify = FALSE)
}

cor_pair <- function(va, vb) {
  stats::cor(t(va), t(vb), method = "spearman", use = "pairwise.complete.obs")
}

compute_repro <- function(Xc, cohort_idx_pairs) {
  G <- nrow(Xc[[1]])
  n_pairs <- length(cohort_idx_pairs)
  U <- matrix(NA_real_, G, n_pairs)
  cor_mats <- vector("list", n_pairs)
  s_self <- matrix(NA_real_, G, n_pairs)

  for (p in seq_len(n_pairs)) {
    a <- cohort_idx_pairs[[p]][1]
    b <- cohort_idx_pairs[[p]][2]
    C <- cor_pair(Xc[[a]], Xc[[b]])
    cor_mats[[p]] <- C
    s <- diag(C)
    s_self[, p] <- s
    U[, p] <- vapply(seq_len(G), function(g) {
      row_vals <- C[g, ]
      bg <- row_vals[-g]
      bg <- bg[is.finite(bg)]
      if (!is.finite(s[g]) || length(bg) == 0) return(NA_real_)
      mean(s[g] > bg, na.rm = TRUE)
    }, numeric(1))
  }

  list(
    ReproScore = rowMeans(U, na.rm = TRUE),
    R_self     = rowMeans(s_self, na.rm = TRUE),
    frac_pairs = rowMeans(U > 0.95, na.rm = TRUE),
    U          = U,
    cor_mats   = cor_mats
  )
}

shuffle_genes_within_cohort <- function(X) {
  lapply(X, function(m) m[sample(nrow(m)), , drop = FALSE])
}

ensure_idr <- function() {
  if (!requireNamespace("idr", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      stop("Package 'idr' required. Install with: BiocManager::install('idr')")
    }
    BiocManager::install("idr", ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(idr))
}

gene_rank_values <- function(df, ccis, ranking = "signed") {
  vec <- setNames(rep(NA_real_, length(ccis)), ccis)
  if (nrow(df) == 0) return(vec)
  agg <- df %>%
    group_by(CCI) %>%
    summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop")
  vals <- agg$LOGFC
  if (ranking == "magnitude") vals <- abs(vals)
  vec[agg$CCI] <- vals
  vec
}

run_pairwise_idr <- function(va, vb, threshold = 0.05) {
  finite <- is.finite(va) & is.finite(vb)
  if (sum(finite) < 10) {
    return(list(pass = character(0), local_idr = numeric(0), ccis = character(0)))
  }
  ccis <- names(va)[finite]
  x <- va[finite]
  y <- vb[finite]
  fit <- tryCatch(
    idr::idr(x, y, mu = 0.5, sigma = 0.1, p = threshold),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(pass = character(0), local_idr = rep(NA_real_, length(ccis)), ccis = ccis))
  }
  local <- fit$idr[, 1]
  names(local) <- ccis
  pass <- names(local)[local <= threshold]
  list(pass = pass, local_idr = local, ccis = ccis)
}
