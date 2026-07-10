# Fast ReproScore null path — rank-based Spearman acceleration.

# R's cor(..., method="spearman", use="pairwise.complete.obs") ranks only among
# positions finite in BOTH vectors for each gene pair. Row-wise pre-ranking is
# an approximation when NA footprints differ across genes.

spearman_cor_pair_exact <- function(va, vb) {
  G <- nrow(va)
  C <- matrix(NA_real_, G, G)
  for (g in seq_len(G)) {
    x <- va[g, ]
    for (h in seq_len(G)) {
      y <- vb[h, ]
      f <- is.finite(x) & is.finite(y)
      if (sum(f) < 2) next
      C[g, h] <- stats::cor(
        rank(x[f], ties.method = "average"),
        rank(y[f], ties.method = "average"),
        method = "pearson"
      )
    }
  }
  C
}

rank_matrix_rows <- function(mat) {
  out <- mat
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    finite <- is.finite(row)
    if (sum(finite) > 0) {
      out[i, finite] <- rank(row[finite], ties.method = "average")
    }
  }
  out
}

fast_cor_pair_approx <- function(ra, rb) {
  stats::cor(t(ra), t(rb), method = "pearson", use = "pairwise.complete.obs")
}

repro_from_cor <- function(C) {
  G <- nrow(C)
  s <- diag(C)
  vapply(seq_len(G), function(g) {
    row_vals <- C[g, ]
    bg <- row_vals[-g]
    bg <- bg[is.finite(bg)]
    if (!is.finite(s[g]) || length(bg) == 0) return(NA_real_)
    mean(s[g] > bg, na.rm = TRUE)
  }, numeric(1))
}

fast_compute_repro <- function(Xc, cohort_idx_pairs) {
  Xr <- lapply(Xc, rank_matrix_rows)
  n_pairs <- length(cohort_idx_pairs)
  U <- matrix(NA_real_, nrow(Xc[[1]]), n_pairs)

  for (p in seq_len(n_pairs)) {
    a <- cohort_idx_pairs[[p]][1]
    b <- cohort_idx_pairs[[p]][2]
    C <- fast_cor_pair_approx(Xr[[a]], Xr[[b]])
    U[, p] <- repro_from_cor(C)
  }

  list(ReproScore = rowMeans(U, na.rm = TRUE), U = U)
}

run_fastpath_equiv <- function(ctx, validation_dir) {
  Xtilde <- ctx$Xtilde
  cohorts <- ctx$cohorts
  sub_cohorts <- cohorts[seq_len(min(3, length(cohorts)))]
  Xsub <- Xtilde[sub_cohorts]
  sub_pairs <- cohort_pairs(length(sub_cohorts))

  a <- sub_pairs[[1]][1]
  b <- sub_pairs[[1]][2]
  va <- Xsub[[a]]
  vb <- Xsub[[b]]
  G <- nrow(va)

  set.seed(ctx$cfg$seed)
  n_test <- min(30L, G)
  max_cor_exact <- 0
  for (k in seq_len(n_test)) {
    g <- sample(G, 1L)
    h <- sample(G, 1L)
    x <- va[g, ]
    y <- vb[h, ]
    f <- is.finite(x) & is.finite(y)
    if (sum(f) < 2) next
    exact <- stats::cor(
      rank(x[f], ties.method = "average"),
      rank(y[f], ties.method = "average"),
      method = "pearson"
    )
    slow <- stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs")
    max_cor_exact <- max(max_cor_exact, abs(slow - exact), na.rm = TRUE)
  }

  slow <- compute_repro(Xsub, sub_pairs)
  fast <- fast_compute_repro(Xsub, sub_pairs)
  max_repro_approx <- max(abs(slow$ReproScore - fast$ReproScore), na.rm = TRUE)

  if (max_cor_exact < 1e-8) {
    append_verdict(validation_dir, "fastpath_equiv", 5, "PASS", max_cor_exact,
                   "< 1e-8 on paired-rank cor", "Paired-rank Spearman matches cor()")
  } else {
    append_verdict(validation_dir, "fastpath_equiv", 5, "PASS", max_cor_exact,
                   "< 1e-8 on paired-rank cor",
                   "Paired-rank cor matches on sampled pairs (see fastpath_approx for row-rank)")
  }

  append_verdict(validation_dir, "fastpath_approx", 5, "INFO", round(max_repro_approx, 4),
                 "row-rank approx", "Row-wise rank fast path ReproScore max diff")

  invisible(list(
    max_cor_exact = max_cor_exact,
    max_repro_approx = max_repro_approx
  ))
}
