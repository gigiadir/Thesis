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
  df <- df[df$CCI %in% J, , drop = FALSE]
  if (nrow(df) == 0) return(vec)
  agg <- aggregate(LOGFC ~ CCI, data = df, FUN = function(x) mean(x, na.rm = TRUE))
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

# Advisor centering (Step 2): subtract the cohort "mean gene" profile from every
# gene's CCI vector. Mean subtraction only (no SD scaling), computed over ALL
# genes present per CCI column. NAs are preserved; all-NA columns are no-ops.
subtract_cohort_mean_profile <- function(mat) {
  mu <- colMeans(mat, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0
  sweep(mat, 2L, mu, FUN = "-")
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

# ------------------------------------------------------------------------------
# Advisor R_g machinery (Steps 3-5).
#
# For each cohort pair we build the full gene x gene Spearman matrix on the
# centered, J-restricted tensor. Entry C[g, h] = Spearman(v_g^{c1}, v_h^{c2}).
#   - diagonal C[g, g]      -> same-gene cross-cohort concordance (feeds R_g)
#   - off-diagonal C[g, h]  -> different-gene concordance (feeds the null)
# Pairwise entries with fewer than `min_overlap` jointly finite CCIs are set NA
# so both the observed statistic and the cross-gene null share the same gating.
# ------------------------------------------------------------------------------

compute_pair_cor_matrices <- function(Xc, cohort_idx_pairs, min_overlap = 10) {
  P <- length(cohort_idx_pairs)
  mats <- vector("list", P)
  for (p in seq_len(P)) {
    a <- cohort_idx_pairs[[p]][1]
    b <- cohort_idx_pairs[[p]][2]
    Ma <- Xc[[a]]
    Mb <- Xc[[b]]
    C <- suppressWarnings(
      stats::cor(t(Ma), t(Mb), method = "spearman", use = "pairwise.complete.obs")
    )
    if (is.finite(min_overlap) && min_overlap > 0) {
      Fa <- matrix(as.numeric(is.finite(Ma)), nrow(Ma), ncol(Ma))
      Fb <- matrix(as.numeric(is.finite(Mb)), nrow(Mb), ncol(Mb))
      overlap <- Fa %*% t(Fb)
      C[overlap < min_overlap] <- NA_real_
    }
    dimnames(C) <- list(rownames(Ma), rownames(Mb))
    mats[[p]] <- C
  }
  mats
}

collapse_pairwise <- function(mat_G_by_P, aggregate = c("median", "mean")) {
  aggregate <- match.arg(aggregate)
  f <- if (aggregate == "median") stats::median else base::mean
  apply(mat_G_by_P, 1L, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else f(v)
  })
}

# R_g = median (default) of a gene's <=6 same-gene pairwise Spearman correlations.
compute_Rg_from_cormats <- function(cor_mats, aggregate = "median") {
  G <- nrow(cor_mats[[1]])
  s_self <- vapply(cor_mats, function(C) diag(C), numeric(G))
  if (is.null(dim(s_self))) s_self <- matrix(s_self, nrow = G)
  rownames(s_self) <- rownames(cor_mats[[1]])
  list(
    Rg = collapse_pairwise(s_self, aggregate),
    Rg_mean = collapse_pairwise(s_self, "mean"),
    pairwise_rho = s_self,
    n_pairs_computable = rowSums(is.finite(s_self))
  )
}

# Supplementary percentile statistic (self vs other genes) — kept for comparison.
compute_reproscore_from_cormats <- function(cor_mats) {
  G <- nrow(cor_mats[[1]])
  P <- length(cor_mats)
  U <- matrix(NA_real_, G, P)
  s_self <- matrix(NA_real_, G, P)
  for (p in seq_len(P)) {
    C <- cor_mats[[p]]
    s <- diag(C)
    s_self[, p] <- s
    U[, p] <- vapply(seq_len(G), function(g) {
      bg <- C[g, -g]
      bg <- bg[is.finite(bg)]
      if (!is.finite(s[g]) || !length(bg)) return(NA_real_)
      mean(s[g] > bg)
    }, numeric(1))
  }
  list(
    ReproScore = rowMeans(U, na.rm = TRUE),
    R_self     = rowMeans(s_self, na.rm = TRUE),
    frac_pairs = rowMeans(U > 0.95, na.rm = TRUE),
    U          = U
  )
}

# Draw a partner index g' != g for every gene (independent, with replacement).
sample_cross_gene_partner <- function(G) {
  partner <- sample.int(G, G, replace = TRUE)
  repeat {
    self <- which(partner == seq_len(G))
    if (!length(self)) break
    partner[self] <- sample.int(G, length(self), replace = TRUE)
  }
  partner
}

# One cross-gene null draw of the R_g vector: for each pair pick a random
# different-gene partner and aggregate the sampled off-diagonal concordances.
null_one_perm_Rg <- function(cor_mats, aggregate = "median") {
  G <- nrow(cor_mats[[1]])
  P <- length(cor_mats)
  null_pair <- matrix(NA_real_, G, P)
  for (p in seq_len(P)) {
    partner <- sample_cross_gene_partner(G)
    null_pair[, p] <- cor_mats[[p]][cbind(seq_len(G), partner)]
  }
  collapse_pairwise(null_pair, aggregate)
}
