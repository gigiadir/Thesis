# Stage 4 — per-gene ReproScore via paired cross-cohort Spearman percentile.

.cohort_pairs <- function(n_cohorts) {
  combn(n_cohorts, 2, simplify = FALSE)
}

.cor_pair <- function(va, vb) {
  stats::cor(t(va), t(vb), method = "spearman", use = "pairwise.complete.obs")
}

.compute_repro <- function(Xc, cohort_idx_pairs) {
  G <- nrow(Xc[[1]])
  n_pairs <- length(cohort_idx_pairs)
  U <- matrix(NA_real_, G, n_pairs)
  cor_mats <- vector("list", n_pairs)
  s_self <- matrix(NA_real_, G, n_pairs)

  for (p in seq_len(n_pairs)) {
    a <- cohort_idx_pairs[[p]][1]
    b <- cohort_idx_pairs[[p]][2]
    C <- .cor_pair(Xc[[a]], Xc[[b]])
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

run_stage_04_reproscore <- function(atlas_env) {
  output_dir <- atlas_env$output_dir
  cohorts <- atlas_env$cohorts
  gene_universe <- atlas_env$gene_universe
  Xtilde <- atlas_env$Xtilde

  message("Stage 4: ReproScore")

  idx_pairs <- .cohort_pairs(length(cohorts))
  pair_labels <- vapply(idx_pairs, function(pr) {
    paste(cohorts[pr[1]], cohorts[pr[2]], sep = "_vs_")
  }, character(1))

  repro <- .compute_repro(Xtilde, idx_pairs)
  names(repro$cor_mats) <- pair_labels

  repro_df <- data.frame(
    gene = gene_universe,
    ReproScore = repro$ReproScore,
    R_self = repro$R_self,
    frac_pairs = repro$frac_pairs,
    stringsAsFactors = FALSE
  )

  if (isTRUE(atlas_env$cfg$eb_shrink)) {
    mu <- mean(repro_df$R_self, na.rm = TRUE)
    repro_df$R_self_shrunk <- 0.5 * repro_df$R_self + 0.5 * mu
  }

  readr::write_tsv(repro_df, file.path(output_dir, "results", "repro_scores.tsv"))
  saveRDS(repro$cor_mats, file.path(output_dir, "results", "cor_pair_matrices.rds"))

  atlas_env$repro <- repro
  atlas_env$repro_df <- repro_df
  atlas_env$cohort_pairs <- idx_pairs
  atlas_env$pair_labels <- pair_labels
  saveRDS(atlas_env, file.path(output_dir, "results", "stage04_atlas_env.rds"))
  message(sprintf("  ReproScore range: [%.3f, %.3f]", min(repro_df$ReproScore, na.rm = TRUE), max(repro_df$ReproScore, na.rm = TRUE)))

  invisible(atlas_env)
}
