# Stage 5 — gene-shuffle nulls and global permutation test.

.shuffle_genes_within_cohort <- function(X) {
  lapply(X, function(m) m[sample(nrow(m)), , drop = FALSE])
}

run_stage_05_nulls <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  gene_universe <- atlas_env$gene_universe
  Xtilde <- atlas_env$Xtilde
  cohort_idx_pairs <- atlas_env$cohort_pairs
  obs <- atlas_env$repro

  n_perm <- cfg$n_perm
  message(sprintf("Stage 5: nulls (n_perm = %d)", n_perm))

  null_repro_scores <- matrix(NA_real_, length(gene_universe), n_perm)
  null_global <- numeric(n_perm)

  obs_global <- mean(obs$ReproScore, na.rm = TRUE) - 0.5

  run_one_perm <- function(b) {
    Xshuf <- .shuffle_genes_within_cohort(Xtilde)
    null_repro <- .compute_repro(Xshuf, cohort_idx_pairs)
    c(null_repro$ReproScore, mean(null_repro$ReproScore, na.rm = TRUE) - 0.5)
  }

  n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  if (.Platform$OS.type == "unix" && n_cores > 1 && n_perm >= 50) {
    message(sprintf("  using parallel::mclapply (%d cores)", n_cores))
    perm_results <- parallel::mclapply(seq_len(n_perm), run_one_perm, mc.cores = n_cores)
  } else {
    perm_results <- lapply(seq_len(n_perm), function(b) {
      if (b %% 100 == 0 || b == 1) message(sprintf("  permutation %d / %d", b, n_perm))
      run_one_perm(b)
    })
  }

  for (b in seq_len(n_perm)) {
    null_repro_scores[, b] <- perm_results[[b]][seq_len(length(gene_universe))]
    null_global[b] <- perm_results[[b]][length(gene_universe) + 1]
  }

  rownames(null_repro_scores) <- gene_universe
  p_emp <- rowMeans(null_repro_scores >= obs$ReproScore, na.rm = TRUE)
  shuffle_fdr <- p.adjust(p_emp, method = "BH")

  repro_df <- atlas_env$repro_df
  repro_df$shuffle_p <- p_emp
  repro_df$shuffle_FDR <- shuffle_fdr
  repro_df$shuffle_FDR_low_power <- TRUE

  global_p <- mean(null_global >= obs_global)

  saveRDS(null_repro_scores, file.path(output_dir, "results", "null_reproscore_matrix.rds"))

  global_lines <- c(
    "=== Global permutation null ===",
    paste("timestamp:", Sys.time()),
    paste("n_perm:", n_perm),
    paste("obs_statistic (mean(ReproScore) - 0.5):", obs_global),
    paste("empirical_p_one_sided:", global_p),
    paste("atlas_provenance_note:", cfg$atlas_provenance_note),
    if (!is.null(cfg$perm_block)) paste("perm_block:", paste(cfg$perm_block, collapse = ",")) else "perm_block: null (no blocking)"
  )
  writeLines(global_lines, file.path(output_dir, "results", "global_null.txt"))

  readr::write_tsv(repro_df, file.path(output_dir, "results", "repro_scores_with_nulls.tsv"))

  atlas_env$repro_df <- repro_df
  atlas_env$null_repro_scores <- null_repro_scores
  atlas_env$global_null <- list(
    obs_statistic = obs_global,
    empirical_p = global_p,
    null_distribution = null_global
  )
  saveRDS(atlas_env, file.path(output_dir, "results", "stage05_atlas_env.rds"))
  message(sprintf("  global null p = %.4f", global_p))

  invisible(atlas_env)
}
