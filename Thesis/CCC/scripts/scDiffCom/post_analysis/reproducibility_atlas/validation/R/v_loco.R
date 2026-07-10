# LOCO validation — leave-one-cohort-out stability.

rank_cor_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2) return(NA_real_)
  stats::cor(x[ok], y[ok], method = "spearman")
}

run_v_loco <- function(ctx, validation_dir) {
  Xtilde <- ctx$Xtilde
  cohorts <- ctx$cohorts
  genes <- ctx$repro_df$gene
  full_scores <- ctx$repro_df$ReproScore
  names(full_scores) <- genes

  loco_results <- list()
  loco_scores <- list()

  for (drop_c in cohorts) {
    keep <- setdiff(cohorts, drop_c)
    Xsub <- Xtilde[keep]
    idx <- cohort_pairs(length(keep))
    repro <- compute_repro(Xsub, idx)
    scores <- repro$ReproScore
    names(scores) <- ctx$gene_universe
    loco_scores[[drop_c]] <- scores

    aligned <- scores[genes]
    rho <- rank_cor_spearman(full_scores, aligned)
    loco_results[[length(loco_results) + 1]] <- data.frame(
      dropped_cohort = drop_c,
      n_cohorts = length(keep),
      rank_rho_vs_full = rho,
      stringsAsFactors = FALSE
    )
  }
  loco_df <- do.call(rbind, loco_results)

  loco_names <- names(loco_scores)
  pair_rhos <- list()
  for (i in seq_along(loco_names)) {
    for (j in seq(i + 1L, length(loco_names))) {
      s1 <- as.numeric(loco_scores[[loco_names[i]]][genes])
      s2 <- as.numeric(loco_scores[[loco_names[j]]][genes])
      rho <- rank_cor_spearman(s1, s2)
      pair_rhos[[length(pair_rhos) + 1]] <- data.frame(
        loco_a = loco_names[i], loco_b = loco_names[j], rank_rho = rho,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(pair_rhos) > 0) {
    pair_df <- do.call(rbind, pair_rhos)
    loco_df <- rbind(loco_df, data.frame(
      dropped_cohort = paste0(pair_df$loco_a, "_vs_", pair_df$loco_b),
      n_cohorts = NA_integer_,
      rank_rho_vs_full = pair_df$rank_rho,
      stringsAsFactors = FALSE
    ))
  }
  readr::write_tsv(loco_df, file.path(validation_dir, "results", "loco_rank_cor.tsv"))

  top_k <- 50L
  top_sets <- lapply(loco_scores, function(s) {
    names(sort(s, decreasing = TRUE))[seq_len(min(top_k, length(s)))]
  })
  robust_core <- Reduce(intersect, top_sets)
  fragile <- character(0)
  for (drop_c in names(top_sets)) {
    others <- setdiff(names(top_sets), drop_c)
    only_with <- setdiff(top_sets[[drop_c]], unlist(top_sets[others]))
    fragile <- c(fragile, only_with)
  }
  fragile <- unique(fragile)

  robust_df <- data.frame(
    gene = robust_core,
    in_all_loco_top50 = TRUE,
    stringsAsFactors = FALSE
  )
  if (length(fragile) > 0) {
    fragile_df <- data.frame(
      gene = fragile,
      in_all_loco_top50 = FALSE,
      stringsAsFactors = FALSE
    )
    robust_df <- rbind(robust_df, fragile_df)
  }
  readr::write_tsv(robust_df, file.path(validation_dir, "results", "loco_robust_core.tsv"))

  min_rho <- min(loco_df$rank_rho_vs_full[loco_df$n_cohorts == 3], na.rm = TRUE)
  if (is.finite(min_rho) && min_rho > 0.8) {
    append_verdict(validation_dir, "loco_stability", "capstone", "PASS", round(min_rho, 4),
                   "all LOCO rho > 0.8", "Reproducibility not driven by one cohort")
  } else {
    worst <- loco_df[which.min(loco_df$rank_rho_vs_full[loco_df$n_cohorts == 3]), ]
    append_verdict(validation_dir, "loco_stability", "capstone", "FAIL", round(min_rho, 4),
                   "all LOCO rho > 0.8",
                   paste("Cohort-driven signal; worst drop:", worst$dropped_cohort))
  }

  ctx
}
