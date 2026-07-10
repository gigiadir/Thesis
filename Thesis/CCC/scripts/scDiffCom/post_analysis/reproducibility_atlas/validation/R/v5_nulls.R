# Stage 5 validation — nulls.

run_v5_nulls <- function(ctx, validation_dir) {
  run_fastpath_equiv(ctx, validation_dir)

  if (!isTRUE(ctx$stage5_complete)) {
    n_done <- if (!is.null(ctx$null_checkpoint)) ctx$null_checkpoint$completed else 0
    for (cid in c("null_centered", "na_exchangeable", "pval_dist", "global_test")) {
      append_verdict(validation_dir, cid, 5, "INFO", n_done,
                     paste0("n_perm=", ctx$cfg$n_perm),
                     "Stage 5 incomplete; skipping until nulls finish")
    }
    return(ctx)
  }

  null_repro <- ctx$null_repro
  if (is.null(null_repro) && !is.null(ctx$stage05_env)) {
    null_repro <- ctx$stage05_env$null_repro_scores
  }
  if (is.null(null_repro)) {
    stop("Stage 5 marked complete but null_reproscore_matrix missing")
  }

  null_mean <- mean(null_repro, na.rm = TRUE)
  null_center_df <- data.frame(
    null_mean = null_mean,
    null_sd = stats::sd(null_repro, na.rm = TRUE),
    n_perm = ncol(null_repro),
    n_genes = nrow(null_repro)
  )
  readr::write_tsv(null_center_df, file.path(validation_dir, "results", "null_center.tsv"))

  png(file.path(validation_dir, "results", "null_center_hist.png"),
      width = 600, height = 500)
  hist(as.numeric(null_repro), breaks = 50, main = "Null ReproScore distribution",
       xlab = "Null ReproScore", col = "grey70", border = "white")
  abline(v = 0.5, col = "red", lwd = 2)
  abline(v = null_mean, col = "blue", lwd = 2, lty = 2)
  dev.off()

  if (null_mean >= 0.48 && null_mean <= 0.52) {
    append_verdict(validation_dir, "null_centered", 5, "PASS", round(null_mean, 4),
                   "[0.48, 0.52]", "Null well-formed, centered at 0.5")
  } else {
    append_verdict(validation_dir, "null_centered", 5, "FAIL", round(null_mean, 4),
                   "[0.48, 0.52]", "Shuffle or score bug; check gene shuffle independence")
  }

  scheme_a_mean <- null_mean
  n_scheme_b <- 5L
  set.seed(ctx$cfg$seed)
  scheme_b_scores <- replicate(n_scheme_b, {
    Xshuf <- shuffle_cci_within_gene(ctx$Xtilde)
    repro <- compute_repro(Xshuf, ctx$cohort_pairs_idx)
    mean(repro$ReproScore, na.rm = TRUE)
  })
  scheme_b_mean <- mean(scheme_b_scores)

  scheme_df <- data.frame(
    scheme = c("A_gene_shuffle", "B_cci_within_gene"),
    null_mean = c(scheme_a_mean, scheme_b_mean),
    n_perm = c(ncol(null_repro), n_scheme_b)
  )
  readr::write_tsv(scheme_df, file.path(validation_dir, "results", "null_scheme_compare.tsv"))

  delta <- abs(scheme_a_mean - scheme_b_mean)
  if (delta > 0.02) {
    append_verdict(validation_dir, "na_exchangeable", 5, "WARN", round(delta, 4),
                   "|delta| <= 0.02", "NA footprint may leak into null")
  } else {
    append_verdict(validation_dir, "na_exchangeable", 5, "PASS", round(delta, 4),
                   "|delta| <= 0.02", "Null robust to NA structure")
  }

  repro_nulls <- ctx$repro_with_nulls
  if (is.null(repro_nulls) && !is.null(ctx$stage05_env)) {
    repro_nulls <- ctx$stage05_env$repro_df
  }
  if (is.null(repro_nulls) || !"shuffle_p" %in% names(repro_nulls)) {
    obs <- ctx$repro_df$ReproScore
    names(obs) <- ctx$repro_df$gene
    p_emp <- vapply(seq_len(nrow(null_repro)), function(g) {
      mean(null_repro[g, ] >= obs[g], na.rm = TRUE)
    }, numeric(1))
    repro_nulls <- ctx$repro_df
    repro_nulls$shuffle_p <- p_emp[match(repro_nulls$gene, ctx$gene_universe)]
    repro_nulls$shuffle_FDR <- p.adjust(repro_nulls$shuffle_p, method = "BH")
  }
  if (!is.null(repro_nulls) && "shuffle_p" %in% names(repro_nulls)) {
    png(file.path(validation_dir, "results", "pval_hist.png"),
        width = 600, height = 500)
    hist(repro_nulls$shuffle_p, breaks = 40, main = "Gene-shuffle p-value distribution",
         xlab = "shuffle_p", col = "steelblue", border = "white")
    dev.off()

    frac_low <- mean(repro_nulls$shuffle_p < 0.05, na.rm = TRUE)
    frac_high <- mean(repro_nulls$shuffle_p > 0.95, na.rm = TRUE)
    if (frac_high > 0.9) {
      append_verdict(validation_dir, "pval_dist", 5, "INFO", round(frac_high, 3),
                     "uniform + spike at 0", "Consistent with no signal")
    } else {
      append_verdict(validation_dir, "pval_dist", 5, "PASS", round(frac_low, 3),
                     "uniform + spike at 0", "Some real signal on null background")
    }
  }

  global_null <- NULL
  if (!is.null(ctx$stage05_env) && !is.null(ctx$stage05_env$global_null)) {
    global_null <- ctx$stage05_env$global_null
  } else if (file.exists(ctx$global_null_txt)) {
    lines <- readLines(ctx$global_null_txt)
    obs_line <- grep("obs_statistic", lines, value = TRUE)
    p_line <- grep("empirical_p", lines, value = TRUE)
    if (length(obs_line) > 0) {
      global_null <- list(
        obs_statistic = as.numeric(sub(".*: ", "", obs_line[1])),
        empirical_p = as.numeric(sub(".*: ", "", p_line[1]))
      )
    }
  }
  if (!is.null(ctx$null_checkpoint) && !is.null(ctx$null_checkpoint$null_global)) {
    null_dist <- ctx$null_checkpoint$null_global
    if (is.null(global_null)) {
      global_null <- list(null_distribution = null_dist)
    } else {
      global_null$null_distribution <- null_dist
    }
  }
  if (!is.null(ctx$stage05_env) && !is.null(ctx$stage05_env$global_null$null_distribution)) {
    global_null <- ctx$stage05_env$global_null
  }

  if (!is.null(global_null)) {
    obs <- global_null$obs_statistic
    null_dist <- global_null$null_distribution
    null_dist <- null_dist[is.finite(null_dist)]
    emp_p <- global_null$empirical_p
    if (is.null(emp_p) && length(null_dist) > 0) {
      emp_p <- mean(null_dist >= obs)
    }

    png(file.path(validation_dir, "results", "global_null_figure.png"),
        width = 700, height = 500)
    hist(null_dist, breaks = 40, main = "Global permutation null",
         xlab = "Null statistic (mean ReproScore - 0.5)", col = "grey70", border = "white")
    abline(v = obs, col = "red", lwd = 2)
    legend("topright", legend = c(sprintf("obs=%.4f", obs), sprintf("p=%.4f", emp_p)),
           col = c("red", NA), lty = c(1, NA), bty = "n")
    dev.off()

    if (is.finite(emp_p) && emp_p < 0.05) {
      append_verdict(validation_dir, "global_test", 5, "PASS", round(emp_p, 4),
                     "p < 0.05", "Atlas-level reproducibility exists")
    } else {
      append_verdict(validation_dir, "global_test", 5, "FAIL", round(emp_p, 4),
                     "p < 0.05", "No atlas-level signal; return to Stage 1/4")
    }
  } else {
    append_verdict(validation_dir, "global_test", 5, "INFO", NA,
                   "p < 0.05", "global_null not available")
  }

  ctx
}
