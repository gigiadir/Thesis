# Stage 3 — cohort-centering: remove cohort main effect per CCI.

.center_cohort_matrix <- function(mat, min_genes = 10, eps = 1e-8) {
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

run_stage_03_center <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  cohorts <- atlas_env$cohorts
  gene_universe <- atlas_env$gene_universe
  X <- atlas_env$X
  min_genes <- cfg$min_genes_per_cci

  message("Stage 3: cohort-centering")

  Xtilde <- lapply(cohorts, function(ds) {
    res <- .center_cohort_matrix(X[[ds]], min_genes = min_genes)
    attr(res$mat, "zero_sd_dropped") <- res$zero_sd
    res$mat
  })
  names(Xtilde) <- cohorts

  all_zero_sd <- unique(unlist(lapply(Xtilde, function(m) attr(m, "zero_sd_dropped"))))
  if (length(all_zero_sd) > 0) {
    message(sprintf("  CCIs with zero sd (kept as NA-centered): %d", length(all_zero_sd)))
  }

  set.seed(cfg$seed)
  hk_genes <- sample(gene_universe, min(20, length(gene_universe)))

  baseline_tbl <- purrr::imap_dfr(cohorts, function(ds, idx) {
    mat <- X[[ds]]
    hk <- mat[hk_genes, , drop = FALSE]
    all_g <- mat
    m_hk <- colMeans(hk, na.rm = TRUE)
    m_all <- colMeans(all_g, na.rm = TRUE)
    finite <- is.finite(m_hk) & is.finite(m_all)
    rho <- if (sum(finite) >= 3) {
      stats::cor(m_hk[finite], m_all[finite], method = "spearman")
    } else {
      NA_real_
    }
    data.frame(
      cohort = ds,
      n_hk_genes = nrow(hk),
      spearman_rho = rho,
      flag_low_rho = is.finite(rho) && rho < 0.7,
      stringsAsFactors = FALSE
    )
  })

  readr::write_tsv(baseline_tbl, file.path(output_dir, "results", "centering_baseline_check.tsv"))
  saveRDS(Xtilde, file.path(output_dir, "results", "Xtilde.rds"))

  atlas_env$Xtilde <- Xtilde
  atlas_env$centering_baseline <- baseline_tbl
  saveRDS(atlas_env, file.path(output_dir, "results", "stage03_atlas_env.rds"))
  message("  Saved Xtilde.rds")

  invisible(atlas_env)
}
