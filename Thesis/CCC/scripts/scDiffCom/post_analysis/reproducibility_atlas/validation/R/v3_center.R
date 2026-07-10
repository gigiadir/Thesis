# Stage 3 validation — centering.

run_v3_center <- function(ctx, validation_dir) {
  X <- ctx$X
  Xtilde <- ctx$Xtilde
  cohorts <- ctx$cohorts

  stacked_before <- stack_cohort_matrix(X, cohorts)
  stacked_after <- stack_cohort_matrix(Xtilde, cohorts)

  sil_before <- pca_silhouette(stacked_before$mat, stacked_before$cohort_label)
  sil_after <- pca_silhouette(stacked_after$mat, stacked_after$cohort_label)

  png(file.path(validation_dir, "results", "batch_before_after.png"),
      width = 1000, height = 450)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (label in c("Before (raw X)", "After (Xtilde)")) {
    stacked <- if (label == "Before (raw X)") stacked_before else stacked_after
    m <- stacked$mat
    finite_rows <- rowSums(is.finite(m)) > ncol(m) * 0.1
    m <- m[finite_rows, , drop = FALSE]
    cl <- stacked$cohort_label[finite_rows]
    m[is.na(m)] <- 0
    pca <- prcomp(m, center = TRUE, scale. = TRUE)
    cols <- rainbow(length(cohorts))[as.integer(factor(cl))]
    plot(pca$x[, 1], pca$x[, 2], col = cols, pch = 16, cex = 0.5,
         main = label, xlab = "PC1", ylab = "PC2")
    legend("topright", legend = cohorts, col = rainbow(length(cohorts)), pch = 16, cex = 0.7)
  }
  dev.off()

  if (is.finite(sil_after) && sil_after > 0.2) {
    status <- "FAIL"
    note <- "Centering insufficient; batch inherited downstream"
  } else if (is.finite(sil_before) && is.finite(sil_after) &&
             sil_before > sil_after) {
    status <- "PASS"
    note <- "Silhouette dropped after centering"
  } else {
    status <- "WARN"
    note <- "Batch separation change ambiguous"
  }
  append_verdict(validation_dir, "batch_collapse", 3, status,
                 sprintf("before=%.3f,after=%.3f", sil_before, sil_after),
                 "after silhouette < 0.2", note)

  cohort_means <- lapply(cohorts, function(ds) {
    colMeans(Xtilde[[ds]], na.rm = TRUE)
  })
  names(cohort_means) <- cohorts
  n_c <- length(cohorts)
  cor_mat <- matrix(NA_real_, n_c, n_c, dimnames = list(cohorts, cohorts))
  for (i in seq_len(n_c)) {
    for (j in seq_len(n_c)) {
      a <- cohort_means[[i]]
      b <- cohort_means[[j]]
      finite <- is.finite(a) & is.finite(b)
      cor_mat[i, j] <- if (sum(finite) >= 3) cor(a[finite], b[finite], method = "spearman") else NA
    }
  }
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("cohort_a", "cohort_b", "spearman_rho")
  readr::write_tsv(cor_df, file.path(validation_dir, "results", "cohort_mean_cor.tsv"))

  mean_cor <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  append_verdict(validation_dir, "cohort_means", 3, "INFO", round(mean_cor, 3),
                 "descriptive",
                 if (mean_cor > 0.5) "High cohort-mean correlation; batch mostly additive"
                 else "Low cohort-mean correlation; centering load-bearing")

  max_mean_dev <- 0
  max_sd_dev <- 0
  min_genes <- ctx$cfg$min_genes_per_cci
  for (ds in cohorts) {
    mat <- Xtilde[[ds]]
    for (j in seq_len(ncol(mat))) {
      vals <- mat[, j]
      finite <- is.finite(vals)
      if (sum(finite) < min_genes) next
      mu <- mean(vals[finite])
      sigma <- stats::sd(vals[finite])
      if (!is.finite(sigma) || sigma == 0) next
      max_mean_dev <- max(max_mean_dev, abs(mu))
      if (sigma < 0.9 || sigma > 1.1) {
        max_sd_dev <- max(max_sd_dev, max(abs(sigma - 0.9), abs(sigma - 1.1)))
      }
    }
  }

  if (max_mean_dev > 1e-6 || max_sd_dev > 0) {
    append_verdict(validation_dir, "center_moments", 3, "FAIL",
                   sprintf("max|mean|=%.2e,max_sd_dev=%.3f", max_mean_dev, max_sd_dev),
                   "mean~0, sd in [0.9,1.1]", "NA-handling bug in per-column moments")
  } else {
    append_verdict(validation_dir, "center_moments", 3, "PASS",
                   sprintf("max|mean|=%.2e", max_mean_dev),
                   "mean~0, sd in [0.9,1.1]", "Post-centering moments OK")
  }

  ctx
}
