# Stage 2 validation — tensor integrity.

run_v2_tensor <- function(ctx, validation_dir) {
  X <- ctx$X
  cohorts <- ctx$cohorts
  gene_universe <- ctx$gene_universe
  malignant_by_cohort <- ctx$malignant_by_cohort
  J <- ctx$J

  na_frac_mat <- sapply(cohorts, function(ds) {
    rowMeans(is.na(X[[ds]]))
  })
  rownames(na_frac_mat) <- gene_universe

  na_gene_cohort <- sapply(cohorts, function(ds) rowMeans(is.na(X[[ds]])))
  low_cov <- gene_universe[apply(na_gene_cohort, 1, max) > 0.85]
  if (length(low_cov) > 0) {
    readr::write_tsv(
      data.frame(gene = low_cov, stringsAsFactors = FALSE),
      file.path(validation_dir, "results", "low_coverage_genes.tsv")
    )
  }

  cohort_mean_na <- colMeans(na_gene_cohort)
  worst_cohort <- cohorts[which.max(cohort_mean_na)]
  worst_na <- max(cohort_mean_na)

  if (worst_na > 0.85) {
    append_verdict(validation_dir, "na_map", 2, "WARN", round(worst_na, 3),
                   "mean NA < 0.85 per cohort", paste("Worst:", worst_cohort))
  } else {
    append_verdict(validation_dir, "na_map", 2, "INFO", round(worst_na, 3),
                   "mean NA < 0.85 per cohort", "NA fractions within expected range")
  }

  png(file.path(validation_dir, "results", "na_fraction_heatmap.png"),
      width = 800, height = 900)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(na_frac_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = FALSE, main = "Gene x cohort NA fraction")
  } else {
    image(t(na_frac_mat[seq(1, nrow(na_frac_mat), length.out = min(200, nrow(na_frac_mat))), ]),
          main = "Gene x cohort NA fraction (subsampled)", xlab = "Cohort", ylab = "Gene (subsampled)")
  }
  dev.off()

  iqrs <- vapply(cohorts, function(ds) {
    vals <- as.numeric(X[[ds]])
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) return(NA_real_)
    stats::IQR(vals)
  }, numeric(1))
  iqr_ratio <- max(iqrs, na.rm = TRUE) / min(iqrs, na.rm = TRUE)

  png(file.path(validation_dir, "results", "logfc_scale_by_cohort.png"),
      width = 700, height = 500)
  boxplot(
    lapply(cohorts, function(ds) as.numeric(X[[ds]])[is.finite(X[[ds]])]),
    names = cohorts, main = "Raw LOGFC spread by cohort", ylab = "LOGFC", las = 2
  )
  dev.off()

  if (iqr_ratio > 2.5) {
    append_verdict(validation_dir, "logfc_scale", 2, "INFO", round(iqr_ratio, 3),
                   "descriptive", "Centering sigma-division is load-bearing")
  } else {
    append_verdict(validation_dir, "logfc_scale", 2, "INFO", round(iqr_ratio, 3),
                   "descriptive", "Cohort LOGFC scales similar")
  }

  total_cells <- 0L
  collapsed_cells <- 0L
  for (ds in cohorts) {
    gene_list <- malignant_by_cohort[[ds]]
    for (g in names(gene_list)) {
      df <- gene_list[[g]]
      if (is.null(df) || nrow(df) == 0) next
      df <- df[df$CCI %in% J, , drop = FALSE]
      if (nrow(df) == 0) next
      tab <- table(df$CCI)
      total_cells <- total_cells + length(tab)
      collapsed_cells <- collapsed_cells + sum(tab > 1)
    }
  }
  collapse_frac <- if (total_cells > 0) collapsed_cells / total_cells else 0

  if (collapse_frac > 0.1) {
    append_verdict(validation_dir, "dup_collapse", 2, "WARN", round(collapse_frac, 4),
                   "<0.10", "Many (gene,CCI) cells built from >1 averaged row")
  } else {
    append_verdict(validation_dir, "dup_collapse", 2, "INFO", round(collapse_frac, 4),
                   "<0.10", "Duplicate collapse fraction acceptable")
  }

  ctx
}
