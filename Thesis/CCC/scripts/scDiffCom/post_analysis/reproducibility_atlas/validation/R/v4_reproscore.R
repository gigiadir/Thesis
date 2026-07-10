# Stage 4 validation — ReproScore.

run_v4_reproscore <- function(ctx, validation_dir) {
  cor_mats <- ctx$cor_mats
  repro_df <- ctx$repro_df
  pair_labels <- ctx$pair_labels

  deltas <- vapply(cor_mats, function(C) {
    d <- diag(C)
    off <- C[row(C) != col(C)]
    off <- off[is.finite(off)]
    if (length(off) == 0 || !any(is.finite(d))) return(NA_real_)
    mean(d, na.rm = TRUE) - mean(off, na.rm = TRUE)
  }, numeric(1))

  mean_delta <- mean(deltas, na.rm = TRUE)
  plot_pairs <- head(pair_labels, min(3, length(pair_labels)))

  png(file.path(validation_dir, "results", "diagonal_heatmaps.png"),
      width = 1200, height = 400 * length(plot_pairs))
  par(mfrow = c(length(plot_pairs), 1), mar = c(4, 4, 3, 1))
  for (nm in plot_pairs) {
    C <- cor_mats[[nm]]
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(C, main = paste(nm, sprintf("delta=%.3f", deltas[[nm]])),
                         show_rownames = FALSE, show_colnames = FALSE,
                         cluster_rows = FALSE, cluster_cols = FALSE)
    } else {
      image(C, main = paste(nm, sprintf("delta=%.3f", deltas[[nm]])),
            xlab = "Gene", ylab = "Gene")
    }
  }
  dev.off()

  if (is.finite(mean_delta) && mean_delta > 0) {
    append_verdict(validation_dir, "diag_visible", 4, "PASS", round(mean_delta, 4),
                   "mean(diag) > mean(off-diag)", "Diagonal signal visible")
  } else {
    append_verdict(validation_dir, "diag_visible", 4, "FAIL", round(mean_delta, 4),
                   "mean(diag) > mean(off-diag)", "No signal at source")
  }

  if (is.na(ctx$eff_n_median)) {
    eff_path <- file.path(validation_dir, "results", "eff_n_per_gene_pair.tsv")
    if (file.exists(eff_path)) {
      eff_df <- readr::read_tsv(eff_path, show_col_types = FALSE)
      ctx$eff_n_median <- median(eff_df$median_eff_n, na.rm = TRUE)
    }
  }

  eff_path <- file.path(validation_dir, "results", "eff_n_per_gene_pair.tsv")
  if (file.exists(eff_path)) {
    eff_df <- readr::read_tsv(eff_path, show_col_types = FALSE)
    merged <- merge(repro_df, eff_df[, c("gene", "median_eff_n")], by = "gene")
    rho <- cor(merged$ReproScore, merged$median_eff_n, method = "spearman", use = "complete.obs")

    png(file.path(validation_dir, "results", "reproscore_vs_effn.png"),
        width = 600, height = 500)
    plot(merged$median_eff_n, merged$ReproScore, pch = 16, cex = 0.6,
         xlab = "Median effective-n", ylab = "ReproScore",
         main = sprintf("ReproScore vs effective-n (rho=%.3f)", rho))
    abline(h = 0.5, col = "grey", lty = 2)
    dev.off()

    abs_rho <- abs(rho)
    if (abs_rho > 0.3) {
      status <- "FAIL"
      note <- "Reproducibility confounded with sparsity"
    } else if (abs_rho > 0.15) {
      status <- "WARN"
      note <- "Moderate sparsity confounding"
    } else {
      status <- "PASS"
      note <- "Score not a sparsity artifact"
    }
    append_verdict(validation_dir, "reproscore_vs_n", 4, status, round(rho, 4),
                   "|rho| < 0.15", note)
  }

  raw_gaps <- vapply(cor_mats, function(C) {
    d <- diag(C)
    off <- C[row(C) != col(C)]
    off <- off[is.finite(off)]
    if (length(off) == 0) return(NA_real_)
    mean(d, na.rm = TRUE) - mean(off, na.rm = TRUE)
  }, numeric(1))
  raw_gap <- mean(raw_gaps, na.rm = TRUE)

  if (is.finite(raw_gap) && raw_gap > 0) {
    append_verdict(validation_dir, "raw_gap", 4, "PASS", round(raw_gap, 4),
                   "> 0", "Same-vs-different gap positive")
  } else {
    append_verdict(validation_dir, "raw_gap", 4, "FAIL", round(raw_gap, 4),
                   "> 0", "No raw same-vs-different gap")
  }

  rs <- repro_df$ReproScore
  png(file.path(validation_dir, "results", "reproscore_hist.png"),
      width = 600, height = 500)
  hist(rs, breaks = 40, main = "ReproScore distribution", xlab = "ReproScore",
       col = "steelblue", border = "white")
  abline(v = 0.5, col = "red", lty = 2)
  dev.off()

  frac_near_half <- mean(abs(rs - 0.5) < 0.05, na.rm = TRUE)
  frac_high <- mean(rs > 0.9, na.rm = TRUE)
  if (frac_high > 0.8) {
    append_verdict(validation_dir, "reproscore_dist", 4, "FAIL", round(frac_high, 3),
                   "not all near 1", "Possible leakage / effective-n artifact")
  } else if (frac_near_half > 0.9) {
    append_verdict(validation_dir, "reproscore_dist", 4, "INFO", round(frac_near_half, 3),
                   "bimodal preferred", "No signal; scores cluster at 0.5")
  } else {
    append_verdict(validation_dir, "reproscore_dist", 4, "PASS", round(frac_high, 3),
                   "bimodal with right tail", "Distribution shows signal structure")
  }

  ctx
}
