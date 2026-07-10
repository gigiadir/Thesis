# Stage 1 validation — vocabulary and effective-n.

run_v1_vocab <- function(ctx, validation_dir) {
  X <- ctx$X
  cohorts <- ctx$cohorts
  gene_universe <- ctx$gene_universe
  J <- ctx$J

  eff_mat <- compute_eff_n(X, cohorts, ctx$cohort_pairs_idx, gene_universe)
  eff_df <- as.data.frame(eff_mat)
  eff_df$gene <- rownames(eff_mat)
  eff_df$median_eff_n <- apply(eff_mat, 1, median, na.rm = TRUE)
  readr::write_tsv(
    eff_df[, c("gene", colnames(eff_mat), "median_eff_n")],
    file.path(validation_dir, "results", "eff_n_per_gene_pair.tsv")
  )

  median_eff_n <- median(eff_df$median_eff_n, na.rm = TRUE)
  ctx$eff_n_median <- median_eff_n

  if (median_eff_n < 15) {
    status <- "FAIL"
    note <- "Spearman unstable; ReproScore largely luck"
  } else if (median_eff_n <= 30) {
    status <- "WARN"
    note <- "Marginal effective sample size"
  } else {
    status <- "PASS"
    note <- "Adequate pairwise-complete CCI counts"
  }
  append_verdict(validation_dir, "eff_n", 1, status, median_eff_n,
                 ">30 (PASS), 15-30 (WARN), <15 (FAIL)", note)

  png(file.path(validation_dir, "results", "eff_n_median_hist.png"),
      width = 700, height = 500)
  hist(eff_df$median_eff_n, breaks = 30, main = "Per-gene median effective-n",
       xlab = "Median pairwise-complete CCI count", col = "steelblue", border = "white")
  abline(v = median_eff_n, col = "red", lwd = 2, lty = 2)
  dev.off()

  cci_meta <- parse_cci(J)
  ct_pair <- paste(cci_meta$EMITTER_CELLTYPE, cci_meta$RECEIVER_CELLTYPE, sep = "->")
  ct_tab <- sort(table(ct_pair), decreasing = TRUE)
  lri_tab <- sort(table(cci_meta$LRI), decreasing = TRUE)

  comp_df <- data.frame(
    type = c(rep("celltype_pair", length(ct_tab)), rep("LRI", min(20, length(lri_tab)))),
    label = c(names(ct_tab), names(lri_tab)[seq_len(min(20, length(lri_tab)))]),
    count = c(as.integer(ct_tab), as.integer(lri_tab)[seq_len(min(20, length(lri_tab)))]),
    fraction = c(as.numeric(ct_tab) / length(J),
               as.integer(lri_tab)[seq_len(min(20, length(lri_tab)))] / length(J)),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(comp_df, file.path(validation_dir, "results", "vocab_composition.tsv"))

  max_frac <- max(as.numeric(ct_tab)) / length(J)
  if (max_frac > 0.6) {
    append_verdict(validation_dir, "vocab_comp", 1, "WARN", round(max_frac, 3),
                   "<0.60", "Single celltype-pair dominates J")
  } else {
    append_verdict(validation_dir, "vocab_comp", 1, "INFO", round(max_frac, 3),
                   "characterization", "Vocabulary spread across celltype pairs")
  }

  ct_df <- data.frame(celltype_pair = names(ct_tab), count = as.integer(ct_tab),
                      fraction = as.numeric(ct_tab) / length(J), stringsAsFactors = FALSE)
  png(file.path(validation_dir, "results", "vocab_composition_bar.png"),
      width = 900, height = 500)
  par(mar = c(10, 4, 3, 1))
  barplot(head(ct_df$count, 15), names.arg = head(ct_df$celltype_pair, 15),
          las = 2, main = "Top celltype-pair CCI counts", ylab = "Count", col = "steelblue")
  dev.off()

  ctx
}
