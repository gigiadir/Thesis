# Diagnostics — support loss recap, off-diagonal blocks, ReproScore vs null disagreement.

run_diagnostics <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  diag_dir <- file.path(output_dir, "results", "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  message("Diagnostics")

  # 1. Support-loss recap
  vocab_report <- readr::read_tsv(
    file.path(output_dir, "results", "vocab_report.tsv"),
    show_col_types = FALSE
  )
  writeLines(
    c(
      "=== Vocabulary support-loss recap ===",
      capture.output(print(vocab_report))
    ),
    file.path(diag_dir, "support_loss_recap.txt")
  )

  # 2. Off-diagonal block detection via mean cross-cohort gene correlation
  cor_mats <- atlas_env$repro$cor_mats
  mean_cor <- Reduce(`+`, cor_mats) / length(cor_mats)
  diag(mean_cor) <- NA

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    png(file.path(diag_dir, "mean_gene_gene_cor_heatmap.png"), width = 900, height = 800)
    pheatmap::pheatmap(
      mean_cor,
      main = "Mean cross-cohort gene x gene Spearman",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      show_colnames = FALSE
    )
    dev.off()
  }

  hc <- stats::hclust(stats::as.dist(1 - mean_cor), method = "average")
  clusters <- stats::cutree(hc, k = min(8, max(2, floor(nrow(mean_cor) / 20))))
  module_tbl <- data.frame(
    gene = rownames(mean_cor),
    module = clusters,
    ReproScore = atlas_env$repro_df$ReproScore[match(rownames(mean_cor), atlas_env$repro_df$gene)],
    stringsAsFactors = FALSE
  )
  readr::write_tsv(module_tbl, file.path(diag_dir, "gene_modules.tsv"))

  module_summary <- module_tbl %>%
    group_by(module) %>%
    summarise(
      n_genes = n(),
      mean_ReproScore = mean(ReproScore, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_tsv(module_summary, file.path(diag_dir, "module_summary.tsv"))

  # Optional module-level ReproScore on centroids
  Xtilde <- atlas_env$Xtilde
  cohorts <- atlas_env$cohorts
  idx_pairs <- atlas_env$cohort_pairs
  module_ids <- sort(unique(clusters))

  centroid_mats <- lapply(cohorts, function(ds) {
    mat <- Xtilde[[ds]]
    cm <- t(vapply(module_ids, function(m) {
      genes_m <- module_tbl$gene[module_tbl$module == m]
      colMeans(mat[genes_m, , drop = FALSE], na.rm = TRUE)
    }, numeric(ncol(mat))))
    rownames(cm) <- paste0("module_", module_ids)
    cm
  })
  names(centroid_mats) <- cohorts

  module_repro <- .compute_repro(centroid_mats, idx_pairs)
  module_repro_df <- data.frame(
    module = paste0("module_", module_ids),
    ReproScore = module_repro$ReproScore,
    R_self = module_repro$R_self,
    stringsAsFactors = FALSE
  )
  readr::write_tsv(module_repro_df, file.path(diag_dir, "module_reproscore.tsv"))

  # 3. ReproScore vs shuffle-null disagreement
  disagree <- atlas_env$repro_df %>%
    mutate(
      category = case_when(
        ReproScore > cfg$reproscore_threshold & shuffle_FDR >= cfg$fdr_threshold ~
          "high_ReproScore_high_shuffle_p_distinctive",
        ReproScore <= cfg$reproscore_threshold & shuffle_FDR < cfg$fdr_threshold ~
          "low_ReproScore_sig_null_investigate",
        TRUE ~ "other"
      )
    ) %>%
    left_join(module_tbl %>% select(gene, module), by = "gene")

  readr::write_tsv(disagree, file.path(diag_dir, "reproscore_null_disagreement.tsv"))

  atlas_env$diagnostics <- list(
    module_tbl = module_tbl,
    module_summary = module_summary,
    disagree = disagree
  )
  saveRDS(atlas_env, file.path(output_dir, "results", "stage_diag_atlas_env.rds"))
  message("  Diagnostics written to results/diagnostics/")

  invisible(atlas_env)
}
