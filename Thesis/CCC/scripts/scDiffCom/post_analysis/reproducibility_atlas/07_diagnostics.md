# Diagnostics




``` r
output_dir <- atlas_env$output_dir
diag_dir <- file.path(output_dir, "results", "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

message("Diagnostics")
```

```
## Diagnostics
```

``` r
vocab_report <- readr::read_tsv(
  file.path(output_dir, "results", "vocab_report.tsv"), show_col_types = FALSE
)
writeLines(
  c("=== Vocabulary support-loss recap ===", capture.output(print(vocab_report))),
  file.path(diag_dir, "support_loss_recap.txt")
)

cor_mats <- atlas_env$repro$cor_mats
mean_cor <- Reduce(`+`, cor_mats) / length(cor_mats)
diag(mean_cor) <- NA

if (requireNamespace("pheatmap", quietly = TRUE)) {
  png(file.path(diag_dir, "mean_gene_gene_cor_heatmap.png"), width = 900, height = 800)
  pheatmap::pheatmap(
    mean_cor, main = "Mean cross-cohort gene x gene Spearman",
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = FALSE, show_colnames = FALSE
  )
  dev.off()
}
```

```
## png 
##   3
```

``` r
hc <- hclust(as.dist(1 - mean_cor), method = "average")
clusters <- cutree(hc, k = min(8, max(2, floor(nrow(mean_cor) / 20))))
rd <- atlas_env$repro_df
module_tbl <- data.frame(
  gene = rownames(mean_cor), module = clusters,
  Rg = rd$Rg[match(rownames(mean_cor), rd$gene)],
  ReproScore = rd$ReproScore[match(rownames(mean_cor), rd$gene)],
  stringsAsFactors = FALSE
)
readr::write_tsv(module_tbl, file.path(diag_dir, "gene_modules.tsv"))

module_summary <- module_tbl %>%
  group_by(module) %>%
  summarise(
    n_genes = n(),
    mean_Rg = mean(Rg, na.rm = TRUE),
    mean_ReproScore = mean(ReproScore, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_tsv(module_summary, file.path(diag_dir, "module_summary.tsv"))

# Pair-coverage reporting: how many cohort pairs contribute to each gene's R_g.
coverage_tbl <- rd %>%
  dplyr::count(n_pairs_computable, name = "n_genes") %>%
  arrange(n_pairs_computable)
readr::write_tsv(coverage_tbl, file.path(diag_dir, "pair_coverage.tsv"))

# EVT calibration method breakdown (empirical vs GPD vs fallback).
if (!is.null(rd$evt_method)) {
  evt_method_tbl <- rd %>% dplyr::count(evt_method, name = "n_genes")
  readr::write_tsv(evt_method_tbl, file.path(diag_dir, "evt_method_breakdown.tsv"))
}

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

module_min_overlap <- if (!is.null(cfg$min_cci_overlap)) cfg$min_cci_overlap else 10L
module_agg <- if (!is.null(cfg$rg_aggregate)) cfg$rg_aggregate else "median"
module_cor_mats <- compute_pair_cor_matrices(centroid_mats, idx_pairs, min_overlap = module_min_overlap)
module_rg <- compute_Rg_from_cormats(module_cor_mats, aggregate = module_agg)
readr::write_tsv(
  data.frame(
    module = paste0("module_", module_ids),
    Rg = module_rg$Rg,
    Rg_mean = module_rg$Rg_mean,
    n_pairs_computable = module_rg$n_pairs_computable,
    stringsAsFactors = FALSE
  ),
  file.path(diag_dir, "module_reproscore.tsv")
)

fdr_col <- if (!is.null(atlas_env$repro_df$evt_FDR)) "evt_FDR" else "shuffle_FDR"
disagree <- atlas_env$repro_df %>%
  mutate(category = case_when(
    .data[[fdr_col]] < cfg$fdr_threshold ~ "atlas_member_fdr_sig",
    .data[[fdr_col]] >= cfg$fdr_threshold & Rg > median(Rg, na.rm = TRUE) ~
      "high_Rg_not_fdr_sig",
    .data[[fdr_col]] >= cfg$fdr_threshold & Rg <= median(Rg, na.rm = TRUE) ~
      "low_Rg_not_fdr_sig",
    TRUE ~ "other"
  )) %>%
  left_join(module_tbl %>% select(gene, module), by = "gene")
readr::write_tsv(disagree, file.path(diag_dir, "reproscore_null_disagreement.tsv"))

atlas_env$diagnostics <- list(module_tbl = module_tbl, module_summary = module_summary, disagree = disagree)
saveRDS(atlas_env, file.path(output_dir, "results", "stage_diag_atlas_env.rds"))
module_summary
```

```
## # A tibble: 8 × 4
##   module n_genes mean_Rg mean_ReproScore
##    <int>   <int>   <dbl>           <dbl>
## 1      1     148  0.227            0.671
## 2      2     133  0.224            0.696
## 3      3      54  0.152            0.633
## 4      4      35  0.132            0.612
## 5      5       4 -0.150            0.398
## 6      6      14  0.0935           0.582
## 7      7       2 -0.0380           0.428
## 8      8       2 -0.111            0.398
```
