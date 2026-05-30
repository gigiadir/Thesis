library(umap)

build_lri_logfc_mat <- function(malignant_list, min_genes = MIN_GENES_PER_CCI) {
  lri_list <- lapply(names(malignant_list), function(gene) {
    malignant_list[[gene]] %>%
      select(LRI, LOGFC) %>%
      group_by(LRI) %>%
      summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop") %>%
      mutate(Gene = gene)
  })

  all_data <- bind_rows(lri_list)
  mat_long <- all_data %>%
    pivot_wider(names_from = Gene, values_from = LOGFC)

  mat <- as.matrix(mat_long[, -1])
  rownames(mat) <- mat_long$LRI
  message("  Matrix dimensions (LRIs x Genes): ", nrow(mat), " x ", ncol(mat))

  n_non_na <- rowSums(!is.na(mat))
  mat <- mat[n_non_na >= min_genes, , drop = FALSE]
  message("  After filtering sparse LRIs: ", nrow(mat), " x ", ncol(mat))
  mat
}

.plot_lri_umap_helper_path <- function() {
  pa_dir <- if (exists("POST_ANALYSIS_DIR", inherits = TRUE)) {
    get("POST_ANALYSIS_DIR", inherits = TRUE)
  } else {
    normalizePath(file.path(getwd(), ".."), winslash = "/")
  }
  file.path(pa_dir, "within_cohort", "R", "04_within_cohort_helpers.R")
}

.ensure_within_cohort_helpers <- function() {
  helper_path <- .plot_lri_umap_helper_path()
  if (!exists("plot_gene_umap_lri", mode = "function", inherits = TRUE) &&
      file.exists(helper_path)) {
    source(helper_path, local = FALSE)
  }
}

plot_lri_heatmap_umap <- function(malignant_list, dataset_label,
                                  k_clusters = 5,
                                  min_genes  = MIN_GENES_PER_CCI,
                                  run_umap   = TRUE,
                                  output_dir = OUTPUT_DIR) {
  message("\n=== GENES x LRI: ", dataset_label, " ===")

  mat     <- build_lri_logfc_mat(malignant_list, min_genes)
  sim_mat <- cosine_sim_na(mat)
  diag(sim_mat) <- NA

  gene_order <- sort(rownames(sim_mat))
  sim_mat    <- sim_mat[gene_order, gene_order]
  n_genes    <- nrow(sim_mat)

  sim_mat_clust <- sim_mat
  diag(sim_mat_clust) <- 1
  sim_mat_clust[is.na(sim_mat_clust)] <- 0
  dist_for_clust <- as.dist(1 - sim_mat_clust)

  heatmap_file <- file.path(output_dir, paste0(dataset_label, "_lri_heatmap.png"))
  png(heatmap_file,
      width  = max(4000, n_genes * 20),
      height = max(4000, n_genes * 20),
      res    = 300)
  pheatmap(
    sim_mat,
    color        = colorRampPalette(c("red", "white", "blue"))(100),
    breaks       = seq(-1, 1, length.out = 101),
    na_col       = "grey80",
    fontsize_row = max(4, 200 / n_genes),
    fontsize_col = max(4, 200 / n_genes),
    angle_col    = "45",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = dist_for_clust,
    clustering_distance_cols = dist_for_clust,
    main         = paste0(dataset_label, " – Gene-Gene Cosine Similarity (LRI)")
  )
  dev.off()
  message("  Heatmap saved to ", heatmap_file)

  if (run_umap) {
    .ensure_within_cohort_helpers()
    if (exists("plot_gene_umap_lri", mode = "function", inherits = TRUE)) {
      plot_gene_umap_lri(
        mat, sim_mat, dataset_label, output_dir,
        k_clusters = k_clusters
      )
    } else {
      message("  UMAP skipped — within_cohort helpers not found.")
    }
  }

  invisible(list(mat = mat, sim_mat = sim_mat))
}
