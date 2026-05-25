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

plot_lri_heatmap_umap <- function(malignant_list, dataset_label,
                                  k_clusters = 5,
                                  min_genes  = MIN_GENES_PER_CCI) {
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

  # Heatmap
  heatmap_file <- file.path(OUTPUT_DIR, paste0(dataset_label, "_lri_heatmap.png"))
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
# 
#   # UMAP
#   mat_genes         <- t(mat)
#   mat_genes_imputed <- mat_genes
#   mat_genes_imputed[is.na(mat_genes_imputed)] <- 0
# 
#   dist_mat   <- as.dist(1 - sim_mat)
#   hclust_res <- hclust(dist_mat, method = "complete")
#   clusters   <- cutree(hclust_res, k = k_clusters)
# 
#   set.seed(42)
#   umap_config             <- umap.defaults
#   umap_config$n_neighbors <- min(15, nrow(mat_genes_imputed) - 1)
#   umap_config$min_dist    <- 0.1
#   umap_config$metric      <- "cosine"
# 
#   umap_result <- umap(mat_genes_imputed, config = umap_config)
# 
#   umap_df <- data.frame(
#     Gene    = rownames(mat_genes_imputed),
#     UMAP1   = umap_result$layout[, 1],
#     UMAP2   = umap_result$layout[, 2],
#     n_LRIs  = rowSums(!is.na(mat_genes)),
#     Cluster = factor(clusters[rownames(mat_genes_imputed)])
#   )
# 
#   p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = Cluster, label = Gene)) +
#     geom_point(size = 3, alpha = 0.9) +
#     geom_text_repel(size = 2.8, max.overlaps = 40,
#                     segment.size = 0.2, segment.alpha = 0.5, show.legend = FALSE) +
#     scale_colour_brewer(palette = "Set1", name = "Cluster") +
#     labs(
#       title    = paste0(dataset_label, " – Gene UMAP (LRI LOGFC space)"),
#       subtitle = paste0("Hierarchical clustering (k = ", k_clusters, ", cosine distance)"),
#       x = "UMAP 1", y = "UMAP 2"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(plot.title = element_text(face = "bold"))
# 
#   umap_file <- file.path(OUTPUT_DIR, paste0(dataset_label, "_lri_gene_umap.png"))
#   ggsave(umap_file, plot = p_umap, width = 10, height = 8, dpi = 300)
#   message("  UMAP saved to ", umap_file)
# 
#   cluster_df <- umap_df %>%
#     select(Gene, Cluster, n_LRIs) %>%
#     arrange(Cluster, Gene)
#   print(cluster_df)
# 
#   cluster_summary <- cluster_df %>%
#     group_by(Cluster) %>%
#     summarise(n_genes = n(), genes = paste(sort(Gene), collapse = ", ")) %>%
#     arrange(Cluster)
#   print(cluster_summary)
# 
#   invisible(list(umap_df = umap_df, clusters = clusters))
}
