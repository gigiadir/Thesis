# Within-cohort differential communication phenotype helpers

is_unknown_celltype <- function(ct) {
  grepl("unknown|other|equivocal", ct, ignore.case = TRUE)
}

filter_unknown_celltypes_df <- function(cci_df) {
  cci_df %>%
    filter(
      !is_unknown_celltype(EMITTER_CELLTYPE),
      !is_unknown_celltype(RECEIVER_CELLTYPE)
    )
}

filter_unknown_celltypes_list <- function(malignant_list) {
  lapply(malignant_list, filter_unknown_celltypes_df)
}

plot_goi_boxplot <- function(malignant_list,
                             gene,
                             dataset_label,
                             dataset_short,
                             output_dir,
                             malignant_label = "H&N Cancer",
                             top_pairs_to_label = 8,
                             plot_width = 18,
                             plot_height = 10) {
  if (!gene %in% names(malignant_list)) {
    stop("Gene ", gene, " not found in malignant list for ", dataset_label)
  }

  plot_df <- malignant_list[[gene]] %>%
    filter_unknown_celltypes_df() %>%
    mutate(ER_PAIR = paste(EMITTER_CELLTYPE, "\u2192", RECEIVER_CELLTYPE))

  pairs_by_spread <- plot_df %>%
    group_by(ER_PAIR) %>%
    summarise(
      logfc_min = min(LOGFC),
      logfc_max = max(LOGFC),
      spread = logfc_max - logfc_min,
      .groups = "drop"
    ) %>%
    slice_max(order_by = spread, n = top_pairs_to_label, with_ties = FALSE) %>%
    pull(ER_PAIR)

  label_df <- bind_rows(
    plot_df %>%
      filter(ER_PAIR %in% pairs_by_spread) %>%
      group_by(ER_PAIR) %>%
      slice_min(order_by = LOGFC, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(extreme = "min"),
    plot_df %>%
      filter(ER_PAIR %in% pairs_by_spread) %>%
      group_by(ER_PAIR) %>%
      slice_max(order_by = LOGFC, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(extreme = "max")
  ) %>%
    distinct(LRI, ER_PAIR, LOGFC, .keep_all = TRUE)

  n_interactions <- nrow(plot_df)
  label_nudge <- 0.35
  y_vals <- c(
    plot_df$LOGFC,
    label_df$LOGFC + label_nudge,
    label_df$LOGFC - label_nudge
  )
  y_lim <- range(y_vals, na.rm = TRUE)
  y_pad <- max(0.08, diff(y_lim) * 0.04)

  p <- ggplot(plot_df, aes(x = ER_PAIR, y = LOGFC, fill = ER_PAIR, color = ER_PAIR)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.4, color = "black", width = 0.55) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.6) +
    geom_text_repel(
      data = dplyr::filter(label_df, extreme == "max"),
      aes(x = ER_PAIR, y = LOGFC, label = LRI),
      size = 2.4,
      direction = "y",
      nudge_y = label_nudge,
      max.overlaps = 30,
      segment.size = 0.2,
      segment.alpha = 0.5,
      box.padding = 0.45,
      point.padding = 0.5,
      min.segment.length = 0.1,
      force = 1.5,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_text_repel(
      data = dplyr::filter(label_df, extreme == "min"),
      aes(x = ER_PAIR, y = LOGFC, label = LRI),
      size = 2.4,
      direction = "y",
      nudge_y = -label_nudge,
      max.overlaps = 30,
      segment.size = 0.2,
      segment.alpha = 0.5,
      box.padding = 0.45,
      point.padding = 0.5,
      min.segment.length = 0.1,
      force = 1.5,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      clip = "off"
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(gene, " In ", malignant_label, ": L-R LogFC across Emitter-Receiver Pairs"),
      subtitle = paste0(
        "Dataset: ", dataset_short, " (", dataset_label, ")",
        "  |  ", n_interactions, " interactions"
      ),
      x = "Cell-Type Interaction (Emitter \u2192 Receiver)",
      y = "LogFC"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      plot.margin = margin(15, 15, 15, 15)
    )

  out_file <- file.path(output_dir, paste0(gene, "_", dataset_short, "_CCI_boxplot.png"))
  ggsave(out_file, plot = p, width = plot_width, height = plot_height, dpi = 300, bg = "white")
  message("  Boxplot saved to ", out_file, " (", nrow(label_df), " min/max LRI labels)")

  label_csv <- file.path(output_dir, paste0(gene, "_", dataset_short, "_CCI_labeled_points.csv"))
  write.csv(
    label_df %>%
      select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, ER_PAIR, LOGFC, REGULATION, extreme) %>%
      arrange(ER_PAIR, LOGFC),
    label_csv,
    row.names = FALSE
  )
  message("  Labeled points saved to ", label_csv)

  invisible(list(plot = p, labels = label_df))
}

plot_gene_umap_lri <- function(mat,
                               sim_mat,
                               dataset_label,
                               output_dir,
                               k_clusters = 5,
                               umap_seed  = 42) {
  mat_genes <- t(mat)
  mat_genes_imputed <- mat_genes
  mat_genes_imputed[is.na(mat_genes_imputed)] <- 0

  sim_mat_clust <- sim_mat
  diag(sim_mat_clust) <- 1
  sim_mat_clust[is.na(sim_mat_clust)] <- 0
  dist_mat <- as.dist(1 - sim_mat_clust)
  hclust_res <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hclust_res, k = k_clusters)

  set.seed(umap_seed)
  umap_config <- umap.defaults
  umap_config$n_neighbors <- min(15, nrow(mat_genes_imputed) - 1)
  umap_config$min_dist <- 0.1
  umap_config$metric <- "cosine"

  umap_result <- umap(mat_genes_imputed, config = umap_config)

  umap_df <- data.frame(
    Gene    = rownames(mat_genes_imputed),
    UMAP1   = umap_result$layout[, 1],
    UMAP2   = umap_result$layout[, 2],
    n_LRIs  = rowSums(!is.na(mat_genes)),
    Cluster = factor(clusters[rownames(mat_genes_imputed)])
  )

  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = Cluster, label = Gene)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      size = 2.8,
      max.overlaps = 40,
      segment.size = 0.2,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    scale_colour_brewer(palette = "Set1", name = "Cluster") +
    labs(
      title = paste0(dataset_label, " \u2013 Gene UMAP (LRI / CCI LOGFC space)"),
      subtitle = paste0("Hierarchical clustering (k = ", k_clusters, ", cosine distance)"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  umap_file <- file.path(output_dir, paste0(dataset_label, "_gene_umap_lri.png"))
  ggsave(umap_file, plot = p_umap, width = 10, height = 8, dpi = 300)
  message("  UMAP saved to ", umap_file)

  cluster_file <- file.path(output_dir, paste0(dataset_label, "_gene_umap_clusters.csv"))
  write.csv(
    umap_df %>% select(Gene, Cluster, n_LRIs) %>% arrange(Cluster, Gene),
    cluster_file,
    row.names = FALSE
  )
  message("  Cluster table saved to ", cluster_file)

  invisible(list(umap_df = umap_df, plot = p_umap, clusters = clusters))
}

plot_gene_gene_lri_heatmap <- function(sim_mat,
                                       dataset_label,
                                       output_dir) {
  # Off-diagonal NAs → 0 for display; diagonal → NA (self-similarity not shown)
  sim_plot <- sim_mat
  sim_plot[is.na(sim_plot)] <- 0
  diag(sim_plot) <- NA

  sim_clust <- sim_mat
  diag(sim_clust) <- 1
  sim_clust[is.na(sim_clust)] <- 0

  n_genes <- nrow(sim_plot)
  font_sz <- max(4, min(7, 180 / n_genes))
  dist_clust <- as.dist(1 - sim_clust)

  heatmap_file <- file.path(output_dir, paste0(dataset_label, "_gene_gene_lri_heatmap.png"))
  px <- max(3000, n_genes * 22)
  png(heatmap_file, width = px, height = px, res = 300)
  pheatmap(
    sim_plot,
    color = colorRampPalette(c("#FF0000", "#000000", "#00FF00"))(100),
    breaks = seq(-1, 1, length.out = 101),
    na_col = "grey90",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = dist_clust,
    clustering_distance_cols = dist_clust,
    clustering_method = "complete",
    fontsize_row = font_sz,
    labels_col = rep("", ncol(sim_plot)),
    treeheight_col = 0,
    border_color = NA,
    legend = TRUE,
    main = paste0(dataset_label, " \u2013 Gene-Gene Cosine Similarity (LRI LOGFC)")
  )
  dev.off()
  message("  Heatmap saved to ", heatmap_file)
  invisible(sim_plot)
}

plot_gene_gene_lri_heatmap_umap <- function(malignant_list,
                                            dataset_label,
                                            output_dir,
                                            min_genes  = MIN_GENES_PER_CCI,
                                            k_clusters = 5,
                                            umap_seed  = 42) {
  message("\n=== GENES x LRI (within-cohort): ", dataset_label, " ===")

  mat <- build_lri_logfc_mat(malignant_list, min_genes = min_genes)
  sim_mat <- cosine_sim_na(mat)
  diag(sim_mat) <- NA

  plot_gene_gene_lri_heatmap(sim_mat, dataset_label, output_dir)
  umap_res <- plot_gene_umap_lri(
    mat, sim_mat, dataset_label, output_dir,
    k_clusters = k_clusters, umap_seed = umap_seed
  )

  cluster_summary <- umap_res$umap_df %>%
    group_by(Cluster) %>%
    summarise(
      n_genes = n(),
      genes = paste(sort(Gene), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(Cluster)
  print(cluster_summary)

  invisible(list(mat = mat, sim_mat = sim_mat, umap = umap_res))
}

compute_fisher_for_gene <- function(sc_obj,
                                    gene_name,
                                    LRI_GO_BP,
                                    malignant_type = MALIGNANT_CELLTYPE,
                                    min_go_size    = 5,
                                    min_detected   = 2) {
  detected <- filter.scDiffCom.cci_table_detected.for.malignant(
    sc_obj,
    malignant_celltype = malignant_type
  )

  raw <- sc_obj@cci_table_raw %>%
    filter(
      EMITTER_CELLTYPE %in% malignant_type |
        RECEIVER_CELLTYPE %in% malignant_type
    )

  if (nrow(detected) == 0) {
    return(NULL)
  }

  raw_go <- raw %>%
    inner_join(LRI_GO_BP, by = "LRI", relationship = "many-to-many")

  detected_go <- detected %>%
    inner_join(LRI_GO_BP, by = "LRI", relationship = "many-to-many")

  if (nrow(detected_go) == 0) {
    return(NULL)
  }

  raw_counts <- raw_go %>%
    dplyr::count(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, GO_NAME, name = "n_raw_in_go")

  raw_total_er <- raw_go %>%
    dplyr::count(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, name = "n_raw_total")

  det_counts <- detected_go %>%
    dplyr::count(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, GO_NAME, REGULATION,
          name = "n_det_in_go_reg")

  det_total_er_reg <- detected_go %>%
    dplyr::count(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, REGULATION, name = "n_det_total_reg")

  contexts <- det_counts %>%
    filter(n_det_in_go_reg >= min_detected) %>%
    left_join(raw_counts, by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "GO_NAME")) %>%
    left_join(raw_total_er, by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")) %>%
    left_join(det_total_er_reg,
              by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "REGULATION")) %>%
    filter(!is.na(n_raw_in_go), n_raw_in_go >= min_go_size)

  if (nrow(contexts) == 0) {
    return(NULL)
  }

  results <- vector("list", nrow(contexts))

  for (i in seq_len(nrow(contexts))) {
    a <- contexts$n_det_in_go_reg[i]
    b <- contexts$n_det_total_reg[i] - a
    c <- contexts$n_raw_in_go[i] - a
    d <- contexts$n_raw_total[i] - contexts$n_raw_in_go[i] - b

    if (any(c(a, b, c, d) < 0)) {
      next
    }

    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))

    results[[i]] <- data.frame(
      Gene = gene_name,
      GO_NAME = contexts$GO_NAME[i],
      EMITTER_CELLTYPE = contexts$EMITTER_CELLTYPE[i],
      RECEIVER_CELLTYPE = contexts$RECEIVER_CELLTYPE[i],
      REGULATION = contexts$REGULATION[i],
      odds_ratio = as.numeric(unname(ft$estimate)),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(results)
}

run_fisher_enrichment <- function(scDiffCom_list,
                                  dataset_label,
                                  output_dir,
                                  LRI_GO_BP,
                                  malignant_type = MALIGNANT_CELLTYPE,
                                  min_go_size    = 5,
                                  min_detected   = 2,
                                  save_every     = 10) {
  checkpoint_file <- file.path(output_dir, paste0(dataset_label, "_fisher_checkpoint.rds"))
  final_file <- file.path(output_dir, paste0(dataset_label, "_fisher_results_final.rds"))

  gene_names <- names(scDiffCom_list)
  n_genes <- length(gene_names)

  if (file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    completed <- checkpoint$completed
    results_so_far <- checkpoint$results
    message(sprintf("Resuming Fisher from checkpoint: %d / %d genes done",
                    length(completed), n_genes))
  } else {
    completed <- character(0)
    results_so_far <- list()
  }

  genes_todo <- setdiff(gene_names, completed)
  message(sprintf("Fisher enrichment: %d genes remaining (of %d)", length(genes_todo), n_genes))

  for (idx in seq_along(genes_todo)) {
    g <- genes_todo[idx]

    res <- tryCatch(
      compute_fisher_for_gene(
        sc_obj = scDiffCom_list[[g]],
        gene_name = g,
        LRI_GO_BP = LRI_GO_BP,
        malignant_type = malignant_type,
        min_go_size = min_go_size,
        min_detected = min_detected
      ),
      error = function(e) {
        message(sprintf("  [!] Error on gene %s: %s", g, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(res) && nrow(res) > 0) {
      results_so_far[[g]] <- res
    }
    completed <- c(completed, g)

    if (idx %% 10 == 0 || idx == length(genes_todo)) {
      message(sprintf("  Progress: %d / %d genes  (%d total rows so far)",
                      idx, length(genes_todo),
                      sum(vapply(results_so_far, nrow, integer(1)))))
    }

    if (idx %% save_every == 0) {
      saveRDS(list(completed = completed, results = results_so_far), checkpoint_file)
      message(sprintf("  [checkpoint saved at gene %d]", idx))
    }
  }

  all_results <- bind_rows(results_so_far)

  if (nrow(all_results) == 0) {
    message("No Fisher results produced — writing empty RDS.")
    all_results <- tibble(
      Gene = character(), GO_NAME = character(),
      EMITTER_CELLTYPE = character(), RECEIVER_CELLTYPE = character(),
      REGULATION = character(), odds_ratio = numeric(),
      p_value = numeric(), bh_p_val = numeric()
    )
  } else {
    all_results <- all_results %>%
      group_by(EMITTER_CELLTYPE, RECEIVER_CELLTYPE, GO_NAME, REGULATION) %>%
      mutate(bh_p_val = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>%
      select(Gene, GO_NAME, EMITTER_CELLTYPE, RECEIVER_CELLTYPE,
             REGULATION, odds_ratio, p_value, bh_p_val) %>%
      arrange(bh_p_val)
  }

  saveRDS(all_results, final_file)
  message(sprintf("Fisher results: %d rows written to %s", nrow(all_results), final_file))
  invisible(all_results)
}

plot_fisher_volcano <- function(fisher_res,
                                selected_genes,
                                dataset_label,
                                output_dir,
                                p_adj_cutoff = 0.05,
                                top_label_n  = 7,
                                bold_go_terms = character(0)) {
  selected_genes <- intersect(selected_genes, unique(fisher_res$Gene))
  if (length(selected_genes) == 0) {
    warning("No selected genes found in Fisher results — skipping volcano plots.")
    return(invisible(NULL))
  }

  fisher_res_sig <- fisher_res %>%
    filter(
      bh_p_val < p_adj_cutoff,
      is.finite(odds_ratio),
      is.finite(bh_p_val),
      REGULATION == "UP"
    ) %>%
    mutate(is_sig = bh_p_val < p_adj_cutoff)

  message(sprintf(
    "Significant UP results: %d / %d rows (%.1f%%)",
    nrow(fisher_res_sig), nrow(fisher_res),
    100 * nrow(fisher_res_sig) / max(nrow(fisher_res), 1)
  ))

  volcano_df <- fisher_res_sig %>%
    filter(Gene %in% selected_genes) %>%
    mutate(
      ER_PAIR = paste(EMITTER_CELLTYPE, "\u2192", RECEIVER_CELLTYPE),
      log2_OR = log2(odds_ratio),
      neg_log_p = -log10(bh_p_val),
      label_text = ifelse(is_sig, GO_NAME, "")
    ) %>%
    group_by(Gene) %>%
    mutate(
      rank_p = rank(-neg_log_p, ties.method = "first"),
      label_text = ifelse(is_sig & rank_p <= top_label_n, GO_NAME, "")
    ) %>%
    ungroup() %>%
    mutate(fontface = ifelse(GO_NAME %in% bold_go_terms, "bold", "plain"))

  if (nrow(volcano_df) == 0) {
    warning("No significant UP rows for selected genes — skipping volcano plots.")
    return(invisible(NULL))
  }

  all_pairs <- unique(volcano_df$ER_PAIR)
  n_pairs <- length(all_pairs)
  er_pal <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(min(8, n_pairs), "Dark2"))(n_pairs),
    all_pairs
  )

  p_volcano <- ggplot(volcano_df, aes(x = log2_OR, y = neg_log_p, colour = ER_PAIR)) +
    geom_hline(
      yintercept = -log10(p_adj_cutoff),
      linetype = "dashed", colour = "grey50", linewidth = 0.4
    ) +
    geom_point(size = 2.2, alpha = 0.85) +
    geom_text_repel(
      data = filter(volcano_df, label_text != ""),
      aes(label = label_text, colour = ER_PAIR, fontface = fontface),
      size = 2.4,
      max.overlaps = 20,
      segment.size = 0.3,
      segment.alpha = 0.5,
      show.legend = FALSE,
      box.padding = 0.3
    ) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = er_pal) +
    scale_x_continuous(
      name = expression(log[2] ~ "(Odds Ratio)"),
      labels = number_format(accuracy = 0.1)
    ) +
    scale_y_continuous(name = expression(-log[10] ~ "(BH p-value)")) +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      strip.background = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(colour = "grey95")
    ) +
    guides(colour = guide_legend(
      title = "Emitter \u2192 Receiver",
      nrow = ceiling(n_pairs / 3),
      override.aes = list(size = 3, alpha = 1)
    )) +
    labs(
      title = paste0("GO-BP Enrichment Volcano \u2014 ", dataset_label),
      subtitle = sprintf("BH p < %.2f  |  top %d GO labels per gene", p_adj_cutoff, top_label_n)
    )

  file_stub <- paste(selected_genes, collapse = "_")
  multi_file <- file.path(output_dir, paste0(file_stub, "_volcano.png"))
  ggsave(multi_file, plot = p_volcano, width = 8, height = 11, dpi = 300)
  message("  Multi-gene volcano saved to ", multi_file)

  for (current_gene in selected_genes) {
    gene_data <- volcano_df %>% filter(Gene == current_gene)
    if (nrow(gene_data) == 0) {
      next
    }

    p_individual <- ggplot(gene_data, aes(x = log2_OR, y = neg_log_p, colour = ER_PAIR)) +
      geom_hline(
        yintercept = -log10(p_adj_cutoff),
        linetype = "dashed", colour = "grey50", linewidth = 0.4
      ) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_text_repel(
        data = filter(gene_data, label_text != ""),
        aes(label = label_text, fontface = fontface),
        size = 3, max.overlaps = 20, segment.size = 0.3
      ) +
      scale_colour_manual(values = er_pal) +
      labs(
        title = paste(current_gene, "\u2014", dataset_label),
        x = expression(log[2] ~ "(Odds Ratio)"),
        y = expression(-log[10] ~ "(BH p-value)")
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.major = element_line(colour = "grey95")
      )

    ggsave(
      file.path(output_dir, paste0(current_gene, "_volcano_no_legend.png")),
      plot = p_individual,
      width = 7,
      height = 7,
      dpi = 300
    )
  }

  invisible(list(multi = p_volcano, data = volcano_df))
}
