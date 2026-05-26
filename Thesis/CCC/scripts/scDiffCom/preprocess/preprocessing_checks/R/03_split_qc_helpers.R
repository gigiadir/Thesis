# Split QC helpers for dispersion (06) and profile PCA (07).
# Requires config: LABEL_NUM, ZSCORE_SPLITS_DIR, RANKGENES_DIR, RESIDUAL_SPLITS_DIR,
# MIN_PATIENTS, MAX_PATIENT_NA_FRAC, N_PC, LABEL_LOW_COR_QUANTILE, N_LABEL_GENES.

row_dispersion <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    return(matrixStats::rowSds(mat, na.rm = TRUE))
  }
  apply(mat, 1, stats::sd, na.rm = TRUE)
}

encode_splits <- function(splits_chr) {
  out <- matrix(NA_real_, nrow = nrow(splits_chr), ncol = ncol(splits_chr),
                dimnames = dimnames(splits_chr))
  for (lev in names(LABEL_NUM)) {
    out[splits_chr == lev] <- LABEL_NUM[[lev]]
  }
  storage.mode(out) <- "double"
  out
}

row_dispersion_expr <- function(mat) {
  row_dispersion(mat)
}

row_dispersion_split <- function(splits_chr) {
  row_dispersion(encode_splits(splits_chr))
}

align_matrices <- function(expr_mat, splits_mat) {
  genes <- intersect(rownames(expr_mat), rownames(splits_mat))
  patients <- intersect(colnames(expr_mat), colnames(splits_mat))
  if (length(genes) == 0L || length(patients) == 0L) {
    return(NULL)
  }
  list(
    expr = expr_mat[genes, patients, drop = FALSE],
    splits = splits_mat[genes, patients, drop = FALSE]
  )
}

compare_dispersion <- function(sd_expr, sd_split) {
  n <- length(sd_expr)
  if (is.null(names(sd_expr))) names(sd_expr) <- seq_len(n)
  valid <- is.finite(sd_expr) & is.finite(sd_split) & sd_expr > 0
  data.frame(
    gene = names(sd_expr)[valid],
    sd_expr = unname(sd_expr[valid]),
    sd_split = unname(sd_split[valid]),
    ratio = unname(sd_split[valid] / sd_expr[valid]),
    stringsAsFactors = FALSE
  )
}

load_splits_from_grouped_dir <- function(base_dir, dataset_name) {
  ds_dir <- file.path(path.expand(base_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    return(NULL)
  }

  files <- list.files(
    ds_dir,
    pattern = paste0("_", dataset_name, "_grouped\\.rds$"),
    full.names = TRUE
  )
  if (length(files) == 0L) {
    return(NULL)
  }

  patients <- NULL
  split_rows <- list()

  for (f in files) {
    df <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(df) || !is.data.frame(df)) next

    exp_col <- grep("_exp$", colnames(df), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(exp_col)) next

    gene <- sub(paste0("_", dataset_name, "_grouped\\.rds$"), "", basename(f))
    if (is.null(patients)) {
      patients <- as.character(df$patient_id)
    }
    vec <- setNames(as.character(df[[exp_col]]), as.character(df$patient_id))
    split_rows[[gene]] <- vec[patients]
  }

  if (length(split_rows) == 0L) {
    return(NULL)
  }

  genes_found <- names(split_rows)
  splits <- matrix(NA_character_, nrow = length(genes_found), ncol = length(patients),
                   dimnames = list(genes_found, patients))
  for (g in genes_found) {
    splits[g, ] <- split_rows[[g]]
  }
  splits
}

load_splits_rankgenes <- function(dataset_name, rankgenes_dir, expr_mat) {
  ds_dir <- file.path(path.expand(rankgenes_dir), dataset_name)
  all_splits_path <- file.path(ds_dir, paste0(dataset_name, ALL_SPLITS_SUFFIX))

  if (file.exists(all_splits_path)) {
    message("  Loading RankGenes splits: ", all_splits_path)
    return(readRDS(all_splits_path))
  }

  splits <- tryCatch(
    load_splits_from_rankgenes_dir(rankgenes_dir, dataset_name),
    error = function(e) NULL
  )
  if (!is.null(splits)) {
    message("  Loaded RankGenes splits from grouped .rds files in ", ds_dir)
    return(splits)
  }

  message("  RankGenes splits not found; building from pseudobulk (column-wise rank tertiles).")
  rank_mat <- compute_rank_matrix(expr_mat)
  build_patient_splits_matrix(expr_mat, rank_mat = rank_mat)
}

load_splits_residual <- function(dataset_name, residual_dir, expr_mat, n_pc = 2L) {
  ds_dir <- file.path(path.expand(residual_dir), dataset_name)
  all_splits_path <- file.path(ds_dir, paste0(dataset_name, ALL_SPLITS_SUFFIX))

  if (file.exists(all_splits_path)) {
    message("  Loading Residual splits: ", all_splits_path)
    return(readRDS(all_splits_path))
  }

  splits <- tryCatch(
    load_splits_from_residual_dir(residual_dir, dataset_name),
    error = function(e) NULL
  )
  if (!is.null(splits)) {
    message("  Loaded Residual splits from grouped .rds files in ", ds_dir)
    return(splits)
  }

  message("  Residual splits not found; building from pseudobulk (PC-regressed tertiles).")
  build_residual_splits_matrix(expr_mat, n_pc = n_pc)$splits
}

load_splits_for_source <- function(source, dataset_name, expr_mat) {
  switch(
    source,
    patient_zscore = load_splits_from_grouped_dir(ZSCORE_SPLITS_DIR, dataset_name),
    rankgenes = load_splits_rankgenes(dataset_name, RANKGENES_DIR, expr_mat),
    residual = load_splits_residual(dataset_name, RESIDUAL_SPLITS_DIR, expr_mat),
    stop("Unknown split source: ", source)
  )
}

plot_dispersion_scatter <- function(cmp_df, ds, source, source_label) {
  n_genes <- nrow(cmp_df)
  pct_split_lt <- 100 * mean(cmp_df$sd_split < cmp_df$sd_expr)
  pct_split_gt <- 100 * mean(cmp_df$sd_split > cmp_df$sd_expr)
  pct_split_zero <- 100 * mean(cmp_df$sd_split == 0)

  subtitle <- sprintf(
    "%d genes | split SD < expr SD: %.1f%% | split SD > expr SD: %.1f%% | split SD = 0: %.1f%%",
    n_genes, pct_split_lt, pct_split_gt, pct_split_zero
  )

  ggplot2::ggplot(cmp_df, ggplot2::aes(x = sd_expr, y = sd_split)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(alpha = 0.35, size = 1.2, color = "#4E79A7") +
    ggplot2::scale_x_continuous(trans = "log1p") +
    ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::labs(
      title = paste0(ds, ": Gene dispersion — expression vs split"),
      subtitle = subtitle,
      x = "Pseudobulk row SD (expression space)",
      y = "Encoded split row SD (split space)",
      caption = source_label
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, color = "grey40"))
}

run_dispersion_stats <- function(cmp_df) {
  wt <- stats::wilcox.test(cmp_df$sd_split, cmp_df$sd_expr, paired = TRUE, alternative = "two.sided")
  wt_greater <- stats::wilcox.test(cmp_df$sd_split, cmp_df$sd_expr, paired = TRUE, alternative = "greater")
  wt_less <- stats::wilcox.test(cmp_df$sd_split, cmp_df$sd_expr, paired = TRUE, alternative = "less")
  ct <- suppressWarnings(stats::cor.test(cmp_df$sd_expr, cmp_df$sd_split, method = "spearman"))

  list(
    n_genes = nrow(cmp_df),
    median_sd_expr = median(cmp_df$sd_expr),
    median_sd_split = median(cmp_df$sd_split),
    median_ratio = median(cmp_df$ratio),
    wilcox_p = wt$p.value,
    wilcox_p_split_gt_expr = wt_greater$p.value,
    wilcox_p_split_lt_expr = wt_less$p.value,
    spearman_rho = unname(ct$estimate),
    spearman_p = ct$p.value
  )
}

row_zscore_matrix <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    rm <- matrixStats::rowMeans2(mat, na.rm = TRUE)
    rsd <- matrixStats::rowSds(mat, na.rm = TRUE)
    rsd[!is.finite(rsd) | rsd == 0] <- NA_real_
    z <- sweep(mat, 1, rm, FUN = "-")
    z <- sweep(z, 1, rsd, FUN = "/")
    return(z)
  }
  t(scale(t(mat)))
}

genes_with_usable_rows <- function(mat, min_patients = MIN_PATIENTS) {
  apply(mat, 1, function(r) {
    ok <- is.finite(r)
    sum(ok) >= min_patients && stats::sd(r[ok]) > 0
  })
}

drop_sparse_patients <- function(mat, max_na_frac = MAX_PATIENT_NA_FRAC,
                               min_patients = MIN_PATIENTS, label = "matrix") {
  if (ncol(mat) == 0L) return(mat)

  na_frac <- colMeans(!is.finite(mat))
  keep <- na_frac < max_na_frac
  dropped <- colnames(mat)[!keep]

  if (length(dropped) > 0L) {
    message(sprintf(
      "  Dropped %d patient(s) from %s (NA fraction >= %.0f%%): %s",
      length(dropped), label, 100 * max_na_frac,
      paste(dropped, collapse = ", ")
    ))
  }

  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < min_patients) {
    warning(
      "Only ", ncol(mat), " patient(s) left in ", label,
      " after NA filtering (need >= ", min_patients, ")."
    )
  }
  mat
}

drop_sparse_patients_pair <- function(z_expr, z_split,
                                    max_na_frac = MAX_PATIENT_NA_FRAC,
                                    min_patients = MIN_PATIENTS) {
  patients <- intersect(colnames(z_expr), colnames(z_split))
  if (length(patients) == 0L) {
    stop("No shared patients between expression and split matrices.")
  }

  na_expr <- colMeans(!is.finite(z_expr[, patients, drop = FALSE]))
  na_split <- colMeans(!is.finite(z_split[, patients, drop = FALSE]))
  keep <- (na_expr < max_na_frac) & (na_split < max_na_frac)
  dropped <- patients[!keep]

  if (length(dropped) > 0L) {
    message(sprintf(
      "  Dropped %d patient(s) from paired PCA (NA in expr and/or split): %s",
      length(dropped), paste(dropped, collapse = ", ")
    ))
  }

  patients <- patients[keep]
  if (length(patients) < min_patients) {
    stop(
      "Fewer than ", min_patients, " patients after paired NA filtering (",
      length(patients), " left)."
    )
  }

  list(
    expr = z_expr[, patients, drop = FALSE],
    split = z_split[, patients, drop = FALSE],
    patients = patients
  )
}

prepare_matrix_for_pca <- function(z_mat, max_na_frac = MAX_PATIENT_NA_FRAC,
                                   min_patients = MIN_PATIENTS,
                                   impute_residual_na = TRUE, label = "PCA") {
  z <- drop_sparse_patients(z_mat, max_na_frac, min_patients, label)

  keep <- genes_with_usable_rows(z, min_patients)
  z <- z[keep, , drop = FALSE]

  gene_ok <- apply(z, 1, function(r) all(is.finite(r)) && stats::sd(r) > 0)
  z_complete <- z[gene_ok, , drop = FALSE]

  if (nrow(z_complete) >= 3L) {
    return(list(z = z_complete, imputed = FALSE, n_genes = nrow(z_complete)))
  }

  if (!impute_residual_na) {
    stop("Fewer than 3 genes with complete rows after patient filtering for ", label, ".")
  }

  n_imputed <- sum(!is.finite(z))
  z[!is.finite(z)] <- 0
  gene_ok <- apply(z, 1, function(r) stats::sd(r) > 0)
  z <- z[gene_ok, , drop = FALSE]

  if (nrow(z) < 3L) {
    stop("Fewer than 3 genes with usable z-scores for ", label, " after imputation.")
  }

  if (n_imputed > 0L) {
    message(sprintf(
      "  %s: imputed %d NA cell(s) with 0 (neutral z-score); %d genes for PCA.",
      label, n_imputed, nrow(z)
    ))
  }

  list(z = z, imputed = TRUE, n_genes = nrow(z))
}

pca_genes <- function(z_mat, n_pc = N_PC, label = "PCA") {
  prep <- prepare_matrix_for_pca(z_mat, label = label)
  z <- prep$z

  pc <- stats::prcomp(z, center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pc$x))
  var_expl <- (pc$sdev^2) / sum(pc$sdev^2)
  list(
    coords = pc$x[, seq_len(n_pc), drop = FALSE],
    var_expl = var_expl,
    rotation = pc$rotation[, seq_len(n_pc), drop = FALSE],
    n_genes = nrow(z),
    n_patients = ncol(z),
    imputed = prep$imputed
  )
}

row_profile_correlation <- function(z_expr, z_split) {
  genes <- intersect(rownames(z_expr), rownames(z_split))
  vapply(genes, function(g) {
    xe <- z_expr[g, ]
    xs <- z_split[g, ]
    ok <- is.finite(xe) & is.finite(xs)
    if (sum(ok) < 3L) return(NA_real_)
    suppressWarnings(stats::cor(xe[ok], xs[ok], method = "pearson"))
  }, numeric(1))
}

build_profile_pca_df <- function(pca_res, space_label) {
  df <- as.data.frame(pca_res$coords)
  colnames(df) <- paste0("PC", seq_len(ncol(df)))
  df$gene <- rownames(pca_res$coords)
  df$space <- space_label
  df$var_pc1 <- round(100 * pca_res$var_expl[1], 1)
  if (length(pca_res$var_expl) >= 2L) {
    df$var_pc2 <- round(100 * pca_res$var_expl[2], 1)
  } else {
    df$var_pc2 <- NA_real_
  }
  df
}

plot_profile_pca <- function(plot_df, df_expr, df_split, ds, source_label, cor_col = "profile_cor") {
  low_cor_thr <- stats::quantile(plot_df[[cor_col]], probs = LABEL_LOW_COR_QUANTILE, na.rm = TRUE)
  plot_df <- plot_df %>%
    dplyr::mutate(
      low_split_informativeness = !is.na(.data[[cor_col]]) &
        .data[[cor_col]] <= low_cor_thr
    )

  label_df <- plot_df %>%
    dplyr::filter(.data$space == "Expression", .data$low_split_informativeness) %>%
    dplyr::arrange(.data[[cor_col]]) %>%
    dplyr::slice_head(n = N_LABEL_GENES)

  p_expr <- ggplot2::ggplot(
    dplyr::filter(plot_df, .data$space == "Expression"),
    ggplot2::aes(x = PC1, y = PC2, color = .data[[cor_col]])
  ) +
    ggplot2::geom_point(alpha = 0.45, size = 1.4) +
    ggplot2::geom_text(
      data = label_df, ggplot2::aes(label = gene),
      size = 2.5, vjust = -0.6, check_overlap = TRUE, show.legend = FALSE
    ) +
    ggplot2::scale_color_viridis_c(option = "magma", name = "Expr–split\ncorrelation", na.value = "grey80") +
    ggplot2::labs(
      title = "Expression space",
      subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%% variance", df_expr$var_pc1[1], df_expr$var_pc2[1]),
      x = "PC1", y = "PC2"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")

  p_split <- ggplot2::ggplot(
    dplyr::filter(plot_df, .data$space == "Split"),
    ggplot2::aes(x = PC1, y = PC2, color = .data[[cor_col]])
  ) +
    ggplot2::geom_point(alpha = 0.45, size = 1.4) +
    ggplot2::scale_color_viridis_c(option = "magma", name = "Expr–split\ncorrelation", na.value = "grey80") +
    ggplot2::labs(
      title = "Split space",
      subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%% variance", df_split$var_pc1[1], df_split$var_pc2[1]),
      x = "PC1", y = "PC2"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")

  (p_expr | p_split) +
    patchwork::plot_annotation(
      title = paste0(ds, ": Gene patient-wise profiles (row z-scored)"),
      subtitle = paste0(
        source_label, " | dark labels = bottom ", 100 * LABEL_LOW_COR_QUANTILE,
        "% expr–split correlation (least informative splits)"
      ),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
    )
}

plot_profile_cor_hist <- function(metrics_df, ds, source_label) {
  ggplot2::ggplot(metrics_df, ggplot2::aes(x = profile_cor)) +
    ggplot2::geom_histogram(bins = 40, fill = "#4E79A7", color = "white", linewidth = 0.2) +
    ggplot2::geom_vline(
      xintercept = stats::quantile(metrics_df$profile_cor, LABEL_LOW_COR_QUANTILE, na.rm = TRUE),
      linetype = "dashed", color = "#E15759"
    ) +
    ggplot2::labs(
      title = paste0(ds, ": Expression vs split patient-profile correlation"),
      subtitle = source_label,
      x = "Per-gene correlation (row z expr vs row z split)",
      y = "Number of genes"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
