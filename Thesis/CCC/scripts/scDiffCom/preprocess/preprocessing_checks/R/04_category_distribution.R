# Category distribution QC (section 04).
# Requires: EXPR_LEVELS, EXPR_PAL, SPLIT_BASE_DIRS, split_source_labels.

read_grouped_file <- function(rds_path, dataset_name) {
  df <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(df) || !is.data.frame(df) || ncol(df) != 3) return(NULL)

  exp_col <- grep("_exp$", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(exp_col)) return(NULL)

  gene <- sub(paste0("_", dataset_name, "_grouped\\.rds$"), "", basename(rds_path))

  df %>%
    dplyr::select(patient_id, dplyr::all_of(exp_col)) %>%
    dplyr::rename(expression_level = dplyr::all_of(exp_col)) %>%
    dplyr::mutate(
      gene = gene,
      expression_level = as.character(expression_level)
    )
}

read_grouped_splits_long <- function(base_dir, dataset_name) {
  ds_dir <- file.path(path.expand(base_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    return(NULL)
  }

  rds_files <- list.files(
    ds_dir,
    pattern = paste0("_", dataset_name, "_grouped\\.rds$"),
    full.names = TRUE
  )
  if (length(rds_files) == 0L) {
    return(NULL)
  }

  agg <- rds_files %>%
    lapply(read_grouped_file, dataset_name = dataset_name) %>%
    Filter(Negate(is.null), .) %>%
    dplyr::bind_rows()

  if (nrow(agg) == 0L) {
    return(NULL)
  }
  agg
}

summarize_category_counts <- function(long_df, expr_levels = EXPR_LEVELS) {
  long_df %>%
    dplyr::filter(!is.na(patient_id), !is.na(expression_level)) %>%
    dplyr::mutate(
      patient_id = as.character(patient_id),
      expression_level = factor(expression_level, levels = expr_levels)
    ) %>%
    dplyr::count(patient_id, expression_level, name = "n") %>%
    tidyr::complete(patient_id, expression_level, fill = list(n = 0))
}

plot_split_category_per_patient <- function(count_df, ds, source_label, expr_pal = EXPR_PAL) {
  ggplot2::ggplot(count_df, ggplot2::aes(x = patient_id, y = n, fill = expression_level)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
      width = 0.8,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_fill_manual(values = expr_pal, drop = FALSE) +
    ggplot2::labs(
      title = paste0(ds, ": Split category distribution per patient"),
      subtitle = source_label,
      x = "patient_id",
      y = "Number of analyzed genes",
      fill = "Split\ncategory"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

plot_split_category_comparison <- function(count_df, ds, expr_pal = EXPR_PAL) {
  ggplot2::ggplot(count_df, ggplot2::aes(x = patient_id, y = n, fill = expression_level)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
      width = 0.8,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::facet_wrap(~ split_source, ncol = 1, scales = "free_x") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_fill_manual(values = expr_pal, drop = FALSE) +
    ggplot2::labs(
      title = paste0(ds, ": Split category distribution — method comparison"),
      x = "patient_id",
      y = "Number of analyzed genes",
      fill = "Split\ncategory"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
