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

comparison_facet_dims <- function(n_panels, ncol = 2L) {
  ncol <- max(1L, as.integer(ncol))
  nrow <- ceiling(n_panels / ncol)
  list(width = max(14, 7 * ncol), height = max(6, 5.5 * nrow))
}

comparison_panel_spacing <- function() {
  x_cm <- if (exists("CATEGORY_QC_COMPARISON_PANEL_SPACING_X", inherits = TRUE)) {
    CATEGORY_QC_COMPARISON_PANEL_SPACING_X
  } else {
    2.5
  }
  y_cm <- if (exists("CATEGORY_QC_COMPARISON_PANEL_SPACING_Y", inherits = TRUE)) {
    CATEGORY_QC_COMPARISON_PANEL_SPACING_Y
  } else {
    1.5
  }
  list(
    x = grid::unit(x_cm, "cm"),
    y = grid::unit(y_cm, "cm")
  )
}

collect_category_counts_by_source <- function(ds, sources = split_sources) {
  out <- list()
  for (src in sources) {
    base_dir <- SPLIT_BASE_DIRS[[src]]
    source_label <- split_source_labels[[src]]
    message("  Source: ", src, " (", source_label, ")")

    long_df <- read_grouped_splits_long(base_dir, ds)
    if (is.null(long_df)) {
      message("  No grouped .rds under ", file.path(base_dir, ds), " — skipping.")
      next
    }

    out[[src]] <- summarize_category_counts(long_df) %>%
      dplyr::mutate(split_source = factor(source_label, levels = source_label))
  }
  out
}

save_split_category_comparison <- function(ds,
                                           comparison_counts,
                                           out_dir = CATEGORY_QC_DIR,
                                           facet_ncol = CATEGORY_QC_COMPARISON_FACET_NCOL,
                                           expr_pal = EXPR_PAL) {
  if (length(comparison_counts) < 2L) {
    message("  Need >= 2 split methods for comparison — skipping ", ds, ".")
    return(invisible(NULL))
  }

  combined <- dplyr::bind_rows(comparison_counts)
  combined$split_source <- factor(
    combined$split_source,
    levels = vapply(comparison_counts, function(x) as.character(x$split_source[1]), character(1))
  )

  p_cmp <- plot_split_category_comparison(
    combined, ds, expr_pal = expr_pal, facet_ncol = facet_ncol
  )
  print(p_cmp)

  dims <- comparison_facet_dims(length(comparison_counts), ncol = facet_ncol)
  cmp_png <- file.path(out_dir, paste0(ds, "_split_category_comparison_facets.png"))
  ggplot2::ggsave(cmp_png, p_cmp, width = dims$width, height = dims$height, dpi = 300)
  message("  Saved comparison: ", cmp_png)
  invisible(cmp_png)
}

plot_split_category_comparison <- function(count_df,
                                           ds,
                                           expr_pal = EXPR_PAL,
                                           facet_ncol = 2L,
                                           panel_spacing = comparison_panel_spacing()) {
  n_facets <- length(unique(count_df$split_source))
  ncol_use <- min(as.integer(facet_ncol), max(1L, n_facets))

  ggplot2::ggplot(count_df, ggplot2::aes(x = patient_id, y = n, fill = expression_level)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.9, preserve = "single"),
      width = 0.8,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::facet_wrap(~ split_source, ncol = ncol_use, scales = "free_x") +
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
      panel.spacing.x = panel_spacing$x,
      panel.spacing.y = panel_spacing$y,
      strip.text = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(b = 6)
      ),
      strip.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )
}
