#!/usr/bin/env Rscript

# HNSCC Atlas publication-ready diagnostic figure generator
# - Figure A: Data overview and composition
# - Figure B: QC assessment across projects
# - Figure C: UMAP integration diagnostics

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(scales)
})

input_rdata <- path.expand("~/scObjects/HNSCC_Atlas.RData")
output_dir <- path.expand("~/Thesis/CCC/outputs/plots/paper")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

warning_log <- character(0)

log_warn <- function(msg) {
  warning(msg, call. = FALSE)
  warning_log <<- c(warning_log, msg)
}

resolve_column <- function(df, candidates, required = FALSE, default_name = NULL) {
  cols <- colnames(df)
  found <- candidates[candidates %in% cols]
  if (length(found) > 0) {
    return(found[[1]])
  }
  if (isTRUE(required)) {
    log_warn(
      paste0(
        "Missing metadata column; expected one of: ",
        paste(candidates, collapse = ", "),
        if (!is.null(default_name)) paste0(". Using fallback label '", default_name, "'.") else "."
      )
    )
  }
  NULL
}

coalesce_unknown <- function(x, unknown_label = "Unknown") {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- unknown_label
  x
}

first_available <- function(named_values, fallback = NA_real_) {
  vals <- named_values[!is.na(named_values)]
  if (length(vals) == 0) fallback else vals[[1]]
}

extract_seurat_object <- function(rdata_path) {
  if (!file.exists(rdata_path)) {
    stop("Input .RData file does not exist: ", rdata_path)
  }

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  obj_names <- ls(e)
  objs <- mget(obj_names, envir = e)
  seurat_names <- names(Filter(function(x) inherits(x, "Seurat"), objs))

  if (length(seurat_names) == 0) {
    stop("No Seurat object found in file: ", rdata_path)
  }

  selected_name <- if ("HNSCC_Atlas" %in% seurat_names) "HNSCC_Atlas" else seurat_names[[1]]
  message("Using Seurat object: ", selected_name)
  objs[[selected_name]]
}

build_table_panel <- function(df, title, subtitle = NULL) {
  if (nrow(df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No data available", size = 4) +
        theme_void() +
        labs(title = title, subtitle = subtitle)
    )
  }

  tbl <- df %>%
    mutate(.row = row_number()) %>%
    pivot_longer(cols = -.row, names_to = "field", values_to = "value") %>%
    mutate(
      field = factor(field, levels = colnames(df)),
      value = as.character(value)
    )

  ggplot(tbl, aes(x = field, y = .row)) +
    geom_tile(fill = "grey97", color = "grey85", linewidth = 0.2) +
    geom_text(aes(label = value), size = 3) +
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.ticks = element_blank()
    )
}

placeholder_panel <- function(title, body_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = body_text, size = 4) +
    theme_void() +
    labs(title = title)
}

save_figure_dual <- function(plot_obj, base_name, width = 14, height = 8) {
  pdf_path <- file.path(output_dir, paste0(base_name, ".pdf"))
  png_path <- file.path(output_dir, paste0(base_name, ".png"))

  ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height, units = "in")
  ggsave(filename = png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300)
  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

prepare_counts_summary <- function(meta, project_col) {
  pre_candidates <- c("pre_qc_cells", "n_cells_pre_qc", "cells_before_qc", "preQC_cells", "Pre_QC_Cells")
  post_candidates <- c("post_qc_cells", "n_cells_post_qc", "cells_after_qc", "postQC_cells", "Post_QC_Cells")

  pre_col <- resolve_column(meta, pre_candidates)
  post_col <- resolve_column(meta, post_candidates)

  by_project <- meta %>%
    transmute(project = .data[[project_col]])

  observed_post <- by_project %>%
    count(project, name = "post_qc_cells_observed")

  if (is.null(pre_col)) {
    log_warn("Pre-QC cell counts were not found; creating NA placeholders by project.")
  }

  if (is.null(post_col)) {
    log_warn("Post-QC metadata counts were not found; using observed cell counts as post-QC totals.")
  }

  supplied_counts <- meta %>%
    transmute(
      project = .data[[project_col]],
      pre_qc_cells = if (!is.null(pre_col)) suppressWarnings(as.numeric(.data[[pre_col]])) else NA_real_,
      post_qc_cells_meta = if (!is.null(post_col)) suppressWarnings(as.numeric(.data[[post_col]])) else NA_real_
    ) %>%
    group_by(project) %>%
    summarise(
      pre_qc_cells = first_available(pre_qc_cells, NA_real_),
      post_qc_cells_meta = first_available(post_qc_cells_meta, NA_real_),
      .groups = "drop"
    )

  observed_post %>%
    left_join(supplied_counts, by = "project") %>%
    mutate(
      post_qc_cells = if_else(is.na(post_qc_cells_meta), post_qc_cells_observed, post_qc_cells_meta),
      qc_retention = if_else(!is.na(pre_qc_cells) & pre_qc_cells > 0, post_qc_cells / pre_qc_cells, NA_real_)
    ) %>%
    transmute(
      Project = as.character(project),
      pre_qc_cells = round(pre_qc_cells),
      post_qc_cells = round(post_qc_cells),
      qc_retention = percent(qc_retention, accuracy = 0.1)
    )
}

collect_umap_reduction <- function(seurat_obj) {
  available <- names(seurat_obj@reductions)
  preferred <- c("umap", "UMAP", "integrated_umap", "harmony_umap")
  found <- preferred[preferred %in% available]
  if (length(found) > 0) return(found[[1]])
  NULL
}

# ---------------------- Data Setup ----------------------
seurat_obj <- extract_seurat_object(input_rdata)
meta <- seurat_obj@meta.data %>%
  tibble::as_tibble(rownames = "cell_id")

project_col <- resolve_column(meta, c("Project", "project", "Dataset", "dataset", "orig.ident"), required = TRUE, default_name = "Project")
patient_col <- resolve_column(meta, c("Patient", "patient", "patient_id", "Sample", "sample", "ident"), required = TRUE, default_name = "Patient")
hpv_col <- resolve_column(meta, c("HPV", "HPV_status", "hpv_status", "HPVStatus"), required = FALSE)
tn_col <- resolve_column(meta, c("TumorNormal", "Tumor_Normal", "Tumor/Normal", "tissue_status", "Status"), required = FALSE)
cell_label_col <- resolve_column(meta, c("Cell_Labels", "CellLabels", "Cell_Type", "CellType", "celltype"), required = FALSE)

if (is.null(project_col)) {
  meta$Project <- "Unknown"
  project_col <- "Project"
}
if (is.null(patient_col)) {
  meta$Patient <- "Unknown"
  patient_col <- "Patient"
}
if (is.null(hpv_col)) {
  meta$HPV <- "Unknown"
  hpv_col <- "HPV"
}
if (is.null(tn_col)) {
  meta$Tumor_Normal <- "Unknown"
  tn_col <- "Tumor_Normal"
}

meta[[project_col]] <- coalesce_unknown(meta[[project_col]])
meta[[patient_col]] <- coalesce_unknown(meta[[patient_col]])
meta[[hpv_col]] <- coalesce_unknown(meta[[hpv_col]])
meta[[tn_col]] <- coalesce_unknown(meta[[tn_col]])

# ---------------------- Figure A ----------------------
metadata_summary <- meta %>%
  transmute(
    Project = .data[[project_col]],
    Patient = .data[[patient_col]],
    HPV = .data[[hpv_col]],
    Tumor_Normal = .data[[tn_col]]
  ) %>%
  count(Project, Patient, HPV, Tumor_Normal, name = "n_cells_post_qc") %>%
  arrange(desc(n_cells_post_qc))

project_level_summary <- metadata_summary %>%
  group_by(Project) %>%
  summarise(
    Patients = n_distinct(Patient),
    HPV_statuses = paste(sort(unique(HPV)), collapse = ", "),
    TumorNormal_statuses = paste(sort(unique(Tumor_Normal)), collapse = ", "),
    post_qc_cells = sum(n_cells_post_qc),
    .groups = "drop"
  ) %>%
  arrange(desc(post_qc_cells))

counts_summary <- prepare_counts_summary(meta, project_col) %>%
  arrange(desc(post_qc_cells))

overview_panel <- build_table_panel(
  project_level_summary %>% select(Project, Patients, HPV_statuses, TumorNormal_statuses, post_qc_cells),
  title = "Metadata Summary by Project",
  subtitle = "Project / Patient / HPV / Tumor-Normal composition"
)

qc_counts_panel <- build_table_panel(
  counts_summary,
  title = "Pre- and Post-QC Cell Counts by Project",
  subtitle = "Pre-QC values are placeholders when unavailable in metadata"
)

if (!is.null(cell_label_col)) {
  composition_df <- meta %>%
    transmute(Project = .data[[project_col]], Cell_Labels = .data[[cell_label_col]]) %>%
    mutate(Cell_Labels = coalesce_unknown(Cell_Labels)) %>%
    count(Project, Cell_Labels, name = "n_cells") %>%
    group_by(Project) %>%
    mutate(frac = n_cells / sum(n_cells)) %>%
    ungroup()

  composition_panel <- ggplot(composition_df, aes(x = Project, y = frac, fill = Cell_Labels)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Cell-Type Composition by Project",
      x = "Project",
      y = "Cell fraction",
      fill = "Cell_Labels"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
} else {
  log_warn("Cell label column not found. Figure A composition panel will be a placeholder.")
  composition_panel <- placeholder_panel(
    title = "Cell-Type Composition by Project",
    body_text = "Missing cell-type metadata column\n(expected: Cell_Labels or alias)"
  )
}

figure_a <- (overview_panel / qc_counts_panel) | composition_panel +
  plot_layout(widths = c(1.2, 1), heights = c(1, 1)) +
  plot_annotation(title = "Figure A: Data Overview and Composition")

save_figure_dual(figure_a, "figureA_overview_composition", width = 16, height = 9.5)

# ---------------------- Figure B ----------------------
qc_metric_candidates <- list(
  nCount_RNA = c("nCount_RNA", "nCount", "RNA_count"),
  nFeature_RNA = c("nFeature_RNA", "nFeature", "RNA_feature"),
  percent_mt = c("percent.mt", "percent_mt", "pct_mt", "mito_percent")
)

qc_col_map <- lapply(qc_metric_candidates, function(cands) resolve_column(meta, cands))
qc_available <- names(qc_col_map)[vapply(qc_col_map, Negate(is.null), logical(1))]

if (length(qc_available) == 0) {
  log_warn("No QC metric columns were found. Figure B will be a placeholder.")
  figure_b <- placeholder_panel(
    title = "Figure B: QC Assessment",
    body_text = "No QC columns found\n(expected nCount_RNA / nFeature_RNA / percent.mt)"
  )
} else {
  qc_long <- lapply(qc_available, function(metric_name) {
    col_name <- qc_col_map[[metric_name]]
    tibble(
      Project = meta[[project_col]],
      metric = metric_name,
      value = suppressWarnings(as.numeric(meta[[col_name]]))
    )
  }) %>%
    bind_rows() %>%
    filter(!is.na(value)) %>%
    mutate(metric = factor(metric, levels = qc_available))

  qc_global_panel <- ggplot(qc_long, aes(x = metric, y = value, fill = metric)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white", color = "grey20") +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(title = "Global QC distributions", x = NULL, y = "Value") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  qc_project_panel <- ggplot(qc_long, aes(x = metric, y = value, fill = metric)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.2) +
    facet_wrap(~ Project, scales = "free_y") +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(title = "QC distributions faceted by Project", x = NULL, y = "Value") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

  counts_for_qc <- counts_summary %>%
    mutate(
      pre_qc_cells_numeric = suppressWarnings(as.numeric(pre_qc_cells)),
      post_qc_cells_numeric = suppressWarnings(as.numeric(post_qc_cells))
    ) %>%
    select(Project, pre_qc_cells_numeric, post_qc_cells_numeric) %>%
    pivot_longer(
      cols = c(pre_qc_cells_numeric, post_qc_cells_numeric),
      names_to = "stage",
      values_to = "cells"
    ) %>%
    mutate(
      stage = recode(
        stage,
        pre_qc_cells_numeric = "Pre-QC",
        post_qc_cells_numeric = "Post-QC"
      )
    )

  qc_count_panel <- ggplot(counts_for_qc, aes(x = Project, y = cells, fill = stage)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.68, color = "white", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(is.na(cells), "NA", comma(round(cells)))),
      position = position_dodge(width = 0.72),
      vjust = -0.3,
      size = 3
    ) +
    labs(title = "Pre-/Post-QC count overlay by Project", x = "Project", y = "Cells", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  figure_b <- (qc_global_panel / qc_count_panel) | qc_project_panel +
    plot_layout(widths = c(1, 1.4), heights = c(1, 0.8)) +
    plot_annotation(title = "Figure B: Quality Control Assessment")
}

save_figure_dual(figure_b, "figureB_qc_diagnostics", width = 16, height = 10)

# ---------------------- Figure C ----------------------
umap_reduction <- collect_umap_reduction(seurat_obj)
if (is.null(umap_reduction)) {
  log_warn("No UMAP reduction found in Seurat object. Figure C will be a placeholder.")
  figure_c <- placeholder_panel(
    title = "Figure C: UMAP Embedding and Integration Diagnostics",
    body_text = "No UMAP reduction available\n(expected reduction name: umap or alias)"
  )
} else if (is.null(cell_label_col)) {
  log_warn("Cell label column not found. Figure C main panel uses unlabeled fallback.")
  p_main <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    pt.size = 0.15,
    cols = "grey55"
  ) +
    labs(
      title = "Global UMAP (unlabeled fallback)",
      subtitle = paste0("Reduction: ", umap_reduction)
    ) +
    theme_minimal(base_size = 11)

  p_proj <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    group.by = project_col,
    split.by = project_col,
    pt.size = 0.15,
    ncol = 3
  ) +
    labs(title = "UMAP split by Project") +
    theme_minimal(base_size = 10)

  p_patient <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    group.by = patient_col,
    split.by = patient_col,
    pt.size = 0.1,
    ncol = 4
  ) +
    labs(title = "UMAP split by Patient") +
    theme_minimal(base_size = 9)

  figure_c <- (p_main | p_proj) / p_patient +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = "Figure C: UMAP Embedding and Integration Diagnostics")
} else {
  p_main <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    group.by = cell_label_col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.18,
    raster = FALSE
  ) +
    labs(
      title = "Global integrated UMAP colored by Cell_Labels",
      subtitle = paste0("Reduction: ", umap_reduction)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_minimal(base_size = 11)

  p_proj <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    group.by = cell_label_col,
    split.by = project_col,
    pt.size = 0.12,
    ncol = 3,
    raster = FALSE
  ) +
    labs(title = "Companion UMAP split by Project") +
    theme_minimal(base_size = 10)

  p_patient <- DimPlot(
    object = seurat_obj,
    reduction = umap_reduction,
    group.by = cell_label_col,
    split.by = patient_col,
    pt.size = 0.1,
    ncol = 4,
    raster = FALSE
  ) +
    labs(title = "Companion UMAP split by Patient") +
    theme_minimal(base_size = 9)

  figure_c <- (p_main | p_proj) / p_patient +
    plot_layout(heights = c(1, 1.25)) +
    plot_annotation(title = "Figure C: UMAP Embedding and Integration Diagnostics")
}

save_figure_dual(figure_c, "figureC_umap_integration", width = 17, height = 11)

message("All requested figure exports completed.")
if (length(warning_log) > 0) {
  message("Warnings captured (handled gracefully):")
  for (warn_msg in unique(warning_log)) {
    message(" - ", warn_msg)
  }
}
