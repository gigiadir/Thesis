#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(scales)
})

output_dir <- path.expand("~/Thesis/CCC/outputs/plots/paper")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

warning_log <- character(0)
saved_files <- character(0)

log_warn <- function(msg) {
  warning(msg, call. = FALSE)
  warning_log <<- c(warning_log, msg)
}

resolve_column <- function(df, candidates) {
  cols <- colnames(df)
  found <- candidates[candidates %in% cols]
  if (length(found) > 0) found[[1]] else NULL
}

coalesce_unknown <- function(x, unknown_label = "Unknown") {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- unknown_label
  x
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
  saved_files <<- c(saved_files, pdf_path, png_path)
  message("Saved: ", pdf_path)
  message("Saved: ", png_path)
}

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

load_seurat_pool <- function(paths) {
  pool <- list()
  for (p in paths) {
    path <- path.expand(p)
    if (!file.exists(path)) {
      next
    }
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    obj_names <- ls(e)
    objs <- mget(obj_names, envir = e)
    for (nm in names(objs)) {
      obj <- objs[[nm]]
      if (inherits(obj, "Seurat")) {
        pool[[nm]] <- obj
      }
    }
  }
  pool
}

pick_seurat_object <- function(pool, target_name) {
  if (length(pool) == 0) return(NULL)
  if (target_name %in% names(pool)) return(pool[[target_name]])

  target_norm <- normalize_name(target_name)
  name_norm <- normalize_name(names(pool))

  idx_exact <- which(name_norm == target_norm)
  if (length(idx_exact) > 0) return(pool[[idx_exact[[1]]]])

  idx_partial <- which(grepl(target_norm, name_norm) | grepl(name_norm, target_norm))
  if (length(idx_partial) > 0) return(pool[[idx_partial[[1]]]])

  NULL
}

collect_umap_reduction <- function(seurat_obj) {
  available <- names(seurat_obj@reductions)
  preferred <- c("umap", "UMAP", "integrated_umap", "harmony_umap")
  found <- preferred[preferred %in% available]
  if (length(found) > 0) return(found[[1]])
  NULL
}

build_qc_violin <- function(seurat_obj, dataset_label, group_col, group_label) {
  # If percent.mt is missing, compute it directly from mitochondrial gene prefix.
  if (!any(c("percent.mt", "percent_mt", "pct_mt", "mito_percent") %in% colnames(seurat_obj@meta.data))) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    message(dataset_label, ": computed percent.mt with PercentageFeatureSet(pattern = '^MT-').")
  }

  meta <- seurat_obj@meta.data %>% tibble::as_tibble(rownames = "cell_id")
  if (is.null(group_col) || !group_col %in% colnames(meta)) {
    return(
      placeholder_panel(
        paste0(dataset_label, ": QC violin"),
        paste0("Missing grouping column for ", group_label)
      )
    )
  }

  metric_map <- list(
    nFeature_RNA = c("nFeature_RNA", "nFeature", "RNA_feature"),
    nCount_RNA = c("nCount_RNA", "nCount", "RNA_count"),
    percent.mt = c("percent.mt", "percent_mt", "pct_mt", "mito_percent")
  )
  metric_cols <- lapply(metric_map, function(cands) resolve_column(meta, cands))
  missing_metrics <- names(metric_cols)[vapply(metric_cols, is.null, logical(1))]
  if (length(missing_metrics) > 0) {
    log_warn(
      paste0(
        dataset_label, ": Missing QC metric(s): ",
        paste(missing_metrics, collapse = ", ")
      )
    )
  }
  available_metrics <- names(metric_cols)[vapply(metric_cols, Negate(is.null), logical(1))]

  if (length(available_metrics) == 0) {
    return(
      placeholder_panel(
        paste0(dataset_label, ": QC violin"),
        "No QC metric columns found"
      )
    )
  }

  plot_df <- lapply(available_metrics, function(metric_name) {
    tibble(
      group_value = coalesce_unknown(meta[[group_col]]),
      metric = metric_name,
      value = suppressWarnings(as.numeric(meta[[metric_cols[[metric_name]]]]))
    )
  }) %>%
    bind_rows() %>%
    filter(!is.na(value))

  ggplot(plot_df, aes(x = group_value, y = value, fill = group_value)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85, linewidth = 0.2) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_y_continuous(labels = label_number(big.mark = ",")) +
    labs(
      title = paste0(dataset_label, ": QC metrics by ", group_label),
      x = group_label,
      y = "Value"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

build_umap_panel <- function(seurat_obj, dataset_label, color_col, panel_title) {
  meta <- seurat_obj@meta.data
  reduction_name <- collect_umap_reduction(seurat_obj)

  if (is.null(reduction_name)) {
    log_warn(paste0(dataset_label, ": No UMAP reduction found."))
    return(
      placeholder_panel(
        paste0(dataset_label, ": ", panel_title),
        "UMAP reduction is missing"
      )
    )
  }
  if (is.null(color_col) || !color_col %in% colnames(meta)) {
    log_warn(paste0(dataset_label, ": Missing metadata column for UMAP panel '", panel_title, "'."))
    return(
      placeholder_panel(
        paste0(dataset_label, ": ", panel_title),
        "Required grouping column is missing"
      )
    )
  }

  DimPlot(
    object = seurat_obj,
    reduction = reduction_name,
    group.by = color_col,
    pt.size = 0.15,
    raster = FALSE
  ) +
    labs(title = panel_title, subtitle = paste0("Colored by ", color_col)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}

extract_dataset <- function(pool, dataset_name) {
  obj <- pick_seurat_object(pool, dataset_name)
  if (is.null(obj)) {
    log_warn(paste0("Could not locate Seurat object for dataset '", dataset_name, "'."))
    return(NULL)
  }
  message("Using Seurat object for ", dataset_name, ": ", dataset_name)
  obj
}

resolve_dataset_columns <- function(seurat_obj) {
  meta <- seurat_obj@meta.data
  list(
    cell_type = resolve_column(meta, c("Cell_Type", "Cell_Labels", "CellType", "celltype", "cell_type")),
    project = resolve_column(meta, c("Project", "project", "Dataset", "dataset", "orig.ident")),
    patient = resolve_column(meta, c("Patient", "patient", "patient_id", "Sample", "sample", "ident"))
  )
}

# ---------------------- Load datasets ----------------------
candidate_paths <- c(
  "~/scObjects/HNSCC_Atlas.RData",
  "~/scObjects/Kurten_HNSC.RData"
)
seurat_pool <- load_seurat_pool(candidate_paths)
if (length(seurat_pool) == 0) {
  stop("No Seurat objects found in candidate paths: ", paste(candidate_paths, collapse = ", "))
}

hnscc_obj <- extract_dataset(seurat_pool, "HNSCC_Atlas")
if (is.null(hnscc_obj)) {
  hnscc_obj <- extract_dataset(seurat_pool, "HNSCC.Atlas")
}
kurten_obj <- extract_dataset(seurat_pool, "Kurten_HNSC")

if (is.null(hnscc_obj) || is.null(kurten_obj)) {
  stop("Could not resolve both required datasets (HNSCC_Atlas and Kurten_HNSC).")
}

hnscc_cols <- resolve_dataset_columns(hnscc_obj)
kurten_cols <- resolve_dataset_columns(kurten_obj)

# ---------------------- HNSCC_Atlas ----------------------
hnscc_qc <- build_qc_violin(
  seurat_obj = hnscc_obj,
  dataset_label = "HNSCC_Atlas",
  group_col = hnscc_cols$cell_type,
  group_label = "Cell_Type"
)

hnscc_umap_celltype <- build_umap_panel(
  seurat_obj = hnscc_obj,
  dataset_label = "HNSCC_Atlas",
  color_col = hnscc_cols$cell_type,
  panel_title = "HNSCC_Atlas UMAP by Cell_Type"
)
hnscc_umap_project <- build_umap_panel(
  seurat_obj = hnscc_obj,
  dataset_label = "HNSCC_Atlas",
  color_col = hnscc_cols$project,
  panel_title = "HNSCC_Atlas UMAP by Project"
)
hnscc_umap_combined <- hnscc_umap_celltype | hnscc_umap_project

save_figure_dual(
  hnscc_qc,
  "hnscc_atlas_qc_violin_by_cell_type",
  width = 16,
  height = 8.5
)
save_figure_dual(
  hnscc_umap_combined,
  "hnscc_atlas_umap_by_celltype_project",
  width = 18,
  height = 8
)

# ---------------------- Kurten_HNSC ----------------------
kurten_qc <- build_qc_violin(
  seurat_obj = kurten_obj,
  dataset_label = "Kurten_HNSC",
  group_col = kurten_cols$patient,
  group_label = "Patient"
)

kurten_umap_celltype <- build_umap_panel(
  seurat_obj = kurten_obj,
  dataset_label = "Kurten_HNSC",
  color_col = kurten_cols$cell_type,
  panel_title = "Kurten_HNSC UMAP by Cell_Type"
)
kurten_umap_project <- build_umap_panel(
  seurat_obj = kurten_obj,
  dataset_label = "Kurten_HNSC",
  color_col = kurten_cols$project,
  panel_title = "Kurten_HNSC UMAP by Project"
)
kurten_umap_patient <- build_umap_panel(
  seurat_obj = kurten_obj,
  dataset_label = "Kurten_HNSC",
  color_col = kurten_cols$patient,
  panel_title = "Kurten_HNSC UMAP by Patient"
)
kurten_umap_combined <- kurten_umap_celltype | kurten_umap_project | kurten_umap_patient

save_figure_dual(
  kurten_qc,
  "kurten_hnsc_qc_violin_by_patient",
  width = 16,
  height = 8.5
)
save_figure_dual(
  kurten_umap_combined,
  "kurten_hnsc_umap_by_celltype_project_patient",
  width = 24,
  height = 8
)

message("All requested figures completed.")
if (length(saved_files) > 0) {
  message("Generated files:")
  for (f in saved_files) message(" - ", f)
}
if (length(warning_log) > 0) {
  message("Warnings captured:")
  for (warn_msg in unique(warning_log)) message(" - ", warn_msg)
}
