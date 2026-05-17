#!/usr/bin/env Rscript

# Build pseudobulk matrices from Seurat objects (malignant cells, detection-filtered genes).
#
# For each job in PSEUDOBULK_CONFIG:
#   (1) seurat_path   — .rds or .RData/.rda
#   (2) patient_col   — e.g. "Patient", "Sample", "patient_id"
#   (3) cell_type_column + malignant_cell_type — e.g. Cell_Type + "Tumor"
#
# Genes kept: detected (counts > 0) in > MIN_DETECT_FRAC_IN_MALIGNANT of malignant cells.
# Matrix values: mean normalized expression (data layer) per patient within malignant cells.
# Rows = genes, columns = patients.

suppressPackageStartupMessages({
  library(purrr)
  library(Seurat)
  library(Matrix)
})

SAVE_MATRICES <- TRUE
SAVE_PER_DATASET_RDS <- TRUE
SAVE_ALL_AS_RDATA <- TRUE
OUT_DIR <- path.expand("~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix")
MIN_DETECT_FRAC_IN_MALIGNANT <- 0.05

PSEUDOBULK_CONFIG <- list(
  list(
    seurat_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Kurten_HNSC.RData",
    dataset_name = "Kurten_HNSC",
    patient_col = "Patient",
    cell_type_column = "Cell_Type",
    malignant_cell_type = "Tumor"
  ),
  list(
    seurat_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Bill_HNSC.RData",
    dataset_name = "Bill_HNSC",
    patient_col = "Patient",
    cell_type_column = "Cell_Type",
    malignant_cell_type = "Tumor"
  ),
  list(
    seurat_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Puram_HNSC.RData",
    dataset_name = "Puram_HNSC",
    patient_col = "Patient",
    cell_type_column = "Cell_Type",
    malignant_cell_type = "Tumor"
  ),
  list(
    seurat_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Choi_HNSC.RData",
    dataset_name = "Choi_HNSC",
    patient_col = "Patient",
    cell_type_column = "Cell_Type",
    malignant_cell_type = "Tumor"
  )
  # list(
  #   seurat_path = "/path/to/cohort.RDS",
  #   dataset_name = "Some_cohort",
  #   patient_col = "ident",
  #   cell_type_column = "Consensus_Cell_Type",
  #   malignant_cell_type = c("Epithelial")
  # )
)

load_seurat_from_path <- function(ds_path) {
  ds_path <- path.expand(ds_path)
  if (!file.exists(ds_path)) {
    stop("Seurat file not found: ", ds_path)
  }
  ext <- tolower(tools::file_ext(ds_path))
  if (ext == "rds") {
    obj <- readRDS(ds_path)
  } else if (ext %in% c("rdata", "rda")) {
    if (requireNamespace("miceadds", quietly = TRUE)) {
      obj <- miceadds::load.Rdata(ds_path, "obj")
    } else {
      env <- new.env(parent = emptyenv())
      loaded <- load(ds_path, envir = env)
      seurat_idx <- vapply(
        loaded,
        function(nm) inherits(get(nm, envir = env), "Seurat"),
        logical(1)
      )
      if (!any(seurat_idx)) {
        stop("No Seurat object found in ", ds_path)
      }
      obj <- get(loaded[which(seurat_idx)[1L]], envir = env)
    }
  } else {
    stop("Unsupported file format: ", ext, ". Use .rds, .RData, or .rda")
  }
  if (!inherits(obj, "Seurat")) {
    stop("File does not contain a Seurat object: ", ds_path)
  }
  obj
}

get_assay_layer_or_slot <- function(obj, layer, slot) {
  assay <- Seurat::DefaultAssay(obj)
  tryCatch(
    SeuratObject::LayerData(obj, layer = layer, assay = assay),
    error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = slot)
    }
  )
}

build_maligexpr_patient_matrix <- function(
    obj,
    dataset_name,
    patient_col,
    cell_type_column,
    malignant_cell_type,
    min_detect_frac = MIN_DETECT_FRAC_IN_MALIGNANT
) {
  if (!patient_col %in% colnames(obj@meta.data)) {
    stop("Column '", patient_col, "' not found in metadata for ", dataset_name)
  }
  if (!cell_type_column %in% colnames(obj@meta.data)) {
    stop("Cell-type column '", cell_type_column, "' not found for ", dataset_name)
  }

  malig_vals <- as.character(malignant_cell_type)
  malig_cell_ids <- rownames(obj@meta.data)[
    as.character(obj@meta.data[[cell_type_column]]) %in% malig_vals
  ]
  malig_cell_ids <- intersect(malig_cell_ids, colnames(obj))
  if (length(malig_cell_ids) == 0L) {
    stop(
      "No malignant cells for ", dataset_name, " (",
      cell_type_column, " %in% ", paste(malig_vals, collapse = ", "), ")"
    )
  }

  obj_malignant <- subset(x = obj, cells = malig_cell_ids)
  malig_cells <- colnames(obj_malignant)
  meta_m <- obj_malignant@meta.data
  meta_full <- obj@meta.data

  counts_mat <- get_assay_layer_or_slot(obj_malignant, layer = "counts", slot = "counts")
  data_mat <- get_assay_layer_or_slot(obj_malignant, layer = "data", slot = "data")

  C <- counts_mat[, malig_cells, drop = FALSE]
  n_m <- ncol(C)
  detect_frac <- Matrix::rowSums(C > 0) / n_m
  genes_keep <- rownames(C)[detect_frac > min_detect_frac]
  genes_keep <- intersect(genes_keep, rownames(data_mat))
  if (length(genes_keep) == 0L) {
    stop(
      "No genes with detection fraction > ", min_detect_frac,
      " in malignant cells for ", dataset_name
    )
  }

  D <- data_mat[genes_keep, malig_cells, drop = FALSE]
  pat_by_cell <- setNames(as.character(meta_m[[patient_col]]), rownames(meta_m))
  patients <- sort(unique(c(
    as.character(meta_full[[patient_col]]),
    as.character(pat_by_cell)
  )))
  patients <- patients[!is.na(patients) & patients != ""]

  mat <- matrix(
    NA_real_,
    nrow = length(genes_keep),
    ncol = length(patients),
    dimnames = list(genes_keep, patients)
  )

  for (p in patients) {
    cc <- names(pat_by_cell)[pat_by_cell == p]
    if (length(cc) == 0L) next
    mat[, p] <- Matrix::rowMeans(D[, cc, drop = FALSE], na.rm = TRUE)
  }

  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}

write_matrix_dims <- function(mat, path) {
  n_patients <- ncol(mat)
  n_genes <- nrow(mat)
  lines <- c(
    paste0("n_patients\t", n_patients),
    paste0("n_genes\t", n_genes),
    "layout\trows = genes (detected in > min fraction of malignant cells), columns = patients; values = mean data layer in malignant cells"
  )
  writeLines(lines, path)
  message(
    "Dimensions: ", n_patients, " patients (columns) x ", n_genes,
    " genes (rows); wrote ", path
  )
}

normalize_pseudobulk_job <- function(job) {
  required <- c(
    "seurat_path", "dataset_name", "patient_col",
    "cell_type_column", "malignant_cell_type"
  )
  missing <- setdiff(required, names(job))
  if (length(missing) > 0L) {
    stop("PSEUDOBULK_CONFIG entry missing field(s): ", paste(missing, collapse = ", "))
  }
  job$malignant_cell_type <- as.character(job$malignant_cell_type)
  job
}

build_pseudobulk_from_config <- function(job) {
  job <- normalize_pseudobulk_job(job)
  ds <- job$dataset_name
  message("Loading ", ds, " from ", job$seurat_path)
  obj <- load_seurat_from_path(job$seurat_path)
  build_maligexpr_patient_matrix(
    obj = obj,
    dataset_name = ds,
    patient_col = job$patient_col,
    cell_type_column = job$cell_type_column,
    malignant_cell_type = job$malignant_cell_type,
    min_detect_frac = MIN_DETECT_FRAC_IN_MALIGNANT
  )
}

run_pseudobulk_pipeline <- function(
    config = PSEUDOBULK_CONFIG,
    out_dir = OUT_DIR,
    save_matrices = SAVE_MATRICES,
    save_per_dataset_rds = SAVE_PER_DATASET_RDS,
    save_all_as_rdata = SAVE_ALL_AS_RDATA,
    min_detect_frac = MIN_DETECT_FRAC_IN_MALIGNANT
) {
  if (isTRUE(save_matrices)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  expr_mats <- purrr::map(config, function(job) {
    ds <- job$dataset_name
    if (is.null(ds) || length(ds) == 0L || !nzchar(ds)) {
      ds <- tools::file_path_sans_ext(basename(job$seurat_path))
    }
    tryCatch(
      {
        job <- normalize_pseudobulk_job(job)
        obj <- load_seurat_from_path(job$seurat_path)
        build_maligexpr_patient_matrix(
          obj = obj,
          dataset_name = ds,
          patient_col = job$patient_col,
          cell_type_column = job$cell_type_column,
          malignant_cell_type = job$malignant_cell_type,
          min_detect_frac = min_detect_frac
        )
      },
      error = function(e) {
        message("Failed to build pseudobulk for ", ds, ": ", conditionMessage(e))
        NULL
      }
    )
  })
  names(expr_mats) <- vapply(config, function(j) j$dataset_name, character(1))
  expr_mats <- purrr::compact(expr_mats)

  message(
    "Built pseudobulk matrices for ", length(expr_mats), " dataset(s): ",
    paste(names(expr_mats), collapse = ", ")
  )

  if (isTRUE(save_matrices) && length(expr_mats) > 0) {
    if (isTRUE(save_per_dataset_rds)) {
      purrr::iwalk(expr_mats, function(mat, ds) {
        rds_path <- file.path(out_dir, paste0(ds, "_pseudobulk_matrix.rds"))
        saveRDS(mat, file = rds_path)
        dims_path <- file.path(out_dir, paste0(ds, "_pseudobulk_matrix.rds_dims.txt"))
        write_matrix_dims(mat, dims_path)
      })
      message("Saved per-dataset matrices and dimension summaries to: ", out_dir)
    }

    if (isTRUE(save_all_as_rdata)) {
      save(expr_mats, file = file.path(out_dir, "expr_mats_gene_by_patient_mean_expr.RData"))
      message("Saved all matrices as .RData to: ", out_dir)
    }
  }

  invisible(expr_mats)
}

if (sys.nframe() == 0L) {
  run_pseudobulk_pipeline()
}
