#!/usr/bin/env Rscript

# Build pseudobulk matrices from Seurat objects (malignant cells).
#
# Genes kept: union of
#   (1) detected (counts > 0) in > min_detect_frac of malignant cells, and
#   (2) force-included genes (default: scDiffCom panel + optional --gene_list / config).
#
# Matrix values: mean normalized expression (data layer) per patient within malignant cells.
# Rows = genes, columns = patients.
#
# CLI examples:
#   Rscript createPseudobulkMatrix.R
#   Rscript createPseudobulkMatrix.R --dataset_name Kurten_HNSC
#   Rscript createPseudobulkMatrix.R --gene_list ~/genes.txt --no_include_panel

suppressPackageStartupMessages({
  library(purrr)
  library(Seurat)
  library(Matrix)
  library(optparse)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(script_dir, "rankGenesSplitUtils.R"))

SAVE_MATRICES <- TRUE
SAVE_PER_DATASET_RDS <- TRUE
SAVE_ALL_AS_RDATA <- TRUE
OUT_DIR <- path.expand("~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix")
MIN_DETECT_FRAC_IN_MALIGNANT <- 0.05
GENE_SOURCES_SUFFIX <- "_pseudobulk_gene_sources.tsv"

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

resolve_force_include_genes <- function(
    force_include_genes = NULL,
    gene_list_path = NULL,
    include_panel = TRUE,
    job_force_include_genes = NULL
) {
  parts <- list()
  if (isTRUE(include_panel)) {
    parts <- c(parts, list(SCDIFFCOM_GENE_PANEL))
  }
  if (!is.null(gene_list_path) && nzchar(gene_list_path)) {
    parts <- c(parts, list(load_gene_list(gene_list_path)))
  }
  if (!is.null(force_include_genes) && length(force_include_genes) > 0L) {
    parts <- c(parts, list(as.character(force_include_genes)))
  }
  if (!is.null(job_force_include_genes) && length(job_force_include_genes) > 0L) {
    parts <- c(parts, list(as.character(job_force_include_genes)))
  }
  if (length(parts) == 0L) {
    return(character())
  }
  unique(unlist(parts, use.names = FALSE))
}

build_gene_source_table <- function(genes_keep, detected, force_in_assay) {
  detected <- intersect(genes_keep, detected)
  force_only <- setdiff(force_in_assay, detected)
  both <- intersect(detected, force_in_assay)
  rows <- list(
    data.frame(gene = both, source = "both", stringsAsFactors = FALSE),
    data.frame(gene = setdiff(detected, both), source = "detected", stringsAsFactors = FALSE),
    data.frame(gene = force_only, source = "force_include", stringsAsFactors = FALSE)
  )
  rows <- rows[vapply(rows, nrow, integer(1)) > 0L]
  if (length(rows) == 0L) {
    return(data.frame(gene = character(), source = character(), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$gene), , drop = FALSE]
}

build_maligexpr_patient_matrix <- function(
    obj,
    dataset_name,
    patient_col,
    cell_type_column,
    malignant_cell_type,
    min_detect_frac = MIN_DETECT_FRAC_IN_MALIGNANT,
    force_include_genes = character(),
    include_panel = TRUE,
    gene_list_path = NULL,
    job_force_include_genes = NULL
) {
  if (!patient_col %in% colnames(obj@meta.data)) {
    stop("Column '", patient_col, "' not found in metadata for ", dataset_name)
  }
  if (!cell_type_column %in% colnames(obj@meta.data)) {
    stop("Cell-type column '", cell_type_column, "' not found for ", dataset_name)
  }

  force_genes <- resolve_force_include_genes(
    force_include_genes = force_include_genes,
    gene_list_path = gene_list_path,
    include_panel = include_panel,
    job_force_include_genes = job_force_include_genes
  )

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
  detected <- rownames(C)[detect_frac > min_detect_frac]
  detected <- intersect(detected, rownames(data_mat))

  force_in_assay <- intersect(force_genes, rownames(data_mat))
  missing_force <- setdiff(force_genes, force_in_assay)
  if (length(missing_force) > 0L) {
    warning(
      dataset_name, ": ", length(missing_force),
      " force-include gene(s) not in assay and skipped: ",
      paste(head(missing_force, 20L), collapse = ", "),
      if (length(missing_force) > 20L) " ..." else "",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  genes_keep <- union(detected, force_in_assay)
  if (length(genes_keep) == 0L) {
    stop(
      "No genes to keep for ", dataset_name,
      " (detection > ", min_detect_frac, " or force-include)."
    )
  }

  n_force_only <- length(setdiff(force_in_assay, detected))
  message(
    dataset_name, ": ", length(detected), " detected + ",
    n_force_only, " force-only = ", length(genes_keep), " genes total"
  )

  gene_sources <- build_gene_source_table(genes_keep, detected, force_in_assay)

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

  list(matrix = mat, gene_sources = gene_sources)
}

write_matrix_dims <- function(mat, path, n_force_only = NA_integer_) {
  n_patients <- ncol(mat)
  n_genes <- nrow(mat)
  layout_note <- paste(
    "rows = genes (detected > min fraction OR force-included),",
    "columns = patients; values = mean data layer in malignant cells"
  )
  lines <- c(
    paste0("n_patients\t", n_patients),
    paste0("n_genes\t", n_genes),
    paste0("layout\t", layout_note)
  )
  if (!is.na(n_force_only)) {
    lines <- c(lines, paste0("n_force_only\t", n_force_only))
  }
  writeLines(lines, path)
  message(
    "Dimensions: ", n_patients, " patients (columns) x ", n_genes,
    " genes (rows); wrote ", path
  )
}

write_gene_sources_tsv <- function(gene_sources, path) {
  utils::write.table(
    gene_sources,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  message("Wrote gene source table: ", path)
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

run_pseudobulk_pipeline <- function(
    config = PSEUDOBULK_CONFIG,
    out_dir = OUT_DIR,
    save_matrices = SAVE_MATRICES,
    save_per_dataset_rds = SAVE_PER_DATASET_RDS,
    save_all_as_rdata = SAVE_ALL_AS_RDATA,
    min_detect_frac = MIN_DETECT_FRAC_IN_MALIGNANT,
    force_include_genes = NULL,
    gene_list_path = NULL,
    include_panel = TRUE,
    dataset_name_filter = NULL
) {
  if (!is.null(dataset_name_filter) && nzchar(dataset_name_filter)) {
    config <- config[vapply(config, function(j) j$dataset_name == dataset_name_filter, logical(1))]
    if (length(config) == 0L) {
      stop("No PSEUDOBULK_CONFIG entry for dataset_name: ", dataset_name_filter)
    }
  }

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
        job_include_panel <- if (!is.null(job$include_panel)) {
          isTRUE(job$include_panel)
        } else {
          include_panel
        }
        obj <- load_seurat_from_path(job$seurat_path)
        res <- build_maligexpr_patient_matrix(
          obj = obj,
          dataset_name = ds,
          patient_col = job$patient_col,
          cell_type_column = job$cell_type_column,
          malignant_cell_type = job$malignant_cell_type,
          min_detect_frac = min_detect_frac,
          force_include_genes = force_include_genes,
          gene_list_path = gene_list_path,
          include_panel = job_include_panel,
          job_force_include_genes = job$force_include_genes
        )
        attr(res$matrix, "gene_sources") <- res$gene_sources
        attr(res$matrix, "dataset_name") <- ds
        res$matrix
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
        gs <- attr(mat, "gene_sources", exact = TRUE)
        if (!is.null(gs)) {
          gs_path <- file.path(out_dir, paste0(ds, GENE_SOURCES_SUFFIX))
          write_gene_sources_tsv(gs, gs_path)
        }
        n_force_only <- NA_integer_
        if (!is.null(gs)) {
          n_force_only <- sum(gs$source == "force_include")
        }
        dims_path <- file.path(out_dir, paste0(ds, "_pseudobulk_matrix.rds_dims.txt"))
        write_matrix_dims(mat, dims_path, n_force_only = n_force_only)
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

parse_cli_and_run <- function() {
  option_list <- list(
    make_option("--dataset_name", type = "character", default = NULL,
                help = "Run only this dataset from PSEUDOBULK_CONFIG [optional]"),
    make_option("--gene_list", type = "character", default = NULL,
                help = "Extra genes to force-include (.txt or .rds) [optional]"),
    make_option("--no_include_panel", action = "store_true", default = FALSE,
                help = "Do not union scDiffCom gene panel (default: panel is included)"),
    make_option("--min_detect_frac", type = "double", default = MIN_DETECT_FRAC_IN_MALIGNANT,
                help = "Detection fraction threshold for expressed genes [default %default]"),
    make_option("--out_dir", type = "character", default = OUT_DIR,
                help = "Output directory [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  include_panel <- !isTRUE(opt$no_include_panel)

  run_pseudobulk_pipeline(
    out_dir = path.expand(opt$out_dir),
    min_detect_frac = opt$min_detect_frac,
    gene_list_path = opt$gene_list,
    include_panel = include_panel,
    dataset_name_filter = opt$dataset_name
  )
}

if (sys.nframe() == 0L) {
  parse_cli_and_run()
}
