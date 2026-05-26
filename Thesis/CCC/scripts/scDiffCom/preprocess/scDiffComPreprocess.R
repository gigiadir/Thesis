#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(miceadds)

# --- 1. CLI Setup ---
option_list <- list(
  make_option("--gene", type="character", default=NULL, help="Gene symbol", metavar="character"),
  make_option("--dataset_path", type="character", default=NULL, help="Path to Seurat object file (.rds or .RData/.rda)", metavar="character"),
  make_option("--column_id", type="character", default=NULL, help="Specific metadata column for patient IDs [optional]", metavar="character"),
  make_option("--cell_type_col", type="character", default="cell_type", help="Metadata column for cell type annotations [default %default]", metavar="character"),
  make_option("--low_q", type="double", default=1/3, help="Lower quantile [default %default]"),
  make_option("--high_q", type="double", default=2/3, help="Upper quantile [default %default]"),
  make_option("--output_path", type="character", default="/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/scDiffCom/Seurats", help="Path to put output files", metavar="character")
)
opt = optparse::parse_args(optparse::OptionParser(option_list=option_list))

# Uncomment to run interactively:
# opt <- list(
#   gene         = "AXL",
#   dataset_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/HNSCC.Atlas.RData",
#   column_id    = "Patient",
#   cell_type_col = "Cell_Type",
#   low_q        = 1/3,
#   high_q       = 2/3,
#   output_path  = "./outputs"
# )

# --- 2. Validation ---
if (is.null(opt$dataset_path) || !file.exists(opt$dataset_path)) {
  stop(paste("Error: Dataset path is invalid or file does not exist:", opt$dataset_path))
}

if (opt$low_q >= opt$high_q) {
  stop(sprintf("Error: --low_q (%f) must be strictly less than --high_q (%f)", opt$low_q, opt$high_q))
}

if (any(c(opt$low_q, opt$high_q) < 0) || any(c(opt$low_q, opt$high_q) > 1)) {
  stop("Error: Quantiles must be between 0 and 1.")
}

if (is.null(opt$gene)) {
  stop("Error: --gene parameter is required.")
}

# --- 3. Globals ---
SEURAT_OBJECTS_OUTPUT_DIR <- opt$output_path

MALIGNANT_CELLS <- c("Epithelial", "Malignant", "Tumor", "Cancer")
PATIENT_COLUMN_OPTIONS <- c("sample", "sampleID", "orig.ident")
MIN_MALIGNANT_CELLS <- 50

# --- 4. Helper Section ---
helpers <- list(
  add_malignant_gene_group = function(
    seurat_obj,
    gene = opt$gene,
    cell_type_col = opt$cell_type_col, 
    assay = "RNA",
    layer = "data",
    quantiles = c(opt$low_q, opt$high_q)
  ) {
    expr_matrix <- GetAssayData(seurat_obj, assay = assay, layer = layer)
    if (!gene %in% rownames(expr_matrix)) stop(paste("Gene", gene, "not found in object."))
    
    cell_metadata <- seurat_obj@meta.data
    cell_metadata$gene_expr_value <- expr_matrix[gene, ]
    
    # --- Patient Column Logic ---
    patient_col <- NULL
    if (!is.null(opt$column_id)) {
      if (opt$column_id %in% colnames(cell_metadata)) {
        patient_col <- opt$column_id
        message("Using user-specified patient column: ", patient_col)
      } else {
        stop(paste("Error: Specified --column_id", opt$column_id, "not found in metadata."))
      }
    } else {
      patient_col <- Filter(function(x) x %in% colnames(cell_metadata), PATIENT_COLUMN_OPTIONS)[1]
      if (is.null(patient_col)) stop("No valid patient column found in metadata. Please provide one via --column_id.")
      message("Auto-detected patient column: ", patient_col)
    }
    
    cell_metadata$patient_id_temp <- cell_metadata[[patient_col]]
    
    # --- Filtering and Aggregation ---
    # Validate the cell_type_col (either default or user-provided)
    if (!cell_type_col %in% colnames(cell_metadata)) {
      stop(paste("Error: Cell type column", cell_type_col, "not found in metadata."))
    }
    
    malignant_metadata <- cell_metadata[cell_metadata[[cell_type_col]] %in% MALIGNANT_CELLS, ]
    if (nrow(malignant_metadata) == 0) stop("No cells match MALIGNANT_CELLS labels.")

    malignant_cell_counts <- table(malignant_metadata$patient_id_temp)
    valid_patients <- names(malignant_cell_counts[malignant_cell_counts > MIN_MALIGNANT_CELLS])
    malignant_metadata <- malignant_metadata[malignant_metadata$patient_id_temp %in% valid_patients, ]
    if (nrow(malignant_metadata) == 0) stop(paste("No patients have more than", MIN_MALIGNANT_CELLS, "malignant cells."))
    message(sprintf("Retained %d patients with > %d malignant cells.", length(valid_patients), MIN_MALIGNANT_CELLS))

    patient_means <- aggregate(gene_expr_value ~ patient_id_temp, data = malignant_metadata, FUN = mean)
    colnames(patient_means) <- c("patient_id", "mean_expr")
    
    qs <- quantile(patient_means$mean_expr, probs = quantiles, na.rm = TRUE)
    
    group_col_name <- paste0(gene, "_EXP")
    patient_means[[group_col_name]] <- with(patient_means, ifelse(
      mean_expr <= qs[1], "LOW",
      ifelse(mean_expr >= qs[2], "HIGH", "MID")
    ))
    
    patient_group_map <- setNames(patient_means[[group_col_name]], patient_means$patient_id)
    seurat_obj@meta.data[[group_col_name]] <- patient_group_map[seurat_obj@meta.data[[patient_col]]]
    
    return(list(
      seurat_obj = seurat_obj,
      patient_summary = patient_means,
      group_col = group_col_name
    ))
  }
)

# --- 5. Main Execution ---
if (!dir.exists(SEURAT_OBJECTS_OUTPUT_DIR)) dir.create(SEURAT_OBJECTS_OUTPUT_DIR, recursive = TRUE)
ds_name = tools::file_path_sans_ext(basename(opt$dataset_path))
out_filename <- paste0(opt$gene, "_", ds_name, "_grouped.rds")
out_dir <- file.path(SEURAT_OBJECTS_OUTPUT_DIR, ds_name)
out_path <- file.path(out_dir, out_filename)

if (!dir.exists(out_dir)) {
  message("Creating directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}

if (file.exists(out_path)) {
  message("SUCCESS: Output file already exists. Skipping processing to save time.")
  message("Path: ", out_path)
  # Quit gracefully with exit code 0 (success)
  quit(save = "no", status = 0)
}

message("Loading Seurat object...")
ext <- tolower(tools::file_ext(opt$dataset_path))
if (ext == "rds") {
  obj <- readRDS(opt$dataset_path)
} else if (ext %in% c("rdata", "rda")) {
  miceadds::load.Rdata(opt$dataset_path, "obj")
} else {
  stop(paste("Unsupported file format:", ext, ". Supported formats: .rds, .RData, .rda"))
}

message(paste("Processing gene:", opt$gene))
results <- helpers$add_malignant_gene_group(seurat_obj = obj)

# --- Summary Printing ---
cat("\n--- Grouping Summary ---\n")
group_col <- results$group_col
summary_df <- results$patient_summary

counts <- table(summary_df[[group_col]])
print(counts)

for (grp in c("HIGH", "LOW")) {
  pts <- summary_df$patient_id[summary_df[[group_col]] == grp]
  cat(sprintf("\n%s Patients (%d):\n", grp, length(pts)))
  cat(paste(pts, collapse = ", "), "\n")
}
cat("------------------------\n\n")

# Save Output
saveRDS(summary_df, out_path)
message("Analysis complete. Saved to: ", out_path)