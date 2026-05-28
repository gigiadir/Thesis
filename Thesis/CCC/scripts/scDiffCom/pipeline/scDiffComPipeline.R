#!/usr/bin/env Rscript
options(future.globals.maxSize = 8000 * 1024^2)

# source("/gpfs0/bgu-ofircohen/group/groupRprofile") # Load shared functions and libraries

# === Diagnostics (remove after debugging) ===
cat("=== Job Info ===\n")
cat("Hostname:", system("hostname", intern = TRUE), "\n")
cat("Date:", format(Sys.time()), "\n")
cat("R version:", R.version.string, "\n")
cat("Rscript path:", system("which Rscript", intern = TRUE), "\n")
cat(".libPaths():\n")
cat(paste(" ", .libPaths(), collapse = "\n"), "\n")
cat("================\n")
# =============================================


library(optparse)
library(Seurat)
library(scDiffCom)
library(dplyr)
library(miceadds)


# --- 1. CLI Setup ---
option_list <- list(
  make_option("--gene", type="character", default=NULL, help="Gene symbol", metavar="character"),
  make_option("--dataset_path", type="character", default=NULL, help="Path to RAW Seurat RDS file", metavar="character"),
  make_option("--patient_summary_path", type="character", default=NULL, help="Path to the patient_summary RDS (output from the previous script)", metavar="character"),
  make_option("--column_id", type="character", default=NULL, help="The metadata column in Seurat that matches patient_id in summary (e.g. sample, orig.ident)", metavar="character"),
  make_option("--output_dir", type="character", default=NULL, help="Base output directory for scDiffCom results", metavar="character"),
  make_option("--cell_type_col", type="character", default="cell_type", help="Metadata column for cell type [default %default]", metavar="character"),
  make_option("--iterations", type="integer", default=10000, help="scDiffCom iterations [default %default]")
)
opt = optparse::parse_args(optparse::OptionParser(option_list=option_list))

# opt <- list(
#   gene         = "CHD4",                               # The gene you want to test
#   dataset_path = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Kurten_HNSC.RData",              # Local path to your RDS file
#   patient_summary_path = "/gpfs0/bgu-ofircohen/users/gigiadir/CCC-PreProcess/results-RankGenes/Kurten_HNSC/CHD4_Kurten_HNSC_grouped.rds",
#   column_id    = "Patient",                                 # Change to "orig.ident" or similar if needed
#   cell_type_col = "Cell_Type",
#   iterations = 10000
# )

# --- 2. Validation ---
if (is.null(opt$dataset_path) || !file.exists(opt$dataset_path)) stop("Error: dataset_path is invalid.")
if (is.null(opt$patient_summary_path) || !file.exists(opt$patient_summary_path)) stop("Error: patient_summary_path is invalid.")
if (is.null(opt$gene)) stop("Error: --gene parameter is required.")
if (is.null(opt$column_id)) stop("Error: --column_id is required to map patients to cells.")
if (is.null(opt$output_dir) || opt$output_dir == "") stop("Error: --output_dir is required and cannot be empty.")

GENE_EXP_COL <- paste0(opt$gene, "_EXP")

# --- 3. Helper Section ---
helpers <- list(
  
  create_seurat_obj = function(ds_path, ps_path, gene_col, patient_col_in_seurat) {
    message("Loading raw Seurat object...")
    ext <- tolower(tools::file_ext(opt$dataset_path))
    if (ext == "rds") {
      obj <- readRDS(opt$dataset_path)
    } else if (ext %in% c("rdata", "rda")) {
      miceadds::load.Rdata(opt$dataset_path, "obj")
    } else {
      stop(paste("Unsupported file format:", ext, ". Supported formats: .rds, .RData, .rda"))
    }
    
    message("Loading patient summary (preprocess output)...")
    ps <- readRDS(ps_path) # This is the summary_df from your previous script
    
    if (!patient_col_in_seurat %in% colnames(obj@meta.data)) {
      stop(sprintf("Error: Column '%s' not found in Seurat metadata.", patient_col_in_seurat))
    }
    
    message(sprintf("Merging %s classifications...", gene_col))
    
    # Map patient IDs to their HIGH/LOW/MID status
    # In your preprocess script, the summary has columns 'patient_id' and the gene_EXP column
    mapping_vec <- setNames(ps[[gene_col]], ps$patient_id)
    
    # Apply mapping to the Seurat metadata based on the user-provided column_id
    obj@meta.data[[gene_col]] <- mapping_vec[as.character(obj@meta.data[[patient_col_in_seurat]])]
    
    return(obj)
  },
  
  get_valid_cell_types = function(seurat_obj) {
    valid_cell_types <- seurat_obj@meta.data %>%
      filter(!!sym(opt$cell_type_col) != "" &
               !!sym(GENE_EXP_COL) %in% c("HIGH", "LOW")) %>%
      group_by(!!sym(opt$cell_type_col)) %>%
      filter(all(c("HIGH", "LOW") %in% !!sym(GENE_EXP_COL))) %>%
      pull(!!sym(opt$cell_type_col)) %>%
      unique()
    
    message(sprintf(
      "Found %d valid cell types containing both HIGH and LOW patients: %s", 
      length(valid_cell_types), 
      paste(valid_cell_types, collapse = ", ")
    ))
    return(valid_cell_types)
  },
  
  run_scdiffcom_for_seurat = function(seurat_obj) {
    valid_cts <- helpers$get_valid_cell_types(seurat_obj)
    if (length(valid_cts) == 0) stop("No valid cell types found for differential analysis.")
    
    # Filter to keep only HIGH/LOW and valid cell types
    seurat_sub <- subset(
      seurat_obj, 
      subset = !!sym(opt$cell_type_col) %in% valid_cts & !!sym(GENE_EXP_COL) %in% c("HIGH", "LOW")
    )
    
    scDiffCom_object <- run_interaction_analysis(
      seurat_object = seurat_sub,
      LRI_species = "human",
      seurat_celltype_id = opt$cell_type_col,
      seurat_condition_id = list(
        column_name = GENE_EXP_COL,
        cond1_name = "LOW",
        cond2_name = "HIGH"
      ),
      iterations = opt$iterations
    )
    return(scDiffCom_object)
  }
)

# --- 4. Main Execution ---
OUTPUT_DIR <- opt$output_dir
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
dataset_name <- tools::file_path_sans_ext(basename(opt$dataset_path))
dataset_output_dir <- file.path(OUTPUT_DIR, dataset_name)
if (!dir.exists(dataset_output_dir)) dir.create(dataset_output_dir, recursive = TRUE)

tryCatch({
  # 1. Create the merged object using the helper
  seurat_obj <- helpers$create_seurat_obj(
    ds_path = opt$dataset_path, 
    ps_path = opt$patient_summary_path, 
    gene_col = GENE_EXP_COL,
    patient_col_in_seurat = opt$column_id
  )
  
  # 2. Log basic stats
  n_high <- sum(seurat_obj@meta.data[[GENE_EXP_COL]] == "HIGH", na.rm = TRUE)
  n_low  <- sum(seurat_obj@meta.data[[GENE_EXP_COL]] == "LOW", na.rm = TRUE)
  message(sprintf("Total cells for analysis -> HIGH: %d, LOW: %d", n_high, n_low))
  
  # 3. Run scDiffCom
  message("Starting scDiffCom analysis...")
  result_obj <- helpers$run_scdiffcom_for_seurat(seurat_obj)
  
  # 4. Save
  out_path <- file.path(dataset_output_dir, paste0(opt$gene, "_", dataset_name, "_scDiffCom.rds"))
  
  saveRDS(result_obj, out_path)
  message("✓ Success! Saved to: ", out_path)
  
}, error = function(e) {
  message("✗ Error during processing: ", e$message)
})