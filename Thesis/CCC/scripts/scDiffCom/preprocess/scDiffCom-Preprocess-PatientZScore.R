#!/usr/bin/env Rscript

# Patient-level z-score tertiles from pseudobulk matrices (createPseudobulkMatrix.R).
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R --dataset_name Choi_HNSC
#   2. Rscript scDiffCom-Preprocess-PatientZScore.R --dataset_name Choi_HNSC

suppressPackageStartupMessages({
  library(optparse)
  library(purrr)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(script_dir, "patientSplitUtils.R"))

option_list <- list(
  make_option("--dataset_name", type = "character", default = NULL,
              help = "Dataset directory name (e.g. Kurten_HNSC)", metavar = "character"),
  make_option("--pseudobulk_dir", type = "character",
              default = "~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix",
              help = paste0(
                "Directory with pseudobulk .rds files (*",
                PSEUDOBULK_MATRIX_SUFFIX,
                ") [default %default]"
              ),
              metavar = "character"),
  make_option("--output_base", type = "character",
              default = "~/CCC-PreProcess/results-Patient-ZScore",
              help = "Base output directory [default %default]",
              metavar = "character"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing per-gene .rds files [default %default]"),
  make_option("--gene_list", type = "character", default = NULL,
              help = "Only write grouped .rds for these genes (.txt/.rds) [optional]",
              metavar = "character")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || opt$dataset_name == "") {
  stop("Error: --dataset_name is required.")
}

pseudobulk_dir <- path.expand(opt$pseudobulk_dir)
output_base <- path.expand(opt$output_base)
in_path <- resolve_pseudobulk_path(opt$dataset_name, pseudobulk_dir)

out_dir <- file.path(output_base, opt$dataset_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading gene×patient pseudobulk matrix: ", in_path)
mat <- load_gene_patient_matrix(in_path)

message("Computing row-wise Z-scores (vectorized). Genes=", nrow(mat), " Patients=", ncol(mat))

z_mat <- NULL
if (requireNamespace("matrixStats", quietly = TRUE)) {
  rm <- matrixStats::rowMeans2(mat, na.rm = TRUE)
  rsd <- matrixStats::rowSds(mat, na.rm = TRUE)
  rsd[is.na(rsd) | rsd == 0] <- NA_real_
  z_mat <- sweep(mat, 1, rm, FUN = "-")
  z_mat <- sweep(z_mat, 1, rsd, FUN = "/")
} else {
  message("Note: package 'matrixStats' not found; falling back to t(scale(t(mat))).")
  z_mat <- t(scale(t(mat)))
}

stopifnot(identical(dim(z_mat), dim(mat)))

patients <- colnames(mat)
genes <- rownames(mat)

splits_mat <- matrix(NA_character_, nrow = length(genes), ncol = length(patients),
                     dimnames = list(genes, patients))
for (i in seq_along(genes)) {
  splits_mat[i, ] <- assign_tertiles(z_mat[i, ])
}

splits_path <- file.path(out_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))
saveRDS(splits_mat, splits_path)
message("Saved all patient splits: ", splits_path)

write_one_gene <- function(i) {
  gene <- genes[[i]]
  out_path <- file.path(out_dir, paste0(gene, "_", opt$dataset_name, "_grouped.rds"))
  if (!opt$overwrite && file.exists(out_path)) return(invisible(NULL))

  df <- data.frame(
    patient_id = patients,
    mean_expr = as.numeric(mat[i, ]),
    stringsAsFactors = FALSE
  )
  exp_col <- paste0(gene, "_exp")
  df[[exp_col]] <- splits_mat[i, ]

  df <- df[, c("patient_id", "mean_expr", exp_col)]
  saveRDS(df, out_path)
  invisible(NULL)
}

gene_idx <- gene_indices_to_write(genes, opt$gene_list)
if (!is.null(opt$gene_list) && nzchar(opt$gene_list)) {
  message("Writing grouped .rds for ", length(gene_idx), " gene(s) from --gene_list")
} else {
  message("Writing per-gene grouped .rds files to: ", out_dir)
}
invisible(purrr::walk(gene_idx, write_one_gene))
message("Done.")
