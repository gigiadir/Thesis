#!/usr/bin/env Rscript

# Residual-based patient splits from pseudobulk matrices (PC-regressed expression).
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-Residual.R --dataset_name Kurten_HNSC --n_pc 2
#
# PCA uses patients with complete pseudobulk across all genes (prcomp on t(mat),
# scale.=TRUE). Per-gene residuals regress expression on PC1..PCn; tertiles are
# rank-based (top third = HIGH, bottom third = LOW).
#
# R libs (if optparse missing from plain Rscript):
#   export R_LIBS_SITE=/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0
#   or source groupRprofile in ~/.Rprofile

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
  make_option("--input_dir", type = "character",
              default = "~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix",
              help = paste0(
                "Directory with pseudobulk .rds files (*",
                PSEUDOBULK_MATRIX_SUFFIX,
                ") [default %default]"
              ),
              metavar = "character"),
  make_option("--output_base", type = "character",
              default = "~/CCC-PreProcess/results-Residual",
              help = "Base output directory [default %default]",
              metavar = "character"),
  make_option("--n_pc", type = "integer", default = 2L,
              help = "Number of patient PCs to regress out (1 or 2) [default %default]"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing per-gene .rds files [default %default]"),
  make_option("--gene_list", type = "character", default = NULL,
              help = "Only write grouped .rds for these genes (.txt/.rds); splits still use full matrix [optional]",
              metavar = "character")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || opt$dataset_name == "") {
  stop("Error: --dataset_name is required.")
}

n_pc <- as.integer(opt$n_pc)
if (!n_pc %in% c(1L, 2L)) {
  stop("Error: --n_pc must be 1 or 2. Got: ", n_pc)
}

input_dir <- path.expand(opt$input_dir)
output_base <- path.expand(opt$output_base)
in_path <- resolve_pseudobulk_path(opt$dataset_name, input_dir)

if (!grepl(PSEUDOBULK_MATRIX_SUFFIX, in_path, fixed = TRUE)) {
  warning(
    "Using legacy gene×patient matrix (not pseudobulk): ", in_path,
    call. = FALSE, immediate. = TRUE
  )
}

out_dir <- file.path(output_base, opt$dataset_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading gene×patient pseudobulk matrix: ", in_path)
mat <- load_gene_patient_matrix(in_path)

patients <- colnames(mat)
genes <- rownames(mat)

message("Computing patient PCs (n_pc=", n_pc, ") and per-gene residuals ...")
res <- build_residual_splits_matrix(mat, n_pc = n_pc)
splits_df <- res$splits
resid_mat <- res$residuals

splits_path <- file.path(out_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))
saveRDS(splits_df, splits_path)
message("Saved all patient splits: ", splits_path)

resid_path <- file.path(out_dir, paste0(opt$dataset_name, RESIDUAL_MATRIX_SUFFIX))
saveRDS(resid_mat, resid_path)
message("Saved residual matrix: ", resid_path)

write_one_gene <- function(i) {
  gene <- genes[[i]]
  out_path <- file.path(out_dir, paste0(gene, "_", opt$dataset_name, "_grouped.rds"))
  if (!opt$overwrite && file.exists(out_path)) return(invisible(NULL))

  df <- data.frame(
    patient_id = patients,
    mean_expr = as.numeric(mat[i, ]),
    stringsAsFactors = FALSE
  )
  exp_col <- paste0(gene, "_EXP")
  df[[exp_col]] <- splits_df[i, ]

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
