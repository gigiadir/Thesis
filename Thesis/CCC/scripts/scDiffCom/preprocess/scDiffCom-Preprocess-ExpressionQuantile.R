#!/usr/bin/env Rscript

# Expression quantile splits from pseudobulk (same rule as scDiffComPreprocess.R).
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-ExpressionQuantile.R --dataset_name Kurten_HNSC
#   3. Rscript compareExpressionSplitEquivalence.R --dataset_name Kurten_HNSC

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
              default = "~/CCC-PreProcess/results-ExpressionQuantile",
              help = "Base output directory [default %default]",
              metavar = "character"),
  make_option("--low_q", type = "double", default = 1/3,
              help = "Lower quantile for LOW group [default %default]"),
  make_option("--high_q", type = "double", default = 2/3,
              help = "Upper quantile for HIGH group [default %default]"),
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
if (opt$low_q >= opt$high_q) {
  stop("Error: --low_q must be strictly less than --high_q.")
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

message(
  "Building expression quantile splits (low_q=", opt$low_q,
  ", high_q=", opt$high_q, ") ..."
)
splits_df <- build_expression_quantile_splits_matrix(
  mat,
  low_q = opt$low_q,
  high_q = opt$high_q
)

splits_path <- file.path(out_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))
saveRDS(splits_df, splits_path)
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
