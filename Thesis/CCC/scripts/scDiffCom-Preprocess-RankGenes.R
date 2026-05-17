#!/usr/bin/env Rscript

# Rank genes per patient using pseudobulk matrices from createPseudobulkMatrix.R
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name Kurten_HNSC
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
source(file.path(script_dir, "rankGenesSplitUtils.R"))

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
              default = "~/CCC-PreProcess/results-RankGenes",
              help = "Base output directory [default %default]",
              metavar = "character"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Overwrite existing per-gene .rds files [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || opt$dataset_name == "") {
  stop("Error: --dataset_name is required.")
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

message("Ranking genes within each patient from pseudobulk expression (column-wise ranks).")
rank_mat <- compute_rank_matrix(mat)
stopifnot(identical(dim(rank_mat), dim(mat)))

patients <- colnames(mat)
genes <- rownames(mat)

message("Building patient tertile splits for all genes ...")
splits_df <- build_patient_splits_matrix(mat, rank_mat = rank_mat)

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

message("Writing per-gene grouped .rds files to: ", out_dir)
invisible(purrr::walk(seq_along(genes), write_one_gene))
message("Done.")
