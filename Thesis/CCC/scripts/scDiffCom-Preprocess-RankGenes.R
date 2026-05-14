#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
})

option_list <- list(
  make_option("--dataset_name", type = "character", default = NULL,
              help = "Dataset directory name (e.g. Kurten_HNSC)", metavar = "character"),
  make_option("--input_dir", type = "character",
              default = "~/Thesis/CCC/outputs/RData_objects/Datasets_metrics",
              help = "Directory containing the compiled gene×patient matrix .rds [default %default]",
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

in_candidates <- c(
  file.path(input_dir, paste0(opt$dataset_name, "_gene_by_patient_means_exp.rds")),
  file.path(input_dir, paste0(opt$dataset_name, "_gene_by_patient_mean_expr.rds"))
)
in_path <- in_candidates[file.exists(in_candidates)][1]

if (is.na(in_path) || is.null(in_path)) {
  stop(
    "Error: Could not find input .rds for dataset '", opt$dataset_name, "'. Tried:\n",
    paste(" -", in_candidates, collapse = "\n")
  )
}

out_dir <- file.path(output_base, opt$dataset_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading gene×patient matrix: ", in_path)
x <- readRDS(in_path)

mat <- if (is.data.frame(x)) {
  m <- as.matrix(x)
  if (!is.null(rownames(x))) rownames(m) <- rownames(x)
  m
} else if (is.matrix(x)) {
  x
} else {
  stop("Input must be a matrix or data.frame. Got: ", paste(class(x), collapse = ", "))
}

if (is.null(rownames(mat)) || anyNA(rownames(mat)) || any(rownames(mat) == "")) {
  stop("Input matrix must have non-empty rownames (gene symbols).")
}
if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(colnames(mat) == "")) {
  stop("Input matrix must have non-empty colnames (patient_ids).")
}

storage.mode(mat) <- "double"

message("Ranking genes within each patient (column-wise ranks).")
rank_mat <- NULL
if (requireNamespace("matrixStats", quietly = TRUE) &&
    exists("colRanks", where = asNamespace("matrixStats"), inherits = FALSE)) {
  # Highest expression gets the highest rank.
  rank_mat <- matrixStats::colRanks(
    mat,
    ties.method = "average",
    preserveShape = TRUE,
    na.last = "keep"
  )
} else {
  message("Note: matrixStats::colRanks not available; falling back to apply(..., 2, rank).")
  rank_mat <- apply(
    mat,
    2,
    function(v) rank(v, ties.method = "average", na.last = "keep")
  )
}

stopifnot(identical(dim(rank_mat), dim(mat)))

patients <- colnames(mat)
genes <- rownames(mat)

assign_tertiles <- function(scores) {
  ok <- !is.na(scores)
  n_ok <- sum(ok)
  out <- rep(NA_character_, length(scores))

  if (n_ok == 0) return(out)
  if (length(unique(scores[ok])) == 1L) {
    out[ok] <- "MID"
    return(out)
  }

  r <- rank(scores[ok], ties.method = "average")
  out_ok <- ifelse(
    r <= n_ok / 3, "LOW",
    ifelse(r > 2 * n_ok / 3, "HIGH", "MID")
  )
  out[ok] <- out_ok
  out
}

write_one_gene <- function(i) {
  gene <- genes[[i]]
  out_path <- file.path(out_dir, paste0(gene, "_", opt$dataset_name, "_grouped.rds"))
  if (!opt$overwrite && file.exists(out_path)) return(invisible(NULL))

  mean_expr <- mat[i, ]
  groups <- assign_tertiles(rank_mat[i, ])

  df <- data.frame(
    patient_id = patients,
    mean_expr = as.numeric(mean_expr),
    stringsAsFactors = FALSE
  )
  exp_col <- paste0(gene, "_EXP")
  df[[exp_col]] <- groups

  df <- df[, c("patient_id", "mean_expr", exp_col)]
  saveRDS(df, out_path)
  invisible(NULL)
}

message("Writing per-gene grouped .rds files to: ", out_dir)
invisible(purrr::walk(seq_along(genes), write_one_gene))
message("Done.")

