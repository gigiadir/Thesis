#!/usr/bin/env Rscript

# Compare legacy scDiffComPreprocess.R quantile splits vs pseudobulk ExpressionQuantile splits.
#
# Run after:
#   scDiffCom-Preprocess-ExpressionQuantile.R
#   and legacy grouped .rds under CCC-PreProcess/results/
#
# Example:
#   Rscript compareExpressionSplitEquivalence.R --dataset_name Kurten_HNSC
#   Rscript compareExpressionSplitEquivalence.R --dataset_name all

suppressPackageStartupMessages({
  library(optparse)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(script_dir, "patientSplitUtils.R"))

HNSC_DATASETS <- c("Puram_HNSC", "Choi_HNSC", "Kurten_HNSC", "Bill_HNSC")

option_list <- list(
  make_option("--dataset_name", type = "character", default = NULL,
              help = "Dataset name (e.g. Kurten_HNSC) or 'all' for HNSC cohorts", metavar = "character"),
  make_option("--legacy_dir", type = "character",
              default = "~/CCC-PreProcess/results",
              help = "Legacy scDiffComPreprocess output [default %default]"),
  make_option("--new_dir", type = "character",
              default = "~/CCC-PreProcess/results-ExpressionQuantile",
              help = "ExpressionQuantile output [default %default]"),
  make_option("--pseudobulk_dir", type = "character",
              default = "~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix",
              help = "Pseudobulk matrix directory (for mean_expr diagnostics) [default %default]"),
  make_option("--output_dir", type = "character",
              default = "~/Thesis/CCC/outputs/QC/split_equivalence",
              help = "Directory for summary CSV files [default %default]"),
  make_option("--genes", type = "character", default = NULL,
              help = "Comma-separated genes or path to gene list; default = genes with legacy grouped .rds",
              metavar = "character"),
  make_option("--low_q", type = "double", default = 1/3,
              help = "Quantile used by preprocess (for reference) [default %default]"),
  make_option("--high_q", type = "double", default = 2/3,
              help = "Quantile used by preprocess (for reference) [default %default]"),
  make_option("--export_discordant", action = "store_true", default = FALSE,
              help = "Write per-gene discordant patient lists [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || opt$dataset_name == "") {
  stop("Error: --dataset_name is required (or use 'all').")
}

legacy_dir <- path.expand(opt$legacy_dir)
new_dir <- path.expand(opt$new_dir)
pseudobulk_dir <- path.expand(opt$pseudobulk_dir)
output_dir <- path.expand(opt$output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- if (tolower(opt$dataset_name) == "all") {
  HNSC_DATASETS
} else {
  opt$dataset_name
}

resolve_gene_list <- function(dataset_name, legacy_ds_dir) {
  if (!is.null(opt$genes) && nzchar(opt$genes)) {
    g <- opt$genes
    if (file.exists(path.expand(g))) {
      return(load_gene_list(g))
    }
    return(strsplit(g, ",", fixed = TRUE)[[1]])
  }
  files <- list.files(
    legacy_ds_dir,
    pattern = paste0("_", dataset_name, "_grouped\\.rds$"),
    full.names = FALSE
  )
  if (length(files) == 0L) {
    return(character())
  }
  sub(paste0("_", dataset_name, "_grouped\\.rds$"), "", files)
}

compare_one_gene <- function(gene, dataset_name, new_splits_mat, pb_mat = NULL) {
  legacy_path <- file.path(legacy_dir, dataset_name,
                           paste0(gene, "_", dataset_name, "_grouped.rds"))
  legacy <- load_one_grouped_gene(legacy_path)
  if (is.null(legacy)) {
    return(NULL)
  }

  if (!gene %in% rownames(new_splits_mat)) {
    return(NULL)
  }

  new_splits <- setNames(
    as.character(new_splits_mat[gene, ]),
    colnames(new_splits_mat)
  )

  patients_legacy <- legacy$patients
  patients_new <- names(new_splits)
  shared <- intersect(patients_legacy, patients_new)

  if (length(shared) == 0L) {
    return(data.frame(
      dataset = dataset_name,
      gene = gene,
      n_patients_legacy = length(patients_legacy),
      n_patients_new = length(patients_new),
      n_shared = 0L,
      n_agree = 0L,
      n_discordant = 0L,
      pct_agreement = NA_real_,
      max_abs_mean_expr_diff = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  leg_lab <- legacy$splits[shared]
  new_lab <- new_splits[shared]
  agree <- leg_lab == new_lab
  n_agree <- sum(agree, na.rm = TRUE)
  n_discordant <- sum(!agree, na.rm = TRUE)
  pct_agreement <- 100 * n_agree / length(shared)

  max_diff <- NA_real_
  if (!is.null(legacy$mean_expr)) {
    leg_expr <- legacy$mean_expr[shared]
    if (!is.null(pb_mat) && gene %in% rownames(pb_mat)) {
      pb_expr <- as.numeric(pb_mat[gene, shared])
      names(pb_expr) <- shared
      diffs <- abs(leg_expr - pb_expr)
      max_diff <- max(diffs[is.finite(diffs)], na.rm = TRUE)
      if (!is.finite(max_diff)) max_diff <- NA_real_
    } else {
      new_path <- file.path(new_dir, dataset_name,
                            paste0(gene, "_", dataset_name, "_grouped.rds"))
      new_g <- load_one_grouped_gene(new_path)
      if (!is.null(new_g) && !is.null(new_g$mean_expr)) {
        diffs <- abs(leg_expr - new_g$mean_expr[shared])
        max_diff <- max(diffs[is.finite(diffs)], na.rm = TRUE)
        if (!is.finite(max_diff)) max_diff <- NA_real_
      }
    }
  }

  if (isTRUE(opt$export_discordant) && n_discordant > 0L) {
    disc <- shared[!agree]
    disc_file <- file.path(
      output_dir,
      paste0(dataset_name, "_", gene, "_discordant_patients.txt")
    )
    writeLines(disc, disc_file)
  }

  data.frame(
    dataset = dataset_name,
    gene = gene,
    n_patients_legacy = length(patients_legacy),
    n_patients_new = length(patients_new),
    n_shared = length(shared),
    n_agree = n_agree,
    n_discordant = n_discordant,
    pct_agreement = pct_agreement,
    max_abs_mean_expr_diff = max_diff,
    stringsAsFactors = FALSE
  )
}

compare_dataset <- function(dataset_name) {
  legacy_ds_dir <- file.path(legacy_dir, dataset_name)
  if (!dir.exists(legacy_ds_dir)) {
    warning("Legacy directory not found: ", legacy_ds_dir, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }

  new_splits_mat <- load_splits_from_grouped_dir(new_dir, dataset_name)
  if (is.null(new_splits_mat)) {
    warning(
      "New ExpressionQuantile splits not found for ", dataset_name,
      ". Run scDiffCom-Preprocess-ExpressionQuantile.R first.",
      call. = FALSE, immediate. = TRUE
    )
    return(NULL)
  }

  pb_mat <- NULL
  pb_path <- tryCatch(
    resolve_pseudobulk_path(dataset_name, pseudobulk_dir),
    error = function(e) NULL
  )
  if (!is.null(pb_path) && file.exists(pb_path)) {
    pb_mat <- load_gene_patient_matrix(pb_path)
  }

  genes <- resolve_gene_list(dataset_name, legacy_ds_dir)
  genes <- intersect(genes, rownames(new_splits_mat))
  if (length(genes) == 0L) {
    warning("No genes to compare for ", dataset_name, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }

  rows <- lapply(genes, function(g) {
    compare_one_gene(g, dataset_name, new_splits_mat, pb_mat = pb_mat)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(NULL)
  }
  do.call(rbind, rows)
}

all_rows <- lapply(datasets, compare_dataset)
all_rows <- all_rows[!vapply(all_rows, is.null, logical(1))]

if (length(all_rows) == 0L) {
  stop("No comparison results produced. Check legacy and ExpressionQuantile output paths.")
}

summary_df <- do.call(rbind, all_rows)
summary_path <- file.path(output_dir, "expression_quantile_equivalence_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)
message("Wrote ", summary_path)

by_ds <- do.call(rbind, lapply(split(summary_df, summary_df$dataset), function(df) {
  data.frame(
    dataset = df$dataset[[1]],
    n_genes = nrow(df),
    mean_pct_agreement = mean(df$pct_agreement, na.rm = TRUE),
    median_pct_agreement = stats::median(df$pct_agreement, na.rm = TRUE),
    n_genes_100pct = sum(df$pct_agreement >= 100 - 1e-9, na.rm = TRUE),
    total_discordant = sum(df$n_discordant),
    stringsAsFactors = FALSE
  )
}))
by_ds_path <- file.path(output_dir, "expression_quantile_equivalence_by_dataset.csv")
utils::write.csv(by_ds, by_ds_path, row.names = FALSE)
message("Wrote ", by_ds_path)

for (i in seq_len(nrow(by_ds))) {
  message(
    by_ds$dataset[[i]], ": ",
    by_ds$n_genes_100pct[[i]], "/", by_ds$n_genes[[i]],
    " genes with 100% agreement on shared patients",
    " (mean agreement ", sprintf("%.1f", by_ds$mean_pct_agreement[[i]]), "%)"
  )
}

message("Done.")
