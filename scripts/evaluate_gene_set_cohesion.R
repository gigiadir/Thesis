#!/usr/bin/env Rscript

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

source(file.path(script_dir, "lib", "cohesion_metrics.R"))

option_list <- list(
  make_option("--reference_rds", type = "character", default = NULL,
              help = "Path to .rds list: set_name -> gene vector", metavar = "character"),
  make_option("--method_name", action = "append", type = "character", default = NULL,
              help = "Method label (repeat per method): --method_name RankGenes", metavar = "character"),
  make_option("--method_input", action = "append", type = "character", default = NULL,
              help = "Path to method matrix .rds/.csv/.tsv (repeat per method)", metavar = "character"),
  make_option("--method_mode", action = "append", type = "character", default = NULL,
              help = "feature_matrix or distance_matrix (repeat per method)", metavar = "character"),
  make_option("--method_similarity", action = "append", type = "character", default = NULL,
              help = "TRUE/FALSE for distance_matrix interpretation (repeat per method)", metavar = "character"),
  make_option("--feature_distance", type = "character", default = "correlation",
              help = "Distance for feature matrices: euclidean|correlation [default %default]", metavar = "character"),
  make_option("--min_set_size", type = "integer", default = 3L,
              help = "Minimum genes per set after intersection [default %default]", metavar = "integer"),
  make_option("--max_positive_pairs", type = "integer", default = 50000L,
              help = "Global cap on positive pairs [default %default]", metavar = "integer"),
  make_option("--max_positive_pairs_per_set", type = "integer", default = 10000L,
              help = "Per-set cap on positive pairs [default %default]", metavar = "integer"),
  make_option("--negative_ratio", type = "double", default = 1,
              help = "Negative-pair count = positive_count * ratio [default %default]", metavar = "double"),
  make_option("--bootstrap_iters", type = "integer", default = 300L,
              help = "Bootstrap iterations for CI [default %default]", metavar = "integer"),
  make_option("--permutation_iters", type = "integer", default = 0L,
              help = "Permutation iterations for p-value (0=skip) [default %default]", metavar = "integer"),
  make_option("--seed", type = "integer", default = 1L,
              help = "Seed for reproducibility [default %default]", metavar = "integer"),
  make_option("--output_dir", type = "character", default = file.path(script_dir, "outputs", "gene_set_cohesion"),
              help = "Output directory [default %default]", metavar = "character"),
  make_option("--run_sanity_check", action = "store_true", default = FALSE,
              help = "Run synthetic sanity scenarios and write sanity outputs [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$reference_rds) || !nzchar(opt$reference_rds)) {
  stop("--reference_rds is required.")
}
if (is.null(opt$method_name) || is.null(opt$method_input) || is.null(opt$method_mode)) {
  stop("Provide repeated --method_name, --method_input, --method_mode values.")
}
if (!(length(opt$method_name) == length(opt$method_input) &&
      length(opt$method_input) == length(opt$method_mode))) {
  stop("--method_name, --method_input, --method_mode must have equal lengths.")
}

if (is.null(opt$method_similarity)) {
  opt$method_similarity <- rep("FALSE", length(opt$method_name))
}
if (length(opt$method_similarity) == 1L && length(opt$method_name) > 1L) {
  opt$method_similarity <- rep(opt$method_similarity, length(opt$method_name))
}
if (length(opt$method_similarity) != length(opt$method_name)) {
  stop("--method_similarity must be length 1 or match number of methods.")
}

read_input_matrix <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    return(readRDS(path))
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    m <- utils::read.table(path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE)
    return(as.matrix(m))
  }
  stop("Unsupported input extension for: ", path, ". Use .rds/.csv/.tsv/.txt")
}

rank_summary <- function(df) {
  out <- df
  out <- out[order(out$separation_ratio, -out$auc, out$cliffs_delta), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  out
}

reference_sets_all <- readRDS(path.expand(opt$reference_rds))
dir.create(path.expand(opt$output_dir), recursive = TRUE, showWarnings = FALSE)

method_rows <- list()
set_level_rows <- list()
coverage_rows <- list()
excluded_rows <- list()
pair_tables <- list()

for (i in seq_along(opt$method_name)) {
  method_name <- opt$method_name[[i]]
  message("Evaluating method: ", method_name)

  raw_matrix <- read_input_matrix(path.expand(opt$method_input[[i]]))
  mode_i <- tolower(opt$method_mode[[i]])
  sim_i <- as.logical(toupper(opt$method_similarity[[i]]))
  if (!mode_i %in% c("feature_matrix", "distance_matrix")) {
    stop("method_mode must be feature_matrix or distance_matrix; got ", mode_i)
  }

  distance_matrix <- prepare_distance_matrix(
    raw_matrix,
    input_mode = mode_i,
    similarity_input = isTRUE(sim_i),
    feature_distance = opt$feature_distance
  )

  prep <- prepare_reference_sets(
    reference_sets = reference_sets_all,
    genes_available = rownames(distance_matrix),
    min_set_size = opt$min_set_size
  )
  reference_sets <- prep$sets
  if (length(reference_sets) == 0L) {
    stop("No reference sets left after filtering for method ", method_name)
  }

  gene_pool <- unique(unlist(reference_sets, use.names = FALSE))
  pair_df <- build_pair_table(
    reference_sets = reference_sets,
    all_genes = gene_pool,
    max_positive_pairs = opt$max_positive_pairs,
    max_positive_pairs_per_set = opt$max_positive_pairs_per_set,
    negative_ratio = opt$negative_ratio,
    seed = opt$seed + i
  )

  eval_res <- compute_method_metrics(
    pair_df = pair_df,
    distance_matrix = distance_matrix,
    bootstrap_iters = opt$bootstrap_iters,
    permutation_iters = opt$permutation_iters,
    seed = opt$seed + i
  )

  summary_row <- eval_res$summary
  summary_row$method <- method_name
  summary_row$input_path <- opt$method_input[[i]]
  summary_row$input_mode <- mode_i
  summary_row$input_similarity <- isTRUE(sim_i)
  summary_row$n_reference_sets <- length(reference_sets)
  summary_row$n_reference_genes <- length(gene_pool)
  method_rows[[i]] <- summary_row

  set_diag <- compute_set_level_diagnostics(reference_sets, distance_matrix)
  set_diag$method <- rep(method_name, nrow(set_diag))
  set_level_rows[[i]] <- set_diag

  cov_df <- prep$coverage
  cov_df$method <- rep(method_name, nrow(cov_df))
  coverage_rows[[i]] <- cov_df

  drop_df <- prep$dropped
  drop_df$method <- rep(method_name, nrow(drop_df))
  excluded_rows[[i]] <- drop_df

  pair_tables[[method_name]] <- eval_res$pair_distances
}

summary_df <- do.call(rbind, method_rows)
summary_ranked <- rank_summary(summary_df)
set_level_df <- do.call(rbind, set_level_rows)
coverage_df <- do.call(rbind, coverage_rows)
excluded_df <- do.call(rbind, excluded_rows)

summary_csv <- file.path(path.expand(opt$output_dir), "method_summary.csv")
summary_rds <- file.path(path.expand(opt$output_dir), "method_summary.rds")
set_level_csv <- file.path(path.expand(opt$output_dir), "set_level_diagnostics.csv")
coverage_csv <- file.path(path.expand(opt$output_dir), "set_coverage.csv")
excluded_csv <- file.path(path.expand(opt$output_dir), "excluded_sets.csv")
pairs_rds <- file.path(path.expand(opt$output_dir), "pair_distances_by_method.rds")

utils::write.csv(summary_ranked, summary_csv, row.names = FALSE)
saveRDS(summary_ranked, summary_rds)
utils::write.csv(set_level_df, set_level_csv, row.names = FALSE)
utils::write.csv(coverage_df, coverage_csv, row.names = FALSE)
utils::write.csv(excluded_df, excluded_csv, row.names = FALSE)
saveRDS(pair_tables, pairs_rds)

if (isTRUE(opt$run_sanity_check)) {
  message("Running synthetic sanity checks...")
  syn_sets <- list(setA = paste0("g", 1:5), setB = paste0("g", 6:10))
  genes <- unlist(syn_sets, use.names = FALSE)

  clustered <- matrix(rnorm(10 * 5, sd = 0.2), nrow = 10, ncol = 5, dimnames = list(genes, paste0("f", 1:5)))
  clustered[syn_sets$setA, ] <- clustered[syn_sets$setA, ] + 3
  clustered[syn_sets$setB, ] <- clustered[syn_sets$setB, ] - 3
  random_mat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5, dimnames = list(genes, paste0("f", 1:5)))

  eval_one_sanity <- function(feat_mat, scenario_name, seed) {
    dmat <- prepare_distance_matrix(feat_mat, input_mode = "feature_matrix", feature_distance = "euclidean")
    pairs <- build_pair_table(
      reference_sets = syn_sets,
      all_genes = genes,
      max_positive_pairs = NULL,
      max_positive_pairs_per_set = NULL,
      negative_ratio = 1,
      seed = seed
    )
    res <- compute_method_metrics(
      pair_df = pairs,
      distance_matrix = dmat,
      bootstrap_iters = 200L,
      permutation_iters = 200L,
      seed = seed
    )
    out <- res$summary
    out$scenario <- scenario_name
    out
  }

  sanity_df <- rbind(
    eval_one_sanity(clustered, "clustered_same_set", opt$seed + 1000L),
    eval_one_sanity(random_mat, "random_features", opt$seed + 2000L)
  )
  utils::write.csv(sanity_df, file.path(path.expand(opt$output_dir), "sanity_check_summary.csv"), row.names = FALSE)
}

message("Done. Wrote:")
message("  - ", summary_csv)
message("  - ", set_level_csv)
message("  - ", coverage_csv)
message("  - ", excluded_csv)
message("  - ", pairs_rds)
