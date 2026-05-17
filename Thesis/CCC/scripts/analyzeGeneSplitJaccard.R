#!/usr/bin/env Rscript

# Compare patient LOW/MID/HIGH splits across genes via Jaccard similarity.
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name Kurten_HNSC
#   3. Rscript analyzeGeneSplitJaccard.R --dataset_name Kurten_HNSC --mode panel
#
# Requires: optparse, pheatmap
# R libraries (if packages missing from plain Rscript):
#   export R_LIBS_SITE=/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0
#   or uncomment source("/gpfs0/bgu-ofircohen/group/groupRprofile") in ~/.Rprofile
#
# Modes:
#   panel — first 100 genes from panel list: Jaccard matrix + readable heatmap
#   all   — all genes: duplicate-split clusters; full matrix only if <= max_genes_full_matrix

PANEL_MAX_GENES <- 123L

suppressPackageStartupMessages({
  library(optparse)
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with install.packages(\"pheatmap\").")
  }
  library(pheatmap)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(script_dir, "rankGenesSplitUtils.R"))

DEFAULT_GENE_PANEL_RDS <- path.expand(
  "~/Thesis/CCC/outputs/RData_objects/GenesLists/Complexes.Oncogenes.OncoKB.Cosmic.NCG.rds"
)

option_list <- list(
  make_option("--dataset_name", type = "character", default = NULL,
              help = "Dataset name (e.g. Kurten_HNSC)", metavar = "character"),
  make_option("--mode", type = "character", default = "panel",
              help = "panel or all [default %default]", metavar = "character"),
  make_option("--rankgenes_dir", type = "character",
              default = "~/CCC-PreProcess/results-RankGenes",
              help = "RankGenes output base directory [default %default]",
              metavar = "character"),
  make_option("--pseudobulk_dir", type = "character",
              default = "~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix",
              help = "Pseudobulk matrix directory [default %default]",
              metavar = "character"),
  make_option("--gene_list", type = "character", default = NULL,
              help = "Optional .rds or .txt gene list (default: scDiffCom panel)",
              metavar = "character"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory [default ~/Thesis/CCC/outputs/split_similarity/{dataset}]",
              metavar = "character"),
  make_option("--no_cluster", action = "store_true", default = FALSE,
              help = "Skip hierarchical clustering and heatmap reordering"),
  make_option("--jaccard_threshold", type = "double", default = 1.0,
              help = "Threshold for duplicate/near-duplicate groups [default %default]"),
  make_option("--max_genes_full_matrix", type = "integer", default = 2500L,
              help = "In 'all' mode, skip full pairwise matrix above this G [default %default]"),
  make_option("--top_pairs", type = "integer", default = 100L,
              help = "Number of top similar gene pairs to export [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || !nzchar(opt$dataset_name)) {
  stop("--dataset_name is required.")
}

mode <- tolower(opt$mode)
if (!mode %in% c("panel", "all")) {
  stop("--mode must be 'panel' or 'all'.")
}

rankgenes_dir <- path.expand(opt$rankgenes_dir)
pseudobulk_dir <- path.expand(opt$pseudobulk_dir)
out_dir <- if (!is.null(opt$output_dir) && nzchar(opt$output_dir)) {
  path.expand(opt$output_dir)
} else {
  path.expand(file.path("~/Thesis/CCC/outputs/split_similarity", opt$dataset_name))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_panel_genes <- function() {
  if (file.exists(DEFAULT_GENE_PANEL_RDS)) {
    return(unique(as.character(readRDS(DEFAULT_GENE_PANEL_RDS))))
  }
  SCDIFFCOM_GENE_PANEL
}

load_splits_for_analysis <- function() {
  ds_dir <- file.path(rankgenes_dir, opt$dataset_name)
  all_splits_path <- file.path(ds_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))

  panel_genes <- if (!is.null(opt$gene_list) && nzchar(opt$gene_list)) {
    load_gene_list(opt$gene_list)
  } else {
    load_panel_genes()
  }

  if (mode == "panel" && length(panel_genes) > PANEL_MAX_GENES) {
    message(
      "Panel mode: using first ", PANEL_MAX_GENES, " of ",
      length(panel_genes), " genes (list order preserved)."
    )
    panel_genes <- panel_genes[seq_len(PANEL_MAX_GENES)]
  }

  if (file.exists(all_splits_path)) {
    message("Loading cached all-patient splits: ", all_splits_path)
    all_splits <- readRDS(all_splits_path)
    if (mode == "panel") {
      keep <- intersect(panel_genes, rownames(all_splits))
      if (length(keep) == 0L) {
        stop("No panel genes found in cached splits at ", all_splits_path)
      }
      missing <- setdiff(panel_genes, keep)
      if (length(missing) > 0L) {
        warning(
          "Panel genes missing from cache: ", length(missing),
          call. = FALSE, immediate. = TRUE
        )
      }
      return(all_splits[intersect(panel_genes, keep), , drop = FALSE])
    }
    return(all_splits)
  }

  if (mode == "panel") {
    message("Loading splits for ", length(panel_genes), " panel genes from RankGenes ...")
    splits <- load_splits_from_rankgenes_dir(
      rankgenes_dir, opt$dataset_name, genes = panel_genes
    )
    missing <- setdiff(panel_genes, rownames(splits))
    if (length(missing) > 0L) {
      warning(
        "Missing grouped .rds for ", length(missing), " gene(s); continuing with ",
        nrow(splits), " genes.", call. = FALSE, immediate. = TRUE
      )
    }
    return(splits[intersect(panel_genes, rownames(splits)), , drop = FALSE])
  }

  message("Recomputing splits from pseudobulk (all genes) ...")
  in_path <- resolve_pseudobulk_path(opt$dataset_name, pseudobulk_dir)
  mat <- load_gene_patient_matrix(in_path)
  rank_mat <- compute_rank_matrix(mat)
  build_patient_splits_matrix(mat, rank_mat = rank_mat)
}

splits <- load_splits_for_analysis()
message("Split matrix: ", nrow(splits), " genes x ", ncol(splits), " patients")

L <- splits_to_integer_matrix(splits)
G <- nrow(L)

# Exact duplicate groups (hash-based, fast for all mode)
keys <- apply(L, 1, paste, collapse = ",")
hash_groups <- split(rownames(splits), keys)
exact_dup <- hash_groups[lengths(hash_groups) > 1L]
message("Exact duplicate split groups (J=1): ", length(exact_dup))

dup_groups <- if (opt$jaccard_threshold >= 1) {
  exact_dup
} else {
  find_duplicate_split_groups(splits, threshold = opt$jaccard_threshold)
}

dup_tsv <- duplicate_groups_to_tsv(dup_groups)
write.table(
  dup_tsv,
  file = file.path(out_dir, "split_duplicate_groups.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

compute_full_matrix <- (mode == "panel") || (G <= opt$max_genes_full_matrix)
J <- NULL
dist_mat <- NULL

if (compute_full_matrix) {
  message("Computing pairwise Jaccard similarity (G=", G, ") ...")
  J <- compute_jaccard_similarity_matrix(L)
  dist_mat <- 1 - J
  diag(dist_mat) <- 0

  saveRDS(J, file.path(out_dir, "jaccard_similarity.rds"))
  saveRDS(dist_mat, file.path(out_dir, "jaccard_distance.rds"))

  pairs_df <- top_similar_pairs(J, n = opt$top_pairs)
  write.csv(
    pairs_df,
    file.path(out_dir, "split_agreement_summary.csv"),
    row.names = FALSE
  )

  do_cluster <- !isTRUE(opt$no_cluster) && G >= 2L
  hc <- NULL
  if (do_cluster) {
    d <- as.dist(dist_mat)
    d[is.na(d)] <- 1
    hc <- hclust(d, method = "average")
    saveRDS(hc, file.path(out_dir, "gene_clusters.rds"))
  }

  D_plot <- dist_mat
  diag(D_plot) <- NA_real_

  heatmap_path <- file.path(out_dir, "jaccard_heatmap.png")
  label_fs <- max(4, 200 / nrow(D_plot))

  png(heatmap_path, width = 4000, height = 4000, res = 300)
  pheatmap(
    D_plot,
    color           = colorRampPalette(c("#d73027", "white", "#4575b4"))(100),
    breaks = seq(0, 1, length.out = 101),
    na_col = "#E8E8E8",
    cluster_rows = if (isTRUE(opt$no_cluster)) FALSE else (if (is.null(hc)) TRUE else hc),
    cluster_cols = if (isTRUE(opt$no_cluster)) FALSE else (if (is.null(hc)) TRUE else hc),
    fontsize_row = label_fs,
    fontsize_col = label_fs,
    angle_col = 45,
    main = paste0(
      opt$dataset_name, " – Gene-Gene Jaccard Split Distance (n=", G, ")"
    )
  )
  dev.off()
  message("Wrote heatmap (pheatmap, distance): ", heatmap_path)
} else {
  message(
    "Skipping full Jaccard matrix (G=", G, " > max_genes_full_matrix=",
    opt$max_genes_full_matrix, "). See split_duplicate_groups.tsv"
  )
  # Top pairs only among representatives of exact duplicate groups
  reps <- vapply(exact_dup, `[[`, character(1), 1L)
  if (length(reps) >= 2L) {
    splits_reps <- splits[reps, , drop = FALSE]
    Lr <- splits_to_integer_matrix(splits_reps)
    Jr <- compute_jaccard_similarity_matrix(Lr)
    pairs_df <- top_similar_pairs(Jr, n = opt$top_pairs)
    pairs_df$note <- "representatives_of_exact_duplicate_groups"
    write.csv(
      pairs_df,
      file.path(out_dir, "split_agreement_summary.csv"),
      row.names = FALSE
    )
  }
}

# Save splits used for reproducibility
saveRDS(splits, file.path(out_dir, paste0(opt$dataset_name, "_splits_used.rds")))
if (mode == "panel") {
  writeLines(rownames(splits), file.path(out_dir, "panel_genes_used.txt"))
}

message("Done. Outputs in: ", out_dir)
