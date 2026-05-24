#!/usr/bin/env Rscript

# Compare patient LOW/MID/HIGH splits across genes via Jaccard similarity.
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name Kurten_HNSC
#   3. Rscript analyzeGeneSplitJaccard.R --dataset_name Kurten_HNSC --mode panel
#
# Requires: optparse, ggplot2, ggrepel
# R libraries (if packages missing from plain Rscript):
#   export R_LIBS_SITE=/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0
#   or uncomment source("/gpfs0/bgu-ofircohen/group/groupRprofile") in ~/.Rprofile
#
# Modes:
#   panel — 123-gene scDiffCom panel: Jaccard matrix + MDS plot
#   all   — all genes: duplicate-split clusters; full matrix only if <= max_genes_full_matrix

suppressPackageStartupMessages({
  library(optparse)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages(\"ggplot2\").")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with install.packages(\"ggrepel\").")
  }
  library(ggplot2)
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
              help = "Optional .rds or .txt gene list (default: 123-gene scDiffCom panel)",
              metavar = "character"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory [default ~/Thesis/CCC/outputs/split_similarity/{dataset}]",
              metavar = "character"),
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

load_splits_for_analysis <- function() {
  ds_dir <- file.path(rankgenes_dir, opt$dataset_name)
  all_splits_path <- file.path(ds_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))

  panel_genes <- if (!is.null(opt$gene_list) && nzchar(opt$gene_list)) {
    load_gene_list(opt$gene_list)
  } else {
    SCDIFFCOM_GENE_PANEL
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

  if (G >= 2L) {
    d <- as.dist(dist_mat)
    d[is.na(d)] <- 1
    hc <- hclust(d, method = "average")
    saveRDS(hc, file.path(out_dir, "gene_clusters.rds"))

    k_mds <- min(2L, G - 1L)
    if (k_mds >= 1L) {
      mds_coords <- cmdscale(d, k = k_mds)
      if (k_mds == 1L) {
        mds_coords <- cbind(mds_coords, 0)
      }
      mds_df <- data.frame(
        gene = rownames(J),
        MDS1 = mds_coords[, 1],
        MDS2 = mds_coords[, 2],
        stringsAsFactors = FALSE
      )
      saveRDS(mds_df, file.path(out_dir, "jaccard_mds_coords.rds"))

      p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
        geom_point(size = 2, alpha = 0.85, colour = "#4575b4") +
        ggrepel::geom_text_repel(
          aes(label = gene),
          size = 2.8,
          max.overlaps = Inf,
          box.padding = 0.3,
          segment.size = 0.2,
          segment.alpha = 0.4
        ) +
        theme_classic(base_size = 13) +
        labs(
          title = paste0(
            opt$dataset_name, " – Gene split agreement (MDS, n=", G, ")"
          ),
          subtitle = paste(
            "2D classical MDS on Jaccard distance (1 − similarity);",
            "closer points = more similar splits"
          )
        )

      mds_path <- file.path(out_dir, "jaccard_mds.png")
      ggsave(mds_path, plot = p_mds, width = 12, height = 9, dpi = 300)
      message("Wrote MDS plot: ", mds_path)
    }
  }
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
