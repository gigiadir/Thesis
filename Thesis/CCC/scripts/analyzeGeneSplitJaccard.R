#!/usr/bin/env Rscript

# Compare patient LOW/MID/HIGH splits across genes via Jaccard similarity.
#
# Run order:
#   1. Rscript createPseudobulkMatrix.R
#   2. Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name Kurten_HNSC
#   3. Rscript analyzeGeneSplitJaccard.R --dataset_name Kurten_HNSC --mode panel
#
# Requires: optparse, pheatmap, ggplot2, ggrepel, ggforce
# R libraries (if packages missing from plain Rscript):
#   export R_LIBS_SITE=/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0
#   or uncomment source("/gpfs0/bgu-ofircohen/group/groupRprofile") in ~/.Rprofile
#
# Modes:
#   panel — 123-gene scDiffCom panel: Jaccard matrix + heatmap + MDS + split-profile PCA
#   all   — all genes: duplicate-split clusters; full matrix only if <= max_genes_full_matrix

suppressPackageStartupMessages({
  library(optparse)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages(\"ggplot2\").")
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with install.packages(\"pheatmap\").")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install with install.packages(\"ggrepel\").")
  }
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required. Install with install.packages(\"ggforce\").")
  }
  library(pheatmap)
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

row_zscore_matrix <- function(mat) {
  z <- mat
  storage.mode(z) <- "double"
  for (i in seq_len(nrow(z))) {
    v <- z[i, ]
    ok <- is.finite(v) & (v != 0)
    if (sum(ok) < 2L) next
    mu <- mean(v[ok])
    sdv <- stats::sd(v[ok])
    if (!is.finite(sdv) || sdv < 1e-8) next
    z[i, ok] <- (v[ok] - mu) / sdv
    z[i, !ok] <- 0
  }
  z
}

prepare_split_profile_matrix <- function(L) {
  keep_pat <- colSums(L > 0L) > 0L
  L <- L[, keep_pat, drop = FALSE]
  keep_genes <- rowSums(L > 0L) >= 2L
  if (!any(keep_genes)) {
    stop("No genes with >= 2 non-missing split assignments for PCA.")
  }
  L <- L[keep_genes, , drop = FALSE]
  z <- row_zscore_matrix(L)
  z
}

run_gene_pca <- function(z_mat, n_pc = 2L) {
  pc <- stats::prcomp(z_mat, center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pc$x), nrow(pc$x))
  var_expl <- (pc$sdev^2) / sum(pc$sdev^2)
  list(
    coords = pc$x[, seq_len(n_pc), drop = FALSE],
    var_expl = var_expl,
    genes = rownames(z_mat)
  )
}

genes_in_dup_groups <- function(dup_groups) {
  if (length(dup_groups) == 0L) return(character())
  unique(unlist(dup_groups, use.names = FALSE))
}

assign_dup_group_ids <- function(dup_groups, all_genes) {
  out <- rep(NA_integer_, length(all_genes))
  names(out) <- all_genes
  for (i in seq_along(dup_groups)) {
    members <- dup_groups[[i]]
    out[members] <- i
  }
  out
}

mean_jaccard_per_gene <- function(J) {
  genes <- rownames(J)
  stats::setNames(
    vapply(seq_along(genes), function(i) {
      row_j <- J[i, ]
      row_j <- row_j[!is.na(row_j)]
      row_j <- row_j[names(row_j) != genes[[i]]]
      if (length(row_j) == 0L) return(NA_real_)
      mean(row_j)
    }, numeric(1)),
    genes
  )
}

build_jaccard_edges <- function(J, min_j) {
  genes <- rownames(J)
  g <- length(genes)
  if (g < 2L) return(data.frame())
  rows <- list()
  k <- 1L
  for (i in seq_len(g - 1L)) {
    for (j in (i + 1L):g) {
      sim <- J[i, j]
      if (!is.na(sim) && sim >= min_j) {
        rows[[k]] <- data.frame(
          gene_a = genes[[i]],
          gene_b = genes[[j]],
          jaccard = sim,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      gene_a = character(),
      gene_b = character(),
      jaccard = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

build_embedding_df <- function(coords, x_col, y_col, genes, J, dup_groups) {
  mean_j <- mean_jaccard_per_gene(J)
  dup_genes <- genes_in_dup_groups(dup_groups)
  dup_ids <- assign_dup_group_ids(dup_groups, genes)

  data.frame(
    gene = genes,
    x = coords[, x_col],
    y = coords[, y_col],
    mean_jaccard = mean_j[genes],
    in_dup_group = genes %in% dup_genes,
    dup_group_id = dup_ids[genes],
    stringsAsFactors = FALSE
  )
}

build_edge_plot <- function(edge_df, embed_df) {
  if (nrow(edge_df) == 0L) {
    return(NULL)
  }
  pos_a <- embed_df[, c("gene", "x", "y")]
  colnames(pos_a) <- c("gene_a", "x1", "y1")
  pos_b <- embed_df[, c("gene", "x", "y")]
  colnames(pos_b) <- c("gene_b", "x2", "y2")
  edge_plot <- merge(edge_df, pos_a, by = "gene_a", all.x = TRUE)
  edge_plot <- merge(edge_plot, pos_b, by = "gene_b", all.x = TRUE)
  edge_plot <- edge_plot[
    is.finite(edge_plot$x1) & is.finite(edge_plot$y1) &
      is.finite(edge_plot$x2) & is.finite(edge_plot$y2),
    ,
    drop = FALSE
  ]
  if (nrow(edge_plot) == 0L) {
    return(NULL)
  }
  edge_plot
}

prepare_dup_mark_layers <- function(dup_df) {
  if (nrow(dup_df) == 0L) {
    return(list(hull = NULL, circle = NULL))
  }
  split_dup <- split(dup_df, dup_df$dup_group_id)
  hull_parts <- lapply(split_dup, function(df) {
    if (nrow(df) < 3L) {
      return(NULL)
    }
    if (nrow(unique(df[, c("x", "y")])) < 3L) {
      return(NULL)
    }
    df
  })
  hull_df <- do.call(rbind, hull_parts)
  if (!is.null(hull_df) && nrow(hull_df) == 0L) {
    hull_df <- NULL
  }

  circle_parts <- lapply(split_dup, function(df) {
    if (nrow(df) < 2L) {
      return(NULL)
    }
    if (nrow(unique(df[, c("x", "y")])) >= 3L) {
      return(NULL)
    }
    data.frame(
      dup_group_id = df$dup_group_id[[1]],
      x = mean(df$x),
      y = mean(df$y),
      stringsAsFactors = FALSE
    )
  })
  circle_df <- do.call(rbind, circle_parts)
  if (!is.null(circle_df) && nrow(circle_df) == 0L) {
    circle_df <- NULL
  }

  list(hull = hull_df, circle = circle_df)
}

plot_gene_embedding <- function(embed_df, edge_df, dup_groups, title, subtitle,
                              xlab, ylab, out_path) {
  embed_df <- embed_df[is.finite(embed_df$x) & is.finite(embed_df$y), , drop = FALSE]
  if (nrow(embed_df) == 0L) {
    warning("No finite embedding coordinates; skipping plot: ", out_path)
    return(invisible(NULL))
  }

  dup_df <- embed_df[embed_df$in_dup_group & !is.na(embed_df$dup_group_id), , drop = FALSE]
  dup_layers <- prepare_dup_mark_layers(dup_df)
  edge_plot <- build_edge_plot(edge_df, embed_df)

  p <- ggplot(embed_df, aes(x = x, y = y)) +
    theme_classic(base_size = 13)

  if (!is.null(edge_plot) && nrow(edge_plot) > 0L) {
    p <- p + geom_segment(
      data = edge_plot,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      inherit.aes = FALSE,
      colour = "grey75",
      linewidth = 0.25,
      alpha = 0.6
    )
  }

  if (!is.null(dup_layers$hull) && nrow(dup_layers$hull) > 0L) {
    p <- p + ggforce::geom_mark_hull(
      data = dup_layers$hull,
      aes(x = x, y = y, group = dup_group_id),
      fill = NA,
      colour = "grey40",
      linewidth = 0.4,
      expand = grid::unit(3, "mm"),
      inherit.aes = FALSE
    )
  }

  if (!is.null(dup_layers$circle) && nrow(dup_layers$circle) > 0L) {
    p <- p + ggforce::geom_mark_circle(
      data = dup_layers$circle,
      aes(x = x, y = y, group = dup_group_id),
      fill = NA,
      colour = "grey40",
      linewidth = 0.4,
      radius = grid::unit(4, "mm"),
      inherit.aes = FALSE
    )
  }

  p <- p +
    geom_point(size = 2, alpha = 0.85, colour = "#4575b4") +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab)

  p <- p + ggrepel::geom_text_repel(
    aes(label = gene),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.3,
    segment.size = 0.2,
    segment.alpha = 0.4
  )

  ggsave(out_path, plot = p, width = 12, height = 9, dpi = 300)
}

option_list <- list(
  make_option("--dataset_name", type = "character", default = NULL,
              help = "Dataset name (e.g. Kurten_HNSC)", metavar = "character"),
  make_option("--mode", type = "character", default = "panel",
              help = "panel or all [default %default]", metavar = "character"),
  make_option("--split_source", type = "character", default = "rankgenes",
              help = "Split method: rankgenes or residual [default %default]",
              metavar = "character"),
  make_option("--rankgenes_dir", type = "character",
              default = "~/CCC-PreProcess/results-RankGenes",
              help = "RankGenes output base directory [default %default]",
              metavar = "character"),
  make_option("--residual_dir", type = "character",
              default = "~/CCC-PreProcess/results-Residual",
              help = "Residual-split output base directory [default %default]",
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
  make_option("--no_cluster", action = "store_true", default = FALSE,
              help = "Skip hierarchical clustering and heatmap reordering [default %default]"),
  make_option("--jaccard_threshold", type = "double", default = 1.0,
              help = "Threshold for duplicate/near-duplicate groups [default %default]"),
  make_option("--max_genes_full_matrix", type = "integer", default = 2500L,
              help = "In 'all' mode, skip full pairwise matrix above this G [default %default]"),
  make_option("--top_pairs", type = "integer", default = 100L,
              help = "Number of top similar gene pairs to export [default %default]"),
  make_option("--edge_jaccard_min", type = "double", default = 0.95,
              help = "Draw edges between gene pairs with Jaccard >= this [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$dataset_name) || !nzchar(opt$dataset_name)) {
  stop("--dataset_name is required.")
}

mode <- tolower(opt$mode)
if (!mode %in% c("panel", "all")) {
  stop("--mode must be 'panel' or 'all'.")
}

split_source <- tolower(opt$split_source)
if (!split_source %in% c("rankgenes", "residual")) {
  stop("--split_source must be 'rankgenes' or 'residual'.")
}

rankgenes_dir <- path.expand(opt$rankgenes_dir)
residual_dir <- path.expand(opt$residual_dir)
splits_dir <- if (split_source == "residual") residual_dir else rankgenes_dir
pseudobulk_dir <- path.expand(opt$pseudobulk_dir)
out_dir <- if (!is.null(opt$output_dir) && nzchar(opt$output_dir)) {
  path.expand(opt$output_dir)
} else {
  path.expand(file.path("~/Thesis/CCC/outputs/split_similarity", opt$dataset_name))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_splits_for_analysis <- function() {
  ds_dir <- file.path(splits_dir, opt$dataset_name)
  all_splits_path <- file.path(ds_dir, paste0(opt$dataset_name, ALL_SPLITS_SUFFIX))
  load_from_dir <- if (split_source == "residual") {
    load_splits_from_residual_dir
  } else {
    load_splits_from_rankgenes_dir
  }
  source_label <- if (split_source == "residual") "Residual" else "RankGenes"

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
    message("Loading splits for ", length(panel_genes), " panel genes from ", source_label, " ...")
    splits <- load_from_dir(splits_dir, opt$dataset_name, genes = panel_genes)
    missing <- setdiff(panel_genes, rownames(splits))
    if (length(missing) > 0L) {
      warning(
        "Missing grouped .rds for ", length(missing), " gene(s); continuing with ",
        nrow(splits), " genes.", call. = FALSE, immediate. = TRUE
      )
    }
    return(splits[intersect(panel_genes, rownames(splits)), , drop = FALSE])
  }

  message("Recomputing splits from pseudobulk (all genes, ", split_source, ") ...")
  in_path <- resolve_pseudobulk_path(opt$dataset_name, pseudobulk_dir)
  mat <- load_gene_patient_matrix(in_path)
  if (split_source == "residual") {
    build_residual_splits_matrix(mat, n_pc = 2L)$splits
  } else {
    rank_mat <- compute_rank_matrix(mat)
    build_patient_splits_matrix(mat, rank_mat = rank_mat)
  }
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

    D_plot <- dist_mat
    diag(D_plot) <- NA_real_
    heatmap_path <- file.path(out_dir, "jaccard_heatmap.png")
    label_fs <- max(4, 200 / nrow(D_plot))
    png(heatmap_path, width = 4000, height = 4000, res = 300)
    pheatmap(
      D_plot,
      color = colorRampPalette(c("#d73027", "white", "#4575b4"))(100),
      breaks = seq(0, 1, length.out = 101),
      na_col = "#E8E8E8",
      cluster_rows = if (isTRUE(opt$no_cluster)) FALSE else hc,
      cluster_cols = if (isTRUE(opt$no_cluster)) FALSE else hc,
      fontsize_row = label_fs,
      fontsize_col = label_fs,
      angle_col = 45,
      main = paste0(
        opt$dataset_name, " \u2013 Gene-Gene Jaccard Split Distance (n=", G, ")"
      )
    )
    dev.off()
    message("Wrote heatmap: ", heatmap_path)

    genes <- rownames(J)
    edge_df <- build_jaccard_edges(J, opt$edge_jaccard_min)

    k_mds <- min(2L, G - 1L)
    mds_stress <- NA_real_
    if (k_mds >= 1L) {
      mds_fit <- stats::cmdscale(d, k = k_mds, eig = TRUE)
      mds_coords <- mds_fit$points
      if (k_mds == 1L) {
        mds_coords <- cbind(mds_coords, 0)
      }
      colnames(mds_coords) <- c("MDS1", "MDS2")
      if (!is.null(mds_fit$GOF) && length(mds_fit$GOF) >= 2L) {
        mds_stress <- round(mds_fit$GOF[[2]], 4)
      }

      mds_embed <- build_embedding_df(
        mds_coords, 1L, 2L, genes, J, dup_groups
      )
      mds_out <- mds_embed
      colnames(mds_out)[colnames(mds_out) == "x"] <- "MDS1"
      colnames(mds_out)[colnames(mds_out) == "y"] <- "MDS2"
      saveRDS(mds_out, file.path(out_dir, "jaccard_mds_coords.rds"))

      mds_subtitle <- paste(
        "Classical MDS on Jaccard distance (1 \u2212 similarity);",
        "closer = more similar splits"
      )
      if (!is.na(mds_stress)) {
        mds_subtitle <- paste(mds_subtitle, "| GOF stress =", mds_stress)
      }

      plot_gene_embedding(
        embed_df = mds_embed,
        edge_df = edge_df,
        dup_groups = dup_groups,
        title = paste0(
          opt$dataset_name, " \u2013 Gene split agreement (MDS, n=", G, ")"
        ),
        subtitle = mds_subtitle,
        xlab = "MDS1",
        ylab = "MDS2",
        out_path = file.path(out_dir, "jaccard_mds.png")
      )
      message("Wrote MDS plot: ", file.path(out_dir, "jaccard_mds.png"))
    }

    z_mat <- prepare_split_profile_matrix(L)
    pca_res <- run_gene_pca(z_mat, n_pc = 2L)
    pca_coords <- pca_res$coords
    colnames(pca_coords) <- c("PC1", "PC2")
    var_pc1 <- round(100 * pca_res$var_expl[1], 1)
    var_pc2 <- if (length(pca_res$var_expl) >= 2L) {
      round(100 * pca_res$var_expl[2], 1)
    } else {
      NA_real_
    }

    pca_embed <- build_embedding_df(
      pca_coords, 1L, 2L, pca_res$genes, J, dup_groups
    )
    pca_out <- pca_embed
    colnames(pca_out)[colnames(pca_out) == "x"] <- "PC1"
    colnames(pca_out)[colnames(pca_out) == "y"] <- "PC2"
    pca_out$var_pc1 <- var_pc1
    pca_out$var_pc2 <- var_pc2
    saveRDS(pca_out, file.path(out_dir, "jaccard_pca_coords.rds"))

    plot_gene_embedding(
      embed_df = pca_embed,
      edge_df = edge_df,
      dup_groups = dup_groups,
      title = paste0(
        opt$dataset_name, " \u2013 Gene split profiles (PCA, n=", nrow(pca_coords), ")"
      ),
      subtitle = paste(
        "Row z-scored patient LOW/MID/HIGH profiles;",
        sprintf("PC1 %.1f%% | PC2 %.1f%% variance", var_pc1, var_pc2)
      ),
      xlab = sprintf("PC1 (%.1f%% var)", var_pc1),
      ylab = sprintf("PC2 (%.1f%% var)", var_pc2),
      out_path = file.path(out_dir, "jaccard_pca.png")
    )
    message("Wrote PCA plot: ", file.path(out_dir, "jaccard_pca.png"))
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
