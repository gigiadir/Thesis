# build_gene_cci_sets: for each gene in a .malignant list, return CCIs ranked by
# |LOGFC|. top_n = NULL returns all interactions (used by Average Overlap); an
# integer truncates to the top-N (used by Jaccard).
build_gene_cci_sets <- function(malignant_list, top_n = TOP_N_CCI) {
  map(malignant_list, function(df) {
    df <- df %>%
      as.data.frame() %>%
      arrange(-abs(LOGFC))
    if (!is.null(top_n)) df <- slice_head(df, n = top_n)
    pull(df, CCI) %>% unique()
  })
}

build_cci_logfc_mat <- function(malignant_list, min_genes = MIN_GENES_PER_CCI) {
  cci_list <- lapply(names(malignant_list), function(gene) {
    malignant_list[[gene]] %>%
      select(CCI, LOGFC) %>%
      group_by(CCI) %>%
      summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop") %>%
      mutate(Gene = gene)
  })

  all_data <- bind_rows(cci_list)
  mat_long <- all_data %>%
    pivot_wider(names_from = Gene, values_from = LOGFC)

  mat <- as.matrix(mat_long[, -1])
  rownames(mat) <- mat_long$CCI
  message("  Matrix dimensions (CCIs x Genes): ", nrow(mat), " x ", ncol(mat))

  n_non_na <- rowSums(!is.na(mat))
  mat <- mat[n_non_na >= min_genes, , drop = FALSE]
  message("  After filtering sparse CCIs: ", nrow(mat), " x ", ncol(mat))
  mat
}

# Build Gene×Dataset × CCI logFC matrix (NA where CCI absent for that gene×dataset).
build_gene_dataset_cci_feat_mat <- function(cci_mats, cancer_map) {
  all_ccis <- sort(unique(unlist(lapply(cci_mats, rownames))))
  gd_rows <- purrr::imap(cci_mats, function(mat, ds) {
    gene_names <- colnames(mat)
    lapply(gene_names, function(g) {
      vec <- setNames(rep(NA_real_, length(all_ccis)), all_ccis)
      present <- rownames(mat)[!is.na(mat[, g])]
      vec[present] <- mat[present, g]
      list(
        label   = paste0(g, "_", ds),
        gene    = g,
        dataset = ds,
        cancer  = cancer_map[[ds]],
        vec     = vec
      )
    })
  }) %>% purrr::flatten()

  non_empty <- sapply(gd_rows, function(x) any(!is.na(x$vec)))
  gd_rows   <- gd_rows[non_empty]

  feat_mat <- do.call(rbind, lapply(gd_rows, `[[`, "vec"))
  rownames(feat_mat) <- vapply(gd_rows, `[[`, "", "label")

  list(
    feat_mat    = feat_mat,
    all_ccis    = all_ccis,
    labels_gd   = vapply(gd_rows, `[[`, "", "label"),
    genes_gd    = vapply(gd_rows, `[[`, "", "gene"),
    datasets_gd = vapply(gd_rows, `[[`, "", "dataset"),
    cancer_gd   = vapply(gd_rows, `[[`, "", "cancer")
  )
}

# HVG CCIs by cross Gene×Dataset logFC variance (matches post-analysis-hvg-ccis).
select_hvg_ccis <- function(malignant_by_ds, cancer_map,
                            top_pct = 0.05, min_obs = 10) {
  cci_mats <- lapply(malignant_by_ds, build_cci_logfc_mat)
  feat_mat <- build_gene_dataset_cci_feat_mat(cci_mats, cancer_map)$feat_mat
  n_obs_per_cci <- colSums(!is.na(feat_mat))
  cci_var <- apply(feat_mat, 2, function(x) var(x, na.rm = TRUE))
  cci_var[n_obs_per_cci < min_obs] <- NA_real_
  eligible <- !is.na(cci_var)
  if (!any(eligible)) {
    stop("No CCIs passed min_obs filter for HVG selection.")
  }
  var_cutoff <- quantile(cci_var[eligible], probs = 1 - top_pct, na.rm = TRUE)
  hvg <- names(cci_var)[eligible & cci_var >= var_cutoff]
  list(
    hvg         = hvg,
    n_eligible  = sum(eligible),
    var_cutoff  = var_cutoff,
    cci_var     = cci_var
  )
}

restrict_ccis_to_hvg <- function(ccis, hvg_ccis) {
  ccis[ccis %in% hvg_ccis]
}

# jaccard_dist: scalar Jaccard distance between two character-vector CCI sets.
jaccard_dist <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  union <- length(union(set_a, set_b))
  if (union == 0) return(NA_real_)
  1 - inter / union
}

# build_gene_dataset_jaccard: N×N Jaccard distance matrix + pairwise CCI intersections
# for one gene across N datasets.
# Returns a list:
#   $dist_mat      — N×N numeric Jaccard distance matrix
#   $intersections — N×N list-matrix of shared CCI character vectors (diagonal = own CCIs)
build_gene_dataset_jaccard <- function(gene, all_sets, labels = DATASET_LABELS) {
  sets    <- map(all_sets, function(ds_sets) ds_sets[[gene]])
  present <- map_lgl(sets, Negate(is.null))
  if (!all(present)) {
    warning(sprintf("Gene '%s' absent in: %s",
                    gene, paste(labels[!present], collapse = ", ")))
  }
  n            <- length(labels)
  dist_mat     <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(labels, labels))
  intersections <- matrix(vector("list", n * n), nrow = n, ncol = n,
                          dimnames = list(labels, labels))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) {
      dist_mat[i, j]      <- NA_real_
      intersections[[i, j]] <- sets[[i]]
    } else if (!is.null(sets[[i]]) && !is.null(sets[[j]])) {
      dist_mat[i, j]      <- jaccard_dist(sets[[i]], sets[[j]])
      intersections[[i, j]] <- intersect(sets[[i]], sets[[j]])
    }
  }
  list(dist_mat = dist_mat, intersections = intersections)
}


# pad_list / dist_fun_superranker_overlap: Average Overlap distance helpers.
# pad_list fills a ranked list to target_len with unique placeholder tokens so
# that SuperRanker::average_overlap receives two same-length rank matrices.
pad_list <- function(lst, target_len) {
  if (length(lst) < target_len) {
    padding <- paste0("MISSING_", seq_len(target_len - length(lst)), "__")
    return(c(lst, padding))
  }
  return(lst[1:target_len])
}

dist_fun_superranker_overlap <- function(list1, list2, k = NULL) {
  if (is.null(k)) k <- max(length(list1), length(list2))
  list1     <- pad_list(list1, k)
  list2     <- pad_list(list2, k)
  all_items <- unique(c(list1, list2))
  ranks1    <- match(all_items, list1);  ranks1[is.na(ranks1)] <- k + 1
  ranks2    <- match(all_items, list2);  ranks2[is.na(ranks2)] <- k + 1
  rank_matrix    <- cbind(ranks1, ranks2)
  overlap_values <- SuperRanker::average_overlap(rank_matrix)
  similarity     <- mean(overlap_values[1:min(k, length(overlap_values))], na.rm = TRUE)
  1 - similarity
}

# build_gene_dataset_ao: N×N Average Overlap distance matrix for one gene across
# N datasets. Uses full ranked CCI lists (no top-N truncation).
build_gene_dataset_ao <- function(gene, all_ranked_sets, labels = DATASET_LABELS) {
  sets    <- map(all_ranked_sets, function(ds_sets) ds_sets[[gene]])
  present <- map_lgl(sets, Negate(is.null))
  if (!all(present)) {
    warning(sprintf("Gene '%s' absent in: %s",
                    gene, paste(labels[!present], collapse = ", ")))
  }
  n        <- length(labels)
  dist_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(labels, labels))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) {
      dist_mat[i, j] <- NA_real_
    } else if (!is.null(sets[[i]]) && !is.null(sets[[j]])) {
      dist_mat[i, j] <- dist_fun_superranker_overlap(sets[[i]], sets[[j]])
    }
  }
  dist_mat
}
