sample_pair_indices <- function(n_items, n_pairs) {
  if (n_items < 2L || n_pairs <= 0L) {
    return(matrix(integer(0), nrow = 0L, ncol = 2L))
  }

  total_pairs <- (n_items * (n_items - 1L)) %/% 2L

  if (n_pairs >= total_pairs) {
    return(t(utils::combn(n_items, 2L)))
  }

  seen <- new.env(parent = emptyenv())
  out <- matrix(0L, nrow = n_pairs, ncol = 2L)
  n_collected <- 0L
  max_attempts <- max(2000L, n_pairs * 20L)
  attempts <- 0L

  while (n_collected < n_pairs && attempts < max_attempts) {
    attempts <- attempts + 1L
    ij <- sort(sample.int(n_items, size = 2L, replace = FALSE))
    key <- paste(ij[[1]], ij[[2]], sep = "|")
    if (!exists(key, envir = seen, inherits = FALSE)) {
      assign(key, TRUE, envir = seen)
      n_collected <- n_collected + 1L
      out[n_collected, ] <- ij
    }
  }

  if (n_collected == 0L) {
    return(matrix(integer(0), nrow = 0L, ncol = 2L))
  }
  out[seq_len(n_collected), , drop = FALSE]
}

pair_key <- function(a, b) {
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

prepare_reference_sets <- function(reference_sets, genes_available, min_set_size = 3L) {
  if (!is.list(reference_sets) || length(reference_sets) == 0L) {
    stop("reference_sets must be a non-empty list: set_name -> character genes.")
  }

  if (is.null(names(reference_sets))) {
    names(reference_sets) <- paste0("set_", seq_along(reference_sets))
  }

  cleaned <- lapply(reference_sets, function(g) {
    unique(as.character(g[!is.na(g) & nzchar(g)]))
  })

  filtered <- lapply(cleaned, intersect, y = genes_available)
  sizes_before <- lengths(cleaned)
  sizes_after <- lengths(filtered)
  keep <- sizes_after >= min_set_size

  list(
    sets = filtered[keep],
    dropped = data.frame(
      set_name = names(filtered)[!keep],
      original_size = unname(sizes_before[!keep]),
      intersected_size = unname(sizes_after[!keep]),
      stringsAsFactors = FALSE
    ),
    coverage = data.frame(
      set_name = names(filtered),
      original_size = unname(sizes_before),
      intersected_size = unname(sizes_after),
      coverage_fraction = ifelse(sizes_before > 0L, sizes_after / sizes_before, NA_real_),
      stringsAsFactors = FALSE
    )
  )
}

build_gene_to_sets <- function(reference_sets) {
  g2s <- list()
  for (set_name in names(reference_sets)) {
    genes <- reference_sets[[set_name]]
    for (g in genes) {
      if (is.null(g2s[[g]])) {
        g2s[[g]] <- set_name
      } else {
        g2s[[g]] <- c(g2s[[g]], set_name)
      }
    }
  }
  g2s
}

build_pair_table <- function(reference_sets,
                             all_genes,
                             max_positive_pairs = NULL,
                             max_positive_pairs_per_set = NULL,
                             negative_ratio = 1,
                             seed = 1L) {
  set.seed(seed)
  genes <- unique(all_genes)
  if (length(genes) < 2L) {
    stop("Need at least 2 genes to build pair table.")
  }

  positive_env <- new.env(parent = emptyenv())
  positive_rows <- list()
  row_id <- 0L

  for (set_name in names(reference_sets)) {
    members <- intersect(reference_sets[[set_name]], genes)
    if (length(members) < 2L) next

    n_possible <- (length(members) * (length(members) - 1L)) %/% 2L
    n_take <- n_possible
    if (!is.null(max_positive_pairs_per_set) && max_positive_pairs_per_set > 0L) {
      n_take <- min(n_take, as.integer(max_positive_pairs_per_set))
    }

    ij <- sample_pair_indices(length(members), n_take)
    if (nrow(ij) == 0L) next

    for (k in seq_len(nrow(ij))) {
      a <- members[ij[k, 1L]]
      b <- members[ij[k, 2L]]
      key <- pair_key(a, b)
      if (!exists(key, envir = positive_env, inherits = FALSE)) {
        assign(key, set_name, envir = positive_env)
        row_id <- row_id + 1L
        positive_rows[[row_id]] <- data.frame(
          gene_a = min(a, b),
          gene_b = max(a, b),
          label = 1L,
          origin_set = set_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(positive_rows) == 0L) {
    stop("No within-set pairs found after filtering sets/genes.")
  }

  positive_df <- do.call(rbind, positive_rows)
  if (!is.null(max_positive_pairs) && max_positive_pairs > 0L && nrow(positive_df) > max_positive_pairs) {
    idx <- sample.int(nrow(positive_df), size = max_positive_pairs, replace = FALSE)
    positive_df <- positive_df[idx, , drop = FALSE]
  }

  n_negative_target <- max(1L, as.integer(ceiling(nrow(positive_df) * negative_ratio)))
  selected_env <- new.env(parent = emptyenv())
  for (k in seq_len(nrow(positive_df))) {
    assign(pair_key(positive_df$gene_a[[k]], positive_df$gene_b[[k]]), TRUE, envir = selected_env)
  }

  negative_rows <- list()
  n_neg <- 0L
  max_attempts <- max(10000L, n_negative_target * 30L)
  attempts <- 0L
  gene_to_sets <- build_gene_to_sets(reference_sets)

  while (n_neg < n_negative_target && attempts < max_attempts) {
    attempts <- attempts + 1L
    ab <- sample(genes, size = 2L, replace = FALSE)
    a <- ab[[1]]
    b <- ab[[2]]
    key <- pair_key(a, b)

    if (exists(key, envir = selected_env, inherits = FALSE)) {
      next
    }
    sets_a <- gene_to_sets[[a]]
    sets_b <- gene_to_sets[[b]]
    if (!is.null(sets_a) && !is.null(sets_b) && length(intersect(sets_a, sets_b)) > 0L) {
      next
    }

    n_neg <- n_neg + 1L
    assign(key, TRUE, envir = selected_env)
    negative_rows[[n_neg]] <- data.frame(
      gene_a = min(a, b),
      gene_b = max(a, b),
      label = 0L,
      origin_set = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  if (n_neg < n_negative_target) {
    warning(
      "Could only sample ", n_neg, " negative pairs out of requested ", n_negative_target,
      call. = FALSE, immediate. = TRUE
    )
  }

  negative_df <- if (length(negative_rows) > 0L) {
    do.call(rbind, negative_rows)
  } else {
    data.frame(gene_a = character(), gene_b = character(), label = integer(), origin_set = character())
  }

  rbind(positive_df, negative_df)
}

prepare_distance_matrix <- function(input_matrix,
                                    input_mode = c("feature_matrix", "distance_matrix"),
                                    similarity_input = FALSE,
                                    feature_distance = c("euclidean", "correlation")) {
  input_mode <- match.arg(input_mode)
  feature_distance <- match.arg(feature_distance)

  if (is.null(rownames(input_matrix)) || any(!nzchar(rownames(input_matrix)))) {
    stop("Input matrix must have non-empty row names as gene IDs.")
  }
  if (anyDuplicated(rownames(input_matrix)) > 0L) {
    stop("Input matrix has duplicated row names (gene IDs).")
  }

  if (input_mode == "distance_matrix") {
    m <- as.matrix(input_matrix)
    if (nrow(m) != ncol(m)) {
      stop("distance_matrix mode requires a square matrix.")
    }
    if (is.null(colnames(m))) {
      colnames(m) <- rownames(m)
    }
    if (!identical(sort(rownames(m)), sort(colnames(m)))) {
      stop("distance_matrix mode requires matching row/column gene names.")
    }
    m <- m[rownames(m), rownames(m), drop = FALSE]
    d <- if (isTRUE(similarity_input)) 1 - m else m
    diag(d) <- 0
    return(d)
  }

  x <- as.matrix(input_matrix)
  storage.mode(x) <- "double"
  if (feature_distance == "euclidean") {
    d <- as.matrix(stats::dist(x, method = "euclidean"))
  } else {
    cor_mat <- stats::cor(t(x), use = "pairwise.complete.obs", method = "pearson")
    cor_mat[!is.finite(cor_mat)] <- 0
    d <- 1 - cor_mat
  }
  diag(d) <- 0
  d
}

extract_pair_distances <- function(pair_df, distance_matrix) {
  if (nrow(pair_df) == 0L) return(numeric())
  idx_a <- match(pair_df$gene_a, rownames(distance_matrix))
  idx_b <- match(pair_df$gene_b, colnames(distance_matrix))
  ok <- !is.na(idx_a) & !is.na(idx_b)
  out <- rep(NA_real_, nrow(pair_df))
  out[ok] <- distance_matrix[cbind(idx_a[ok], idx_b[ok])]
  out
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  m <- length(x)
  n <- length(y)
  if (m == 0L || n == 0L) return(NA_real_)
  r <- rank(c(x, y), ties.method = "average")
  u <- sum(r[seq_len(m)]) - (m * (m + 1)) / 2
  (2 * u) / (m * n) - 1
}

auc_from_distance <- function(within_dist, between_dist) {
  pos_scores <- -within_dist[is.finite(within_dist)]
  neg_scores <- -between_dist[is.finite(between_dist)]
  m <- length(pos_scores)
  n <- length(neg_scores)
  if (m == 0L || n == 0L) return(NA_real_)
  r <- rank(c(pos_scores, neg_scores), ties.method = "average")
  u <- sum(r[seq_len(m)]) - (m * (m + 1)) / 2
  u / (m * n)
}

bootstrap_ci <- function(within_dist, between_dist, n_boot = 500L, seed = 1L) {
  set.seed(seed)
  x <- within_dist[is.finite(within_dist)]
  y <- between_dist[is.finite(between_dist)]
  if (length(x) == 0L || length(y) == 0L || n_boot < 10L) {
    return(list(
      ratio = c(NA_real_, NA_real_),
      within_mean = c(NA_real_, NA_real_),
      between_mean = c(NA_real_, NA_real_)
    ))
  }

  ratio <- numeric(n_boot)
  wx <- numeric(n_boot)
  by <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    sx <- sample(x, size = length(x), replace = TRUE)
    sy <- sample(y, size = length(y), replace = TRUE)
    mx <- mean(sx)
    my <- mean(sy)
    ratio[[b]] <- if (isTRUE(all.equal(my, 0))) NA_real_ else mx / my
    wx[[b]] <- mx
    by[[b]] <- my
  }

  list(
    ratio = as.numeric(stats::quantile(ratio, probs = c(0.025, 0.975), na.rm = TRUE)),
    within_mean = as.numeric(stats::quantile(wx, probs = c(0.025, 0.975), na.rm = TRUE)),
    between_mean = as.numeric(stats::quantile(by, probs = c(0.025, 0.975), na.rm = TRUE))
  )
}

permutation_p_value <- function(within_dist, between_dist, n_perm = 1000L, seed = 1L) {
  x <- within_dist[is.finite(within_dist)]
  y <- between_dist[is.finite(between_dist)]
  if (length(x) == 0L || length(y) == 0L || n_perm < 10L) return(NA_real_)
  observed <- mean(x) / mean(y)
  pooled <- c(x, y)
  nx <- length(x)

  set.seed(seed)
  perm_ratio <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    idx <- sample.int(length(pooled), size = nx, replace = FALSE)
    x_perm <- pooled[idx]
    y_perm <- pooled[-idx]
    perm_ratio[[i]] <- mean(x_perm) / mean(y_perm)
  }

  (sum(perm_ratio <= observed) + 1) / (length(perm_ratio) + 1)
}

compute_method_metrics <- function(pair_df,
                                   distance_matrix,
                                   bootstrap_iters = 500L,
                                   permutation_iters = 0L,
                                   seed = 1L) {
  d <- extract_pair_distances(pair_df, distance_matrix)
  pair_df$distance <- d
  pair_df <- pair_df[is.finite(pair_df$distance), , drop = FALSE]

  within <- pair_df$distance[pair_df$label == 1L]
  between <- pair_df$distance[pair_df$label == 0L]

  mean_within <- mean(within)
  mean_between <- mean(between)
  ratio <- if (isTRUE(all.equal(mean_between, 0))) NA_real_ else mean_within / mean_between
  ci <- bootstrap_ci(within, between, n_boot = bootstrap_iters, seed = seed)
  auc <- auc_from_distance(within, between)
  delta <- cliffs_delta(within, between)
  p_perm <- if (permutation_iters > 0L) {
    permutation_p_value(within, between, n_perm = permutation_iters, seed = seed + 1L)
  } else {
    NA_real_
  }

  summary_row <- data.frame(
    n_pairs_within = sum(pair_df$label == 1L),
    n_pairs_between = sum(pair_df$label == 0L),
    mean_within = mean_within,
    median_within = stats::median(within),
    mean_between = mean_between,
    median_between = stats::median(between),
    separation_ratio = ratio,
    ratio_ci_low = ci$ratio[[1]],
    ratio_ci_high = ci$ratio[[2]],
    mean_within_ci_low = ci$within_mean[[1]],
    mean_within_ci_high = ci$within_mean[[2]],
    mean_between_ci_low = ci$between_mean[[1]],
    mean_between_ci_high = ci$between_mean[[2]],
    cliffs_delta = delta,
    auc = auc,
    permutation_p = p_perm,
    stringsAsFactors = FALSE
  )

  list(summary = summary_row, pair_distances = pair_df)
}

compute_set_level_diagnostics <- function(reference_sets, distance_matrix) {
  rows <- lapply(names(reference_sets), function(set_name) {
    genes <- intersect(reference_sets[[set_name]], rownames(distance_matrix))
    if (length(genes) < 2L) return(NULL)
    pairs <- utils::combn(genes, 2L)
    dvals <- distance_matrix[cbind(match(pairs[1, ], rownames(distance_matrix)),
                                   match(pairs[2, ], colnames(distance_matrix)))]
    data.frame(
      set_name = set_name,
      n_genes = length(genes),
      n_pairs = length(dvals),
      mean_distance = mean(dvals, na.rm = TRUE),
      median_distance = stats::median(dvals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame(
      set_name = character(),
      n_genes = integer(),
      n_pairs = integer(),
      mean_distance = numeric(),
      median_distance = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}
