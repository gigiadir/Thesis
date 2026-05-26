# Shared helpers for RankGenes, Residual, ExpressionQuantile, and split-similarity analysis.
# Sourced by scDiffCom-Preprocess-RankGenes.R, scDiffCom-Preprocess-Residual.R,
# scDiffCom-Preprocess-ExpressionQuantile.R, compareExpressionSplitEquivalence.R,
# and analyzeGeneSplitJaccard.R.

PSEUDOBULK_MATRIX_SUFFIX <- "_pseudobulk_matrix.rds"
ALL_SPLITS_SUFFIX <- "_all_patient_splits.rds"
RESIDUAL_MATRIX_SUFFIX <- "_residual_matrix.rds"

# scDiffCom panel (keep in sync with Main-scDiffComPipeline.R)
SCDIFFCOM_GENE_PANEL <- c(
  "ABI1", "ACTB", "ACTG1", "APH1A", "ARID1A", "ARID1B", "ARID2", "AXL",
  "BARD1", "BAZ1A", "BCL6", "BLM", "BPTF", "BRCA1", "BRIP1", "BUB1B",
  "CARM1", "CCNC", "CD28", "CDC27", "CDH1", "CDH11", "CDH2", "CDK8",
  "CDKN2A", "CHD4", "COL1A1", "COL1A2", "CREBBP", "CSF1", "CSF1R", "CSF2",
  "CTLA4", "CTNNB1", "CUL3", "CUL4A", "CUL7", "CYFIP1", "DNMT3B", "EGFR",
  "EP300", "ERBB2", "ERCC2", "ERCC3", "ESR1", "FANCA", "FANCC", "FANCE",
  "FANCF", "FANCG", "FANCL", "GATA3", "GNA11", "GNAQ", "GNB1", "GPS2",
  "HDAC1", "HDAC2", "HDAC4", "HDAC7", "HLA-A", "HLA-B", "HLA-C", "IGF1",
  "IGF2", "LDB1", "LMO2", "LTB", "MAPK1", "MAPK3", "MDM2", "MED1", "MED12",
  "MKI67", "MLH1", "MNAT1", "NBN", "NCOA3", "NCOR1", "NCOR2", "NCSTN",
  "NDC80", "NDUFB9", "NPM1", "NUF2", "NUP98", "PARP1", "PBRM1", "PHC2",
  "PIGA", "POLR2A", "PSMB2", "RABEP1", "RAD21", "RAD50", "RAD51B", "RAD51C",
  "SDC4", "SIN3A", "SKP2", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCB1",
  "SMARCD1", "SMARCE1", "SMC1A", "STAG1", "STAG2", "STAT1", "STAT2",
  "TBL1XR1", "TCEB1", "THBS1", "THRAP3", "TP53", "TRAK1", "VEGFA", "VHL",
  "VWF", "XRCC1", "XRCC2", "YY1"
)

LABEL_LEVELS <- c("LOW", "MID", "HIGH")

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

# Quantile-based LOW/MID/HIGH (matches scDiffComPreprocess.R).
assign_quantile_groups <- function(scores, low_q = 1/3, high_q = 2/3) {
  if (low_q >= high_q) {
    stop("low_q must be strictly less than high_q. Got: ", low_q, ", ", high_q)
  }
  out <- rep(NA_character_, length(scores))
  ok <- is.finite(scores)
  if (!any(ok)) return(out)

  qs <- stats::quantile(scores[ok], probs = c(low_q, high_q), na.rm = TRUE)
  out[ok] <- ifelse(
    scores[ok] <= qs[[1L]], "LOW",
    ifelse(scores[ok] >= qs[[2L]], "HIGH", "MID")
  )
  out
}

build_expression_quantile_splits_matrix <- function(mat, low_q = 1/3, high_q = 2/3) {
  genes <- rownames(mat)
  patients <- colnames(mat)
  splits <- matrix(NA_character_, nrow = length(genes), ncol = length(patients),
                   dimnames = list(genes, patients))

  for (i in seq_along(genes)) {
    splits[i, ] <- assign_quantile_groups(mat[i, ], low_q = low_q, high_q = high_q)
  }
  splits
}

load_gene_patient_matrix <- function(path) {
  x <- readRDS(path)
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
  mat
}

resolve_pseudobulk_path <- function(dataset_name, input_dir) {
  matrix_name <- paste0(dataset_name, PSEUDOBULK_MATRIX_SUFFIX)
  in_candidates <- c(
    file.path(input_dir, dataset_name, matrix_name),
    file.path(input_dir, matrix_name),
    file.path(input_dir, dataset_name, paste0(dataset_name, "_gene_by_patient_means_exp.rds")),
    file.path(input_dir, paste0(dataset_name, "_gene_by_patient_means_exp.rds")),
    file.path(input_dir, dataset_name, paste0(dataset_name, "_gene_by_patient_mean_expr.rds")),
    file.path(input_dir, paste0(dataset_name, "_gene_by_patient_mean_expr.rds"))
  )
  in_path <- in_candidates[file.exists(in_candidates)][1]
  if (is.na(in_path) || is.null(in_path)) {
    stop(
      "Could not find pseudobulk matrix for '", dataset_name, "'. Tried:\n",
      paste(" -", in_candidates, collapse = "\n")
    )
  }
  in_path
}

compute_rank_matrix <- function(mat) {
  if (requireNamespace("matrixStats", quietly = TRUE) &&
      exists("colRanks", where = asNamespace("matrixStats"), inherits = FALSE)) {
    return(matrixStats::colRanks(
      mat,
      ties.method = "average",
      preserveShape = TRUE,
      na.last = "keep"
    ))
  }
  apply(mat, 2, function(v) rank(v, ties.method = "average", na.last = "keep"))
}

build_patient_splits_matrix <- function(mat, rank_mat = NULL) {
  if (is.null(rank_mat)) {
    rank_mat <- compute_rank_matrix(mat)
  }
  stopifnot(identical(dim(rank_mat), dim(mat)))

  genes <- rownames(mat)
  patients <- colnames(mat)
  splits <- matrix(NA_character_, nrow = length(genes), ncol = length(patients),
                   dimnames = list(genes, patients))

  for (i in seq_along(genes)) {
    splits[i, ] <- assign_tertiles(rank_mat[i, ])
  }
  splits
}

compute_patient_pcs <- function(mat, n_pc = 2L) {
  n_pc <- as.integer(n_pc)
  if (!n_pc %in% c(1L, 2L)) {
    stop("n_pc must be 1 or 2. Got: ", n_pc)
  }

  patients <- colnames(mat)
  t_mat <- t(mat)
  complete <- stats::complete.cases(t_mat)
  n_complete <- sum(complete)

  if (n_complete < n_pc + 2L) {
    stop(
      "Need at least ", n_pc + 2L,
      " patients with complete pseudobulk rows for PCA; got ", n_complete
    )
  }

  t_complete <- t_mat[complete, , drop = FALSE]
  gene_sd <- apply(t_complete, 2, stats::sd, na.rm = TRUE)
  variable_genes <- is.finite(gene_sd) & gene_sd > 0
  n_variable <- sum(variable_genes)
  if (n_variable < n_pc + 2L) {
    stop(
      "Need at least ", n_pc + 2L,
      " genes with nonzero variance across complete patients for PCA; got ",
      n_variable
    )
  }
  if (n_variable < ncol(t_complete)) {
    message(
      "PCA: excluded ", ncol(t_complete) - n_variable,
      " zero-variance gene(s) across complete patients."
    )
  }

  pc_fit <- stats::prcomp(
    t_complete[, variable_genes, drop = FALSE],
    center = TRUE,
    scale. = TRUE
  )
  n_pc_use <- min(n_pc, ncol(pc_fit$x), nrow(pc_fit$x))

  pcs <- matrix(NA_real_, nrow = length(patients), ncol = n_pc_use,
                dimnames = list(patients, paste0("PC", seq_len(n_pc_use))))
  pc_cols <- seq_len(n_pc_use)
  pcs[complete, pc_cols] <- pc_fit$x[, pc_cols, drop = FALSE]

  if (n_complete < length(patients)) {
    message(
      "PCA used ", n_complete, " of ", length(patients),
      " patients (complete cases only)."
    )
  }

  pcs
}

compute_residual_matrix <- function(mat, pcs) {
  genes <- rownames(mat)
  patients <- colnames(mat)
  resid <- matrix(NA_real_, nrow = length(genes), ncol = length(patients),
                  dimnames = list(genes, patients))

  if (is.null(rownames(pcs))) rownames(pcs) <- patients
  if (!identical(rownames(pcs), patients)) {
    stop("pcs rownames must match colnames(mat).")
  }

  pc_ok <- stats::complete.cases(pcs)
  if (sum(pc_ok) < 3L) {
    return(resid)
  }

  P <- cbind(Intercept = 1, pcs[pc_ok, , drop = FALSE])
  Q <- qr.Q(qr(P))
  Y <- mat[, pc_ok, drop = FALSE]

  gene_complete <- apply(Y, 1, function(r) all(is.finite(r)))
  if (any(gene_complete)) {
    Yc <- Y[gene_complete, , drop = FALSE]
    Yt <- t(Yc)
    Rc <- Yt - Q %*% crossprod(Q, Yt)
    resid[which(gene_complete), pc_ok] <- t(Rc)
  }

  incomplete <- which(!gene_complete)
  if (length(incomplete) > 0L) {
    pcs_ok <- pcs[pc_ok, , drop = FALSE]
    for (ii in incomplete) {
      y <- mat[ii, pc_ok, drop = FALSE]
      ok <- is.finite(y)
      if (sum(ok) < 3L) next
      fit <- stats::lm(y[ok] ~ pcs_ok[ok, , drop = FALSE])
      r <- stats::residuals(fit)
      resid[ii, pc_ok][ok] <- r
    }
  }

  resid
}

build_residual_splits_matrix <- function(mat, n_pc = 2L, pcs = NULL) {
  if (is.null(pcs)) {
    pcs <- compute_patient_pcs(mat, n_pc = n_pc)
  }

  resid_mat <- compute_residual_matrix(mat, pcs)
  stopifnot(identical(dim(resid_mat), dim(mat)))

  genes <- rownames(mat)
  patients <- colnames(mat)
  splits <- matrix(NA_character_, nrow = length(genes), ncol = length(patients),
                   dimnames = list(genes, patients))

  for (i in seq_along(genes)) {
    splits[i, ] <- assign_tertiles(resid_mat[i, ])
  }

  list(splits = splits, residuals = resid_mat, pcs = pcs)
}

splits_to_integer_matrix <- function(splits) {
  L <- matrix(0L, nrow = nrow(splits), ncol = ncol(splits),
              dimnames = dimnames(splits))
  for (lev in seq_along(LABEL_LEVELS)) {
    L[splits == LABEL_LEVELS[[lev]]] <- as.integer(lev)
  }
  L
}

gene_indices_to_write <- function(all_genes, gene_list_path = NULL) {
  if (is.null(gene_list_path) || !nzchar(gene_list_path)) {
    return(seq_along(all_genes))
  }
  wanted <- load_gene_list(gene_list_path)
  idx <- match(wanted, all_genes)
  missing <- wanted[is.na(idx)]
  if (length(missing) > 0L) {
    warning(
      length(missing), " gene(s) from --gene_list not in pseudobulk matrix: ",
      paste(head(missing, 15L), collapse = ", "),
      if (length(missing) > 15L) " ..." else "",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) {
    stop("No genes from --gene_list found in pseudobulk matrix.")
  }
  unique(idx)
}

load_gene_list <- function(path = NULL) {
  if (is.null(path) || !nzchar(path)) {
    return(SCDIFFCOM_GENE_PANEL)
  }
  path <- path.expand(path)
  if (!file.exists(path)) {
    stop("Gene list file not found: ", path)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds" || ext == "rda") {
    g <- readRDS(path)
  } else {
    g <- scan(path, what = character(), quiet = TRUE)
  }
  if (is.data.frame(g)) {
    g <- g[[1]]
  }
  unique(as.character(g))
}

load_splits_from_grouped_dir <- function(base_dir, dataset_name, genes = NULL) {
  ds_dir <- file.path(path.expand(base_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    return(NULL)
  }

  all_splits_path <- file.path(ds_dir, paste0(dataset_name, ALL_SPLITS_SUFFIX))
  if (is.null(genes) && file.exists(all_splits_path)) {
    return(readRDS(all_splits_path))
  }

  if (is.null(genes)) {
    files <- list.files(ds_dir, pattern = paste0("_", dataset_name, "_grouped\\.rds$"),
                        full.names = TRUE)
    if (length(files) == 0L) {
      return(NULL)
    }
    genes <- sub(paste0("_", dataset_name, "_grouped\\.rds$"), "",
                 basename(files))
  }

  patients <- NULL
  split_rows <- list()

  for (gene in genes) {
    f <- file.path(ds_dir, paste0(gene, "_", dataset_name, "_grouped.rds"))
    if (!file.exists(f)) next
    df <- readRDS(f)
    exp_col <- grep("_EXP$", colnames(df), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(exp_col)) next
    if (is.null(patients)) {
      patients <- as.character(df$patient_id)
    }
    vec <- setNames(as.character(df[[exp_col]]), as.character(df$patient_id))
    split_rows[[gene]] <- vec[patients]
  }

  if (length(split_rows) == 0L) {
    return(NULL)
  }

  genes_found <- names(split_rows)
  splits <- matrix(NA_character_, nrow = length(genes_found), ncol = length(patients),
                   dimnames = list(genes_found, patients))
  for (g in genes_found) {
    splits[g, ] <- split_rows[[g]]
  }
  splits
}

load_one_grouped_gene <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(df) || !is.data.frame(df)) {
    return(NULL)
  }
  exp_col <- grep("_EXP$", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(exp_col) || !"patient_id" %in% colnames(df)) {
    return(NULL)
  }
  patients <- as.character(df$patient_id)
  mean_expr <- if ("mean_expr" %in% colnames(df)) {
    setNames(as.numeric(df$mean_expr), patients)
  } else {
    NULL
  }
  splits <- setNames(as.character(df[[exp_col]]), patients)
  list(patients = patients, mean_expr = mean_expr, splits = splits)
}

load_splits_from_rankgenes_dir <- function(rankgenes_dir, dataset_name, genes = NULL) {
  ds_dir <- file.path(path.expand(rankgenes_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    stop("RankGenes directory not found: ", ds_dir)
  }

  splits <- load_splits_from_grouped_dir(rankgenes_dir, dataset_name, genes = genes)
  if (is.null(splits)) {
    stop("No splits loaded from ", ds_dir)
  }
  splits
}

load_splits_from_residual_dir <- function(residual_dir, dataset_name, genes = NULL) {
  ds_dir <- file.path(path.expand(residual_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    stop("Residual splits directory not found: ", ds_dir)
  }

  all_splits_path <- file.path(ds_dir, paste0(dataset_name, ALL_SPLITS_SUFFIX))
  if (is.null(genes) && file.exists(all_splits_path)) {
    return(readRDS(all_splits_path))
  }

  if (is.null(genes)) {
    files <- list.files(ds_dir, pattern = paste0("_", dataset_name, "_grouped\\.rds$"),
                        full.names = TRUE)
    if (length(files) == 0L) {
      stop("No grouped .rds files in ", ds_dir)
    }
    genes <- sub(paste0("_", dataset_name, "_grouped\\.rds$"), "",
                 basename(files))
  }

  patients <- NULL
  split_rows <- list()

  for (gene in genes) {
    f <- file.path(ds_dir, paste0(gene, "_", dataset_name, "_grouped.rds"))
    if (!file.exists(f)) next
    df <- readRDS(f)
    exp_col <- grep("_EXP$", colnames(df), value = TRUE)[1]
    if (is.na(exp_col)) next
    if (is.null(patients)) {
      patients <- as.character(df$patient_id)
    }
    vec <- setNames(as.character(df[[exp_col]]), as.character(df$patient_id))
    split_rows[[gene]] <- vec[patients]
  }

  if (length(split_rows) == 0L) {
    stop("No splits loaded from ", ds_dir)
  }

  genes_found <- names(split_rows)
  splits <- matrix(NA_character_, nrow = length(genes_found), ncol = length(patients),
                   dimnames = list(genes_found, patients))
  for (g in genes_found) {
    splits[g, ] <- split_rows[[g]]
  }
  splits
}

compute_jaccard_similarity_matrix <- function(L) {
  G <- nrow(L)
  J <- matrix(NA_real_, G, G, dimnames = list(rownames(L), rownames(L)))
  diag(J) <- 1

  for (i in seq_len(G - 1L)) {
    ai <- L[i, ]
    for (j in (i + 1L):G) {
      aj <- L[j, ]
      valid <- (ai > 0L) & (aj > 0L)
      union <- (ai > 0L) | (aj > 0L)
      n_union <- sum(union)
      if (n_union == 0L) {
        sim <- NA_real_
      } else {
        sim <- sum(valid & (ai == aj)) / n_union
      }
      J[i, j] <- J[j, i] <- sim
    }
  }
  J
}

find_duplicate_split_groups <- function(splits, threshold = 1) {
  L <- splits_to_integer_matrix(splits)
  keys <- apply(L, 1, paste, collapse = ",")
  grp <- split(rownames(splits), keys)
  grp <- grp[lengths(grp) >= 1L]

  if (threshold >= 1) {
    dup <- grp[lengths(grp) > 1L]
    return(dup)
  }

  # Near-duplicates: cluster by Jaccard >= threshold (greedy on representatives)
  genes <- rownames(splits)
  J <- compute_jaccard_similarity_matrix(L)
  used <- rep(FALSE, length(genes))
  names(used) <- genes
  clusters <- list()

  for (g in genes) {
    if (used[[g]]) next
    members <- genes[!used & J[g, ] >= threshold]
    members <- members[!is.na(members)]
    if (length(members) == 0L) members <- g
    used[members] <- TRUE
    clusters[[length(clusters) + 1L]] <- members
  }
  clusters
}

duplicate_groups_to_tsv <- function(duplicate_groups) {
  rows <- lapply(seq_along(duplicate_groups), function(i) {
    members <- duplicate_groups[[i]]
    data.frame(
      group_id = i,
      representative = members[[1]],
      n_genes = length(members),
      members = paste(members, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

top_similar_pairs <- function(J, n = 50L) {
  genes <- rownames(J)
  G <- length(genes)
  if (G < 2L) return(data.frame())
  pairs <- list()
  k <- 1L
  for (i in seq_len(G - 1L)) {
    for (j in (i + 1L):G) {
      pairs[[k]] <- c(genes[[i]], genes[[j]], J[i, j])
      k <- k + 1L
    }
  }
  df <- as.data.frame(do.call(rbind, pairs), stringsAsFactors = FALSE)
  colnames(df) <- c("gene_a", "gene_b", "jaccard_similarity")
  df$jaccard_similarity <- as.numeric(df$jaccard_similarity)
  df <- df[order(-df$jaccard_similarity), ]
  head(df, n)
}
