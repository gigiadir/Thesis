# Shared helpers for RankGenes and split-similarity analysis.
# Sourced by scDiffCom-Preprocess-RankGenes.R and analyzeGeneSplitJaccard.R.

PSEUDOBULK_MATRIX_SUFFIX <- "_pseudobulk_matrix.rds"
ALL_SPLITS_SUFFIX <- "_all_patient_splits.rds"

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
  in_candidates <- c(
    file.path(input_dir, paste0(dataset_name, PSEUDOBULK_MATRIX_SUFFIX)),
    file.path(input_dir, paste0(dataset_name, "_gene_by_patient_means_exp.rds")),
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

splits_to_integer_matrix <- function(splits) {
  L <- matrix(0L, nrow = nrow(splits), ncol = ncol(splits),
              dimnames = dimnames(splits))
  for (lev in seq_along(LABEL_LEVELS)) {
    L[splits == LABEL_LEVELS[[lev]]] <- as.integer(lev)
  }
  L
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

load_splits_from_rankgenes_dir <- function(rankgenes_dir, dataset_name, genes = NULL) {
  ds_dir <- file.path(path.expand(rankgenes_dir), dataset_name)
  if (!dir.exists(ds_dir)) {
    stop("RankGenes directory not found: ", ds_dir)
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
