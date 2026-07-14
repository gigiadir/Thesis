# Shared helpers for reproducibility atlas validation QC.

STAR_CHECKS <- c(
  "eff_n", "batch_collapse", "diag_visible", "Rg_vs_n", "raw_gap",
  "null_centered", "na_exchangeable", "global_test", "loco_stability",
  "fastpath_equiv"
)

init_verdict_file <- function(validation_dir) {
  path <- file.path(validation_dir, "results", "VERDICT.tsv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    writeLines(
      "check_id\tstage\tstatus\tvalue\texpected\tnote",
      path
    )
  }
  invisible(path)
}

append_verdict <- function(validation_dir, check_id, stage, status, value = NA,
                           expected = "", note = "") {
  path <- file.path(validation_dir, "results", "VERDICT.tsv")
  val_str <- if (is.na(value)) "" else as.character(value)
  line <- paste(
    check_id, stage, status, val_str, expected, note,
    sep = "\t"
  )
  cat(line, "\n", file = path, append = TRUE)
  invisible(line)
}

read_verdict <- function(validation_dir) {
  path <- file.path(validation_dir, "results", "VERDICT.tsv")
  if (!file.exists(path)) {
    return(data.frame(
      check_id = character(), stage = character(), status = character(),
      value = character(), expected = character(), note = character(),
      stringsAsFactors = FALSE
    ))
  }
  readr::read_tsv(path, show_col_types = FALSE)
}

parse_cci <- function(J) {
  m <- regexec("^(.+?)_(.+?)_(.+)$", J)
  parts <- regmatches(J, m)
  data.frame(
    CCI = J,
    EMITTER_CELLTYPE = vapply(parts, function(x) if (length(x) >= 4) x[2] else NA_character_, character(1)),
    RECEIVER_CELLTYPE = vapply(parts, function(x) if (length(x) >= 4) x[3] else NA_character_, character(1)),
    LRI = vapply(parts, function(x) if (length(x) >= 4) x[4] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

safe_read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_tsv(path, show_col_types = FALSE)
}

load_validation_context <- function(output_dir,
                                    nulls_dir = NULL,
                                    atlas_dir = NULL) {
  output_dir <- path.expand(output_dir)
  if (is.null(nulls_dir)) {
    batch_cfg_path <- file.path(atlas_dir, "config_batch.yml")
    if (file.exists(batch_cfg_path)) {
      batch_cfg <- yaml::read_yaml(batch_cfg_path)
      nulls_dir <- path.expand(batch_cfg$output_dir)
    } else {
      nulls_dir <- output_dir
    }
  } else {
    nulls_dir <- path.expand(nulls_dir)
  }

  cfg_path <- file.path(atlas_dir, "config.yml")
  cfg <- yaml::read_yaml(cfg_path)
  cfg$output_dir <- output_dir
  if (!is.null(cfg$seed)) set.seed(cfg$seed)

  results_dir <- file.path(output_dir, "results")
  nulls_results_dir <- file.path(nulls_dir, "results")

  read_or_stop <- function(fname, label) {
    path <- file.path(results_dir, fname)
    if (!file.exists(path)) stop("Missing ", label, ": ", path)
    if (grepl("\\.tsv$", fname)) safe_read_tsv(path) else safe_read_rds(path)
  }

  X <- read_or_stop("X.rds", "tensor X")
  Xtilde <- read_or_stop("Xtilde.rds", "centered tensor Xtilde")
  J <- read_or_stop("J.rds", "vocabulary J")
  na_mask <- read_or_stop("na_mask.rds", "na_mask")
  repro_df <- read_or_stop("repro_scores.tsv", "repro_scores")
  cor_mats <- read_or_stop("cor_pair_matrices.rds", "cor_pair_matrices")
  gene_universe <- read_or_stop("gene_universe.tsv", "gene_universe")$gene

  stage04_path <- file.path(results_dir, "stage04_atlas_env.rds")
  stage04_env <- safe_read_rds(stage04_path)
  if (is.null(stage04_env)) {
    stop("Missing stage04_atlas_env.rds — required for dup_collapse and repro internals")
  }

  cohorts <- stage04_env$cohorts
  cohort_pairs_idx <- stage04_env$cohort_pairs
  pair_labels <- stage04_env$pair_labels
  repro <- stage04_env$repro
  malignant_by_cohort <- stage04_env$malignant_by_cohort
  cci_by_cohort <- stage04_env$cci_by_cohort

  stage05_env <- safe_read_rds(file.path(nulls_results_dir, "stage05_atlas_env.rds"))
  null_repro <- safe_read_rds(file.path(nulls_results_dir, "null_reproscore_matrix.rds"))
  repro_with_nulls <- safe_read_tsv(file.path(nulls_results_dir, "repro_scores_with_nulls.tsv"))
  global_null_txt <- file.path(nulls_results_dir, "global_null.txt")
  null_checkpoint <- safe_read_rds(file.path(nulls_results_dir, "null_perm_checkpoint.rds"))

  if (!is.null(null_checkpoint) && (is.null(null_repro) || ncol(null_repro) < cfg$n_perm)) {
    if (null_checkpoint$completed >= cfg$n_perm) {
      # New schema stores null_Rg; fall back to the legacy field name.
      null_repro <- if (!is.null(null_checkpoint$null_Rg)) {
        null_checkpoint$null_Rg
      } else {
        null_checkpoint$null_repro_scores
      }
    }
  }

  stage5_complete <- !is.null(stage05_env) ||
    (!is.null(null_repro) && ncol(null_repro) >= cfg$n_perm)

  list(
    cfg = cfg,
    output_dir = output_dir,
    nulls_dir = nulls_dir,
    results_dir = results_dir,
    nulls_results_dir = nulls_results_dir,
    atlas_dir = atlas_dir,
    X = X,
    Xtilde = Xtilde,
    J = J,
    na_mask = na_mask,
    repro_df = repro_df,
    cor_mats = cor_mats,
    gene_universe = gene_universe,
    cohorts = cohorts,
    cohort_pairs_idx = cohort_pairs_idx,
    pair_labels = pair_labels,
    repro = repro,
    malignant_by_cohort = malignant_by_cohort,
    cci_by_cohort = cci_by_cohort,
    stage04_env = stage04_env,
    stage05_env = stage05_env,
    null_repro = null_repro,
    repro_with_nulls = repro_with_nulls,
    global_null_txt = global_null_txt,
    null_checkpoint = null_checkpoint,
    stage5_complete = stage5_complete,
    eff_n_median = NA_real_
  )
}

compute_eff_n <- function(X, cohorts, cohort_pairs_idx, gene_universe) {
  n_pairs <- length(cohort_pairs_idx)
  pair_names <- vapply(cohort_pairs_idx, function(pr) {
    paste(cohorts[pr[1]], cohorts[pr[2]], sep = "_vs_")
  }, character(1))
  eff_mat <- matrix(NA_integer_, length(gene_universe), n_pairs,
                    dimnames = list(gene_universe, pair_names))
  for (p in seq_len(n_pairs)) {
    a <- cohort_pairs_idx[[p]][1]
    b <- cohort_pairs_idx[[p]][2]
    ca <- cohorts[a]
    cb <- cohorts[b]
    for (g in seq_along(gene_universe)) {
      x <- X[[ca]][g, ]
      y <- X[[cb]][g, ]
      eff_mat[g, p] <- sum(is.finite(x) & is.finite(y))
    }
  }
  eff_mat
}

stack_cohort_matrix <- function(X, cohorts) {
  mats <- lapply(cohorts, function(ds) X[[ds]])
  combined <- do.call(rbind, mats)
  cohort_label <- rep(cohorts, each = nrow(X[[1]]))
  list(mat = combined, cohort_label = cohort_label)
}

pca_silhouette <- function(mat, cohort_label, n_pc = 10) {
  finite_rows <- rowSums(is.finite(mat)) > ncol(mat) * 0.1
  m <- mat[finite_rows, , drop = FALSE]
  cl <- cohort_label[finite_rows]
  m[is.na(m)] <- 0
  if (nrow(m) < 10 || length(unique(cl)) < 2) return(NA_real_)
  pca <- prcomp(m, center = TRUE, scale. = TRUE)
  k <- min(n_pc, ncol(pca$x))
  d <- dist(pca$x[, seq_len(k), drop = FALSE])
  sil <- cluster::silhouette(as.integer(factor(cl)), d)
  mean(sil[, "sil_width"])
}

shuffle_cci_within_gene <- function(X) {
  lapply(X, function(m) {
    out <- m
    for (g in seq_len(nrow(m))) {
      row <- m[g, ]
      finite <- is.finite(row)
      if (sum(finite) > 1) {
        out[g, finite] <- sample(row[finite])
      }
    }
    out
  })
}

report_gate <- function(validation_dir, gate_name, checks = NULL) {
  verdict <- read_verdict(validation_dir)
  if (!is.null(checks)) {
    verdict <- verdict[verdict$check_id %in% checks, , drop = FALSE]
  }
  message("\n=== ", gate_name, " ===")
  if (nrow(verdict) == 0) {
    message("(no verdict rows)")
    return(invisible(verdict))
  }
  for (i in seq_len(nrow(verdict))) {
    row <- verdict[i, ]
    message(sprintf("  [%s] %s: %s — %s (value=%s)",
                    row$status, row$check_id, row$note, row$expected, row$value))
  }
  invisible(verdict)
}

gate_a_should_stop <- function(validation_dir) {
  verdict <- read_verdict(validation_dir)
  eff <- verdict[verdict$check_id == "eff_n", , drop = FALSE]
  nrow(eff) > 0 && eff$status[1] == "FAIL"
}

gate_b_decision <- function(validation_dir) {
  verdict <- read_verdict(validation_dir)
  get_status <- function(id) {
    rows <- verdict[verdict$check_id == id, , drop = FALSE]
    if (nrow(rows) == 0) return(NA_character_)
    rows$status[nrow(rows)]
  }
  list(
    null_centered = get_status("null_centered"),
    raw_gap = get_status("raw_gap"),
    diag_visible = get_status("diag_visible"),
    Rg_vs_n = get_status("Rg_vs_n"),
    global_test = get_status("global_test")
  )
}
