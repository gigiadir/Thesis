# Stage 2 — extract logFC tensor X: list of 4 matrices [G x J].

.gene_logfc_vector <- function(df, J) {
  vec <- setNames(rep(NA_real_, length(J)), J)
  if (nrow(df) == 0) return(vec)
  agg <- df %>%
    group_by(CCI) %>%
    summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop")
  present <- agg$CCI
  vec[present] <- agg$LOGFC
  vec
}

.assert_tensor_alignment <- function(X, gene_universe, J, cohorts) {
  stopifnot(length(X) == length(cohorts))
  for (i in seq_along(cohorts)) {
    stopifnot(identical(rownames(X[[i]]), gene_universe))
    stopifnot(identical(colnames(X[[i]]), J))
    stopifnot(nrow(X[[i]]) == length(gene_universe))
    stopifnot(ncol(X[[i]]) == length(J))
  }
  if (i > 1) {
    for (i in seq(2, length(cohorts))) {
      stopifnot(identical(rownames(X[[1]]), rownames(X[[i]])))
      stopifnot(identical(colnames(X[[1]]), colnames(X[[i]])))
    }
  }
}

run_stage_02_tensor <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  cohorts <- atlas_env$cohorts
  gene_universe <- atlas_env$gene_universe
  malignant_by_cohort <- atlas_env$malignant_by_cohort
  J <- atlas_env$J
  na_mask <- atlas_env$na_mask

  message("Stage 2: build logFC tensor X [G x J]")

  X <- lapply(cohorts, function(ds) {
    gene_list <- malignant_by_cohort[[ds]]
    mat <- vapply(
      gene_universe,
      function(g) .gene_logfc_vector(gene_list[[g]], J),
      numeric(length(J))
    )
    mat <- t(mat)
    rownames(mat) <- gene_universe
    colnames(mat) <- J
    mat
  })
  names(X) <- cohorts

  .assert_tensor_alignment(X, gene_universe, J, cohorts)

  na_frac <- vapply(X, function(m) mean(is.na(m)), numeric(1))
  na_tbl <- data.frame(
    cohort = cohorts,
    na_fraction = na_frac,
    n_genes = nrow(X[[1]]),
    n_ccis = ncol(X[[1]]),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(na_tbl, file.path(output_dir, "results", "tensor_na_fraction.tsv"))
  for (i in seq_along(cohorts)) {
    message(sprintf("  %s: NA fraction = %.3f", cohorts[[i]], na_frac[[i]]))
  }

  saveRDS(X, file.path(output_dir, "results", "X.rds"))
  saveRDS(na_mask, file.path(output_dir, "results", "na_mask.rds"))

  atlas_env$X <- X
  saveRDS(atlas_env, file.path(output_dir, "results", "stage02_atlas_env.rds"))
  message("  Saved X.rds")

  invisible(atlas_env)
}
