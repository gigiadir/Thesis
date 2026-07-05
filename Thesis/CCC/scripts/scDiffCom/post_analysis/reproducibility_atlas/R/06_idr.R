# Stage 6 — pairwise IDR on DE CCI rankings per gene.

.ensure_idr <- function() {
  if (!requireNamespace("idr", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      stop("Package 'idr' required. Install with: BiocManager::install('idr')")
    }
    BiocManager::install("idr", ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(idr))
}

.gene_rank_values <- function(df, ccis, ranking = "signed") {
  vec <- setNames(rep(NA_real_, length(ccis)), ccis)
  if (nrow(df) == 0) return(vec)
  agg <- df %>%
    group_by(CCI) %>%
    summarise(LOGFC = mean(LOGFC, na.rm = TRUE), .groups = "drop")
  vals <- agg$LOGFC
  if (ranking == "magnitude") vals <- abs(vals)
  vec[agg$CCI] <- vals
  vec
}

.run_pairwise_idr <- function(va, vb, threshold = 0.05) {
  finite <- is.finite(va) & is.finite(vb)
  if (sum(finite) < 10) {
    return(list(pass = character(0), local_idr = numeric(0), ccis = character(0)))
  }
  ccis <- names(va)[finite]
  x <- va[finite]
  y <- vb[finite]
  fit <- tryCatch(
    idr::idr(x, y, mu = 0.5, sigma = 0.1, p = threshold),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(pass = character(0), local_idr = rep(NA_real_, length(ccis)), ccis = ccis))
  }
  local <- fit$idr[, 1]
  names(local) <- ccis
  pass <- names(local)[local <= threshold]
  list(pass = pass, local_idr = local, ccis = ccis)
}

run_stage_06_idr <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  cohorts <- atlas_env$cohorts
  gene_universe <- atlas_env$gene_universe
  malignant_by_cohort <- atlas_env$malignant_by_cohort
  idx_pairs <- atlas_env$cohort_pairs
  pair_labels <- atlas_env$pair_labels
  ranking <- cfg$idr_ranking
  threshold <- cfg$idr_threshold
  n_pairs <- length(idx_pairs)

  message("Stage 6: IDR (ranking = ", ranking, ")")
  .ensure_idr()

  idr_long <- list()
  idr_summary <- vector("list", length(gene_universe))
  names(idr_summary) <- gene_universe

  for (g in gene_universe) {
    ccis_g <- sort(unique(unlist(lapply(cohorts, function(ds) {
      malignant_by_cohort[[ds]][[g]]$CCI
    }))))
    if (length(ccis_g) < 10) {
      idr_summary[[g]] <- list(
        idr_pass_fraction = NA_real_,
        idr_pass_CCIs = character(0)
      )
      next
    }

    vecs <- lapply(cohorts, function(ds) {
      .gene_rank_values(malignant_by_cohort[[ds]][[g]], ccis = ccis_g, ranking = ranking)
    })
    names(vecs) <- cohorts

    pair_pass <- vector("list", n_pairs)
    names(pair_pass) <- pair_labels
    cci_pair_pass <- setNames(
      rep(0L, length(ccis_g)),
      ccis_g
    )

    for (p in seq_len(n_pairs)) {
      a <- idx_pairs[[p]][1]
      b <- idx_pairs[[p]][2]
      res <- .run_pairwise_idr(vecs[[a]], vecs[[b]], threshold = threshold)
      pair_pass[[p]] <- res$pass
      if (length(res$pass) > 0) {
        cci_pair_pass[res$pass] <- cci_pair_pass[res$pass] + 1L
      }
      if (length(res$ccis) > 0) {
        idr_long[[length(idr_long) + 1]] <- data.frame(
          gene = g,
          pair = pair_labels[[p]],
          CCI = res$ccis,
          local_idr = as.numeric(res$local_idr[res$ccis]),
          pass = res$ccis %in% res$pass,
          stringsAsFactors = FALSE
        )
      }
    }

    pass_frac_by_cci <- cci_pair_pass / n_pairs
    idr_pass_CCIs <- names(pass_frac_by_cci)[pass_frac_by_cci >= 0.5]
    idr_pass_fraction <- if (length(ccis_g) > 0) {
      mean(pass_frac_by_cci >= 0.5)
    } else {
      NA_real_
    }

    idr_summary[[g]] <- list(
      idr_pass_fraction = idr_pass_fraction,
      idr_pass_CCIs = idr_pass_CCIs
    )
  }

  idr_df <- data.frame(
    gene = gene_universe,
    idr_pass_fraction = vapply(idr_summary, function(x) x$idr_pass_fraction, numeric(1)),
    idr_pass_CCIs = vapply(idr_summary, function(x) paste(x$idr_pass_CCIs, collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )

  idr_long_df <- if (length(idr_long) > 0) {
    dplyr::bind_rows(idr_long)
  } else {
    data.frame(
      gene = character(), pair = character(), CCI = character(),
      local_idr = numeric(), pass = logical()
    )
  }

  readr::write_tsv(idr_df, file.path(output_dir, "results", "idr_summary.tsv"))
  readr::write_csv(idr_long_df, file.path(output_dir, "results", "atlas_cci_long.csv"))

  atlas_env$idr_df <- idr_df
  atlas_env$idr_long_df <- idr_long_df
  saveRDS(atlas_env, file.path(output_dir, "results", "stage06_atlas_env.rds"))
  message(sprintf("  IDR complete for %d genes", nrow(idr_df)))

  invisible(atlas_env)
}
