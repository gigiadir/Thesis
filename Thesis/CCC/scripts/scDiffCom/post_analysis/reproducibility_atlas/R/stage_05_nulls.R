# Stage 5 — cross-gene pairing null + EVT/GPD calibration of R_g.
#
# Null model (advisor Step 5): for each gene g and cohort pair (c1, c2), pair
# g's cohort-c1 vector against a RANDOM different gene g' in cohort c2 and take
# the same median-across-pairs collapse used for the observed R_g. Because
# C[g, g'] (the off-diagonal of the Stage 4 gene x gene Spearman matrix) already
# IS Spearman(v_g^{c1}, v_{g'}^{c2}), each null draw is a cheap resampling of
# off-diagonal entries — no correlation recomputation required.
#
# Calibration (advisor Step 6, Knijnenburg 2009): per-gene hybrid empirical/GPD
# p-value (see stage_05_evt.R). The identical treatment is applied to the global
# statistic (mean R_g across genes).

resolve_n_cores <- function(requested = NULL, max_cores = 50L) {
  max_cores <- as.integer(max_cores)
  if (length(max_cores) != 1L || !is.finite(max_cores) || max_cores < 1L) {
    stop("max_cores must be a positive integer")
  }

  if (is.null(requested)) {
    env_cores <- Sys.getenv("MC_CORES", unset = NA_character_)
    if (!is.na(env_cores) && nzchar(env_cores)) {
      requested <- suppressWarnings(as.integer(env_cores))
    }
  }

  if (is.null(requested) || length(requested) != 1L || !is.finite(requested)) {
    requested <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  }

  as.integer(max(1L, min(max_cores, requested)))
}

.null_perm_checkpoint_path <- function(output_dir) {
  file.path(output_dir, "results", "null_perm_checkpoint.rds")
}

load_null_perm_checkpoint <- function(output_dir, n_perm, resume = TRUE) {
  if (!isTRUE(resume)) {
    return(NULL)
  }

  path <- .null_perm_checkpoint_path(output_dir)
  if (!file.exists(path)) {
    return(NULL)
  }

  ckpt <- readRDS(path)
  if (!identical(as.integer(ckpt$n_perm), as.integer(n_perm))) {
    message("Ignoring null_perm_checkpoint.rds (n_perm mismatch: ",
            ckpt$n_perm, " vs ", n_perm, ")")
    return(NULL)
  }

  ckpt
}

save_null_perm_checkpoint <- function(output_dir,
                                      null_Rg,
                                      null_global,
                                      completed,
                                      n_perm,
                                      seed) {
  path <- .null_perm_checkpoint_path(output_dir)
  saveRDS(
    list(
      completed = as.integer(completed),
      null_Rg = null_Rg,
      null_global = null_global,
      n_perm = as.integer(n_perm),
      seed = seed,
      timestamp = Sys.time()
    ),
    path
  )
  invisible(path)
}

delete_null_perm_checkpoint <- function(output_dir) {
  path <- .null_perm_checkpoint_path(output_dir)
  if (file.exists(path)) {
    file.remove(path)
    message("Removed ", basename(path))
  }
  invisible(path)
}

# Generate the cross-gene null R_g matrix (G x n_perm), resumable + parallel.
run_cross_gene_null <- function(cor_mats,
                                gene_universe,
                                n_perm,
                                n_cores,
                                aggregate = "median",
                                resume = TRUE,
                                checkpoint_every = 50L,
                                output_dir,
                                seed = NULL) {
  n_perm <- as.integer(n_perm)
  checkpoint_every <- max(1L, as.integer(checkpoint_every))
  n_genes <- length(gene_universe)

  null_Rg <- matrix(NA_real_, n_genes, n_perm)
  null_global <- rep(NA_real_, n_perm)
  completed <- 0L

  ckpt <- load_null_perm_checkpoint(output_dir, n_perm, resume = resume)
  if (!is.null(ckpt)) {
    completed <- min(ckpt$completed, n_perm)
    null_Rg[, seq_len(completed)] <- ckpt$null_Rg[, seq_len(completed)]
    null_global[seq_len(completed)] <- ckpt$null_global[seq_len(completed)]
    message(sprintf(
      "Resuming from permutation %d / %d (checkpoint %s)",
      completed + 1L, n_perm, ckpt$timestamp
    ))
  }

  if (completed >= n_perm) {
    return(list(null_Rg = null_Rg, null_global = null_global, completed = completed))
  }

  if (!is.null(seed)) set.seed(seed)

  run_one_perm <- function(b) {
    v <- null_one_perm_Rg(cor_mats, aggregate = aggregate)
    c(v, mean(v, na.rm = TRUE))
  }

  remaining <- seq.int(completed + 1L, n_perm)
  use_parallel <- .Platform$OS.type == "unix" && n_cores > 1L && length(remaining) >= 50L

  if (use_parallel) {
    message(sprintf(
      "  parallel::mclapply (%d cores, %d perms remaining, chunk size %d)",
      n_cores, length(remaining), checkpoint_every
    ))
    chunk_starts <- seq(remaining[[1]], n_perm, by = checkpoint_every)
    for (chunk_start in chunk_starts) {
      chunk_end <- min(chunk_start + checkpoint_every - 1L, n_perm)
      chunk <- seq.int(chunk_start, chunk_end)
      perm_results <- parallel::mclapply(
        chunk,
        function(b) {
          if (b == chunk_start) message(sprintf("  permutation %d / %d", b, n_perm))
          run_one_perm(b)
        },
        mc.cores = n_cores
      )
      for (i in seq_along(chunk)) {
        b <- chunk[[i]]
        null_Rg[, b] <- perm_results[[i]][seq_len(n_genes)]
        null_global[b] <- perm_results[[i]][n_genes + 1L]
      }
      completed <- chunk_end
      save_null_perm_checkpoint(output_dir, null_Rg, null_global, completed, n_perm, seed)
    }
  } else {
    for (b in remaining) {
      if (b %% 100L == 0L || b == remaining[[1]]) {
        message(sprintf("  permutation %d / %d", b, n_perm))
      }
      result <- run_one_perm(b)
      null_Rg[, b] <- result[seq_len(n_genes)]
      null_global[b] <- result[n_genes + 1L]
      completed <- b
      if (completed %% checkpoint_every == 0L || completed == n_perm) {
        save_null_perm_checkpoint(output_dir, null_Rg, null_global, completed, n_perm, seed)
      }
    }
  }

  list(null_Rg = null_Rg, null_global = null_global, completed = completed)
}

run_stage_05_nulls <- function(atlas_env,
                               n_perm = NULL,
                               n_cores = NULL,
                               resume = TRUE,
                               checkpoint_every = NULL,
                               max_cores = 50L,
                               cfg = NULL) {
  output_dir <- atlas_env$output_dir
  gene_universe <- atlas_env$gene_universe
  obs <- atlas_env$repro
  cor_mats <- obs$cor_mats
  aggregate <- if (!is.null(obs$aggregate)) obs$aggregate else "median"

  if (is.null(cor_mats)) {
    stop("atlas_env$repro$cor_mats missing — rerun Stage 4 with the R_g implementation.")
  }
  if (is.null(obs$Rg)) {
    stop("atlas_env$repro$Rg missing — rerun Stage 4 with the R_g implementation.")
  }

  if (is.null(cfg) && !is.null(atlas_env$cfg)) cfg <- atlas_env$cfg

  if (is.null(n_perm)) {
    n_perm <- if (!is.null(cfg) && !is.null(cfg$n_perm)) cfg$n_perm else 1000L
  }
  n_perm <- as.integer(n_perm)

  if (is.null(checkpoint_every)) {
    checkpoint_every <- if (!is.null(cfg) && !is.null(cfg$null_checkpoint_every)) {
      cfg$null_checkpoint_every
    } else {
      50L
    }
  }

  if (is.null(max_cores) && !is.null(cfg) && !is.null(cfg$null_max_cores)) {
    max_cores <- cfg$null_max_cores
  }

  n_cores <- resolve_n_cores(requested = n_cores, max_cores = max_cores)
  seed <- if (!is.null(cfg) && !is.null(cfg$seed)) cfg$seed else NULL

  evt_n_exc <- if (!is.null(cfg) && !is.null(cfg$evt_n_exc)) as.integer(cfg$evt_n_exc) else 250L
  evt_gof_alpha <- if (!is.null(cfg) && !is.null(cfg$evt_gof_alpha)) cfg$evt_gof_alpha else 0.05
  evt_exceedance_min <- if (!is.null(cfg) && !is.null(cfg$evt_exceedance_min)) {
    as.integer(cfg$evt_exceedance_min)
  } else {
    10L
  }
  evt_gof_boot <- if (!is.null(cfg) && !is.null(cfg$evt_gof_boot)) as.integer(cfg$evt_gof_boot) else 199L

  message(sprintf("Stage 5: cross-gene null (n_perm = %d, n_cores = %d, aggregate = %s, resume = %s)",
                  n_perm, n_cores, aggregate, resume))

  perm_out <- run_cross_gene_null(
    cor_mats = cor_mats,
    gene_universe = gene_universe,
    n_perm = n_perm,
    n_cores = n_cores,
    aggregate = aggregate,
    resume = resume,
    checkpoint_every = checkpoint_every,
    output_dir = output_dir,
    seed = seed
  )

  null_Rg <- perm_out$null_Rg
  null_global <- perm_out$null_global
  rownames(null_Rg) <- gene_universe

  obs_Rg <- obs$Rg

  # Per-gene EVT/GPD calibration (empirical in the bulk, GPD in the tail).
  message("Stage 5: per-gene EVT/GPD calibration")
  evt_tbl <- calibrate_evt_genes(
    obs_vec = obs_Rg,
    null_mat = null_Rg,
    n_exc = evt_n_exc,
    gof_alpha = evt_gof_alpha,
    exceedance_min = evt_exceedance_min,
    gof_boot = evt_gof_boot,
    seed = seed,
    n_cores = n_cores
  )

  evt_FDR <- p.adjust(evt_tbl$evt_p, method = "BH")

  repro_df <- atlas_env$repro_df
  repro_df$empirical_p <- evt_tbl$empirical_p
  repro_df$evt_p <- evt_tbl$evt_p
  repro_df$evt_method <- evt_tbl$evt_method
  repro_df$evt_n_exc <- evt_tbl$evt_n_exc
  repro_df$evt_xi <- evt_tbl$evt_xi
  repro_df$evt_gof_p <- evt_tbl$evt_gof_p
  repro_df$empirical_FDR <- p.adjust(evt_tbl$empirical_p, method = "BH")
  repro_df$evt_FDR <- evt_FDR
  # Legacy aliases so Stage 6 / diagnostics / validation keep resolving.
  repro_df$shuffle_p <- evt_tbl$empirical_p
  repro_df$shuffle_FDR <- evt_FDR
  repro_df$shuffle_FDR_low_power <- repro_df$n_pairs_computable < length(cor_mats)

  # Global statistic: mean R_g across genes, calibrated the same way.
  obs_global <- mean(obs_Rg, na.rm = TRUE)
  global_cal <- calibrate_evt_gpd(
    obs = obs_global, null_vec = null_global,
    n_exc = evt_n_exc, gof_alpha = evt_gof_alpha,
    exceedance_min = evt_exceedance_min, gof_boot = evt_gof_boot,
    seed = seed
  )
  global_p_emp <- (sum(null_global >= obs_global, na.rm = TRUE) + 1) /
    (sum(is.finite(null_global)) + 1)

  results_dir <- file.path(output_dir, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(null_Rg, file.path(results_dir, "null_reproscore_matrix.rds"))

  provenance_note <- if (!is.null(cfg) && !is.null(cfg$atlas_provenance_note)) {
    cfg$atlas_provenance_note
  } else {
    ""
  }

  writeLines(
    c(
      "=== Global cross-gene null (mean R_g) ===",
      paste("timestamp:", Sys.time()),
      paste("n_perm:", n_perm),
      paste("n_cores:", n_cores),
      paste("aggregate:", aggregate),
      paste("obs_statistic_mean_Rg:", obs_global),
      paste("null_mean:", mean(null_global, na.rm = TRUE)),
      paste("null_sd:", stats::sd(null_global, na.rm = TRUE)),
      paste("empirical_p_one_sided:", global_p_emp),
      paste("evt_p:", global_cal$p),
      paste("evt_method:", global_cal$method),
      paste("evt_gof_p:", global_cal$gof_p),
      paste("atlas_provenance_note:", provenance_note)
    ),
    file.path(results_dir, "global_null.txt")
  )

  readr::write_tsv(repro_df, file.path(results_dir, "repro_scores_with_nulls.tsv"))

  delete_null_perm_checkpoint(output_dir)

  atlas_env$repro_df <- repro_df
  atlas_env$null_Rg <- null_Rg
  atlas_env$evt <- evt_tbl
  atlas_env$global_null <- list(
    obs_statistic = obs_global,
    empirical_p = global_p_emp,
    evt_p = global_cal$p,
    evt = global_cal,
    null_distribution = null_global
  )

  save.atlas.checkpoint(atlas_env, 5, results_dir = results_dir)

  message("Global empirical p: ", signif(global_p_emp, 4),
          " | global EVT p: ", signif(global_cal$p, 4),
          " (", global_cal$method, ")")

  invisible(list(
    atlas_env = atlas_env,
    global_p = global_p_emp,
    global_evt_p = global_cal$p,
    repro_df = repro_df
  ))
}
