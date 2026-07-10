# Stage 5 — gene-shuffle nulls and global permutation test.

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
                                      null_repro_scores,
                                      null_global,
                                      completed,
                                      n_perm,
                                      seed) {
  path <- .null_perm_checkpoint_path(output_dir)
  saveRDS(
    list(
      completed = as.integer(completed),
      null_repro_scores = null_repro_scores,
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

run_null_permutations <- function(Xtilde,
                                  cohort_idx_pairs,
                                  gene_universe,
                                  n_perm,
                                  n_cores,
                                  resume = TRUE,
                                  checkpoint_every = 50L,
                                  output_dir,
                                  seed = NULL) {
  n_perm <- as.integer(n_perm)
  checkpoint_every <- max(1L, as.integer(checkpoint_every))
  n_genes <- length(gene_universe)

  null_repro_scores <- matrix(NA_real_, n_genes, n_perm)
  null_global <- rep(NA_real_, n_perm)
  completed <- 0L

  ckpt <- load_null_perm_checkpoint(output_dir, n_perm, resume = resume)
  if (!is.null(ckpt)) {
    completed <- min(ckpt$completed, n_perm)
    null_repro_scores[, seq_len(completed)] <-
      ckpt$null_repro_scores[, seq_len(completed)]
    null_global[seq_len(completed)] <- ckpt$null_global[seq_len(completed)]
    message(sprintf(
      "Resuming from permutation %d / %d (checkpoint %s)",
      completed + 1L, n_perm, ckpt$timestamp
    ))
  }

  if (completed >= n_perm) {
    return(list(
      null_repro_scores = null_repro_scores,
      null_global = null_global,
      completed = completed
    ))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  run_one_perm <- function(b) {
    Xshuf <- shuffle_genes_within_cohort(Xtilde)
    null_repro <- compute_repro(Xshuf, cohort_idx_pairs)
    c(null_repro$ReproScore, mean(null_repro$ReproScore, na.rm = TRUE) - 0.5)
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
          if (b == chunk_start || (b - chunk_start) %% 100L == 0L) {
            message(sprintf("  permutation %d / %d", b, n_perm))
          }
          run_one_perm(b)
        },
        mc.cores = n_cores
      )
      for (i in seq_along(chunk)) {
        b <- chunk[[i]]
        null_repro_scores[, b] <- perm_results[[i]][seq_len(n_genes)]
        null_global[b] <- perm_results[[i]][n_genes + 1L]
      }
      completed <- chunk_end
      save_null_perm_checkpoint(
        output_dir, null_repro_scores, null_global, completed, n_perm, seed
      )
    }
  } else {
    for (b in remaining) {
      if (b %% 100L == 0L || b == remaining[[1]]) {
        message(sprintf("  permutation %d / %d", b, n_perm))
      }
      result <- run_one_perm(b)
      null_repro_scores[, b] <- result[seq_len(n_genes)]
      null_global[b] <- result[n_genes + 1L]
      completed <- b
      if (completed %% checkpoint_every == 0L || completed == n_perm) {
        save_null_perm_checkpoint(
          output_dir, null_repro_scores, null_global, completed, n_perm, seed
        )
      }
    }
  }

  list(
    null_repro_scores = null_repro_scores,
    null_global = null_global,
    completed = completed
  )
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
  Xtilde <- atlas_env$Xtilde
  cohort_idx_pairs <- atlas_env$cohort_pairs
  obs <- atlas_env$repro

  if (is.null(cfg) && !is.null(atlas_env$cfg)) {
    cfg <- atlas_env$cfg
  }

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

  message(sprintf("Stage 5: nulls (n_perm = %d, n_cores = %d, resume = %s)",
                  n_perm, n_cores, resume))

  perm_out <- run_null_permutations(
    Xtilde = Xtilde,
    cohort_idx_pairs = cohort_idx_pairs,
    gene_universe = gene_universe,
    n_perm = n_perm,
    n_cores = n_cores,
    resume = resume,
    checkpoint_every = checkpoint_every,
    output_dir = output_dir,
    seed = seed
  )

  null_repro_scores <- perm_out$null_repro_scores
  null_global <- perm_out$null_global
  rownames(null_repro_scores) <- gene_universe

  obs_global <- mean(obs$ReproScore, na.rm = TRUE) - 0.5
  p_emp <- rowMeans(null_repro_scores >= obs$ReproScore, na.rm = TRUE)

  repro_df <- atlas_env$repro_df
  repro_df$shuffle_p <- p_emp
  repro_df$shuffle_FDR <- p.adjust(p_emp, method = "BH")
  repro_df$shuffle_FDR_low_power <- TRUE

  global_p <- mean(null_global >= obs_global)

  results_dir <- file.path(output_dir, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(null_repro_scores, file.path(results_dir, "null_reproscore_matrix.rds"))

  provenance_note <- if (!is.null(cfg) && !is.null(cfg$atlas_provenance_note)) {
    cfg$atlas_provenance_note
  } else {
    ""
  }

  writeLines(
    c(
      "=== Global permutation null ===",
      paste("timestamp:", Sys.time()),
      paste("n_perm:", n_perm),
      paste("n_cores:", n_cores),
      paste("obs_statistic:", obs_global),
      paste("empirical_p_one_sided:", global_p),
      paste("atlas_provenance_note:", provenance_note)
    ),
    file.path(results_dir, "global_null.txt")
  )

  readr::write_tsv(repro_df, file.path(results_dir, "repro_scores_with_nulls.tsv"))

  delete_null_perm_checkpoint(output_dir)

  atlas_env$repro_df <- repro_df
  atlas_env$null_repro_scores <- null_repro_scores
  atlas_env$global_null <- list(
    obs_statistic = obs_global,
    empirical_p = global_p,
    null_distribution = null_global
  )

  save.atlas.checkpoint(atlas_env, 5, results_dir = results_dir)

  invisible(list(
    atlas_env = atlas_env,
    global_p = global_p,
    repro_df = repro_df
  ))
}
