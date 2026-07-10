#!/usr/bin/env Rscript
# Resume validation from stage 5 tail through biology (skip stages 1-4).

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
VALIDATION_DIR <- normalizePath(dirname(script_path), winslash = "/")
ATLAS_DIR <- normalizePath(file.path(VALIDATION_DIR, ".."), winslash = "/")

conda_lib <- Sys.getenv("R_LIBS_SITE", unset = NA_character_)
if (!is.na(conda_lib) && nzchar(conda_lib) && dir.exists(conda_lib)) {
  .libPaths(conda_lib)
} else {
  env_prefix <- path.expand("~/.conda/envs/scDiffComPipeline_env/lib/R/library")
  if (dir.exists(env_prefix)) .libPaths(env_prefix)
}

source(file.path(ATLAS_DIR, "R", "atlas_helpers.R"))
source(file.path(VALIDATION_DIR, "R", "validation_helpers.R"))
source(file.path(VALIDATION_DIR, "R", "fast_reproscore_null.R"))
source(file.path(VALIDATION_DIR, "R", "v5_nulls.R"))
source(file.path(VALIDATION_DIR, "R", "v_loco.R"))
source(file.path(VALIDATION_DIR, "R", "v_biology.R"))

output_dir <- path.expand(
  "/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/scDiffCom/reproducibility_atlas/v1-v2"
)
nulls_dir <- path.expand(
  "/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/scDiffCom/reproducibility_atlas/v1-v2"
)

init_verdict_file(VALIDATION_DIR)
verdict <- read_verdict(VALIDATION_DIR)
if (nrow(verdict) > 0) {
  verdict <- verdict[!duplicated(verdict$check_id, fromLast = TRUE), , drop = FALSE]
  readr::write_tsv(verdict, file.path(VALIDATION_DIR, "results", "VERDICT.tsv"))
}

message("=== Validation tail (stage 5 remainder + LOCO + biology) ===")
message("Started: ", Sys.time())

ctx <- load_validation_context(
  output_dir = output_dir,
  nulls_dir = nulls_dir,
  atlas_dir = ATLAS_DIR
)
message("Stage 5 complete: ", ctx$stage5_complete)

existing <- read_verdict(VALIDATION_DIR)$check_id
need <- function(id) !(id %in% existing)

if (need("fastpath_equiv")) {
  message("--- fastpath_equiv ---")
  run_fastpath_equiv(ctx, VALIDATION_DIR)
}

if (need("null_centered")) {
  message("--- null_centered (full v5) ---")
  ctx <- run_v5_nulls(ctx, VALIDATION_DIR)
} else {
  message("--- stage 5 tail checks ---")
  null_repro <- ctx$null_repro
  if (is.null(null_repro) && !is.null(ctx$stage05_env)) {
    null_repro <- ctx$stage05_env$null_repro_scores
  }
  null_mean <- mean(null_repro, na.rm = TRUE)

  if (need("na_exchangeable")) {
    message("  na_exchangeable ...")
    n_scheme_b <- 50L
    set.seed(ctx$cfg$seed)
    scheme_b_scores <- replicate(n_scheme_b, {
      Xshuf <- shuffle_cci_within_gene(ctx$Xtilde)
      repro <- compute_repro(Xshuf, ctx$cohort_pairs_idx)
      mean(repro$ReproScore, na.rm = TRUE)
    })
    scheme_b_mean <- mean(scheme_b_scores)
    scheme_df <- data.frame(
      scheme = c("A_gene_shuffle", "B_cci_within_gene"),
      null_mean = c(null_mean, scheme_b_mean),
      n_perm = c(ncol(null_repro), n_scheme_b)
    )
    readr::write_tsv(scheme_df, file.path(VALIDATION_DIR, "results", "null_scheme_compare.tsv"))
    delta <- abs(null_mean - scheme_b_mean)
    if (delta > 0.02) {
      append_verdict(VALIDATION_DIR, "na_exchangeable", 5, "WARN", round(delta, 4),
                     "|delta| <= 0.02", "NA footprint may leak into null")
    } else {
      append_verdict(VALIDATION_DIR, "na_exchangeable", 5, "PASS", round(delta, 4),
                     "|delta| <= 0.02", "Null robust to NA structure")
    }
  }

  if (need("pval_dist")) {
    message("  pval_dist ...")
    repro_nulls <- ctx$repro_with_nulls
    if (is.null(repro_nulls) && !is.null(ctx$stage05_env)) {
      repro_nulls <- ctx$stage05_env$repro_df
    }
    if (!is.null(repro_nulls) && "shuffle_p" %in% names(repro_nulls)) {
      png(file.path(VALIDATION_DIR, "results", "pval_hist.png"), width = 600, height = 500)
      hist(repro_nulls$shuffle_p, breaks = 40, main = "Gene-shuffle p-value distribution",
           xlab = "shuffle_p", col = "steelblue", border = "white")
      dev.off()
      frac_low <- mean(repro_nulls$shuffle_p < 0.05, na.rm = TRUE)
      frac_high <- mean(repro_nulls$shuffle_p > 0.95, na.rm = TRUE)
      if (frac_high > 0.9) {
        append_verdict(VALIDATION_DIR, "pval_dist", 5, "INFO", round(frac_high, 3),
                       "uniform + spike at 0", "Consistent with no signal")
      } else {
        append_verdict(VALIDATION_DIR, "pval_dist", 5, "PASS", round(frac_low, 3),
                       "uniform + spike at 0", "Some real signal on null background")
      }
    }
  }

  if (need("global_test")) {
    message("  global_test ...")
    global_null <- NULL
    if (!is.null(ctx$stage05_env) && !is.null(ctx$stage05_env$global_null)) {
      global_null <- ctx$stage05_env$global_null
    } else if (file.exists(ctx$global_null_txt)) {
      lines <- readLines(ctx$global_null_txt)
      obs_line <- grep("obs_statistic", lines, value = TRUE)
      p_line <- grep("empirical_p", lines, value = TRUE)
      if (length(obs_line) > 0) {
        global_null <- list(
          obs_statistic = as.numeric(sub(".*: ", "", obs_line[1])),
          empirical_p = as.numeric(sub(".*: ", "", p_line[1]))
        )
      }
    }
    if (!is.null(global_null)) {
      obs <- global_null$obs_statistic
      null_dist <- global_null$null_distribution
      if (!is.null(null_dist)) {
        null_dist <- null_dist[is.finite(null_dist)]
      }
      emp_p <- global_null$empirical_p
      if (is.null(emp_p) && length(null_dist) > 0) {
        emp_p <- mean(null_dist >= obs)
      }
      if (!is.null(null_dist) && length(null_dist) > 0) {
        png(file.path(VALIDATION_DIR, "results", "global_null_figure.png"),
            width = 700, height = 500)
        hist(null_dist, breaks = 40, main = "Global permutation null",
             xlab = "Null statistic (mean ReproScore - 0.5)", col = "grey70", border = "white")
        abline(v = obs, col = "red", lwd = 2)
        legend("topright", legend = c(sprintf("obs=%.4f", obs), sprintf("p=%.4f", emp_p)),
               col = c("red", NA), lty = c(1, NA), bty = "n")
        dev.off()
      }
      if (is.finite(emp_p) && emp_p < 0.05) {
        append_verdict(VALIDATION_DIR, "global_test", 5, "PASS", round(emp_p, 4),
                       "p < 0.05", "Atlas-level reproducibility exists")
      } else {
        append_verdict(VALIDATION_DIR, "global_test", 5, "FAIL", round(emp_p, 4),
                       "p < 0.05", "No atlas-level signal; return to Stage 1/4")
      }
    } else {
      append_verdict(VALIDATION_DIR, "global_test", 5, "INFO", NA,
                     "p < 0.05", "global_null not available")
    }
  }
}

if (need("loco_stability")) {
  message("--- LOCO ---")
  ctx <- run_v_loco(ctx, VALIDATION_DIR)
}

if (need("controls") || need("enrichment")) {
  message("--- biology ---")
  ctx <- run_v_biology(ctx, VALIDATION_DIR, controls = NULL)
}

report_gate(VALIDATION_DIR, "FINAL VERDICT")
verdict <- read_verdict(VALIDATION_DIR)
message("\nSummary:")
print(table(verdict$status))

key_artifacts <- c(
  "batch_before_after.png", "diagonal_heatmaps.png",
  "global_null_figure.png", "loco_rank_cor.tsv"
)
message("\nKey credibility artifacts:")
for (a in key_artifacts) {
  path <- file.path(VALIDATION_DIR, "results", a)
  message(sprintf("  %s: %s", a, if (file.exists(path)) "OK" else "MISSING"))
}

message("\nFinished: ", Sys.time())
