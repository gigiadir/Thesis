#!/usr/bin/env Rscript
# Recover Stage 5 outputs from a completed null_perm_checkpoint (no re-permutation).
#
#   Rscript recover_stage05_outputs.R
#   Rscript recover_stage05_outputs.R --config config_batch.yml

parse_recover_args <- function(args) {
  opts <- list(config = NULL)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--config", "-c")) {
      i <- i + 1L
      opts$config <- args[[i]]
    } else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript recover_stage05_outputs.R [options]",
        "",
        "Options:",
        "  --config, -c PATH   Config YAML (default: config_batch.yml in atlas dir)",
        "  --help, -h          Show this help",
        sep = "\n"
      ))
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg, " (use --help)")
    }
    i <- i + 1L
  }
  opts
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
ATLAS_DIR <- normalizePath(dirname(script_path), winslash = "/")

# Prefer conda env libraries over site-wide R 4.6 packages (ABI mismatch).
conda_lib <- Sys.getenv("R_LIBS_SITE", unset = NA_character_)
if (!is.na(conda_lib) && nzchar(conda_lib) && dir.exists(conda_lib)) {
  .libPaths(conda_lib)
} else {
  env_prefix <- path.expand("~/.conda/envs/scDiffComPipeline_env/lib/R/library")
  if (dir.exists(env_prefix)) .libPaths(env_prefix)
}

source(file.path(ATLAS_DIR, "R/atlas_helpers.R"))
source(file.path(ATLAS_DIR, "R/stage_05_nulls.R"))

opts <- parse_recover_args(commandArgs(trailingOnly = TRUE))

config_path <- if (!is.null(opts$config)) {
  normalizePath(opts$config, winslash = "/", mustWork = TRUE)
} else {
  file.path(ATLAS_DIR, "config_batch.yml")
}

if (!file.exists(config_path)) {
  stop("Config not found: ", config_path)
}

cfg <- yaml::read_yaml(config_path)
cfg$base_results_dir <- path.expand(cfg$base_results_dir)
cfg$output_dir <- path.expand(cfg$output_dir)
if (nzchar(Sys.getenv("OUTPUT_DIR", unset = ""))) {
  cfg$output_dir <- path.expand(Sys.getenv("OUTPUT_DIR"))
}

results_dir <- file.path(cfg$output_dir, "results")
checkpoint_path <- file.path(results_dir, "null_perm_checkpoint.rds")
stage04_path <- file.path(results_dir, "stage04_atlas_env.rds")

if (!file.exists(stage04_path)) {
  stop("Missing stage04 checkpoint: ", stage04_path)
}
if (!file.exists(checkpoint_path)) {
  stop("Missing null perm checkpoint: ", checkpoint_path,
       " — cannot recover without completed permutations.")
}

ckpt <- readRDS(checkpoint_path)
if (!is.finite(ckpt$completed) || ckpt$completed < ckpt$n_perm) {
  stop("Null perm checkpoint incomplete: ", ckpt$completed, " / ", ckpt$n_perm)
}

message("=== Stage 5 output recovery (no re-permutation) ===")
message("Started: ", Sys.time())
message("Output dir: ", cfg$output_dir)
message("Null perms completed: ", ckpt$completed, " / ", ckpt$n_perm)
message(".libPaths(): ", paste(.libPaths(), collapse = ", "))

atlas_env <- readRDS(stage04_path)
atlas_env$output_dir <- cfg$output_dir
atlas_env$cfg <- cfg

out <- run_stage_05_nulls(
  atlas_env = atlas_env,
  n_perm = ckpt$n_perm,
  n_cores = 1L,
  resume = TRUE,
  checkpoint_every = cfg$null_checkpoint_every,
  max_cores = cfg$null_max_cores,
  cfg = cfg
)

message("Global empirical p: ", out$global_p)
message("Finished: ", Sys.time())
message("Wrote: ", file.path(results_dir, "repro_scores_with_nulls.tsv"))
message("Wrote: ", file.path(results_dir, "stage05_atlas_env.rds"))
