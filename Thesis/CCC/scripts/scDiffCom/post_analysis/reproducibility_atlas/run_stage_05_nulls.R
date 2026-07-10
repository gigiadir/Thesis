#!/usr/bin/env Rscript
# Stage 5 batch driver — gene-shuffle nulls (resumable, parallel).
#
#   Rscript run_stage_05_nulls.R
#   Rscript run_stage_05_nulls.R --n-perm 100 --ncores 4
#   Rscript run_stage_05_nulls.R --no-resume
#   Rscript run_stage_05_nulls.R --config config.yml

parse_cli_args <- function(args) {
  opts <- list(
    config = NULL,
    n_perm = NULL,
    n_cores = NULL,
    resume = TRUE,
    checkpoint_every = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--config", "-c")) {
      i <- i + 1L
      opts$config <- args[[i]]
    } else if (arg %in% c("--n-perm", "-n")) {
      i <- i + 1L
      opts$n_perm <- as.integer(args[[i]])
    } else if (arg %in% c("--ncores", "--cores")) {
      i <- i + 1L
      opts$n_cores <- as.integer(args[[i]])
    } else if (arg == "--checkpoint-every") {
      i <- i + 1L
      opts$checkpoint_every <- as.integer(args[[i]])
    } else if (arg == "--no-resume") {
      opts$resume <- FALSE
    } else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript run_stage_05_nulls.R [options]",
        "",
        "Options:",
        "  --config, -c PATH          Config YAML (default: config.yml in atlas dir)",
        "  --n-perm, -n N             Number of permutations (default: cfg$n_perm)",
        "  --ncores, --cores N        Parallel cores, clamped to cfg$null_max_cores",
        "  --checkpoint-every N       Save resume checkpoint every N perms (default: 50)",
        "  --no-resume                Ignore/delete resume and run from scratch",
        "  --help, -h                 Show this help",
        "",
        "Environment:",
        "  MC_CORES                     Default core count if --ncores omitted",
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

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))

config_path <- if (!is.null(opts$config)) {
  normalizePath(opts$config, winslash = "/", mustWork = TRUE)
} else {
  file.path(ATLAS_DIR, "config.yml")
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

checkpoint_path <- file.path(cfg$output_dir, "results", "stage04_atlas_env.rds")
if (!file.exists(checkpoint_path)) {
  stop("Missing checkpoint: ", checkpoint_path, " — run Stage 4 first.")
}

message("=== Reproducibility Atlas — Stage 5 (nulls) ===")
message("Started: ", Sys.time())
message("Atlas dir: ", ATLAS_DIR)
message("Config: ", config_path)
message("Output dir: ", cfg$output_dir)
message("Checkpoint: ", checkpoint_path)

atlas_env <- readRDS(checkpoint_path)
# Batch jobs must read/write on GPFS; stage04 RDS may carry a Dropbox symlink path.
atlas_env$output_dir <- cfg$output_dir
required_fields <- c("Xtilde", "repro", "cohort_pairs", "repro_df", "gene_universe", "output_dir")
missing_fields <- setdiff(required_fields, names(atlas_env))
if (length(missing_fields) > 0) {
  stop("stage04_atlas_env.rds missing fields: ", paste(missing_fields, collapse = ", "))
}

atlas_env$cfg <- cfg
n_perm <- if (!is.null(opts$n_perm)) opts$n_perm else cfg$n_perm
max_cores <- if (!is.null(cfg$null_max_cores)) cfg$null_max_cores else 50L
n_cores <- resolve_n_cores(requested = opts$n_cores, max_cores = max_cores)
checkpoint_every <- if (!is.null(opts$checkpoint_every)) {
  opts$checkpoint_every
} else if (!is.null(cfg$null_checkpoint_every)) {
  cfg$null_checkpoint_every
} else {
  50L
}

message(sprintf("n_perm = %d | n_cores = %d | resume = %s | checkpoint_every = %d",
                n_perm, n_cores, opts$resume, checkpoint_every))

if (!isTRUE(opts$resume)) {
  delete_null_perm_checkpoint(atlas_env$output_dir)
}

out <- run_stage_05_nulls(
  atlas_env = atlas_env,
  n_perm = n_perm,
  n_cores = n_cores,
  resume = opts$resume,
  checkpoint_every = checkpoint_every,
  max_cores = max_cores,
  cfg = cfg
)

message("Global empirical p: ", out$global_p)
message("Finished: ", Sys.time())
message("Wrote: ", file.path(cfg$output_dir, "results", "stage05_atlas_env.rds"))
