#!/usr/bin/env Rscript
# Reproducibility atlas driver.
#
# Run all stages:
#   Rscript run_atlas.R
#
# Run one stage (loads prior checkpoint automatically):
#   Rscript run_atlas.R 0    # inspect + load
#   Rscript run_atlas.R 1    # CCI vocabulary
#   Rscript run_atlas.R 2    # tensor X
#   Rscript run_atlas.R 3    # centering Xtilde
#   Rscript run_atlas.R 4    # ReproScore
#   Rscript run_atlas.R 5    # nulls
#   Rscript run_atlas.R 6    # IDR
#   Rscript run_atlas.R 7    # atlas table
#   Rscript run_atlas.R diag # diagnostics

ATLAS_DIR <- normalizePath(
  if (length(grep("^--file=", commandArgs(trailingOnly = FALSE))) > 0) {
    dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
  } else {
    "."
  },
  winslash = "/"
)

POST_ANALYSIS_DIR <- normalizePath(file.path(ATLAS_DIR, ".."), winslash = "/")
PROJECT_SCRIPTS   <- normalizePath(file.path(POST_ANALYSIS_DIR, "..", ".."), winslash = "/")

source(file.path(PROJECT_SCRIPTS, "utils/Utils.R"))
source(file.path(POST_ANALYSIS_DIR, "R/00_setup.R"))
source(file.path(POST_ANALYSIS_DIR, "R/01_load_filter.R"))
source(file.path(POST_ANALYSIS_DIR, "config/hnsc_datasets.R"))

stage_files <- c(
  "00_inspect_objects.R",
  "01_cci_vocabulary.R",
  "02_build_tensor.R",
  "03_center.R",
  "04_reproscore.R",
  "05_nulls.R",
  "06_idr.R",
  "07_atlas.R",
  "diagnostics.R"
)
for (f in stage_files) {
  source(file.path(ATLAS_DIR, "R", f))
}

cfg <- yaml::read_yaml(file.path(ATLAS_DIR, "config.yml"))
cfg$base_results_dir <- path.expand(cfg$base_results_dir)
cfg$output_dir       <- path.expand(cfg$output_dir)

BASE_RESULTS_DIR <- cfg$base_results_dir
MALIGNANT_CELLTYPE <- cfg$malignant_celltype
results_dir <- file.path(cfg$output_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(cfg$seed)

stage_arg <- commandArgs(trailingOnly = TRUE)[1]
run_all <- is.na(stage_arg) || stage_arg == "" || stage_arg == "all"
run_diag <- identical(stage_arg, "diag")
stage_num <- if (!run_all && !run_diag) as.integer(stage_arg) else NA_integer_

message("=== Reproducibility Atlas ===")
message("Config: ", file.path(ATLAS_DIR, "config.yml"))
message("BASE_RESULTS_DIR: ", BASE_RESULTS_DIR)
message("Output: ", cfg$output_dir)

.load_stage_env <- function(n) {
  path <- file.path(results_dir, sprintf("stage%02d_atlas_env.rds", n))
  if (!file.exists(path)) {
    stop("Missing checkpoint: ", path, " — run stage ", n, " first.")
  }
  readRDS(path)
}

.run_stage <- function(n, atlas_env) {
  atlas_env <- switch(as.character(n),
    "0" = run_stage_00_inspect(cfg),
    "1" = run_stage_01_vocab(atlas_env),
    "2" = run_stage_02_tensor(atlas_env),
    "3" = run_stage_03_center(atlas_env),
    "4" = run_stage_04_reproscore(atlas_env),
    "5" = run_stage_05_nulls(atlas_env),
    "6" = run_stage_06_idr(atlas_env),
    "7" = run_stage_07_atlas(atlas_env),
    stop("Unknown stage: ", n)
  )
  invisible(atlas_env)
}

if (run_diag) {
  atlas_env <- .load_stage_env(7)
  atlas_env <- run_diagnostics(atlas_env)
  message("Diagnostics complete.")
  quit(save = "no", status = 0)
}

if (run_all) {
  stages <- 0:7
} else {
  if (is.na(stage_num) || stage_num < 0 || stage_num > 7) {
    stop("Stage must be 0–7, 'all', or 'diag'. Got: ", stage_arg)
  }
  stages <- stage_num
  message("Running stage(s): ", paste(stages, collapse = ", "))
}

atlas_env <- NULL
for (s in stages) {
  if (s == 0) {
    stage00_path <- file.path(results_dir, "stage00_atlas_env.rds")
    if (file.exists(stage00_path)) {
      message("Stage 0 checkpoint exists — delete stage00_atlas_env.rds to re-run from scratch.")
      atlas_env <- readRDS(stage00_path)
    } else {
      atlas_env <- run_stage_00_inspect(cfg)
      file.copy(
        file.path(ATLAS_DIR, "config.yml"),
        file.path(results_dir, "config_used.yml"),
        overwrite = TRUE
      )
    }
  } else {
    if (is.null(atlas_env)) {
      atlas_env <- .load_stage_env(s - 1)
    }
    atlas_env <- .run_stage(s, atlas_env)
  }
  message("Stage ", s, " complete.")
}

if (run_all) {
  atlas_env <- run_diagnostics(atlas_env)
  message("=== Atlas pipeline complete ===")
}
