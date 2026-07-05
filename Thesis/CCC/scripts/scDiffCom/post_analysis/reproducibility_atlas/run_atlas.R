#!/usr/bin/env Rscript
# Reproducibility atlas driver.
#
# Run all stages:
#   Rscript run_atlas.R
#
# Run one stage (loads prior checkpoint automatically):
#   Rscript run_atlas.R 0    # inspect + load
#   Rscript run_atlas.R 1    # CCI vocabulary
#   ...
#   Rscript run_atlas.R diag # diagnostics
#
# Interactive: open index.Rmd in RStudio.

ATLAS_DIR <- normalizePath(
  if (length(grep("^--file=", commandArgs(trailingOnly = FALSE))) > 0) {
    dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
  } else {
    "."
  },
  winslash = "/"
)

source(file.path(ATLAS_DIR, "R/00_atlas_setup.R"))
init.atlas.session(ATLAS_DIR)

stage_arg <- commandArgs(trailingOnly = TRUE)[1]
run_all <- is.na(stage_arg) || stage_arg == "" || stage_arg == "all"
run_diag <- identical(stage_arg, "diag")
stage_num <- if (!run_all && !run_diag) as.integer(stage_arg) else NA_integer_

message("=== Reproducibility Atlas ===")
message("Config: ", file.path(ATLAS_DIR, "config.yml"))
message("BASE_RESULTS_DIR: ", BASE_RESULTS_DIR)
message("Output: ", cfg$output_dir)

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
  atlas_env <- load.atlas.checkpoint(7)
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
      atlas_env <- load.atlas.checkpoint(s - 1)
    }
    atlas_env <- .run_stage(s, atlas_env)
  }
  message("Stage ", s, " complete.")
}

if (run_all) {
  atlas_env <- run_diagnostics(atlas_env)
  message("=== Atlas pipeline complete ===")
}
