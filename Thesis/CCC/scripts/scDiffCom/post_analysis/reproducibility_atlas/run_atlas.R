#!/usr/bin/env Rscript
# Reproducibility atlas driver — run one stage at a time.

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
source(file.path(ATLAS_DIR, "R/00_inspect_objects.R"))

cfg <- yaml::read_yaml(file.path(ATLAS_DIR, "config.yml"))
cfg$base_results_dir <- path.expand(cfg$base_results_dir)
cfg$output_dir       <- path.expand(cfg$output_dir)

BASE_RESULTS_DIR <- cfg$base_results_dir
MALIGNANT_CELLTYPE <- cfg$malignant_celltype
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(cfg$seed)

message("=== Reproducibility Atlas ===")
message("Config: ", file.path(ATLAS_DIR, "config.yml"))
message("BASE_RESULTS_DIR: ", BASE_RESULTS_DIR)

atlas_env <- run_stage_00_inspect(cfg)

file.copy(
  file.path(ATLAS_DIR, "config.yml"),
  file.path(cfg$output_dir, "results", "config_used.yml"),
  overwrite = TRUE
)

message("Stage 0 complete.")
