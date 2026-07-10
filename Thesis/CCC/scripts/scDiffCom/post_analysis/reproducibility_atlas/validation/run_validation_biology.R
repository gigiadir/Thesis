#!/usr/bin/env Rscript
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
source(file.path(VALIDATION_DIR, "R", "v_biology.R"))

output_dir <- Sys.getenv("ATLAS_OUTPUT_DIR")
nulls_dir <- Sys.getenv("ATLAS_NULLS_DIR")
if (!nzchar(output_dir)) output_dir <- yaml::read_yaml(file.path(ATLAS_DIR, "config.yml"))$output_dir
if (!nzchar(nulls_dir)) nulls_dir <- output_dir

message("=== Validation biology only ===")
message("Started: ", Sys.time())
ctx <- load_validation_context(output_dir = output_dir, nulls_dir = nulls_dir, atlas_dir = ATLAS_DIR)
ctx <- run_v_biology(ctx, VALIDATION_DIR, controls = NULL)
report_gate(VALIDATION_DIR, "FINAL VERDICT")
verdict <- read_verdict(VALIDATION_DIR)
message("Summary:"); print(table(verdict$status))
message("Finished: ", Sys.time())
