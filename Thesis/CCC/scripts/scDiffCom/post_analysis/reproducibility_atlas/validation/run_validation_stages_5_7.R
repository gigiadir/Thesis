#!/usr/bin/env Rscript
# Complete validation stages 5–6 (nulls + LOCO/biology); append to existing VERDICT.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
VALIDATION_DIR <- normalizePath(dirname(script_path), winslash = "/")
ATLAS_DIR <- normalizePath(file.path(VALIDATION_DIR, ".."), winslash = "/")

.atlas_extra_libs <- c(
  path.expand("~/R/x86_64-redhat-linux-gnu-library/4.6"),
  "/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0"
)
.atlas_extra_libs <- .atlas_extra_libs[dir.exists(.atlas_extra_libs)]
if (length(.atlas_extra_libs)) .libPaths(c(.libPaths(), .atlas_extra_libs))

source(file.path(ATLAS_DIR, "R", "atlas_helpers.R"))
source(file.path(VALIDATION_DIR, "R", "validation_helpers.R"))
source(file.path(VALIDATION_DIR, "R", "fast_reproscore_null.R"))
source(file.path(VALIDATION_DIR, "R", "v5_nulls.R"))
source(file.path(VALIDATION_DIR, "R", "v_loco.R"))
source(file.path(VALIDATION_DIR, "R", "v_biology.R"))

output_dir <- Sys.getenv("ATLAS_OUTPUT_DIR")
nulls_dir <- Sys.getenv("ATLAS_NULLS_DIR")
if (!nzchar(output_dir)) output_dir <- yaml::read_yaml(file.path(ATLAS_DIR, "config.yml"))$output_dir
if (!nzchar(nulls_dir)) nulls_dir <- output_dir

message("=== Validation stages 5 + LOCO/biology ===")
message("Started: ", Sys.time())
ctx <- load_validation_context(output_dir = output_dir, nulls_dir = nulls_dir, atlas_dir = ATLAS_DIR)
ctx <- run_v5_nulls(ctx, VALIDATION_DIR)
report_gate(VALIDATION_DIR, "GATE B (Stages 3-5)")
decision <- gate_b_decision(VALIDATION_DIR)
message("GATE B decision: null=", decision$null_centered, " global=", decision$global_test)
ctx <- run_v_loco(ctx, VALIDATION_DIR)
ctx <- run_v_biology(ctx, VALIDATION_DIR, controls = NULL)
report_gate(VALIDATION_DIR, "FINAL VERDICT")
verdict <- read_verdict(VALIDATION_DIR)
message("Summary:"); print(table(verdict$status))
message("Finished: ", Sys.time())
