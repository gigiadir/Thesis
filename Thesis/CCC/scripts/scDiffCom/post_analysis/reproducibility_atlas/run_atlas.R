#!/usr/bin/env Rscript
# CLI driver — executes section Rmds via knitr (same code as interactive index.Rmd).
#
#   Rscript run_atlas.R       # all sections
#   Rscript run_atlas.R 3     # section 03_center.Rmd only
#   Rscript run_atlas.R diag  # 07_diagnostics.Rmd

ATLAS_DIR <- normalizePath(
  if (any(grepl("^--file=", commandArgs(trailingOnly = FALSE)))) {
    dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
  } else ".",
  winslash = "/"
)

# Library resolution for stages 0-4 (run under system R 4.6):
#   * group R_4.6.0 (default site lib) provides the Bioconductor / DE stack
#     (edgeR, DESeq2, limma, NMF, data.table, tidyverse, ...).
#   * scDiffCom is only built for R 4.5, so append the group R_4.5.0 lib last —
#     it is pure-R and loads under 4.6 while its compiled deps resolve from 4.6.0.
#   * the user lib holds SuperRanker.
.atlas_extra_libs <- c(
  path.expand("~/R/x86_64-redhat-linux-gnu-library/4.6"),
  "/gpfs0/bgu-ofircohen/group/R_packages/R_4.5.0"
)
.atlas_extra_libs <- .atlas_extra_libs[dir.exists(.atlas_extra_libs)]
if (length(.atlas_extra_libs)) .libPaths(c(.libPaths(), .atlas_extra_libs))

source(file.path(ATLAS_DIR, "R/00_atlas_setup.R"))
init.atlas.session(ATLAS_DIR)

section_files <- c(
  "00_inspect.Rmd",
  "01_cci_vocabulary.Rmd",
  "02_build_tensor.Rmd",
  "03_center.Rmd",
  "04_reproscore.Rmd",
  "05_nulls.Rmd",
  "06_atlas.Rmd",
  "07_diagnostics.Rmd"
)

stage_arg <- commandArgs(trailingOnly = TRUE)[1]
run_all <- is.na(stage_arg) || stage_arg == "" || stage_arg == "all"
run_diag <- identical(stage_arg, "diag")

if (run_diag) {
  knit.atlas.section("07_diagnostics.Rmd")
  message("Diagnostics complete.")
  quit(save = "no")
}

if (run_all) {
  to_run <- section_files
} else {
  stage_num <- as.integer(stage_arg)
  if (is.na(stage_num) || stage_num < 0 || stage_num > 6) {
    stop("Stage must be 0–6, 'all', or 'diag'. Got: ", stage_arg)
  }
  to_run <- section_files[stage_num + 1]
  message("Running: ", to_run)
}

message("=== Reproducibility Atlas ===")
for (f in to_run) {
  message("\n--- ", f, " ---")
  knit.atlas.section(f)
}
if (run_all) message("\n=== Atlas pipeline complete ===")
