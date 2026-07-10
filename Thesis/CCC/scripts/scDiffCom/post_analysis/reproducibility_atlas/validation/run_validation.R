#!/usr/bin/env Rscript
# Reproducibility Atlas validation QC driver.
#
#   Rscript run_validation.R --through gate_a
#   Rscript run_validation.R --through gate_b
#   Rscript run_validation.R --through all
#   Rscript run_validation.R --priority star_only

parse_validation_args <- function(args) {
  opts <- list(
    output_dir = NULL,
    nulls_dir = NULL,
    through = "all",
    priority = "all",
    controls_file = NULL,
    fresh = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--output-dir", "-o")) {
      i <- i + 1L
      opts$output_dir <- args[[i]]
    } else if (arg %in% c("--nulls-dir", "-n")) {
      i <- i + 1L
      opts$nulls_dir <- args[[i]]
    } else if (arg == "--through") {
      i <- i + 1L
      opts$through <- args[[i]]
    } else if (arg == "--priority") {
      i <- i + 1L
      opts$priority <- args[[i]]
    } else if (arg == "--controls-file") {
      i <- i + 1L
      opts$controls_file <- args[[i]]
    } else if (arg == "--fresh") {
      opts$fresh <- TRUE
    } else if (arg %in% c("--help", "-h")) {
      cat(paste(
        "Usage: Rscript run_validation.R [options]",
        "",
        "Options:",
        "  --output-dir, -o PATH   Stages 0-4 results root (default: config.yml output_dir)",
        "  --nulls-dir, -n PATH    Stage 5 results root (default: config_batch.yml output_dir)",
        "  --through LEVEL         gate_a | gate_b | all (default: all)",
        "  --priority MODE         all | star_only",
        "  --controls-file PATH    YAML with positive/negative gene lists",
        "  --fresh                 Overwrite VERDICT.tsv",
        "  --help, -h              Show help",
        sep = "\n"
      ))
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
    i <- i + 1L
  }
  opts
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
VALIDATION_DIR <- normalizePath(dirname(script_path), winslash = "/")
ATLAS_DIR <- normalizePath(file.path(VALIDATION_DIR, ".."), winslash = "/")

# Prefer conda env libraries over site-wide R 4.6 packages (ABI mismatch).
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
source(file.path(VALIDATION_DIR, "R", "v1_vocab.R"))
source(file.path(VALIDATION_DIR, "R", "v2_tensor.R"))
source(file.path(VALIDATION_DIR, "R", "v3_center.R"))
source(file.path(VALIDATION_DIR, "R", "v4_reproscore.R"))
source(file.path(VALIDATION_DIR, "R", "v5_nulls.R"))
source(file.path(VALIDATION_DIR, "R", "v_loco.R"))
source(file.path(VALIDATION_DIR, "R", "v_biology.R"))

opts <- parse_validation_args(commandArgs(trailingOnly = TRUE))

cfg <- yaml::read_yaml(file.path(ATLAS_DIR, "config.yml"))
output_dir <- if (!is.null(opts$output_dir)) opts$output_dir else cfg$output_dir

if (opts$fresh) {
  verdict_path <- file.path(VALIDATION_DIR, "results", "VERDICT.tsv")
  if (file.exists(verdict_path)) file.remove(verdict_path)
}

init_verdict_file(VALIDATION_DIR)

message("=== Reproducibility Atlas Validation QC ===")
message("Started: ", Sys.time())
message("Output dir (0-4): ", path.expand(output_dir))
message("Through: ", opts$through)

ctx <- load_validation_context(
  output_dir = output_dir,
  nulls_dir = opts$nulls_dir,
  atlas_dir = ATLAS_DIR
)
message("Nulls dir (5): ", ctx$nulls_dir)
message("Stage 5 complete: ", ctx$stage5_complete)

should_run <- function(check_id, stage_num) {
  if (opts$priority == "star_only" && !(check_id %in% STAR_CHECKS)) {
    return(FALSE)
  }
  TRUE
}

run_stage <- function(name, fn, ...) {
  message("\n--- ", name, " ---")
  fn(ctx, VALIDATION_DIR, ...)
}

if (opts$through == "gate_a") {
  ctx <- run_v1_vocab(ctx, VALIDATION_DIR)
  ctx <- run_v2_tensor(ctx, VALIDATION_DIR)
  report_gate(VALIDATION_DIR, "GATE A (Stages 1-2)")
  if (gate_a_should_stop(VALIDATION_DIR)) {
    message("\nGATE A BLOCKED: eff_n = FAIL. Fix vocabulary before trusting Stage 4/5.")
    quit(save = "no", status = 1)
  }
}

if (opts$through %in% c("gate_b", "all")) {
  if (opts$through == "all") {
    ctx <- run_v1_vocab(ctx, VALIDATION_DIR)
    ctx <- run_v2_tensor(ctx, VALIDATION_DIR)
    report_gate(VALIDATION_DIR, "GATE A (Stages 1-2)")
    if (gate_a_should_stop(VALIDATION_DIR)) {
      message("\nGATE A BLOCKED: eff_n = FAIL. Fix vocabulary before trusting Stage 4/5.")
      quit(save = "no", status = 1)
    }
  }
  ctx <- run_v3_center(ctx, VALIDATION_DIR)
  ctx <- run_v4_reproscore(ctx, VALIDATION_DIR)
  ctx <- run_v5_nulls(ctx, VALIDATION_DIR)

  report_gate(VALIDATION_DIR, "GATE B (Stages 3-5)")
  decision <- gate_b_decision(VALIDATION_DIR)
  message("\nGATE B decision tree:")
  if (!is.na(decision$null_centered) && decision$null_centered == "FAIL") {
    message("  null_centered FAIL -> fix null (5a/5b) and rerun Stage 5")
  } else if (!is.na(decision$null_centered) && decision$null_centered == "PASS" &&
             !is.na(decision$raw_gap) && decision$raw_gap == "FAIL") {
    message("  null OK but no signal upstream -> revisit Stage 1 vocabulary / centering")
  } else if (!is.na(decision$reproscore_vs_n) && decision$reproscore_vs_n == "FAIL") {
    message("  reproscore_vs_n FAIL -> add min-n floor and rerun from Stage 4")
  } else if (!is.na(decision$global_test) && decision$global_test == "PASS") {
    message("  all key checks PASS -> proceed to LOCO + biology")
  } else if (!ctx$stage5_complete) {
    message("  Stage 5 incomplete -> rerun validation --through gate_b when nulls finish")
  }

  if (opts$through == "gate_b") {
    quit(save = "no", status = 0)
  }
}

if (opts$through == "all") {
  ctx <- run_v_loco(ctx, VALIDATION_DIR)

  controls <- NULL
  if (!is.null(opts$controls_file) && file.exists(opts$controls_file)) {
    controls <- yaml::read_yaml(opts$controls_file)
  }
  ctx <- run_v_biology(ctx, VALIDATION_DIR, controls = controls)

  report_gate(VALIDATION_DIR, "FINAL VERDICT")
  verdict <- read_verdict(VALIDATION_DIR)
  message("\nSummary:")
  print(table(verdict$status))

  key_artifacts <- c(
    "batch_before_after.png",
    "diagonal_heatmaps.png",
    "global_null_figure.png",
    "loco_rank_cor.tsv"
  )
  message("\nKey credibility artifacts:")
  for (a in key_artifacts) {
    path <- file.path(VALIDATION_DIR, "results", a)
    message(sprintf("  %s: %s", a, if (file.exists(path)) "OK" else "MISSING"))
  }
}

message("\nFinished: ", Sys.time())
message("Wrote: ", file.path(VALIDATION_DIR, "results", "VERDICT.tsv"))
