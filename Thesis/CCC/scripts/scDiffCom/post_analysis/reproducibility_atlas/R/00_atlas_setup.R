# Session setup for reproducibility atlas (Rmd + run_atlas.R).

init.atlas.session <- function(atlas_dir = NULL) {
  if (is.null(atlas_dir)) {
    atlas_dir <- normalizePath(".", winslash = "/")
  } else {
    atlas_dir <- normalizePath(atlas_dir, winslash = "/")
  }

  post_analysis_dir <- normalizePath(file.path(atlas_dir, ".."), winslash = "/")
  project_scripts   <- normalizePath(file.path(post_analysis_dir, "..", ".."), winslash = "/")

  source(file.path(project_scripts, "utils/Utils.R"))
  source(file.path(post_analysis_dir, "R/00_setup.R"))
  source(file.path(post_analysis_dir, "R/01_load_filter.R"))
  source(file.path(post_analysis_dir, "config/hnsc_datasets.R"))
  source(file.path(atlas_dir, "R/atlas_helpers.R"))

  cfg <- yaml::read_yaml(file.path(atlas_dir, "config.yml"))
  cfg$base_results_dir <- path.expand(cfg$base_results_dir)
  cfg$output_dir       <- path.expand(cfg$output_dir)

  assign("ATLAS_DIR", atlas_dir, envir = .GlobalEnv)
  assign("POST_ANALYSIS_DIR", post_analysis_dir, envir = .GlobalEnv)
  assign("cfg", cfg, envir = .GlobalEnv)
  assign("BASE_RESULTS_DIR", cfg$base_results_dir, envir = .GlobalEnv)
  assign("MALIGNANT_CELLTYPE", cfg$malignant_celltype, envir = .GlobalEnv)
  assign("results_dir", file.path(cfg$output_dir, "results"), envir = .GlobalEnv)
  dir.create(get("results_dir", envir = .GlobalEnv), recursive = TRUE, showWarnings = FALSE)
  set.seed(cfg$seed)

  invisible(list(
    ATLAS_DIR = atlas_dir,
    cfg = cfg,
    results_dir = get("results_dir", envir = .GlobalEnv)
  ))
}

ensure.atlas.setup <- function(atlas_dir = NULL) {
  if (exists("cfg", envir = .GlobalEnv) && exists("results_dir", envir = .GlobalEnv)) {
    return(invisible(TRUE))
  }
  init.atlas.session(atlas_dir)
}

load.atlas.checkpoint <- function(stage_n, results_dir = get("results_dir", envir = .GlobalEnv)) {
  path <- file.path(results_dir, sprintf("stage%02d_atlas_env.rds", stage_n))
  if (!file.exists(path)) {
    stop("Missing checkpoint: ", path, " — run the previous section first.")
  }
  readRDS(path)
}

knit.atlas.section <- function(section_file, atlas_dir = get("ATLAS_DIR", envir = .GlobalEnv)) {
  path <- file.path(atlas_dir, "sections", section_file)
  if (!file.exists(path)) stop("Section not found: ", path)
  knitr::knit(path, envir = .GlobalEnv)
}
