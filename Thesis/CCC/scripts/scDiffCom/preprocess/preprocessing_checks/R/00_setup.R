# Path resolution (set PREPROCESS_CHECKS_DIR in index.Rmd before sourcing).

if (!exists("PREPROCESS_CHECKS_DIR", inherits = FALSE)) {
  PREPROCESS_CHECKS_DIR <- normalizePath(".", winslash = "/")
}
if (!exists("PREPROCESS_DIR", inherits = FALSE)) {
  PREPROCESS_DIR <- normalizePath(file.path(PREPROCESS_CHECKS_DIR, ".."), winslash = "/")
}
if (!exists("PROJECT_SCRIPTS", inherits = FALSE)) {
  PROJECT_SCRIPTS <- normalizePath(
    file.path(PREPROCESS_CHECKS_DIR, "..", "..", ".."),
    winslash = "/"
  )
}

source(file.path(PREPROCESS_DIR, "scDiffComGenePanel.R"), local = FALSE)
source(file.path(PREPROCESS_DIR, "patientSplitUtils.R"), local = FALSE)
