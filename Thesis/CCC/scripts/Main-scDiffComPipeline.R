# Backward-compatible entry point (Phase 2 layout).
# Run from Thesis/CCC/scripts/ (or: Rscript Main-scDiffComPipeline.R).
stub_dir <- {
  fa <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(fa)) {
    dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
  } else {
    normalizePath(".", winslash = "/")
  }
}
source(file.path(stub_dir, "scDiffCom/pipeline/Main-scDiffComPipeline.R"), local = FALSE)
