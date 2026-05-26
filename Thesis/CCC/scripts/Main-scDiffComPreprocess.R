# Backward-compatible entry point (Phase 2 layout).
# Run from Thesis/CCC/scripts/ (or: Rscript Main-scDiffComPreprocess.R).
stub_dir <- {
  fa <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(fa)) {
    dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
  } else {
    normalizePath(".", winslash = "/")
  }
}
source(file.path(stub_dir, "scDiffCom/preprocess/Main-scDiffComPreprocess.R"), local = FALSE)
