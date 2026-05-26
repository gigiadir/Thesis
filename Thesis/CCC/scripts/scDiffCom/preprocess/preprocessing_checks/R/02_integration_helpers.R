# Requires SEURAT_DIR, N_FEATURES, N_PCS, CT_COL from config.

stamp_source <- function(fname, label) {
  message("Loading raw Seurat object...")
  path <- file.path(SEURAT_DIR, fname)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
  } else if (ext %in% c("rdata", "rda")) {
    miceadds::load.Rdata(path, "obj")
  } else {
    stop(paste("Unsupported file format:", ext, ". Supported formats: .rds, .RData, .rda"))
  }
  obj$source <- label
  message("Loaded ", fname, " — ", ncol(obj), " cells")
  obj
}

ensure_umap <- function(obj, label) {
  if ("umap" %in% names(obj@reductions)) {
    message(label, ": UMAP already present — skipping preprocessing.")
    return(obj)
  }
  message(label, ": no UMAP found — running NormalizeData → PCA → UMAP ...")
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = N_FEATURES, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = N_PCS, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:N_PCS, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = 1:N_PCS, verbose = FALSE)
  obj
}

ind_umap_plot <- function(obj, label) {
  col <- if (CT_COL %in% colnames(obj@meta.data)) CT_COL else "seurat_clusters"
  Seurat::DimPlot(obj, reduction = "umap", group.by = col, label = TRUE, repel = TRUE,
                  pt.size = 0.3) +
    ggplot2::ggtitle(label) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = grid::unit(0.35, "cm"))
}
