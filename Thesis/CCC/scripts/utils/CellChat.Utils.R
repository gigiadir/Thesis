remotes::install_github("sqjin/CellChat")

library(Seurat)
library(dplyr)
library(CellChat)

create_cell_chat_object <- function(
    seurat_obj,
    assay        = DefaultAssay(seurat_obj),
    slot         = "data",
    cell_type_col = "cell_type",
    min_cells    = 10
) {
  message("🔹 Extracting assay data (", assay, ":", slot, ") ...")
  data.input <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  message("🔹 Preparing metadata (group = ", cell_type_col, ") ...")
  meta_full <- seurat_obj@meta.data
  if (!cell_type_col %in% colnames(meta_full)) {
    stop("Meta.data must contain a '", cell_type_col, "' column.")
  }
  
  group_vec <- meta_full[[cell_type_col]]
  
  # 1) remove NA cell types
  keep_cells <- !is.na(group_vec)
  if (!any(keep_cells)) {
    stop("No cells with non-NA ", cell_type_col, " labels.")
  }
  
  group_vec <- group_vec[keep_cells]
  data.input <- data.input[, keep_cells, drop = FALSE]
  
  # 2) drop unused factor levels
  group_vec <- droplevels(factor(group_vec))
  
  meta <- data.frame(group = group_vec)
  rownames(meta) <- colnames(data.input)
  
  message("🔹 Creating CellChat object ...")
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  message("🔹 Loading CellChatDB.human ...")
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB
  
  message("🔹 Subsetting database to expressed genes ...")
  cellchat <- subsetData(cellchat)
  
  message("🔹 Identifying over-expressed genes ...")
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  message("🔹 Identifying over-expressed interactions ...")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  message("🔹 Computing communication probabilities ...")
  # safer to use defaults; your version with raw.use=TRUE is fine if your CellChat version supports it
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL)
  
  message("🔹 Filtering communications (min.cells = ", min_cells, ") ...")
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
  message("🔹 Computing pathway-level communication probabilities ...")
  cellchat <- computeCommunProbPathway(cellchat)
  
  message("🔹 Aggregating network and computing centrality ...")
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  message("✅ CellChat object successfully created with ", 
          length(unique(meta$group)), " cell groups (",
          paste(levels(meta$group), collapse = ", "), ").")
  
  cellchat
}


drop_cell_types_from_cell_chat_object <- function(cellchat, drop_types) {
  # Ensure input is valid
  if (!all(drop_types %in% unique(cellchat@idents))) {
    warning("Some specified cell types not found in CellChat object.")
  }
  
  # Keep only cells not in drop_types
  keep_cells <- names(cellchat@idents)[!(cellchat@idents %in% drop_types)]
  
  # Subset CellChat object
  cellchat@idents <- cellchat@idents[keep_cells]
  cellchat@meta <- cellchat@meta[keep_cells, , drop = FALSE]
  cellchat@data.signaling <- cellchat@data.signaling[, keep_cells, drop = FALSE]
  
  # Drop interactions involving the removed cell types
  cellchat@net$weight <- cellchat@net$weight[
    !(rownames(cellchat@net$weight) %in% drop_types),
    !(colnames(cellchat@net$weight) %in% drop_types),
    drop = FALSE
  ]
  
  # Optionally, re-identify group information
  cellchat@group <- factor(cellchat@idents)
  
  return(cellchat)
}


#' Plot CellChat comparison for pathway chunks and save to disk
#'
#' @param cellchat_merged A merged CellChat object with two groups (e.g., high vs low).
#' @param pathway_chunks  A list where each element is a character vector of pathway names.
#' @param dataset_label   A short label for the dataset/cohort (e.g., "HNSCC", "Breast").
#' @param comparison_label A human-readable comparison label for titles and filenames
#'                         (e.g., "AXL high vs low", "TNBC vs Luminal").
#' @param output_base_dir Directory where all outputs go (a subfolder per dataset is created).
#' @param measure         rankNet `measure` argument (default: "weight").
#' @param mode            rankNet `mode` argument (default: "comparison").
#' @param file_format     "png" or "pdf" (default: "png").
#' @param width_px        Image width in pixels (for PNG). Default: 2000.
#' @param height_px       Image height in pixels (for PNG). Default: 1500.
#' @param dpi             DPI for PNG. Default: 300.
#' @param verbose         If TRUE, prints progress.
#'
#' @return A data.frame log with chunk index, n_pathways, and output path.
plot_dccc_chunks <- function(
    cellchat_merged,
    pathway_chunks,
    dataset_label,
    comparison_label,
    output_base_dir = "../outputs/dCCC",
    measure = "weight",
    mode = "comparison",
    file_format = "png",
    width_px = 2000,
    height_px = 1500,
    dpi = 300,
    verbose = TRUE
) {
  stopifnot(is.list(pathway_chunks), length(pathway_chunks) > 0)
  if (!file_format %in% c("png","pdf")) {
    stop("file_format must be 'png' or 'pdf'.")
  }
  
  # Helper to make safe file/dir names
  sanitize <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    gsub("_+", "_", x)
  }
  
  # Folder structure: <output_base_dir>/<dataset_label>/<comparison_label>/
  dataset_dir    <- file.path(output_base_dir, sanitize(dataset_label))
  comparison_dir <- file.path(dataset_dir, sanitize(comparison_label))
  if (!dir.exists(comparison_dir)) dir.create(comparison_dir, recursive = TRUE)
  
  # Title template for plots
  title_text <- sprintf("%s — %s dCCC in signaling pathways", comparison_label, dataset_label)
  
  # Run over chunks
  out_log <- vector("list", length(pathway_chunks))
  for (i in seq_along(pathway_chunks)) {
    pathways_chunk <- pathway_chunks[[i]]
    if (length(pathways_chunk) == 0) {
      if (verbose) message(sprintf("[chunk %02d] empty pathway list — skipping", i))
      out_log[[i]] <- data.frame(
        chunk_index = i, n_pathways = 0, output = NA_character_, status = "skipped_empty",
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Filename
    base_name <- sprintf("Comparison_%s_%s_%02d", sanitize(comparison_label), sanitize(dataset_label), i)
    file_path <- file.path(comparison_dir, paste0(base_name, ".", file_format))
    
    if (verbose) {
      message(sprintf("[chunk %02d] %d pathways -> %s", i, length(pathways_chunk), basename(file_path)))
    }
    
    # Build plot
    p <- try(
      rankNet(
        cellchat_merged,
        mode = mode,
        measure = measure,
        signaling = pathways_chunk
      ) + ggplot2::labs(title = title_text),
      silent = TRUE
    )
    
    if (inherits(p, "try-error")) {
      warning(sprintf("[chunk %02d] rankNet failed: %s", i, as.character(p)))
      out_log[[i]] <- data.frame(
        chunk_index = i, n_pathways = length(pathways_chunk),
        output = NA_character_, status = "rankNet_error", stringsAsFactors = FALSE
      )
      next
    }
    
    # Save
    save_ok <- TRUE
    if (file_format == "png") {
      # Use ggsave for robustness
      tryCatch({
        ggplot2::ggsave(
          filename = file_path, plot = p,
          width = width_px, height = height_px,
          units = "px", dpi = dpi, limitsize = FALSE
        )
      }, error = function(e) {
        warning(sprintf("[chunk %02d] ggsave PNG failed: %s. Falling back to grDevices::png()", i, e$message))
        save_ok <<- FALSE
      })
      
      # Fallback device route if ggsave fails (rare)
      if (!save_ok) {
        grDevices::png(file_path, width = width_px, height = height_px, res = dpi)
        print(p)
        grDevices::dev.off()
      }
    } else if (file_format == "pdf") {
      tryCatch({
        ggplot2::ggsave(filename = file_path, plot = p, device = "pdf")
      }, error = function(e) {
        warning(sprintf("[chunk %02d] ggsave PDF failed: %s", i, e$message))
        file_path <<- NA_character_
      })
    }
    
    out_log[[i]] <- data.frame(
      chunk_index = i,
      n_pathways = length(pathways_chunk),
      output = file_path,
      status = ifelse(is.na(file_path), "save_error", "ok"),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out_log)
}



CellChat.Utils.Get.InfoFlow.Pathway <- function(cellchat, thresh = 0.05) {
  prob <- cellchat@netP$prob   # [source, target, pathway]
  pval <- cellchat@netP$pval
  
  # zero-out non-significant edges (same as CellChat)
  prob[pval > thresh] <- 0
  
  # information flow per pathway = sum over all source/target
  info <- apply(prob, 3, sum)
  return(info)  # named vector, names = pathway_name
}

CellChat.Utils.Get.InfoFlow.LR <- function(cellchat, thresh = 0.05) {
  prob <- cellchat@net$prob   # [source, target, LR_pair]
  pval <- cellchat@net$pval
  
  prob[pval > thresh] <- 0
  info <- apply(prob, 3, sum)
  return(info)  # named vector, names = LR pair IDs
}
