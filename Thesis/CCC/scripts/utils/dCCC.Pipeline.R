# dCCC.Pipeline.R
# Requires:
# - split_seurat_two_groups()   (from Seurat.Utils.R)
# - create_cell_chat_object()   (from CellChat.Utils.R)
# - plot_dccc_chunks()          (from your plotting util)
library(Seurat)
library(dplyr)
library(CellChat)

dCCC.Pipeline <- function(
    seurat_obj,
    # --- splitting args (passed to split_seurat_two_groups) ---
    split_args = list(
      mode = "gene_expression",
      gene = NULL,
      patient_col = "patient_id",
      assay = NULL,  # NULL => DefaultAssay(seurat_obj)
      slot  = "data",
      agg_fun = "mean",
      cutoff_mode = "quantile",
      q_low = 0.5,
      custom_cut = NULL,
      add_meta = TRUE,
      condition_col = NULL,
      group1_values = NULL,
      group2_values = NULL,
      case_insensitive = TRUE,
      per_cell_majority = TRUE,
      majority_threshold = 0.5,
      group1_label = "Group1",
      group2_label = "Group2",
      exclude_ambiguous = TRUE
    ),
    # --- CellChat creation args ---
    cellchat_args = list(
      assay = NULL,  # NULL => DefaultAssay(seurat_obj)
      slot  = "data"
    ),
    # --- labels for merged CellChat object ---
    add_names = c("Group1", "Group2"),
    # --- plotting via plot_dccc_chunks() ---
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
  # ---- sanity ----
  stopifnot(is.list(split_args), is.list(cellchat_args))
  stopifnot(is.list(pathway_chunks), length(pathway_chunks) > 0)
  if (length(add_names) != 2) stop("add_names must have length 2 (names for group1/group2).")
  if (!exists("plot_dccc_chunks")) stop("plot_dccc_chunks() must be defined and sourced.")
  
  # ---- resolve defaults tied to the provided object ----
  if (is.null(split_args$assay))    split_args$assay    <- DefaultAssay(seurat_obj)
  if (is.null(cellchat_args$assay)) cellchat_args$assay <- DefaultAssay(seurat_obj)
  
  # ---- 1) Split Seurat into two groups ----
  split_call <- c(list(seurat_obj = seurat_obj), split_args)
  split_res  <- do.call(split_seurat_two_groups, split_call)
  
  if (is.null(split_res$group1) || is.null(split_res$group2)) {
    stop("split_seurat_two_groups() did not return $group1/$group2.")
  }
  
  # ---- 2) Build CellChat per group ----
  cc1 <- do.call(create_cell_chat_object, c(list(seurat_obj = split_res$group1), cellchat_args))
  cc2 <- do.call(create_cell_chat_object, c(list(seurat_obj = split_res$group2), cellchat_args))
  
  # ---- 3) Merge ----
  merged_cc <- mergeCellChat(c(cc1, cc2), add.names = add_names)
  
  # ---- 4) Plot with plot_dccc_chunks() ----
  log_df <- plot_dccc_chunks(
    cellchat_merged   = merged_cc,
    pathway_chunks    = pathway_chunks,
    dataset_label     = dataset_label,
    comparison_label  = comparison_label,
    output_base_dir   = output_base_dir,
    measure           = measure,
    mode              = mode,
    file_format       = file_format,
    width_px          = width_px,
    height_px         = height_px,
    dpi               = dpi,
    verbose           = verbose
  )
  
  # ---- return useful artifacts ----
  list(
    split = split_res,
    cellchat_group1 = cc1,
    cellchat_group2 = cc2,
    merged_cellchat = merged_cc,
    output_log = log_df,
    output_dir = file.path(
      output_base_dir,
      gsub("[^A-Za-z0-9._-]+", "_", dataset_label),
      gsub("[^A-Za-z0-9._-]+", "_", comparison_label)
    )
  )
}