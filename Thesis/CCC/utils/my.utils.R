remotes::install_github("sqjin/CellChat")

library(Seurat)
library(dplyr)
library(CellChat)

split_by_gene_expression <- function(
    seurat_obj,
    gene,
    patient_col = "patient_id",   # e.g., "patient_id" or "orig.ident"
    assay       = DefaultAssay(seurat_obj),
    slot        = "data",          # "data" (normalized/log1p) or "counts" (raw)
    agg_fun     = c("mean", "median"),  # how to aggregate per patient
    cutoff_mode = c("median", "quantile", "custom"),
    q_low       = 0.5,             # used if cutoff_mode="quantile" (e.g., 0.25)
    custom_cut  = NULL,            # numeric cutoff if cutoff_mode="custom"
    add_meta    = TRUE             # store per-cell gene expression in meta.data
) {
  agg_fun     <- match.arg(agg_fun)
  cutoff_mode <- match.arg(cutoff_mode)
  
  # --- checks ---
  if (!patient_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data.", patient_col))
  }
  if (!(gene %in% rownames(seurat_obj[[assay]]))) {
    stop(sprintf("Gene '%s' not found in assay '%s'.", gene, assay))
  }
  if (!slot %in% c("data","counts","scale.data")) {
    stop("slot must be one of 'data', 'counts', or 'scale.data'.")
  }
  
  # --- pull per-cell expression for the gene ---
  gexpr <- GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, , drop = TRUE]
  
  if (add_meta) {
    seurat_obj[[paste0(gene, "_expr_", slot)]] <- gexpr
  }
  
  meta <- seurat_obj@meta.data %>%
    mutate(.cell = colnames(seurat_obj),
           .patient = .data[[patient_col]],
           .gexpr   = gexpr)
  
  # --- patient-level aggregation ---
  # mean or median expression across all cells of that patient
  patient_summary <- meta %>%
    group_by(.patient) %>%
    summarise(
      n_cells      = n(),
      expr_mean    = mean(.gexpr),
      expr_median  = median(.gexpr),
      pct_expr_pos = mean(.gexpr > 0) * 100,
      .groups = "drop"
    ) %>%
    mutate(expr_agg = if (agg_fun == "mean") expr_mean else expr_median)
  
  # --- choose cutoff ---
  cutoff <- switch(
    cutoff_mode,
    median   = stats::median(patient_summary$expr_agg, na.rm = TRUE),
    quantile = as.numeric(stats::quantile(patient_summary$expr_agg, probs = q_low, na.rm = TRUE)),
    custom   = {
      if (is.null(custom_cut)) stop("Provide custom_cut when cutoff_mode = 'custom'.")
      custom_cut
    }
  )
  
  # --- label patients ---
  # By default: Low = expr_agg <= cutoff ; High = expr_agg > cutoff
  patient_summary <- patient_summary %>%
    mutate(group = ifelse(expr_agg > cutoff, "High", "Low"))
  
  patients_high <- patient_summary %>% filter(group == "High") %>% pull(.patient)
  patients_low  <- patient_summary %>% filter(group == "Low")  %>% pull(.patient)
  
  # --- subset Seurat object into two groups by patients ---
  seu_high <- subset(seurat_obj, cells = rownames(meta %>% filter(.patient %in% patients_high)))
  seu_low  <- subset(seurat_obj, cells = rownames(meta %>% filter(.patient %in% patients_low)))
  
  # --- return ---
  list(
    high      = seu_high,
    low       = seu_low,
    cutoff    = cutoff,
    cutoff_mode = cutoff_mode,
    agg_fun   = agg_fun,
    assay     = assay,
    slot      = slot,
    groups_table = patient_summary %>%
      arrange(desc(expr_agg)) %>%
      rename(patient = .patient,
             expr_patient = expr_agg)
  )
}


create_cell_chat_object <- function(
    seurat_obj,
    assay = DefaultAssay(seurat_obj),
    slot = "data") {
  
  message("🔹 Extracting assay data (", assay, ":", slot, ") ...")
  data.input <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  
  message("🔹 Preparing metadata (group = cell_type) ...")
  if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
    stop("Meta.data must contain a 'cell_type' column.")
  }
  meta <- data.frame(group = seurat_obj@meta.data$cell_type)
  rownames(meta) <- colnames(seurat_obj)
  
  message("🔹 Creating CellChat object ...")
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  message("🔹 Loading CellChatDB.human ...")
  CellChatDB <- CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  message("🔹 Subsetting database to expressed genes ...")
  cellchat <- subsetData(cellchat)
  
  message("🔹 Identifying over-expressed genes ...")
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  message("🔹 Identifying over-expressed interactions ...")
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  message("🔹 Computing communication probabilities (type = triMean) ...")
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL, raw.use = TRUE)
  
  message("🔹 Filtering communications (min.cells = 10) ...")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  message("🔹 Computing pathway-level communication probabilities ...")
  cellchat <- computeCommunProbPathway(cellchat)
  
  message("🔹 Aggregating network and computing centrality ...")
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  message("✅ CellChat object successfully created with ", 
          length(unique(meta$group)), " cell groups.")
  
  # --- return ---
  cellchat
}