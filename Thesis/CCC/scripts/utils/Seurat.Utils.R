library(Seurat)
library(dplyr)



load.and.create.seurat <- function(mtx_path, genes_path, cells_path, metadata_path) {
  curr.data <-ReadMtx(
    mtx=mtx_path,
    features=genes_path,
    cells=cells_path,
    feature.column = 1,
    feature.sep = "\n",
    cell.sep = ",",
    skip.cell = 1
  )
  
  metadata <- read.csv(metadata_path)
  
  cells <- read.csv(cells_path)
  row.names(cells) <- cells$cell_name
  
  seurat <- CreateSeuratObject(
    counts = curr.data,
    meta.data = metadata
  )
  seurat <- AddMetaData(seurat, cells)
  Idents(seurat) <- seurat$cell_type
  
  return(seurat)
}

seurat.pipeline <- function(pbmc, is_tpm = F) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern="^MT-")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  all.genes = rownames(pbmc)
  
  if(is_tpm) {
    pbmc <- SetAssayData(pbmc, slot="data", new.data = log1p(pbmc@assays$RNA$counts))
  } else{
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
  pbmc <- FindNeighbors(pbmc)
  pbmc <- FindClusters(pbmc)
  pbmc <- RunUMAP(pbmc, dims=1:20)
  
  return (pbmc)
}

#################################################

split_by_gene_expression <- function(
    seurat_obj,
    gene,
    patient_col = "patient_id",   # e.g., "patient_id" or "orig.ident"
    assay       = DefaultAssay(seurat_obj),
    slot        = "data",          # "data" (normalized/log1p), "counts", or "scale.data"
    agg_fun     = c("mean", "median"),
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
  
  # --- per-cell expression for the gene ---
  gexpr <- GetAssayData(seurat_obj, assay = assay, slot = slot)[gene, , drop = TRUE]
  
  if (add_meta) {
    seurat_obj[[paste0(gene, "_expr_", slot)]] <- gexpr
  }
  
  meta <- seurat_obj@meta.data %>%
    mutate(.cell = colnames(seurat_obj),
           .patient = .data[[patient_col]],
           .gexpr   = gexpr)
  
  # --- patient-level aggregation ---
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
  patient_summary <- patient_summary %>%
    mutate(group = ifelse(expr_agg > cutoff, "High", "Low"))
  
  patients_high <- patient_summary %>% filter(group == "High") %>% pull(.patient)
  patients_low  <- patient_summary %>% filter(group == "Low")  %>% pull(.patient)
  
  # --- subset Seurat object by patients ---
  seu_high <- subset(seurat_obj, cells = rownames(meta %>% filter(.patient %in% patients_high)))
  seu_low  <- subset(seurat_obj, cells = rownames(meta %>% filter(.patient %in% patients_low)))
  
  # --- return ---
  list(
    group1_label = "High",
    group2_label = "Low",
    group1   = seu_high,
    group2   = seu_low,
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

# ========= 2) Metadata-based, patient-level split (e.g., TNBC vs ER+/Luminal) =========
# Works whether condition is already per-patient or per-cell.
split_by_patient_metadata <- function(
    seurat_obj,
    condition_col,                     # metadata column with subtype/status/etc.
    group1_values,                     # character vector for Group1 membership
    group2_values,                     # character vector for Group2 membership
    patient_col    = "patient_id",
    per_cell_majority = TRUE,          # if condition varies per cell, take majority vote within patient
    majority_threshold = 0.5,          # >50% default; tweak as needed
    case_insensitive = TRUE,
    group1_label = "Group1",
    group2_label = "Group2",
    exclude_ambiguous = TRUE           # drop patients who are ambiguous/unknown
) {
  if (!all(c(condition_col, patient_col) %in% colnames(seurat_obj@meta.data))) {
    stop("condition_col or patient_col not found in meta.data.")
  }
  
  meta <- seurat_obj@meta.data %>%
    mutate(.cell = colnames(seurat_obj),
           .patient = .data[[patient_col]],
           .cond_raw = .data[[condition_col]])
  
  # normalize case if needed
  if (case_insensitive) {
    meta <- meta %>% mutate(.cond = tolower(as.character(.cond_raw)))
    group1_values <- tolower(group1_values)
    group2_values <- tolower(group2_values)
  } else {
    meta <- meta %>% mutate(.cond = as.character(.cond_raw))
  }
  
  # derive patient-level label
  if (per_cell_majority) {
    patient_cond <- meta %>%
      group_by(.patient, .cond) %>%
      summarise(n = n(), .groups = "drop_last") %>%
      mutate(prop = n / sum(n)) %>%
      arrange(.patient, desc(prop), desc(n)) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(
        group = case_when(
          .cond %in% group1_values & prop > majority_threshold ~ group1_label,
          .cond %in% group2_values & prop > majority_threshold ~ group2_label,
          TRUE ~ "Ambiguous/Unknown"
        )
      ) %>%
      select(.patient, majority_value = .cond, majority_prop = prop, group)
  } else {
    # assume per-patient annotation: take first non-NA per patient
    patient_cond <- meta %>%
      group_by(.patient) %>%
      summarise(value = first(na.omit(.cond)), .groups = "drop") %>%
      mutate(
        group = case_when(
          value %in% group1_values ~ group1_label,
          value %in% group2_values ~ group2_label,
          TRUE ~ "Ambiguous/Unknown"
        )
      ) %>%
      rename(majority_value = value) %>%
      mutate(majority_prop = NA_real_)
  }
  
  if (exclude_ambiguous) {
    valid_patients <- patient_cond %>% filter(group %in% c(group1_label, group2_label)) %>% pull(.patient)
    patient_cond <- patient_cond %>% filter(.patient %in% valid_patients)
  }
  
  patients_g1 <- patient_cond %>% filter(group == group1_label) %>% pull(.patient)
  patients_g2 <- patient_cond %>% filter(group == group2_label) %>% pull(.patient)
  
  # subset
  cells <- tibble(.cell = colnames(seurat_obj)) %>%
    bind_cols(seurat_obj@meta.data %>% select(all_of(patient_col))) %>%
    rename(.patient = !!sym(patient_col))
  
  seu_g1 <- subset(seurat_obj, cells = cells %>% filter(.patient %in% patients_g1) %>% pull(.cell))
  seu_g2 <- subset(seurat_obj, cells = cells %>% filter(.patient %in% patients_g2) %>% pull(.cell))
  
  list(
    group1_label = group1_label,
    group2_label = group2_label,
    group1       = seu_g1,
    group2       = seu_g2,
    per_cell_majority = per_cell_majority,
    majority_threshold = majority_threshold,
    groups_table = patient_cond %>%
      arrange(desc(group), .patient)
  )
}

# ========= 3) One front-door that routes to gene/metadata mode =========
# Usage:
# - Gene mode:   split_seurat_two_groups(seu, mode="gene_expression", gene="AXL", ...)
# - Meta mode:   split_seurat_two_groups(seu, mode="metadata", condition_col="Subtype",
#                                        group1_values=c("TNBC"),
#                                        group2_values=c("Luminal A","Luminal B","ER+","HER2"))
split_seurat_two_groups <- function(
    seurat_obj,
    mode = c("gene_expression", "metadata"),
    # gene mode args
    gene = NULL,
    patient_col = "patient_id",
    assay = DefaultAssay(seurat_obj),
    slot = "data",
    agg_fun = c("mean", "median"),
    cutoff_mode = c("median", "quantile", "custom"),
    q_low = 0.5,
    custom_cut = NULL,
    add_meta = TRUE,
    # metadata mode args
    condition_col = NULL,
    group1_values = NULL,
    group2_values = NULL,
    case_insensitive = TRUE,
    per_cell_majority = TRUE,
    majority_threshold = 0.5,
    group1_label = "Group1",
    group2_label = "Group2",
    exclude_ambiguous = TRUE
) {
  mode <- match.arg(mode)
  
  if (mode == "gene_expression") {
    if (is.null(gene)) stop("Provide 'gene' for mode='gene_expression'.")
    return(
      split_by_gene_expression(
        seurat_obj = seurat_obj,
        gene = gene,
        patient_col = patient_col,
        assay = assay,
        slot = slot,
        agg_fun = agg_fun,
        cutoff_mode = cutoff_mode,
        q_low = q_low,
        custom_cut = custom_cut,
        add_meta = add_meta
      )
    )
  }
  
  if (mode == "metadata") {
    if (is.null(condition_col) || is.null(group1_values) || is.null(group2_values)) {
      stop("For mode='metadata', provide 'condition_col', 'group1_values', and 'group2_values'.")
    }
    return(
        split_by_patient_metadata(
        seurat_obj = seurat_obj,
        condition_col = condition_col,
        group1_values = group1_values,
        group2_values = group2_values,
        patient_col = patient_col,
        per_cell_majority = per_cell_majority,
        majority_threshold = majority_threshold,
        case_insensitive = case_insensitive,
        group1_label = group1_label,
        group2_label = group2_label,
        exclude_ambiguous = exclude_ambiguous
      )
    )
  }
  stop("Not supported mode.")
}

# ============================== EXAMPLES ==============================
# 1) AXL high vs low (per-patient mean; 50th percentile splitter)
# axl_split <- split_seurat_two_groups(
#   Kurten.2021.Seurat,
#   mode = "gene_expression",
#   gene = "AXL",
#   patient_col = "orig.ident",
#   assay = "alra",
#   slot = "data",
#   agg_fun = "mean",
#   cutoff_mode = "quantile",
#   q_low = 0.5
# )
# axl_split$group1_label  # "High"
# axl_split$group2_label  # "Low"
# axl_split$groups_table  # patient-level summary (mean/median, etc.)

# 2) Breast cancer: TNBC vs ER+/Luminal (metadata-based)
# tnbc_vs_erlum_split <- split_seurat_two_groups(
#   Breast.Seurat,
#   mode = "metadata",
#   condition_col = "Subtype",   # e.g., values like "TNBC","Luminal A","Luminal B","HER2","ER+"
#   group1_values = c("TNBC"),
#   group2_values = c("Luminal A","Luminal B","ER+","HER2"),
#   group1_label = "TNBC",
#   group2_label = "ERplus_Luminal",
#   patient_col = "patient_id",
#   per_cell_majority = TRUE,    # majority vote if Subtype varies across a patient's cells
#   majority_threshold = 0.5
# )
# tnbc_vs_erlum_split$group1_label
# tnbc_vs_erlum_split$group2_label
# tnbc_vs_erlum_split$groups_table  # shows each patient's assigned group and majority value/prop