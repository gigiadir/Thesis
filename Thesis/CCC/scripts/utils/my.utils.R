add_gene_expression_group <- function(
    seurat_obj,
    gene = "AXL",
    patient_col_options = c("orig.ident"),
    assay = DefaultAssay(seurat_obj),
    layer = "data",
    quantiles = c(1/3, 2/3)
) {
  # --- Input checks ---
  stopifnot(length(quantiles) == 2, quantiles[1] < quantiles[2])
  
  # --- Extract expression matrix ---
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (!gene %in% rownames(expr_matrix)) {
    stop(paste("Gene", gene, "not found in assay", assay))
  }
  
  # --- Extract per-cell expression ---
  gene_expr <- expr_matrix[gene, ]
  cell_metadata <- seurat_obj@meta.data
  
  available_cols <- colnames(cell_metadata)
  patient_col <- NULL
  for(col in patient_col_options) {
    if (col %in% available_cols) {
      patient_col <- col
      cat(sprintf("  Using patient column: %s\n", patient_col))
      break
    }
  }
  
  if (is.null(patient_col)) {
    stop(paste("Metadata must contain a patient identifier column. Available columns: ", available_cols))
  }
  
  # --- Annotate metadata with expression and patient ID ---
  cell_metadata$gene_expr_value <- gene_expr
  cell_metadata$patient_id_temp <- cell_metadata[[patient_col]]
  
  # --- Compute mean expression per patient ---
  patient_means <- aggregate(gene_expr_value ~ patient_id_temp, data = cell_metadata, FUN = mean)
  colnames(patient_means) <- c("patient_id", "mean_expr")
  
  # --- Determine quantile cutoffs ---
  qs <- quantile(patient_means$mean_expr, probs = quantiles, na.rm = TRUE)
  low_cutoff  <- qs[1]
  high_cutoff <- qs[2]
  
  # --- Assign group labels ---
  patient_means$GROUP <- with(patient_means, ifelse(
    mean_expr <= low_cutoff, "LOW",
    ifelse(mean_expr >= high_cutoff, "HIGH", NA)
  ))
  
  # --- Map patient → group ---
  patient_group_map <- setNames(patient_means$GROUP, patient_means$patient_id)
  
  # --- Assign each cell its group ---
  group_vector <- patient_group_map[cell_metadata[[patient_col]]]
  
  # --- Add to Seurat metadata (dynamic column name) ---
  group_col <- paste0(gene, "_EXP")
  seurat_obj@meta.data[[group_col]] <- group_vector  
  # --- Return updated Seurat and summary table ---
  list(
    seurat_obj = seurat_obj,
    patient_summary = setNames(patient_means, c("patient_id", paste0("mean_", gene), group_col))
  )
}

cosine_sim_na <- function(mat) {
  # mat: rows = CCIs, cols = genes
  # pairwise cosine between *columns*
  n <- ncol(mat)
  sim <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      xi <- mat[, i]; xj <- mat[, j]
      ok <- !is.na(xi) & !is.na(xj)
      if (sum(ok) < 2) next
      denom <- sqrt(sum(xi[ok]^2)) * sqrt(sum(xj[ok]^2))
      sim[i, j] <- if (denom == 0) NA else sum(xi[ok] * xj[ok]) / denom
    }
  }
  sim
}
