load.dataset.scDiffComs <- function(dataset_name,
                                    results_dir = BASE_RESULTS_DIR,
                                    do_switch   = FALSE) {
  message("Loading scDiffComs results for dataset: ", dataset_name)
  rds_dir   <- file.path(results_dir, dataset_name)
  rds_files <- list.files(
    rds_dir,
    pattern   = paste0("_", dataset_name, "_scDiffCom\\.rds$"),
    full.names = TRUE
  )

  scDiffCom_list <- list()

  for (rds_path in rds_files) {
    fname        <- basename(rds_path)
    gene         <- str_remove(fname, paste0("_", dataset_name, "_scDiffCom\\.rds$"))
    message("  Loading gene: ", gene)
    scDiffCom_obj <- readRDS(rds_path)
    scDiffCom_list[[gene]] <- scDiffCom_obj
  }

  return(scDiffCom_list)
}

filter.scDiffCom.cci_table_detected.for.malignant <- function(scDiffCom_obj,
                                                               malignant_celltype = MALIGNANT_CELLTYPE) {
  scDiffCom_obj@cci_table_detected %>%
    filter(
      (EMITTER_CELLTYPE %in% malignant_celltype | RECEIVER_CELLTYPE %in% malignant_celltype),
      IS_CCI_DE == TRUE,
      is.finite(LOGFC)
    )
}

filter.scDiffCom.cci_table_detected.for.celltypes <- function(scDiffCom_obj,
                                                               emitter_celltype,
                                                               receiver_celltype) {
  scDiffCom_obj@cci_table_detected %>%
    filter(
      EMITTER_CELLTYPE == emitter_celltype,
      RECEIVER_CELLTYPE == receiver_celltype,
      IS_CCI_DE == TRUE,
      is.finite(LOGFC)
    )
}
