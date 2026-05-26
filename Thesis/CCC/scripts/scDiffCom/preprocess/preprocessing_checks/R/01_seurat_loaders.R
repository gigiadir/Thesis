Load.Seurat.With.Patients.Data <- function(ds_path, ps_path, gene_col, patient_col_in_seurat) {
  message("Loading raw Seurat object...")
  ext <- tolower(tools::file_ext(ds_path))
  if (ext == "rds") {
    obj <- readRDS(ds_path)
  } else if (ext %in% c("rdata", "rda")) {
    miceadds::load.Rdata(ds_path, "obj")
  } else {
    stop(paste("Unsupported file format:", ext, ". Supported formats: .rds, .RData, .rda"))
  }

  message("Loading patient summary (preprocess output)...")
  ps <- readRDS(ps_path)

  if (!patient_col_in_seurat %in% colnames(obj@meta.data)) {
    stop(sprintf("Error: Column '%s' not found in Seurat metadata.", patient_col_in_seurat))
  }

  message(sprintf("Merging %s classifications...", gene_col))

  mapping_vec <- setNames(ps[[gene_col]], ps$patient_id)
  obj@meta.data[[gene_col]] <- mapping_vec[as.character(obj@meta.data[[patient_col_in_seurat]])]

  list(
    seurat = obj,
    patient_summary = as.data.frame(ps)
  )
}

plot_patient_composition <- function(obj, patient_summary, gene_col_exp_name, ds_name,
                                     cell_type_column = "Consensus_Cell_Type",
                                     patient_col_in_seurat = "ident",
                                     use_proportions = FALSE) {
  gene_col <- grep("_EXP$", colnames(patient_summary), value = TRUE)

  meta <- obj@meta.data
  meta$patient <- as.character(meta[[patient_col_in_seurat]])
  meta <- meta[!is.na(meta$patient), ]

  all_patients <- unique(meta$patient)
  ordered_patients <- data.frame(patient = all_patients, stringsAsFactors = FALSE) %>%
    dplyr::left_join(patient_summary[, c("patient_id", "mean_expr")],
                     by = c("patient" = "patient_id")) %>%
    dplyr::arrange(dplyr::desc(mean_expr)) %>%
    dplyr::pull(patient)

  label_map <- setNames(patient_summary[[gene_col]], patient_summary$patient_id)
  x_labels <- setNames(
    ifelse(ordered_patients %in% names(label_map),
           paste0(ordered_patients, "\n", label_map[ordered_patients]),
           ordered_patients),
    ordered_patients
  )

  cell_counts <- meta %>%
    dplyr::filter(!is.na(.data[[cell_type_column]])) %>%
    dplyr::rename(cell_type = dplyr::all_of(cell_type_column)) %>%
    dplyr::count(patient, cell_type) %>%
    dplyr::group_by(patient) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cell_type = factor(cell_type,
                                    levels = c("Epithelial",
                                               sort(setdiff(unique(cell_type), "Epithelial")))))

  y_var <- if (use_proportions) "prop" else "n"
  y_label <- if (use_proportions) "Proportion" else "Cell Count"

  p <- ggplot(cell_counts, aes(x = patient, y = .data[[y_var]], fill = cell_type)) +
    geom_bar(stat = "identity")

  if (!use_proportions) {
    p <- p + geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.4) +
      scale_y_continuous(breaks = function(x) sort(unique(c(pretty(x), 50))))
  }

  p +
    scale_x_discrete(limits = ordered_patients, labels = x_labels) +
    labs(x = "Patient", y = y_label,
         title = paste(ds_name, "Cell Type Distribution:", gene_col_exp_name),
         fill = "Cell Type") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 10),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
}
