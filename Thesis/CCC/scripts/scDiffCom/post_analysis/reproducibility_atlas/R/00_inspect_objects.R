# Stage 0 — discover cohorts/genes, verify scDiffCom structure, build filtered lists.

is_unknown_celltype <- function(ct) {
  grepl("unknown|other|equivocal", ct, ignore.case = TRUE)
}

filter_unknown_celltypes <- function(cci_df) {
  cci_df %>%
    filter(
      !is_unknown_celltype(EMITTER_CELLTYPE),
      !is_unknown_celltype(RECEIVER_CELLTYPE)
    )
}

.genes_from_rds_dir <- function(results_dir, cohort) {
  rds_dir <- file.path(results_dir, cohort)
  files <- list.files(
    rds_dir,
    pattern = paste0("_", cohort, "_scDiffCom\\.rds$"),
    full.names = FALSE
  )
  sort(str_remove(files, paste0("_", cohort, "_scDiffCom\\.rds$")))
}

run_stage_00_inspect <- function(cfg) {
  results_dir <- path.expand(cfg$base_results_dir)
  output_dir  <- path.expand(cfg$output_dir)
  dir.create(file.path(output_dir, "results"), recursive = TRUE, showWarnings = FALSE)

  cohorts <- cfg$cohorts
  message("Stage 0: inspect objects")
  message("  results_dir: ", results_dir)
  message("  output_dir:  ", output_dir)

  genes_by_cohort <- stats::setNames(
    lapply(cohorts, .genes_from_rds_dir, results_dir = results_dir),
    cohorts
  )
  for (ds in cohorts) {
    message(sprintf("  %s: %d genes on disk", ds, length(genes_by_cohort[[ds]])))
  }

  gene_universe <- Reduce(intersect, genes_by_cohort)
  gene_union    <- sort(unique(unlist(genes_by_cohort)))
  message(sprintf("  gene intersection: %d | union: %d", length(gene_universe), length(gene_union)))

  if (length(gene_universe) == 0) {
    stop("Gene intersection is empty — check base_results_dir and cohort names.")
  }

  missing_tbl <- purrr::map_dfr(cohorts, function(ds) {
    missing <- setdiff(gene_universe, genes_by_cohort[[ds]])
    data.frame(
      cohort = ds,
      n_present = length(genes_by_cohort[[ds]]),
      n_in_universe = length(gene_universe),
      n_missing_from_universe = length(missing),
      missing_genes = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  readr::write_tsv(
    data.frame(gene = gene_universe, stringsAsFactors = FALSE),
    file.path(output_dir, "results", "gene_universe.tsv")
  )
  readr::write_tsv(missing_tbl, file.path(output_dir, "results", "genes_missing_per_cohort.tsv"))

  scDiffComs_by_cohort <- stats::setNames(
    lapply(cohorts, load.dataset.scDiffComs, results_dir = results_dir),
    cohorts
  )
  scDiffComs_by_cohort <- lapply(scDiffComs_by_cohort, function(lst) lst[gene_universe])

  malignant_by_cohort <- lapply(scDiffComs_by_cohort, function(lst) {
    lapply(lst, filter.scDiffCom.cci_table_detected.for.malignant)
  })

  if (isTRUE(cfg$filter_unknown_celltypes)) {
    malignant_by_cohort <- lapply(malignant_by_cohort, function(lst) {
      lapply(lst, filter_unknown_celltypes)
    })
  }

  inspect_lines <- c(
    "=== scDiffCom object inspection (Stage 0) ===",
    paste("timestamp:", Sys.time()),
    paste("results_dir:", results_dir),
    paste("cohorts:", paste(cohorts, collapse = ", ")),
    paste("gene_universe:", length(gene_universe), "genes"),
    paste("de_filter: IS_CCI_DE == TRUE (via filter.scDiffCom.cci_table_detected.for.malignant)"),
    paste("filter_unknown_celltypes:", cfg$filter_unknown_celltypes),
    ""
  )

  sample_cohort <- cohorts[[1]]
  sample_gene   <- gene_universe[[1]]
  sample_obj    <- scDiffComs_by_cohort[[sample_cohort]][[sample_gene]]
  sample_df     <- sample_obj@cci_table_detected

  inspect_lines <- c(
    inspect_lines,
    sprintf("Sample object: %s / %s", sample_cohort, sample_gene),
    paste("class:", paste(class(sample_obj), collapse = ", ")),
    paste("slots:", paste(slotNames(sample_obj), collapse = ", ")),
    paste("cci_table_detected columns:", paste(names(sample_df), collapse = ", ")),
    sprintf("cci_table_detected dim: %d rows x %d cols", nrow(sample_df), ncol(sample_df)),
    "",
    "IS_CCI_DE (all detected, tumor-involved sample rows):",
    capture.output(print(table(
      sample_df$IS_CCI_DE[
        sample_df$EMITTER_CELLTYPE %in% MALIGNANT_CELLTYPE |
          sample_df$RECEIVER_CELLTYPE %in% MALIGNANT_CELLTYPE
      ],
      useNA = "ifany"
    ))),
    "",
    "head(cci_table_detected, key columns):",
    capture.output(print(head(sample_df[, c(
      "CCI", "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "LRI",
      "LOGFC", "IS_CCI_DE", "BH_P_VALUE_DE", "REGULATION"
    )], 5))),
    "",
    "DE-filtered malignant rows per cohort (after unknown-type filter):"
  )

  de_counts <- purrr::imap(malignant_by_cohort, function(lst, ds) {
    n_rows <- sum(vapply(lst, nrow, integer(1)))
    sprintf("  %s: %d total DE rows across %d genes", ds, n_rows, length(lst))
  })
  inspect_lines <- c(inspect_lines, de_counts, "")

  writeLines(
    inspect_lines,
    file.path(output_dir, "results", "inspect_report.txt")
  )
  message("Wrote inspect_report.txt")

  atlas_env <- new.env(parent = globalenv())
  atlas_env$cfg                  <- cfg
  atlas_env$results_dir          <- results_dir
  atlas_env$output_dir           <- output_dir
  atlas_env$cohorts              <- cohorts
  atlas_env$gene_universe        <- gene_universe
  atlas_env$genes_by_cohort      <- genes_by_cohort
  atlas_env$scDiffComs_by_cohort <- scDiffComs_by_cohort
  atlas_env$malignant_by_cohort  <- malignant_by_cohort

  saveRDS(
    atlas_env,
    file.path(output_dir, "results", "stage00_atlas_env.rds")
  )
  message("Saved stage00_atlas_env.rds")

  invisible(atlas_env)
}
