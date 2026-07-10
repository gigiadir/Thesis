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

is.unknown.celltype <- function(ct) {
  grepl("unknown|other|equivocal|multi", ct, ignore.case = TRUE)
}

filter.unknown.celltypes <- function(cci_df) {
  cci_df %>%
    filter(
      !is.unknown.celltype(EMITTER_CELLTYPE),
      !is.unknown.celltype(RECEIVER_CELLTYPE)
    )
}

report.filter.malignant <- function(orig, filt, label) {
  n_orig <- sum(vapply(orig, nrow, integer(1)))
  n_filt <- sum(vapply(filt, nrow, integer(1)))
  message(sprintf("  %-35s  %d → %d rows (removed %d)", label, n_orig, n_filt, n_orig - n_filt))
}

genes.from.scDiffCom.results <- function(cohort, results_dir = BASE_RESULTS_DIR) {
  rds_dir <- file.path(results_dir, cohort)
  files <- list.files(
    rds_dir,
    pattern = paste0("_", cohort, "_scDiffCom\\.rds$"),
    full.names = FALSE
  )
  sort(str_remove(files, paste0("_", cohort, "_scDiffCom\\.rds$")))
}

discover.gene.universe.from.cohorts <- function(cohorts,
                                                  results_dir = BASE_RESULTS_DIR,
                                                  verbose = TRUE) {
  results_dir <- path.expand(results_dir)
  genes_by_cohort <- stats::setNames(
    lapply(cohorts, genes.from.scDiffCom.results, results_dir = results_dir),
    cohorts
  )
  if (verbose) {
    for (ds in cohorts) {
      message(sprintf("  %s: %d genes on disk", ds, length(genes_by_cohort[[ds]])))
    }
  }
  gene_universe <- Reduce(intersect, genes_by_cohort)
  gene_union    <- sort(unique(unlist(genes_by_cohort)))
  list(
    genes_by_cohort = genes_by_cohort,
    gene_universe = gene_universe,
    gene_union = gene_union
  )
}

intersect.gene.names.across.cohorts <- function(malignant_by_cohort) {
  Reduce(intersect, lapply(malignant_by_cohort, names))
}

load.malignant.tables.for.cohort <- function(cohort,
                                             results_dir = BASE_RESULTS_DIR,
                                             genes = NULL,
                                             malignant_celltype = MALIGNANT_CELLTYPE,
                                             verbose = TRUE) {
  if (is.null(genes)) {
    scDiffComs <- load.dataset.scDiffComs(cohort, results_dir = results_dir)
    return(lapply(
      scDiffComs,
      filter.scDiffCom.cci_table_detected.for.malignant,
      malignant_celltype = malignant_celltype
    ))
  }

  out <- list()
  for (gene in genes) {
    rds_path <- file.path(
      results_dir, cohort,
      sprintf("%s_%s_scDiffCom.rds", gene, cohort)
    )
    if (!file.exists(rds_path)) {
      if (verbose) warning("Missing: ", rds_path)
      next
    }
    out[[gene]] <- filter.scDiffCom.cci_table_detected.for.malignant(
      readRDS(rds_path),
      malignant_celltype = malignant_celltype
    )
  }
  out
}

# Shared load + DE malignant filter (+ optional unknown-celltype filter) for one or more cohorts.
# Used by post_analysis/sections/01_load_and_filter.Rmd and reproducibility_atlas Stage 0.
load.and.filter.malignant.by.cohorts <- function(cohorts,
                                                   results_dir = BASE_RESULTS_DIR,
                                                   genes = NULL,
                                                   filter_unknown_celltypes = TRUE,
                                                   malignant_celltype = MALIGNANT_CELLTYPE,
                                                   checkpoint_dir = NULL,
                                                   verbose = TRUE) {
  results_dir <- path.expand(results_dir)
  discovered <- discover.gene.universe.from.cohorts(
    cohorts, results_dir = results_dir, verbose = verbose
  )
  genes_by_cohort <- discovered$genes_by_cohort

  if (is.null(genes)) {
    genes <- discovered$gene_universe
  }

  malignant_by_cohort <- stats::setNames(vector("list", length(cohorts)), cohorts)
  malignant_orig_by_cohort <- if (isTRUE(filter_unknown_celltypes)) {
    stats::setNames(vector("list", length(cohorts)), cohorts)
  } else {
    NULL
  }

  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  for (ds in cohorts) {
    ckpt_path <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, paste0(ds, "_malignant.rds"))
    } else {
      NULL
    }

    if (!is.null(ckpt_path) && file.exists(ckpt_path)) {
      if (verbose) message("  Loading checkpoint for ", ds)
      malignant_list <- readRDS(ckpt_path)
    } else {
      if (verbose) {
        message(
          "  Loading + filtering DE tables for ", ds,
          " (", length(genes), " genes)..."
        )
      }
      malignant_list <- load.malignant.tables.for.cohort(
        cohort = ds,
        results_dir = results_dir,
        genes = genes,
        malignant_celltype = malignant_celltype,
        verbose = verbose
      )
      if (!is.null(ckpt_path)) {
        saveRDS(malignant_list, ckpt_path)
        if (verbose) message("  Checkpoint saved: ", ckpt_path)
      }
    }

    if (isTRUE(filter_unknown_celltypes)) {
      malignant_orig_by_cohort[[ds]] <- malignant_list
      malignant_list <- lapply(malignant_list, filter.unknown.celltypes)
    }

    malignant_by_cohort[[ds]] <- malignant_list
  }

  list(
    malignant_by_cohort = malignant_by_cohort,
    malignant_orig_by_cohort = malignant_orig_by_cohort,
    genes_by_cohort = genes_by_cohort,
    gene_universe = genes
  )
}
