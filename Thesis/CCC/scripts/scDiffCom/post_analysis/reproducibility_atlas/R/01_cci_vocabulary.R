# Stage 1 — fixed CCI vocabulary J (intersection or union + NA mask).

.cci_union_per_cohort <- function(malignant_by_cohort) {
  lapply(malignant_by_cohort, function(gene_list) {
    sort(unique(unlist(lapply(gene_list, function(df) df$CCI))))
  })
}

run_stage_01_vocab <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir
  cohorts <- atlas_env$cohorts
  malignant_by_cohort <- atlas_env$malignant_by_cohort

  message("Stage 1: build CCI vocabulary")
  cci_by_cohort <- .cci_union_per_cohort(malignant_by_cohort)

  for (ds in cohorts) {
    message(sprintf("  %s: %d unique DE CCIs (union across genes)", ds, length(cci_by_cohort[[ds]])))
  }

  J_intersection <- Reduce(intersect, cci_by_cohort)
  J_union        <- sort(unique(unlist(cci_by_cohort)))

  vocab_mode <- cfg$vocab_mode
  if (!vocab_mode %in% c("intersection", "union")) {
    stop("vocab_mode must be 'intersection' or 'union', got: ", vocab_mode)
  }

  J <- if (vocab_mode == "intersection") J_intersection else J_union
  if (length(J) == 0) {
    stop("CCI vocabulary J is empty under vocab_mode = ", vocab_mode)
  }

  na_mask <- lapply(cci_by_cohort, function(cci_set) {
    !J %in% cci_set
  })
  names(na_mask) <- cohorts

  support_loss <- 1 - length(J_intersection) / length(J_union)

  vocab_report <- data.frame(
    metric = c(
      "J_intersection", "J_union", "J_selected", "vocab_mode",
      "support_loss_fraction", "n_cohorts"
    ),
    value = c(
      length(J_intersection),
      length(J_union),
      length(J),
      vocab_mode,
      support_loss,
      length(cohorts)
    ),
    stringsAsFactors = FALSE
  )

  per_cohort <- purrr::imap_dfr(cci_by_cohort, function(cci_set, ds) {
    data.frame(
      cohort = ds,
      n_de_ccis = length(cci_set),
      n_in_J = sum(J %in% cci_set),
      n_missing_from_J = sum(!J %in% cci_set),
      stringsAsFactors = FALSE
    )
  })

  readr::write_tsv(vocab_report, file.path(output_dir, "results", "vocab_report.tsv"))
  readr::write_tsv(per_cohort, file.path(output_dir, "results", "vocab_per_cohort.tsv"))
  saveRDS(J, file.path(output_dir, "results", "J.rds"))
  saveRDS(na_mask, file.path(output_dir, "results", "vocab_na_mask.rds"))

  message(sprintf("  |J| = %d (mode: %s); support loss intersection/union = %.1f%%",
                  length(J), vocab_mode, 100 * support_loss))

  atlas_env$J        <- J
  atlas_env$na_mask  <- na_mask
  atlas_env$cci_by_cohort <- cci_by_cohort
  saveRDS(atlas_env, file.path(output_dir, "results", "stage01_atlas_env.rds"))

  invisible(atlas_env)
}
