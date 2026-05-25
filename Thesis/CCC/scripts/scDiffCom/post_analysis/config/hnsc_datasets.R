# H&N cohort configuration — single source for active post-analysis datasets.

HNSC_DATASETS <- c(
  "Kurten_HNSC",
  "Puram_HNSC",
  "Choi_HNSC",
  "Bill_HNSC"
)

CANCER_FILTER <- "H&N"

DATASET_LABELS <- HNSC_DATASETS

ds_short_map <- c(
  Kurten_HNSC = "Kurten",
  Puram_HNSC  = "Puram",
  Choi_HNSC   = "Choi",
  Bill_HNSC   = "Bill"
)

ds_cancer_map <- stats::setNames(rep(CANCER_FILTER, length(HNSC_DATASETS)), HNSC_DATASETS)

source_colors <- c(
  Kurten = "#F781BF",
  Puram  = "#999999",
  Choi   = "#66C2A5",
  Bill   = "#FC8D62"
)

cancer_shapes <- c(`H&N` = 15)

noise_ds_meta <- data.frame(
  ds_full     = HNSC_DATASETS,
  ds_short    = unname(ds_short_map[HNSC_DATASETS]),
  Cancer_Type = CANCER_FILTER,
  stringsAsFactors = FALSE
)

.cci_post_colors <- c(
  Kurten_HNSC = "#F781BF",
  Puram_HNSC  = "#999999",
  Choi_HNSC   = "#66C2A5",
  Bill_HNSC   = "#FC8D62"
)

.hnsc_malignant_orig_catalog <- function() {
  stats::setNames(
    lapply(HNSC_DATASETS, function(ds) get(paste0(ds, ".malignant.orig"), inherits = TRUE)),
    HNSC_DATASETS
  )
}
