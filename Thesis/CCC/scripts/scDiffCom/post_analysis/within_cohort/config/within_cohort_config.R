# Within-cohort phenotypes — edit TARGET_DATASET before knitting index.Rmd

TARGET_DATASET <- "Bill_HNSC"   # Puram_HNSC | Choi_HNSC | Kurten_HNSC

GENE_OF_INT       <- "AXL"
MALIGNANT_TYPE    <- "Tumor"
MIN_GENES_PER_LRI <- 2
K_CLUSTERS        <- 5
UMAP_SEED         <- 42

# AXL boxplot — label min/max LogFC only for top N ER pairs by LogFC spread
BOXPLOT_TOP_PAIRS_TO_LABEL <- 8
BOXPLOT_WIDTH  <- 18
BOXPLOT_HEIGHT <- 10

# Fisher enrichment (disabled in index.Rmd — slow; uncomment section 04 to run)
MIN_GO_SIZE   <- 5
MIN_DETECTED  <- 2
SAVE_EVERY    <- 10

# Volcano plots
SELECTED_GENES <- c("AXL", "ERBB2", "EGFR", "HLA-A", "MKI67")
P_ADJ_CUTOFF   <- 0.05
TOP_LABEL_N    <- 7
BOLD_GO_TERMS  <- c(
  "extracellular matrix organization",
  "angiogenesis",
  "extracellular structure organization"
)

WITHIN_COHORT_OUTPUT_DIR <- file.path(
  path.expand("~/Thesis/CCC/outputs/scDiffCom/plots/within_cohort"),
  TARGET_DATASET
)
