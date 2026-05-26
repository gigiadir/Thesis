# Shared config for preprocessing_checks sections 04, 06, 07.

DATASETS <- c("Puram_HNSC", "Choi_HNSC", "Kurten_HNSC", "Bill_HNSC")

# "patient_zscore" | "rankgenes" | "residual" | "expression_quantile" | "all"
SPLIT_SOURCE <- "all"

split_sources <- switch(
  SPLIT_SOURCE,
  patient_zscore = "patient_zscore",
  rankgenes = "rankgenes",
  residual = "residual",
  expression_quantile = "expression_quantile",
  # facet_wrap fills row-wise (ncol=2): expression | z-score / rankgenes | residual
  all = c("expression_quantile", "patient_zscore", "rankgenes", "residual"),
  stop(
    "SPLIT_SOURCE must be 'patient_zscore', 'rankgenes', 'residual', ",
    "'expression_quantile', or 'all'. Got: ", SPLIT_SOURCE
  )
)

split_source_labels <- c(
  expression_quantile = "Expression tertiles (pseudobulk quantile)",
  patient_zscore = "Patient-ZScore tertiles",
  rankgenes = "RankGenes tertiles (within-patient rank)",
  residual = "Residual tertiles (PC-regressed pseudobulk)"
)

PSEUDOBULK_DIR <- path.expand("~/Thesis/CCC/outputs/RData_objects/pseudobulk_matrix")
ZSCORE_SPLITS_DIR <- path.expand("~/CCC-PreProcess/results-Patient-ZScore")
RANKGENES_DIR <- path.expand("~/CCC-PreProcess/results-RankGenes")
RESIDUAL_SPLITS_DIR <- path.expand("~/CCC-PreProcess/results-Residual")
EXPRESSION_QUANTILE_SPLITS_DIR <- path.expand("~/CCC-PreProcess/results-ExpressionQuantile")

SPLIT_BASE_DIRS <- c(
  patient_zscore = ZSCORE_SPLITS_DIR,
  rankgenes = RANKGENES_DIR,
  residual = RESIDUAL_SPLITS_DIR,
  expression_quantile = EXPRESSION_QUANTILE_SPLITS_DIR
)

CATEGORY_QC_DIR <- path.expand("~/Thesis/CCC/outputs/plots/QC/split_category_distribution")

# Section 04: set RUN_PER_SOURCE=FALSE and RUN_COMPARISON=TRUE to save only facet comparisons.
CATEGORY_QC_RUN_PER_SOURCE <- TRUE
CATEGORY_QC_RUN_COMPARISON <- TRUE
CATEGORY_QC_COMPARISON_FACET_NCOL <- 2L
# Gaps between facet panels in the 2x2 comparison plot (centimeters).
CATEGORY_QC_COMPARISON_PANEL_SPACING_X <- 2.5
CATEGORY_QC_COMPARISON_PANEL_SPACING_Y <- 1.5
DISPERSION_DIR <- path.expand("~/Thesis/CCC/outputs/plots/QC/gene_dispersion")
PROFILE_PCA_DIR <- path.expand("~/Thesis/CCC/outputs/plots/QC/gene_profile_pca")

EXPR_LEVELS <- c("LOW", "MID", "HIGH")
EXPR_PAL <- c(LOW = "#4E79A7", MID = "#F28E2B", HIGH = "#E15759")

LABEL_NUM <- c(LOW = 1, MID = 2, HIGH = 3)

N_PC <- 2L
MIN_PATIENTS <- 5L
MAX_PATIENT_NA_FRAC <- 0.25
LABEL_LOW_COR_QUANTILE <- 0.10
N_LABEL_GENES <- 12L

# Merge + ComBat integration (section 03)
SEURAT_DIR <- path.expand("~/scObjects")
INTEGRATION_OUTPUT_DIR <- path.expand("~/Thesis/CCC/outputs/scDiffCom/integration")
N_FEATURES <- 3000
N_PCS <- 30
CT_COL <- "Cell_Type"
SEURAT_V5 <- utils::packageVersion("Seurat") >= "5.0.0"

CELLTYPE_QC_DIR <- path.expand("~/Thesis/CCC/outputs/scDiffCom/QC")
