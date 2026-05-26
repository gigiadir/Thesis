# scDiffCom preprocessing checks

Split from the monolithic `scDiffCom-Preprocessing-Checks.Rmd` (archived under `archive/`).

## Knit

- **Full report:** knit `index.Rmd` here, or the wrapper at `../scDiffCom-Preprocessing-Checks.Rmd` (set working directory to `preprocess/`).
- **Config:** edit `SPLIT_SOURCE` in [config/hnsc_preprocess_checks.R](config/hnsc_preprocess_checks.R) before knitting sections 04, 06, or 07.

| `SPLIT_SOURCE` | Methods included |
|----------------|------------------|
| `patient_zscore` | Patient-ZScore tertiles only |
| `rankgenes` | RankGenes tertiles only |
| `residual` | Residual tertiles only |
| `expression_quantile` | Expression tertiles only (`~/CCC-PreProcess/results-ExpressionQuantile`) |
| `all` | all four (default; skips missing dirs) |

## Sections

| File | Purpose |
|------|---------|
| `01_load_axl_cohorts.Rmd` | Load H&N Seurat + AXL patient summaries |
| `02_celltype_distribution.Rmd` | Cell-type composition per patient |
| `03_merge_combat_umap.Rmd` | ComBat integration UMAP (4 cohorts) |
| `04_split_category_distribution.Rmd` | LOW/MID/HIGH counts per patient **per split method** + facet comparison |
| `05_ctla4_atlas_sanity.Rmd` | Legacy CTLA4 table (`eval=FALSE`; needs external `helpers`) |
| `06_gene_dispersion.Rmd` | Pseudobulk vs split dispersion |
| `07_gene_profile_pca.Rmd` | Gene profile PCA (expression vs split space) |

## Partial runs

Run the `setup-preprocessing-checks` chunk in `index.Rmd`, then knit individual section files (or run chunks in order):

- Seurat QC: setup → 01 → 02
- Split-category comparison only: in config set `CATEGORY_QC_RUN_PER_SOURCE <- FALSE` and `CATEGORY_QC_RUN_COMPARISON <- TRUE`, then setup → run chunk `split-category-comparison-facets` in `sections/04_split_category_distribution.Rmd` (or the matching child chunk in `index.Rmd`)
- Split-category (all plots): setup → section 04 (requires grouped `.rds` under each split output dir)
- Pseudobulk/split QC: setup → 06 / 07

## Outputs

- Category plots: `~/Thesis/CCC/outputs/plots/QC/split_category_distribution/`
- Dispersion: `~/Thesis/CCC/outputs/plots/QC/gene_dispersion/`
- Profile PCA: `~/Thesis/CCC/outputs/plots/QC/gene_profile_pca/`
- Integration: `~/Thesis/CCC/outputs/scDiffCom/integration/`

## Pipeline order

See [../README.md](../README.md): pseudobulk → RankGenes (and optionally PatientZScore / Residual / ExpressionQuantile) → run this notebook.
