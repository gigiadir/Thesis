# scDiffCom Post-Analysis — Documentation

> **Layout (2025):** Analysis code lives under [`scDiffCom/post_analysis/`](scDiffCom/post_analysis/). Knit [`index.Rmd`](scDiffCom/post_analysis/index.Rmd) or the wrapper [`scDiffCom-PostAnalysis.Rmd`](scDiffCom-PostAnalysis.Rmd). Full monolith backup: `post_analysis/archive/scDiffCom-PostAnalysis.Rmd.monolith`.

Active focus: **four H&N cohorts** (`Kurten_HNSC`, `Puram_HNSC`, `Choi_HNSC`, `Bill_HNSC`). See [`post_analysis/README.md`](scDiffCom/post_analysis/README.md) for run order and adding cohorts.

---

## Module map

| # | File | Former monolith chunk |
|---|------|---------------------|
| Setup | `index.Rmd` + `R/00_setup.R`, `R/01_load_filter.R`, `config/hnsc_datasets.R` | `setup`, `utils` |
| 1 | `sections/01_load_and_filter.Rmd` | `load-data`, `filter-unknown-celltypes` |
| 2 | `sections/02_qc_celltypes.Rmd` | `celltype-noise-diagnostics`, `qc-celltype-distributions` |
| 3 | `sections/03_build_cci_sets.Rmd` + `R/02_cci_helpers.R` | `cci-sets` |
| 4 | `sections/04_goi_cross_cohort.Rmd` | `gene-gene-jaccard` (GOI heatmaps) |
| 5 | `sections/05_global_jaccard_umap.Rmd` | `global-jaccard-umap` |
| 6 | `sections/06_gene_gene_jaccard.Rmd` | `gene-gene-clustering` |
| 7 | `sections/07_lri_per_dataset.Rmd` + `R/03_lri_helpers.R` | `genes-lri` |
| 8 | `sections/08_lri_cross_cohort_umap.Rmd` | `gene-dataset-lri-umap` |
| 9 | `sections/09_hvg_ccis.Rmd` | `post-analysis-hvg-ccis` |
| 10 | `sections/10_limma_batch.Rmd` | `post-analysis-limma-batch` |

### Within-cohort module (single dataset)

| # | File | Outputs |
|---|------|---------|
| Setup | `within_cohort/index.Rmd` + `config/within_cohort_config.R` | Configurable `TARGET_DATASET` |
| W1 | `within_cohort/sections/01_load_dataset.Rmd` | `TARGET.scDiffComs`, `TARGET.malignant` |
| W2 | `within_cohort/sections/02_axl_boxplot.Rmd` | `{GENE}_{dataset}_CCI_boxplot.png` |
| W3 | `within_cohort/sections/03_gene_gene_lri.Rmd` + `within_cohort/R/04_within_cohort_helpers.R` | `{dataset}_gene_gene_lri_heatmap.png`, `{dataset}_gene_umap_lri.png` |
| W4 | `within_cohort/sections/04_gobp_fisher_volcano.Rmd` | `{dataset}_fisher_results_final.rds`, volcano PNGs |

Knit [`within_cohort/index.Rmd`](scDiffCom/post_analysis/within_cohort/index.Rmd) — separate from the cross-cohort driver.

Thesis Methods prose: [post_analysis/METHODS.md](scDiffCom/post_analysis/METHODS.md) (cross-cohort), [within_cohort/METHODS.md](scDiffCom/post_analysis/within_cohort/METHODS.md) (within-cohort).

---

## 1. Analysis globals (`R/00_setup.R`)

Loads packages, sets `BASE_RESULTS_DIR`, `MALIGNANT_CELLTYPE` (`Tumor`), `OUTPUT_DIR`, builds `LRI_GO_BP`.

## 2. Load / filter helpers (`R/01_load_filter.R`)

- `load.dataset.scDiffComs()` — load per-gene RDS from pipeline output
- `filter.scDiffCom.cci_table_detected.for.malignant()` — DE CCIs involving tumor compartment

## 3. Load and filter (`sections/01_load_and_filter.Rmd`)

Loads four H&N `.malignant` lists; filters unknown cell types; intersects genes across cohorts (`Top.N.HN.Genes`).

## 4. CCI helpers (`R/02_cci_helpers.R`)

`build_gene_cci_sets`, Jaccard / Average Overlap distances, `build_gene_dataset_cci_feat_mat`, `select_hvg_ccis`, etc.

## 5. Build CCI sets (`sections/03_build_cci_sets.Rmd`)

`TOP_N_CCI` (500), `GENE_OF_INT` (AXL), `all_top_cci_sets`, `all_cci_sets_ranked`, `gene_cci_top500_summary`.

## 6. Gene-of-interest cross-cohort (`sections/04_goi_cross_cohort.Rmd`)

4×4 Jaccard and AO heatmaps for `GENE_OF_INT` across H&N datasets.

## 7–10. Downstream

- **§5** Global Gene×Dataset UMAP (Jaccard, AO, MNN) — `CANCER_FILTER = "H&N"`
- **§6** Gene×Gene consensus matrices (HVG-restricted CCIs)
- **§7–8** LRI-space heatmaps and cross-cohort UMAP
- **§9–10** HVG CCI selection, limma batch correction, validation UMAPs

---

Detailed prose for Breast/Lung multi-cohort design remains in the archived monolith; not maintained in the active H&N path.
