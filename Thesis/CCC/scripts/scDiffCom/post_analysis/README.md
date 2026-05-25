# scDiffCom Post-Analysis (H&N)

Cross-cohort comparison of scDiffCom results after pipeline jobs complete.

## Quick start

1. Ensure RDS files exist under `BASE_RESULTS_DIR` (see `R/00_setup.R`).
2. In RStudio, open `index.Rmd`, set working directory to this folder (`post_analysis/`), knit.
3. Or knit the legacy wrapper: `scripts/scDiffCom-PostAnalysis.Rmd` (from `scripts/`).

## Run order

| Section | File | Main outputs |
|---------|------|----------------|
| Setup | `index.Rmd` | Libraries, paths, helpers |
| 1 | `sections/01_load_and_filter.Rmd` | `*.malignant`, `Top.N.HN.Genes` |
| 2 | `sections/02_qc_celltypes.Rmd` | `plots/v2/qualityCheck/*.png` |
| 3 | `sections/03_build_cci_sets.Rmd` | `all_top_cci_sets`, `all_cci_sets_ranked` |
| 4 | `sections/04_goi_cross_cohort.Rmd` | `{GENE_OF_INT}_jaccard_umap.png`, `_ao_heatmap.png` |
| 5 | `sections/05_global_jaccard_umap.Rmd` | `global_*_umap_h&n.png` |
| 6 | `sections/06_gene_gene_jaccard.Rmd` | `gene_gene_*_heatmap_h&n.png` |
| 7 | `sections/07_lri_per_dataset.Rmd` | `{dataset}_lri_heatmap.png` |
| 8 | `sections/08_lri_cross_cohort_umap.Rmd` | `gene_dataset_lri_umap_h&n.png` |
| 9 | `sections/09_hvg_ccis.Rmd` | `cci_hvg_variance_ranking_h&n.csv` |
| 10 | `sections/10_limma_batch.Rmd` | `cci_post_analysis_matrices_h&n.rds`, CCI UMAPs |

Default plot directory: `~/Thesis/CCC/outputs/scDiffCom/plots/v5-split-by-rank` (`OUTPUT_DIR` in `R/00_setup.R`).

## Layout

```
post_analysis/
  index.Rmd              # driver — knit this
  config/hnsc_datasets.R # cohort names, colors, cancer filter
  R/                     # sourced helpers (same logic as monolith)
  sections/              # one analysis block per file
  archive/               # monolith backup + multi-cancer notes
```

## Adding a new H&N cohort

1. Run scDiffCom pipeline → RDS under `BASE_RESULTS_DIR/<DatasetName>/`.
2. Edit `config/hnsc_datasets.R` — add to `HNSC_DATASETS`, `ds_short_map`, `source_colors`.
3. Edit `sections/01_load_and_filter.Rmd` — load + filter blocks (copy an existing cohort).
4. Edit `sections/03_build_cci_sets.Rmd` — build `*_top_cci_sets` and list entries.
5. Re-knit `index.Rmd`.

## Multi-cancer (Breast / Lung)

Not maintained in the active path. See `archive/multi_cancer/README.md` and `archive/scDiffCom-PostAnalysis.Rmd.monolith`.

## Related docs

- [scDiffCom-PostAnalysis-Documentation.md](../../scDiffCom-PostAnalysis-Documentation.md) — section-level logic (updated for this layout)
- [scripts/README.md](../../README.md) — full scripts tree
