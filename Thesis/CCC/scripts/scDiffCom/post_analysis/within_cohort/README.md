# Within-Cohort Differential Communication Phenotypes (H&N)

Standalone module for **single-cohort** scDiffCom post-analysis — complementary to the cross-cohort driver in [`../index.Rmd`](../index.Rmd).

## Quick start

1. Ensure per-gene RDS files exist under `~/CCC-scDiffCom/results/split-by-rank-genes-v2/<DatasetName>/`.
2. Edit [`config/within_cohort_config.R`](config/within_cohort_config.R) — set `TARGET_DATASET` (`Kurten_HNSC`, `Puram_HNSC`, `Choi_HNSC`, or `Bill_HNSC`).
3. In RStudio, open `index.Rmd`, set working directory to this folder (`within_cohort/`), knit.

```bash
cd Thesis/CCC/scripts/scDiffCom/post_analysis/within_cohort
Rscript -e "rmarkdown::render('index.Rmd')"
```

## Outputs

Figures and tables are written to:

```
~/Thesis/CCC/outputs/scDiffCom/plots/within_cohort/{TARGET_DATASET}/
```

| File | Description |
|------|-------------|
| `AXL_{dataset}_CCI_boxplot.png` | GOI L-R LogFC across emitter→receiver pairs |
| `{dataset}_gene_gene_lri_heatmap.png` | Gene×gene cosine similarity (LRI LOGFC space) |
| `{dataset}_gene_umap_lri.png` | Gene UMAP with hierarchical clusters |
| `{dataset}_gene_umap_clusters.csv` | Cluster membership table |
| `{dataset}_fisher_results_final.rds` | Fisher GO-BP enrichment (all genes) |
| `{dataset}_fisher_checkpoint.rds` | Resume checkpoint (optional) |
| `{genes}_volcano.png` | Multi-gene Fisher volcano |
| `{gene}_volcano_no_legend.png` | Per-gene volcano |

## Sections

| # | File | Main outputs |
|---|------|--------------|
| Setup | `index.Rmd` | Config, shared helpers |
| 1 | `sections/01_load_dataset.Rmd` | Load one H&N cohort, filter malignant + unknown cell types |
| 2 | `sections/02_axl_boxplot.Rmd` | AXL (or `GENE_OF_INT`) boxplot |
| 3 | `sections/03_gene_gene_lri.Rmd` | Gene-gene heatmap + UMAP |
| 4 | `sections/04_gobp_fisher_volcano.Rmd` | Fisher + volcano (**disabled** — uncomment in `index.Rmd`) |

## Methods documentation

Thesis-ready description of representations, modules, and parameters: [METHODS.md](METHODS.md).

## Config reference

See [`config/within_cohort_config.R`](config/within_cohort_config.R) for `GENE_OF_INT`, Fisher thresholds, volcano gene list, and output paths.
