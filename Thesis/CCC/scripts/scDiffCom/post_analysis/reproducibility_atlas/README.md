# Cross-Cohort CCI Reproducibility Atlas

## Interactive (Rmd)

Open [`index.Rmd`](index.Rmd) in RStudio and run chunks section by section.

Each section can also be opened standalone under [`sections/`](sections/) — it loads the prior stage checkpoint automatically.

## CLI

```bash
cd reproducibility_atlas
Rscript run_atlas.R 0   # one stage
Rscript run_atlas.R     # all stages + diagnostics
```

## Shared load/filter logic

Loading and filtering reuse [`../R/01_load_filter.R`](../R/01_load_filter.R):

- `load.and.filter.malignant.by.cohorts()` — used by Stage 0 and post-analysis `01_load_and_filter.Rmd`
- `discover.gene.universe.from.cohorts()` — gene intersection from RDS filenames
- `filter.unknown.celltypes()` / `report.filter.malignant()`

## Config

Edit [`config.yml`](config.yml) before running. A copy is saved to `results/config_used.yml`.
