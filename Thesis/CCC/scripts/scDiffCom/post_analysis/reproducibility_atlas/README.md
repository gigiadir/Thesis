# Cross-Cohort CCI Reproducibility Atlas

## Layout

| Path | Purpose |
|------|---------|
| [`index.Rmd`](index.Rmd) | Interactive driver — run section chunks in order |
| [`sections/*.Rmd`](sections/) | **Stage code lives here** (inspect, vocab, tensor, …) |
| [`R/atlas_helpers.R`](R/atlas_helpers.R) | Shared functions (ReproScore, IDR, tensor helpers) |
| [`R/00_atlas_setup.R`](R/00_atlas_setup.R) | Session init, checkpoints, `knit.atlas.section()` |
| [`../R/01_load_filter.R`](../R/01_load_filter.R) | Load/filter scDiffCom (shared with post_analysis) |

## Interactive

Open [`index.Rmd`](index.Rmd) or any file under [`sections/`](sections/) in RStudio. Each section contains the full stage logic; standalone sections load the prior checkpoint if `atlas_env` is missing.

## CLI

```bash
Rscript run_atlas.R      # all sections
Rscript run_atlas.R 0    # one section (00_inspect.Rmd)
Rscript run_atlas.R diag
```

## Config

[`config.yml`](config.yml) — copied to `results/config_used.yml` on first run.
