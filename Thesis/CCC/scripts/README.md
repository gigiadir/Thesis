# Thesis/CCC Scripts

| Path | Purpose |
|------|---------|
| [scDiffCom/post_analysis/](scDiffCom/post_analysis/) | **Post-analysis** — cross-cohort H&N comparison (knit `index.Rmd`) |
| [scDiffCom-PostAnalysis.Rmd](scDiffCom-PostAnalysis.Rmd) | Legacy wrapper → `post_analysis/index.Rmd` |
| `Main-scDiffComPipeline.R`, `scDiffComPipeline.R` | Cluster / local scDiffCom runs |
| `Main-scDiffComPreprocess.R`, `scDiffComPreprocess*.R` | Preprocessing before pipeline |
| [utils/](utils/) | Shared R utilities (`Utils.R`, Seurat helpers, etc.) |
| `*Tutorial.Rmd`, `*Playground.Rmd` | Exploratory notebooks |

Pipeline outputs are read from `~/CCC-scDiffCom/results/` (see post-analysis `R/00_setup.R`).
