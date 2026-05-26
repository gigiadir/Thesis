# Thesis/CCC Scripts

| Path | Purpose |
|------|---------|
| [scDiffCom/](scDiffCom/) | **scDiffCom workflow** — preprocess → pipeline → post_analysis |
| [scDiffCom/post_analysis/](scDiffCom/post_analysis/) | Cross-cohort H&N comparison (knit `index.Rmd`) |
| [scDiffCom-PostAnalysis.Rmd](scDiffCom-PostAnalysis.Rmd) | Legacy wrapper → `post_analysis/index.Rmd` |
| [tutorials/](tutorials/) | CellChat, NicheNet, exploratory notebooks |
| [utils/](utils/) | Shared R utilities (`Utils.R`, Seurat helpers, etc.) |
| [SYNC_TO_SCRIPTS.md](SYNC_TO_SCRIPTS.md) | Copy pipeline edits to `~/Scripts` for cluster |

## scDiffCom layout

```
scDiffCom/
  preprocess/   # pseudobulk, rank genes → CCC-PreProcess
  pipeline/     # scDiffCom jobs → CCC-scDiffCom (sync to ~/Scripts)
  post_analysis/  # knit index.Rmd
```

## Root stubs

| Stub | Redirects to |
|------|----------------|
| `Main-scDiffComPipeline.R` | `scDiffCom/pipeline/Main-scDiffComPipeline.R` |
| `Main-scDiffComPreprocess.R` | `scDiffCom/preprocess/Main-scDiffComPreprocess.R` |

Pipeline outputs: `~/CCC-scDiffCom/results/`. Post-analysis plot paths: `scDiffCom/post_analysis/R/00_setup.R`.
