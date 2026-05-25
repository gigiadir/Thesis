# scDiffCom — scripts layout

End-to-end workflow for Head & Neck (and archived multi-cohort) scDiffCom analysis.

## Run order

1. **[preprocess/](preprocess/)** — pseudobulk matrices, patient splits, rank genes  
   Output: `~/CCC-PreProcess/results-RankGenes/`

2. **[pipeline/](pipeline/)** — per-gene scDiffCom jobs (cluster via `~/Scripts`)  
   Output: `~/CCC-scDiffCom/results/split-by-rank-genes/`

3. **[post_analysis/](post_analysis/)** — cross-cohort comparison (knit `index.Rmd`)  
   Output: `~/Thesis/CCC/outputs/scDiffCom/plots/`

## Entry points

| Task | File |
|------|------|
| Generate qsub script (Thesis) | `pipeline/Main-scDiffComPipeline.R` |
| Preprocess job generator | `preprocess/Main-scDiffComPreprocess.R` |
| Post-analysis report | `post_analysis/index.Rmd` |
| Legacy post-analysis wrapper | `../scDiffCom-PostAnalysis.Rmd` |

Cluster runs use copies under `~/Scripts/` — see [../SYNC_TO_SCRIPTS.md](../SYNC_TO_SCRIPTS.md).
