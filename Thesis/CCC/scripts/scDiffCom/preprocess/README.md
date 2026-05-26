# scDiffCom preprocess

Typical order:

1. `Rscript createPseudobulkMatrix.R --dataset_name <cohort>`
2. `Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name <cohort>`
3. Optional: `scDiffCom-Preprocess-Residual.R`, then knit `preprocessing_checks/index.Rmd` or `scDiffCom-Preprocessing-Checks.Rmd` (wrapper)

Batch jobs: `Main-scDiffComPreprocess.R` → `submit_jobs.sh` (written next to `scDiffComPreprocess.R`).

Output: `~/CCC-PreProcess/results-RankGenes/` (used by `../pipeline/`).
