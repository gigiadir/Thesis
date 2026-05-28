# scDiffCom preprocess

Typical order:

1. `Rscript createPseudobulkMatrix.R --dataset_name <cohort>` (patients need `> 50` malignant cells; set one `malignant_cell_type` per dataset: `Epithelial`, `Malignant`, `Tumor`, or `Cancer`)
2. `Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name <cohort>` (production pipeline input)
3. Optional: `scDiffCom-Preprocess-ExpressionQuantile.R`, `scDiffCom-Preprocess-Residual.R`, `scDiffCom-Preprocess-PatientZScore.R`
4. Optional: `compareExpressionSplitEquivalence.R` (legacy vs pseudobulk expression quantile)
5. Optional: knit `preprocessing_checks/index.Rmd` or `scDiffCom-Preprocessing-Checks.Rmd` (wrapper)

Batch jobs: `Main-scDiffComPreprocess.R` → `submit_jobs.sh` (legacy per-gene Seurat quantile splits).

## Split methods

| Script | Output dir | Rule |
|--------|------------|------|
| `scDiffComPreprocess.R` | `~/CCC-PreProcess/results/` | Seurat malignant means → `quantile(1/3, 2/3)` |
| `scDiffCom-Preprocess-ExpressionQuantile.R` | `~/CCC-PreProcess/results-ExpressionQuantile/` | Pseudobulk → same quantile rule as legacy |
| `scDiffCom-Preprocess-RankGenes.R` | `~/CCC-PreProcess/results-RankGenes/` | Pseudobulk → column-wise rank tertiles |
| `scDiffCom-Preprocess-Residual.R` | `~/CCC-PreProcess/results-Residual/` | Pseudobulk → PC-regressed tertiles |
| `scDiffCom-Preprocess-PatientZScore.R` | `~/CCC-PreProcess/results-Patient-ZScore/` | Pseudobulk → row z-score tertiles |

## Expression quantile + legacy equivalence

```bash
Rscript createPseudobulkMatrix.R --dataset_name Kurten_HNSC
Rscript scDiffCom-Preprocess-ExpressionQuantile.R --dataset_name Kurten_HNSC
Rscript compareExpressionSplitEquivalence.R --dataset_name Kurten_HNSC
```

Summaries: `~/Thesis/CCC/outputs/QC/split_equivalence/expression_quantile_equivalence_*.csv`

**When splits should match exactly:** same patients, same `mean_expr`, same `--low_q` / `--high_q`.

**Known reasons for disagreement** (even with identical quantile code):

- Pseudobulk uses one `malignant_cell_type` per dataset; legacy uses all of `MALIGNANT_CELLS` — pick the label that matches your cohort (often `Tumor` for HNSC).
- Rebuild pseudobulk after changing filters (`createPseudobulkMatrix.R` uses `> 50` malignant cells by default, same as legacy).
- Different gene sets (legacy runs one gene per job; pseudobulk is full matrix).

If `max_abs_mean_expr_diff` is near zero but labels differ, check quantile tie-breaking; if diffs are large, filtering or cell-type definitions differ.

## Split-rule concordance (Figure 3.3)

Pairwise agreement of **HIGH / LOW** labels between partition rules, summarized as the mean fraction of concordant extreme-labeled patients across driver-panel genes (default: 123-gene scDiffCom panel).

**Metric (per gene, rule pair A vs B):**

- Comparable patients: both rules assign `LOW` or `HIGH` (exclude `MID` and missing).
- Per-gene score: fraction of comparable patients with the same extreme label.
- Matrix entry: mean of per-gene scores across aligned panel genes.

**Prerequisites:** pseudobulk plus all four preprocess outputs for the cohort:

```bash
Rscript createPseudobulkMatrix.R --dataset_name Kurten_HNSC
Rscript scDiffCom-Preprocess-ExpressionQuantile.R --dataset_name Kurten_HNSC
Rscript scDiffCom-Preprocess-RankGenes.R --dataset_name Kurten_HNSC
Rscript scDiffCom-Preprocess-Residual.R --dataset_name Kurten_HNSC
Rscript scDiffCom-Preprocess-PatientZScore.R --dataset_name Kurten_HNSC
Rscript plotSplitRuleConcordance.R --dataset_name Kurten_HNSC
```

**Outputs:** `~/Thesis/CCC/outputs/split_rule_concordance/{dataset}/`

- `split_rule_concordance_per_gene.csv` — long-form per-gene pairwise scores
- `split_rule_concordance_summary.csv` — mean concordance per rule pair
- `split_rule_concordance_matrix.csv` — symmetric 4×4 matrix
- `split_rule_concordance_heatmap.png` / `.pdf` — publication figure
- `split_rule_concordance_run_metadata.tsv` — run parameters
