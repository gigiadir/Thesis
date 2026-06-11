# Gene-Set Cohesion Benchmark

This benchmark checks whether genes expected to behave similarly (same complex/pathway/set) are closer to each other than unrelated genes under different method outputs (e.g., RankGenes / ExpressionTertiles / Residual).

## Inputs

- Reference sets: `.rds` list with shape `set_name -> character vector of genes`.
  - Example:
    - `list(MHC_I = c("HLA-A", "HLA-B", "HLA-C"), IFN = c("IFIT1", "IFIT3"))`
- One or more method inputs:
  - `feature_matrix` mode: gene-by-feature matrix (genes in rownames, numeric features in columns).
  - `distance_matrix` mode: square gene-by-gene matrix (distance or similarity).

## Main outputs

- `method_summary.csv` / `.rds`: per-method score table and ranking.
- `set_level_diagnostics.csv`: cohesion per biological set per method.
- `set_coverage.csv`: reference-set coverage after gene intersection.
- `excluded_sets.csv`: sets removed by minimum-size filter.
- `pair_distances_by_method.rds`: sampled pair table with labels and distances.
- Optional `sanity_check_summary.csv`: synthetic checks.

## Primary metric

- `separation_ratio = mean(within_set_distance) / mean(between_set_distance)`
  - Lower is better.
  - `< 1` means same-set genes are closer on average.

## Supporting metrics

- `mean_within`, `median_within`
- `mean_between`, `median_between`
- `cliffs_delta` (distributional effect size)
- `auc` (distance-based discrimination: same-set vs different-set pairs)
- Optional `permutation_p`
- Bootstrap confidence intervals for separation ratio and means

## Usage

Run from the `scripts` directory or pass absolute paths.

### Example 1: Compare methods from feature matrices

```bash
Rscript scripts/evaluate_gene_set_cohesion.R \
  --reference_rds data/reference_sets.rds \
  --method_name RankGenes \
  --method_input data/rankgenes_gene_features.rds \
  --method_mode feature_matrix \
  --method_name ExpressionTertiles \
  --method_input data/expression_tertiles_gene_features.rds \
  --method_mode feature_matrix \
  --method_name Residual \
  --method_input data/residual_gene_features.rds \
  --method_mode feature_matrix \
  --feature_distance correlation \
  --min_set_size 3 \
  --output_dir outputs/gene_set_cohesion
```

### Example 2: Compare precomputed gene-gene similarity matrices

```bash
Rscript scripts/evaluate_gene_set_cohesion.R \
  --reference_rds data/reference_sets.rds \
  --method_name RankGenesSplitSpace \
  --method_input data/rankgenes_jaccard_similarity.rds \
  --method_mode distance_matrix \
  --method_similarity TRUE \
  --method_name ResidualSplitSpace \
  --method_input data/residual_jaccard_similarity.rds \
  --method_mode distance_matrix \
  --method_similarity TRUE \
  --output_dir outputs/gene_set_cohesion_split_space
```

## Notes on fair comparison

- Keep the same `reference_rds` and sampling settings across methods.
- Use identical `seed`, `negative_ratio`, and pair caps for reproducible comparisons.
- Check `set_coverage.csv` to ensure one method is not advantaged by higher gene coverage.
- If runtime is high, reduce:
  - `--max_positive_pairs`
  - `--max_positive_pairs_per_set`
  - `--bootstrap_iters`
  - `--permutation_iters`
