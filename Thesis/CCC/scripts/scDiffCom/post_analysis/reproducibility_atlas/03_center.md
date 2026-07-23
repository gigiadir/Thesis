# Stage 3 — Cohort centering




``` r
output_dir <- atlas_env$output_dir
cohorts <- atlas_env$cohorts
gene_universe <- atlas_env$gene_universe
X_full <- atlas_env$X_full
J <- atlas_env$J

message("Stage 3: cohort mean-gene subtraction (full CCIs), then restrict to J")
```

```
## Stage 3: cohort mean-gene subtraction (full CCIs), then restrict to J
```

``` r
# Step 2 (advisor): subtract the cohort "mean gene" profile mu_c[j] from every
# gene's vector, computed over ALL cohort-tested CCIs (unbiased centering). Mean
# subtraction only — no SD scaling (Spearman is scale-invariant per vector, and
# the advisor removes only the cohort-common additive signal).
Xtilde_full <- lapply(cohorts, function(ds) subtract_cohort_mean_profile(X_full[[ds]]))
names(Xtilde_full) <- cohorts

# Restrict the centered tensor to the intersection J for downstream scoring, so
# all cohorts share identical, aligned CCI columns.
Xtilde <- lapply(cohorts, function(ds) Xtilde_full[[ds]][, J, drop = FALSE])
names(Xtilde) <- cohorts
assert_tensor_alignment(Xtilde, gene_universe, J, cohorts)

set.seed(cfg$seed)
hk_genes <- sample(gene_universe, min(75, length(gene_universe)))

centering_baseline <- purrr::imap_dfr(cohorts, function(ds, idx) {
  mat <- X_full[[ds]]
  m_hk <- colMeans(mat[hk_genes, , drop = FALSE], na.rm = TRUE)
  m_all <- colMeans(mat, na.rm = TRUE)
  finite <- is.finite(m_hk) & is.finite(m_all)
  rho <- if (sum(finite) >= 3) cor(m_hk[finite], m_all[finite], method = "spearman") else NA_real_
  data.frame(
    cohort = ds, n_hk_genes = length(hk_genes), spearman_rho = rho,
    flag_low_rho = is.finite(rho) && rho < 0.7, stringsAsFactors = FALSE
  )
})

readr::write_tsv(centering_baseline, file.path(output_dir, "results", "centering_baseline_check.tsv"))
saveRDS(Xtilde_full, file.path(output_dir, "results", "Xtilde_full.rds"))
saveRDS(Xtilde, file.path(output_dir, "results", "Xtilde.rds"))

atlas_env$Xtilde_full <- Xtilde_full
atlas_env$Xtilde <- Xtilde
atlas_env$centering_baseline <- centering_baseline
save.atlas.checkpoint(atlas_env, 3)
```

```
## Saved stage03_atlas_env.rds
```

``` r
centering_baseline
```

```
##        cohort n_hk_genes spearman_rho flag_low_rho
## 1 Kurten_HNSC         75    0.7421550        FALSE
## 2  Puram_HNSC         75    0.8128033        FALSE
## 3   Choi_HNSC         75    0.6347274         TRUE
## 4   Bill_HNSC         75    0.4761342         TRUE
```
