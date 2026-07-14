# Cross-Cohort CCI Reproducibility Atlas — Advisor Handout

**Project:** HNSCC cell–cell interaction (CCI) reproducibility across four atlas subsets  
**Cohorts:** Kurten_HNSC, Puram_HNSC, Choi_HNSC, Bill_HNSC  
**Gene universe:** 392-gene intersection (DE-filtered malignant CCI)  
**CCI vocabulary:** 2,401 interactions (intersection mode)  
**Primary statistic:** **R_g** — per-gene cross-cohort reproducibility score  
**Pipeline version:** `reproducibility_atlas/v1-v2`

---

## Executive summary

We ask whether genes show a **reproducible CCI logFC pattern** across four independent HNSCC scRNA-seq cohorts. For each gene we build a vector of logFC values across 2,401 CCIs in each cohort, remove the cohort-common signal by subtracting the cohort "mean-gene" profile, and measure cross-cohort agreement with Spearman correlation.

The per-gene score is **R_g = the median of the gene's ≤6 same-gene pairwise cross-cohort Spearman correlations** (one per cohort pair). Significance is assessed against a **cross-gene pairing null** (gene *g* in cohort A vs a *different* gene in cohort B), calibrated per gene with **Extreme Value Theory / Generalized Pareto Distribution (GPD) tail extrapolation** (Knijnenburg et al. 2009).

**Headline result:** Global atlas-level reproducibility is significant (mean R_g = **0.195**, empirical **p = 0.001**, EVT global **p ≪ 1e-300**). **41 genes** pass per-gene `evt_FDR < 0.05` (18 strong, R_g ≥ 0.5; 23 moderate, 0.3 ≤ R_g < 0.5). All 392 genes have full 6/6 cohort-pair coverage.

**Main caveats:**
1. All four cohorts are HNSCC atlas subsets — not fully independent biological replicates.
2. **Leave-one-cohort-out (LOCO) now FAILs the strict threshold** (min rank ρ = **0.76** when Kurten_HNSC is dropped, vs the ρ > 0.8 gate). The gene ranking is *moderately* cohort-sensitive; Kurten_HNSC has the largest single-cohort influence. Report the atlas as cross-cohort consistency, and disclose the LOCO sensitivity.

---

## Old method vs new method (what changed and why)

The pipeline previously scored genes with a **percentile ReproScore** against a **gene-label shuffle** null. The advisor's methodology replaces this with a **correlation-scale R_g** and a **cross-gene pairing** null with EVT tail calibration. ReproScore is retained only as a **supplementary diagnostic**.

| Step | Old method | New method (advisor R_g) | Why the change |
|------|-----------|--------------------------|----------------|
| **Centering (Stage 3)** | Per-cohort, per-CCI **z-score** across genes (subtract mean *and* divide by SD) | Per-cohort **mean-gene subtraction** only (subtract the cohort's mean profile across genes per CCI; no SD scaling) | Removes the cohort-common signal without distorting the rank structure Spearman relies on; SD scaling is redundant under rank correlation |
| **Per-gene statistic (Stage 4)** | **ReproScore** = fraction of *other* genes whose cross-match is weaker than gene *g*'s self-match, averaged over 6 pairs (a percentile in [0,1], null ≈ 0.5) | **R_g** = **median** of the gene's ≤6 same-gene pairwise Spearman correlations (a correlation in [−1,1], null ≈ 0) | A direct, interpretable effect size (a correlation) instead of a rank-percentile; median is robust to one discordant cohort pair |
| **Partial coverage** | Implicit | Explicit `n_pairs_computable`; a pair needs ≥ `min_cci_overlap` (10) jointly-observed CCIs, and a gene needs ≥ `atlas_min_pairs` (1) computable pairs | Makes missing-data handling transparent and auditable |
| **Null model (Stage 5)** | **Shuffle gene labels** within each cohort, recompute ReproScore ×1,000 | **Cross-gene pairing:** for each pair, compare gene *g* to a *random different* gene (resample off-diagonal entries of the precomputed correlation matrices) ×1,000 | Directly tests "does the *same* gene agree more than *different* genes?"; reuses the already-computed matrices, so it is far cheaper |
| **Null center** | ≈ **0.5** (percentile scale) | ≈ **0** (correlation scale); observed −2e-04 | Correct baseline for a correlation statistic |
| **p-value / tail** | Empirical p = (#null ≥ obs + 1)/(N + 1); floored at 1/1001 | **Hybrid empirical + GPD**: empirical in the bulk; when observed lands in the extreme tail (< 10 null exceedances) fit a **GPD** to the top exceedances and read the p-value off the fitted tail, with an **Anderson–Darling** goodness-of-fit gate | 1,000 permutations cannot resolve p-values below ~1e-3; GPD extrapolation gives calibrated tail p-values "with fewer permutations" (Knijnenburg 2009) |
| **Atlas membership (Stage 6)** | `shuffle_FDR < 0.05` (52 genes) | `evt_FDR < 0.05 AND n_pairs_computable ≥ 1` (41 genes), tiered by R_g strength | FDR now rests on calibrated tail p-values; strength tier adds an interpretable effect-size layer |

**Bottom line:** the new method reports an **effect size you can defend** (a cross-cohort correlation) with a **null that directly matches the biological question** and a **calibrated tail** for the strongest genes.

---

## Data flow (Stages 0–7)

```
Stage 0  Load & inspect     → 4 cohorts, 392-gene intersection
Stage 1  CCI vocabulary     → J = 2,401 CCIs (intersection); full per-cohort vocab kept for centering
Stage 2  Tensor             → X_full (all CCIs, for centering) + X (J-aligned, for scoring)
Stage 3  Centering          → subtract cohort mean-gene profile on X_full → X̃, subset to J
Stage 4  R_g                → per cohort pair: gene×gene Spearman → R_g = median of self-correlations
Stage 5  Null + EVT         → cross-gene pairing null (1,000×) → empirical_p, evt_p, evt_FDR, global test
Stage 6  Atlas assembly     → atlas_member = evt_FDR < fdr_threshold; strength tiers by R_g
Stage 7  Diagnostics        → modules, heatmaps, coverage, EVT method breakdown
```

---

## Stage 0 — Load and inspect (gene universe)

**Problem:** Each cohort has its own scDiffCom per-gene RDS files. We need a **shared gene list** and **comparable malignant CCI tables** before building cross-cohort tensors.

**Solution:**

1. **Discover genes on disk** per cohort from `base_results_dir` (one RDS per gene × cohort).
2. **Gene universe** = intersection across all four cohorts → **392 genes** (all must be present in every cohort).
3. **Load + filter** malignant CCI rows from each gene's `cci_table_detected`:
   - Tumor involved: `EMITTER_CELLTYPE` or `RECEIVER_CELLTYPE` ∈ {Tumor}
   - DE only: `IS_CCI_DE == TRUE`
   - Finite `LOGFC`
4. **Optional:** drop CCIs involving unknown cell types (`filter_unknown_celltypes`).

**Stage 0 outputs:** `gene_universe.tsv` (392 genes), `genes_missing_per_cohort.tsv` (sanity, 0 missing), `inspect_report.txt`, `stage00_atlas_env.rds`.

**Interpretation:** Stage 0 fixes *which genes* enter the atlas. Everything downstream is conditional on the 392-gene intersection and DE-filtered malignant CCIs.

---

## Stage 1 — CCI vocabulary (dual: full for centering, J for scoring)

**Problem:** Genes do not share the same set of detected CCIs across cohorts. Spearman correlation for scoring requires a **common CCI index J**; but *centering* must see the full CCI signal of each cohort to correctly estimate the cohort-common profile.

**Solution — two tiers:**

1. **Full per-cohort vocabulary** (`cci_by_cohort`): every DE CCI observed in a cohort. Used *only* for Stage 3 centering, so the "mean-gene" profile is estimated on the complete signal.
2. **Intersection J** = CCIs present in **all** cohorts → **2,401** (`vocab_mode: intersection`). Used for all *scoring* (Stages 4–6) so profile vectors are aligned across cohorts.

Support loss and a per-cohort NA mask are recorded.

**Why dual:** Centering on the full vocabulary avoids leaking an intersection-induced bias into the centered values; restricting *scoring* to J keeps cross-cohort correlations well-defined.

**Stage 1 outputs:** `vocab_report.tsv`, `vocab_per_cohort.tsv`, `J.rds`, `vocab_na_mask.rds`.

**Validation:** Median effective-n = **237** pairwise-complete CCIs per gene (`eff_n` PASS). Vocabulary spread across celltype pairs (`vocab_comp` INFO, 0.125).

---

## Stage 2 — Tensor (gene × CCI logFC profiles)

**Problem:** Convert per-gene tables into **fixed-shape numeric tensors** suitable for correlation.

**Solution:** For each cohort **c** build two matrices:

$$X^{\text{full}}_{c}[g,j] = \text{mean}(\text{LOGFC}) \text{ over the cohort's full CCI vocabulary}$$
$$X_{c}[g,j] = \text{same, restricted to the 2,401 intersection CCIs } J$$

- Rows = genes (fixed order from Stage 0); columns = CCIs.
- **NA** when gene *g* has no DE row for CCI *j* in cohort *c*.
- `X_full` (native columns) feeds centering; `X` (J-aligned, 392 × 2,401) feeds scoring/QC.

**Sparsity:** Mean NA fraction ≈ **72%** per gene–cohort (`na_map` INFO) — expected, handled by pairwise-complete Spearman.

**Stage 2 outputs:** `X_full.rds`, `X.rds`, `tensor_na_fraction.tsv`, `na_mask.rds`, `stage02_atlas_env.rds`.

**Validation:** Cohort LOGFC IQR ratio = **1.19** (`logfc_scale` INFO); no duplicate CCI collapse (`dup_collapse` INFO 0).

---

## Stage 3 — Centering (mean-gene subtraction)

**Problem:** All genes in a cohort share a common CCI response component (the cohort's average profile). Left in, it inflates cross-gene similarity and masks gene-specific patterns.

**Solution:** For each cohort **separately**, subtract the **mean-gene profile** — the average across genes for each CCI column — from every gene's vector:

$$\tilde{X}_{c,g,j} = X^{\text{full}}_{c,g,j} - \mu_{c,j}, \qquad \mu_{c,j} = \frac{1}{|\text{genes}|}\sum_{g} X^{\text{full}}_{c,g,j}$$

Centering is computed on **`X_full`**, then the result is **subset to the intersection J** → `Xtilde`.

**What this is NOT:**
- Not a z-score — **mean subtraction only**, no division by SD (rank correlation is scale-invariant).
- Not centering across cohorts or across CCIs within a gene.

**Interpretation:** After centering we compare **gene-specific deviations** from the cohort-common CCI profile.

**Validation:** Post-centering column means ≈ 0 (`center_moments` PASS, max|mean| = 1.8e-16). The cohort-common signal is real and worth removing (`cohort_means` INFO −0.023; `batch_collapse` PASS — silhouette drops after centering).

---

## Stage 4 — Gene × gene correlation and R_g

### The correlation matrix C (per cohort pair)

For one cohort pair (a, b), each gene **g** has a centered CCI profile vector of length 2,401 in each cohort:

$$C_{ab}[g,h] = \text{Spearman}\big(\tilde{v}_a(g),\, \tilde{v}_b(h)\big) \quad\text{(pairwise-complete; NA if } < 10 \text{ shared CCIs)}$$

| Entry | Meaning | Feeds |
|-------|---------|-------|
| **C[g, g]** (diagonal) | Does gene *g*'s CCI pattern in A match **the same gene** in B? | **R_g** |
| **C[g, h], g ≠ h** (off-diagonal) | Does gene *g* in A look like a **different** gene *h* in B? | **the null** |

**Validation:** mean(diagonal) − mean(off-diagonal) = **+0.19** (`diag_visible` PASS; `raw_gap` PASS).

### R_g (per gene) — the primary statistic

Collect gene *g*'s diagonal correlation across the (up to) 6 cohort pairs, then take the **median**:

$$R_g = \text{median}\big\{\, C_{ab}[g,g] : \text{cohort pairs } (a,b) \text{ with } \ge 10 \text{ shared CCIs} \,\big\}$$

- Range [−1, 1]; **higher = more reproducible**. Null ≈ 0.
- `Rg_mean` (mean instead of median) is stored alongside for comparison.
- `n_pairs_computable` records how many of the 6 pairs contributed.

**Full universe (392 genes):** median R_g = **0.207**, mean = **0.195**, range **−0.38 to 0.68**. All 392 genes have **6/6** computable pairs.

### ReproScore (supplementary only)

The old percentile statistic is still computed for comparison: `ReproScore` (median 0.67, range 0.25–0.95), `R_self` (mean diagonal), `frac_pairs` (fraction of pairs with percentile > 0.95). **Not used for membership.**

**Stage 4 outputs:** `repro_scores.tsv`, `pairwise_rho.tsv`, `cor_pair_matrices.rds`.

---

## Stage 5 — Cross-gene null and EVT calibration

**Null (cross-gene pairing):** For each cohort pair, replace gene *g*'s self-correlation `C[g,g]` with its correlation to a **randomly chosen different gene** `C[g,g′]` (resampled from the off-diagonal of the already-computed matrices), then aggregate to a null R_g. Repeat **1,000×**.

**Why this null:** it asks exactly "does the *same* gene agree across cohorts more than *different* genes do?" — and because it resamples precomputed matrix entries, it needs no re-correlation.

**Why the null centers at 0:** unrelated gene pairs have no systematic rank agreement → null R_g ≈ 0. **Validation:** null mean = **−2e-04** (`null_centered` PASS, expected [−0.05, 0.05]); robust to NA structure (`na_exchangeable` PASS).

### EVT / GPD tail calibration (Knijnenburg et al. 2009)

1,000 permutations cannot resolve p-values below ≈ 1e-3, yet the strongest genes sit far in the tail. Per gene:

- **Bulk (≥ 10 null values ≥ observed):** use the empirical p-value (reliable).
- **Tail (< 10 exceedances):** fit a **Generalized Pareto Distribution** to the top exceedances of the null (Pickands–Balkema–de Haan theorem → the tail of almost any distribution is GPD), gated by an **Anderson–Darling** goodness-of-fit test (parametric bootstrap). If GOF rejects, shrink the exceedance count and refit; if no fit is acceptable, fall back to empirical. The p-value is then read off the fitted tail.

This is the "fewer permutations, more accurate P-values" recipe (Knijnenburg et al. 2009, *Bioinformatics* 25(12):i161–i168).

### Global test (headline)

**Statistic:** mean(R_g) over 392 genes = **0.195**  
**Empirical one-sided p:** **0.001** (top of 1,000 nulls; null mean ≈ 0, SD ≈ 0.009)  
**EVT global p:** **≪ 1e-300** (GPD; GOF p = 0.39 → adequate fit)

### Per-gene EVT method usage (Stage 7 diagnostic)

Of 392 genes: **320 empirical**, **72 GPD** (i.e. 72 genes were extreme enough to need tail extrapolation).

**Stage 5 outputs:** `repro_scores_with_nulls.tsv`, `null_reproscore_matrix.rds` (392 × 1,000), `global_null.txt`, `stage05_atlas_env.rds`.

---

## Stage 6 — Atlas assembly (membership rule)

**Membership:**

```
atlas_member := evt_FDR < fdr_threshold (0.05)  AND  n_pairs_computable >= atlas_min_pairs (1)
```

→ **41 genes** (`evt_FDR < 0.05`). Genes are ranked within the atlas by **R_g**.

**Strength tiers** (effect-size layer on top of significance):

| Tier | Rule | n |
|------|------|---|
| strong | atlas_member & R_g ≥ 0.5 | 18 |
| moderate | atlas_member & 0.3 ≤ R_g < 0.5 | 23 |
| weak | atlas_member & R_g < 0.3 | 0 |
| not_member | evt_FDR ≥ 0.05 | 351 |

**Top genes (by R_g among members):** CDC27, SMARCA2, BCL3, NCOA3, TYK2, BIRC3, FANCL, NFKB2, TRAF2, MAPK1, HLA-C, TCEB1 — an NF-κB / DNA-damage / receptor-signaling flavored set.

> `reproscore_threshold = 0.95` in the config is **legacy** and unused for membership.

**Stage 6 outputs:** `atlas.csv`, `atlas.parquet`, `stage06_atlas_env.rds`.

---

## Atlas column dictionary (`atlas.csv`)

One row per gene (392 rows). Columns, in order:

| Column | Type | Meaning |
|--------|------|---------|
| `gene` | chr | Gene symbol |
| **`Rg`** | num | **Primary statistic** — median of the gene's ≤6 same-gene cross-cohort Spearman correlations. Range [−1,1]; higher = more reproducible; null ≈ 0 |
| `Rg_mean` | num | Same as `Rg` but using the **mean** instead of the median (sensitivity check) |
| `n_pairs_computable` | int | How many of the 6 cohort pairs had ≥ 10 shared CCIs and contributed to `Rg` (here: 6 for all genes) |
| `ReproScore` | num | **Supplementary** old percentile statistic (0.5 = null-like; higher = more reproducible). Not used for membership |
| `R_self` | num | Mean diagonal Spearman (raw self-correlation strength across pairs) |
| `frac_pairs` | num | Fraction of pairs whose ReproScore percentile exceeds 0.95 (supplementary) |
| `empirical_p` | num | One-sided permutation p from the cross-gene null: (#null R_g ≥ observed + 1)/(N + 1); floor ≈ 1/1001 |
| **`evt_p`** | num | **Calibrated p-value** — empirical in the bulk, GPD tail-extrapolated in the extreme tail |
| `evt_method` | chr | Which estimator produced `evt_p`: `empirical`, `gpd`, `empirical_smallN`, or `empirical_fallback` |
| `evt_n_exc` | int | Number of tail exceedances used for the GPD fit (NA if empirical) |
| `evt_xi` | num | GPD shape parameter ξ of the fitted tail (NA if empirical). ξ < 0 ⇒ bounded (light) tail |
| `evt_gof_p` | num | Anderson–Darling goodness-of-fit p for the GPD fit; large ⇒ GPD is adequate (NA if empirical) |
| `empirical_FDR` | num | Benjamini–Hochberg FDR over `empirical_p` |
| **`evt_FDR`** | num | **BH FDR over `evt_p`** — the value that drives atlas membership |
| `shuffle_p` | num | Legacy gene-shuffle p (old method), kept for back-comparison |
| `shuffle_FDR` | num | Legacy BH FDR over `shuffle_p` |
| `shuffle_FDR_low_power` | lgl | Legacy flag: per-gene FDR is exploratory with 4 cohorts |
| **`atlas_member`** | lgl | `evt_FDR < 0.05 & n_pairs_computable ≥ 1` (41 TRUE) |
| **`strength_tier`** | chr | `strong` / `moderate` / `weak` / `not_member` (by `Rg` among members) |

---

## Validation summary (QC pipeline)

Automated checks in `validation/` read pipeline outputs and write `validation/results/VERDICT.tsv`. Overall: **22 PASS, 5 WARN, 14 INFO, 1 FAIL**.

### Critical checks

| Check | Result | What it proves |
|-------|--------|----------------|
| `eff_n` | **PASS** (median 237 CCIs/pair) | Enough data for stable Spearman |
| `center_moments` | **PASS** (max abs mean ≈ 2e-16) | Mean-gene subtraction worked (means ≈ 0) |
| `diag_visible` | **PASS** (Δ = 0.19) | Same-gene correlation > different-gene correlation |
| `null_centered` | **PASS** (−2e-04) | Cross-gene null correctly centered at 0 |
| `global_test` | **PASS** (p = 0.001) | Atlas-level signal is real |
| `loco_stability` | **FAIL** (min ρ = 0.76) | Ranking is **moderately cohort-sensitive** (see below) |

### LOCO detail (leave-one-cohort-out rank correlation vs full atlas)

| Dropped cohort | rank ρ vs full |
|----------------|----------------|
| Kurten_HNSC | **0.76** (worst) |
| Choi_HNSC | 0.78 |
| Bill_HNSC | 0.80 |
| Puram_HNSC | 0.81 |

All four are close to the 0.8 gate; **Kurten_HNSC** has the most influence. Interpretation: the atlas is broadly stable but not immune to single-cohort removal — disclose this and consider it when making per-gene claims.

### Caveats (WARN — not blockers)

| Check | Issue |
|-------|-------|
| `Rg_vs_n` | Mild sparsity confounding (ρ = 0.25 between R_g and pairwise-complete counts) |
| `controls` | Default positive/negative controls (VEGFA/ACTB) do not separate strongly — treat biological priors cautiously |

### Key figures for presentation

| File | Shows |
|------|-------|
| `validation/results/batch_before_after.png` | PCA before/after centering |
| `validation/results/diagonal_heatmaps.png` | Diagonal (same-gene) signal in C matrices |
| `validation/results/global_null_figure.png` | Observed mean R_g vs cross-gene null |
| `validation/results/loco_rank_cor.tsv` | Rank stability when dropping each cohort |

---

## Decision tree (is the result trustworthy?)

```
Effective-n adequate?  ──No──► Fix vocabulary (Stage 1)
        │  Yes
        ▼
Null centered at 0?    ──No──► Fix cross-gene null (Stage 5)
        │  Yes
        ▼
Global test significant? ─No──► No atlas signal; revisit upstream
        │  Yes
        ▼
LOCO stable (ρ>0.8)?   ──No──► Cohort-driven; report which cohort ◄── CURRENTLY HERE (Kurten ρ=0.76)
        │  Yes
        ▼
Proceed to atlas interpretation (41 evt_FDR-significant genes)
```

Current status: significant global signal and 41 FDR-significant genes, but the LOCO gate is not fully met — frame per-gene claims with the Kurten sensitivity caveat.

---

## Limitations to disclose

1. **Non-independent cohorts:** All four are HNSCC atlas subsets from related projects.
2. **LOCO sensitivity:** Dropping Kurten_HNSC lowers rank ρ to 0.76 (< 0.8 gate); the ranking is moderately cohort-dependent.
3. **Intersection vocabulary:** 2,401 CCIs — sparse NA structure (median ~72% NA per gene-cohort).
4. **Mild sparsity confounding:** R_g correlates weakly (ρ ≈ 0.25) with the number of usable CCIs.
5. **R_g ≠ biology per se:** High R_g means cross-cohort CCI *pattern* consistency, not known cancer-driver status; default controls did not separate cleanly.

---

## One-paragraph methods blurb (for slides/paper)

> For each of 392 intersection genes we extracted mean malignant CCI logFC across 2,401 shared interactions in four HNSCC cohorts. Per cohort, we subtracted the mean-gene CCI profile to remove the cohort-common signal (mean subtraction; no scaling). For each of six cohort pairs we computed Spearman correlations between genes' CCI profiles, forming a gene×gene matrix whose diagonal measures same-gene cross-cohort agreement. Each gene's reproducibility score R_g is the median of its ≤6 diagonal correlations. Significance was assessed against a cross-gene pairing null (same gene vs a random different gene) over 1,000 permutations, calibrated per gene by GPD tail extrapolation with an Anderson–Darling goodness-of-fit gate (extreme value theory; Knijnenburg et al. 2009). Global reproducibility was significant (mean R_g = 0.195; empirical p = 0.001; EVT p ≪ 1e-300). Forty-one genes passed evt_FDR < 0.05 (18 with R_g ≥ 0.5). Leave-one-cohort-out rank correlations ranged 0.76–0.81, indicating broadly stable but moderately cohort-sensitive rankings.

---

## File locations

| Artifact | Path |
|----------|------|
| Gene universe | `results/gene_universe.tsv` |
| CCI vocabulary | `results/J.rds`, `results/vocab_report.tsv` |
| Tensors | `results/X_full.rds` (centering), `results/X.rds` (scoring) |
| Centered tensor | `results/Xtilde.rds` |
| Pairwise correlations | `results/cor_pair_matrices.rds`, `results/pairwise_rho.tsv` |
| R_g table | `results/repro_scores.tsv` |
| With nulls / FDR | `results/repro_scores_with_nulls.tsv` |
| Global null | `results/global_null.txt` |
| **Atlas** | `results/atlas.csv`, `results/atlas.parquet` |
| Diagnostics | `results/diagnostics/` (modules, coverage, EVT breakdown, heatmap) |
| Validation verdict | `validation/results/VERDICT.tsv` |
| Pipeline README | `README.md` |

---

*Generated for thesis advisor review. Pipeline: `reproducibility_atlas/` (advisor R_g methodology).*
