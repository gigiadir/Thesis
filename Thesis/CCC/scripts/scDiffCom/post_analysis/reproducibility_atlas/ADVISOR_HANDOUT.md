# Cross-Cohort CCI Reproducibility Atlas — Advisor Handout

**Project:** HNSCC cell–cell interaction (CCI) reproducibility across four atlas subsets  
**Cohorts:** Kurten_HNSC, Puram_HNSC, Choi_HNSC, Bill_HNSC  
**Gene universe:** 392-gene intersection (DE-filtered malignant CCI)  
**CCI vocabulary:** 2,401 interactions (intersection mode)  
**Pipeline version:** `reproducibility_atlas/v1-v2`

---

## Executive summary

We ask whether genes show a **reproducible CCI logFC pattern** across four independent HNSCC scRNA-seq cohorts. For each gene, we build a vector of logFC values across 2,401 CCIs in each cohort, compare cross-cohort profile similarity using Spearman correlation, and convert that into a **ReproScore** (a percentile against all cross-gene matches). A **gene-shuffle null** (1,000 permutations) calibrates significance.

**Headline result:** Global atlas-level reproducibility is statistically significant (empirical **p = 0**). **52 genes** pass per-gene `shuffle_FDR < 0.05`. Leave-one-cohort-out analysis shows rankings are stable (min rank ρ = **0.81**). The pipeline passes all critical validation checks (effective sample size, null calibration, diagonal signal).

**Main caveat:** All four cohorts are HNSCC atlas subsets — not fully independent biological replicates. Results should be framed as cross-cohort consistency within a disease context, not independent replication.

---

## Data flow (Stages 0–5)

```
Stage 0  Load & inspect     → 4 cohorts, 392-gene intersection
Stage 1  CCI vocabulary     → J = 2,401 CCIs (intersection)
Stage 2  Tensor X            → X[cohort, gene, CCI] = mean LOGFC
Stage 3  Centering           → X̃ = column-wise z-score per cohort
Stage 4  ReproScore          → gene × gene Spearman → ReproScore
Stage 5  Nulls               → gene-shuffle null, shuffle_FDR, global test
```

---

## Stage 3 — Centering (why and how)

**Problem:** Cohorts may differ in overall LOGFC scale (batch, sample size, DE sensitivity).

**Solution:** For each cohort **separately**, for each CCI column **j**, z-score across the 392 genes:

$$\tilde{X}_{c,g,j} = \frac{X_{c,g,j} - \mu_{c,j}}{\sigma_{c,j}}$$

**What this is NOT:**
- Not centering across cohorts
- Not centering across CCIs within a gene (rows are not z-scored)

**Interpretation:** After centering, we compare **shapes** of CCI profiles (“is gene *g* relatively high in TGF-β CCIs and low in chemokine CCIs?”), not raw magnitudes.

**Validation:** Post-centering column means ≈ 0 and SDs ≈ 1 (`center_moments` PASS).

---

## Stage 4 — Gene × gene correlation and ReproScore

### The correlation matrix C

For one cohort pair (a, b), each gene **g** has a CCI profile vector of length 2,401 in each cohort. We compute:

$$C[g,h] = \text{Spearman}\big(\tilde{v}_a(g),\, \tilde{v}_b(h)\big)$$

| Entry | Meaning |
|-------|---------|
| **C[g, g]** (diagonal) | Does gene *g*'s CCI pattern in cohort A match **the same gene** in cohort B? |
| **C[g, h], g ≠ h** (off-diagonal) | Does gene *g* in A look like a **different** gene *h* in B? |

**High diagonal:** reproducible cross-cohort pattern for that gene.  
**High off-diagonal:** confounding — hard to distinguish self from cross matches.

**Validation:** mean(diagonal) − mean(off-diagonal) = **+0.17** (`diag_visible` PASS).

### The percentile score u_{a,b}(g)

For gene *g* in pair (a,b):

$$u_{a,b}(g) = \frac{1}{|H|}\sum_{h \neq g} \mathbf{1}\big[C[g,g] > C[g,h]\big]$$

= fraction of other genes *h* whose cross-match is **weaker** than *g*'s self-match.

| u_{a,b}(g) | Interpretation |
|------------|----------------|
| ≈ 1.0 | Self-pattern beats almost all cross-gene matches |
| ≈ 0.5 | Indistinguishable from background (no signal) |

### ReproScore (per gene)

$$\text{ReproScore}(g) = \text{mean of } u_{a,b}(g) \text{ over 6 cohort pairs}$$

**Stage 4 outputs** (`repro_scores.tsv`):

| Column | Meaning |
|--------|---------|
| `gene` | Gene symbol |
| `ReproScore` | Main score (0.5 = null-like; higher = more reproducible) |
| `R_self` | Mean diagonal Spearman (raw self-correlation strength) |
| `frac_pairs` | Fraction of pairs with u > 0.95 |

**Full universe (392 genes):** median ReproScore = **0.68**, range 0.26–0.95.

---

## Stage 5 — Gene-shuffle null and significance

**Null:** Shuffle gene labels **within each cohort** (breaks true gene identity while preserving cohort structure). Recompute ReproScore 1,000 times.

**Why null should center at 0.5:** Under shuffle, no gene should systematically beat cross-gene matches → ReproScore ≈ 0.5.

**Validation:** Null mean = **0.4999** (`null_centered` PASS).

### Outputs (`repro_scores_with_nulls.tsv`)

| Column | Meaning |
|--------|---------|
| `shuffle_p` | Empirical p-value: fraction of nulls ≥ observed ReproScore |
| `shuffle_FDR` | BH-adjusted `shuffle_p` |
| `shuffle_FDR_low_power` | Flag: per-gene FDR is exploratory with 4 cohorts |

### Global test (headline)

**Statistic:** mean(ReproScore) − 0.5 = **0.162**  
**Empirical one-sided p:** **0** (observed in extreme right tail of 1,000 nulls)

> **Recommendation:** Lead with the **global permutation test** for the atlas claim. Per-gene FDR is supplementary.

---

## Atlas membership (recommended rule)

**Do not use `reproscore_threshold = 0.95`** — no gene exceeds 0.95 (max = 0.950).

**Recommended (FDR-only):**

```
atlas_member := shuffle_FDR < 0.05
```

→ **52 genes** (ReproScore range 0.82–0.95 among these)

**Optional strength filter:**

```
atlas_member := shuffle_FDR < 0.05  AND  ReproScore > 0.85
```

→ **44 genes** (drops 8 borderline low-ReproScore FDR hits)

**Rank within atlas by ReproScore** (many genes tie at shuffle_p = 0).

| Tier | Rule | n |
|------|------|---|
| Core | FDR < 0.05 & ReproScore ≥ 0.90 | 11 |
| Extended | FDR < 0.05 & 0.85 ≤ ReproScore < 0.90 | 33 |
| Borderline | FDR < 0.05 & ReproScore < 0.85 | 8 |

**Top genes (by ReproScore among FDR-significant):** SMARCA2, HLA-A, ERBB3, CDC27, BCL3, FANCL, TYK2, HLA-C, ARID2, THBS1.

---

## Validation summary (QC pipeline)

Automated checks in `validation/` read pipeline outputs and write `validation/results/VERDICT.tsv`.

### Critical checks — all PASS

| Check | Result | What it proves |
|-------|--------|----------------|
| `eff_n` | PASS (median 237 CCIs/pair) | Enough data for stable Spearman |
| `diag_visible` | PASS (Δ = 0.17) | Self-correlation > cross-correlation |
| `null_centered` | PASS (0.4999) | Null is well-formed |
| `global_test` | PASS (p = 0) | Atlas-level signal is real |
| `loco_stability` | PASS (ρ = 0.81) | Not driven by one cohort |

### Caveats (WARN — not blockers)

| Check | Issue |
|-------|-------|
| `reproscore_vs_n` | Mild sparsity confounding (ρ = 0.24) |
| `batch_collapse` | Cohort batch effect was weak pre-centering |

### Key figures for presentation

| File | Shows |
|------|-------|
| `validation/results/batch_before_after.png` | PCA before/after centering |
| `validation/results/diagonal_heatmaps.png` | Diagonal signal in C matrices |
| `validation/results/global_null_figure.png` | Observed statistic vs null |
| `validation/results/loco_rank_cor.tsv` | Stability when dropping one cohort |

---

## Decision tree (is the result trustworthy?)

```
Effective-n adequate?  ──No──► Fix vocabulary (Stage 1)
        │
       Yes
        ▼
Null centered at 0.5?  ──No──► Fix shuffle logic (Stage 5)
        │
       Yes
        ▼
Global test significant?  ──No──► No atlas signal; revisit upstream
        │
       Yes
        ▼
LOCO stable?  ──No──► Cohort-driven; report which cohort
        │
       Yes
        ▼
Proceed to atlas interpretation (52 FDR-significant genes)
```

---

## Limitations to disclose

1. **Non-independent cohorts:** All four are HNSCC atlas subsets from related projects.
2. **Intersection vocabulary:** 2,401 CCIs — sparse NA structure (median ~72% NA per gene-cohort).
3. **Per-gene FDR:** Underpowered with 4 cohorts; global test is the primary inferential claim.
4. **ReproScore ≠ biology per se:** High score means cross-cohort CCI *pattern* consistency, not necessarily known cancer driver status.
5. **IDR not used for membership:** IDR concordance was weak in validation; membership based on shuffle FDR only.

---

## One-paragraph methods blurb (for slides/paper)

> For each of 392 intersection genes, we extracted mean malignant CCI logFC across 2,401 shared interactions in four HNSCC cohorts. Per-cohort CCI columns were z-scored across genes to remove scale effects. For each of six cohort pairs, we computed Spearman correlations between genes' CCI profiles, yielding a gene×gene matrix. ReproScore measures, for each gene, the fraction of other genes whose cross-cohort cross-match is weaker than the gene's self-match, averaged over pairs. Significance was assessed by 1,000 gene-label shuffle permutations per cohort. Global atlas reproducibility was significant (empirical p = 0). Fifty-two genes passed FDR < 0.05; leave-one-cohort-out rank correlations exceeded 0.8 for all cohort drops.

---

## File locations

| Artifact | Path |
|----------|------|
| ReproScore table | `results/repro_scores.tsv` |
| With nulls / FDR | `results/repro_scores_with_nulls.tsv` |
| Global null | `results/global_null.txt` |
| Validation verdict | `validation/results/VERDICT.tsv` |
| Pipeline README | `README.md` |

---

*Generated for thesis advisor review. Pipeline: `reproducibility_atlas/` on branch `cursor/reproducibility-atlas-stage0`.*
