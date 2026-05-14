# Technical Documentation: `HeadAndNeckAnalysis.Rmd`

**Author:** Adir Gigi | **Date:** 2024-09-29 | **Project:** Thesis — Figure 1 (Bulk)

---

## 1. Executive Summary

This analysis investigates the role of the **AXL receptor tyrosine kinase** in the head and neck squamous cell carcinoma (HNSCC) tumor microenvironment (TME). Using three independent single-cell RNA-seq (scRNA-seq) cohorts, the script pursues two overarching biological questions:

1. **Which cell types within HNSCC tumors most highly express AXL, and which genes co-vary with AXL within those cell types?** This is addressed by computing within-cell-type co-expression correlations and intersecting results across independent datasets to produce robust, cell-type-specific **AXL co-expression signatures**.
2. **In bulk RNA-seq data, do AXL expression signatures in specific cell types (particularly Macrophages and Malignant cells) correlate with immune exclusion or T-cell exhaustion phenotypes?** This links AXL activity to the broader immunosuppressive TME landscape.

The analysis culminates in **Figure 1C**, a patient-level scatter/line plot correlating per-patient AXL signature scores with T-cell abundance, CD4/CD8 ratios, T-cell exhaustion, and PD-1 (PDCD1) expression.

---

## 2. Pipeline Breakdown

### Stage 0 — Environment Setup
- Sources project-level constants (`constants.R`) and utility libraries (`Utils.R`, `Seurat.Utils.R`, `RNA.seq.Utils.R`, `CIBERSORT.R`).
- Loads two pre-built gene set databases: `GSDB.list.of.lists` (curated cell-type and pathway gene sets) and `GSDB.list.inclusive` (broader inclusive gene set collection, including Tirosh T-cell exhaustion signatures).

### Stage 1 — Data Loading & Initial Seurat Construction
Three scRNA-seq datasets are loaded from the **3CA (Tirosh lab) HNSCC collection**:
- **Kürten 2021** — UMI count matrix (`.mtx`)
- **Cillo 2020** — UMI count matrix (`.mtx`)
- **Puram 2017** — TPM matrix (`.mtx`; treated as pre-normalized)

Each is assembled into a `Seurat` object via a custom `load.and.create.seurat()` wrapper that ingests the expression matrix, barcode table, gene table, and pre-annotated metadata (including `cell_type` and `patient` fields).

### Stage 2 — Preprocessing & QC (`seurat.pipeline`)
Each Seurat object passes through a custom `seurat.pipeline()` wrapper. For Puram 2017 the `is_tpm = TRUE` flag is set, indicating the pipeline branches to skip library-size normalization (data is already TPM). The pipeline is presumed to include, at minimum:
- Normalization (log-normalization for UMI datasets)
- Highly variable gene selection
- PCA dimensionality reduction
- UMAP embedding
- Standard Seurat cell clustering

### Stage 3 — Imputation via ALRA
`RunALRA()` (from `SeuratWrappers`) is applied to all three objects, storing imputed values in a dedicated `"alra"` assay slot. **All downstream gene expression queries use the `alra` assay**, addressing the extreme sparsity/dropout characteristic of scRNA-seq data.

### Stage 4 — Initial Visualization (Cell Type & AXL/HAVCR2 Expression)
Produces exploratory plots to confirm data quality and motivate the AXL-centric analysis (detailed in Section 5).

### Stage 5 — AXL Co-Expression Correlation Analysis
The core biological inference step. For each dataset and each cell type, the script:
1. Subsets cells by `cell_type`.
2. Computes pairwise Spearman-rank correlations between AXL and all other genes using `scran::correlatePairs()` on the ALRA-imputed expression matrix.
3. Filters to retain gene pairs with **FDR < 0.1** and takes the **top 1,000 genes** ranked by `rho`.

This is run independently for: **Fibroblast, Macrophage, Malignant, T_cell** (and other cell types present in each dataset).

### Stage 6 — Cross-Dataset Signature Construction
The per-study AXL co-expression gene lists are **intersected across datasets** using `Reduce(intersect, ...)` to distill a consensus, study-replicated signature:

| Cell Type | Datasets Intersected |
|---|---|
| Macrophage | Kürten 2021 ∩ Cillo 2020 |
| Fibroblast | Kürten 2021 ∩ Puram 2017 |
| Malignant | Kürten 2021 ∩ Puram 2017 |
| T_cell | Kürten 2021 ∩ Puram 2017 |

> **Note:** Cillo 2020 is used for Macrophage but not for Fibroblast/Malignant/T_cell — likely because the Cillo dataset lacks those annotations or has insufficient cell numbers for those types.

A **parallel signature track** is constructed for CIBERSORTx deconvolution: using `RunPrestoAll()` (Presto fast Wilcoxon), the top 1,000 markers ranked by AUC are found per cell type per study, then intersected across studies (capped at **50 genes per cell type** for the signature matrix). The Dendritic cell signature is manually overridden to use only the Kürten ∩ Cillo intersection.

### Stage 7 — Bulk RNA-seq Deconvolution (CIBERSORTx)
- The HNSCC bulk RNA-seq dataset (`HNSC.Bulk.Dataset`) is exported as a tab-delimited count matrix.
- The Kürten 2021 signature matrix (pre-computed and loaded) is used as the reference panel.
- `CIBERSORT()` is called to infer cell-type proportions across bulk tumor samples.

### Stage 8 — Figure 1C: Patient-Level Correlation (Malignant & Macrophage)
Two parallel analysis blocks compute per-patient summaries of:
- **AXL signal** in Malignant or Macrophage cells (either raw expression quantiles or `AddModuleScore` using the consensus AXL signature).
- **Immune context** (T-cell, CD8, CD4, exhaustion signature module scores; PD-1 expression in T-cells) aggregated per patient.

These per-patient data frames are merged and visualized as scatter + regression plots.

### Stage 9 — Bulk Correlation Heatmap (Cross-Cell-Type)
A cross-validated heatmap is generated correlating, for each bulk sample:
- AXL signature AUCell scores (from `calc.Signatures.AUC`) for each cell type
- CIBERSORT-inferred cell type proportions

This produces a Pearson correlation matrix visualized as a labeled heatmap.

---

## 3. Key Methodologies

### R Packages

| Package | Role |
|---|---|
| `Seurat` | Single-cell data container, normalization, UMAP, clustering, module scoring, data fetching |
| `SeuratWrappers` | `RunALRA()` — low-rank imputation of scRNA-seq dropout |
| `scran` | `correlatePairs()` — rank-based gene-pair correlation with FDR control |
| `dplyr` | Per-patient data aggregation, filtering, joining |
| `ggplot2` | All visualizations (scatter, violin, line, tile/heatmap) |
| `reshape2` | `melt()` — wide-to-long transformation for multi-metric plotting |
| `VennDiagram` | Three-way Venn diagrams of gene set overlaps |
| Custom `CIBERSORT.R` | Local implementation of the CIBERSORT deconvolution algorithm |

### Statistical & Algorithmic Methods

- **ALRA Imputation:** Adaptively-thresholded Low Rank Approximation. Recovers biologically plausible zero values from the sparse scRNA-seq count matrix by assuming a low-rank underlying signal. Critical here because AXL is a low-to-moderate expression gene prone to dropout.
- **`scran::correlatePairs()`:** Spearman-rank correlation with permutation-based FDR estimation. Appropriate for non-normally distributed scRNA-seq data. Applied within cell-type subsets to avoid confounding by cell-type composition.
- **`Seurat::AddModuleScore()`:** Computes a per-cell composite score for a gene set by comparing average expression of the set against a random background of similar-expression genes (Tirosh et al. method). Used to score AXL and immune signatures.
- **`RunPrestoAll()` (Presto):** Fast Wilcoxon rank-sum test with AUC statistic for one-vs-rest differential expression. Used to find canonical cell-type marker genes.
- **`calc.Signatures.AUC()`:** Custom AUCell-style function ranking bulk samples by expression and computing area under the recovery curve for each signature. Top 20% of ranked genes are used.
- **`CIBERSORT`:** Linear support vector regression-based deconvolution of bulk RNA-seq into cell-type proportions using a reference signature matrix.
- **Pearson correlation** for comparing AXL AUCell scores vs. CIBERSORT proportions across patients.
- **Linear regression smoothing** (`geom_smooth(method = "lm")`) for patient-level scatter plots.

---

## 4. Visualization Catalog

### 4.1 UMAP Plots by Cell Type (3 plots)
| File | Description |
|---|---|
| `Puram2017-UMAP-cell-type.png` | UMAP embedding of Puram 2017 cells, colored by pre-annotated `cell_type`. Axes: UMAP1/UMAP2 (unitless). Signal: confirms cell-type separation; quality check of the Seurat pipeline. |
| `Cillo2020-UMAP-cell-type.png` | Same for Cillo 2020. |
| `Kurten2021-UMAP-cell-type.png` | Same for Kürten 2021. |

### 4.2 AXL Expression Violin Plots (4 plots)
| File | Description |
|---|---|
| `Kurten2021-AXL-by-cell-type.png` | Kürten 2021, **pre-ALRA** (default RNA assay). X: cell type, Y: normalized AXL expression. Signal: baseline AXL expression level per cell type before imputation. |
| `Kurten2021-AXL-by-cell-type-Alra.png` | Kürten 2021, **post-ALRA** (`alra` assay). Signal: effect of imputation on AXL signal recovery, particularly in low-expressing populations. |
| `Cillo2020-AXL-by-cell-type.png` | Cillo 2020, post-ALRA. |
| `Puram2017-AXL-by-cell-type.png` | Puram 2017, post-ALRA. |

### 4.3 HAVCR2 (TIM-3) Violin Plot (1 plot)
| File | Description |
|---|---|
| `Kurten2021-HAVCR2-by-cell-type.png` | Kürten 2021, post-ALRA. X: cell type, Y: HAVCR2 expression. Signal: motivated by a prior pan-cancer bulk analysis showing HAVCR2 is highly correlated with AXL and enriched in Macrophages and CAFs. Validates that co-expression extends to the single-cell level. |

### 4.4 Venn Diagram — T_cell AXL Signature (1 plot)
| File | Description |
|---|---|
| `T_cell_signatures_venn2.png` | Three-way Venn diagram of AXL-correlated genes (FDR < 0.1, top 1,000) within T-cells across all three studies. Signal: quantifies the degree of cross-study consensus in the T-cell AXL co-expression signature. The intersection forms the final T_cell AXL signature. |

### 4.5 Figure 1C — Malignant Cell AXL vs. Immune Metrics (2 plots, in-notebook)
Both plots share the same data; Kürten 2021 patients are the unit of observation.
- **X-axis:** Per-patient 99th-percentile AXL expression in Malignant cells (`quantile_0.99_AXL`).
- **Y-axis:** One of several immune metrics: mean T-cell signature score, mean CD8 score, mean CD4 score, mean T-cell exhaustion score, mean/median/75th-percentile PD-1 expression in T-cells.
- **Plot A:** Scatter + linear regression lines per metric (color-coded).
- **Plot B:** Line plot connecting patients per metric (shows rank ordering).
- **Signal:** Tests whether high AXL in tumor cells co-occurs with an immunosuppressed or exhausted T-cell compartment.

### 4.6 Figure 1C — Macrophage AXL Signature vs. Immune Metrics (2 plots, in-notebook)
Structurally identical to 4.5, but X-axis is:
- **X-axis:** Per-patient mean AXL module score in Macrophage cells (`mean_correlation_score_from_axl_signature`), derived from the consensus Macrophage AXL signature.
- **Signal:** Tests the same immunosuppression hypothesis from the Macrophage compartment rather than Malignant cells.

### 4.7 Cross-Cell-Type Correlation Heatmap (1 plot, in-notebook)
- **X-axis:** CIBERSORT-inferred cell-type proportion (one column per cell type).
- **Y-axis:** AUCell-scored AXL co-expression signature enrichment in bulk samples (one row per cell type).
- **Fill:** Pearson correlation coefficient (blue = negative, red = positive, white = 0).
- **Annotation:** Numeric correlation values printed in each tile.
- **Signal:** The off-diagonal entries reveal whether AXL activity in one cell type (e.g., Macrophage) tracks with the abundance of a different cell type (e.g., T-cell), illuminating cell-cell coordination in the TME. The diagonal tests whether AXL signature enrichment in a cell type correlates with that same cell type's abundance.
- **Caption text in code:** "Based on 469 patients: Pearson correlation between AXL cell-type signature scores and CIBERSORT-inferred cell type score across samples."

---

## 5. Data Flow

```
[scRNA-seq Raw Input]
  Kürten 2021 (.mtx + Cells.csv + Genes.txt + Meta-data.csv)
  Cillo 2020  (.mtx + Cells.csv + Genes.txt + Meta-data.csv)
  Puram 2017  (.mtx TPM + Cells.csv + Genes.txt + Meta-data.csv)
         │
         ▼
[load.and.create.seurat()]  ──→  Seurat objects with cell_type metadata
         │
         ▼
[seurat.pipeline()]  ──→  Normalized + Clustered + UMAP-embedded Seurat objects
  (is_tpm=TRUE for Puram)
         │
         ▼
[RunALRA()]  ──→  Seurat objects with additional "alra" assay slot
         │
         ├──→  [UMAP + Violin plots]  (Stage 4)
         │
         ▼
[correlatePairs() on alra assay, per cell type]
  ──→  AXL-correlated gene lists (FDR < 0.1, top 1000) per {dataset × cell type}
         │
         ▼
[Reduce(intersect) across datasets per cell type]
  ──→  AXL.signature.per.cell.type  (consensus cross-study signatures)
         │
         ├──→  [Venn Diagram]  (Stage 6)
         │
         ▼
[RunPrestoAll() on alra assay, per cell type]
  ──→  Top 1000 markers by AUC per {dataset × cell type}
         │
         ▼
[Reduce(intersect) + truncate to top 50]
  ──→  Cell.Type.Signatures.For.Matrix  (for CIBERSORTx reference panel)
         │
         ▼
[Bulk RNA-seq Input]  HNSC.Bulk.Dataset$counts  (469 patients)
         │
         ├──→  [CIBERSORT()]  ──→  Bulk.Deconvolution  (cell-type proportions per patient)
         │
         ├──→  [calc.Signatures.AUC()]  ──→  AXL AUCell scores per cell type per patient
         │
         ▼
[Pearson correlation matrix]
  ──→  Correlation Heatmap  (AXL AUCell × CIBERSORT proportions)

[Figure 1C — Kürten 2021 patient-level analysis]
  AddModuleScore(AXL signatures)  ──→  Per-cell scores
  FetchData(AXL, PDCD1)          ──→  Raw expression values
  group_by(Patient) + summarize() ──→  Per-patient statistics
  inner_join()                    ──→  Merged patient table
         │
         ▼
  [Scatter + Line plots: AXL signal vs. immune metrics]
```

---

## 6. Technical Observations

### Hard-Coded Thresholds & Parameters

| Parameter | Value | Location | Rationale |
|---|---|---|---|
| `FDR` threshold for `correlatePairs` | `< 0.1` | `find_correlated_pairs_per_cell_type()` | Permissive FDR to retain a broad list before cross-dataset intersection; the intersection step provides the true statistical stringency |
| `top_n` correlated pairs | `1,000` | `find_correlated_pairs_per_cell_type()` | Captures a broad co-expression neighborhood around AXL before cross-study filtering |
| `n_markers` (Presto) | `1,000` | `find.markers.by.cell.type()` | Same philosophy — broad initial pool to be refined by intersection |
| `SIGNATURE_MAX_LENGTH` | `50` | Signature matrix construction | Hard ceiling on genes per cell type in the CIBERSORTx reference panel; prevents overfitting |
| `top.n.percent.of.genes.in.ranking.used` (AUCell) | `20` | `calc.Signatures.AUC()` call | Standard AUCell parameterization — uses top 20% of genes ranked by expression per sample |
| AXL quantile used for Figure 1C X-axis | `0.99` | Figure 1C (Malignant block) | Among multiple quantiles computed (0.95–0.9995), 99th percentile was selected, likely to capture "high-AXL" cell signal while excluding extreme outliers |
| PD-1 quantile summary | `0.75` (in addition to mean/median) | Figure 1C | T-cell PD-1 (PDCD1) is summarized at three levels; 75th percentile captures "high-expressing" T-cells |

### Analytical Design Choices

- **ALRA used for all gene expression queries** — both for the correlation analysis and the violin plots. This is a deliberate choice to recover biologically meaningful signal from dropout. The pre/post ALRA violin plots for AXL serve as an explicit validation of this choice.
- **Dataset-specific intersection strategy for signatures:** The Macrophage signature uses Kürten ∩ Cillo (Puram excluded), while Fibroblast/Malignant/T_cell use Kürten ∩ Puram. This asymmetry is likely driven by the cell types present in each study's metadata.
- **Dendritic signature is manually overridden** after the automated loop, using only Kürten ∩ Cillo. This suggests the automated loop produced an unsatisfactory result (perhaps Puram 2017 had no Dendritic annotation), requiring manual correction.
- **`CIBERSORTx.Construct.Signature.Matrix` is a stub (`#TODO: IMPLEMENT`)** — the signature matrix was loaded from a pre-saved `.RData` object (`Kurten.2021.Signature.Matrix`), indicating this construction step was completed outside the notebook or in a prior session.
- **Figure 1C is implemented twice** — once for Malignant cells (using raw 99th-percentile AXL expression) and once for Macrophage cells (using an AXL module score). This duality tests whether the AXL-immunosuppression relationship is cell-autonomous (tumor-intrinsic) or TME-mediated (macrophage-driven).
- **469 patients** are referenced in the heatmap caption, establishing the bulk cohort size as the TCGA HNSC dataset (most likely TCGA-HNSC, given the `HNSC.Bulk.Dataset` naming convention).

---

*Documentation generated from `HeadAndNeckAnalysis.Rmd` — 2026-05-02.*
