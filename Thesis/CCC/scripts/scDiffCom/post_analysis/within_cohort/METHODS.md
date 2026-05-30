# Within-Cohort Differential Communication Phenotypes

Methods documentation for the within-cohort post-analysis module (dCCC). Suitable for adaptation in a thesis Methods section.

---

## Overview

To characterize how perturbation of individual cancer-associated genes reshapes intercellular signaling **within a single head and neck squamous cell carcinoma (HNSC) cohort**, we developed a within-cohort post-analysis module downstream of **scDiffCom**. For each cohort and each gene in a predefined panel, scDiffCom compares cell–cell communication (CCC) between patients with **high** versus **low** malignant expression of that gene. The within-cohort module aggregates these gene-specific differential CCC profiles into comparable representations and visual summaries, which we refer to as **within-cohort differential communication phenotypes**.

The pipeline is implemented in R (`within_cohort/index.Rmd` and helper functions in `post_analysis/within_cohort/`) and is run **independently per cohort** (e.g. Kurten_HNSC, Bill_HNSC, Puram_HNSC, Choi_HNSC). It is distinct from the cross-cohort post-analysis, which compares phenotypes across datasets after batch correction.

---

## Upstream inputs: scDiffCom per-gene objects

For each cohort \(d\) and each cancer gene \(g\) in the analysis panel, we loaded a scDiffCom result object produced by the upstream pipeline (`{g}_{d}_scDiffCom.rds` from `split-by-rank-genes-v2`). Each object contains:

- **`cci_table_detected`**: differentially detected cell–cell interactions (CCIs), with ligand–receptor identifier **LRI**, emitter and receiver cell types, log-fold change **LOGFC** of the interaction score between HIGH and LOW expression groups, regulation direction (**REGULATION**: UP/DOWN), and flag **`IS_CCI_DE`**.
- **`cci_table_raw`**: the broader set of tested interactions involving the malignant compartment (used for GO enrichment background; Module 4).

Patient stratification into HIGH/LOW groups follows the rank-based preprocessing used in the scDiffCom pipeline (pseudobulk expression tertiles per gene per cohort). The within-cohort module does not re-run scDiffCom; it operates on stored outputs.

---

## Module 1: Data loading and filtering

**Loading.** All per-gene scDiffCom objects for cohort \(d\) were loaded into a list \(\mathcal{G}_d = \{ \text{scDiffCom}^{(g)} : g \in \mathcal{G} \}\), where \(\mathcal{G}\) is the gene panel available for that cohort.

**Malignant-focused differential interactions.** From each `cci_table_detected`, we retained rows satisfying:

1. The malignant compartment is involved: emitter or receiver cell type is **Tumor** (HNSC annotation).
2. `IS_CCI_DE == TRUE` (interaction significantly differential between HIGH and LOW).
3. `LOGFC` is finite.

This yields, for each gene \(g\), a table \(T^{(g)}_d\) of malignant-involved DE interactions.

**Cell-type quality filter.** We excluded interactions where either emitter or receiver matched unknown, equivocal, or other ambiguous labels (case-insensitive pattern matching). Original tables were preserved; filtered tables were used for all downstream steps.

**Output of Module 1:** per-cohort list of filtered tables \(\{ T^{(g)}_d \}_{g \in \mathcal{G}}\), used as the common input to Modules 2–3 (and optionally 4).

---

## Module 2: Gene-of-interest dissection — emitter–receiver resolved dCCC

**Purpose.** To inspect how perturbation of a single gene (default: **AXL**) redistributes differential signaling across **cell-type axes** in one cohort, without collapsing across emitter–receiver pairs.

**Unit of observation.** Each row in \(T^{(\text{AXL})}_d\) is one differential LRI at a specific **emitter → receiver** cell-type pair. We define:

\[
\text{ER\_pair} = (\text{EMITTER\_CELLTYPE}) \rightarrow (\text{RECEIVER\_CELLTYPE})
\]

**Visualization.** For each ER_pair, we plotted the distribution of **LOGFC** across all differential LRIs in that pair (boxplot + jittered points). A horizontal reference line at LOGFC = 0 indicates no change. The subtitle reports cohort identity and **\(N\)**, the total number of interaction points (rows) in the plot.

**LRI annotation.** To highlight extreme interactions without overcrowding the figure, we labeled LRIs only for the **top 8 ER_pairs by LogFC spread** (max(LOGFC) − min(LOGFC) within the pair). For each selected pair, we labeled the LRI achieving the **minimum** and **maximum** LOGFC (up to 16 labels total). Labels were placed with directional text repulsion along the y-axis to avoid overlap with boxplot whiskers. Axis limits were set to the data range plus 4% padding, accounting for label offsets.

**Interpretation.** This module describes a **gene-specific, spatially resolved dCCC phenotype**: which cell-type channels gain or lose ligand–receptor signaling when AXL (or another gene of interest) is perturbed in the malignant compartment.

---

## Module 3: Gene-level dCCC phenotypes in LRI–LOGFC space

**Purpose.** To compare **all genes** in the panel within one cohort by asking: when gene \(g\) is perturbed, which ligand–receptor programs change in a similar way?

### 3.1 Core representation: LRI × gene LOGFC matrix

For each gene \(g\), we collapsed \(T^{(g)}_d\) to one value per LRI by averaging LOGFC over all malignant-involved DE rows sharing that LRI:

\[
M_{l,g} = \frac{1}{|S_{l,g}|} \sum_{r \in S_{l,g}} \text{LOGFC}_r
\]

where \(S_{l,g}\) is the set of rows for ligand–receptor \(l\) and gene \(g\), and \(l\) indexes the **LRI** (ligand–receptor pair, agnostic of emitter/receiver cell type at this step).

Stacking all genes yields matrix **\(M_d \in \mathbb{R}^{L \times G}\)** (rows = LRIs, columns = genes). Entries can be NA when an LRI was not differentially detected for that gene.

**Sparsity filter.** LRIs present in fewer than **2** genes were removed:

\[
\mathcal{L}' = \{ l : \sum_g \mathbf{1}(M_{l,g} \text{ not NA}) \geq 2 \}
\]

This retains comparable rows for similarity analysis while dropping extremely sparse LRIs.

**Interpretation.** Each column of \(M_d\) is the **within-cohort differential communication phenotype of gene \(g\)**: a vector of LRI-level log-fold changes induced by perturbing \(g\) in that tumor microenvironment.

### 3.2 Gene–gene similarity (cosine)

Pairwise **cosine similarity** between columns of \(M_d\) defines gene–gene concordance of dCCC phenotypes. For genes \(g_i, g_j\), using only LRIs where both values are observed:

\[
\text{sim}(g_i, g_j) = \frac{\sum_{l \in \mathcal{L}_{ij}} M_{l,g_i} M_{l,g_j}}{\|\mathbf{m}_{g_i}\|_2 \|\mathbf{m}_{g_j}\|_2}
\]

where \(\mathcal{L}_{ij}\) is the set of LRIs with non-missing entries for both genes. The diagonal (self-similarity) was set to NA for visualization.

**Heatmap.** We clustered genes hierarchically (complete linkage on distance \(1 - \text{sim}\)). The heatmap displays \(\text{sim}(g_i,g_j)\) on \([-1,1]\) (red–black–green). Off-diagonal missing similarities were imputed to 0 for display; **diagonal cells were set to NA** (grey) to avoid emphasizing trivial self-similarity. Column labels were omitted (symmetric matrix); row labels show gene symbols.

**Interpretation.** Clusters in the heatmap identify **groups of cancer genes that induce similar LRI-level dCCC programs** within the cohort.

### 3.3 Low-dimensional embedding and discrete clusters (UMAP)

To visualize gene relationships, we transposed \(M_d\) to obtain gene × LRI matrices. Missing values were imputed to **0** for embedding only (not for cosine similarity). We applied:

1. **Hierarchical clustering** on cosine distance (complete linkage), cutting the tree at **\(k = 5\)** clusters.
2. **UMAP** (cosine metric, `n_neighbors = min(15, G−1)`, `min_dist = 0.1`, seed = 42).

Genes were plotted in UMAP space, colored by hierarchical cluster, with gene symbols labeled. A CSV export records cluster membership and the number of non-NA LRIs per gene.

**Interpretation.** The UMAP provides a **2D summary of within-cohort dCCC phenotype similarity**, complementary to the heatmap: nearby genes perturb overlapping ligand–receptor communication programs.

---

## Module 4 (optional): GO biological process enrichment of dCCC contexts

> **Note:** This module is implemented but **disabled by default** in the active pipeline because Fisher testing over the full gene panel is computationally intensive. Enable by uncommenting `sections/04_gobp_fisher_volcano.Rmd` in `index.Rmd`.

**Purpose.** To test whether UP- or DOWN-regulated DE interactions for gene \(g\) are **enriched** in specific biological contexts defined by **(emitter, receiver, GO biological process)**.

**GO annotation.** LRIs were linked to Gene Ontology biological process (GO-BP) terms via `scDiffCom::LRI_human$LRI_curated_GO`, joined to GO term names.

**Contingency framework.** For each gene \(g\) and each context \((E, R, \text{GO})\) and regulation direction \(R \in \{\text{UP}, \text{DOWN}\}\):

|  | In GO context | Not in GO context |
|--|---------------|-------------------|
| DE, regulation \(R\) | \(a\) | \(b\) |
| Not DE / other | \(c\) | \(d\) |

- **Foreground:** rows in \(T^{(g)}_d\) with the given regulation.
- **Background:** `cci_table_raw` restricted to malignant-involved interactions, annotated with GO-BP.

Contexts required ≥ **5** LRIs in the GO term in the raw background and ≥ **2** DE interactions in the foreground. **Fisher’s exact test** was applied per context; **Benjamini–Hochberg** adjustment was performed across tests grouped by \((E, R, \text{GO}, \text{REGULATION})\).

**Volcano plots.** For selected genes, UP-regulated enrichments were visualized with \(x = \log_2(\text{odds ratio})\), \(y = -\log_{10}(p_{\text{BH}})\), points colored by emitter→receiver pair, and top GO terms labeled per gene.

**Interpretation.** This module links **dCCC phenotypes to interpretable biological processes and cell-type axes**, moving from LRI-level vectors to formal enrichment statistics.

---

## Representations summary

| Representation | Dimensions | What it encodes | Used in |
|----------------|------------|-----------------|---------|
| Filtered DE table \(T^{(g)}_d\) | interactions × metadata | Per-gene malignant dCCC at full CCI resolution (LRI + cell types) | Module 2, Fisher background/foreground |
| LRI × gene matrix \(M_d\) | LRIs × genes | Gene-specific dCCC phenotype (LRI-averaged LOGFC) | Module 3 |
| Gene × gene similarity | genes × genes | Concordance of dCCC phenotypes (cosine) | Heatmap |
| UMAP coordinates | genes × 2 | Nonlinear embedding of \(M_d^\top\) | UMAP plot |
| Hierarchical clusters | genes | Discrete phenotype groups (\(k=5\)) | UMAP color, CSV |
| Fisher results (optional) | tests × contexts | GO-BP enrichment per gene/context | Volcano |

---

## Parameters and software

| Parameter | Value |
|-----------|--------|
| Malignant cell type | Tumor |
| Min genes per LRI (matrix filter) | 2 |
| UMAP clusters (\(k\)) | 5 |
| UMAP `n_neighbors` | min(15, \(G-1\)) |
| UMAP `min_dist` | 0.1 |
| UMAP seed | 42 |
| Boxplot: ER pairs labeled | Top 8 by LogFC spread |
| Boxplot figure size | 18 × 10 in, 300 dpi |
| Fisher: min GO size / min DE (optional) | 5 / 2 |
| BH cutoff (volcano, optional) | 0.05 |

**Software:** R (tidyverse, scDiffCom, ggplot2, ggrepel, pheatmap, umap, gplots).

**Code:** `Thesis/CCC/scripts/scDiffCom/post_analysis/within_cohort/`

**Outputs:** `~/Thesis/CCC/outputs/scDiffCom/plots/within_cohort/{cohort}/`

---

## Suggested thesis subsection structure

1. **Within-cohort differential communication analysis** (overview)
2. **Gene-of-interest dissection** (AXL boxplot; cite Kurten/Bill as examples)
3. **Gene-level phenotype clustering in LRI space** (matrix, heatmap, UMAP)
4. **GO-BP enrichment** (optional paragraph if Fisher is enabled for the final thesis)

---

## Example Results sentence (template)

> In cohort Kurten_HNSC, perturbation of 124 cancer genes yielded 781 LRIs with detectable malignant dCCC across ≥2 genes. Cosine similarity and UMAP embedding identified five gene clusters with distinct LRI-level phenotypes. Dissection of AXL revealed [describe dominant ER pairs] with \(N = 905\) differential interactions across 27 emitter–receiver pairs.

Replace numbers with cohort-specific values from plot subtitles and `{cohort}_gene_umap_clusters.csv`.

---

## Pipeline flow

```
scDiffCom per-gene RDS (HIGH vs LOW)
        ↓
Module 1: Filter (Tumor, IS_CCI_DE, unknown cell types)
        ↓
    ┌───┴───┐
    ↓       ↓
Module 2   Module 3 (+ optional Module 4)
AXL boxplot   LRI×gene matrix → cosine heatmap + UMAP
(ER-resolved) Fisher GO-BP (optional)
```
