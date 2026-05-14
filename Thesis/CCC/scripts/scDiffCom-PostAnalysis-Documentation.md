# scDiffCom-PostAnalysis — Code Documentation

> This document describes the logic, inputs, and outputs for each major section of `scDiffCom-PostAnalysis.Rmd` (lines 1–536). The underlying R code is not modified.

---

## Table of Contents

1. [Analysis Globals & Helpers](#1-analysis-globals--helpers)
2. [Analysis Utils](#2-analysis-utils)
3. [Load Data](#3-load-data)
4. [Shared Config, Helpers & CCI Set Construction](#4-shared-config-helpers--cci-set-construction)
5. [CTLA4 Case Study — 6×6 Jaccard Heatmap](#5-ctla4-case-study--66-jaccard-heatmap-across-all-datasets)
6. [Global Jaccard UMAP — All Genes × All Six Datasets](#6-global-jaccard-umap--all-genes--all-six-datasets)
7. [Gene–Gene Consensus Jaccard Clustering](#7-genegene-consensus-jaccard-clustering)

---

## 1. Analysis Globals & Helpers

### Logic & Objectives

This setup chunk establishes the global runtime environment for all downstream analyses. It loads required R packages, defines project-wide constants, and pre-computes a curated lookup table that links ligand–receptor interactions (LRIs) to biological-process Gene Ontology (GO) terms.

The constants (e.g., `TOP_N_CCI`, `MALIGNANT_CELLTYPE`, `USE_COSINE`) act as single-point-of-control switches — changing a value here propagates through every subsequent chunk without manual edits elsewhere.

The **LRI–GO** table is built by joining `scDiffCom`'s curated LRI catalogue with GO biological-process annotations, enabling optional GO-term enrichment of any CCI set identified later.

### Inputs

- No external files. All data originate from built-in `scDiffCom` package objects:
  - `scDiffCom::LRI_human$LRI_curated_GO`
  - `scDiffCom::gene_ontology_level`

### Outputs

| Variable | Description |
|---|---|
| `BASE_RESULTS_DIR` | Root directory containing per-dataset scDiffCom `.rds` result files |
| `MALIGNANT_CELLTYPE` | Cell-type label treated as the tumour compartment (`"Epithelial"`) |
| `USE_COSINE` | Flag: use cosine distance (vs. correlation) for vector-space comparisons |
| `MIN_GENES_PER_CCI` | Minimum number of genes a CCI must appear in to pass sparse-filter steps |
| `OUTPUT_DIR` | Directory where output PNG figures are saved |
| `LRI_GO_BP` | Data frame linking each LRI to its biological-process GO ID and term name |

---

## 2. Analysis Utils

### Logic & Objectives

This chunk defines two reusable filtering helper functions that are called repeatedly during data loading. Centralising the filters here prevents logic duplication and ensures consistent filtering criteria across all six datasets.

- **`load.dataset.scDiffComs`** — discovers all `.rds` files for a named dataset directory, loads each one, and returns a named list keyed by gene name. This abstraction makes it easy to add new datasets without changing downstream code.
- **`filter.scDiffCom.cci_table_detected.for.malignant`** — retains only CCIs where at least one endpoint (emitter or receiver) is the designated malignant cell type, the interaction is statistically differentially expressed (`IS_CCI_DE == TRUE`), and the log fold-change is a finite number. Filtering to finite `LOGFC` removes interactions that were detected but could not be quantified (e.g., division-by-zero edge cases in scDiffCom).
- **`filter.scDiffCom.cci_table_detected.for.celltypes`** — a more specific variant that filters to a precise emitter/receiver pair; reserved for targeted cell-type pairwise comparisons.

### Inputs

- `dataset_name`: string matching a subdirectory under `BASE_RESULTS_DIR`
- `scDiffCom_obj`: an scDiffCom S4 object (slot `@cci_table_detected` accessed)
- `malignant_celltype`: defaults to the global `MALIGNANT_CELLTYPE` (`"Epithelial"`)

### Outputs

| Object | Description |
|---|---|
| `load.dataset.scDiffComs()` | Named list: `gene → scDiffCom S4 object` |
| `filter...for.malignant()` | Filtered data frame of DE CCIs involving the malignant cell type |
| `filter...for.celltypes()` | Filtered data frame of DE CCIs for a specific emitter–receiver pair |

---

## 3. Load Data

### Logic & Objectives

Six published single-cell RNA-seq tumour datasets are loaded here — three Breast cancer and three Lung cancer cohorts. Using multiple independent datasets from each cancer type is a core design principle: any CCI pattern that appears consistently across datasets is more likely to be a genuine biological signal rather than a dataset-specific artefact.

All six datasets were processed through scDiffCom pipeline step 4 (`_step4`), meaning conditions are already correctly oriented (treatment vs. control direction) and `do_switch = FALSE`.

After loading, each dataset's gene list is immediately filtered to malignant-involving interactions using `filter.scDiffCom.cci_table_detected.for.malignant`, producing the `.malignant` lists used in all downstream analyses.

### Inputs

- `.rds` files on disk under `BASE_RESULTS_DIR`, one per gene per dataset.
  - **Breast datasets**: `Bassez2021_Breast_step4`, `Qian2020_Breast_step4`, `Wu2021_Breast_step4`
  - **Lung datasets**: `Chan2021_Lung_step4`, `Laughney2020_Lung_step4`, `Xing2021_Lung_step4`

### Outputs

| Variable | Description |
|---|---|
| `Bassez2021_Breast_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Bassez Breast) |
| `Qian2020_Breast_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Qian Breast) |
| `Wu2021_Breast_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Wu Breast) |
| `Chan2021_Lung_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Chan Lung) |
| `Laughney2020_Lung_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Laughney Lung) |
| `Xing2021_Lung_step4.scDiffComs` | Named list: gene → raw scDiffCom S4 object (Xing Lung) |
| `*.malignant` (6 lists) | Filtered data frames: only DE CCIs involving the Epithelial/malignant compartment |

---

## 4. Shared Config, Helpers & CCI Set Construction

### Logic & Objectives

This chunk is the analytical workhorse that converts raw scDiffCom results into the core data structure used by every subsequent visualisation: a **set of top-N CCIs per gene per dataset**.

**Why top 500 CCIs?**
Each gene's scDiffCom result may contain hundreds of cell–cell interactions (CCIs) that are differentially expressed. To compare interaction profiles across datasets of differing sizes and cell-type compositions, we truncate to the top 500 interactions ranked by absolute log fold-change (`|LOGFC|`). This has two benefits:
1. It focuses on the most strongly differential interactions, reducing noise from borderline significant CCIs.
2. It creates a fair, fixed-length representation for each gene × dataset unit, making cross-dataset Jaccard comparisons interpretable.

**Jaccard distance** (`jaccard_dist`) measures the dissimilarity between two CCI sets as:

$$d_J(A, B) = 1 - \frac{|A \cap B|}{|A \cup B|}$$

A value of **0** means the two sets are identical (complete conservation of communication patterns); a value of **1** means no overlap at all (entirely divergent communication).

**`build_gene_dataset_jaccard`** assembles an N×N distance matrix for one gene across all N datasets, enabling direct comparison of how conserved that gene's tumour microenvironment communication profile is across cohorts.

**`gene_cci_top500_summary`** is a tidy inspection table (Gene | n_CCIs | CCI, one row per gene–CCI pair per dataset) that allows manual browsing of which specific interactions drive each gene's profile.

### Inputs

- `*.malignant` lists (6 lists, produced in the Load Data chunk)
- Global constants: `TOP_N_CCI` (500), `GENE_OF_INT` (`"CTLA4"`), `GENES_OF_INT`, `DATASET_LABELS`, `SEED`

### Outputs

| Variable | Description |
|---|---|
| `*_top_cci_sets` (6 named lists) | Per-dataset: gene → character vector of top-500 CCI IDs by `\|LOGFC\|` |
| `all_top_cci_sets` | Master named list combining all six `*_top_cci_sets` lists; keyed by full dataset name |
| `gene_cci_top500_summary` | Named list of tidy data frames (one per dataset): Gene, n_CCIs, CCI columns |
| `ds_short_map` | Named character vector mapping full dataset names to short display labels |
| `ds_cancer_map` | Named character vector mapping full dataset names to cancer type (`"Breast"` / `"Lung"`) |
| `source_colors` | Named colour palette for per-dataset UMAP aesthetics |
| `cancer_shapes` | Named shape palette: circle (Breast) vs. triangle (Lung) |

---

## 5. CTLA4 Case Study — 6×6 Jaccard Heatmap Across All Datasets

### Logic & Objectives

CTLA4 is selected as a representative immune-checkpoint gene to illustrate how CCI conservation can be assessed across all six tumour cohorts simultaneously. The 6×6 Jaccard distance matrix captures, for each pair of datasets, how much of CTLA4's top-500 DE interaction repertoire is shared.

This cross-tissue comparison (Breast vs. Lung) tests whether CTLA4-mediated tumour microenvironment communication is cancer-type-specific or universal. Within-cancer-type pairs (e.g., Bassez vs. Qian, both Breast) are expected to show lower distances (greater overlap) than cross-cancer-type pairs (Breast vs. Lung), assuming cancer-type-specific biology dominates. Deviations from this expectation suggest conserved pan-cancer immune mechanisms.

### Inputs

- `all_top_cci_sets`: master CCI set list (from chunk 4)
- `GENE_OF_INT`: `"CTLA4"`
- `TOP_N_CCI`: 500
- `DATASET_LABELS`: the six full dataset name strings

### Outputs

| Variable | Description |
|---|---|
| `ctla4_jaccard_mat` | 6×6 numeric matrix of pairwise Jaccard distances for CTLA4 |
| `matrix_ctla4_jaccard` | Alias of `ctla4_jaccard_mat`, exposed for interactive inspection |
| `ctla4_display_mat` | Display copy with `NA` cells replaced by 1 (for heatmap rendering) |
| `cell_annot` | Character matrix of formatted distance values for heatmap cell labels |
| Inline pheatmap figure | Rendered heatmap in the notebook output |

### How to Interpret — Jaccard Heatmap

The heatmap displays the **Jaccard distance** between every pair of datasets for the CTLA4 gene, computed on the top 500 CCIs (ranked by `|LOGFC|`):

$$d_J = 1 - \frac{|\text{CCI set}_i \cap \text{CCI set}_j|}{|\text{CCI set}_i \cup \text{CCI set}_j|}$$

- **Cell value = 0** (deep blue): the two datasets share an identical set of top CCIs for CTLA4 — complete conservation of the communication pattern.
- **Cell value = 1** (deep red): the two datasets share no CCIs at all — entirely divergent CTLA4-mediated communication.
- **Intermediate values**: partial overlap; e.g., 0.4 means roughly 60% of the union of CCIs is shared.
- **Diagonal**: always 0 (a dataset is identical to itself).
- **`N/A` cells** (displayed as 1 / coloured red): CTLA4 was absent (no detected DE CCIs) in one of the two datasets being compared.
- **Annotation bars** (top and left): colour-code each dataset by cancer type — red = Breast, blue = Lung. Clustering of same-colour rows/columns indicates cancer-type-specific conservation.

---

## 6. Global Jaccard UMAP — All Genes × All Six Datasets

### Logic & Objectives

To understand the **global landscape** of gene-level communication profiles across all datasets simultaneously, we embed every (gene, dataset) pair into 2D using UMAP. This reveals whether genes cluster by their CCI similarity profile, and whether those clusters align with cancer type, dataset of origin, or are gene-specific patterns.

The pipeline has five stages:
1. **Vocabulary construction**: collect every unique CCI that appears in any top-500 set across all six datasets — the global CCI vocabulary.
2. **Binary vectorisation**: for each (gene, dataset) pair, construct a binary presence/absence vector over this vocabulary. A `1` at position *k* means CCI *k* was in that gene's top-500 set for that dataset.
3. **Pairwise Jaccard distances**: compute the full distance matrix between all (gene, dataset) vectors using efficient matrix multiplication (`tcrossprod`).
4. **UMAP embedding**: reduce the high-dimensional Jaccard distance matrix to 2D with `uwot::umap`, using `n_neighbors = 15`, `min_dist = 0.1`.
5. **Visualisation**: scatter plot coloured by dataset source and shaped by cancer type, with gene-name labels.

### Inputs

- `all_top_cci_sets`: master CCI set list (from chunk 4)
- `ds_short_map`, `ds_cancer_map`, `source_colors`, `cancer_shapes`: metadata maps (from chunk 4)
- `GENES_OF_INT`: genes to highlight (`"CTLA4"`)
- `SEED`: 42 (for UMAP reproducibility)
- `K_UMAP`: 15 (UMAP number of neighbours), `MIN_DIST_UMAP`: 0.1

### Outputs

| Variable | Description |
|---|---|
| `bin_mat` | Binary matrix: rows = Gene_Dataset pairs, columns = global CCI vocabulary |
| `jac_dist_mat` | Full pairwise Jaccard distance matrix between all Gene_Dataset pairs |
| `umap_global` | Raw UMAP coordinate matrix (2 columns) |
| `umap_global_df` | Data frame: UMAP1, UMAP2, Label, Gene, Source, Cancer_Type, highlight flag |
| `matrix_umap_input` | Alias of `bin_mat`, exposed for inspection |
| `global_jaccard_umap.png` | Saved figure in `OUTPUT_DIR` |

### How to Interpret — Global Jaccard UMAP

**What each point represents:**
Each point in the UMAP is a **Gene–Dataset pair** (e.g., `CTLA4_Bassez2021`). There is one point per gene per dataset in which that gene had at least one DE CCI involving the malignant compartment.

**The underlying vector:**
Each point is embedded based on a **binary vector** over the global CCI vocabulary. Entry *k* in this vector is `1` if CCI *k* appeared among the top 500 `|LOGFC|`-ranked CCIs for that gene in that dataset, and `0` otherwise. Two points are close in the UMAP if their binary vectors have high Jaccard similarity (i.e., they share many of the same top-ranked CCIs).

**Aesthetic mappings:**

| Aesthetic | Variable | Meaning |
|---|---|---|
| **Color** | Dataset Source | Which of the six cohorts the point comes from (e.g., red = Bassez2021, orange = Qian2020). Points of the same colour come from the same study. |
| **Shape** | Cancer Type | Circle (●) = Breast cancer dataset; Triangle (▲) = Lung cancer dataset. |
| **Label** | Gene name | The gene whose CCI profile is being represented. |

**How to read clusters:**
- Points that cluster together share similar CCI interaction profiles across a given gene.
- If points from the same **cancer type** (same shape) cluster together regardless of dataset, it suggests cancer-type-specific communication programmes.
- If points from the same **dataset** (same colour) cluster together regardless of gene, it suggests technical/cohort effects dominate over biology.
- Genes labelled near each other have similar TME communication fingerprints.

---

## 7. Gene–Gene Consensus Jaccard Clustering

### Logic & Objectives

Where the UMAP explores individual gene × dataset profiles, this chunk asks a complementary question: **which genes share similar overall communication signatures when considered across all six datasets together?**

To answer this, each gene's CCI repertoire is aggregated into a **union set** (all CCIs that appear in at least one dataset's top-500 list for that gene). This consensus representation collapses dataset-level noise and focuses on the broadest communication fingerprint of each gene. The top 50 genes by total union CCI count are then selected to focus on the most communication-active genes.

A **Gene × Gene Jaccard similarity matrix** is computed on these union sets, and **Ward-D2 hierarchical clustering** is applied to identify modules of co-regulated communicating genes — genes in the same module tend to recruit the same CCIs across the tumour microenvironment.

**Why Ward-D2?** Ward's method minimises within-cluster variance at each merge step, producing compact, well-separated modules that are easier to interpret biologically.

**Why the top 50 genes by union CCI count?** Genes with very few CCIs across all datasets contribute little signal and inflate the distance matrix with near-1 values. Restricting to the top 50 most active genes gives a heatmap that is both readable and biologically informative.

### Inputs

- `all_top_cci_sets`: master CCI set list (from chunk 4)
- `GENES_OF_INT`: genes to highlight in the annotation bar
- `TOP_N_GENES`: 50 (number of genes to retain)
- `N_MODULES`: 5 (number of hierarchical clusters to cut)

### Outputs

| Variable | Description |
|---|---|
| `gene_cci_union` | Named list: gene → union of all CCI IDs across all six datasets |
| `top_genes` | Character vector of the top-50 genes by union CCI count |
| `gene_bin_mat` | Binary matrix: rows = top-50 genes, columns = union CCI vocabulary |
| `gene_jac_sim` | 50×50 Jaccard **similarity** matrix (Gene × Gene) |
| `gene_jac_dist` | 50×50 Jaccard **distance** matrix (= 1 − similarity) |
| `matrix_global_gene_similarity` | Alias of `gene_jac_sim`, exposed for inspection |
| `hc_genes` | `hclust` object (Ward-D2) used to order rows/columns in the heatmap |
| `gene_modules` | Named integer vector: gene → module index (1 to `N_MODULES`) |
| `module_summary` | Data frame: Gene, Module label (M1–M5), n_CCIs, sorted by module and CCI count |
| `gene_gene_jaccard_heatmap.png` | High-resolution (5000×5000 px, 300 dpi) heatmap saved to `OUTPUT_DIR` |

### How to Interpret — Gene–Gene Jaccard Similarity Heatmap

**Matrix values:**
Each cell (i, j) shows the **Jaccard similarity** between genes *i* and *j* based on their union CCI sets:

$$\text{sim}_J(i, j) = \frac{|\text{CCI union}_i \cap \text{CCI union}_j|}{|\text{CCI union}_i \cup \text{CCI union}_j|}$$

- **Value = 1** (dark blue): genes *i* and *j* recruit an identical set of CCIs across all datasets — they may operate in the same signalling pathway or share ligand/receptor partners.
- **Value = 0** (white): no CCI is shared between the two genes — entirely distinct communication programmes.

**Row/column ordering:**
Genes are ordered by Ward-D2 hierarchical clustering on the distance matrix (`1 − similarity`), so visually proximate rows/columns are most similar. Colour blocks along the diagonal correspond to the **M1–M5 modules** annotated in the sidebar.

**Annotation bars:**

| Bar | Meaning |
|---|---|
| **Module (M1–M5)** | Hierarchical cluster membership — genes within a module share overlapping CCI profiles and may represent co-regulated communication programmes |
| **Highlight (Yes/No)** | Gold = gene is in `GENES_OF_INT` (e.g., CTLA4); grey = all other genes |

**Module summary table (`module_summary`):**
Lists every top-50 gene with its module assignment and total union CCI count, sorted to make the most communication-active genes within each module immediately visible.
