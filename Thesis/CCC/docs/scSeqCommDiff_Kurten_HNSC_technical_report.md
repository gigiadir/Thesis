# scSeqCommDiff on Kurten_HNSC: Technical Report

**Audience:** Bioinformatics and data engineering  
**Pipeline script:** [`Thesis/CCC/scripts/run_scSeqCommDiff_Kurten_HNSC.R`](../scripts/run_scSeqCommDiff_Kurten_HNSC.R)  
**Default outputs:** `Thesis/CCC/outputs/scSeqCommDiff/Kurten_HNSC/`  
**Scope:** Proof-of-concept (single `Rscript` run; no batch/`qsub` wrapper)

**Primary references**

- Cesaro G, Baruzzo G, Tussardi G, Di Camillo B. *Differential cellular communication inference framework for large-scale single-cell RNA-sequencing data.* NAR Genomics and Bioinformatics. 2025;7(2):lqaf084. [doi:10.1093/nargab/lqaf084](https://doi.org/10.1093/nargab/lqaf084)
- Baruzzo G, Cesaro G, Di Camillo B. *Identify, quantify and characterize cellular communication from single-cell RNA sequencing data with scSeqComm.* Bioinformatics. 2022;38(7):1920–1929. [doi:10.1093/bioinformatics/btac036](https://doi.org/10.1093/bioinformatics/btac036)
- Package tutorial: [https://sysbiobig.gitlab.io/scSeqComm/articles/scseqcomm.html](https://sysbiobig.gitlab.io/scSeqComm/articles/scseqcomm.html)
- CClens: [https://gitlab.com/sysbiobig/cclens](https://gitlab.com/sysbiobig/cclens)

---

## 1. Core Methodology & Calculation

### 1.1 Framework overview

**scSeqCommDiff** extends the base **scSeqComm** framework to compare cell–cell communication between two experimental conditions at two coupled levels:

| Level | Biological question | Mechanistic layer |
|-------|---------------------|-------------------|
| **Intercellular** | Which ligand–receptor (L-R) pairs differ between conditions across sender→receiver cell-type pairs? | L-R expression in `cluster_L` (sender) and `cluster_R` (receiver) |
| **Intracellular** | In receptor-expressing clusters, which pathway-linked target genes differ downstream of receptor activation? | Receptor → TF → target via transcriptional regulatory network (TRN) and pathway priors |

The Kurten_HNSC POC uses `scSeqComm::scSeqComm_differential()` with `scenario = "multi-sample"`, which is the mode intended when **biological replicates (patients)** exist within each condition and intra-condition variance must be preserved.

### 1.2 Intercellular signaling model

For each L-R pair \((L, R)\) and ordered cell-type pair \((C_L, C_R)\):

1. **Per-sample intercellular scores** are computed using the selected scheme (default **`scSeqComm`**: combines expression and cell-type specificity, and supports multi-subunit complexes in the L-R database via comma-separated subunits).

2. In **multi-sample** mode, differential intercellular inference compares the **distribution of sample-level scores** between conditions (e.g. AXL `HIGH` vs `LOW` patients), rather than pooling all cells and treating them as independent replicates.

3. The pipeline sets `alternative_inter = "two.sided"` and `cond_names = c("HIGH", "LOW")`, producing:
   - `S_inter_HIGH`, `S_inter_LOW`: condition-specific aggregate intercellular evidence
   - `pvalue_adj_S_inter`: multiplicity-adjusted significance of the difference (default correction via `padj_method`, typically Bonferroni family)
   - `effsize_S_inter`: standardized magnitude of change (Cohen’s \(d\)-style effect size in multi-sample mode; contrasts with `logFC_S_inter` in multi-condition mode)

**Interpretation:** A row with `ligand = GAS6`, `receptor = AXL`, `cluster_L = Macrophage`, `cluster_R = Tumor` describes differential **Macrophage→Tumor** signaling through GAS6–AXL between AXL-high and AXL-low patient groups. Positive `effsize_S_inter` under `cond_names = c("HIGH", "LOW")` indicates stronger intercellular score in **HIGH** relative to **LOW**.

### 1.3 Intracellular signaling model

Intracellular analysis links receptors to downstream transcriptional consequences:

1. **Prior structure:** `TF_reg_DB` (default: merged TRRUST v2, HTRIdb, RegNetwork high-confidence edges) defines TF → target gene regulons. `R_TF_association` (default: `TF_PPR_KEGG_human`) provides receptor–TF links within KEGG pathways via Personalized PageRank (PPR) scores.

2. **Differential target identification:** For each receptor \(R\) in receiver cluster \(C_R\), target genes in the TRN are tested for differential expression between conditions. In **multi-sample** mode, `DEmethod` defaults to **`pseudo.wilcoxon`**: expression is aggregated to **cluster × sample** pseudobulks (`aggregation_method = "mean"` by default), then Wilcoxon tests are run between condition groups at the pseudobulk level. This explicitly respects **patient as replicate**, unlike cell-level Wilcoxon tests that inflate \(N\).

3. **Intracellular score `S_intra`:** Summarizes evidence that downstream pathway activity differs, conditioned on the receptor–pathway association. Pathway-specific rows include:
   - `pathway`: KEGG/Reactome pathway name
   - `genes`, `up_genes`, `down_genes`: DE targets linked to the receptor in that pathway
   - Directionality of DE follows `cond_names` order when `only.pos = TRUE`; the Kurten script uses `only.pos = FALSE`

### 1.4 Multi-sample vs multi-condition (why Kurten uses multi-sample)

| Aspect | `multi-condition` | `multi-sample` (Kurten POC) |
|--------|-------------------|-----------------------------|
| Unit of replication | Cells pooled within condition | **Patient** (`Sample_ID`) |
| Intercellular test | Permutation over pooled cells (`Nrep` default 1000) | **Sample-level** comparison of intercellular scores |
| Intracellular DE | Cell-level Wilcoxon (default) | **Pseudobulk Wilcoxon** per cluster×sample |
| Requires `Sample_ID` | Optional (ignored) | **Required** |

Kurten_HNSC has **5 HIGH and 5 LOW** patients after excluding `MID` and unannotated cells. Multi-sample mode avoids pseudoreplication from ~40k cells and aligns with the patient-stratified AXL grouping already used in the scDiffCom arm of this project.

### 1.5 Comparison to scDiffCom (project context)

Your existing [`scDiffComPipeline.R`](../scripts/scDiffComPipeline.R) performs permutation-based differential CCI with explicit `seurat_condition_id` and patient-aware structure. scSeqCommDiff differs in:

- Joint **inter + intra** modeling in one `differential_comm` table
- Built-in TRN/pathway intracellular layer (not only L-R logFC)
- Native **CClens** export for interactive exploration

Both tools should be interpreted as complementary: scDiffCom for curated L-R GO/CCI ranking; scSeqCommDiff for pathway-resolved downstream responses and sample-level intercellular statistics.

---

## 2. Data Flow & Output Architecture

### 2.1 Pipeline inputs (Kurten_HNSC)

| Input | Path / column | Role |
|-------|---------------|------|
| Seurat object | `scObjects/Kurten_HNSC.RData` | Normalized `RNA/data` matrix |
| AXL grouping | `CCC-PreProcess/results/Kurten_HNSC/AXL_Kurten_HNSC_grouped.rds` | Patient → `AXL_EXP` (`HIGH`/`LOW`/`MID`) |
| `Cell_ID` | `colnames(gene_expr)` | Cell barcodes |
| `Cluster_ID` | `Cell_Type` | Cell-type labels |
| `Condition_ID` | `AXL_EXP` (HIGH/LOW only in POC) | Experimental contrast |
| `Sample_ID` | `Patient` | Biological replicate |

**Pre-QC in script:** Drops unannotated cells, optional `Multi`, cell types with `< min_cells` (default 30) in either condition, and `Epithelial` (79 cells; fails min_cells in this dataset).

### 2.2 `scSeqComm_differential()` return object

R `list` of length 3:

```r
scSeqCommDiff_res$differential_comm          # primary integrative table
scSeqCommDiff_res$`intercellular signaling`  # per-condition intercellular results
scSeqCommDiff_res$`intracellular signaling`  # intracellular pathway-level results
```

Saved as `Kurten_HNSC_scSeqCommDiff_res.rds`.

### 2.3 `differential_comm` column glossary (multi-sample)

| Column | Type | Meaning |
|--------|------|---------|
| `ligand` | character | Ligand gene (or comma-separated subunits) |
| `receptor` | character | Receptor gene (or subunits) |
| `cluster_L` | character | Sender cell type (ligand-expressing) |
| `cluster_R` | character | Receiver cell type (receptor-expressing) |
| `S_inter_HIGH` | numeric | Intercellular score in HIGH patients |
| `S_inter_LOW` | numeric | Intercellular score in LOW patients |
| `pvalue_adj_S_inter` | numeric | Adjusted p-value for differential intercellular score |
| `effsize_S_inter` | numeric | Effect size (HIGH vs LOW) |
| `S_intra` | numeric | Intracellular differential evidence (pathway-specific rows) |
| `pathway` | character | KEGG/Reactome pathway name |
| `genes` | character | Comma-separated DE target genes |
| `up_genes` | character | Targets up in HIGH vs LOW (order-dependent) |
| `down_genes` | character | Targets down in HIGH vs LOW |

Additional columns may appear depending on package version and scoring options. The full unfiltered table is written to `Kurten_HNSC_differential_comm_full.tsv`.

### 2.4 Post-analysis outputs (Tumor → AXL hierarchy)

Post-analysis does **not** re-run inference; it subsets `differential_comm`:

```
differential_comm (full)
    └── differential_comm_Tumor_involved.tsv   # cluster_L == "Tumor" | cluster_R == "Tumor"
            └── AXL_Tumor_differential_comm_all.tsv
                    └── AXL_Tumor_differential_comm_significant.tsv  # p_adj & |effsize| filters
```

| File | Filter logic |
|------|----------------|
| `differential_comm_Tumor_involved.tsv` | Malignant-centric CCI only |
| `AXL_Tumor_differential_comm_all.tsv` | Tumor subset + `AXL` in ligand or receptor (incl. complexes) |
| `AXL_Tumor_differential_comm_significant.tsv` | Above + `pvalue_adj_S_inter < 0.05` + `\|effsize_S_inter\| ≥ 0.5` |
| `AXL_Tumor_summary.txt` | Human-readable counts and top hits |
| `Kurten_HNSC_CClens_Tumor_diff_input.tsv` | CClens-ready (default: all Tumor-involved rows) |

Auxiliary QC: `qc_celltype_condition_counts.tsv`, `patient_balance.tsv`, `run_manifest.json`, `sessionInfo.txt`.

### 2.5 CClens mapping

CClens accepts a flat table and supports **two-condition differential** uploads ([CClens README](https://gitlab.com/sysbiobig/cclens/-/raw/master/README.md)).

The script adds:

- `sender` ← `cluster_L`
- `receiver` ← `cluster_R`

Required conceptual fields for differential mode:

- `ligand`, `receptor`, `sender`, `receiver`
- Condition-specific intercellular scores: `S_inter_HIGH`, `S_inter_LOW`
- Optional: `pvalue_adj_S_inter`, `effsize_S_inter`, `S_intra`, gene list columns

**Usage:** After `devtools::install_gitlab("sysbiobig/cclens")`, run `cclens::run_cclens()`, upload `Kurten_HNSC_CClens_Tumor_diff_input.tsv`, confirm score column names in the UI, then use filters (ligand/receptor/cell-type sliders; p-value and effect-size thresholds) across tabs: *Tabular Visualization*, *Overall Communication*, *Cell Clusters*, *Ligands and Receptors*.

CLI: `--cclens_source tumor|axl|axl_significant` selects which filtered table is exported.

---

## 3. Application to the HNSC (Kurten) Dataset

### 3.1 Biological and technical context

- **Source object:** Kurten et al. HNSC scRNA-seq (`Kurten_HNSC`, 73,080 cells raw; ~40,799 after HIGH/LOW + annotation QC in a typical prep run).
- **Contrast:** Patient-level **AXL expression in malignant cells** (from prior `scDiffComPreprocess.R` quantile split), not native `HPV+`/`HPV-` metadata. Five patients per arm; three `MID` patients excluded; HN02–HN06 absent from AXL grouping (insufficient malignant cells in preprocess).
- **HPV confounding:** HIGH/LOW groups are not HPV-balanced (e.g. HPV+ enriched in some HIGH patients). Results describe **AXL-stratified** communication, not HPV status per se.
- **CD45 fractions:** Both `CD45p` and `CD45n` cells are retained; `Sample_ID = Patient` keeps paired fractions under the same replicate.

### 3.2 What the Tumor-centric filter yields

Only interactions where **`Tumor` is sender or receiver** are interpreted in downstream tables and this report’s application layer.

| Orientation | Example | Data points retained |
|-------------|---------|---------------------|
| **X → Tumor** | `Macrophage → Tumor`, GAS6–AXL | Ligand from immune/stroma; receptor on tumor; differential `S_inter_*`, `effsize_S_inter`, intracellular DE in `cluster_R = Tumor` |
| **Tumor → X** | `Tumor → Fibroblast`, AXL–? | Tumor-secreted ligand; receiver stroma/immune response |
| **Excluded** | `CD8T → Macrophage` | Dropped from `differential_comm_Tumor_involved.tsv` |

This mirrors the malignant-focused filter in [`scDiffCom-PostAnalysis.Rmd`](../scripts/scDiffCom-PostAnalysis.Rmd) (`EMITTER_CELLTYPE` / `RECEIVER_CELLTYPE` ∈ Tumor).

### 3.3 What the AXL filter yields (within Tumor subset)

Applied **after** the Tumor filter. A row is kept if `AXL` appears in `ligand` or `receptor` (including multi-subunit strings like `GAS6,AXL`).

**Intercellular data points isolated:**

- L-R pairs where AXL is the **receptor** (e.g. GAS6–AXL, PROS1–AXL, TGFB1–AXL) with sender/receiver cell types and `S_inter_HIGH` vs `S_inter_LOW`
- Pairs where AXL is the **ligand** (less common in curated DBs) describing tumor-autocrine or tumor→stroma signaling
- Statistical calls: `pvalue_adj_S_inter`, `effsize_S_inter` for condition contrast

**Intracellular data points isolated:**

- Rows with `receptor == AXL` and `cluster_R == Tumor`: pathway names, `S_intra`, and `up_genes`/`down_genes` reflecting **transcriptional response in tumor cells** downstream of AXL between HIGH and LOW groups
- Pathway examples often include RTK, integrin, or TAM-relevant KEGG terms when GAS6–AXL axis is active

**Significant tier (`AXL_Tumor_differential_comm_significant.tsv`):** Subset with `pvalue_adj_S_inter < 0.05` and `|effsize_S_inter| ≥ 0.5` (tunable via `--p_adj_max`, `--effsize_min`). These are the primary candidates for follow-up and CClens `--cclens_source axl_significant`.

### 3.4 Running the POC

```bash
# Install once
R -e "devtools::install_gitlab('sysbiobig/scseqcomm')"

# Full run (memory-intensive; use subset_genes or bigmatrix if needed)
Rscript Thesis/CCC/scripts/run_scSeqCommDiff_Kurten_HNSC.R --n_cores 4

# Re-post-process existing RDS
Rscript Thesis/CCC/scripts/run_scSeqCommDiff_Kurten_HNSC.R --skip_run
```

**Memory fallbacks:** `--subset_genes` (L-R + TRN genes only), `--bigmatrix`, or interactive node with ≥64 GB RAM recommended for full 45k×40k matrix.

### 3.5 Expected validation checks

1. `patient_balance.tsv`: 5 HIGH, 5 LOW patients.
2. Every row in `differential_comm_Tumor_involved.tsv` has `Tumor` in `cluster_L` or `cluster_R`.
3. `AXL_Tumor_*` files ⊆ Tumor-involved set.
4. Spot-check known TAM axis pair (GAS6–AXL, `Macrophage → Tumor`) in significant export if biology supports it.

---

## Appendix: Script parameter reference

| CLI flag | Default | Notes |
|----------|---------|-------|
| `--seurat_path` | `scObjects/Kurten_HNSC.RData` | |
| `--axl_patient_summary` | `CCC-PreProcess/.../AXL_Kurten_HNSC_grouped.rds` | |
| `--cond_high` / `--cond_low` | HIGH / LOW | `cond_names` order |
| `--min_cells` | 30 | Per cell type per condition |
| `--n_cores` | 4 | Parallel backend `doParallel` |
| `--subset_genes` | off | POC memory relief |
| `--focus_gene` | AXL | Post-analysis highlight |
| `--tumor_label` | Tumor | Post-analysis malignant filter |
| `--cclens_source` | tumor | CClens export scope |

---

*Document version: aligned with `run_scSeqCommDiff_Kurten_HNSC.R` POC implementation.*
