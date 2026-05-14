#!/usr/bin/env Rscript
# completeness_audit.R
#
# Compares every Gene × Dataset pair from a master gene list against the
# files that actually exist under CCC-PreProcess/results, then reports
# and saves the gaps.
#
# Usage (from shell):
#   Rscript ~/Scripts/completeness_audit.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# ─── Configuration ─────────────────────────────────────────────────────────────
GENE_LIST_PATH <- "/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/data/Complexes.Oncogenes.OncoKB.Cosmic.NCG.rds"
RESULTS_ROOT   <- path.expand("~/CCC-PreProcess/results")
OUTPUT_DIR     <- path.expand("~/CCC-PreProcess/audit")

TARGET_DATASETS <- c(
  "Bassez2021_Breast_step4",
  "Qian2020_Breast_step4",
  "Wu2021_Breast_step4"
)

# ─── 1. Load master gene list ──────────────────────────────────────────────────
message("── Step 1: Loading master gene list")

stopifnot(file.exists(GENE_LIST_PATH))
raw <- readRDS(GENE_LIST_PATH)

master_genes <- if (is.data.frame(raw)) {
  # Accept any common gene-column naming convention
  gene_col <- intersect(
    c("gene", "Gene", "GENE", "gene_name", "Gene_name", "symbol", "Symbol", "hgnc_symbol"),
    colnames(raw)
  )[1]
  if (is.na(gene_col)) {
    stop(
      "Cannot identify gene column. Available columns: ",
      paste(colnames(raw), collapse = ", ")
    )
  }
  message("  Using column: '", gene_col, "'")
  unique(raw[[gene_col]])
} else {
  unique(as.character(raw))
}

master_genes <- master_genes[!is.na(master_genes) & nzchar(master_genes)]
message(sprintf("  %d unique genes loaded", length(master_genes)))

# ─── 2. File Discovery ─────────────────────────────────────────────────────────
message("── Step 2: Discovering files")

if (!dir.exists(RESULTS_ROOT)) {
  stop("Results directory not found: ", RESULTS_ROOT)
}

all_files <- list.files(
  path      = RESULTS_ROOT,
  pattern   = "\\.rds$",
  recursive = TRUE,
  full.names = FALSE
)
message(sprintf("  %d .rds files found under %s", length(all_files), RESULTS_ROOT))

# ─── 3. Parse Gene & Dataset from filenames ────────────────────────────────────
# File pattern: {GENE}_{DATASET}_grouped.rds
#
# Gene names can contain underscores, so we anchor on the *known* dataset
# strings rather than splitting naively on '_'.
#
# Regex strategy:
#   ^(.+?)_(<DATASET1>|<DATASET2>|<DATASET3>)_grouped\.rds$
#
# The non-greedy (.+?) expands character-by-character until the regex engine
# finds a position where the literal dataset name follows. Because the dataset
# strings are long and specific, this unambiguously separates gene from dataset
# even when the gene name itself contains underscores.
message("── Step 3: Parsing filenames")

dataset_alt <- paste(TARGET_DATASETS, collapse = "|")
file_regex  <- sprintf("^(.+?)_(%s)_grouped\\.rds$", dataset_alt)

existing <- tibble(path = all_files) %>%
  mutate(
    filename = basename(path),
    gene     = str_match(filename, file_regex)[, 2],
    dataset  = str_match(filename, file_regex)[, 3]
  ) %>%
  filter(!is.na(gene), !is.na(dataset)) %>%
  distinct(gene, dataset)

n_unmatched <- length(all_files) - nrow(existing)
if (n_unmatched > 0) {
  message(sprintf(
    "  WARNING: %d file(s) did not match the expected pattern and were skipped",
    n_unmatched
  ))
}
message(sprintf("  %d valid gene×dataset pairs parsed from disk", nrow(existing)))

# ─── 4. Build Theoretical Complete Set ────────────────────────────────────────
message("── Step 4: Building theoretical complete set")

theoretical <- crossing(
  gene    = master_genes,
  dataset = TARGET_DATASETS
)
message(sprintf(
  "  %d total pairs expected  (%d genes × %d datasets)",
  nrow(theoretical), length(master_genes), length(TARGET_DATASETS)
))

# ─── 5. Gap Analysis ──────────────────────────────────────────────────────────
message("── Step 5: Computing gaps")

missing_pairs <- theoretical %>%
  anti_join(existing, by = c("gene", "dataset")) %>%
  arrange(dataset, gene)

extra_pairs <- existing %>%
  anti_join(theoretical, by = c("gene", "dataset")) %>%
  arrange(dataset, gene)

message(sprintf("  Missing pairs  (in master, not on disk): %d", nrow(missing_pairs)))
message(sprintf("  Extra pairs    (on disk, not in master): %d", nrow(extra_pairs)))

# ─── 6. Summary Table ─────────────────────────────────────────────────────────
message("── Step 6: Summarising")

summary_tbl <- missing_pairs %>%
  count(dataset, name = "n_missing") %>%
  # Ensure every dataset appears even if it has zero missing files
  right_join(tibble(dataset = TARGET_DATASETS), by = "dataset") %>%
  mutate(
    n_missing    = replace_na(n_missing, 0L),
    n_total      = length(master_genes),
    n_present    = n_total - n_missing,
    pct_complete = round(100 * n_present / n_total, 2)
  ) %>%
  arrange(dataset)

cat("\n")
cat("══ Completeness Summary ════════════════════════════════════════════════════\n")
print(summary_tbl, n = Inf)
cat("════════════════════════════════════════════════════════════════════════════\n\n")

# ─── 7. Save Outputs ──────────────────────────────────────────────────────────
message("── Step 7: Saving outputs")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

missing_csv  <- file.path(OUTPUT_DIR, "missing_pairs.csv")
summary_csv  <- file.path(OUTPUT_DIR, "completeness_summary.csv")
missing_rds  <- file.path(OUTPUT_DIR, "missing_pairs.rds")

write.csv(missing_pairs, missing_csv,  row.names = FALSE)
write.csv(summary_tbl,   summary_csv, row.names = FALSE)
saveRDS(missing_pairs,   missing_rds)

if (nrow(extra_pairs) > 0) {
  write.csv(extra_pairs, file.path(OUTPUT_DIR, "extra_pairs.csv"), row.names = FALSE)
  message("  Extra pairs also saved (files on disk not in master gene list)")
}

message(sprintf("  missing_pairs.csv  → %s", missing_csv))
message(sprintf("  missing_pairs.rds  → %s", missing_rds))
message(sprintf("  completeness_summary.csv → %s", summary_csv))
message("Done.")
