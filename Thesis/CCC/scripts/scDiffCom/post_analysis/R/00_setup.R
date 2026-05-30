# Utils.R is sourced from index.Rmd before this file runs.

library(tidyverse) # Loads dplyr, tidyr, ggplot2, purrr, stringr, etc.
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(scDiffCom)
library(ggplot2)
#library(ggrepel)
library(SuperRanker)


# BASE_RESULTS_DIR <- "~/Thesis/CCC/outputs/scDiffCom/scDiffComs/Consensus_Cell_Type"
BASE_RESULTS_DIR <- "~/CCC-scDiffCom/results/split-by-rank-genes-v2"
# MALIGNANT_CELLTYPE <- c("Epithelial", "Malignant", "Tumor", "Cancer")
MALIGNANT_CELLTYPE <- c("Tumor")
USE_COSINE         <- TRUE           # TRUE → cosine distance instead of correlation
MIN_GENES_PER_CCI  <- 2              # CCIs must appear in ≥ N genes (filters sparse rows)
# OUTPUT_DIR         <- "~/Thesis/CCC/outputs/scDiffCom/plots/v2"            # where to save output PNGs (v2 — original, no unknown-type filter)
OUTPUT_DIR         <- "~/Thesis/CCC/outputs/scDiffCom/plots/v6"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

LRI_GO <- scDiffCom::LRI_human$LRI_curated_GO
BP_GO_INFO <- scDiffCom::gene_ontology_level %>%
  filter(ASPECT == "biological_process") %>%
  select(ID, NAME)
LRI_GO_BP <- LRI_GO %>%
  inner_join(BP_GO_INFO, by = c("GO_ID" = "ID")) %>%
  mutate(GO_NAME = .data$NAME) %>%
  select(LRI, GO_ID, GO_NAME)
