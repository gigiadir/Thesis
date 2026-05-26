#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# scSeqCommDiff POC: Kurten_HNSC (AXL HIGH vs LOW, multi-sample)
#
# Prerequisites (install once; use user library — group R_LIBS is read-only):
#   export R_LIBS_USER="$HOME/R/x86_64-redhat-linux-gnu-library/4.5"
#   Rscript -e 'install.packages("devtools", lib=Sys.getenv("R_LIBS_USER"))'
#   Rscript -e 'devtools::install_gitlab("sysbiobig/scseqcomm")'
#   (If already installed: "Skipping install ... SHA1 has not changed" is OK.)
#   (Do not pass lib= to install_gitlab — devtools 2.4+ ignores/rejects it.)
#
# CClens (optional, separate command):
#   Rscript -e 'devtools::install_gitlab("sysbiobig/cclens")'
#
# POC usage (interactive or single run; no batch wrapper):
#   Rscript run_scSeqCommDiff_Kurten_HNSC.R --n_cores 4
#   Rscript run_scSeqCommDiff_Kurten_HNSC.R --skip_run  # prep + post only
#
# Reference: Cesaro et al. NAR Genomics and Bioinformatics 2025, lqaf084
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(jsonlite)
})

# --- CLI -----------------------------------------------------------------------
option_list <- list(
  make_option("--seurat_path", type = "character",
              default = "/gpfs0/bgu-ofircohen/users/gigiadir/scObjects/Kurten_HNSC.RData",
              help = "Path to Kurten_HNSC Seurat .RData"),
  make_option("--axl_patient_summary", type = "character",
              default = "/gpfs0/bgu-ofircohen/users/gigiadir/CCC-PreProcess/results/Kurten_HNSC/AXL_Kurten_HNSC_grouped.rds",
              help = "Patient-level AXL HIGH/LOW/MID assignments"),
  make_option("--output_dir", type = "character",
              default = "/gpfs0/bgu-ofircohen/users/gigiadir/Thesis/CCC/outputs/scSeqCommDiff/Kurten_HNSC",
              help = "Output directory"),
  make_option("--patient_col", type = "character", default = "Patient",
              help = "Seurat metadata column for sample/patient ID"),
  make_option("--cell_type_col", type = "character", default = "Cell_Type",
              help = "Seurat metadata column for cell type (Cluster_ID)"),
  make_option("--condition_col", type = "character", default = "AXL_EXP",
              help = "Metadata column for experimental condition (HIGH/LOW)"),
  make_option("--cond_high", type = "character", default = "HIGH",
              help = "Label for condition 1 (numerator in effect-size interpretation)"),
  make_option("--cond_low", type = "character", default = "LOW",
              help = "Label for condition 2"),
  make_option("--lr_db", type = "character", default = "ConnectomeDB_2020",
              help = "LR database: ConnectomeDB_2020 or Kumar_2018"),
  make_option("--min_cells", type = "integer", default = 30L,
              help = "Minimum cells per cell type per condition (scSeqComm)"),
  make_option("--n_cores", type = "integer", default = 4L,
              help = "Parallel cores for scSeqComm_differential"),
  make_option("--count_thr", type = "double", default = 1,
              help = "Expression threshold for DE gene testing (log-normalized data)"),
  make_option("--drop_multi", action = "store_true", default = FALSE,
              help = "Exclude Cell_Type == 'Multi'"),
  make_option("--subset_genes", action = "store_true", default = FALSE,
              help = "Restrict gene_expr to LR + TRN genes (memory POC fallback)"),
  make_option("--bigmatrix", action = "store_true", default = FALSE,
              help = "Pass bigmatrix=TRUE to scSeqComm_differential"),
  make_option("--p_adj_max", type = "double", default = 0.05,
              help = "Adjusted p-value cutoff for 'significant' AXL export"),
  make_option("--effsize_min", type = "double", default = 0.5,
              help = "Minimum |effsize_S_inter| for 'significant' AXL export"),
  make_option("--focus_gene", type = "character", default = "AXL",
              help = "Gene symbol for downstream highlight filter"),
  make_option("--tumor_label", type = "character", default = "Tumor",
              help = "Cell type label for malignant cluster in post-analysis"),
  make_option("--skip_run", action = "store_true", default = FALSE,
              help = "Skip scSeqComm_differential; load existing RDS from output_dir"),
  make_option("--cclens_source", type = "character", default = "tumor",
              help = "CClens export table: tumor | axl | axl_significant")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Helpers -------------------------------------------------------------------
stopifnot_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

gene_in_lr_field <- function(field, gene) {
  if (is.na(field) || !nzchar(field)) return(FALSE)
  parts <- unlist(strsplit(field, ",", fixed = TRUE))
  trimws(parts) %in% gene
}

row_involves_gene <- function(df, gene, ligand_col = "ligand", receptor_col = "receptor") {
  lig <- df[[ligand_col]]
  rec <- df[[receptor_col]]
  vapply(seq_len(nrow(df)), function(i) {
    gene_in_lr_field(lig[i], gene) || gene_in_lr_field(rec[i], gene)
  }, logical(1))
}

filter_tumor_involved <- function(df, tumor_label) {
  df[df$cluster_L == tumor_label | df$cluster_R == tumor_label, , drop = FALSE]
}

to_cclens_table <- function(df) {
  out <- df
  if ("cluster_L" %in% names(out)) out$sender <- out$cluster_L
  if ("cluster_R" %in% names(out)) out$receiver <- out$cluster_R
  out
}

load_seurat_with_axl <- function(seurat_path, axl_ps_path, patient_col,
                               condition_col, cond_levels, cell_type_col,
                               drop_multi) {
  message("Loading Seurat object: ", seurat_path)
  ext <- tolower(tools::file_ext(seurat_path))
  if (ext == "rds") {
    obj <- readRDS(seurat_path)
  } else if (ext %in% c("rdata", "rda")) {
    env <- new.env()
    load(seurat_path, envir = env)
    nms <- ls(env)
    if ("Kurten_HNSC" %in% nms) {
      obj <- get("Kurten_HNSC", envir = env)
    } else {
      obj <- get(nms[[1]], envir = env)
    }
  } else {
    stop("Unsupported Seurat file extension: ", ext)
  }

  message("Loading AXL patient summary: ", axl_ps_path)
  ps <- readRDS(axl_ps_path)
  if (!all(c("patient_id", condition_col) %in% names(ps))) {
    stop("Patient summary must contain 'patient_id' and ", condition_col)
  }
  mapping <- setNames(ps[[condition_col]], as.character(ps$patient_id))
  if (!patient_col %in% colnames(obj@meta.data)) {
    stop("Column not found in Seurat metadata: ", patient_col)
  }
  obj@meta.data[[condition_col]] <- mapping[as.character(obj@meta.data[[patient_col]])]

  keep <- !is.na(obj@meta.data[[cell_type_col]]) &
    obj@meta.data[[condition_col]] %in% cond_levels
  if (drop_multi) {
    keep <- keep & obj@meta.data[[cell_type_col]] != "Multi"
  }
  message("Subsetting to ", sum(keep), " cells (annotated, ", paste(cond_levels, collapse = "/"), ")")
  obj <- subset(obj, cells = colnames(obj)[keep])
  list(seurat = obj, patient_summary = ps)
}

build_cell_metadata <- function(seurat_obj, patient_col, cell_type_col,
                                condition_col) {
  md <- seurat_obj@meta.data
  data.frame(
    Cell_ID = colnames(seurat_obj),
    Cluster_ID = as.character(md[[cell_type_col]]),
    Condition_ID = as.character(md[[condition_col]]),
    Sample_ID = as.character(md[[patient_col]]),
    stringsAsFactors = FALSE
  )
}

qc_cell_types <- function(cell_metadata, min_cells, cond_levels) {
  ct_counts <- cell_metadata %>%
    group_by(Cluster_ID, Condition_ID) %>%
    summarise(n_cells = n(), .groups = "drop")

  valid <- ct_counts %>%
    group_by(Cluster_ID) %>%
    summarise(
      has_all_conds = all(cond_levels %in% Condition_ID),
      min_n = min(n_cells),
      .groups = "drop"
    ) %>%
    filter(has_all_conds, min_n >= min_cells)

  dropped <- setdiff(unique(cell_metadata$Cluster_ID), valid$Cluster_ID)
  if (length(dropped) > 0) {
    message("Dropping cell types failing min_cells (", min_cells, ") in both conditions: ",
            paste(dropped, collapse = ", "))
  }
  valid$Cluster_ID
}

write_manifest <- function(path, obj) {
  write(jsonlite::toJSON(obj, pretty = TRUE, auto_unbox = TRUE), path)
  invisible(obj)
}

# --- Paths ---------------------------------------------------------------------
stopifnot_dir(opt$output_dir)
paths <- list(
  rds = file.path(opt$output_dir, "Kurten_HNSC_scSeqCommDiff_res.rds"),
  diff_full = file.path(opt$output_dir, "Kurten_HNSC_differential_comm_full.tsv"),
  tumor_tsv = file.path(opt$output_dir, "differential_comm_Tumor_involved.tsv"),
  axl_all = file.path(opt$output_dir, paste0(opt$focus_gene, "_Tumor_differential_comm_all.tsv")),
  axl_sig = file.path(opt$output_dir, paste0(opt$focus_gene, "_Tumor_differential_comm_significant.tsv")),
  axl_summary = file.path(opt$output_dir, paste0(opt$focus_gene, "_Tumor_summary.txt")),
  cclens = file.path(opt$output_dir, "Kurten_HNSC_CClens_Tumor_diff_input.tsv"),
  manifest = file.path(opt$output_dir, "run_manifest.json"),
  session = file.path(opt$output_dir, "sessionInfo.txt"),
  qc_counts = file.path(opt$output_dir, "qc_celltype_condition_counts.tsv"),
  patient_balance = file.path(opt$output_dir, "patient_balance.tsv")
)

cond_levels <- c(opt$cond_high, opt$cond_low)
cond_names <- c(opt$cond_high, opt$cond_low)

# --- Load data -----------------------------------------------------------------
loaded <- load_seurat_with_axl(
  seurat_path = opt$seurat_path,
  axl_ps_path = opt$axl_patient_summary,
  patient_col = opt$patient_col,
  condition_col = opt$condition_col,
  cond_levels = cond_levels,
  cell_type_col = opt$cell_type_col,
  drop_multi = opt$drop_multi
)
seurat_obj <- loaded$seurat

cell_metadata <- build_cell_metadata(
  seurat_obj, opt$patient_col, opt$cell_type_col, opt$condition_col
)

valid_cell_types <- qc_cell_types(cell_metadata, opt$min_cells, cond_levels)
cell_metadata <- cell_metadata[cell_metadata$Cluster_ID %in% valid_cell_types, , drop = FALSE]
seurat_obj <- subset(seurat_obj, cells = cell_metadata$Cell_ID)

message("Final cells: ", nrow(cell_metadata), "; cell types: ", length(valid_cell_types))

write.table(
  cell_metadata %>% count(Cluster_ID, Condition_ID, Sample_ID, name = "n_cells"),
  paths$qc_counts, sep = "\t", row.names = FALSE, quote = FALSE
)

patient_balance <- cell_metadata %>%
  distinct(Sample_ID, Condition_ID) %>%
  count(Condition_ID, name = "n_patients")
write.table(patient_balance, paths$patient_balance, sep = "\t", row.names = FALSE, quote = FALSE)
message("Patients per condition:")
print(patient_balance)

# --- Expression matrix ---------------------------------------------------------
message("Extracting normalized expression (RNA/data) as dgCMatrix ...")
gene_expr <- SeuratObject::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
if (!inherits(gene_expr, "dgCMatrix")) {
  gene_expr <- as(gene_expr, "dgCMatrix")
}
cell_metadata <- cell_metadata[match(colnames(gene_expr), cell_metadata$Cell_ID), , drop = FALSE]
stopifnot(identical(cell_metadata$Cell_ID, colnames(gene_expr)))

# --- Run scSeqCommDiff ---------------------------------------------------------
if (!opt$skip_run) {
  if (!requireNamespace("scSeqComm", quietly = TRUE)) {
    stop(
      "Package 'scSeqComm' is not installed. Install with:\n",
      "  devtools::install_gitlab('sysbiobig/scseqcomm')"
    )
  }
  suppressPackageStartupMessages({
    library(scSeqComm)
    if (opt$n_cores > 1L) {
      if (!requireNamespace("doRNG", quietly = TRUE)) {
        stop("Install doRNG for parallel scSeqComm runs.")
      }
      library(doRNG)
    }
  })

  LR_db <- switch(
    opt$lr_db,
    ConnectomeDB_2020 = {
      data(LR_pairs_ConnectomeDB_2020, package = "scSeqComm", envir = environment())
      LR_pairs_ConnectomeDB_2020
    },
    Kumar_2018 = {
      data(LR_pairs_Kumar_2018, package = "scSeqComm", envir = environment())
      LR_pairs_Kumar_2018
    },
    stop("Unknown --lr_db: ", opt$lr_db)
  )

  data(TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High, package = "scSeqComm", envir = environment())
  data(TF_PPR_KEGG_human, package = "scSeqComm", envir = environment())
  TF_TG_db <- TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High
  TF_PPR <- TF_PPR_KEGG_human

  if (isTRUE(opt$subset_genes)) {
    lr_genes <- unique(c(LR_db$ligand, LR_db$receptor))
    lr_genes <- unique(unlist(strsplit(lr_genes, ",", fixed = TRUE)))
    trn_genes <- unique(unlist(TF_TG_db))
    keep_genes <- intersect(rownames(gene_expr), unique(c(lr_genes, trn_genes, opt$focus_gene)))
    message("Gene subsetting: ", length(keep_genes), " genes")
    gene_expr <- gene_expr[keep_genes, , drop = FALSE]
  }

  message("Running scSeqComm_differential (multi-sample, cond: ",
          paste(cond_names, collapse = " vs "), ") ...")
  scSeqCommDiff_res <- scSeqComm::scSeqComm_differential(
    gene_expr = gene_expr,
    cell_metadata = cell_metadata,
    scenario = "multi-sample",
    cond_names = cond_names,
    inter_signaling = TRUE,
    intra_signaling = TRUE,
    LR_pairs_DB = LR_db,
    inter_scores = "scSeqComm",
    TF_reg_DB = TF_TG_db,
    R_TF_association = TF_PPR,
    DEmethod = NULL,
    alternative_inter = "two.sided",
    only.pos = FALSE,
    N_cores = opt$n_cores,
    backend = "doParallel",
    bigmatrix = opt$bigmatrix,
    min_cells = opt$min_cells,
    count_thr = opt$count_thr
  )

  message("Saving full result to ", paths$rds)
  saveRDS(scSeqCommDiff_res, paths$rds)
} else {
  if (!file.exists(paths$rds)) {
    stop("--skip_run requires existing RDS at: ", paths$rds)
  }
  message("Loading existing result: ", paths$rds)
  scSeqCommDiff_res <- readRDS(paths$rds)
}

# --- Post-analysis -------------------------------------------------------------
diff_comm <- scSeqCommDiff_res$differential_comm
if (is.null(diff_comm) || nrow(diff_comm) == 0L) {
  stop("differential_comm is empty.")
}

message("Writing full differential_comm (", nrow(diff_comm), " rows)")
write.table(diff_comm, paths$diff_full, sep = "\t", row.names = FALSE, quote = FALSE)

tumor_comm <- filter_tumor_involved(diff_comm, opt$tumor_label)
message("Tumor-involved interactions: ", nrow(tumor_comm), " / ", nrow(diff_comm))
write.table(tumor_comm, paths$tumor_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

gene_mask <- row_involves_gene(tumor_comm, opt$focus_gene)
axl_comm <- tumor_comm[gene_mask, , drop = FALSE]
message(opt$focus_gene, "-involved (Tumor subset): ", nrow(axl_comm), " rows")
write.table(axl_comm, paths$axl_all, sep = "\t", row.names = FALSE, quote = FALSE)

p_col <- "pvalue_adj_S_inter"
es_col <- "effsize_S_inter"
sig_ok <- rep(TRUE, nrow(axl_comm))
if (p_col %in% names(axl_comm)) {
  sig_ok <- sig_ok & !is.na(axl_comm[[p_col]]) & axl_comm[[p_col]] < opt$p_adj_max
}
if (es_col %in% names(axl_comm)) {
  sig_ok <- sig_ok & !is.na(axl_comm[[es_col]]) & abs(axl_comm[[es_col]]) >= opt$effsize_min
}
axl_sig <- axl_comm[sig_ok, , drop = FALSE]
write.table(axl_sig, paths$axl_sig, sep = "\t", row.names = FALSE, quote = FALSE)

summary_lines <- c(
  paste0("scSeqCommDiff Kurten_HNSC POC summary"),
  paste0("Date: ", Sys.time()),
  paste0("Condition contrast: ", opt$cond_high, " vs ", opt$cond_low, " (multi-sample)"),
  paste0("Cells analyzed: ", nrow(cell_metadata)),
  paste0("differential_comm rows (full): ", nrow(diff_comm)),
  paste0("Tumor-involved rows: ", nrow(tumor_comm)),
  paste0(opt$focus_gene, "-involved (Tumor subset): ", nrow(axl_comm)),
  paste0(opt$focus_gene, " significant (p_adj<", opt$p_adj_max, ", |effsize|>=", opt$effsize_min, "): ", nrow(axl_sig)),
  "",
  "Top ", opt$focus_gene, " significant intercellular (by |effsize_S_inter|):"
)
if (nrow(axl_sig) > 0L && all(c("ligand", "receptor", "cluster_L", "cluster_R", es_col) %in% names(axl_sig))) {
  top <- axl_sig[order(-abs(axl_sig[[es_col]])), ]
  top <- head(top[, c("ligand", "receptor", "cluster_L", "cluster_R", p_col, es_col)], 15)
  summary_lines <- c(summary_lines, capture.output(print(top, row.names = FALSE)))
} else {
  summary_lines <- c(summary_lines, "(none)")
}
writeLines(summary_lines, paths$axl_summary)

# --- CClens export -------------------------------------------------------------
cclens_df <- switch(
  opt$cclens_source,
  tumor = tumor_comm,
  axl = axl_comm,
  axl_significant = axl_sig,
  stop("Unknown --cclens_source: ", opt$cclens_source)
)
cclens_out <- to_cclens_table(cclens_df)
write.table(cclens_out, paths$cclens, sep = "\t", row.names = FALSE, quote = FALSE)
message("CClens table (source=", opt$cclens_source, "): ", paths$cclens)

# --- Manifest & session --------------------------------------------------------
manifest <- list(
  seurat_path = normalizePath(opt$seurat_path, mustWork = TRUE),
  axl_patient_summary = normalizePath(opt$axl_patient_summary, mustWork = TRUE),
  output_dir = normalizePath(opt$output_dir, mustWork = TRUE),
  scenario = "multi-sample",
  cond_names = cond_names,
  n_cells = nrow(cell_metadata),
  n_cell_types = length(valid_cell_types),
  n_patients_high = sum(patient_balance$Condition_ID == opt$cond_high),
  n_patients_low = sum(patient_balance$Condition_ID == opt$cond_low),
  differential_comm_rows = nrow(diff_comm),
  tumor_involved_rows = nrow(tumor_comm),
  focus_gene = opt$focus_gene,
  focus_gene_rows = nrow(axl_comm),
  focus_gene_significant_rows = nrow(axl_sig),
  lr_db = opt$lr_db,
  min_cells = opt$min_cells,
  skip_run = opt$skip_run
)
write_manifest(paths$manifest, manifest)
writeLines(capture.output(sessionInfo()), paths$session)

message("Done. Outputs in: ", opt$output_dir)
