#!/usr/bin/env Rscript

# Figure 3.3: Partition-rule label concordance (H&N).
# Mean fraction of patients with concordant HIGH or LOW labels between rule pairs,
# averaged across driver-panel genes.
#
# Requires preprocess outputs for all four rules on the same cohort.
#
# Example:
#   Rscript plotSplitRuleConcordance.R --dataset_name Kurten_HNSC
#
# Requires: optparse, ggplot2

suppressPackageStartupMessages({
  library(optparse)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages(\"ggplot2\").")
  }
  library(ggplot2)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(".", winslash = "/")
}
source(file.path(script_dir, "scDiffComGenePanel.R"))
source(file.path(script_dir, "patientSplitUtils.R"))

EXTREME_LABELS <- c("LOW", "HIGH")

RULE_SPECS <- list(
  ExpressionQuantile = list(
    dir_opt = "expression_quantile_dir",
    default = "~/CCC-PreProcess/results-ExpressionQuantile"
  ),
  RankGenes = list(
    dir_opt = "rankgenes_dir",
    default = "~/CCC-PreProcess/results-RankGenes"
  ),
  Residual = list(
    dir_opt = "residual_dir",
    default = "~/CCC-PreProcess/results-Residual"
  ),
  PatientZScore = list(
    dir_opt = "patient_zscore_dir",
    default = "~/CCC-PreProcess/results-Patient-ZScore"
  )
)

option_list <- list(
  make_option("--dataset_name", type = "character", default = "Kurten_HNSC",
              help = "Dataset / cohort name [default %default]", metavar = "character"),
  make_option("--expression_quantile_dir", type = "character",
              default = RULE_SPECS$ExpressionQuantile$default,
              help = "ExpressionQuantile output base [default %default]"),
  make_option("--rankgenes_dir", type = "character",
              default = RULE_SPECS$RankGenes$default,
              help = "RankGenes output base [default %default]"),
  make_option("--residual_dir", type = "character",
              default = RULE_SPECS$Residual$default,
              help = "Residual output base [default %default]"),
  make_option("--patient_zscore_dir", type = "character",
              default = RULE_SPECS$PatientZScore$default,
              help = "PatientZScore output base [default %default]"),
  make_option("--gene_list", type = "character", default = NULL,
              help = "Optional .rds or .txt gene list (default: scDiffCom driver panel)",
              metavar = "character"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory [default ~/Thesis/CCC/outputs/split_rule_concordance/{dataset}]",
              metavar = "character"),
  make_option("--skip_plot", action = "store_true", default = FALSE,
              help = "Write tables only, skip figure [default %default]")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

dataset_name <- opt$dataset_name
if (is.null(dataset_name) || !nzchar(dataset_name)) {
  stop("--dataset_name is required.")
}

rule_dirs <- list(
  ExpressionQuantile = path.expand(opt$expression_quantile_dir),
  RankGenes = path.expand(opt$rankgenes_dir),
  Residual = path.expand(opt$residual_dir),
  PatientZScore = path.expand(opt$patient_zscore_dir)
)

out_dir <- if (!is.null(opt$output_dir) && nzchar(opt$output_dir)) {
  path.expand(opt$output_dir)
} else {
  path.expand(file.path(
    "~/Thesis/CCC/outputs/split_rule_concordance",
    dataset_name
  ))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_genes <- if (!is.null(opt$gene_list) && nzchar(opt$gene_list)) {
  load_gene_list(opt$gene_list)
} else {
  SCDIFFCOM_GENE_PANEL
}

load_rule_splits <- function(rule_name, base_dir, dataset_name, genes) {
  ds_dir <- file.path(base_dir, dataset_name)
  if (!dir.exists(ds_dir)) {
    stop("Directory not found for rule ", rule_name, ": ", ds_dir)
  }
  splits <- load_splits_from_grouped_dir(base_dir, dataset_name, genes = genes)
  if (is.null(splits) || nrow(splits) == 0L) {
    stop("No splits loaded for rule ", rule_name, " in ", ds_dir)
  }
  splits
}

align_splits_list <- function(splits_list) {
  rule_names <- names(splits_list)
  genes <- Reduce(intersect, lapply(splits_list, rownames))
  if (length(genes) == 0L) {
    stop("No genes shared across all split rules.")
  }

  patients <- Reduce(intersect, lapply(splits_list, colnames))
  if (length(patients) == 0L) {
    stop("No patients shared across all split rules.")
  }

  aligned <- lapply(splits_list, function(m) {
    m[genes, patients, drop = FALSE]
  })
  list(splits = aligned, genes = genes, patients = patients)
}

gene_pair_extreme_concordance <- function(lab_a, lab_b) {
  extreme_a <- lab_a %in% EXTREME_LABELS
  extreme_b <- lab_b %in% EXTREME_LABELS
  comparable <- extreme_a & extreme_b
  n_comp <- sum(comparable, na.rm = TRUE)
  if (n_comp == 0L) {
    return(list(
      concordance = NA_real_,
      n_comparable = 0L,
      n_concordant = 0L
    ))
  }
  la <- lab_a[comparable]
  lb <- lab_b[comparable]
  n_conc <- sum(la == lb, na.rm = TRUE)
  list(
    concordance = n_conc / n_comp,
    n_comparable = n_comp,
    n_concordant = n_conc
  )
}

compute_per_gene_concordance <- function(aligned_splits, rule_a, rule_b) {
  genes <- rownames(aligned_splits[[rule_a]])
  mat_a <- aligned_splits[[rule_a]]
  mat_b <- aligned_splits[[rule_b]]

  rows <- lapply(genes, function(g) {
    res <- gene_pair_extreme_concordance(mat_a[g, ], mat_b[g, ])
    data.frame(
      gene = g,
      rule_a = rule_a,
      rule_b = rule_b,
      concordance = res$concordance,
      n_comparable = res$n_comparable,
      n_concordant = res$n_concordant,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_pairwise <- function(per_gene_df, rule_names) {
  pairs <- combn(rule_names, 2, simplify = FALSE)
  pair_rows <- lapply(pairs, function(p) {
    sub <- per_gene_df[per_gene_df$rule_a == p[[1]] & per_gene_df$rule_b == p[[2]], ]
    ok <- is.finite(sub$concordance)
    data.frame(
      rule_a = p[[1]],
      rule_b = p[[2]],
      mean_concordance = mean(sub$concordance[ok], na.rm = TRUE),
      n_genes = nrow(sub),
      n_genes_used = sum(ok),
      mean_n_comparable = mean(sub$n_comparable[ok], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pair_rows)
}

build_summary_matrix <- function(summary_df, rule_names) {
  n <- length(rule_names)
  M <- matrix(NA_real_, n, n, dimnames = list(rule_names, rule_names))
  diag(M) <- 1

  for (i in seq_len(nrow(summary_df))) {
    a <- summary_df$rule_a[[i]]
    b <- summary_df$rule_b[[i]]
    v <- summary_df$mean_concordance[[i]]
    M[a, b] <- v
    M[b, a] <- v
  }
  M
}

run_sanity_checks <- function(M, summary_df, per_gene_df) {
  rule_names <- rownames(M)
  n <- length(rule_names)

  if (!all(is.finite(diag(M)))) {
    stop("Sanity check failed: diagonal is not all 1.")
  }
  if (abs(max(diag(M)) - 1) > 1e-9 || abs(min(diag(M)) - 1) > 1e-9) {
    stop("Sanity check failed: diagonal values are not 1.")
  }

  off <- M[upper.tri(M)]
  off <- off[is.finite(off)]
  if (length(off) > 0L && (min(off) < -1e-9 || max(off) > 1 + 1e-9)) {
    stop("Sanity check failed: off-diagonal concordance outside [0, 1].")
  }

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      if (!is.finite(M[i, j]) || !is.finite(M[j, i])) next
      if (abs(M[i, j] - M[j, i]) > 1e-9) {
        stop(
          "Sanity check failed: matrix not symmetric (",
          rule_names[[i]], " vs ", rule_names[[j]], ")."
        )
      }
    }
  }

  # Spot-check first gene, first off-diagonal pair
  if (n >= 2L && nrow(per_gene_df) > 0L) {
    g0 <- per_gene_df$gene[[1]]
    ra <- rule_names[[1]]
    rb <- rule_names[[2]]
    sub <- per_gene_df[
      per_gene_df$gene == g0 &
        per_gene_df$rule_a == ra &
        per_gene_df$rule_b == rb,
    ]
    if (nrow(sub) == 1L && is.finite(sub$concordance[[1]])) {
      la <- aligned[[ra]][g0, ]
      lb <- aligned[[rb]][g0, ]
      manual <- gene_pair_extreme_concordance(la, lb)
      if (abs(manual$concordance - sub$concordance[[1]]) > 1e-9) {
        stop("Sanity check failed: manual spot-check mismatch for gene ", g0)
      }
    }
  }

  message("Sanity checks passed (symmetry, diagonal=1, range, spot-check).")
  invisible(TRUE)
}

plot_concordance_heatmap <- function(M, dataset_name, out_png, out_pdf) {
  rule_names <- rownames(M)
  M_plot <- M
  diag(M_plot) <- NA_real_
  long <- expand.grid(
    rule_a = factor(rule_names, levels = rule_names),
    rule_b = factor(rule_names, levels = rev(rule_names)),
    stringsAsFactors = TRUE
  )
  long$concordance <- mapply(
    function(a, b) M_plot[as.character(a), as.character(b)],
    as.character(long$rule_a),
    as.character(long$rule_b)
  )
  low_bound <- 0.5
  high_bound <- 1
  subtitle_text <- paste(
    dataset_name,
    "mean fraction of patients with concordant HIGH/LOW labels",
    "between rule pairs, averaged across driver-panel genes",
    sep = "\n"
  )

  p <- ggplot(long, aes(x = rule_a, y = rule_b, fill = concordance)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = ifelse(is.finite(concordance), sprintf("%.2f", concordance), "")),
      color = "black",
      size = 4
    ) +
    scale_fill_gradient(
      limits = c(low_bound, high_bound),
      low = "#e9f0f8",
      high = "#2166ac",
      name = "Mean\nconcordance",
      na.value = "grey90"
    ) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title.position = "plot",
      plot.subtitle = element_text(lineheight = 1.05, margin = margin(b = 8)),
      plot.margin = margin(t = 8, r = 20, b = 8, l = 8)
    ) +
    labs(
      title = "Partition-rule label concordance (H&N)",
      subtitle = subtitle_text,
      x = "Rule",
      y = "Rule"
    )

  ggsave(out_png, plot = p, width = 9, height = 7.2, dpi = 300)
  ggsave(out_pdf, plot = p, width = 9, height = 7.2)
  invisible(p)
}

write_run_metadata <- function(path, opt, rule_dirs, genes, patients, n_genes_used) {
  meta <- data.frame(
    key = c(
      "dataset_name",
      "n_panel_genes_requested",
      "n_genes_aligned",
      "n_patients_aligned",
      "n_genes_with_any_pair",
      "expression_quantile_dir",
      "rankgenes_dir",
      "residual_dir",
      "patient_zscore_dir",
      "gene_list",
      "metric",
      "timestamp"
    ),
    value = c(
      dataset_name,
      as.character(length(panel_genes)),
      as.character(length(genes)),
      as.character(length(patients)),
      as.character(n_genes_used),
      rule_dirs$ExpressionQuantile,
      rule_dirs$RankGenes,
      rule_dirs$Residual,
      rule_dirs$PatientZScore,
      if (!is.null(opt$gene_list) && nzchar(opt$gene_list)) opt$gene_list else "SCDIFFCOM_GENE_PANEL",
      "HIGH/LOW concordance on comparable extreme-labeled patients; mean across genes",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(meta, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- main ---

message("Loading splits for ", dataset_name, " (", length(panel_genes), " panel genes) ...")
splits_list <- lapply(names(rule_dirs), function(rn) {
  message("  ", rn, " ...")
  load_rule_splits(rn, rule_dirs[[rn]], dataset_name, panel_genes)
})
names(splits_list) <- names(rule_dirs)

aligned_result <- align_splits_list(splits_list)
aligned <- aligned_result$splits
genes <- aligned_result$genes
patients <- aligned_result$patients

missing_genes <- setdiff(panel_genes, genes)
if (length(missing_genes) > 0L) {
  warning(
    length(missing_genes), " panel gene(s) missing from one or more rules; using ",
    length(genes), " aligned genes.",
    call. = FALSE, immediate. = TRUE
  )
}

message(
  "Aligned: ", length(genes), " genes x ", length(patients), " patients across ",
  length(rule_dirs), " rules."
)

rule_names <- names(aligned)
pair_list <- combn(rule_names, 2, simplify = FALSE)
per_gene_parts <- lapply(pair_list, function(p) {
  compute_per_gene_concordance(aligned, p[[1]], p[[2]])
})
per_gene_df <- do.call(rbind, per_gene_parts)

summary_df <- summarize_pairwise(per_gene_df, rule_names)
M <- build_summary_matrix(summary_df, rule_names)

run_sanity_checks(M, summary_df, per_gene_df)

per_gene_path <- file.path(out_dir, "split_rule_concordance_per_gene.csv")
summary_path <- file.path(out_dir, "split_rule_concordance_summary.csv")
matrix_path <- file.path(out_dir, "split_rule_concordance_matrix.csv")
meta_path <- file.path(out_dir, "split_rule_concordance_run_metadata.tsv")

utils::write.csv(per_gene_df, per_gene_path, row.names = FALSE)
utils::write.csv(summary_df, summary_path, row.names = FALSE)
utils::write.csv(
  data.frame(rule = rownames(M), M, check.names = FALSE),
  matrix_path,
  row.names = FALSE
)

n_genes_used <- min(summary_df$n_genes_used, na.rm = TRUE)
write_run_metadata(meta_path, opt, rule_dirs, genes, patients, n_genes_used)

message("Wrote ", per_gene_path)
message("Wrote ", summary_path)
message("Wrote ", matrix_path)
message("Wrote ", meta_path)

if (!isTRUE(opt$skip_plot)) {
  png_path <- file.path(out_dir, "split_rule_concordance_heatmap.png")
  pdf_path <- file.path(out_dir, "split_rule_concordance_heatmap.pdf")
  plot_concordance_heatmap(M, dataset_name, png_path, pdf_path)
  message("Wrote ", png_path)
  message("Wrote ", pdf_path)
}

message("Pairwise mean concordance:")
for (i in seq_len(nrow(summary_df))) {
  message(
    "  ", summary_df$rule_a[[i]], " vs ", summary_df$rule_b[[i]], ": ",
    sprintf("%.3f", summary_df$mean_concordance[[i]]),
    " (n_genes_used=", summary_df$n_genes_used[[i]], ")"
  )
}

message("Done. Outputs in: ", out_dir)
