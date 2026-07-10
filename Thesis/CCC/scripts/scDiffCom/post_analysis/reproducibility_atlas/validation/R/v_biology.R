# Biological validation — IDR concordance, controls, enrichment.

compute_idr_summary <- function(ctx) {
  malignant_by_cohort <- ctx$malignant_by_cohort
  cohorts <- ctx$cohorts
  gene_universe <- ctx$gene_universe
  idx_pairs <- ctx$cohort_pairs_idx
  pair_labels <- ctx$pair_labels
  ranking <- ctx$cfg$idr_ranking
  threshold <- ctx$cfg$idr_threshold
  n_pairs <- length(idx_pairs)

  ensure_idr()

  idr_pass_fraction <- vapply(gene_universe, function(g) {
    ccis_g <- sort(unique(unlist(lapply(cohorts, function(ds) {
      malignant_by_cohort[[ds]][[g]]$CCI
    }))))
    if (length(ccis_g) < 10) return(NA_real_)

    vecs <- lapply(cohorts, function(ds) {
      gene_rank_values(malignant_by_cohort[[ds]][[g]], ccis = ccis_g, ranking = ranking)
    })
    cci_pair_pass <- setNames(rep(0L, length(ccis_g)), ccis_g)
    for (p in seq_len(n_pairs)) {
      a <- idx_pairs[[p]][1]
      b <- idx_pairs[[p]][2]
      res <- run_pairwise_idr(vecs[[a]], vecs[[b]], threshold = threshold)
      if (length(res$pass) > 0) cci_pair_pass[res$pass] <- cci_pair_pass[res$pass] + 1L
    }
    pass_frac_by_cci <- cci_pair_pass / n_pairs
    mean(pass_frac_by_cci >= 0.5)
  }, numeric(1))

  data.frame(
    gene = gene_universe,
    idr_pass_fraction = idr_pass_fraction,
    stringsAsFactors = FALSE
  )
}

DEFAULT_POSITIVE_CONTROLS <- c(
  "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
  "VEGFA", "FGF2", "PDGFA", "PDGFB", "CXCL12", "CCL2", "IL6"
)

DEFAULT_NEGATIVE_CONTROLS <- c(
  "GAPDH", "ACTB", "B2M", "RPL13A", "RPLP0"
)

run_v_biology <- function(ctx, validation_dir, controls = NULL) {
  repro_df <- ctx$repro_df
  gene_universe <- ctx$gene_universe

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' required for IDR validation")
  }
  suppressPackageStartupMessages(library(dplyr))

  idr_df <- compute_idr_summary(ctx)
  merged <- merge(repro_df, idr_df, by = "gene")
  rho <- cor(merged$ReproScore, merged$idr_pass_fraction,
             method = "spearman", use = "complete.obs")

  png(file.path(validation_dir, "results", "reproscore_vs_idr.png"),
      width = 600, height = 500)
  plot(merged$idr_pass_fraction, merged$ReproScore, pch = 16, cex = 0.6,
       xlab = "IDR pass fraction", ylab = "ReproScore",
       main = sprintf("ReproScore vs IDR (rho=%.3f)", rho))
  abline(h = 0.5, col = "grey", lty = 2)
  dev.off()

  merged$repro_rank <- rank(-merged$ReproScore, ties.method = "average")
  merged$idr_rank <- rank(-merged$idr_pass_fraction, ties.method = "average")
  merged$rank_diff <- abs(merged$repro_rank - merged$idr_rank)
  disagree <- merged[order(-merged$rank_diff), ][seq_len(min(20, nrow(merged))), ]
  readr::write_tsv(disagree[, c("gene", "ReproScore", "idr_pass_fraction", "repro_rank", "idr_rank")],
                   file.path(validation_dir, "results", "method_disagreement.tsv"))

  if (is.finite(rho) && rho > 0.2) {
    append_verdict(validation_dir, "reproscore_vs_idr", 6, "PASS", round(rho, 4),
                   "positive concordance", "Two methods agree")
  } else if (is.finite(rho) && rho > 0) {
    append_verdict(validation_dir, "reproscore_vs_idr", 6, "WARN", round(rho, 4),
                   "positive concordance", "Weak concordance; methods may measure different things")
  } else {
    append_verdict(validation_dir, "reproscore_vs_idr", 6, "WARN", round(rho, 4),
                   "positive concordance", "Near-zero concordance")
  }

  if (is.null(controls)) {
    pos <- intersect(DEFAULT_POSITIVE_CONTROLS, gene_universe)
    neg <- intersect(DEFAULT_NEGATIVE_CONTROLS, gene_universe)
    controls <- list(positive = pos, negative = neg)
    message("Using default control genes (positive: ", paste(pos, collapse = ", "),
            "; negative: ", paste(neg, collapse = ", "), ")")
    message("Override with --controls-file to customize after user review.")
  }

  all_ctrl <- c(controls$positive, controls$negative)
  ctrl_df <- repro_df[repro_df$gene %in% all_ctrl, ]
  ctrl_df$control_type <- ifelse(ctrl_df$gene %in% controls$positive, "positive", "negative")
  ctrl_df$rank <- rank(-ctrl_df$ReproScore, ties.method = "average")
  ctrl_df <- ctrl_df[order(ctrl_df$rank), ]
  readr::write_tsv(ctrl_df, file.path(validation_dir, "results", "control_ranks.tsv"))

  pos_ranks <- ctrl_df$rank[ctrl_df$control_type == "positive"]
  neg_ranks <- ctrl_df$rank[ctrl_df$control_type == "negative"]
  n_genes <- nrow(repro_df)
  pos_ok <- length(pos_ranks) == 0 || mean(pos_ranks) < n_genes * 0.5
  neg_ok <- length(neg_ranks) == 0 || mean(neg_ranks) > n_genes * 0.5

  if (pos_ok && neg_ok && length(pos_ranks) > 0) {
    append_verdict(validation_dir, "controls", 7, "PASS",
                   sprintf("pos_mean_rank=%.0f", mean(pos_ranks)),
                   "positives high, negatives low", "Canonical controls rank as expected")
  } else {
    append_verdict(validation_dir, "controls", 7, "WARN",
                   sprintf("pos_mean_rank=%.0f,neg_mean_rank=%.0f",
                           if (length(pos_ranks)) mean(pos_ranks) else NA,
                           if (length(neg_ranks)) mean(neg_ranks) else NA),
                   "positives high, negatives low",
                   "Known biology not ranking up; trust pipeline less")
  }

  atlas_genes <- repro_df$gene[repro_df$ReproScore > ctx$cfg$reproscore_threshold]
  if (requireNamespace("scDiffCom", quietly = TRUE)) {
    lri_go <- scDiffCom::LRI_human$LRI_curated_GO
    if (!is.null(lri_go) && nrow(lri_go) > 0) {
      sig_genes <- unique(unlist(strsplit(lri_go$LRI, ":")))
      atlas_in_sig <- sum(atlas_genes %in% sig_genes)
      bg_in_sig <- sum(gene_universe %in% sig_genes)
      n_atlas <- length(atlas_genes)
      n_bg <- length(gene_universe)
      n_sig <- length(sig_genes)
      expected <- n_atlas * bg_in_sig / n_bg
      p_val <- stats::phyper(atlas_in_sig - 1, bg_in_sig, n_bg - bg_in_sig, n_atlas, lower.tail = FALSE)

      enrich_df <- data.frame(
        term = "LRI_curated_signaling",
        atlas_n = n_atlas,
        atlas_in_term = atlas_in_sig,
        background_in_term = bg_in_sig,
        background_n = n_bg,
        expected = expected,
        p_value = p_val,
        stringsAsFactors = FALSE
      )
      readr::write_tsv(enrich_df, file.path(validation_dir, "results", "atlas_enrichment.tsv"))

      if (p_val < 0.05) {
        append_verdict(validation_dir, "enrichment", 7, "PASS", round(p_val, 4),
                       "p < 0.05", "Atlas enriched for LRI signaling annotations")
      } else {
        append_verdict(validation_dir, "enrichment", 7, "INFO", round(p_val, 4),
                       "p < 0.05", "No significant LRI enrichment; atlas may still be valid")
      }
    } else {
      append_verdict(validation_dir, "enrichment", 7, "INFO", NA,
                     "p < 0.05", "LRI_curated_GO not available")
    }
  } else {
    append_verdict(validation_dir, "enrichment", 7, "INFO", NA,
                   "p < 0.05", "scDiffCom package not available for enrichment")
  }

  ctx
}
