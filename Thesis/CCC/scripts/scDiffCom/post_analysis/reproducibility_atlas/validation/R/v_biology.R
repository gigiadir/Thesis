# Biological validation — controls and enrichment.

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

  atlas_genes <- {
    fdr_thr <- if (!is.null(ctx$cfg$fdr_threshold)) ctx$cfg$fdr_threshold else 0.05
    repro <- if (!is.null(ctx$repro_with_nulls) && "shuffle_FDR" %in% names(ctx$repro_with_nulls)) {
      ctx$repro_with_nulls
    } else {
      repro_df
    }
    if ("shuffle_FDR" %in% names(repro)) {
      repro$gene[repro$shuffle_FDR < fdr_thr]
    } else {
      repro$gene[repro$ReproScore > ctx$cfg$reproscore_threshold]
    }
  }
  if (requireNamespace("scDiffCom", quietly = TRUE)) {
    lri_go <- scDiffCom::LRI_human$LRI_curated_GO
    if (!is.null(lri_go) && nrow(lri_go) > 0) {
      sig_genes <- unique(unlist(strsplit(lri_go$LRI, ":")))
      atlas_in_sig <- sum(atlas_genes %in% sig_genes)
      bg_in_sig <- sum(gene_universe %in% sig_genes)
      n_atlas <- length(atlas_genes)
      n_bg <- length(gene_universe)
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
