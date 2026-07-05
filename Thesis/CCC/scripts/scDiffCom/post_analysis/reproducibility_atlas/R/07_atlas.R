# Stage 7 — assemble final atlas tables.

run_stage_07_atlas <- function(atlas_env) {
  cfg <- atlas_env$cfg
  output_dir <- atlas_env$output_dir

  message("Stage 7: assemble atlas")

  atlas <- atlas_env$repro_df %>%
    left_join(atlas_env$idr_df, by = "gene") %>%
    mutate(
      atlas_member = ReproScore > cfg$reproscore_threshold &
        shuffle_FDR < cfg$fdr_threshold &
        !is.na(idr_pass_fraction)
    )

  readr::write_csv(atlas, file.path(output_dir, "results", "atlas.csv"))

  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(atlas, file.path(output_dir, "results", "atlas.parquet"))
  }

  n_members <- sum(atlas$atlas_member, na.rm = TRUE)
  message(sprintf("  atlas.csv: %d genes, %d atlas members", nrow(atlas), n_members))

  atlas_env$atlas <- atlas
  saveRDS(atlas_env, file.path(output_dir, "results", "stage07_atlas_env.rds"))

  invisible(atlas_env)
}
