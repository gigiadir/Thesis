# Cross-Cohort CCI Reproducibility Atlas

## Layout

| Path | Purpose |
|------|---------|
| [`index.Rmd`](index.Rmd) | Interactive driver — run section chunks in order |
| [`sections/*.Rmd`](sections/) | **Stage code lives here** (inspect, vocab, tensor, …) |
| [`R/atlas_helpers.R`](R/atlas_helpers.R) | Shared functions (ReproScore, tensor helpers) |
| [`R/stage_05_nulls.R`](R/stage_05_nulls.R) | Stage 5 null permutations (shared by Rmd + batch script) |
| [`R/00_atlas_setup.R`](R/00_atlas_setup.R) | Session init, checkpoints, `knit.atlas.section()` |
| [`run_stage_05_nulls.R`](run_stage_05_nulls.R) | Standalone Stage 5 batch driver (resumable) |
| [`submit_stage05_nulls.sh`](submit_stage05_nulls.sh) | SGE `qsub` wrapper for remote jobs |
| [`../R/01_load_filter.R`](../R/01_load_filter.R) | Load/filter scDiffCom (shared with post_analysis) |

## Interactive

Open [`index.Rmd`](index.Rmd) or any file under [`sections/`](sections/) in RStudio. Each section contains the full stage logic; standalone sections load the prior checkpoint if `atlas_env` is missing.

## CLI

```bash
Rscript run_atlas.R      # all sections
Rscript run_atlas.R 0    # one section (00_inspect.Rmd)
Rscript run_atlas.R 5    # Stage 5 via knitr (same core logic as batch script)
Rscript run_atlas.R diag
```

## Stage 5 batch (remote job)

Stage 5 (`n_perm=1000`) is slow. Use the standalone script — slim init (no knitr/scDiffCom), resumable checkpoints, parallel `mclapply` (max 50 cores).

**Prerequisite:** `results/stage04_atlas_env.rds` must exist.

### Direct Rscript

```bash
cd ~/Thesis/CCC/scripts/scDiffCom/post_analysis/reproducibility_atlas

# Full run (resume enabled by default)
Rscript run_stage_05_nulls.R

# Quick test
Rscript run_stage_05_nulls.R --n-perm 100 --ncores 4

# Resume after job kill (automatic if null_perm_checkpoint.rds exists)
Rscript run_stage_05_nulls.R

# Force restart
Rscript run_stage_05_nulls.R --no-resume
```

### SGE job (recommended)

```bash
cd ~/Thesis/CCC/scripts/scDiffCom/post_analysis/reproducibility_atlas
bash submit_stage05_nulls.sh

# More cores (match NCORES to -pe shared N)
NCORES=16 bash submit_stage05_nulls.sh

# Override parallel environment (cluster default: shared; use none to omit -pe)
SGE_PE=shared NCORES=16 bash submit_stage05_nulls.sh
SGE_PE=none bash submit_stage05_nulls.sh
```

Logs: `$HOME/CCC-scDiffCom/results/reproducibility_atlas/logs/atlas_nulls.o*` / `.e*`

**Resource guidance:**
- Parallel env: `SGE_PE=shared` (N1GE on bhn1095); list PEs with `qconf -spl`
- Cores: up to `null_max_cores` (default 50); set `NCORES` and match `-pe shared N`
- Memory: ~8–16 GB
- Wall time: hours for 1000 perms; resume checkpoints saved every `null_checkpoint_every` (default 50)

### Stage 5 outputs

| File | Content |
|------|---------|
| `repro_scores_with_nulls.tsv` | `shuffle_p`, `shuffle_FDR` added |
| `null_reproscore_matrix.rds` | G × n_perm null matrix |
| `global_null.txt` | Global permutation p-value |
| `stage05_atlas_env.rds` | Checkpoint for Stage 6 |
| `null_perm_checkpoint.rds` | Temporary resume file (deleted on success) |

After Stage 5: `Rscript run_atlas.R 6` or open `sections/06_atlas.Rmd`.

### Sync Stage 5 outputs (batch → interactive)

Stage 5 batch uses [`config_batch.yml`](config_batch.yml) (`output_dir` on GPFS). Interactive sections use [`config.yml`](config.yml) (`output_dir` under `~/Thesis/CCC/outputs/...`). After a batch run, copy these four files into the Thesis `results/` directory:

```bash
SRC=~/CCC-scDiffCom/results/reproducibility_atlas/v1-v2/results
DST=~/Thesis/CCC/outputs/scDiffCom/reproducibility_atlas/v1-v2/results
cp "$SRC"/{stage05_atlas_env.rds,repro_scores_with_nulls.tsv,null_reproscore_matrix.rds,global_null.txt} "$DST"/
```

Then patch `output_dir` inside `stage05_atlas_env.rds` to match `config.yml` (batch checkpoints store the GPFS path):

```bash
Rscript -e 'p <- path.expand("~/Thesis/CCC/outputs/scDiffCom/reproducibility_atlas/v1-v2/results/stage05_atlas_env.rds"); e <- readRDS(p); o <- path.expand("~/Thesis/CCC/outputs/scDiffCom/reproducibility_atlas/v1-v2"); e$output_dir <- o; if (!is.null(e$cfg)) e$cfg$output_dir <- o; saveRDS(e, p)'
```

Alternatively, if `null_perm_checkpoint.rds` still exists on GPFS, [`recover_stage05_outputs.R`](recover_stage05_outputs.R) can rebuild outputs locally — but a straight `cp` is faster when the batch run already finished.

## Config

[`config.yml`](config.yml) — copied to `results/config_used.yml` on first run.

Stage 5 keys: `n_perm`, `null_checkpoint_every`, `null_max_cores`, `seed`.

Stage 6 atlas membership: `shuffle_FDR < fdr_threshold` (default 0.05). ReproScore ranks members within the atlas.
