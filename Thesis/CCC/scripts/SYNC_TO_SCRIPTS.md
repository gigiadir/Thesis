# Sync Thesis scripts to ~/Scripts (cluster)

Compute nodes often cannot see `Thesis/CCC/scripts` (symlink/Dropbox). Pipeline jobs are submitted from **GPFS** paths under `~/Scripts/`.

## What to sync after editing Thesis

| Thesis file | Copy to |
|-------------|---------|
| `scDiffCom/pipeline/scDiffComPipeline.R` | `~/Scripts/scDiffComPipeline.R` |
| `~/Scripts/run_with_scDiffComPipeline_env.sh` | Update separately if wrapper changes |

```bash
cp -f ~/Thesis/CCC/scripts/scDiffCom/pipeline/scDiffComPipeline.R \
      ~/Scripts/scDiffComPipeline.R
```

## Regenerating `submit_scDiffCom_jobs.sh`

`Main-scDiffComPipeline.R` writes `submit_scDiffCom_jobs.sh` next to `dirname(script_path)`.

- **Cluster workflow (default):** `script_path` points to `~/Scripts/scDiffComPipeline.R` → generated script is `~/Scripts/submit_scDiffCom_jobs.sh`.
- **Thesis-only test:** temporarily point `script_path` at `scDiffCom/pipeline/scDiffComPipeline.R` to generate under `pipeline/`, then copy to `~/Scripts/` if needed.

The committed `submit_scDiffCom_jobs.sh` in git is **ignored** (regenerate locally). See `.gitignore`.

## What not to sync

- `scDiffCom/post_analysis/` — interactive / report knitting
- `scDiffCom/preprocess/` — unless you run preprocess on cluster via qsub
- `tutorials/` — exploratory only
- `utils/` — sourced from Thesis paths in notebooks
