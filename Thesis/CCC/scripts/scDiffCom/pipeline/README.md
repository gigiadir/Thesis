# scDiffCom pipeline

| File | Role |
|------|------|
| `Main-scDiffComPipeline.R` | Build `submit_scDiffCom_jobs.sh` qsub commands |
| `scDiffComPipeline.R` | Per-gene scDiffCom runner (sync copy to `~/Scripts/`) |
| `run_scSeqCommDiff_Kurten_HNSC.R` | scSeqCommDiff POC (Kurten HNSC) |

`submit_scDiffCom_jobs.sh` is **gitignored** — regenerate after editing the gene/dataset list in `Main-scDiffComPipeline.R`. See [../../SYNC_TO_SCRIPTS.md](../../SYNC_TO_SCRIPTS.md).
