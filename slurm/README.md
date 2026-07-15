# SLURM batch scripts (stage 1, 3, 4)

Run the memory-/compute-heavy pipeline stages as batch jobs instead of an interactive
RStudio session. Stage 2 (PCA) is a standalone diagnostic and is left out of the chain;
run it interactively if you want it.

## One-time setup

1. Fill in `config.site.R` at the repo root (paths for the delivery) — the jobs read it
   automatically via the `.methyl-root` marker, so **submit from the repo root**.
2. In `slurm/env.sh`, set `R_MODULE` to your cluster's R module (e.g. `R/4.5.3`), or leave
   it empty if `Rscript` is already on the batch PATH.
3. Pick `NPARTS` (stage-3 CpG chunks). Default is 5. If you change it, change it in **two**
   places: `NPARTS` in `slurm/env.sh` (or `export METHYL_NPARTS=...`) **and** `--array=1-N`
   in `slurm/stage3.sbatch`. `stage3.sbatch` asserts they agree and fails fast if not.

## Run the whole chain

```bash
bash slurm/submit.sh              # stage1 -> stage3 array -> stage4, wired with afterok
# or with overrides:
METHYL_NPARTS=10 R_MODULE=R/4.5.3 bash slurm/submit.sh   # (also set --array=1-10 in stage3.sbatch)
```

## Or submit stages individually

```bash
sbatch slurm/stage1.sbatch
sbatch --dependency=afterok:<stage1_jobid> slurm/stage3.sbatch
sbatch --dependency=afterok:<stage3_jobid> slurm/stage4_merge.sbatch
```

## Resources (tune per node)

| Stage | Default `--mem` | `--time` | Notes |
|-------|-----------------|----------|-------|
| 1     | 96G  | 12h | read + noob + dasen at ~1600 samples; detectionP chunked (`METHYL_DETP_CHUNK`, default 200) |
| 3     | 48G/task | 24h | each array task loads the full dasen betas (~24G) then does its CpG slice |
| 4     | 32G  | 2h  | merge NPARTS chunks -> `output/B.adjusted.platebatches.txt` |

Logs land in `slurm/logs/`. Track with `squeue -u $USER`.
