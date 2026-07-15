#!/usr/bin/env bash
#SBATCH --job-name=cats-meth-s3
#SBATCH --output=slurm/logs/stage3_%A_%a.out
#SBATCH --error=slurm/logs/stage3_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --array=1-5
# Stage 3: one array task per CpG chunk. Each task loads the full dasen betas (~24G) then
# processes its slice, so size --mem for the load. --array size must equal NPARTS (env.sh).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source slurm/env.sh

: "${SLURM_ARRAY_TASK_ID:?run via sbatch --array=1-N, or set SLURM_ARRAY_TASK_ID}"
if [ "${SLURM_ARRAY_TASK_COUNT:-$NPARTS}" != "$NPARTS" ]; then
    echo "stage3: --array size (${SLURM_ARRAY_TASK_COUNT:-unset}) != NPARTS ($NPARTS)" >&2; exit 1
fi

echo "stage3 chunk ${SLURM_ARRAY_TASK_ID}/${NPARTS} start $(date) on $(hostname)"
Rscript 3.methylation_adjust.chunked.R --nparts "$NPARTS" --part "$SLURM_ARRAY_TASK_ID"
echo "stage3 chunk ${SLURM_ARRAY_TASK_ID}/${NPARTS} done  $(date)"
