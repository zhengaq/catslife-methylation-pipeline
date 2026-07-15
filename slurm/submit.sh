#!/usr/bin/env bash
# Submit stage 1 -> stage 3 (array) -> stage 4 as an afterok dependency chain, from the repo root.
#   bash slurm/submit.sh
#   METHYL_NPARTS=10 R_MODULE=R/4.5.3 bash slurm/submit.sh   # also set --array=1-10 in stage3.sbatch
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p slurm/logs

EXPORTS="ALL"
[ -n "${METHYL_NPARTS:-}" ] && EXPORTS="$EXPORTS,METHYL_NPARTS=$METHYL_NPARTS"
[ -n "${R_MODULE:-}" ]      && EXPORTS="$EXPORTS,R_MODULE=$R_MODULE"

jid1=$(sbatch --parsable --export="$EXPORTS" slurm/stage1.sbatch)
echo "stage1 submitted: $jid1"
jid3=$(sbatch --parsable --export="$EXPORTS" --dependency=afterok:"$jid1" slurm/stage3.sbatch)
echo "stage3 submitted: $jid3 (after $jid1)"
jid4=$(sbatch --parsable --export="$EXPORTS" --dependency=afterok:"$jid3" slurm/stage4_merge.sbatch)
echo "stage4 submitted: $jid4 (after $jid3)"
echo "Track: squeue -u \$USER"
