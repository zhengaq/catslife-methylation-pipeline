#!/usr/bin/env bash
#SBATCH --job-name=cats-meth-s1
#SBATCH --output=slurm/logs/stage1_%j.out
#SBATCH --error=slurm/logs/stage1_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
# Stage 1: read IDATs -> detP QC -> noob -> dasen. Submit from the repo root.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source slurm/env.sh

# Skip the unused stage-1 checkpoints by default (disk); override with METHYL_SAVE_INTERMEDIATES=TRUE.
export METHYL_SAVE_INTERMEDIATES="${METHYL_SAVE_INTERMEDIATES:-FALSE}"

echo "stage1 start $(date) on $(hostname); DETP_CHUNK=${METHYL_DETP_CHUNK:-200} SAVE_INTERMEDIATES=$METHYL_SAVE_INTERMEDIATES"
Rscript 1.methylation_minfi.R
echo "stage1 done  $(date)"
