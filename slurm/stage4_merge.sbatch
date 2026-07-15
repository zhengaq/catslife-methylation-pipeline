#!/usr/bin/env bash
#SBATCH --job-name=cats-meth-s4
#SBATCH --output=slurm/logs/stage4_%j.out
#SBATCH --error=slurm/logs/stage4_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
# Stage 4: merge the NPARTS chunks -> output/B.adjusted.platebatches.txt. Run after stage 3.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source slurm/env.sh

echo "stage4 merge (NPARTS=$NPARTS) start $(date) on $(hostname)"
Rscript 4.methylation_merge.chunked.R
echo "stage4 done  $(date)"
