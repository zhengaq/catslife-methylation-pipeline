#!/usr/bin/env bash
# Shared settings for the SLURM batch scripts (sourced, not executed).
#   R_MODULE — module that puts Rscript on PATH (e.g. R/4.5.3); empty if already on PATH.
#   NPARTS   — stage-3 CpG chunks; must equal stage3.sbatch's --array size (asserted there).
NPARTS="${METHYL_NPARTS:-5}"
export METHYL_NPARTS="$NPARTS"

R_MODULE="${R_MODULE:-}"
if [ -n "$R_MODULE" ] && command -v module >/dev/null 2>&1; then module load "$R_MODULE"; fi
command -v Rscript >/dev/null 2>&1 || { echo "env.sh: Rscript not on PATH — set R_MODULE" >&2; exit 1; }
