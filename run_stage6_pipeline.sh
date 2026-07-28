#!/usr/bin/env bash
#
# run_stage6_pipeline.sh — run the stage-6 sensitivity & validity checks with checkpoint/resume:
#   1. validity  stage6/validity_clocks.R   (clock-table validity/sensitivity; needs stage-5 mAge_clocks.csv)
#   2. pca       stage6/pca_sex_batch.R      (sex-chromosome / batch-structure PCA; needs stage-1 betas + minfi)
#
# Checkpoint/resume: each check that finishes writes logs/.ckpt/<check>.done. A re-run skips any
# check whose marker exists, so after a failure you just launch again and it resumes at the failed
# check. Unlike the stage-5 pipeline the checks are INDEPENDENT and not fail-fast: a failed check is
# logged and the pipeline moves on (only successful checks checkpoint), so a broken dependency for
# one check never blocks the other. Each check is streamed to the console and to
# logs/stage6_pipeline_<check>.log; the run exits non-zero if any check failed.
#
# Paths (inputs/outputs) are NOT set here; each check resolves them through config.R / config.site.R
# and fails loud on its own if an input is missing.
#
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT="${RSCRIPT:-Rscript}"   # overridable (e.g. a specific Rscript, or a stub in tests)

usage() {
  cat >&2 <<'USAGE'
run_stage6_pipeline.sh — run the stage-6 sensitivity & validity checks, with per-check
checkpoint/resume and logging.

  ./run_stage6_pipeline.sh            run/resume: skip completed checks, run the rest
  ./run_stage6_pipeline.sh --status   show done/pending per check, then exit
  ./run_stage6_pipeline.sh --from C   redo from check C (validity|pca)
  ./run_stage6_pipeline.sh --force    ignore all checkpoints, redo everything

Detached (the PCA check loads the full beta matrix), logging to a file, surviving logout:
  setsid nohup ./run_stage6_pipeline.sh > stage6_pipeline.log 2>&1 < /dev/null &
USAGE
}

FORCE=0; FROM=""; STATUS=0
while [ $# -gt 0 ]; do
  case "$1" in
    --force)   FORCE=1 ;;
    --from)    shift; FROM="${1:-}" ;;
    --status)  STATUS=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

LOGDIR="${METHYL_LOGS_DIR:-logs}"   # export METHYL_LOGS_DIR=<work/logs> to keep logs out of the repo
CKPT="$LOGDIR/.ckpt"
mkdir -p "$CKPT"

# ordered checks: "name|Rscript path"
STEPS=(
  "validity|stage6/validity_clocks.R"
  "pca|stage6/pca_sex_batch.R"
)
step_names() { local s; for s in "${STEPS[@]}"; do printf '%s ' "${s%%|*}"; done; }
ts() { date '+%Y-%m-%d %H:%M:%S'; }

if [ "$STATUS" -eq 1 ]; then
  for s in "${STEPS[@]}"; do
    IFS='|' read -r name script <<< "$s"
    [ -f "$CKPT/$name.done" ] && st="done" || st="pending"
    printf "%-10s %-32s %s\n" "$name" "$script" "$st"
  done
  exit 0
fi

# --force clears every checkpoint; --from clears the named check and all later ones.
if [ "$FORCE" -eq 1 ]; then rm -f "$CKPT"/*.done; fi
if [ -n "$FROM" ]; then
  case " $(step_names) " in
    *" $FROM "*) : ;;
    *) echo "unknown --from check: $FROM (valid: $(step_names))" >&2; exit 2 ;;
  esac
  seen=0
  for s in "${STEPS[@]}"; do
    name="${s%%|*}"
    [ "$name" = "$FROM" ] && seen=1
    [ "$seen" -eq 1 ] && rm -f "$CKPT/$name.done"
  done
fi

echo "[$(ts)] stage-6 pipeline start (checks: $(step_names))"
failed=0
for s in "${STEPS[@]}"; do
  IFS='|' read -r name script <<< "$s"
  if [ -f "$CKPT/$name.done" ]; then
    echo "[$(ts)] skip  $name  (checkpoint present)"
    continue
  fi
  echo "[$(ts)] run   $name  ($script)"
  log="$LOGDIR/stage6_pipeline_${name}.log"
  # process substitution (not a pipe) so $? is the Rscript exit code, not tee's
  if "$RSCRIPT" "$script" > >(tee "$log") 2>&1; then
    touch "$CKPT/$name.done"
    echo "[$(ts)] ok    $name"
  else
    rc=$?
    echo "[$(ts)] FAIL  $name  (exit $rc) — see $log. Fix the cause and re-run to resume here." >&2
    failed=$((failed + 1))
  fi
done
if [ "$failed" -gt 0 ]; then
  echo "[$(ts)] stage-6 pipeline finished with $failed failed check(s); re-run to resume them" >&2
  exit 1
fi
echo "[$(ts)] stage-6 pipeline complete"
