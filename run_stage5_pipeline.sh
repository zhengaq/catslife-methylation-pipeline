#!/usr/bin/env bash
#
# run_stage5_pipeline.sh — orchestrate the full stage-5 chain end to end:
#   1. person_table  scripts/build/build_person_table.R   (admin .sav + sample list -> CLEAN_ID_FILE)
#   2. dyads         scripts/build/catslife_id_dyads.R     (person table -> DYADS_FILE)
#   3. phenotype     scripts/build/build_phenotype_file.R  (the ID bridge -> PhenotypeFile.csv)
#   4. stage5        run_stage5.R                          (clocks + descriptive/validation)
#
# Checkpoint/resume: each step that finishes writes logs/.ckpt/<step>.done. A re-run skips any
# step whose marker exists, so after a failure you just launch again and it resumes at the failed
# step. Every step is fail-fast (a non-zero exit stops the pipeline) and is streamed to both the
# console and logs/stage5_pipeline_<step>.log.
#
# Paths (inputs/outputs) are NOT set here; every step resolves them through config.R / config.site.R
# and fails loud on its own if an input is missing.
#
set -uo pipefail
cd "$(dirname "$0")"
RSCRIPT="${RSCRIPT:-Rscript}"   # overridable (e.g. a specific Rscript, or a stub in tests)

usage() {
  cat >&2 <<'USAGE'
run_stage5_pipeline.sh — build the ID bridge (person table -> dyads -> phenotype) and run stage 5,
with per-step checkpoint/resume and logging.

  ./run_stage5_pipeline.sh            run/resume: skip completed steps, run the rest
  ./run_stage5_pipeline.sh --status   show done/pending per step, then exit
  ./run_stage5_pipeline.sh --from S   redo from step S (person_table|dyads|phenotype|stage5)
  ./run_stage5_pipeline.sh --force    ignore all checkpoints, redo everything

Detached (stage 5 is long), logging to a file, surviving logout:
  setsid nohup ./run_stage5_pipeline.sh > stage5_pipeline.log 2>&1 < /dev/null &
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

LOGDIR=logs
CKPT="$LOGDIR/.ckpt"
mkdir -p "$CKPT"

# ordered steps: "name|Rscript path"
STEPS=(
  "person_table|scripts/build/build_person_table.R"
  "dyads|scripts/build/catslife_id_dyads.R"
  "phenotype|scripts/build/build_phenotype_file.R"
  "stage5|run_stage5.R"
)
step_names() { local s; for s in "${STEPS[@]}"; do printf '%s ' "${s%%|*}"; done; }
ts() { date '+%Y-%m-%d %H:%M:%S'; }

if [ "$STATUS" -eq 1 ]; then
  for s in "${STEPS[@]}"; do
    IFS='|' read -r name script <<< "$s"
    [ -f "$CKPT/$name.done" ] && st="done" || st="pending"
    printf "%-14s %-42s %s\n" "$name" "$script" "$st"
  done
  exit 0
fi

# --force clears every checkpoint; --from clears the named step and all later ones.
if [ "$FORCE" -eq 1 ]; then rm -f "$CKPT"/*.done; fi
if [ -n "$FROM" ]; then
  case " $(step_names) " in
    *" $FROM "*) : ;;
    *) echo "unknown --from step: $FROM (valid: $(step_names))" >&2; exit 2 ;;
  esac
  seen=0
  for s in "${STEPS[@]}"; do
    name="${s%%|*}"
    [ "$name" = "$FROM" ] && seen=1
    [ "$seen" -eq 1 ] && rm -f "$CKPT/$name.done"
  done
fi

echo "[$(ts)] stage-5 pipeline start (steps: $(step_names))"
for s in "${STEPS[@]}"; do
  IFS='|' read -r name script <<< "$s"
  if [ -f "$CKPT/$name.done" ]; then
    echo "[$(ts)] skip  $name  (checkpoint present)"
    continue
  fi
  echo "[$(ts)] run   $name  ($script)"
  log="$LOGDIR/stage5_pipeline_${name}.log"
  # process substitution (not a pipe) so $? is the Rscript exit code, not tee's
  if "$RSCRIPT" "$script" > >(tee "$log") 2>&1; then
    touch "$CKPT/$name.done"
    echo "[$(ts)] ok    $name"
  else
    rc=$?
    echo "[$(ts)] FAIL  $name  (exit $rc) — see $log. Fix the cause and re-run to resume here." >&2
    exit "$rc"
  fi
done
echo "[$(ts)] stage-5 pipeline complete"
