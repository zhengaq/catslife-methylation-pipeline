#!/usr/bin/env bash
#
# reorganize_stage2.sh — refinements applied AFTER reorganize.sh (stage 1) has run and the
# tree is already split into ../data/ + ../work/{derived,intermediate,results,logs,_review}.
#
# It assumes stage 1 left things here (from reorganize.sh):
#   work/results/mAge_clocks.csv, work/results/mAge_clocks_adjusted.csv
#   work/results/reports/          (current figures: PCA_sex_batch.pdf, sex_qc.csv, ...)
#   work/_review/CATSLife_methylation_clocks_QC_report.md
#   work/_review/output/           (older run: violin_*.png, PCA_plots.pdf, sample_missingness.txt,
#                                    plus stale duplicates of the current figures)
#
# Refinements (your three points on the aftermath):
#   1. clock matrices are analysis-ready data     -> work/derived/
#   2. the QC report is a deliverable             -> work/results/
#   3. the older run's unique figures are report  -> work/results/reports/
#
# Every operation is a MOVE (reversible). Nothing is hard-deleted. DRY-RUN by default;
# pass --apply to execute. Absent sources are skipped and reported, so if your tree differs
# from the assumptions above a dry run shows exactly what would (not) move — check it first.
#
set -uo pipefail

APPLY=0
BASE_OVERRIDE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1 ;;
    --base)  shift; BASE_OVERRIDE="${1:-}" ;;
    -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

REPO="$(cd "$(dirname "$0")" && pwd)"
[ -f "$REPO/config.R" ] || { echo "refusing: $REPO has no config.R (run from the repo root)" >&2; exit 1; }
BASE="${BASE_OVERRIDE:-$(dirname "$REPO")}"
WORK="$BASE/work"
[ -d "$WORK" ] || { echo "refusing: $WORK not found (run reorganize.sh --apply first, or set --base)" >&2; exit 1; }

DERIVED="$WORK/derived"
RESULTS="$WORK/results"
RESREPORTS="$RESULTS/reports"
REVIEW="$WORK/_review"

echo "base:    $BASE"
[ "$APPLY" -eq 1 ] && echo "mode:    APPLY (moving files)" || echo "mode:    DRY RUN (no changes; pass --apply to execute)"
echo

moved=0; skipped=0

# mv2 <src-abs> <dest-dir>
mv2() {
  local src="$1" dest="$2"
  if [ ! -e "$src" ]; then printf '  skip (absent)  %s\n' "${src#$BASE/}"; skipped=$((skipped+1)); return 0; fi
  if [ "$APPLY" -eq 1 ]; then
    mkdir -p "$dest"; mv "$src" "$dest/"
    printf '  moved  %-42s -> %s/\n' "${src#$WORK/}" "${dest#$BASE/}"
  else
    printf '  DRY    %-42s -> %s/\n' "${src#$WORK/}" "${dest#$BASE/}"
  fi
  moved=$((moved+1))
}
section() { echo; echo "== $1 =="; }

# ---- 1. clock matrices -> work/derived/ ---------------------------------
section "clock matrices: results/ -> derived/"
mv2 "$RESULTS/mAge_clocks.csv"          "$DERIVED"
mv2 "$RESULTS/mAge_clocks_adjusted.csv" "$DERIVED"

# ---- 2. QC report -> work/results/ --------------------------------------
section "QC report: _review/ -> results/"
mv2 "$REVIEW/CATSLife_methylation_clocks_QC_report.md" "$RESULTS"

# ---- 3. older run's unique figures -> work/results/reports/ -------------
# The violins, the stage-2 PCA_plots.pdf, and the per-sample missingness table live only in
# the older output/ run; move them in with the current figures. The three figures that also
# exist (current) in results/reports/ are left behind in _review/output/ as duplicates.
section "older figures: _review/output/ -> results/reports/"
if [ -d "$REVIEW/output" ]; then
  for f in "$REVIEW"/output/violin_*.png; do [ -e "$f" ] && mv2 "$f" "$RESREPORTS"; done
  mv2 "$REVIEW/output/PCA_plots.pdf"        "$RESREPORTS"
  mv2 "$REVIEW/output/sample_missingness.txt" "$RESREPORTS"
else
  echo "  skip (absent)  _review/output/"
fi

# ---- summary + next steps ------------------------------------------------
echo
echo "== summary =="
echo "  $moved item(s) $([ "$APPLY" -eq 1 ] && echo moved || echo "would move"); $skipped absent/skipped."
if [ "$APPLY" -eq 1 ]; then
  cat <<EOF

Done (stage 2).
  work/derived/   now also holds mAge_clocks{,_adjusted}.csv
  work/results/   now holds the QC report
  work/results/reports/  now holds the violins + PCA_plots.pdf + sample_missingness.txt

Left for you:
  - work/_review/output/ retains only the 3 figures that duplicate results/reports/
    (Sample_missingness_tissue.jpg, minfi_QC.missingness_detP.pdf, sex_qc.csv). Confirm the
    results/reports/ copies are the ones you want, then delete _review/output/ and _review/.
  - The violins predate the final sex-exclusion clock re-run; regenerate them with stage 5's
    report step if you need post-exclusion versions.
EOF
else
  echo "  re-run with --apply to perform the moves."
fi
