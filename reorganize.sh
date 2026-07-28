#!/usr/bin/env bash
#
# reorganize.sh — restructure the CATSLife methylation working directory so the git
# checkout holds only code, and all data lives in sibling data/ (raw inputs) and
# work/ (derived / intermediate / results / logs) trees.
#
#   REPO/            the git checkout this script sits in (code only, after this runs)
#   ../data/         raw inputs (read-only)
#   ../work/derived/       built ID bridge: person table, dyads, phenotype, cell proportions
#   ../work/intermediate/  stage 1-4 regenerable artifacts (*.RDat, B.*, testMpca, rank_corr)
#   ../work/results/       deliverables: mAge_clocks*, tables/, sensitivity/, reports/(figures)
#   ../work/logs/          run logs + orchestrator checkpoints
#   ../work/_review/       stale / duplicate / session cruft — inspect, then delete yourself
#
# Every operation is a MOVE (reversible on the same filesystem) or a quarantine into
# _review/. Nothing is hard-deleted. DRY-RUN by default; pass --apply to execute.
#
# STAGE 1 of the reorganization: it relocates the original flat repo-root layout into
# data/ + work/. The follow-up refinements (clocks -> derived/, report -> results/, the
# older figures -> results/reports/) are in reorganize_stage2.sh, which runs AFTER this.
#
set -uo pipefail

# ---- args ----------------------------------------------------------------
APPLY=0
BASE_OVERRIDE=""
while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1 ;;
    --base)  shift; BASE_OVERRIDE="${1:-}" ;;
    -h|--help)
      sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

# ---- locate repo + targets ----------------------------------------------
REPO="$(cd "$(dirname "$0")" && pwd)"
[ -f "$REPO/config.R" ] || { echo "refusing: $REPO has no config.R (run from the repo root)" >&2; exit 1; }
BASE="${BASE_OVERRIDE:-$(dirname "$REPO")}"
[ "$BASE" != "$REPO" ] || { echo "refusing: target base must differ from the repo" >&2; exit 1; }

DATA="$BASE/data"
DERIVED="$BASE/work/derived"
INTER="$BASE/work/intermediate"
RESULTS="$BASE/work/results"
RESREPORTS="$RESULTS/reports"
LOGS="$BASE/work/logs"
REVIEW="$BASE/work/_review"

echo "repo:    $REPO"
echo "base:    $BASE"
[ "$APPLY" -eq 1 ] && echo "mode:    APPLY (moving files)" || echo "mode:    DRY RUN (no changes; pass --apply to execute)"
echo

moved=0; skipped=0

# move <destdir> <relpath...> — move each existing REPO/<relpath> into <destdir>
move() {
  local dest="$1"; shift
  local rel src
  for rel in "$@"; do
    src="$REPO/$rel"
    if [ ! -e "$src" ]; then skipped=$((skipped+1)); continue; fi
    if [ "$APPLY" -eq 1 ]; then
      mkdir -p "$dest"; mv "$src" "$dest/"
      printf '  moved  %-46s -> %s/\n' "$rel" "${dest#$BASE/}"
    else
      printf '  DRY    %-46s -> %s/\n' "$rel" "${dest#$BASE/}"
    fi
    moved=$((moved+1))
  done
}

# move_contents <destdir> <reldir> — move everything inside REPO/<reldir> into <destdir>,
# then drop the emptied dir (so <reldir> is not re-nested one level down under <destdir>).
move_contents() {
  local dest="$1" rel="$2" src="$REPO/$2"
  if [ ! -d "$src" ]; then skipped=$((skipped+1)); return 0; fi
  if [ "$APPLY" -eq 1 ]; then
    mkdir -p "$dest"
    find "$src" -mindepth 1 -maxdepth 1 -exec mv {} "$dest/" \;
    rmdir "$src" 2>/dev/null || true
    printf '  merged %-46s -> %s/\n' "$rel/*" "${dest#$BASE/}"
  else
    printf '  DRY    %-46s -> %s/\n' "$rel/* (contents)" "${dest#$BASE/}"
  fi
  moved=$((moved+1))
}

section() { echo; echo "== $1 =="; }

# ---- RAW inputs -> data/ -------------------------------------------------
section "raw inputs -> data/"
move "$DATA" \
  "cats12_admin_withLabAge1_10Feb2026.sav" \
  "Buffy Coat DNA Methylation Sample List.xlsx"

# ---- DERIVED (built ID bridge) -> work/derived/ --------------------------
section "derived id-bridge -> work/derived/"
move "$DERIVED" \
  "catslife_person_table.sav" \
  "catslife_dyads.csv" \
  "PhenotypeFile.csv" \
  "reports/cell_proportions.blood.saliva.txt"

# ---- INTERMEDIATE (regenerable stage 1-4) -> work/intermediate/ ----------
# rgSetflt is the retired F_RGFLT checkpoint -> quarantine it BEFORE the RDat sweep.
section "stale checkpoint -> work/_review/"
move "$REVIEW" "methylation_data_rgSetflt.RDat"

section "intermediates -> work/intermediate/"
move "$INTER" \
  "dasen_betas.RDat" \
  "methylation_data_raw.RDat" \
  "methylation_data_detP.RDat" \
  "methylation_data_noob.RDat" \
  "methylation_data_noobflt.RDat" \
  "methylation_data_dasen.RDat" \
  "B.adjusted.platebatches.txt" \
  "testMpca.RDat" \
  "rank_corr.rds" \
  "methylation_data_detP.failedprobe.txt" \
  "methylation_data_detP.failedsamp.txt" \
  "reports/dasen_Mpca_pca.RDat" \
  "reports/B.residualized.blood.1.txt" \
  "reports/B.residualized.blood.2.txt" \
  "reports/B.residualized.blood.3.txt" \
  "reports/B.residualized.blood.4.txt" \
  "reports/B.residualized.blood.5.txt" \
  "reports/B.adjusted.regression.blood.1.txt" \
  "reports/B.adjusted.regression.blood.2.txt" \
  "reports/B.adjusted.regression.blood.3.txt" \
  "reports/B.adjusted.regression.blood.4.txt" \
  "reports/B.adjusted.regression.blood.5.txt"

# ---- RESULTS -> work/results/ -------------------------------------------
section "results -> work/results/"
move "$RESULTS" \
  "mAge_clocks.csv" \
  "mAge_clocks_adjusted.csv" \
  "tables" \
  "sensitivity"
# report figures (leave reports/CATSLife_..._QC_report.md tracked in the repo)
move "$RESREPORTS" \
  "reports/PCA_sex_batch.pdf" \
  "reports/pca_sex_batch_summary.csv" \
  "reports/Sample_missingness_tissue.jpg" \
  "reports/minfi_QC.missingness_detP.pdf" \
  "reports/sex_qc.csv"

# ---- LOGS + checkpoints -> work/logs/ -----------------------------------
section "logs + checkpoints -> work/logs/"
move "$LOGS" \
  "stage1.log" "stage2.log" "stage3.log" "stage4.log" "stage5.log" \
  "stage5_pipeline.log"
move_contents "$LOGS" "logs"   # flatten logs/* (incl .ckpt) into work/logs/, not work/logs/logs/

# ---- QUARANTINE: stale duplicate run + session cruft -> work/_review/ -----
# output/ is an older run whose figures duplicate reports/ (and holds the only violins);
# quarantined rather than merged, so you can salvage the violins or regenerate them.
# The root-level QC report is a loose duplicate of the git-tracked reports/ copy (which
# arrives/updates via `git pull`); quarantine the loose one so the tracked copy is canonical.
section "stale duplicates + session cruft -> work/_review/ (delete after inspection)"
move "$REVIEW" \
  "output" \
  "CATSLife_methylation_clocks_QC_report.md" \
  "Rplots.pdf" \
  ".Rhistory" \
  ".Rproj.user" \
  "output.txt" \
  "directory-tree.txt"

# ---- summary + next steps ------------------------------------------------
echo
echo "== summary =="
echo "  $moved item(s) $([ "$APPLY" -eq 1 ] && echo moved || echo "would move"); $skipped absent/skipped."
if [ "$APPLY" -eq 1 ]; then
  cat <<EOF

Done (stage 1). The repo now holds code only. Next, run reorganize_stage2.sh to apply the
refinements: move the clock matrices to work/derived/, the QC report to work/results/, and
the older figures out of _review into work/results/reports/.
Your finished results are readable now under: $RESULTS
EOF
else
  echo "  re-run with --apply to perform the moves."
fi
