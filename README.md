# CATSLife Methylation Pipeline

Sequential Illumina EPIC preprocessing and epigenetic-age clock computation for
CATSLife DNA methylation data. Tissues: blood (PBMC + Buffy Coat) and Saliva,
plus Cell Line QC controls. Supports both **EPIC v1** (legacy) and **EPIC
v2.0** arrays via `ARRAY_VERSION` in `config.R`.

## Pipeline stages

| # | Script | Does |
|---|--------|------|
| 1 | `1.methylation_minfi.R` | Read IDATs (minfi) -> detection-p QC (sample + probe missingness) -> `mapToGenome` + drop SNP loci -> `noob` background correction -> `dasen` normalization (wateRmelon) -> betas / M-values + QC PDFs |
| 2 | `2.methylation_pca.R` | PCA on M-values (optionally a random CpG subset), colored by tissue |
| 3 | `3.methylation_adjust.chunked.R` | Estimate cell-type proportions (EpiDISH for blood, BeadSorted/`estimateLC` for saliva) -> residualize betas on cell proportions -> residualize on plate batch. Chunked via `--part`/`--nparts` for parallel (e.g. SLURM array) runs |
| 4 | `4.methylation_merge.chunked.R` | Merge the per-chunk blood + saliva outputs -> `B.adjusted.platebatches.txt` |
| 5 | `run_stage5.R` (sources `stage5/`) | ~17 epigenetic clocks via `dnaMethyAge` (a primary pass on stage 1's unadjusted betas, plus a second comparison pass on stage 4's adjusted betas), then descriptive stats and clock-vs-age validation |

Stage 5 needs a person-level phenotype file (age, sex, family). Build it by
running `scripts/build/build_person_table.R` (merges the individual-admin `.sav` with
the sample list into the person table `CLEAN_ID_FILE`), then
`scripts/build/catslife_id_dyads.R` and `scripts/build/build_phenotype_file.R`, before
`run_stage5.R` (see below).

## The ID bridge (array IDs <-> person IDs)

Stages 1-4 key on the array id, which is a de-identified `random_id`, whereas stage 5 keys on person/family ids (i.e., `aid`/`pfamid`).
`scripts/build/build_person_table.R` merges the admin file with the sample list, so the person table carries
`random_id`; `scripts/build/build_phenotype_file.R` then joins the sheet on `random_id`.
The `_<wave>` suffix gives the `Wave`: one row **per sample**, with `Age` from the
admin `LabAge` (wave 2) or `LabAge1` (wave 1) per row.

- keeps intentional duplicate pairs (`DUPS_FILE`) as a grouped consistency check
  rather than dropping them,
- excludes known-problem samples (`PROBLEM_HISTORY_FILE`),
- flags (not drops) likely cross-wave longitudinal resamples (`IBD_FILE`), and
- writes a three-source sex QC report (`output/reports/sex_qc.csv`).

## Environment

1. Confirm R version (expect 4.5.x):
   ```bash
   R --version
   ```
2. Install the CRAN + Bioconductor packages this QC + clock-measurement pipeline
   needs, pinned to the Bioconductor 3.22 release so versions match what it was
   built and tested against:
   ```r
   install.packages("BiocManager")
   BiocManager::install(version = "3.22", update = FALSE, ask = FALSE)

   BiocManager::install(c(
     "remotes", "optparse", "data.table", "tidyverse",
     "corrplot", "haven", "readxl",
     "minfi", "wateRmelon", "EpiDISH", "ExperimentHub",
     "FlowSorted.Blood.EPIC", "BeadSorted.Saliva.EPIC",
     "IlluminaHumanMethylationEPICmanifest",
     "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
     "IlluminaHumanMethylationEPICv2manifest",
     "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
     "impute", "RPMM"
   ), update = FALSE, ask = FALSE)
   ```
3. Install the two GitHub-only packages, pinned to the same commits this
   pipeline was tested against:
   ```r
   remotes::install_github("hhhh5/ewastools@bb1a67dfa737c9bf3166bd366df2979c6228c724", upgrade = "never")
   remotes::install_github("yiluyucheng/dnaMethyAge@0d40c9bb02c7d9a4d91ce01d150a9c0b12edf1f4", upgrade = "never")
   ```
4. Sanity-check the install:
   ```r
   for (p in c("minfi", "wateRmelon", "EpiDISH", "FlowSorted.Blood.EPIC",
               "BeadSorted.Saliva.EPIC", "ewastools", "dnaMethyAge")) {
     suppressMessages(library(p, character.only = TRUE)); cat("OK ", p, "\n")
   }
   ```


### Run

1. On a new delivery wave, run `Rscript scripts/build/inspect_2026_delivery.R` first.
   It prints the sample sheet / admin file / pedigree / QC-file structure
   (column names, row counts) and a crosswalk dry-run match rate, without
   loading any participant data wholesale, since exact filenames vary between
   waves.
2. `cp config.site.example.R config.site.R` and edit `config.site.R`.
   Set input root(s), the writable output root (`METHYL_ANALYSIS_DIR`), the 
   IDAT directory + sample sheet, the pedigree file, the individual-admin 
   file (`METHYL_ADMIN_FILE`), and the sample list (`METHYL_SAMPLE_LIST_FILE`, 
   the `random_id`<->`nidaid` crosswalk). 
3. Confirm the mapping resolves and inputs/outputs exist:
   ```r
   source("config.R"); describe_paths()
   validate_paths("stage1"); validate_paths("person_table")
   validate_paths("phenotype_bridge"); validate_paths("stage5")
   ```
4. Run in order:
   ```bash
   Rscript 1.methylation_minfi.R
   Rscript 2.methylation_pca.R
   Rscript 3.methylation_adjust.chunked.R --nparts 5 --part 1   # repeat for parts 1..5
   Rscript 4.methylation_merge.chunked.R

   # ID bridge: person table -> dyads -> phenotype
   Rscript scripts/build/build_person_table.R       # ADMIN_FILE + SAMPLE_LIST_FILE -> CLEAN_ID_FILE
   Rscript scripts/build/catslife_id_dyads.R
   Rscript scripts/build/build_phenotype_file.R
   Rscript run_stage5.R
   ```
   `build_phenotype_file.R` fails loud (rather than silently guessing) if any
   non-control sample's `random_id` does not resolve to a person.
5. Set `METHYL_ARRAY_VERSION=v1` if the raw IDATs are legacy EPIC v1 arrays
   (the default is `v2`).

## Optional tools

- `scripts/build/inspect_2026_delivery.R` - read-only recon + crosswalk dry-run for a
  new delivery wave (run this first, see step 1 above).
- `scripts/build/crosscheck_genomestudio_betas.R` - compares our computed dasen betas
  against GenomeStudio's own `Methylation_Profile.txt` for a probe/sample
  spot-check. An independent sanity check, not part of the normal run.

## Additional notes

- Existing sex discrepancy. 4 participants with identified mis-match on stated vs. genotyped; Any additional mismatch (i.e., admin file's self-report vs. pedigree file) will be flagged in 
  `output/reports/sex_qc.csv`.
- Age is the admin file's `LabAge` (wave 2) or `LabAge1` (wave 1), assigned by the
  sample's `Wave` (the sheet's `_2` suffix). Samples with no age for their wave
  resolve normally but have no age-acceleration value.
- The array->person crosswalk is `random_id` -> `nidaid` (via the sample list) -> the admin
  person
