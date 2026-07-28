### config.site.example.R — environment / server path profile (TEMPLATE).
###
### Copy this to  config.site.R  (gitignored, never committed) and set the paths. config.R loads config.site.R automatically and reads these
### via the METHYL_* env vars, so NO analysis script needs editing. Only set what
### differs from the in-repo defaults; delete the lines you don't need.
###
### After editing, check the mapping from an R session at the repo root:
###   source("config.R"); describe_paths(); validate_paths("stage1")

Sys.setenv(
  ## --- where the code lives (so the working directory doesn't matter) ---
  METHYL_PROJECT_DIR  = "/home/you/catslife-methylation",

  ## --- roots: inputs (often read-only) vs writable outputs ---
  METHYL_DATA_DIR     = "/secure/catslife/inputs",       # default root for the input files below
  METHYL_ANALYSIS_DIR = "/scratch/you/catslife/work",    # base for all outputs (must be writable)

  ## --- output sub-trees (optional): separate analysis-ready data / regenerable
  ## artifacts / findings / logs. Omit these four to keep everything flat under
  ## METHYL_ANALYSIS_DIR (the historical layout). ---
  METHYL_DERIVED_DIR      = "/scratch/you/catslife/work/derived",       # person table, dyads, phenotype, cell props, mAge_clocks*
  METHYL_INTERMEDIATE_DIR = "/scratch/you/catslife/work/intermediate",  # *.RDat, adjusted betas, rank_corr (regenerable)
  METHYL_RESULTS_DIR      = "/scratch/you/catslife/work/results",       # tables/, sensitivity/, reports/(figures)
  METHYL_LOGS_DIR         = "/scratch/you/catslife/work/logs",          # run logs + orchestrator checkpoints

  ## --- stage 1 (raw arrays) ---
  ## IDATs live one dir per Sentrix barcode under IDAT_DIR (Released_Data/Data/<barcode>/...) —
  ## point IDAT_DIR at that Data/ root. The sample sheet is a SIBLING of Data/, not nested in it.
  METHYL_IDAT_DIR     = "/secure/catslife/Released_Data/Data",
  METHYL_SAMPLE_SHEET = "/secure/catslife/Released_Data/GenomeStudio_Project_and_Files/SampleSheet_S_Reynolds_Smolen_CognitiveAging_2_EPIC.csv",
  ## Pedigree/SIF file (.xlsx: Family/Individual/Father/Mother/Sex/Subject_ID/Population).
  METHYL_ID_KEY       = "/secure/catslife/Released_Data/Sample_Information/S_Reynolds_Smolen_SIF.xlsx",

  ## --- person table (admin .sav + sample-list crosswalk) ---
  ## build_person_table.R merges ADMIN_FILE (person source) with SAMPLE_LIST_FILE
  ## (the random_id<->nidaid crosswalk, .xlsx) on nidaid and produces CLEAN_ID_FILE —
  ## a derived data product, so point it under DERIVED_DIR.
  METHYL_ADMIN_FILE       = "/secure/catslife/individual_admin.sav",
  METHYL_SAMPLE_LIST_FILE = "/secure/catslife/Buffy Coat DNA Methylation Sample List.xlsx",
  METHYL_CLEAN_ID_FILE    = "/scratch/you/catslife/work/derived/catslife_person_table.sav",

  ## --- stage 5 (clocks + additional analysis): point to your filenames ---
  METHYL_BETAS_FILE      = "/secure/catslife/BetasFile.csv",
  METHYL_PHENOTYPE_FILE  = "/scratch/you/catslife/work/derived/PhenotypeFile.csv",  # built by build_phenotype_file.R
  METHYL_HORVATH_CPGS    = "/secure/catslife/HorvathCpGsFile.csv",
  METHYL_TWIN_PHENO_FILE = "/secure/catslife/MethylationAges_PhenotypeInfo_File.sav"
)

## Notes:
##  - CLEAN_ID_FILE is produced by build_person_table.R from ADMIN_FILE — run that first.
##  - DYADS_FILE is produced under DERIVED_DIR; TWIN_RS_FILE under ANALYSIS_DIR — no need to set them.
##  - Sys.setenv here is visible to R only. For run_stage5_pipeline.sh / run_stage6_pipeline.sh
##    to write logs under work/logs, also export it in your shell:
##      export METHYL_LOGS_DIR=/scratch/you/catslife/work/logs
