### config.site.example.R: template for the site path profile.
###
### Copy it to config.site.R (gitignored) and set the paths. config.R loads config.site.R
### automatically and reads these METHYL_* variables, so no analysis script needs editing.
### Set only what differs from the defaults in config.R and delete the other lines.
###
### Then check the mapping from an R session at the repo root:
###   source("config.R"); describe_paths(); validate_paths("stage1")

Sys.setenv(
  ## Where the code lives (so the working directory doesn't matter)
  METHYL_PROJECT_DIR  = "/home/you/catslife-methylation",

  ## Roots: inputs (often read-only) and writable outputs
  METHYL_DATA_DIR     = "/secure/catslife/inputs",       # default root for the input files below
  METHYL_ANALYSIS_DIR = "/scratch/you/catslife/work",    # base for all outputs (must be writable)

  ## Optional output sub-trees for analysis-ready data, regenerable artifacts, findings and
  ## logs. Omit them to keep all outputs in one directory under METHYL_ANALYSIS_DIR.
  METHYL_DERIVED_DIR      = "/scratch/you/catslife/work/derived",       # person table, dyads, phenotype, cell props, mAge_clocks*
  METHYL_INTERMEDIATE_DIR = "/scratch/you/catslife/work/intermediate",  # *.RDat, adjusted betas, rank_corr (regenerable)
  METHYL_RESULTS_DIR      = "/scratch/you/catslife/work/results",       # tables/, sensitivity/, reports/(figures)
  METHYL_LOGS_DIR         = "/scratch/you/catslife/work/logs",          # run logs + orchestrator checkpoints

  ## Stage 1 (raw arrays). IDATs live one directory per Sentrix barcode
  ## (Released_Data/Data/<barcode>/...); point IDAT_DIR at that Data/ directory. The sample
  ## sheet sits next to Data/, not inside it.
  METHYL_IDAT_DIR     = "/secure/catslife/Released_Data/Data",
  METHYL_SAMPLE_SHEET = "/secure/catslife/Released_Data/GenomeStudio_Project_and_Files/SampleSheet_S_Reynolds_Smolen_CognitiveAging_2_EPIC.csv",
  ## Pedigree/SIF file (.xlsx: Family/Individual/Father/Mother/Sex/Subject_ID/Population).
  METHYL_ID_KEY       = "/secure/catslife/Released_Data/Sample_Information/S_Reynolds_Smolen_SIF.xlsx",

  ## Person table. build_person_table.R merges ADMIN_FILE with SAMPLE_LIST_FILE (the
  ## random_id <-> nidaid crosswalk, .xlsx) on nidaid and writes CLEAN_ID_FILE, a derived
  ## file that belongs under DERIVED_DIR.
  METHYL_ADMIN_FILE       = "/secure/catslife/individual_admin.sav",
  METHYL_SAMPLE_LIST_FILE = "/secure/catslife/Buffy Coat DNA Methylation Sample List.xlsx",
  METHYL_CLEAN_ID_FILE    = "/scratch/you/catslife/work/derived/catslife_person_table.sav",

  ## Sample-label corrections (see config.R, "Sample-identity corrections")
  METHYL_SAMPLE_SWAPS_FILE = "/secure/catslife/C12_Methylation_Sample_Swaps.csv",

  ## Stage 5: the phenotype file written by build_phenotype_file.R
  METHYL_PHENOTYPE_FILE  = "/scratch/you/catslife/work/derived/PhenotypeFile.csv"
)

## Notes:
##  - DYADS_FILE is written under DERIVED_DIR; no need to set it.
##  - Use absolute paths: a relative path resolves against whichever directory R starts in.
##  - Every name needs the METHYL_ prefix; config.R stops on a bare name or a <placeholder>.
