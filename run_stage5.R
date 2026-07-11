#!/usr/bin/env Rscript
### run_stage5.R — orchestrate the core stage-5 scope: epigenetic clock computation +
### descriptive stats / validation (population -> reliability -> longitudinal -> report).
### Inferential analyses (LME, twin-corr, meta-analysis, ADCE heritability) are a separate,
### explicitly-invoked step, out of the core QC scope: see run_additional_analysis.R.
source("config.R")
source("stage5/helpers.R")
source("stage5/population.R")
source("stage5/reliability.R")     # duplicate technical reliability of the clocks (DUPS_FILE)
source("stage5/longitudinal.R")    # wave-1 ∩ wave-2 within-person stability (descriptive QC)
source("stage5/report.R")
cat("run_stage5 complete\n")
