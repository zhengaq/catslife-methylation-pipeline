#!/usr/bin/env Rscript
### run_stage5.R: epigenetic clock computation and descriptive statistics / validation
### (population -> reliability -> longitudinal -> report).
source("config.R")
source("stage5/helpers.R")
source("stage5/population.R")
source("stage5/reliability.R")     # technical reliability of the clocks across duplicate aliquots
source("stage5/longitudinal.R")    # within-person stability across waves
source("stage5/report.R")
cat("run_stage5 complete\n")
