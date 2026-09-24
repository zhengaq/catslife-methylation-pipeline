#!/usr/bin/env Rscript
### run_stage6.R: sensitivity and validity checks on the stage 1-5 output. The checks are
### independent; one that fails (e.g. a missing input or package) is reported and the
### others still run.
###   - validity_clocks.R : clock-vs-age validity under non-independence, DNAmTL identity,
###                         ID resolution, clock-NA propagation (needs stage-5 mAge_clocks.csv)
###   - pca_sex_batch.R   : sex-chromosome / batch structure PCA (needs stage-1 betas + minfi)
source("config.R")

run_check <- function(name, path) {
    cat("\n----- stage 6:", name, "-----\n")
    ok <- tryCatch({ source(path, local = new.env()); TRUE },
                   error = function(e) { cat("SKIP/FAIL", name, ":", conditionMessage(e), "\n"); FALSE })
    invisible(ok)
}

run_check("clock validity/sensitivity", "stage6/validity_clocks.R")
run_check("sex-chromosome / batch PCA", "stage6/pca_sex_batch.R")
cat("\nrun_stage6 complete\n")
