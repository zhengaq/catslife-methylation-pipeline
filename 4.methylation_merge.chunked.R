#!/usr/bin/env Rscript
### Merge the per-chunk batch-adjusted blood + saliva outputs from stage 3 into
### one matrix (rows = CpGs, cols = samples) -> output/B.adjusted.platebatches.txt
### Streamed chunk-by-chunk (fwrite append), so the full cohort matrix is never held in memory.
source("config.R")
suppressMessages(library(data.table))

nchunks  <- NPARTS
out      <- file.path(ANALYSIS_DIR, "B.adjusted.platebatches.txt")
ref_cols <- NULL

for (chunk in 1:nchunks) {
    ds <- fread(file.path(REPORT_DIR, paste0("B.adjusted.regression.saliva.", chunk, ".txt")), header = TRUE)
    db <- fread(file.path(REPORT_DIR, paste0("B.adjusted.regression.blood.",  chunk, ".txt")), header = TRUE)

    ### Defensive: drop any trailing all-NA "V#" column
    for (nm in grep("^V[0-9]+$", colnames(ds), value = TRUE)) if (all(is.na(ds[[nm]]))) ds[, (nm) := NULL]
    for (nm in grep("^V[0-9]+$", colnames(db), value = TRUE)) if (all(is.na(db[[nm]]))) db[, (nm) := NULL]

    ### Merge blood + saliva for this chunk on the shared CpG column
    if (all(ds[[1]] == db[[1]])) {
        d.tmp <- cbind(db, ds[, -1])
    } else {
        d.tmp <- merge(ds, db)
    }

    ### Every chunk must carry the same columns in the same order (the appended rows rely on it)
    if (is.null(ref_cols)) {
        ref_cols <- colnames(d.tmp)
    } else if (!all(colnames(d.tmp) == ref_cols)) {
        cat("Colnames of chunk", chunk, "DO NOT MATCH chunk 1. Exiting now.\n", date(), "\n")
        quit(save = 'no')
    }

    ### Header only on the first chunk (append=FALSE also truncates any stale file); append the rest
    fwrite(d.tmp, out, quote = FALSE, sep = "\t", row.names = FALSE,
           col.names = (chunk == 1), append = (chunk > 1))
    rm(ds, db, d.tmp); gc()
    cat("Completed chunk", chunk, "\n")
}
cat("Wrote", out, "\n", date(), "\n")
