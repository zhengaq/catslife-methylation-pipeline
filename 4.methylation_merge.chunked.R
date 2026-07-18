#!/usr/bin/env Rscript
### Merge the per-chunk batch-adjusted blood + saliva outputs from stage 3 into
### one matrix (rows = CpGs, cols = samples) -> output/B.adjusted.platebatches.txt
source("config.R")
suppressMessages(library(data.table))

nchunks <- NPARTS

d <- NULL
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

    if (is.null(d)) {
        d <- d.tmp
    } else if (all(colnames(d) == colnames(d.tmp))) {
        d <- rbind(d, d.tmp)
    } else {
        cat("Colnames of d & d.tmp DO NOT MATCH.\nchunk:", chunk, "Exiting now.\n", date(), "\n")
        quit(save = 'no')
    }
    cat("Completed chunk", chunk, "\n")
}
cat("Completed merging.\n", date(), "\n")

out <- file.path(ANALYSIS_DIR, "B.adjusted.platebatches.txt")
fwrite(d, out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
cat("Wrote", out, "\n")
