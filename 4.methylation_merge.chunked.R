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
    ### Read whichever per-tissue adjusted files stage 3 wrote (a single-tissue wave writes one),
    ### and merge the present tissues on the shared CpG column.
    tissue_files <- c(blood  = file.path(REPORT_DIR, paste0("B.adjusted.regression.blood.",  chunk, ".txt")),
                      saliva = file.path(REPORT_DIR, paste0("B.adjusted.regression.saliva.", chunk, ".txt")))
    tissue_files <- tissue_files[file.exists(tissue_files)]
    if (!length(tissue_files))
        stop("stage 4: no per-tissue adjusted files for chunk ", chunk, " under ", REPORT_DIR)

    parts <- lapply(tissue_files, function(f) {
        d <- fread(f, header = TRUE)
        ### Defensive: drop any trailing all-NA "V#" column
        for (nm in grep("^V[0-9]+$", colnames(d), value = TRUE)) if (all(is.na(d[[nm]]))) d[, (nm) := NULL]
        d
    })
    ### One tissue -> passthrough; two -> cbind on matching CpG order, else merge
    d.tmp <- Reduce(function(x, y) if (all(x[[1]] == y[[1]])) cbind(x, y[, -1]) else merge(x, y), parts)

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
    rm(parts, d.tmp); gc()
    cat("Completed chunk", chunk, "\n")
}
cat("Wrote", out, "\n", date(), "\n")
