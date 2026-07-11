### stage5/reliability.R — per-clock technical reliability across the intentional-
### duplicate aliquots (DUPS_FILE) that build_phenotype_file.R retains and tags with
### a shared DupGroupID. Same-DNA re-runs, so this measures assay + pipeline
### reliability (not biology). Writes output/tables/clock_duplicate_reliability.csv;
### no duplicates -> an empty, well-formed table (not an error).
source("config.R"); source("stage5/helpers.R")

TABLES_DIR <- file.path(ANALYSIS_DIR, "tables"); dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
OUT      <- file.path(TABLES_DIR, "clock_duplicate_reliability.csv")
REL_COLS <- c("clock", "n_groups", "n_samples", "icc", "retest_r", "mean_abs_diff", "repeatability_coef")

m  <- read.csv(file.path(ANALYSIS_DIR, "mAge_clocks.csv"))
ph <- read.csv(PHENOTYPE_FILE)
clock_cols <- intersect(LME_CLOCKS, names(m))

write_reliability <- function(df) {
  write.csv(df, OUT, row.names = FALSE)
  cat("reliability: wrote", basename(OUT), "-", nrow(df), "clock row(s)\n")
}
empty_tab <- function(msg) {
  cat("reliability:", msg, "- no duplicate reliability to compute\n")
  write_reliability(setNames(data.frame(matrix(nrow = 0, ncol = length(REL_COLS))), REL_COLS))
}

## population.R drops DupGroupID from mAge_clocks.csv, so recover it from the
## phenotype by a Sample join.
if (!"DupGroupID" %in% names(ph)) {
  empty_tab("PHENOTYPE_FILE has no DupGroupID column (a cohort with no intentional-duplicate pairs)")
} else {
  d   <- merge(m, unique(ph[, c("Sample", "DupGroupID")]), by = "Sample", all.x = TRUE)
  d   <- d[!is.na(d$DupGroupID), , drop = FALSE]
  gsz <- table(d$DupGroupID)
  d   <- d[d$DupGroupID %in% names(gsz)[gsz >= 2L], , drop = FALSE]
  if (nrow(d) == 0L || !length(clock_cols)) {
    empty_tab("no DupGroupID group has >= 2 retained samples")
  } else {
    cat("reliability: ", length(unique(d$DupGroupID)), " duplicate group(s), ", nrow(d),
        " samples, ", length(clock_cols), " clocks\n", sep = "")
    res <- do.call(rbind, lapply(clock_cols, function(cl) dup_reliability_one_clock(d, cl)))
    write_reliability(res)
  }
}
