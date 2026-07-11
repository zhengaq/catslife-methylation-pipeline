### stage5/longitudinal.R — per-clock within-person wave-1 ∩ wave-2 measurement
### stability among LongitudinalGroupIDs (IBD-flagged cross-wave resamples) sampled
### in BOTH waves within a tissue: cross-wave correlation, mean Δ epigenetic vs
### chronological age, and their ratio. Same person years apart, so this is
### biological stability, distinct from reliability.R's same-DNA technical
### replicates. Descriptive only. No both-wave overlap -> an empty, well-formed
### table (not an error). Writes output/tables/clock_wave_overlap.csv.
source("config.R"); source("stage5/helpers.R")

TABLES_DIR <- file.path(ANALYSIS_DIR, "tables"); dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
OUT      <- file.path(TABLES_DIR, "clock_wave_overlap.csv")
OUT_COLS <- c("clock", "tissue", "n_pairs", "n_individuals", "retest_r",
              "mean_delta_epi", "sd_delta_epi", "mean_delta_chrono", "epi_per_chrono_year")

m  <- read.csv(file.path(ANALYSIS_DIR, "mAge_clocks.csv"))
ph <- read.csv(PHENOTYPE_FILE)
clock_cols <- intersect(LME_CLOCKS, names(m))

write_overlap <- function(df) {
  write.csv(df, OUT, row.names = FALSE)
  cat("longitudinal: wrote", basename(OUT), "-", nrow(df), "clock x tissue row(s)\n")
}
empty_tab <- function(msg) {
  cat("longitudinal:", msg, "- no both-wave overlap to describe\n")
  write_overlap(setNames(data.frame(matrix(nrow = 0, ncol = length(OUT_COLS))), OUT_COLS))
}

## population.R drops LongitudinalGroupID/Wave from mAge_clocks.csv, so recover
## them by a Sample join. Both columns are added by the phenotype bridge.
if (!all(c("LongitudinalGroupID", "Wave") %in% names(ph)) || !length(clock_cols)) {
  empty_tab("PHENOTYPE_FILE lacks LongitudinalGroupID/Wave (a single-wave cohort)")
} else {
  d   <- merge(m, unique(ph[, c("Sample", "LongitudinalGroupID", "Wave")]), by = "Sample", all.x = TRUE)
  res <- do.call(rbind, lapply(clock_cols, function(cl) wave_overlap_one_clock(d, cl)))
  if (is.null(res) || !nrow(res)) {
    empty_tab("no LongitudinalGroupID is sampled in both wave 1 and wave 2")
  } else {
    write_overlap(res)
    cat("longitudinal:", sum(res$n_pairs), "cross-wave pair-observation(s) over",
        length(unique(res$tissue)), "tissue(s)\n")
  }
}
