### stage5/population.R — clock computation + cross-tissue rank-correlation
### bootstrap. Writes mAge_clocks.csv, rank_corr.rds, and (when
### B.adjusted.platebatches.txt exists) mAge_clocks_adjusted.csv.
source("config.R"); source("stage5/helpers.R")
suppressMessages({
  library(data.table); library(dplyr); library(tibble); library(dnaMethyAge)
})

## X-prefix beta colnames explicitly: as.data.frame() does NOT apply R's
## check.names "X"-mangling to matrix colnames (unlike data.frame()), and the
## clock join below expects the "X"-prefixed sample ids.
x_prefix_cols <- function(m) { colnames(m) <- paste0("X", colnames(m)); m }

## Fraction of Horvath2013's CpGs present in the beta matrix — a broadly-covered proxy
## for whether the full clock loop will find its probes (catches a wrong array version,
## missing v2 probe-id canonicalization, or wrong orientation). Prints the rate and
## returns TRUE/FALSE; the caller decides whether a miss is fatal (Pass 1) or a skip
## (Pass 2, the optional additive pass).
assert_clock_cpg_coverage <- function(Betas, threshold = CLOCK_CPG_COVERAGE_MIN) {
  data("HorvathS2013", package = "dnaMethyAge", envir = environment())
  probes <- setdiff(coefs$Probe, "Intercept")
  rate <- mean(probes %in% rownames(Betas))
  cat("Horvath2013 clock CpG match-rate against this beta matrix:", round(rate, 3), "\n")
  invisible(rate >= threshold)
}

## Shared clock-computation loop (used for both the unadjusted and adjusted passes).
clock_specs <- tibble::tribble(
  ~key,           ~mage,                        ~accel,
  "HannumG2013",  "Hannum_mAge",                "Hannum_Age_Acceleration",
  "HorvathS2013", "Horvath_mAge",               "Horvath_Age_Acceleration",
  "HorvathS2018", "Horvath2_mAge",              "Horvath2_Age_Acceleration",
  "ZhangQ2019",   "ZhangQ_mAge",                "ZhangQ_Age_Acceleration",
  "DunedinPACE",  "Dunedin_Pace",               "Dunedin_Pace_Acceleration",
  "LevineM2018",  "PhenoAge_mAge",              "PhenoAge_Acceleration",
  "YangZ2016",    "epiTOC_mitoticdivisions",    "epiTOC_Acceleration",
  "epiTOC2",      "epiTOC2_mitoticdivisions",   "epiTOC2_Acceleration",
  "PCGrimAge",    "PCGrimAge_mAge",             "PCGrimAge_Acceleration",
  "ZhangY2017",   "ZhangY_mAge",                "ZhangY_Acceleration",
  "LuA2019",      "LuA_mAge",                   "LuA_Acceleration",
  "ShirebyG2020", "ShirebyG2020_mAge",          "ShirebyG2020_Acceleration",
  "McEwenL2019",  "PedBE_mAge",                 "PedBE_Acceleration",
  "LuA2023p2",    "PanM2_mAge",                 "PanM2_Acceleration",
  "LuA2023p3",    "PanM3_mAge",                 "PanM3_Acceleration"
)
compute_clocks <- function(Betas, info2) {
  ## methyAge() can crash (rowMeans on a degenerate/empty subset) when a clock's
  ## probes are ~100% missing from Betas; NA-fill that one clock rather than abort
  ## the whole pass.
  clock_tab <- function(key, mage, accel) {
    a <- tryCatch(as.data.frame(methyAge(Betas, age_info = info2, clock = key)),
                  error = function(e) {
                    warning("compute_clocks: ", key, " crashed (likely ~0% probe coverage) - ",
                            "filling ", mage, "/", accel, " with NA: ", conditionMessage(e))
                    data.frame(Sample = info2$Sample, mAge = NA_real_, Age_Acceleration = NA_real_)
                  })
    names(a)[names(a) == "mAge"] <- mage
    names(a)[names(a) == "Age_Acceleration"] <- accel
    a[, c("Sample", mage, accel)]
  }
  merged <- info2
  for (i in seq_len(nrow(clock_specs))) {
    s <- clock_specs[i, ]
    merged <- left_join(merged, clock_tab(s$key, s$mage, s$accel), by = "Sample")
  }
  invisible(tryCatch(methyAge(Betas, age_info = info2, clock = "BernabeuE2023c"),
                     error = function(e) NULL))   # cAge: author note "does not run"; unused
  ## Undo this function's own "X" prefix — unconditional, unlike strip_x_prefix()
  ## (which only strips a digit-following X), since it removes exactly what was added.
  merged$Sample <- sub("^X", "", merged$Sample)
  merged
}

## ------------------------------------------------------- Phenotype table ----
Demographics <- read.csv(PHENOTYPE_FILE)
if ("DNASource" %in% names(Demographics) && !"DNA_Source" %in% names(Demographics))
  Demographics$DNA_Source <- Demographics$DNASource
## Canonicalize defensively and drop Cell Line QC controls — the clocks are trained
## on human tissue, not cell lines.
Demographics$DNA_Source <- canonicalize_dna_source(Demographics$DNA_Source)
Demographics <- Demographics[Demographics$DNA_Source != "Cell_Line", ]
if (!"IndividualID" %in% names(Demographics))
  stop("PHENOTYPE_FILE missing IndividualID — build it with scripts/build/build_phenotype_file.R.")
info <- Demographics[, intersect(c("Sample", "IndividualID", "Age", "Sex", "FamilyID", "DNA_Source"),
                                 names(Demographics))]
info$Sample <- paste0("X", info$Sample); info2 <- info

## ------------------------------------------------------ Pass 1: unadjusted ----
## Published clocks are fixed-weight predictors trained on normalized input, so
## they run on stage-1's dasen betas, not the cell/plate-adjusted stage-4 betas.
dv    <- load_one(F_DASENB)
Betas <- dv$b; rm(dv); gc()                    ### the M-values in dv are unused here; free them before the beta copy
Betas <- x_prefix_cols(as.data.frame(Betas))
if (!assert_clock_cpg_coverage(Betas))
  stop("population.R: too few clock CpGs in the beta matrix — check ARRAY_VERSION / ",
       "probe-id canonicalization before trusting clock output.")
merged <- compute_clocks(Betas, info2)
rm(Betas); gc()   ### Pass 1 betas consumed; free before the rank bootstrap and the Pass 2 adjusted betas
write.csv(merged, file.path(ANALYSIS_DIR, "mAge_clocks.csv"), row.names = FALSE)
cat("clocks: wrote mAge_clocks.csv -", ncol(merged), "cols,", nrow(merged), "samples\n")

## ------------------------------------------ Rank-correlation bootstrap ----
## Cross-tissue consistency of each clock: bootstrap the Spearman rank correlation
## across tissues (+ age in PBMC). Needs >=2 tissues present — check upfront rather
## than run 100 iterations x 15 clocks guaranteed to yield an all-NA matrix.
tissues_present <- intersect(unique(merged$DNA_Source), DNA_SOURCES)
if (length(tissues_present) < 2) {
  cat("rank: skipping cross-tissue rank-correlation bootstrap — insufficient tissue coverage: only ",
      length(tissues_present), " tissue(s) present in this cohort (",
      paste(tissues_present, collapse = ", "), "); a cross-tissue comparison needs >=2. ",
      "This reflects the cohort's tissue scope (e.g. a single-tissue delivery wave), ",
      "not a data-quality problem.\n", sep = "")
  na_result <- list(mean = NA_real_, median = NA_real_, min = NA_real_, max = NA_real_, fisher_r = NA_real_)
  rank_results <- setNames(lapply(LME_CLOCKS, function(cl) na_result), LME_CLOCKS)
} else {
  rank_results <- setNames(lapply(LME_CLOCKS, function(cl) rank_corr_one_clock(merged, cl)), LME_CLOCKS)
}
saveRDS(rank_results, file.path(ANALYSIS_DIR, "rank_corr.rds"))
cat("rank: wrote rank_corr.rds (", length(rank_results), "clocks )\n")

## --------------------------------------------------------- Pass 2: adjusted ----
## Additive, comparison-only second pass on stage 4's cell/plate-adjusted betas.
## NOT wired into report.R/inferential — mAge_clocks.csv (Pass 1) stays the sole
## downstream input.
if (file.exists(ADJUSTED_BETAS_FILE)) {
  BetasAdj <- fread(ADJUSTED_BETAS_FILE) %>% column_to_rownames(var = "CpG")
  BetasAdj <- x_prefix_cols(BetasAdj)
  if (assert_clock_cpg_coverage(BetasAdj)) {
    mergedAdj <- compute_clocks(BetasAdj, info2)
    write.csv(mergedAdj, file.path(ANALYSIS_DIR, "mAge_clocks_adjusted.csv"), row.names = FALSE)
    cat("clocks (adjusted): wrote mAge_clocks_adjusted.csv -", ncol(mergedAdj), "cols,",
        nrow(mergedAdj), "samples\n")
  } else {
    cat("population.R: skipping the adjusted-betas clock pass — insufficient clock-CpG coverage\n")
  }
} else {
  cat("population.R: ADJUSTED_BETAS_FILE not found (", ADJUSTED_BETAS_FILE,
      ") — skipping the adjusted-betas clock pass\n")
}
