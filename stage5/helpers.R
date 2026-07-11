### stage5/helpers.R — constants + shared functions for the core stage-5 modules.
suppressMessages({
  library(dplyr); library(tidyr); library(purrr)
})

## ---- Constants ------------------------------------------------------------
DNA_SOURCES <- c("PBMC", "Buffy_Coat", "Saliva")

## Clock mAge columns in the merged clock table (population.R).
LME_CLOCKS <- c("Dunedin_Pace", "Hannum_mAge", "Horvath_mAge", "Horvath2_mAge",
                "ZhangQ_mAge", "PhenoAge_mAge", "epiTOC_mitoticdivisions",
                "epiTOC2_mitoticdivisions", "PCGrimAge_mAge", "ZhangY_mAge",
                "LuA_mAge", "ShirebyG2020_mAge", "PedBE_mAge", "PanM2_mAge", "PanM3_mAge")

## ---- Statistical helpers --------------------------------------------------
## Fisher z-transform; the clamp keeps r = +/-1 (e.g. a correlation matrix's
## diagonal) from producing +/-Inf.
z_fisher <- function(x, eps = 1e-6) atanh(pmin(pmax(x, -1 + eps), 1 - eps))

## Family-structure code -> zygosity string (1=Adopted, 2=BiologicalSiblings,
## 3=DZ, 4=MZ; matches catslife_id_dyads.R's famtype_pair).
zygosity_from_famtype <- function(x) {
  dplyr::case_when(
    x == 1 ~ "Adopted", x == 2 ~ "BiologicalSiblings",
    x == 3 ~ "DZ",       x == 4 ~ "MZ", TRUE ~ NA_character_)
}

## ---- Per-clock cross-tissue rank bootstrap --------------------------------
## One member per family; rank the clock within each tissue and rank Age within
## PBMC; Spearman-correlate the 4 rank vectors; repeat n_iter times and aggregate.
rank_corr_one_clock <- function(df, clock_col, n_iter = 100, seed = 123) {
  set.seed(seed)
  cor_results <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    fam <- df %>% group_by(FamilyID) %>% slice_sample(n = 1) %>%
      select(FamilyID, IndividualID) %>% ungroup()
    one <- df %>% semi_join(fam, by = c("FamilyID", "IndividualID"))
    sal <- subset(one, DNA_Source == "Saliva")
    bc  <- subset(one, DNA_Source == "Buffy_Coat")
    pb  <- subset(one, DNA_Source == "PBMC")
    sal$Rank_clock <- rank(sal[[clock_col]], ties.method = "average")
    bc$Rank_clock  <- rank(bc[[clock_col]],  ties.method = "average")
    pb$Rank_clock  <- rank(pb[[clock_col]],  ties.method = "average")
    pb$Rank_Age    <- rank(pb$Age,           ties.method = "average")
    rd <- merge(sal[, c("IndividualID", "Rank_clock")],
                bc[,  c("IndividualID", "Rank_clock")],
                by = "IndividualID", suffixes = c("_Saliva", "_BC"))
    rd <- merge(rd, pb[, c("IndividualID", "Rank_clock", "Rank_Age")], by = "IndividualID")
    colnames(rd)[4] <- "Rank_clock_PBMC"
    cor_results[[i]] <- cor(rd[, 2:5], method = "spearman")
  }
  arr <- simplify2array(cor_results)
  list(mean     = apply(arr, 1:2, mean),
       median   = apply(arr, 1:2, median),
       min      = apply(arr, 1:2, min),
       max      = apply(arr, 1:2, max),
       fisher_r = tanh(apply(z_fisher(arr), 1:2, mean)))
}

## ---- Duplicate technical reliability (reliability.R) ----------------------
## Duplicate aliquots (DUPS_FILE) are the same DNA re-run and therefore
## EXCHANGEABLE (no rater/occasion order), so the reliability model is the
## one-way random-effects ICC(1,1), not the two-way ICC(2,1). Unbalanced group
## sizes use the n0 (average-group-size) correction; returns NA when fewer than
## two groups carry >= 2 finite values.
icc_oneway <- function(x, g) {
  ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- as.character(g[ok])
  tb <- table(g); x <- x[g %in% names(tb)[tb >= 2L]]; g <- g[g %in% names(tb)[tb >= 2L]]
  k <- length(unique(g)); N <- length(x)
  if (k < 2L || N <= k) return(NA_real_)
  grand <- mean(x); gmean <- tapply(x, g, mean); ns <- tapply(x, g, length)
  MSB <- sum(ns * (gmean - grand)^2) / (k - 1L)     # gmean/ns share tapply's level order
  MSW <- sum((x - gmean[g])^2) / (N - k)
  n0  <- (N - sum(ns^2) / N) / (k - 1L)             # average group size (unbalanced correction)
  denom <- MSB + (n0 - 1) * MSW
  if (denom == 0) return(NA_real_)
  (MSB - MSW) / denom
}

## Per-clock reliability: ICC(1,1) (primary); test-retest Pearson r + mean
## absolute difference over the 2-member groups (member order within a pair is
## arbitrary for exchangeable aliquots, so r is secondary); and the Bland-Altman
## repeatability coefficient 2.77 * within-group SD.
dup_reliability_one_clock <- function(d, clock_col, group_col = "DupGroupID") {
  x <- d[[clock_col]]; g <- as.character(d[[group_col]])
  ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- g[ok]
  tb <- table(g); keep <- names(tb)[tb >= 2L]
  x <- x[g %in% keep]; g <- g[g %in% keep]
  if (!length(x))
    return(data.frame(clock = clock_col, n_groups = 0L, n_samples = 0L, icc = NA_real_,
                      retest_r = NA_real_, mean_abs_diff = NA_real_, repeatability_coef = NA_real_))
  k <- length(unique(g)); N <- length(x); gmean <- tapply(x, g, mean)
  MSW <- if (N > k) sum((x - gmean[g])^2) / (N - k) else NA_real_
  pr <- names(tb)[tb == 2L]                          # order-based diagnostics on exact pairs only
  r <- NA_real_; mad <- NA_real_
  if (length(pr)) {
    a <- vapply(pr, function(p) x[g == p][1], numeric(1))
    b <- vapply(pr, function(p) x[g == p][2], numeric(1))
    mad <- mean(abs(a - b))
    if (length(pr) >= 3L && sd(a) > 0 && sd(b) > 0) r <- suppressWarnings(cor(a, b))
  }
  data.frame(clock = clock_col, n_groups = k, n_samples = N, icc = icc_oneway(x, g),
             retest_r = r, mean_abs_diff = mad,
             repeatability_coef = if (is.na(MSW)) NA_real_ else 2.77 * sqrt(MSW))
}

## ---- Wave-1 ∩ wave-2 longitudinal descriptives (longitudinal.R) -----------
## Per clock, within-person change across waves among LongitudinalGroupIDs
## sampled in BOTH waves within the SAME tissue: one row per tissue with the
## within-person wave-1<->wave-2 correlation, mean Δ epigenetic + chronological
## age, and their ratio (epigenetic years per calendar year).
wave_overlap_one_clock <- function(d, clock_col) {
  keep <- is.finite(d[[clock_col]]) & !is.na(d$LongitudinalGroupID) & d$Wave %in% c(1, 2)
  b <- d[keep, c("LongitudinalGroupID", "DNA_Source", "Wave", "Age", clock_col)]
  if (!nrow(b)) return(NULL)
  names(b)[names(b) == clock_col] <- "val"
  ## collapse any repeat sample within a (group, tissue, wave) to its mean, then pair waves
  b  <- aggregate(cbind(val, Age) ~ LongitudinalGroupID + DNA_Source + Wave, data = b, FUN = mean)
  mg <- merge(b[b$Wave == 1, ], b[b$Wave == 2, ],
              by = c("LongitudinalGroupID", "DNA_Source"), suffixes = c("_1", "_2"))
  if (!nrow(mg)) return(NULL)
  do.call(rbind, lapply(split(mg, mg$DNA_Source), function(s) {
    d_epi <- s$val_2 - s$val_1; d_chr <- s$Age_2 - s$Age_1; mchr <- mean(d_chr)
    r <- if (nrow(s) >= 3L && sd(s$val_1) > 0 && sd(s$val_2) > 0)
      suppressWarnings(cor(s$val_1, s$val_2)) else NA_real_
    data.frame(clock = clock_col, tissue = s$DNA_Source[1], n_pairs = nrow(s),
               n_individuals = length(unique(s$LongitudinalGroupID)), retest_r = r,
               mean_delta_epi = mean(d_epi), sd_delta_epi = sd(d_epi), mean_delta_chrono = mchr,
               epi_per_chrono_year = if (is.finite(mchr) && abs(mchr) > 1e-8) mean(d_epi) / mchr else NA_real_,
               row.names = NULL)
  }))
}
