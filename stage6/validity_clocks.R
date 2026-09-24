### stage6/validity_clocks.R: sensitivity and validity checks on the clock table
### (DERIVED_DIR/mAge_clocks.csv). Base R only; no betas needed.
###   - clock-vs-age correlations recomputed on one sample per person and one per family,
###     to test for inflation from co-twin / repeated-wave non-independence;
###   - DNAmTL (dnaMethyAge clock "LuA2019") estimates telomere length in kb, so its
###     negative correlation with age is expected;
###   - ID-bridge resolution and family/wave structure (unresolved ids, individuals, families);
###   - clock-NA counts, with all-NA rows traced to stage-1 QC drops (no betas) or an
###     identity flag (clock_excluded).
### Outputs: SENS_DIR/{validity_independence,clock_na,id_resolution,allNA_samples}.csv

source("config.R")

CLOCKS_FILE <- file.path(DERIVED_DIR, "mAge_clocks.csv")
if (!file.exists(CLOCKS_FILE))
    stop("stage6/validity_clocks.R: missing ", CLOCKS_FILE, "; run stage 5 first.")
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

m      <- read.csv(CLOCKS_FILE, check.names = FALSE, stringsAsFactors = FALSE)
idcol  <- if ("aid"    %in% names(m)) "aid"    else "IndividualID"
famcol <- if ("pfamid" %in% names(m)) "pfamid" else "FamilyID"
## Clock value columns (the mAge / pace / mitotic-division outputs), excluding the paired
## *_Acceleration and meta columns; follows the naming in stage5/population.R::clock_specs.
clocks <- grep("(_mAge|_mitoticdivisions|Dunedin_Pace)$", names(m), value = TRUE)
n      <- nrow(m)

r_age <- function(d, cl) {
    x <- suppressWarnings(as.numeric(d[[cl]])); y <- suppressWarnings(as.numeric(d$Age))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(x[ok], y[ok]))
}
one_per <- function(d, key) {                 # one random row per key level (drops NA-key rows)
    k <- d[[key]]; d <- d[!is.na(k) & k != "", , drop = FALSE]
    s <- split(seq_len(nrow(d)), d[[key]])
    d[vapply(s, function(ix) if (length(ix) == 1) ix else sample(ix, 1), integer(1)), , drop = FALSE]
}

cat("=== stage 6: clock validity/sensitivity on", basename(CLOCKS_FILE),
    "(", n, "samples,", length(clocks), "clocks ) ===\n")

## ID resolution + non-independence structure ----
n_id_na  <- sum(is.na(m[[idcol]])  | m[[idcol]]  == "")
n_fam_na <- sum(is.na(m[[famcol]]) | m[[famcol]] == "")
n_indiv  <- length(unique(m[[idcol]][ !(is.na(m[[idcol]])  | m[[idcol]]  == "") ]))
n_fam    <- length(unique(m[[famcol]][!(is.na(m[[famcol]]) | m[[famcol]] == "") ]))
persons  <- m[!duplicated(m[[idcol]]), c(idcol, famcol)]
write.csv(data.frame(metric = c("samples","unresolved_id","unresolved_fam","distinct_individuals","distinct_families"),
                     value  = c(n, n_id_na, n_fam_na, n_indiv, n_fam)),
          file.path(SENS_DIR, "id_resolution.csv"), row.names = FALSE)
cat(sprintf("ID: %d samples | unresolved %s=%d %s=%d | %d individuals | %d families\n",
            n, idcol, n_id_na, famcol, n_fam_na, n_indiv, n_fam))
cat("samples-per-individual:"); print(table(as.integer(table(m[[idcol]]))))

## clock-NA propagation + stage-1 QC cross-reference ----
cm     <- vapply(clocks, function(cl) is.finite(suppressWarnings(as.numeric(m[[cl]]))), logical(n))
na_per <- colSums(!cm)
allna  <- rowSums(!cm) == length(clocks)
x2 <- data.frame(clock = clocks, n_NA = as.integer(na_per), pct_NA = round(100 * na_per / n, 2))
write.csv(x2[order(-x2$n_NA), ], file.path(SENS_DIR, "clock_na.csv"), row.names = FALSE)
bad <- m[allna, intersect(c("Sample", idcol, famcol, "Age", "Sex", "Identity_flag", "clock_excluded"), names(m))]
## A row is all-NA either because its identity is flagged (clock_excluded == TRUE) or because
## the sample was dropped at stage-1 QC and has no betas. The stage-1 sample-missingness
## report, when present, confirms the second case.
excl_col <- if ("clock_excluded" %in% names(m)) m$clock_excluded[allna] %in% TRUE else rep(FALSE, nrow(bad))
miss_f   <- file.path(REPORT_DIR, "sample_missingness.txt")
if (nrow(bad) && file.exists(miss_f)) {
    miss <- read.delim(miss_f, stringsAsFactors = FALSE)
    miss$key <- paste0(miss$Sentrix_ID, "_", miss$Sentrix_Position)
    bad$stage1_missingness <- round(miss$samp.miss.prop[match(bad$Sample, miss$key)], 4)
    bad$qc_dropped         <- bad$stage1_missingness > SAMPLE_MISSINGNESS
}
write.csv(bad, file.path(SENS_DIR, "allNA_samples.csv"), row.names = FALSE)
cat(sprintf("clock-NA: %d/%d samples fully complete; %d lose >=1 clock; %d all-NA (%d identity-flagged",
            sum(rowSums(!cm) == 0), n, sum(rowSums(!cm) >= 1), sum(allna), sum(excl_col)))
if (!is.null(bad$qc_dropped))
    cat(sprintf(", %d stage-1 QC drops", sum(bad$qc_dropped %in% TRUE & !excl_col)))
cat(")\n")

## age-validity under one-per-person / one-per-family ----
set.seed(123)
d_person <- one_per(m, idcol)
d_family <- one_per(d_person, famcol)
validity <- data.frame(
    clock        = clocks,
    mean         = round(vapply(clocks, function(cl) mean(suppressWarnings(as.numeric(m[[cl]])), na.rm = TRUE), numeric(1)), 3),
    sd           = round(vapply(clocks, function(cl) sd(suppressWarnings(as.numeric(m[[cl]])),   na.rm = TRUE), numeric(1)), 3),
    r_all        = round(vapply(clocks, function(cl) r_age(m,        cl), numeric(1)), 3),
    r_person     = round(vapply(clocks, function(cl) r_age(d_person, cl), numeric(1)), 3),
    r_family     = round(vapply(clocks, function(cl) r_age(d_family, cl), numeric(1)), 3),
    n_all        = vapply(clocks, function(cl) sum(is.finite(suppressWarnings(as.numeric(m[[cl]]))) & is.finite(m$Age)), integer(1)),
    row.names    = NULL)
validity$shift_family <- round(validity$r_family - validity$r_all, 3)
validity$note <- ifelse(validity$clock == "LuA_mAge",
                        "DNAmTL (LuA2019): telomere length in kb; negative r vs age is expected", "")
write.csv(validity, file.path(SENS_DIR, "validity_independence.csv"), row.names = FALSE)
cat(sprintf("validity: N all=%d one-per-person=%d one-per-family=%d | max |shift(all->family)|=%.3f\n",
            n, nrow(d_person), nrow(d_family), max(abs(validity$shift_family), na.rm = TRUE)))
if ("LuA_mAge" %in% clocks)
    cat(sprintf("note: LuA_mAge is DNAmTL (telomere length, range [%.2f, %.2f] kb), not an age estimate\n",
                min(suppressWarnings(as.numeric(m$LuA_mAge)), na.rm = TRUE),
                max(suppressWarnings(as.numeric(m$LuA_mAge)), na.rm = TRUE)))
cat("stage6/validity_clocks: wrote 4 tables to", SENS_DIR, "\n")
