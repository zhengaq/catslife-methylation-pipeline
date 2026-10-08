### stage6/sample_swaps_check.R: checks SAMPLE_SWAPS_FILE against the phenotype file and the clock
### table, without the betas.
###   - Audit (stops on a failure, after writing the report): every relabel and identity flag in the
###     phenotype file is one SAMPLE_SWAPS_FILE lists, applied as written; each relabeled array
###     resolves to a person; and mAge_clocks.csv carries the phenotype file's random_id and identity
###     flag for every sample, so a clock table computed before a change to the swaps file fails.
###   - Clock-age evidence (reported, never stops): the age clocks in SWAP_AGE_CLOCKS, calibrated
###     against chronological age on the arrays SAMPLE_SWAPS_FILE does not list, give each array a
###     consensus age estimate. A relabel is supported when that estimate fits the corrected person's
###     age and not the sheet person's (within SWAP_AGE_Z robust SDs of the calibration error), and
###     contradicted in the opposite case. The clocks see age only: a swap between two people of
###     similar age (co-twins, same-age siblings) fits both labels, and a cross-sex relabel is
###     tested through age like any other.
### Base R + haven; reads SAMPLE_SWAPS_FILE, PHENOTYPE_FILE, CLEAN_ID_FILE and
### DERIVED_DIR/mAge_clocks.csv.
### Output: REPORT_DIR/sample_swaps_check.csv (one row per SAMPLE_SWAPS_FILE row)

source("config.R")
suppressMessages(library(haven))

CLOCKS_FILE <- file.path(DERIVED_DIR, "mAge_clocks.csv")
for (f in c(SAMPLE_SWAPS_FILE, PHENOTYPE_FILE, CLEAN_ID_FILE, CLOCKS_FILE))
    if (!file.exists(f)) stop("stage6/sample_swaps_check.R: missing ", f)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

## A DNA_Source is calibrated only with at least MIN_CAL_N untouched arrays, and a clock enters its
## consensus only when it correlates with chronological age at least MIN_CAL_R among them.
MIN_CAL_N <- 30
MIN_CAL_R <- 0.5

## inputs ----
sw <- read_sample_swaps(SAMPLE_SWAPS_FILE)
ph <- read.csv(PHENOTYPE_FILE, colClasses = "character", na.strings = c("", "NA"))
need <- c("Sample", "Subject_ID", "Subject_ID_sheet", "IndividualID", "DNA_Source", "Age", "Identity_flag")
if (length(miss <- setdiff(need, names(ph))))
    stop("stage6/sample_swaps_check.R: ", PHENOTYPE_FILE, " lacks column(s) ", paste(miss, collapse = ", "),
         "; rebuild it with scripts/build/build_phenotype_file.R.")
ph$Age <- as.numeric(ph$Age)
ph$Identity_flag <- ph$Identity_flag %in% "TRUE"
ck <- read.csv(CLOCKS_FILE, check.names = FALSE, stringsAsFactors = FALSE)
if (length(miss <- setdiff(c("Sample", "random_id", "Identity_flag", "clock_excluded"), names(ck))))
    stop("stage6/sample_swaps_check.R: ", CLOCKS_FILE, " lacks column(s) ", paste(miss, collapse = ", "),
         "; rerun stage 5 on the phenotype file.")
pt <- read_sav(CLEAN_ID_FILE)
pt <- data.frame(random_id = as.integer(pt$random_id), age = as.numeric(pt$age),
                 age_w1 = as.numeric(pt$age_w1), nsex = as.numeric(pt$nsex))
pt <- pt[!is.na(pt$random_id), ]

## Age and admin sex of the person a label names; the age is the one for the label's wave
## (LabAge for wave 2, LabAge1 for wave 1), as in the phenotype file.
person_age <- function(label) {
    i <- match(subject_base_id(label), pt$random_id)
    ifelse(subject_wave(label) == 2, pt$age[i], pt$age_w1[i])
}
person_sex <- function(label) {
    s <- pt$nsex[match(subject_base_id(label), pt$random_id)]
    ifelse(s %in% 1, "M", ifelse(s %in% 0, "F", NA_character_))
}

## Audit ----
key    <- sheet_label(ph$Subject_ID_sheet)
row_of <- match(key, sw$from)                     # the swaps-file row naming each array, or NA
act    <- sw$action[row_of]
problems <- character()
fail <- function(bad, msg, what)
    if (any(bad)) problems <<- c(problems, paste0(msg, ": ", paste(unique(what[bad]), collapse = ", ")))
fail(is.na(row_of) & ph$Subject_ID != ph$Subject_ID_sheet,
     "the phenotype file relabels array(s) that SAMPLE_SWAPS_FILE does not list",
     paste(ph$Subject_ID_sheet, "->", ph$Subject_ID))
fail(act %in% "relabel" & ph$Subject_ID != sw$to[row_of],
     "relabel(s) not applied as SAMPLE_SWAPS_FILE writes them",
     paste0(ph$Subject_ID_sheet, " -> ", ph$Subject_ID, " (file: ", sw$to[row_of], ")"))
fail(act %in% "flag" & ph$Subject_ID != ph$Subject_ID_sheet,
     "identity-flagged array(s) relabeled", ph$Subject_ID_sheet)
fail(ph$Identity_flag != (act %in% "flag"),
     "Identity_flag differs from the flags in SAMPLE_SWAPS_FILE", ph$Subject_ID_sheet)
fail(act %in% "relabel" & is.na(ph$IndividualID),
     "relabeled array(s) without a person", ph$Subject_ID)

ci  <- match(ck$Sample, ph$Sample)
rid <- suppressWarnings(as.integer(ck$random_id))
fail(is.na(ci), "mAge_clocks.csv sample(s) absent from the phenotype file", ck$Sample)
in_ph <- !is.na(ci)
stale <- in_ph & (is.na(rid) | rid != subject_base_id(ph$Subject_ID)[ci] |
                  as.logical(ck$Identity_flag) != ph$Identity_flag[ci])
fail(stale %in% TRUE, "mAge_clocks.csv disagrees with the phenotype file on random_id or Identity_flag for",
     ck$Sample)
if (length(problems)) problems <- c(problems, paste(
    "rebuild the phenotype file from the current SAMPLE_SWAPS_FILE",
    "(scripts/build/build_phenotype_file.R) and rerun stage 5"))

## Clock-age consensus ----
## Each clock is calibrated within a DNA_Source by a linear fit of the clock on chronological age
## over the untouched arrays (not in SAMPLE_SWAPS_FILE, clocks not excluded), and inverted to an age
## estimate; the consensus is the median over the clocks that track age. Its error on the untouched
## arrays (median, MAD) scales every comparison below.
clocks <- intersect(SWAP_AGE_CLOCKS, names(ck))
if (!length(clocks))
    stop("stage6/sample_swaps_check.R: none of SWAP_AGE_CLOCKS (", paste(SWAP_AGE_CLOCKS, collapse = ", "),
         ") is a column of ", CLOCKS_FILE)
cx <- ck[in_ph, c("Sample", "clock_excluded", clocks)]
cx$DNA_Source <- ph$DNA_Source[ci[in_ph]]
cx$Age        <- ph$Age[ci[in_ph]]
cx$cal        <- is.na(row_of[ci[in_ph]]) & !(cx$clock_excluded %in% TRUE) & is.finite(cx$Age)
cx$age_hat    <- NA_real_
calib <- data.frame(DNA_Source = sort(unique(cx$DNA_Source)), n = NA_integer_, clocks = NA_character_,
                    centre = NA_real_, scale = NA_real_, why = NA_character_, stringsAsFactors = FALSE)
for (j in seq_len(nrow(calib))) {
    g <- cx$DNA_Source %in% calib$DNA_Source[j]; cal <- g & cx$cal
    calib$n[j] <- sum(cal)
    if (sum(cal) < MIN_CAL_N) { calib$why[j] <- paste("fewer than", MIN_CAL_N, "untouched arrays"); next }
    est <- vapply(clocks, function(k) {
        y  <- suppressWarnings(as.numeric(cx[[k]])); ok <- cal & is.finite(y)
        r  <- if (sum(ok) >= MIN_CAL_N) suppressWarnings(stats::cor(y[ok], cx$Age[ok])) else NA
        if (!isTRUE(r >= MIN_CAL_R)) return(rep(NA_real_, nrow(cx)))
        b <- stats::coef(stats::lm(y[ok] ~ cx$Age[ok]))
        (y - b[[1]]) / b[[2]]
    }, numeric(nrow(cx)))
    est <- matrix(est, nrow = nrow(cx), dimnames = list(NULL, clocks))
    use <- colSums(is.finite(est)) > 0
    if (!any(use)) { calib$why[j] <- "no clock tracks chronological age"; next }
    cx$age_hat[g] <- apply(est[g, use, drop = FALSE], 1, function(v)
        if (any(is.finite(v))) stats::median(v, na.rm = TRUE) else NA_real_)
    e <- (cx$age_hat - cx$Age)[cal]
    calib$clocks[j] <- paste(clocks[use], collapse = " ")
    calib$centre[j] <- stats::median(e, na.rm = TRUE)
    calib$scale[j]  <- stats::mad(e, na.rm = TRUE)
}
for (j in seq_len(nrow(calib))) {
    if (is.na(calib$why[j])) cat(sprintf(
        "calibration %s: %d untouched arrays | clocks %s | error median %+.2f y, robust SD %.2f y | the clocks tell two labels apart only when their ages differ by more than %.1f y\n",
        calib$DNA_Source[j], calib$n[j], calib$clocks[j], calib$centre[j], calib$scale[j],
        SWAP_AGE_Z * calib$scale[j]))
    else cat(sprintf("calibration %s: no clock-age evidence (%s)\n", calib$DNA_Source[j], calib$why[j]))
}

## Report: one row per SAMPLE_SWAPS_FILE row ----
ai <- match(sw$from, key)                         # the array each row names, NA if not in the phenotype file
xi <- match(ph$Sample[ai], cx$Sample)
cb <- match(cx$DNA_Source[xi], calib$DNA_Source)
age_sheet <- person_age(sw$from); age_corr <- person_age(sw$to)
age_hat   <- cx$age_hat[xi]
z_sheet   <- (age_hat - age_sheet - calib$centre[cb]) / calib$scale[cb]
z_corr    <- (age_hat - age_corr  - calib$centre[cb]) / calib$scale[cb]
sex_sheet <- person_sex(sw$from); sex_corr <- person_sex(sw$to)
why <- ifelse(is.na(ai), "not in the phenotype file",
       ifelse(is.na(xi) | cx$clock_excluded[xi] %in% TRUE, "no clock values",
       ifelse(!is.na(calib$why[cb]), calib$why[cb],
       ifelse(!is.finite(age_hat), "no clock values",
       ifelse(!is.finite(age_sheet), "no age for the sheet person",
       ifelse(!is.finite(age_corr), "no age for the corrected person", NA_character_))))))
fit_c <- abs(z_corr) <= SWAP_AGE_Z; fit_s <- abs(z_sheet) <= SWAP_AGE_Z
evidence <- ifelse(!is.na(why), paste("not testable:", why),
            ifelse(sw$action == "flag", ifelse(fit_c, "label fits", "label does not fit"),
            ifelse(fit_c & !fit_s, "supports relabel",
            ifelse(!fit_c & fit_s, "contradicts relabel",
            ifelse(fit_c, "both labels fit", "neither label fits")))))
out <- data.frame(
    Sample = ph$Sample[ai], Subject_ID_sheet = sw$from, Subject_ID = sw$to, action = sw$action,
    notes = sw$notes, in_phenotype = !is.na(ai),
    relabel_type = ifelse(sw$action == "flag", "flag",
                   ifelse(is.na(sex_sheet) | is.na(sex_corr), "sex unknown",
                   ifelse(sex_sheet != sex_corr, "cross-sex", "same-sex"))),
    Sex_sheet_person = sex_sheet, Sex_corrected_person = sex_corr, DNA_Source = cx$DNA_Source[xi],
    Age_sheet_person = round(age_sheet, 2), Age_corrected_person = round(age_corr, 2),
    age_clocks = round(age_hat, 2), z_sheet = round(z_sheet, 2), z_corrected = round(z_corr, 2),
    clock_evidence = evidence, stringsAsFactors = FALSE)
write.csv(out, file.path(REPORT_DIR, "sample_swaps_check.csv"), row.names = FALSE)

rl <- out$action == "relabel"
cat(sprintf("SAMPLE_SWAPS_FILE: %d relabel(s) (%d cross-sex, %d same-sex), %d flag(s); %d row(s) not in the phenotype file\n",
            sum(rl), sum(out$relabel_type == "cross-sex"), sum(out$relabel_type == "same-sex"),
            sum(!rl), sum(!out$in_phenotype)))
tb <- table(out$clock_evidence)
cat("clock-age evidence:\n", sprintf("  %3d  %s\n", as.integer(tb), names(tb)), sep = "")
for (i in which(out$clock_evidence %in% c("contradicts relabel", "neither label fits", "label does not fit")))
    cat(sprintf("WARNING: %s for %s -> %s (%s): clock age %.1f, sheet person %.1f, corrected person %.1f\n",
                out$clock_evidence[i], out$Subject_ID_sheet[i], out$Subject_ID[i], out$Sample[i],
                out$age_clocks[i], out$Age_sheet_person[i], out$Age_corrected_person[i]))
cat("stage6/sample_swaps_check: wrote sample_swaps_check.csv to", REPORT_DIR, "\n")

if (length(problems))
    stop("stage6/sample_swaps_check.R: the phenotype file or clock table does not match SAMPLE_SWAPS_FILE:\n  - ",
         paste(problems, collapse = "\n  - "))
cat("audit: the phenotype file and mAge_clocks.csv match SAMPLE_SWAPS_FILE\n")
