### config.R: central configuration for the CATSLife methylation pipeline, sourced by
### every stage. describe_paths() prints the resolved paths; validate_paths() checks them.

## Project root + site profile ----
.find_root <- function() {
    r <- Sys.getenv("METHYL_PROJECT_DIR", "")
    if (nzchar(r)) return(normalizePath(r, mustWork = FALSE))
    d <- normalizePath(getwd(), mustWork = FALSE)
    while (!file.exists(file.path(d, ".methyl-root")) && dirname(d) != d) d <- dirname(d)
    if (file.exists(file.path(d, ".methyl-root"))) d else normalizePath(getwd(), mustWork = FALSE)
}
.root <- .find_root()
.site <- Sys.getenv("METHYL_SITE_CONFIG", file.path(.root, "config.site.R"))
if (file.exists(.site)) {
    message("config: loading site profile ", .site)
    .env0 <- Sys.getenv()
    source(.site, local = TRUE)
    .env1 <- Sys.getenv()
    .set  <- names(.env1)[!names(.env1) %in% names(.env0) | .env1[names(.env1)] != .env0[names(.env1)]]
    ## config.R reads only METHYL_-prefixed variables, so a bare name (CLEAN_ID_FILE) set by the
    ## profile would be ignored and the default used instead.
    .src   <- readLines(file.path(.root, "config.R"), warn = FALSE)
    .known <- sub("^METHYL_", "", unique(unlist(regmatches(.src, gregexpr("METHYL_[A-Z0-9_]+", .src)))))
    .bare  <- intersect(.set, .known)
    if (length(.bare))
        stop("config: ", .site, " sets ", paste(.bare, collapse = ", "), " without the METHYL_ prefix; ",
             "rename to ", paste0("METHYL_", .bare, collapse = ", "), call. = FALSE)
    .ph <- .set[grepl("<[^>]+>", .env1[.set])]
    if (length(.ph))
        stop("config: ", .site, " leaves a <placeholder> in ", paste(.ph, collapse = ", "), call. = FALSE)
    rm(.env0, .env1, .set, .src, .known, .bare, .ph)
}

## Run mode ----
## ARRAY_VERSION: "v2" (EPIC v2.0, default) or "v1" (legacy EPIC v1). Env: METHYL_ARRAY_VERSION=v1
ARRAY_VERSION <- match.arg(Sys.getenv("METHYL_ARRAY_VERSION", "v2"), c("v1", "v2"))

## SAVE_INTERMEDIATES=FALSE skips stage 1's optional .RDat checkpoints (raw/detP/noob/noobflt/
## dasen; only dasen_betas.RDat is read downstream), saving tens of GB but disabling resume.
SAVE_INTERMEDIATES <- !(toupper(Sys.getenv("METHYL_SAVE_INTERMEDIATES", "TRUE")) %in% c("FALSE", "0", "NO"))

## METHYL_RESUME=TRUE makes stage 1 start from the latest checkpoint a SAVE_INTERMEDIATES=TRUE run
## left behind: F_DASEN (only beta/M extraction remains), else F_NOOBFLT (re-run dasen only), else
## F_RAW/F_DETP (skip the IDAT read / detection p-values). Default FALSE (fresh run).
RESUME <- toupper(Sys.getenv("METHYL_RESUME", "FALSE")) %in% c("TRUE", "1", "YES")

## DASEN_STREAM=TRUE (default) normalizes with dasen_stream() below, which gives the same output
## as wateRmelon::dasen at a fraction of the memory. METHYL_DASEN_STREAM=FALSE uses stock
## wateRmelon::dasen, which needs >100GB at cohort scale; use it only on small data.
DASEN_STREAM <- !(toupper(Sys.getenv("METHYL_DASEN_STREAM", "TRUE")) %in% c("FALSE", "0", "NO"))

## Directories ----
PROJECT_DIR  <- Sys.getenv("METHYL_PROJECT_DIR",  .root)
DATA_DIR     <- Sys.getenv("METHYL_DATA_DIR",     file.path(PROJECT_DIR, "data"))
ANALYSIS_DIR <- Sys.getenv("METHYL_ANALYSIS_DIR", file.path(PROJECT_DIR, "output"))
## Output sub-trees. Each defaults to ANALYSIS_DIR (one flat directory); a site profile can point
## each at its own directory to separate analysis-ready data (derived), regenerable stage 1-4
## artifacts (intermediate), and findings (results). REPORT_DIR/TABLES_DIR/SENS_DIR sit under results.
DERIVED_DIR      <- Sys.getenv("METHYL_DERIVED_DIR",      ANALYSIS_DIR)
INTERMEDIATE_DIR <- Sys.getenv("METHYL_INTERMEDIATE_DIR", ANALYSIS_DIR)
RESULTS_DIR      <- Sys.getenv("METHYL_RESULTS_DIR",      ANALYSIS_DIR)
LOGS_DIR         <- Sys.getenv("METHYL_LOGS_DIR",         file.path(ANALYSIS_DIR, "logs"))  # run logs + checkpoints
REPORT_DIR       <- file.path(RESULTS_DIR, "reports")
TABLES_DIR       <- file.path(RESULTS_DIR, "tables")
SENS_DIR         <- file.path(RESULTS_DIR, "sensitivity")
## Outputs must not land in the code checkout itself (a subdirectory such as output/ is fine).
for (.d in c(ANALYSIS_DIR, DERIVED_DIR, INTERMEDIATE_DIR, RESULTS_DIR, LOGS_DIR))
    if (normalizePath(.d, mustWork = FALSE) == normalizePath(PROJECT_DIR, mustWork = FALSE))
        stop("config: an output directory resolves to the code checkout (", PROJECT_DIR, "); set ",
             "METHYL_ANALYSIS_DIR or the METHYL_*_DIR sub-trees in the site profile", call. = FALSE)
for (.d in c(DERIVED_DIR, INTERMEDIATE_DIR, REPORT_DIR, TABLES_DIR, SENS_DIR))
    dir.create(.d, recursive = TRUE, showWarnings = FALSE)

## Inputs ----
## IDATs live one directory per Sentrix barcode under IDAT_DIR
## (Released_Data/Data/<barcode>/<barcode>_R##C##_{Grn,Red}.idat). See load_targets().
IDAT_DIR        <- Sys.getenv("METHYL_IDAT_DIR",     DATA_DIR)
SAMPLE_SHEET    <- Sys.getenv("METHYL_SAMPLE_SHEET", file.path(DATA_DIR, "sample_sheet_ids.csv"))
## Pedigree/SIF file (.xlsx: Family/Individual/Father/Mother/Sex/Subject_ID/Population).
## Used for Father/Mother + sex cross-check
ID_KEY          <- Sys.getenv("METHYL_ID_KEY",       file.path(DATA_DIR, "SIF.xlsx"))

## Intermediate outputs ----
F_RAW     <- file.path(INTERMEDIATE_DIR, "methylation_data_raw.RDat")
F_DETP    <- file.path(INTERMEDIATE_DIR, "methylation_data_detP.RDat")
F_NOOB    <- file.path(INTERMEDIATE_DIR, "methylation_data_noob.RDat")
F_NOOBFLT <- file.path(INTERMEDIATE_DIR, "methylation_data_noobflt.RDat")
F_DASEN   <- file.path(INTERMEDIATE_DIR, "methylation_data_dasen.RDat")
F_DASENB  <- file.path(INTERMEDIATE_DIR, "dasen_betas.RDat")

## QC parameters ----
DETP_THRESHOLD     <- 0.05  # detP < this  => probe call is "detected"
SAMPLE_MISSINGNESS <- 0.01  # drop sample if > this fraction of probes undetected
PROBE_MISSINGNESS  <- 0.01  # drop probe  if > this fraction of samples undetected

## Stage 2 (PCA) ----
## Default: PCA the top PCA_NCPG most-variable CpGs (standard for a structure/QC PCA;
## see the README). METHYL_PCA_SUBSET=FALSE PCAs all CpGs (heavy at cohort scale).
PCA_SUBSET <- toupper(Sys.getenv("METHYL_PCA_SUBSET", "TRUE")) %in% c("TRUE","1","YES")
PCA_NCPG   <- 5000

## Stage 3/4 (chunking) ----
NPARTS <- as.integer(Sys.getenv("METHYL_NPARTS", "5"))  # CpG chunks for stages 3 & 4
## 0 = residualize all CpGs; >0 caps the per-CpG residualization loop (a quick partial run).
RESID_CPG_LIMIT <- as.integer(Sys.getenv("METHYL_RESID_CPG_LIMIT", "0"))

## Stage 5 (clocks) inputs ----
## Min fraction of a reference clock's CpGs that must be present in the beta matrix
## before trusting clock output (population.R::assert_clock_cpg_coverage).
CLOCK_CPG_COVERAGE_MIN <- as.numeric(Sys.getenv("METHYL_CLOCK_CPG_COVERAGE_MIN", "0.90"))
PHENOTYPE_FILE <- Sys.getenv("METHYL_PHENOTYPE_FILE", file.path(DERIVED_DIR, "PhenotypeFile.csv"))
## Stage 4's cell/plate-adjusted betas; read by stage5/population.R's adjusted clock pass.
ADJUSTED_BETAS_FILE <- Sys.getenv("METHYL_ADJUSTED_BETAS_FILE",
                                  file.path(INTERMEDIATE_DIR, "B.adjusted.platebatches.txt"))
## Stage 3's cell-proportion table; an analysis-ready covariate source for PHENOTYPE_FILE.
CELL_PROPORTIONS_FILE <- Sys.getenv("METHYL_CELL_PROPORTIONS_FILE",
                                    file.path(DERIVED_DIR, "cell_proportions.blood.saliva.txt"))

## QC inputs (keyed on Subject_ID), applied in build_phenotype_file.R:
##  - DUPS_FILE: intentional duplicate pairs; both aliquots are retained and tagged with a
##    shared DupGroupID, whose agreement stage5/reliability.R quantifies.
##  - PROBLEM_HISTORY_FILE (xlsx): subjects with "problem remains at release" = yes excluded.
##  - IBD_FILE: PLINK --genome; DUPLICATED=Yes & EXPECTED=No rows flagged (not dropped).
DUPS_FILE            <- Sys.getenv("METHYL_DUPS_FILE",            file.path(DATA_DIR, "DUPS.csv"))
PROBLEM_HISTORY_FILE <- Sys.getenv("METHYL_PROBLEM_HISTORY_FILE", file.path(DATA_DIR, "Project_Problem_History.xlsx"))
IBD_FILE             <- Sys.getenv("METHYL_IBD_FILE",             file.path(DATA_DIR, "IBD_results.csv"))
## SAMPLE_SWAPS_FILE (label corrections) is defined with its helpers below.

## CATSLife person table ----
## ADMIN_FILE: the individual-admin .sav (aid/pfamid/nsex/famtype/adopted/LabAge/nidaid/...).
## SAMPLE_LIST_FILE: the per-wave array<->person crosswalk (random_id <-> nidaid + CATSLife
## wave). build_person_table.R merges the two (on nidaid) into CLEAN_ID_FILE, so the sheet's
## PI Provided Subject ID (= random_id) reaches the person world (aid/pfamid).
ADMIN_FILE       <- Sys.getenv("METHYL_ADMIN_FILE",       file.path(DATA_DIR, "individual_admin.sav"))
SAMPLE_LIST_FILE <- Sys.getenv("METHYL_SAMPLE_LIST_FILE", file.path(DATA_DIR, "sample_list.xlsx"))
## CLEAN_ID_FILE is written by build_person_table.R; DYADS_FILE by catslife_id_dyads.R.
CLEAN_ID_FILE    <- Sys.getenv("METHYL_CLEAN_ID_FILE",    file.path(DERIVED_DIR, "catslife_person_table.sav"))
DYADS_FILE       <- Sys.getenv("METHYL_DYADS_FILE",       file.path(DERIVED_DIR, "catslife_dyads.csv"))

## Vocabulary canonicalization ----
## Canonicalize tissue names ("Buffy Coat" -> "Buffy_Coat") as the sample sheet is read,
## so downstream code sees one form. Stops on an unrecognized value.
DNA_SOURCE_MAP <- c("Buffy Coat" = "Buffy_Coat", "Buffy_Coat" = "Buffy_Coat",
                     "Cell Line"  = "Cell_Line",  "Cell_Line"  = "Cell_Line",
                     "PBMC" = "PBMC", "Saliva" = "Saliva")
canonicalize_dna_source <- function(x) {
    out <- unname(DNA_SOURCE_MAP[x])
    bad <- is.na(out) & !is.na(x)
    if (any(bad))
        stop("canonicalize_dna_source: unrecognized DNA_Source value(s): ",
             paste(unique(x[bad]), collapse = ", "))
    out
}

## Digit-leading IDs (Sentrix barcodes) get an "X" prepended by R's check.names when
## used as column names; strip that leading X (R only adds it before a digit).
strip_x_prefix <- function(x) sub("^X(?=[0-9])", "", x, perl = TRUE)

## Strip a trailing "_<wave>" suffix ("14254_2" -> "14254") to recover the base Subject_ID.
## Narrow (only a trailing "_<digits>") so it won't touch a "_2D" duplicate-aliquot marker.
strip_wave_suffix <- function(x) sub("_[0-9]+$", "", x)

## Parse a Subject_ID to its integer base key (the array-facing random_id): strip the
## dup-aliquot "D" marker (so both members of a duplicate pair resolve to the same person)
## and the wave suffix, then require a pure-integer residue. Anything else (e.g. a control)
## is an error instead of a silent NA.
subject_base_id <- function(x) {
    base <- strip_wave_suffix(sub("_[0-9]*D$", "", x))
    bad  <- !grepl("^[0-9]+$", base)
    if (any(bad))
        stop("subject_base_id: non-numeric Subject_ID base(s): ",
             paste(unique(x[bad]), collapse = ", "))
    as.integer(base)
}

## Parse the wave from a Subject_ID's "_<digits>" suffix (optionally with a "D" dup marker):
## "14254" -> 1, "14254_2" -> 2, "557_2D" -> 2, "11747_D" -> 1 (no digits -> wave 1).
subject_wave <- function(x) {
    w <- sub("^.*_([0-9]+)D?$", "\\1", x)
    ifelse(w == x, 1L, suppressWarnings(as.integer(w)))
}

## Classify an IBD-flagged genetic-duplicate pair (DUPLICATED=Yes & EXPECTED=No) against the
## person table (columns random_id, pfamid, ZygGroup): "cross_wave" = same wave-stripped
## Subject_ID base (one participant resampled across waves); "mz" = different persons in the
## same pfamid with ZygGroup==1 (MZ co-twins, genetically identical by design); "unexpected" =
## anything else (a likely sample swap/mislabel needing manual review). Vectorized over s1/s2.
classify_ibd_pair <- function(s1, s2, person) {
    r1 <- match(subject_base_id(s1), person$random_id)
    r2 <- match(subject_base_id(s2), person$random_id)
    same_base <- strip_wave_suffix(s1) == strip_wave_suffix(s2)
    is_mz <- !is.na(r1) & !is.na(r2) & person$pfamid[r1] == person$pfamid[r2] &
             person$ZygGroup[r1] %in% 1 & person$ZygGroup[r2] %in% 1
    ifelse(same_base, "cross_wave", ifelse(is_mz, "mz", "unexpected"))
}

## Sample-identity corrections ----
## SAMPLE_SWAPS_FILE lists sample-sheet labels (PI Provided Subject ID) found to be wrong, in
## columns "Incorrect Random ID" and "Correct Random ID" (a swap is two rows, one per sample).
## Each row is one of:
##   relabel: the sample labeled <incorrect> is really <correct>;
##   exclude: <correct> is UNKNOWN_RANDOM_ID (the sample's person is unknown); the sample is dropped;
##   flag:    <incorrect> == <correct> (identity doubtful); the sample is kept with Identity_flag,
##            and its clocks are NA-filled (clock_excluded) while EXCLUDE_IDENTITY_FLAGGED is TRUE.
## A cohort with no corrections supplies the file with its header only.
SAMPLE_SWAPS_FILE <- Sys.getenv("METHYL_SAMPLE_SWAPS_FILE", file.path(DATA_DIR, "sample_swaps.csv"))
UNKNOWN_RANDOM_ID <- 99999L
EXCLUDE_IDENTITY_FLAGGED <- !(toupper(Sys.getenv("METHYL_EXCLUDE_IDENTITY_FLAGGED", "TRUE")) %in% c("FALSE", "0", "NO"))

## Read SAMPLE_SWAPS_FILE into data.frame(from, to, action, notes). Header names are
## whitespace-trimmed; blank rows are ignored.
read_sample_swaps <- function(path = SAMPLE_SWAPS_FILE) {
    s <- read.csv(path, check.names = FALSE, colClasses = "character", na.strings = character(0))
    names(s) <- trimws(names(s))
    need <- c("Incorrect Random ID", "Correct Random ID")
    if (!all(need %in% names(s)))
        stop("read_sample_swaps: ", path, " needs columns ", paste0('"', need, '"', collapse = " and "))
    from  <- trimws(s[["Incorrect Random ID"]]); to <- trimws(s[["Correct Random ID"]])
    notes <- if ("Notes" %in% names(s)) trimws(s[["Notes"]]) else rep("", length(from))
    keep  <- nzchar(from) | nzchar(to)
    from  <- from[keep]; to <- to[keep]; notes <- notes[keep]
    if (any(!nzchar(from) | !nzchar(to)))
        stop("read_sample_swaps: a row gives only one of the incorrect/correct ids")
    if (anyDuplicated(from))
        stop("read_sample_swaps: id(s) listed more than once as incorrect: ",
             paste(unique(from[duplicated(from)]), collapse = ", "))
    action <- ifelse(from == to, "flag",
                     ifelse(strip_wave_suffix(to) == as.character(UNKNOWN_RANDOM_ID), "exclude", "relabel"))
    if (anyDuplicated(to[action == "relabel"]))
        stop("read_sample_swaps: two samples relabeled to the same id")
    data.frame(from = from, to = to, action = action, notes = notes, stringsAsFactors = FALSE)
}

## Apply read_sample_swaps() output to the sample sheet's labels `sid`. Returns
## list(subject_id = corrected labels, action = per-sample action or NA). Every listed label
## must occur exactly once on the sheet. A duplicate-aliquot label ("<id>D") of a listed id
## stops the build, because the file does not say whether the aliquot needs the same correction.
apply_sample_swaps <- function(sid, swaps) {
    n <- vapply(swaps$from, function(f) sum(sid == f, na.rm = TRUE), integer(1))
    if (any(n != 1))
        stop("apply_sample_swaps: listed id(s) not found exactly once on the sample sheet: ",
             paste0(swaps$from[n != 1], " (", n[n != 1], "x)", collapse = ", "))
    aliquot <- sid[!sid %in% swaps$from & sub("_?D$", "", sid) %in% swaps$from]
    if (length(aliquot))
        stop("apply_sample_swaps: duplicate aliquot(s) of a listed id; list them in the swaps file too: ",
             paste(aliquot, collapse = ", "))
    i      <- match(sid, swaps$from)
    action <- swaps$action[i]
    out    <- ifelse(action %in% "relabel", swaps$to[i], sid)
    clash  <- intersect(out[duplicated(out)], swaps$to[swaps$action == "relabel"])
    if (length(clash))
        stop("apply_sample_swaps: relabeling leaves more than one sample labeled ",
             paste(clash, collapse = ", "), "; is the other half of a swap missing?")
    list(subject_id = out, action = action)
}

## EPIC v2 gives some replicate probes an id suffix ("cg#######_TC21"); clock and
## cell-type references key on the bare "cg########" id, so strip the suffix. Where
## two rows collapse to one bare id, keep the lower-missingness row. No-op for v1.
canonicalize_v2_probe_ids <- function(betas, array_version = ARRAY_VERSION) {
    if (array_version != "v2") return(betas)
    ids  <- rownames(betas)
    bare <- sub("_[A-Za-z0-9]+$", "", ids)
    if (identical(bare, ids)) return(betas)  # nothing suffixed; no-op
    miss <- rowSums(is.na(betas))
    ord  <- order(bare, miss)                 # lowest missingness first per bare id
    betas <- betas[ord, , drop = FALSE]; bare <- bare[ord]
    keep_idx_v2 <- !duplicated(bare)           # first (=lowest-missingness) occurrence
    betas <- betas[keep_idx_v2, , drop = FALSE]
    rownames(betas) <- bare[keep_idx_v2]
    betas
}

## Helpers ----
## Load an .RDat holding a single object, returning it regardless of its name.
load_one <- function(path) {
    e  <- new.env(parent = emptyenv())
    nm <- load(path, envir = e)
    e[[nm[1]]]
}

## Safe negative indexing: x[-drop] selects NOTHING when drop is empty
## (seq_len(n)[-integer(0)] == integer(0)). Returns the indices to KEEP, so
## x[keep_idx(nrow(x), rm), ] works whether or not anything is dropped.
keep_idx <- function(total, drop) if (length(drop)) seq_len(total)[-drop] else seq_len(total)

## Detection p-values in sample-sized batches (exact: each sample is scored vs its own control
## background), to cap the ~60GB peak of an all-at-once call at cohort scale. minfi::-qualified
## since ewastools also exports detectionP. Batch size via METHYL_DETP_CHUNK; chunk >= ncol = one call.
detectionP_chunked <- function(rgSet, chunk = as.integer(Sys.getenv("METHYL_DETP_CHUNK", "200"))) {
    n <- ncol(rgSet)
    if (is.na(chunk) || chunk < 1L || chunk >= n) return(minfi::detectionP(rgSet))
    idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
    cat("detectionP_chunked:", n, "samples in", length(idx), "batch(es) of up to", chunk, "\n")
    parts <- lapply(seq_along(idx), function(k) {
        cat("  detectionP batch", k, "/", length(idx), "\n"); utils::flush.console()
        minfi::detectionP(rgSet[, idx[[k]], drop = FALSE])
    })
    do.call(cbind, parts)[, colnames(rgSet), drop = FALSE]   # realign to original column order
}

## preprocessNoob in sample-batches. noob (dyeMethod "single", the default) is per-sample, so
## processing sample-subsets and reassembling is exact. The Meth/Unmeth matrices are preallocated
## and filled in place, so peak memory stays near one full MethylSet (~25GB plus the input) where
## an all-at-once call needs >70GB at ~1600 samples. Batch size via METHYL_NOOB_CHUNK. (dasen is
## cross-sample and cannot be batched this way; dasen_stream() below streams it instead.)
preprocessNoob_chunked <- function(rgSet, chunk = as.integer(Sys.getenv("METHYL_NOOB_CHUNK", "200"))) {
    n <- ncol(rgSet)
    if (is.na(chunk) || chunk < 1L || chunk >= n) return(minfi::preprocessNoob(rgSet, verbose = TRUE))
    idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
    cat("preprocessNoob_chunked:", n, "samples in", length(idx), "batch(es) of up to", chunk, "\n")
    M <- U <- NULL; pmeth <- ""
    for (k in seq_along(idx)) {
        cat("  noob batch", k, "/", length(idx), "\n"); utils::flush.console()
        ms <- minfi::preprocessNoob(rgSet[, idx[[k]], drop = FALSE])
        if (is.null(M)) {
            M <- matrix(NA_real_, nrow(ms), n, dimnames = list(rownames(ms), colnames(rgSet)))
            U <- M
            pmeth <- minfi::preprocessMethod(ms)
        }
        stopifnot(identical(rownames(ms), rownames(M)))
        M[, idx[[k]]] <- minfi::getMeth(ms)
        U[, idx[[k]]] <- minfi::getUnmeth(ms)
        rm(ms); gc()
    }
    minfi::MethylSet(Meth = M, Unmeth = U, annotation = minfi::annotation(rgSet),
                     preprocessMethod = pmeth)
}

## Streaming dasen: a drop-in for wateRmelon::dasen(MethylSet) with the same output to machine
## precision, peaking at ~40-60GB where stock dasen needs >100GB at ~1600 samples. Stock memory
## goes to limma::normalizeQuantiles, which holds the input submatrix, a full sorted copy and a
## full rank matrix at once. Quantile normalization needs none of them: build the reference
## distribution one column at a time (pass 1), then map each column onto it in place (pass 2).
##
## qn_stream: limma::normalizeQuantiles(A, ties=TRUE), streamed. refcols selects the columns that
## define the reference distribution (all columns = cohort average, the stock behavior; a subset,
## e.g. wave 1, anchors the reference on it). Every column is mapped onto that reference. Ties are
## resolved with rank() + approx() interpolation (ties.method="average"), as limma's ties=TRUE
## path does; assigning by sorted position would give different values on tied data.
qn_stream <- function(A, refcols) {
    n1 <- nrow(A); i <- (0:(n1 - 1)) / (n1 - 1)
    acc <- numeric(n1)
    for (j in refcols) acc <- acc + sort.int(A[, j], method = "quick")   # pass 1: reference
    m <- acc / length(refcols)
    for (j in seq_len(ncol(A))) {                                        # pass 2: apply in place
        r <- rank(A[, j])                                               # ties.method="average"
        A[, j] <- approx(i, m, (r - 1) / (n1 - 1), ties = list("ordered", mean))$y
    }
    A
}

## dfsfit: per-sample background offset applied to Type I probes, optionally smoothed across
## samples by an lm on Sentrix row/column (roco). wateRmelon:::dfs2 (the internal stock dasen
## uses) gives one scalar per sample, so the lm runs over a length-n vector and is cheap. roco is
## parsed as stock dasen parses it. If the position model cannot be fit (e.g. colnames without a
## Sentrix position), smoothing is skipped with a message, as stock dasen effectively does.
dfsfit_stream <- function(mn, onetwo, roco) {
    mdf <- vapply(seq_len(ncol(mn)), function(j) wateRmelon:::dfs2(mn[, j], onetwo), numeric(1))
    if (!is.null(roco)) {
        scol <- as.numeric(substr(roco, 6, 6)); srow <- as.numeric(substr(roco, 3, 3))
        fit <- try(lm(mdf ~ srow + scol), silent = TRUE)
        if (!inherits(fit, "try-error")) mdf <- fit$fitted.values
        else message("dfsfit_stream: Sentrix position model failed, skipping roco smoothing")
    }
    isI <- onetwo == "I"
    mn[isI, ] <- mn[isI, ] - matrix(rep(mdf, sum(isI)), byrow = TRUE, nrow = sum(isI))
    mn
}

## dasen_stream: dfsfit (roco on Meth, none on Unmeth), then quantile-normalize each channel x
## probe type. reference = NULL uses the cohort average (all samples), as stock dasen does. Passing
## reference = <the wave-1 columns> (integer indices or a logical mask) estimates the reference
## from those columns and still maps every sample onto it. Returns a minfi MethylSet; getBeta()/
## getM()/betas() apply offset 100 by default, which equals dasen's default fudge=100, so no fudge
## is applied here. wateRmelon:::got is the probe design-type accessor stock dasen uses.
dasen_stream <- function(mset, reference = NULL) {
    mns <- minfi::getMeth(mset); uns <- minfi::getUnmeth(mset)
    onetwo <- wateRmelon:::got(mset)
    if (anyNA(mns) || anyNA(uns))
        stop("dasen_stream: NA intensities; the streaming fast path assumes complete data")
    refcols <- if (is.null(reference)) seq_len(ncol(mns))
               else if (is.logical(reference)) which(reference) else reference
    roco <- substring(colnames(mns), regexpr("R0[1-9]C0[1-9]", colnames(mns)))
    mns <- dfsfit_stream(mns, onetwo, roco = roco)
    uns <- dfsfit_stream(uns, onetwo, roco = NULL); gc()
    for (t in c("I", "II")) {
        r <- onetwo == t
        mns[r, ] <- qn_stream(mns[r, , drop = FALSE], refcols)
        uns[r, ] <- qn_stream(uns[r, , drop = FALSE], refcols); gc()
    }
    minfi::MethylSet(Meth = mns, Unmeth = uns, annotation = minfi::annotation(mset),
                     preprocessMethod = minfi::preprocessMethod(mset))
}

## GenomeStudio sample sheets carry a [Header]/[Manifests]/[Data] preamble; the table
## starts after "[Data]". A flat sheet with no marker parses from the top (skip=0).
## check.names=FALSE preserves headers with spaces ("PI Provided Subject ID").
read_sample_sheet <- function(path) {
    preamble <- readLines(path, n = 200L, warn = FALSE)
    marker   <- grep("^\\[Data\\]", preamble)
    skip     <- if (length(marker)) marker[1] else 0L
    read.csv(path, skip = skip, check.names = FALSE, stringsAsFactors = FALSE)
}

## Sample-sheet / targets data frame.
load_targets <- function() {
    t <- read_sample_sheet(SAMPLE_SHEET)
    ## IDATs live one directory per Sentrix barcode (Released_Data/Data/<barcode>/...);
    ## barcode = Sample_Group before the first "_".
    shape_ok <- grepl("^[0-9]+_R[0-9]+C[0-9]+$", t$Sample_Group)
    if (any(!shape_ok))
        stop("load_targets: Sample_Group not in <barcode>_R##C## form for ",
             sum(!shape_ok), " row(s): ",
             paste(utils::head(t$Sample_Group[!shape_ok], 5), collapse = ", "))

    barcode     <- sub("_.*$", "", t$Sample_Group)
    t$Basename  <- file.path(IDAT_DIR, barcode, t$Sample_Group)
    missing_dir <- !dir.exists(file.path(IDAT_DIR, unique(barcode)))
    if (any(missing_dir))
        stop("load_targets: barcode subdirectory missing under IDAT_DIR (", IDAT_DIR, "): ",
             paste(unique(barcode)[missing_dir], collapse = ", "))
    t
}

## Raw RGChannelSet from the IDATs named in the sample sheet.
load_raw_rgSet <- function() {
    minfi::read.metharray.exp(targets = load_targets(), force = TRUE)
}

## Path bridge: see + validate the logical -> physical mapping ----
## Registry of the logical paths, each tagged role (input/output/root) and the
## entry point that needs it. Drives describe_paths() and validate_paths().
.path_registry <- function() {
    row <- function(name, path, role, stage) data.frame(name = name, role = role,
                                                         stage = stage, path = path,
                                                         stringsAsFactors = FALSE)
    rbind(
        row("PROJECT_DIR",     PROJECT_DIR,     "root",   "all"),
        row("DATA_DIR",        DATA_DIR,        "input",  "all"),
        row("ANALYSIS_DIR",    ANALYSIS_DIR,    "output", "all"),
        row("DERIVED_DIR",     DERIVED_DIR,     "output", "all"),
        row("INTERMEDIATE_DIR", INTERMEDIATE_DIR, "output", "all"),
        row("RESULTS_DIR",     RESULTS_DIR,     "output", "all"),
        row("REPORT_DIR",      REPORT_DIR,      "output", "all"),
        row("IDAT_DIR",        IDAT_DIR,        "input",  "stage1"),
        row("SAMPLE_SHEET",    SAMPLE_SHEET,    "input",  "stage1"),
        row("ID_KEY",          ID_KEY,          "input",  "stage1"),
        row("ADMIN_FILE",       ADMIN_FILE,       "input",  "person_table"),
        row("SAMPLE_LIST_FILE", SAMPLE_LIST_FILE, "input",  "person_table"),
        row("CLEAN_ID_FILE",    CLEAN_ID_FILE,    "output", "person_table"),
        row("ID_KEY",          ID_KEY,          "input",  "phenotype_bridge"),
        row("SAMPLE_SHEET",    SAMPLE_SHEET,    "input",  "phenotype_bridge"),
        row("CLEAN_ID_FILE",   CLEAN_ID_FILE,   "input",  "phenotype_bridge"),
        row("DYADS_FILE",      DYADS_FILE,      "input",  "phenotype_bridge"),
        row("CELL_PROPORTIONS_FILE", CELL_PROPORTIONS_FILE, "input", "phenotype_bridge"),
        row("DUPS_FILE",       DUPS_FILE,       "input",  "phenotype_bridge"),
        row("PROBLEM_HISTORY_FILE", PROBLEM_HISTORY_FILE, "input", "phenotype_bridge"),
        row("IBD_FILE",        IBD_FILE,        "input",  "phenotype_bridge"),
        row("SAMPLE_SWAPS_FILE", SAMPLE_SWAPS_FILE, "input", "phenotype_bridge"),
        row("PHENOTYPE_FILE",  PHENOTYPE_FILE,  "output", "phenotype_bridge"),
        ## ADJUSTED_BETAS_FILE is optional: population.R skips the adjusted pass without it.
        row("PHENOTYPE_FILE",  PHENOTYPE_FILE,  "input",  "stage5"))
}

.path_status <- function(path, role) {
    if (role %in% c("output", "root")) {
        if (!dir.exists(path)) return(if (role == "root") "MISSING" else "absent (made at run)")
        if (file.access(path, 2) == 0) "writable" else "NOT writable"
    } else {
        if (!file.exists(path)) return("MISSING")
        if (file.access(path, 4) == 0) "readable" else "NOT readable"
    }
}

## Print the resolved logical -> physical mapping with status. Run on the server
## to confirm the scripts will read/write where you expect.
describe_paths <- function() {
    reg <- .path_registry()
    reg$status <- mapply(.path_status, reg$path, reg$role)
    cat("CATSLife methylation paths\n")
    cat("  project root:", PROJECT_DIR, "\n")
    cat("  site profile:", if (file.exists(.site)) .site else "(none)", "\n\n")
    w <- max(nchar(reg$name))
    for (i in seq_len(nrow(reg)))
        cat(sprintf("  %-*s  %-6s %-6s  %-20s  %s\n", w, reg$name[i], reg$role[i],
                    reg$stage[i], reg$status[i], reg$path[i]))
    invisible(reg)
}

## Pre-run check for one entry point (a stage in the registry, e.g. "stage1", "person_table",
## "phenotype_bridge", "stage5") or "all": its inputs must exist and ANALYSIS_DIR be writable.
validate_paths <- function(stage = "all") {
    reg <- .path_registry()
    stages <- setdiff(unique(reg$stage), "all")
    if (!stage %in% c("all", stages))
        stop("validate_paths: unknown stage \"", stage, "\"; one of: ",
             paste(c("all", stages), collapse = ", "), call. = FALSE)
    problems <- character(0)
    dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
    if (file.access(ANALYSIS_DIR, 2) != 0)
        problems <- c(problems, paste("ANALYSIS_DIR not writable:", ANALYSIS_DIR))
    want <- if (stage == "all") reg$stage else c(stage, "all")
    need <- reg[reg$role == "input" & reg$stage %in% want, ]
    for (i in seq_len(nrow(need)))
        if (!file.exists(need$path[i]))
            problems <- c(problems, paste0("missing input ", need$name[i], ": ", need$path[i]))
    if (length(problems))
        stop("validate_paths(\"", stage, "\") failed:\n  - ",
             paste(problems, collapse = "\n  - "),
             "\nRun describe_paths() and set the paths in config.site.R.", call. = FALSE)
    cat("validate_paths(\"", stage, "\"): OK\n", sep = "")
    invisible(TRUE)
}

## Fail-loud validation for build_phenotype_file.R: every non-control sample must resolve
## to a person (IndividualID = admin aid) and family (FamilyID = pfamid); no duplicate
## samples; recognized DNA_Source; the cohort overlaps the dyad table; no excluded
## Subject_ID slipped through. Prints the crosswalk match rate; stops with an actionable
## message. `excluded_subject_ids` = Subject_IDs the caller already tried to exclude.
validate_phenotype_bridge <- function(pheno, dyads, excluded_subject_ids = character(0)) {
    problems <- character(0)
    if (length(excluded_subject_ids) && "Subject_ID" %in% names(pheno)) {
        slipped <- intersect(excluded_subject_ids, pheno$Subject_ID)
        if (length(slipped))
            problems <- c(problems, paste0("excluded Subject_ID(s) still present in the assembled ",
                "phenotype file (DUPS_FILE/PROBLEM_HISTORY_FILE exclusion didn't take): ",
                paste(slipped, collapse = ", ")))
    }
    ## A miss => NA IndividualID => the sheet random_id <-> person crosswalk failed.
    resolved <- sum(!is.na(pheno$IndividualID)); total <- nrow(pheno)
    cat(sprintf("validate_phenotype_bridge(): crosswalk match rate %d/%d (%.1f%%)\n",
                resolved, total, if (total) 100 * resolved / total else 0))
    if (resolved < total) {
        unresolved <- if ("Subject_ID" %in% names(pheno)) pheno$Subject_ID[is.na(pheno$IndividualID)]
                      else pheno$Sample[is.na(pheno$IndividualID)]
        problems <- c(problems, paste0(total - resolved,
            " sample(s) whose random_id did not resolve to a person (sheet random_id <-> person crosswalk failed): ",
            paste(unresolved, collapse = ", ")))
    }
    ## Age can be legitimately NA (LabAge/LabAge1 are sparse for some persons), so a
    ## resolved-but-NA Age is a note, not an error: those samples have no age acceleration.
    age.na <- "Age" %in% names(pheno) & is.na(pheno$Age) & !is.na(pheno$IndividualID)
    if (any(age.na))
        cat("validate_phenotype_bridge(): note -", sum(age.na),
            "resolved sample(s) have no LabAge (age acceleration not computable for them)\n")
    fam.na <- is.na(pheno$FamilyID) & !is.na(pheno$IndividualID)
    if (any(fam.na))
        problems <- c(problems, paste0(sum(fam.na),
            " sample(s) resolved an aid but no FamilyID (pfamid missing in the person table): ",
            paste(pheno$Sample[fam.na], collapse = ", ")))
    dup <- duplicated(pheno$Sample)
    if (any(dup))
        problems <- c(problems, paste0("duplicate Sample id(s): ",
            paste(unique(pheno$Sample[dup]), collapse = ", ")))
    bad_dna <- setdiff(unique(pheno$DNA_Source), DNA_SOURCES)
    if (length(bad_dna))
        problems <- c(problems, paste0("unrecognized DNA_Source value(s) after canonicalization: ",
            paste(bad_dna, collapse = ", ")))
    have_ids <- unique(stats::na.omit(pheno$IndividualID))
    dyad_ids <- unique(c(dyads$aid.x, dyads$aid.y))
    if (nrow(dyads) > 0 && length(intersect(have_ids, dyad_ids)) == 0)
        problems <- c(problems, paste0("none of the ", length(have_ids), " phenotype-file individuals ",
            "appear in DYADS_FILE's aid.x/aid.y; the random_id <-> person crosswalk is very likely ",
            "keyed wrong (heritability/twin-corr would silently degrade to no family structure)"))
    if (length(problems))
        stop("validate_phenotype_bridge() failed:\n  - ", paste(problems, collapse = "\n  - "),
             "\nCheck the sample sheet's PI Provided Subject ID (random_id) against SAMPLE_LIST_FILE ",
             "(inspect via scripts/build/inspect_2026_delivery.R's crosswalk dry-run).", call. = FALSE)
    cat("validate_phenotype_bridge(): OK -", nrow(pheno), "samples,", length(have_ids), "individuals\n")
    invisible(TRUE)
}

## Test profile ----
## METHYL_TEST_PROFILE names an R file sourced last, so it can override the constants and
## functions above (e.g. to inject a synthetic RGChannelSet). Leave it unset for real runs.
.test_profile <- Sys.getenv("METHYL_TEST_PROFILE", "")
if (nzchar(.test_profile)) {
    if (!file.exists(.test_profile) && file.exists(file.path(.root, .test_profile)))
        .test_profile <- file.path(.root, .test_profile)
    if (!file.exists(.test_profile))
        stop("METHYL_TEST_PROFILE set but not found: ", .test_profile)
    message("config: loading test profile ", .test_profile)
    source(.test_profile)
}
