### config.R — central configuration for the CATSLife methylation pipeline,
### sourced by every stage. describe_paths() and
### validate_paths() show + check the resolved mapping.

## ---- Project root + site profile ------------------------------------------
.find_root <- function() {
    r <- Sys.getenv("METHYL_PROJECT_DIR", "")
    if (nzchar(r)) return(normalizePath(r, mustWork = FALSE))
    d <- normalizePath(getwd(), mustWork = FALSE)
    while (!file.exists(file.path(d, ".methyl-root")) && dirname(d) != d) d <- dirname(d)
    if (file.exists(file.path(d, ".methyl-root"))) d else normalizePath(getwd(), mustWork = FALSE)
}
.root <- .find_root()
.site <- Sys.getenv("METHYL_SITE_CONFIG", file.path(.root, "config.site.R"))
if (file.exists(.site)) { message("config: loading site profile ", .site); source(.site, local = TRUE) }

## ---- Run mode -------------------------------------------------------------
## ARRAY_VERSION: "v2" (EPIC v2.0, default) or "v1" (legacy EPIC v1). Env: METHYL_ARRAY_VERSION=v1
ARRAY_VERSION <- match.arg(Sys.getenv("METHYL_ARRAY_VERSION", "v2"), c("v1", "v2"))

## SAVE_INTERMEDIATES=FALSE skips stage 1's optional .RDat checkpoints (raw/detP/noob/dasen and
## the noobflt resume bundle; only dasen_betas.RDat is read downstream), saving tens of GB.
## F_NOOBFLT is the resume point: with it saved, METHYL_RESUME=TRUE re-runs only dasen. Default TRUE.
SAVE_INTERMEDIATES <- !(toupper(Sys.getenv("METHYL_SAVE_INTERMEDIATES", "TRUE")) %in% c("FALSE", "0", "NO"))

## METHYL_RESUME=TRUE makes stage 1 reload the saved raw/detP checkpoints (from a prior
## SAVE_INTERMEDIATES=TRUE run) instead of re-reading IDATs / recomputing detP — handy for
## resuming after a crash without repeating the slow read. Default FALSE (fresh run).
RESUME <- toupper(Sys.getenv("METHYL_RESUME", "FALSE")) %in% c("TRUE", "1", "YES")

## DASEN_STREAM=TRUE (default) makes stage 1 normalize with the streaming dasen_stream() below
## (memory-safe, output identical to wateRmelon::dasen — see its header). METHYL_DASEN_STREAM=FALSE
## restores stock wateRmelon::dasen, which is cross-sample and needs >100GB at cohort scale — for
## A/B verification on small data only, never for the full delivery.
DASEN_STREAM <- !(toupper(Sys.getenv("METHYL_DASEN_STREAM", "TRUE")) %in% c("FALSE", "0", "NO"))

## ---- Directories ----------------------------------------------------------
PROJECT_DIR  <- Sys.getenv("METHYL_PROJECT_DIR",  .root)
DATA_DIR     <- Sys.getenv("METHYL_DATA_DIR",     file.path(PROJECT_DIR, "data"))
ANALYSIS_DIR <- Sys.getenv("METHYL_ANALYSIS_DIR", file.path(PROJECT_DIR, "output"))
REPORT_DIR   <- file.path(ANALYSIS_DIR, "reports")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- Inputs ---------------------------------------------------------------
## IDATs live one directory per Sentrix barcode under IDAT_DIR
## (Released_Data/Data/<barcode>/<barcode>_R##C##_{Grn,Red}.idat). See load_targets().
IDAT_DIR        <- Sys.getenv("METHYL_IDAT_DIR",     DATA_DIR)
SAMPLE_SHEET    <- Sys.getenv("METHYL_SAMPLE_SHEET", file.path(DATA_DIR, "sample_sheet_ids.csv"))
## Pedigree/SIF file (.xlsx: Family/Individual/Father/Mother/Sex/Subject_ID/Population).
## Used for Father/Mother + sex cross-check
ID_KEY          <- Sys.getenv("METHYL_ID_KEY",       file.path(DATA_DIR, "SIF.xlsx"))

## ---- Intermediate outputs -------------------------------------------------
F_RAW     <- file.path(ANALYSIS_DIR, "methylation_data_raw.RDat")
F_DETP    <- file.path(ANALYSIS_DIR, "methylation_data_detP.RDat")
F_RGFLT   <- file.path(ANALYSIS_DIR, "methylation_data_rgSetflt.RDat")
F_NOOB    <- file.path(ANALYSIS_DIR, "methylation_data_noob.RDat")
F_NOOBFLT <- file.path(ANALYSIS_DIR, "methylation_data_noobflt.RDat")
F_DASEN   <- file.path(ANALYSIS_DIR, "methylation_data_dasen.RDat")
F_DASENB  <- file.path(ANALYSIS_DIR, "dasen_betas.RDat")

## ---- QC parameters --------------------------------------------------------
DETP_THRESHOLD     <- 0.05  # detP < this  => probe call is "detected"
SAMPLE_MISSINGNESS <- 0.01  # drop sample if > this fraction of probes undetected
PROBE_MISSINGNESS  <- 0.01  # drop probe  if > this fraction of samples undetected

## ---- Stage 2 (PCA) --------------------------------------------------------
## Default: PCA the top PCA_NCPG most-variable CpGs (standard for a structure/QC PCA;
## see README.runtime.md). METHYL_PCA_SUBSET=FALSE PCAs all CpGs (heavy at cohort scale).
PCA_SUBSET <- toupper(Sys.getenv("METHYL_PCA_SUBSET", "TRUE")) %in% c("TRUE","1","YES")
PCA_NCPG   <- 5000

## ---- Stage 3/4 (chunking) -------------------------------------------------
NPARTS <- as.integer(Sys.getenv("METHYL_NPARTS", "5"))  # CpG chunks for stages 3 & 4
## 0 = residualize all CpGs; >0 caps the per-CpG residualization loop (a quick partial run).
RESID_CPG_LIMIT <- as.integer(Sys.getenv("METHYL_RESID_CPG_LIMIT", "0"))
## ---- Stage 5 (clocks) inputs ----------------------------------------------
## Min fraction of a reference clock's CpGs that must be present in the beta matrix
## before trusting clock output (population.R::assert_clock_cpg_coverage).
CLOCK_CPG_COVERAGE_MIN <- as.numeric(Sys.getenv("METHYL_CLOCK_CPG_COVERAGE_MIN", "0.90"))
BETAS_FILE     <- Sys.getenv("METHYL_BETAS_FILE",     file.path(DATA_DIR, "BetasFile.csv"))
PHENOTYPE_FILE <- Sys.getenv("METHYL_PHENOTYPE_FILE", file.path(DATA_DIR, "PhenotypeFile.csv"))
## HORVATH_CPGS: declared for the fallback monolith only; no stage5/ module reads it.
HORVATH_CPGS   <- Sys.getenv("METHYL_HORVATH_CPGS",   file.path(DATA_DIR, "HorvathCpGsFile.csv"))
## Stage 4's cell/plate-adjusted betas; read by stage5/population.R's adjusted clock pass.
ADJUSTED_BETAS_FILE <- Sys.getenv("METHYL_ADJUSTED_BETAS_FILE",
                                  file.path(ANALYSIS_DIR, "B.adjusted.platebatches.txt"))
## Stage 3's cell-proportion report; source of PHENOTYPE_FILE's cell_* covariates.
CELL_PROPORTIONS_FILE <- Sys.getenv("METHYL_CELL_PROPORTIONS_FILE",
                                    file.path(REPORT_DIR, "cell_proportions.blood.saliva.txt"))

## QC inputs (keyed on Subject_ID), applied in build_phenotype_file.R:
##  - DUPS_FILE: intentional duplicate pairs; BOTH aliquots are retained and tagged with a
##    shared DupGroupID (a consistency check; stage5/reliability.R quantifies it), not dropped.
##  - PROBLEM_HISTORY_FILE (xlsx): subjects with "problem remains at release" = yes excluded.
##  - IBD_FILE: PLINK --genome; DUPLICATED=Yes & EXPECTED=No rows flagged (not dropped).
DUPS_FILE            <- Sys.getenv("METHYL_DUPS_FILE",            file.path(DATA_DIR, "DUPS.csv"))
PROBLEM_HISTORY_FILE <- Sys.getenv("METHYL_PROBLEM_HISTORY_FILE", file.path(DATA_DIR, "Project_Problem_History.xlsx"))
IBD_FILE             <- Sys.getenv("METHYL_IBD_FILE",             file.path(DATA_DIR, "IBD_results.csv"))

## ---- CATSLife person table ------------------------------------------------
## ADMIN_FILE: the individual-admin .sav (aid/pfamid/nsex/famtype/adopted/LabAge/nidaid/...).
## SAMPLE_LIST_FILE: the per-wave array<->person crosswalk (random_id <-> nidaid + CATSLife
## wave). build_person_table.R merges the two (on nidaid) into CLEAN_ID_FILE, so the sheet's
## PI Provided Subject ID (= random_id) reaches the person world (aid/pfamid).
ADMIN_FILE       <- Sys.getenv("METHYL_ADMIN_FILE",       file.path(DATA_DIR, "individual_admin.sav"))
SAMPLE_LIST_FILE <- Sys.getenv("METHYL_SAMPLE_LIST_FILE", file.path(DATA_DIR, "sample_list.xlsx"))

## ---- Additional-analysis inputs (run_additional_analysis.R) ----------------
## CLEAN_ID_FILE is produced by build_person_table.R.
CLEAN_ID_FILE   <- Sys.getenv("METHYL_CLEAN_ID_FILE",   file.path(DATA_DIR, "CATSLife_pseudo_id.sav"))
DYADS_FILE      <- Sys.getenv("METHYL_DYADS_FILE",      file.path(ANALYSIS_DIR, "catslife_dyads.csv"))
TWIN_PHENO_FILE <- Sys.getenv("METHYL_TWIN_PHENO_FILE", file.path(ANALYSIS_DIR, "synth_twin_pheno.sav"))
TWIN_RS_FILE    <- Sys.getenv("METHYL_TWIN_RS_FILE",    file.path(ANALYSIS_DIR, "synth_twin_rs.xlsx"))

## ---- Vocabulary canonicalization --------------------------------------------
## Canonicalize tissue names ("Buffy Coat" -> "Buffy_Coat") once as the sample sheet
## enters memory, so downstream sees one form. Fails loud on an unrecognized value.
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
## dup-aliquot "D" marker (aliquots are retained per the DUPS consistency check, so both
## members resolve to the same person) and the wave suffix, then require a pure-integer
## residue — fail loud otherwise (e.g. a control) rather than coercing to NA.
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

## ---- Helpers --------------------------------------------------------------
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
## processing sample-subsets and reassembling is exact. Preallocate the Meth/Unmeth matrices and
## fill column-blocks in place (no cbind double-hold), so peak stays ~= one full MethylSet (~25GB
## + the input) instead of the >70GB an all-at-once call needs at ~1600 samples. NOTE: dasen is
## cross-sample, so it is NOT chunkable by sample the way noob is — but it IS streamable; see
## dasen_stream() below for the memory-safe equivalent. Batch via METHYL_NOOB_CHUNK.
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

## Streaming/low-memory dasen — a drop-in for wateRmelon::dasen(MethylSet) whose output is
## IDENTICAL to stock (verified to machine precision, beta max|diff| ~3e-16; guarded by
## test/test_dasen_stream.R) but whose peak is ~40-60GB instead of the >100GB stock needs at
## cohort scale (1642 samples). Stock's memory blows up inside limma::normalizeQuantiles, which
## holds the input submatrix PLUS a full sorted copy PLUS a full rank matrix. Quantile
## normalization needs none of those: build the reference by streaming (pass 1), then map each
## column onto it in place (pass 2). See methylation/troubleshoot/dasen-memory-and-normalization.md.
##
## qn_stream: limma::normalizeQuantiles(A, ties=TRUE), streamed. refcols selects which columns
## define the reference distribution (all columns = cohort average = stock default; a subset =
## wave-1 anchor). Every column is still mapped onto that reference. The apply step is
## rank()+approx() interpolation (ties.method="average"), matching limma's ties=TRUE path exactly
## — naive sorted-position assignment would NOT reproduce stock.
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

## dfsfit: per-sample background offset applied to Type I probes, plus dfsfit's optional
## cross-sample Sentrix row/col (roco) lm smoothing of the per-sample offset scalars. dfs2() gives
## one scalar per sample (streamable); the lm runs over that length-n vector, so it is cheap and
## kept exactly as stock. roco is extracted the same way stock dasen extracts it; the lm is wrapped
## in try() so a degenerate position model (e.g. non-Sentrix colnames) skips smoothing rather than
## erroring — which is what stock effectively does on such data too. wateRmelon:::dfs2 is the same
## internal stock dasen calls, so the background offset stays byte-for-byte faithful.
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

## dasen_stream: dfsfit (roco on Meth, none on Unmeth) then quantile-normalize each channel x
## probe-type. reference = NULL => cohort average (all samples; == stock default). To use a wave-1
## anchor instead, pass
## reference = <the wave-1 columns> (integer indices or a logical mask) — the reference is then
## estimated from those columns while every sample is still mapped onto it. Returns a minfi
## MethylSet; downstream getBeta()/getM()/betas() use offset 100 by default, matching dasen's
## default fudge=100 (this is why the returned MethylSet reproduces stock betas without applying
## fudge here). wateRmelon:::got is the same design-type accessor stock dasen uses.
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
    ## IDATs live one directory per Sentrix barcode (Released_Data/Data/<barcode>/...), NOT
    ## flat; barcode = Sample_Group before the first "_".
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

## ---- Path bridge: see + validate the logical -> physical mapping -----------
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
        row("PHENOTYPE_FILE",  PHENOTYPE_FILE,  "output", "phenotype_bridge"),
        row("PHENOTYPE_FILE",  PHENOTYPE_FILE,  "input",  "stage5"),
        ## ADJUSTED_BETAS_FILE intentionally NOT required: population.R's adjusted pass
        ## skips gracefully when it's absent.
        row("CLEAN_ID_FILE",   CLEAN_ID_FILE,   "input",  "additional_analysis"),
        row("TWIN_PHENO_FILE", TWIN_PHENO_FILE, "input",  "additional_analysis"),
        row("DYADS_FILE",      DYADS_FILE,      "output", "additional_analysis"),
        row("TWIN_RS_FILE",    TWIN_RS_FILE,    "output", "additional_analysis"))
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

## Fail-fast pre-run check. stage in {"stage1","phenotype_bridge","stage5",
## "additional_analysis","all"}: the required inputs must be readable and
## ANALYSIS_DIR writable.
validate_paths <- function(stage = c("all", "stage1", "person_table", "phenotype_bridge", "stage5", "additional_analysis")) {
    stage <- match.arg(stage)
    reg <- .path_registry()
    problems <- character(0)
    dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
    if (file.access(ANALYSIS_DIR, 2) != 0)
        problems <- c(problems, paste("ANALYSIS_DIR not writable:", ANALYSIS_DIR))
    want <- switch(stage,
                   all = c("stage1", "person_table", "phenotype_bridge", "stage5", "additional_analysis", "all"),
                   stage1 = c("stage1", "all"),
                   person_table = c("person_table", "all"),
                   phenotype_bridge = c("phenotype_bridge", "all"),
                   stage5 = c("stage5", "all"),
                   additional_analysis = c("additional_analysis", "all"))
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
            " sample(s) resolved an aid but no FamilyID (pfamid missing in the person table — should never happen): ",
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
            "appear in DYADS_FILE's aid.x/aid.y — the random_id <-> person crosswalk is very likely ",
            "keyed wrong (heritability/twin-corr would silently degrade to no family structure)"))
    if (length(problems))
        stop("validate_phenotype_bridge() failed:\n  - ", paste(problems, collapse = "\n  - "),
             "\nCheck the sample sheet's PI Provided Subject ID (random_id) against SAMPLE_LIST_FILE ",
             "(inspect via scripts/build/inspect_2026_delivery.R's crosswalk dry-run).", call. = FALSE)
    cat("validate_phenotype_bridge(): OK -", nrow(pheno), "samples,", length(have_ids), "individuals\n")
    invisible(TRUE)
}

## ---- Dev-only test profile (smoke tests) — never shipped -------------------
## Sourced LAST so it can override the constants and functions above (e.g. inject a
## synthetic RGChannelSet). Unset in real runs; the profile lives under test/ and is
## excluded from the public runtime.
.test_profile <- Sys.getenv("METHYL_TEST_PROFILE", "")
if (nzchar(.test_profile)) {
    if (!file.exists(.test_profile) && file.exists(file.path(.root, .test_profile)))
        .test_profile <- file.path(.root, .test_profile)
    if (!file.exists(.test_profile))
        stop("METHYL_TEST_PROFILE set but not found: ", .test_profile)
    message("config: loading test profile ", .test_profile)
    source(.test_profile)
}
