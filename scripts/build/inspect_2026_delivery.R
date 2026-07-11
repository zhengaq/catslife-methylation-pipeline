#!/usr/bin/env Rscript
### Usage: Rscript scripts/build/inspect_2026_delivery.R [--full-rows]
###   --full-rows   also print a few example DATA rows — only use against
###                 de-identified examples, never participant files
suppressMessages({ library(readxl); library(haven) })
source("config.R")

args      <- commandArgs(trailingOnly = TRUE)
show_rows <- "--full-rows" %in% args
n_rows    <- if (show_rows) 3L else 0L

root <- Sys.getenv("METHYL_2026_RECON_DIR", "data/methylation_2026")
p <- function(...) file.path(root, ...)

section <- function(title) cat("\n====", title, "====\n")

## ---- 1. Sample sheet (GenomeStudio bracket-preamble format) --------------
section("Sample sheet")
sheet_path <- Sys.getenv("METHYL_RECON_SAMPLE_SHEET", p("example_SampleSheet_2_EPIC.csv"))
if (file.exists(sheet_path)) {
    raw <- readLines(sheet_path, n = 15L, warn = FALSE)
    cat("First 15 raw lines (preamble check):\n"); cat(paste0("  ", raw), sep = "\n")
    d <- read_sample_sheet(sheet_path)
    cat("\nParsed: ", nrow(d), "rows x", ncol(d), "columns\n")
    cat("Columns:", paste(names(d), collapse = " | "), "\n")
    if (show_rows) print(utils::head(d, n_rows))
} else cat("(not found:", sheet_path, ")\n")

## ---- 2. Pedigree / SIF file (xlsx) ----------------------------------------
section("Pedigree / SIF file")
sif_path <- Sys.getenv("METHYL_RECON_SIF", p("example_SIF.xlsx"))
if (file.exists(sif_path)) {
    for (s in excel_sheets(sif_path)) {
        d <- read_excel(sif_path, sheet = s)
        cat("Sheet '", s, "': ", nrow(d), " rows x ", ncol(d), " columns\n", sep = "")
        cat("Columns:", paste(names(d), collapse = " | "), "\n")
        if (show_rows) print(as.data.frame(utils::head(d, n_rows)))
    }
} else cat("(not found:", sif_path, ")\n")

## ---- 2b. CATSLife individual-admin file (person source-of-truth) ----------
## nidaid is the join key to the sample list (-> CLEAN_ID_FILE).
section("Individual-admin file (.sav)")
admin_path <- Sys.getenv("METHYL_RECON_ADMIN", p("Example_individual_admin.sav"))
adm <- NULL
if (file.exists(admin_path)) {
    adm <- read_sav(admin_path)
    cat(nrow(adm), "rows x", ncol(adm), "columns\n")
    cat("Columns:", paste(names(adm), collapse = " | "), "\n")
    if ("aid" %in% names(adm)) cat("aid: ", sum(is.na(adm$aid)), " NA, ", length(unique(adm$aid)), " unique\n", sep = "")
    if ("nidaid" %in% names(adm)) cat("nidaid (join key): ", length(unique(adm$nidaid)), " unique\n", sep = "")
    if ("LabAge" %in% names(adm)) cat("LabAge (wave-2 age): ", sum(is.na(adm$LabAge)), " NA\n", sep = "")
    if ("LabAge1" %in% names(adm)) cat("LabAge1 (wave-1 age): ", sum(is.na(adm$LabAge1)), " NA\n", sep = "")
    if (show_rows) print(as.data.frame(utils::head(adm, n_rows)))
} else cat("(not found:", admin_path, ")\n")

## ---- 2c. Sample list: the array<->person crosswalk (random_id <-> nidaid + wave) --
section("Sample list (random_id <-> nidaid + wave)")
sl_path <- Sys.getenv("METHYL_RECON_SAMPLE_LIST", p("Buffy Coat DNA Methylation Sample List.xlsx"))
sl <- NULL
if (file.exists(sl_path)) {
    sl <- read_excel(sl_path)
    cat(nrow(sl), "rows x", ncol(sl), "columns\n")
    cat("Columns:", paste(names(sl), collapse = " | "), "\n")
    if ("random_id" %in% names(sl)) cat("random_id: ", length(unique(sl$random_id)), " unique\n", sep = "")
    if ("nidaid" %in% names(sl)) cat("nidaid: ", length(unique(sl$nidaid)), " unique\n", sep = "")
    if ("CATSLife" %in% names(sl)) { cat("wave (CATSLife):\n"); print(table(sl$CATSLife)) }
    if (show_rows) print(as.data.frame(utils::head(sl, n_rows)))
} else cat("(not found:", sl_path, ")\n")

## ---- 2d. Crosswalk dry-run: sheet random_id -> sample list -> nidaid -> admin --
## Previews what validate_phenotype_bridge() enforces, per hop, so a low rate
## localizes the break.
section("Crosswalk dry-run (sheet random_id -> [sample list] nidaid -> [admin] person)")
if (file.exists(sheet_path)) {
    subj <- read_sample_sheet(sheet_path)[["PI Provided Subject ID"]]
    subj <- subj[!is.na(subj) & !grepl("^METHYL", subj)]                  # drop controls
    rid  <- suppressWarnings(as.integer(strip_wave_suffix(sub("_[0-9]*D$", "", subj))))
    cat(sprintf("%d non-control sheet sample(s); %d are wave 2 (\"_2\" suffix)\n",
                length(subj), sum(subject_wave(subj) == 2)))
    if (!is.null(sl) && "random_id" %in% names(sl)) {
        in_sl <- rid %in% suppressWarnings(as.integer(sl$random_id))
        cat(sprintf("  hop 1  sheet random_id in the sample list: %d/%d (%.1f%%)\n",
                    sum(in_sl, na.rm = TRUE), length(rid), 100 * sum(in_sl, na.rm = TRUE) / length(rid)))
    } else cat("  (hop 1 skipped: sample list not available)\n")
    if (!is.null(sl) && !is.null(adm) && "nidaid" %in% names(sl) && "nidaid" %in% names(adm)) {
        in_adm <- as.character(sl$nidaid) %in% as.character(adm$nidaid)
        cat(sprintf("  hop 2  sample-list nidaid in the admin file: %d/%d (%.1f%%)\n",
                    sum(in_adm), length(in_adm), 100 * sum(in_adm) / length(in_adm)))
    } else cat("  (hop 2 skipped: sample list or admin not available)\n")
} else cat("(skipped: sample sheet not available)\n")

## ---- 3. Samples_Table (QC/signal summary, NOT an ID bridge) --------------
section("Samples_Table (QC summary)")
st_path <- Sys.getenv("METHYL_RECON_SAMPLES_TABLE", p("example_Samples_Table.csv"))
if (file.exists(st_path)) {
    d <- read.csv(st_path, check.names = FALSE, nrows = if (show_rows) n_rows else 0)
    cat("Columns:", paste(names(d), collapse = " | "), "\n")
    if (show_rows) print(d)
} else cat("(not found:", st_path, ")\n")

## ---- 4. DUPS file (intentional QC duplicate pairs) -----------------------
section("DUPS file")
dups_path <- Sys.getenv("METHYL_RECON_DUPS", p("example_Sample_DUPS.csv"))
if (file.exists(dups_path)) {
    d <- read.csv(dups_path, check.names = FALSE)
    cat(nrow(d), "row(s). Columns:", paste(names(d), collapse = " | "), "\n")
    if (show_rows) print(utils::head(d, n_rows))
} else cat("(not found:", dups_path, ")\n")

## ---- 5. IBD relatedness file (PLINK --genome output) ----------------------
section("IBD relatedness file")
ibd_path <- Sys.getenv("METHYL_RECON_IBD", p("example_IBD_unexp_dups_exp.csv"))
if (file.exists(ibd_path)) {
    d <- read.csv(ibd_path, check.names = FALSE)
    cat(nrow(d), "row(s). Columns:", paste(names(d), collapse = " | "), "\n")
    if ("DUPLICATED" %in% names(d) && "EXPECTED" %in% names(d))
        print(table(DUPLICATED = d$DUPLICATED, EXPECTED = d$EXPECTED))
} else cat("(not found:", ibd_path, ")\n")

## ---- 6. Project problem history --------------------------------------------
section("Project problem history")
prob_path <- Sys.getenv("METHYL_RECON_PROBLEM_HISTORY", p("example_project_problem.xlsx"))
if (file.exists(prob_path)) {
    for (s in excel_sheets(prob_path)) {
        d <- read_excel(prob_path, sheet = s)
        cat("Sheet '", s, "': ", nrow(d), " rows x ", ncol(d), " columns\n", sep = "")
        cat("Columns:", paste(names(d), collapse = " | "), "\n")
    }
} else cat("(not found:", prob_path, ")\n")

## ---- 7. Methylation_Profile.txt (wide GenomeStudio Final Report) ---------
## Header only — the file can be tens of thousands of columns wide, so never load it wholesale.
section("Methylation_Profile.txt (header only, never loaded wholesale)")
profile_path <- Sys.getenv("METHYL_RECON_PROFILE", p("Methylation_Profile_headers_rows.txt"))
if (file.exists(profile_path)) {
    header <- strsplit(readLines(profile_path, n = 1L, warn = FALSE), "\t")[[1]]
    cat("Total columns:", length(header), "\n")
    cat("First ~50 annotation columns (stops before the per-sample block):\n")
    cat(paste0("  ", utils::head(header, 50)), sep = "\n")
} else cat("(not found:", profile_path, ")\n")

## ---- 8. Sentrix barcode subdirectories vs sample sheet's Sentrix_ID ------
## Self-skips (rather than false-matching unrelated subdirectories) unless IDAT_DIR
## looks like a Released_Data/Data/ root with digit-named barcode subdirs.
section("Barcode subdirectory cross-check (IDAT_DIR)")
if (dir.exists(IDAT_DIR) && file.exists(sheet_path)) {
    all_dirs <- list.dirs(IDAT_DIR, full.names = FALSE, recursive = FALSE)
    have_dirs <- grep("^[0-9]+$", all_dirs, value = TRUE)  # only digit-only (barcode-shaped) names
    sheet_barcodes <- unique(sub("_.*$", "", read_sample_sheet(sheet_path)$Sample_Group))
    if (!length(have_dirs)) {
        cat("(skipped: no digit-named subdirectories under IDAT_DIR (", IDAT_DIR, ") — ",
            "point METHYL_IDAT_DIR at Released_Data/Data/ to run this check against a delivery)\n", sep = "")
    } else {
        missing <- setdiff(sheet_barcodes, have_dirs)
        cat(length(have_dirs), "barcode dir(s) under IDAT_DIR;", length(sheet_barcodes),
            "distinct barcode(s) in the sample sheet.\n")
        if (length(missing)) cat("MISSING subdirectories for:", paste(missing, collapse = ", "), "\n")
        else cat("All sample-sheet barcodes have a matching subdirectory.\n")
    }
} else cat("(skipped: IDAT_DIR or sample sheet not available)\n")

cat("\nDone.\n")
