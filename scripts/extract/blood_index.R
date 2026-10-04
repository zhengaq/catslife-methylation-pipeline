#!/usr/bin/env Rscript
### scripts/extract/blood_index.R: bridge a list of Buffy Coat lab IDs to array positions
### (Sample_Group) and extract their cell-type proportions. Not part of the pipeline run:
### it reads a completed run's outputs and writes only to [out dir].
###
### lab ID = SampleLog on SAMPLE_LIST_FILE -> nidaid -> all of that person's waves -> array,
### joined on PHENOTYPE_FILE's swap-corrected Subject_ID so a swapped array goes to the person
### whose DNA it holds. One row per (person, wave, array); rows that cannot be bridged are kept
### with the reason in `status`.
###
### Usage: Rscript scripts/extract/blood_index.R <lab-ID file> [out dir]
###   <lab-ID file>  one column headed "lab ID"
###   [out dir]      default: the lab-ID file's directory
### Writes blood_index.csv and cell_proportions.blood_subset.csv.
source("config.R")
suppressMessages({ library(dplyr); library(readr); library(readxl) })

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript scripts/extract/blood_index.R <lab-ID file> [out dir]")
request_file <- args[1]
out_dir      <- if (length(args) > 1) args[2] else dirname(request_file)

CELL_TYPES <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")

## A copy cut off mid-row would silently lose samples.
for (f in c(request_file, SAMPLE_LIST_FILE, SAMPLE_SWAPS_FILE, PHENOTYPE_FILE, CELL_PROPORTIONS_FILE))
    if (!file.exists(f)) stop("blood_index: missing ", f)
for (f in c(PHENOTYPE_FILE, CELL_PROPORTIONS_FILE)) {
    con <- file(f, "rb"); seek(con, file.size(f) - 1); last <- readBin(con, "raw", 1); close(con)
    if (last != as.raw(10)) stop("blood_index: ", f, " does not end in a newline; the copy is truncated")
}

req <- read_tsv(request_file, col_types = cols(`lab ID` = "i"))
if (anyDuplicated(req$`lab ID`))
    stop("blood_index: duplicate lab IDs in ", request_file, ": ",
         paste(unique(req$`lab ID`[duplicated(req$`lab ID`)]), collapse = ", "))

## The list is edited by hand: ids may be typed as text, so read everything as text.
sl <- read_excel(SAMPLE_LIST_FILE, col_types = "text", na = c("", "NA"), trim_ws = TRUE)
need <- c("nidaid", "SampleLog", "CATSLife", "random_id")
if (!all(need %in% names(sl)))
    stop("blood_index: ", SAMPLE_LIST_FILE, " needs columns ", paste(need, collapse = ", "))
whole <- function(x, what) {
    if (any(bad <- !is.na(x) & !grepl("^[0-9]+(\\.0+)?$", x)))
        stop("blood_index: ", what, " value(s) that are not whole numbers: ", paste(unique(x[bad]), collapse = ", "))
    as.integer(as.numeric(x))
}
sl <- sl %>%
    filter(!is.na(SampleLog)) %>%
    transmute(nidaid, `lab ID` = whole(SampleLog, "SampleLog"), random_id = whole(random_id, "random_id"),
              wave = sub("^CATSLife ([12])$", "CATSLife\\1", CATSLife))
bad <- setdiff(unique(sl$wave), c("CATSLife1", "CATSLife2"))
if (length(bad)) stop("blood_index: unrecognized CATSLife value(s): ", paste(bad, collapse = ", "))
if (anyDuplicated(sl$`lab ID`))
    stop("blood_index: duplicate SampleLog(s) in ", SAMPLE_LIST_FILE, ": ",
         paste(unique(sl$`lab ID`[duplicated(sl$`lab ID`)]), collapse = ", "))

hit      <- left_join(req, sl, by = "lab ID", relationship = "one-to-one")
unlisted <- hit %>% filter(is.na(random_id)) %>% mutate(status = "lab ID not on the sample list")
no_nid   <- hit %>% filter(!is.na(random_id), is.na(nidaid)) %>%
    mutate(status = "no nidaid on the sample list")

people <- unique(na.omit(hit$nidaid))
rows <- sl %>% filter(nidaid %in% people) %>%
    mutate(requested = `lab ID` %in% req$`lab ID`) %>%
    group_by(nidaid) %>%
    mutate(wave_coverage = if (n_distinct(wave) == 2) "both" else paste(first(wave), "only")) %>%
    ungroup()

swaps <- read_sample_swaps(SAMPLE_SWAPS_FILE)
ph <- read_csv(PHENOTYPE_FILE, col_types = cols(.default = "c")) %>%
    filter(DNA_Source == "Buffy_Coat")
relabeled <- ph$Subject_ID != ph$Subject_ID_sheet
unknown <- setdiff(ph$Subject_ID_sheet[relabeled], swaps$from[swaps$action == "relabel"])
if (length(unknown))
    stop("blood_index: PHENOTYPE_FILE relabels ", paste(unknown, collapse = ", "),
         ", which SAMPLE_SWAPS_FILE does not list; the phenotype file was built from another swaps file")
ph <- ph %>%
    transmute(random_id = subject_base_id(Subject_ID), wave = paste0("CATSLife", Wave),
              Sample_Group = Sample, Subject_ID_sheet,
              swap = ifelse(relabeled, swaps$notes[match(Subject_ID_sheet, swaps$from)], NA),
              aliquot = grepl("D$", Subject_ID), identity_flag = Identity_flag == "TRUE")

cells <- read_tsv(CELL_PROPORTIONS_FILE, show_col_types = FALSE)

status_of <- function(sample_group, swap, from, identity_flag, aliquot, in_cells) {
    if (is.na(sample_group)) return("no array in phenotype file")
    notes <- c(if (!is.na(swap)) paste0("relabeled from ", from, " (", swap, ")"),
               if (identity_flag) "identity flag: fingerprint discrepancy",
               if (aliquot) "technical duplicate aliquot",
               if (!in_cells) "no cell-type proportions for this array")
    if (length(notes)) paste(notes, collapse = "; ") else "ok"
}

out <- rows %>%
    left_join(ph, by = c("random_id", "wave"), relationship = "one-to-many") %>%
    mutate(status = mapply(status_of, Sample_Group, swap, Subject_ID_sheet, identity_flag, aliquot,
                           Sample_Group %in% cells$Sample_Group)) %>%
    bind_rows(mutate(no_nid, requested = TRUE), mutate(unlisted, requested = TRUE)) %>%
    arrange(nidaid, wave, Sample_Group) %>%
    select(nidaid, `lab ID`, wave, Sample_Group, random_id, requested, wave_coverage, status)

sub <- out %>%
    filter(!is.na(Sample_Group)) %>%
    inner_join(select(cells, Sample_Group, all_of(CELL_TYPES)), by = "Sample_Group") %>%
    select(nidaid, `lab ID`, wave, Sample_Group, wave_coverage, status, all_of(CELL_TYPES))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(out, file.path(out_dir, "blood_index.csv"), na = "")
write_csv(sub, file.path(out_dir, "cell_proportions.blood_subset.csv"), na = "")

cat("blood_index:", nrow(req), "lab IDs,", length(people), "persons ->", nrow(out), "index rows,",
    nrow(sub), "cell-proportion rows in", out_dir, "\n")
print(as.data.frame(count(out, wave, status)), row.names = FALSE)
