#!/usr/bin/env Rscript
### scripts/build/build_phenotype_file.R: build PHENOTYPE_FILE, linking array ids (the
### sheet's PI Provided Subject ID = random_id) to person ids (aid/pfamid) through the
### random_id the person table carries. One row per sample with a Wave column (from the
### "_<n>" suffix); Age is the admin LabAge (wave 2) or LabAge1 (wave 1) for that row. The
### SIF supplies only Father/Mother and a sex cross-check. Label corrections from
### SAMPLE_SWAPS_FILE are applied first. Intentional DUPS pairs are kept and grouped.
### validate_phenotype_bridge() stops the build if a random_id doesn't resolve.
###
### Run after: stage 1 (sample sheet), build_person_table.R, catslife_id_dyads.R, stage 3.
source("config.R"); source("stage5/helpers.R")
suppressMessages({
    library(dplyr); library(tidyr); library(readr); library(readxl); library(haven); library(purrr)
})
validate_paths("phenotype_bridge")

## Two sexes disagree only when both are known, so the empty Sex_geno column never
## raises a flag.
disagree <- function(a, b) !is.na(a) & !is.na(b) & a != b

## 1. Sample sheet: Subject_ID, DNA_Source, drop controls ----
sheet <- read_sample_sheet(SAMPLE_SHEET) %>%
    mutate(DNA_Source = canonicalize_dna_source(DNA_Source),
           Subject_ID = `PI Provided Subject ID`) %>%
    filter(DNA_Source != "Cell_Line", !grepl("^METHYL", Subject_ID))
if (!"SIF_Sex" %in% names(sheet)) sheet$SIF_Sex <- NA_character_

## 1b. Sample-identity corrections (SAMPLE_SWAPS_FILE) ----
## Relabel swapped samples before any identity join. Subject_ID_sheet keeps the sheet's
## label, which the vendor-keyed QC files (DUPS, PROBLEM_HISTORY) still refer to. The
## sheet's SIF_Sex describes the labeled person, so it moves with the corrected label.
swaps <- read_sample_swaps(SAMPLE_SWAPS_FILE)
fix   <- apply_sample_swaps(sheet$Subject_ID, swaps)
sheet <- sheet %>% mutate(Subject_ID_sheet = Subject_ID, Subject_ID = fix$subject_id,
                          identity_action = fix$action)
relabeled <- which(sheet$identity_action %in% "relabel")
sheet$SIF_Sex[relabeled] <- sheet$SIF_Sex[match(sheet$Subject_ID[relabeled], sheet$Subject_ID_sheet)]
pwalk(filter(sheet, !is.na(identity_action)) %>% select(Subject_ID_sheet, Subject_ID, identity_action),
      function(Subject_ID_sheet, Subject_ID, identity_action)
          cat("build_phenotype_file: identity", identity_action, Subject_ID_sheet,
              if (identity_action == "relabel") paste("->", Subject_ID) else "", "\n"))
excluded_subject_ids <- character()

## 2. DUPS: retain both aliquots, tag a shared DupGroupID ----
## Intentional technical replicates, kept so their consistency can be checked. Both
## members (canonical and "_D" re-run) get DupGroupID = the canonical id.
dup_group <- if (file.exists(DUPS_FILE)) {
    dups <- read_csv(DUPS_FILE, col_types = cols(.default = "c"))
    tibble(Subject_ID = c(dups[["Subject ID 1"]], dups[["Subject ID 2"]]),
           DupGroupID = c(dups[["Subject ID 1"]], dups[["Subject ID 1"]])) %>%
        distinct(Subject_ID, .keep_all = TRUE)
} else tibble(Subject_ID = character(), DupGroupID = character())
sheet <- left_join(sheet, dup_group, by = c("Subject_ID_sheet" = "Subject_ID"))
cat("build_phenotype_file: retained", sum(!is.na(sheet$DupGroupID)),
    "intentional-duplicate sample(s) in", n_distinct(sheet$DupGroupID, na.rm = TRUE), "group(s)\n")

## 3. PROBLEM_HISTORY: exclude subjects with a known unresolved issue ----
if (file.exists(PROBLEM_HISTORY_FILE)) {
    flagged <- read_excel(PROBLEM_HISTORY_FILE) %>%
        filter(tolower(trimws(`Does a problem remain at release?`)) %in% c("yes", "y"))
    if (nrow(flagged)) {
        pwalk(list(flagged[["Subject ID"]], flagged[["Problem Description"]]),
              ~ cat("build_phenotype_file: excluding Subject_ID", .x, "-", .y, "\n"))
        drop <- sheet$Subject_ID_sheet %in% flagged[["Subject ID"]]
        excluded_subject_ids <- c(excluded_subject_ids, sheet$Subject_ID[drop])
        sheet <- sheet[!drop, ]
    }
}

## Person table (used by the IBD MZ check and the crosswalk below) ----
person <- read_sav(CLEAN_ID_FILE) %>%
    mutate(aid = as.integer(aid), random_id = as.integer(random_id)) %>%
    filter(!is.na(random_id))                       # only sampled persons carry a random_id

## Samples whose random_id is in the sample list but has no confirmed nidaid yet ("pending",
## e.g. a newly assigned id) have no person to join; leave them out, by name, until the sample
## list is completed. Any other random_id that fails to resolve still stops the build below.
sl_pending <- read_sample_list() %>% filter(is.na(nidaid)) %>% pull(random_id)
sl_pending <- setdiff(sl_pending, person$random_id)
pend <- subject_base_id(sheet$Subject_ID) %in% sl_pending
walk(sheet$Subject_ID[pend], ~ cat("build_phenotype_file: excluding Subject_ID", .x,
                                   "- its random_id has no confirmed nidaid in SAMPLE_LIST_FILE\n"))
excluded_subject_ids <- c(excluded_subject_ids, sheet$Subject_ID[pend])
sheet <- sheet[!pend, ]

## 4. IBD: flag cross-wave resamples and unexpected duplicates ----
## DUPLICATED=Yes & EXPECTED=No = a genetic duplicate not in DUPS_FILE. classify_ibd_pair()
## (config.R) splits these against the person table: "cross_wave" = one participant resampled
## across waves (shared LongitudinalGroupID); "mz" = an MZ co-twin pair (same pfamid,
## ZygGroup=1), expected; "unexpected" = a likely swap/mislabel, reported for review.
## No sample is dropped here.
long_group <- tibble(Subject_ID = character(), LongitudinalGroupID = character())
if (file.exists(IBD_FILE)) {
    id_map <- select(sheet, Sample_ID, Subject_ID)
    ibd <- read_csv(IBD_FILE, show_col_types = FALSE) %>%
        filter(DUPLICATED == "Yes", EXPECTED == "No") %>%
        left_join(id_map, by = c("IID1" = "Sample_ID")) %>% rename(s1 = Subject_ID) %>%
        left_join(id_map, by = c("IID2" = "Sample_ID")) %>% rename(s2 = Subject_ID) %>%
        filter(!is.na(s1), !is.na(s2)) %>%
        mutate(pair_class = classify_ibd_pair(s1, s2, person))
    pwalk(filter(ibd, pair_class == "cross_wave"), function(s1, s2, PI_HAT, ...)
        cat("build_phenotype_file: flagged", s1, "<->", s2, "as a likely cross-wave resample (PI_HAT=", PI_HAT, ")\n"))
    pwalk(filter(ibd, pair_class == "mz"), function(s1, s2, PI_HAT, ...)
        cat("build_phenotype_file: flagged", s1, "<->", s2, "as an expected MZ co-twin pair (same pfamid, ZygGroup=1, PI_HAT=", PI_HAT, ")\n"))
    pwalk(filter(ibd, pair_class == "unexpected"), function(s1, s2, ...)
        cat("build_phenotype_file: WARNING - unexpected duplicate", s1, "<->", s2,
            "is neither a cross-wave resample nor an MZ co-twin pair - needs manual review\n"))
    pairs <- filter(ibd, pair_class == "cross_wave")
    long_group <- bind_rows(
        transmute(pairs, Subject_ID = s1, LongitudinalGroupID = strip_wave_suffix(s1)),
        transmute(pairs, Subject_ID = s2, LongitudinalGroupID = strip_wave_suffix(s1))) %>%
        distinct(Subject_ID, .keep_all = TRUE)
}
sheet <- left_join(sheet, long_group, by = "Subject_ID")

## 5. Crosswalk: sheet random_id + Wave -> person table ----
## subject_base_id() gives the array-facing random_id and subject_wave() the wave; the
## person table carries random_id, so one join brings in identity (aid), family (pfamid),
## sex (nsex) and both waves' ages.
bridge <- sheet %>%
    mutate(random_id = subject_base_id(Subject_ID), Wave = subject_wave(Subject_ID)) %>%
    left_join(select(person, random_id, aid, pfamid, nsex, age, age_w1, famtype),
              by = "random_id", relationship = "many-to-one") %>%
    mutate(IndividualID = aid, FamilyID = pfamid)

## 6. SIF: Father/Mother + pedigree sex ----
## SIF Individual is the vendor Subject_ID with underscores removed, not the aid, so the
## SIF is not used for identity. Founder/parent rows carry Subject_ID = NA; drop them.
sif <- read_excel(ID_KEY) %>%
    filter(!is.na(Subject_ID)) %>%
    mutate(Subject_ID = as.character(Subject_ID)) %>%
    select(Subject_ID, Father, Mother, Sex_ped_num = Sex)
bridge <- left_join(bridge, sif, by = "Subject_ID", relationship = "many-to-one")
bad_sex <- setdiff(unique(bridge$Sex_ped_num), c(1, 2, NA))
if (length(bad_sex))
    stop("build_phenotype_file: SIF Sex value(s) outside {1,2,NA} PED coding: ",
         paste(bad_sex, collapse = ", "))

## 7. Sex QC: genotype vs admin (primary) vs pedigree ----
## nsex is the primary source; Sex_geno is reserved for a genotype-based call and is NA
## for now. Disagreements are flagged, not dropped, and written to sex_qc.csv together
## with the sheet's SIF_Sex as a further cross-check.
bridge <- bridge %>% mutate(
    Sex       = case_when(nsex == 1 ~ "M", nsex == 0 ~ "F", TRUE ~ NA_character_),
    Sex_ped   = case_when(Sex_ped_num == 1 ~ "M", Sex_ped_num == 2 ~ "F", TRUE ~ NA_character_),
    Sex_geno  = NA_character_,
    Sex_sheet = if_else(SIF_Sex %in% c("M", "F"), SIF_Sex, NA_character_),
    Sex_flag  = disagree(Sex, Sex_ped) | disagree(Sex, Sex_geno))
bridge %>%
    transmute(Sample = Sample_Group, Subject_ID, IndividualID,
              Sex_admin = Sex, Sex_ped, Sex_geno, Sex_sheet, Sex_flag) %>%
    write_csv(file.path(REPORT_DIR, "sex_qc.csv"))
cat("build_phenotype_file: sex QC -", sum(bridge$Sex_flag, na.rm = TRUE),
    "sample(s) with an admin-vs-pedigree/genotype disagreement (flagged, not dropped); wrote sex_qc.csv\n")

## 8. DUPS consistency: paired aliquots must resolve to one person/sex ----
dup_check <- bridge %>%
    filter(!is.na(DupGroupID)) %>%
    group_by(DupGroupID) %>%
    summarise(discordant = n_distinct(IndividualID) > 1 | n_distinct(Sex) > 1, .groups = "drop")
bridge <- bridge %>%
    left_join(dup_check, by = "DupGroupID") %>%
    mutate(Dup_flag = coalesce(discordant, FALSE)) %>%
    select(-discordant)
cat("build_phenotype_file: DUPS consistency -", sum(dup_check$discordant),
    "of", nrow(dup_check), "duplicate group(s) discordant on person/sex\n")

## 9. FamilyType/Zygosity from the dyad table (per-aid representative) ----
dyads <- read_csv(DYADS_FILE, show_col_types = FALSE)
aid_famtype <- dyads %>%
    pivot_longer(c(aid.x, aid.y), values_to = "aid") %>%
    mutate(aid = as.integer(aid)) %>%
    distinct(aid, .keep_all = TRUE) %>%
    select(aid, FamilyType = famtype_pair)
bridge <- bridge %>%
    left_join(aid_famtype, by = c("IndividualID" = "aid"), relationship = "many-to-one") %>%
    mutate(Zygosity = zygosity_from_famtype(FamilyType))

## 10. Cell proportions (stage 3 output) ----
if (file.exists(CELL_PROPORTIONS_FILE)) {
    non_cell <- c("Sample_Group", "DNA Source", "Population", "SIF_Sex", "Core_Lab_ID",
                  "Study_ID", "Family_Relationship", "Family", "DNA_Source", "Array")
    cellprop <- read_tsv(CELL_PROPORTIONS_FILE, show_col_types = FALSE)
    cell_cols <- setdiff(names(cellprop), non_cell)
    cellprop <- cellprop %>%
        rename_with(~ paste0("cell_", .x), all_of(cell_cols)) %>%
        select(Sample_Group, starts_with("cell_"))
    bridge <- left_join(bridge, cellprop, by = "Sample_Group", relationship = "many-to-one")
} else cat("build_phenotype_file: CELL_PROPORTIONS_FILE not found; no cell covariates\n")

## 11. Assemble, validate, write ----
out <- bridge %>% transmute(
    Sample = Sample_Group, Subject_ID, IndividualID, FamilyID, Wave, FamilyType, Zygosity,
    DNA_Source, Age = if_else(Wave == 2, age, age_w1),   # LabAge (wave-2) / LabAge1 (wave-1)
    Sex, Sex_geno, Sex_flag, Identity_flag = identity_action %in% "flag", Subject_ID_sheet,
    DupGroupID, Dup_flag, Sample_Plate,
    LongitudinalGroupID = coalesce(LongitudinalGroupID, Subject_ID))
cell_out <- select(bridge, starts_with("cell_"))
if (ncol(cell_out)) out <- bind_cols(out, cell_out)

validate_phenotype_bridge(out, dyads, excluded_subject_ids = excluded_subject_ids)
write_csv(out, PHENOTYPE_FILE)
cat("build_phenotype_file: wrote", PHENOTYPE_FILE, "-", nrow(out), "samples,",
    n_distinct(out$IndividualID, na.rm = TRUE), "individuals\n")
