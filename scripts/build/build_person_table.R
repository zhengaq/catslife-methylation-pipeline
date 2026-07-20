#!/usr/bin/env Rscript
### scripts/build/build_person_table.R — derive the CATSLife person table (CLEAN_ID_FILE) by
### merging the individual-admin .sav (ADMIN_FILE) with the array<->person crosswalk
### (SAMPLE_LIST_FILE: random_id <-> nidaid) on nidaid. The result carries random_id (the
### sheet's join key) alongside aid/pfamid, so build_phenotype_file.R can reach the person
### world in one hop. Feeds catslife_id_dyads.R and build_phenotype_file.R. Run FIRST.
source("config.R")
suppressMessages({ library(dplyr); library(haven); library(readxl) })
validate_paths("person_table")

admin <- read_sav(ADMIN_FILE)
if (!"ZygGroup" %in% names(admin)) admin$ZygGroup <- NA_real_   # zygosity index (1=MZ); tolerate admins that lack it

## Sample list = the array<->person crosswalk; dedupe to one (nidaid, random_id) per person.
sample_list <- read_excel(SAMPLE_LIST_FILE) %>%
    transmute(nidaid = as.character(nidaid), random_id = as.integer(random_id)) %>%
    distinct(nidaid, random_id)

## Every sampled nidaid must be present in the admin file — fail loud, not silent NA.
missing <- setdiff(sample_list$nidaid, as.character(admin$nidaid))
if (length(missing))
    stop("build_person_table: ", length(missing), " sample-list nidaid(s) absent from the admin file: ",
         paste(utils::head(missing, 10), collapse = ", "))

## nsex (0=F,1=M) is the authoritative self-reported sex. Fail loud on an out-of-range
## code rather than coercing to NA, so a malformed admin file stops here.
bad_nsex <- setdiff(unique(admin$nsex), c(0, 1, NA))
if (length(bad_nsex))
    stop("build_person_table: nsex value(s) outside the expected {0,1,NA} coding: ",
         paste(bad_nsex, collapse = ", "))

person <- admin %>%
    transmute(
        project,
        ## aid is NA-free & unique; the numeric ID column is NA for un-sampled persons,
        ## so key everything on as.integer(aid), never ID.
        aid       = as.integer(aid),
        nidaid    = as.character(nidaid),
        pfamid    = as.integer(pfamid),
        age       = LabAge,              # wave-2 age (also the person-level age for stage 5)
        age_w1    = LabAge1,             # wave-1 age (mapped to wave-1 rows in the bridge)
        nsex,                            # for the bridge's sex cross-check
        female    = case_when(nsex == 0 ~ 1, nsex == 1 ~ 0, TRUE ~ NA_real_),
        white     = case_when(racecat_ORIG == 5 ~ 1, racecat_ORIG %in% c(9, NA) ~ NA_real_, TRUE ~ 0),
        hispanic  = case_when(hispanic_ORIG == "Y" ~ 1, hispanic_ORIG == "N" ~ 0, TRUE ~ NA_real_),
        adopted,
        famtype,
        ZygGroup  = as.integer(ZygGroup)   # zygosity index (1=MZ); auto-classifies IBD duplicates
    ) %>%
    ## attach random_id (NA for un-sampled persons — kept so dyads retain full families)
    left_join(sample_list, by = "nidaid", relationship = "many-to-one") %>%
    group_by(pfamid) %>%
    mutate(household_composition = case_when(
        all(adopted == 0, na.rm = TRUE) ~ "all_biological",
        all(adopted == 1, na.rm = TRUE) ~ "all_adopted",
        any(adopted == 0, na.rm = TRUE) & any(adopted == 1, na.rm = TRUE) ~ "mixed_bio_adoptive",
        TRUE ~ NA_character_)) %>%
    ungroup()

if (anyNA(person$aid) || anyDuplicated(person$aid))
    stop("build_person_table: aid must be unique and non-NA after as.integer(aid)")
if (anyNA(person$pfamid))
    stop("build_person_table: pfamid has NA(s) — every person must have a family id")

write_sav(person, CLEAN_ID_FILE)
cat("build_person_table: wrote", CLEAN_ID_FILE, "-", nrow(person), "persons,",
    sum(!is.na(person$random_id)), "with a random_id,", length(unique(person$pfamid)), "families\n")
