#!/usr/bin/env Rscript
### Form sibling dyads from the cleaned CATSLife ID table and assign the pair-level
### relationship `famtype_pair`. Runs on CLEAN_ID_FILE. For mixed CAP households the
### pair type depends on both members' adoption status.
### famtype_pair: 1=adoptive, 2=biological sibling, 3=DZ, 4=MZ.
source("config.R")
suppressMessages({ library(dplyr); library(haven) })

CATSLife_ID <- read_sav(CLEAN_ID_FILE)

CATSLife_ID_wide <- CATSLife_ID %>%
  left_join(CATSLife_ID,
            by = join_by(pfamid, famtype, household_composition, project),
            suffix = c(".x", ".y"), relationship = "many-to-many") %>%
  filter(aid.x < aid.y) %>%
  mutate(
    # Retain famtype (TRUE ~ famtype) except for the special CAP-sibling recodes;
    # the leading guard keeps DZ/MZ twins from being collapsed to "biological sibling".
    famtype_pair = case_when(
      famtype %in% c(3, 4)            ~ as.numeric(famtype),  # twins retain zygosity (3=DZ, 4=MZ)
      adopted.x == 0 & adopted.y == 0 ~ 2,                     # non-twin biological siblings
      adopted.x == 1 | adopted.y == 1 ~ 1,                     # adoptive relationship
      TRUE                            ~ as.numeric(famtype)
    )
  )

dir.create(dirname(DYADS_FILE), recursive = TRUE, showWarnings = FALSE)
write.csv(CATSLife_ID_wide, DYADS_FILE, row.names = FALSE)
cat("Wrote", DYADS_FILE, "(", nrow(CATSLife_ID_wide), "dyads )\n")
print(table(project = CATSLife_ID_wide$project, famtype_pair = CATSLife_ID_wide$famtype_pair))