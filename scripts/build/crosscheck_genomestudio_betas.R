#!/usr/bin/env Rscript
### scripts/build/crosscheck_genomestudio_betas.R — optional sanity check (not part
### of the normal pipeline): compares this pipeline's minfi->noob->dasen betas
### (F_DASENB) against GenomeStudio's Methylation_Profile.txt AVG_Beta values, for a
### small probes x samples sample.
###
suppressMessages({ library(data.table) })
source("config.R")

N_SAMPLES <- as.integer(Sys.getenv("METHYL_CROSSCHECK_N_SAMPLES", "5"))
N_PROBES  <- as.integer(Sys.getenv("METHYL_CROSSCHECK_N_PROBES",  "200"))
profile_path <- Sys.getenv("METHYL_PROFILE_FILE",
                           file.path(DATA_DIR, "methylation_2026", "Methylation_Profile_headers_rows.txt"))

if (!file.exists(F_DASENB))
    stop("crosscheck_genomestudio_betas: ", F_DASENB, " not found — run stage 1 first.")
if (!file.exists(profile_path))
    stop("crosscheck_genomestudio_betas: Methylation_Profile.txt not found (", profile_path,
         ") — set METHYL_PROFILE_FILE.")

dasen <- load_one(F_DASENB)
betas <- canonicalize_v2_probe_ids(dasen$b)   # bare cg-id rownames, matching Methylation_Profile's NAME column

## ---- Map our Sample_Group (betas' colnames) -> Sample_Name (profile's column prefix) ----
sheet <- read_sample_sheet(SAMPLE_SHEET)
if (!all(c("Sample_Group", "Sample_Name") %in% names(sheet)))
    stop("crosscheck_genomestudio_betas: SAMPLE_SHEET has no Sample_Name column — ",
         "can't map to Methylation_Profile.txt's per-sample column prefixes.")
group_to_name <- setNames(sheet$Sample_Name, sheet$Sample_Group)

header <- strsplit(readLines(profile_path, n = 1L, warn = FALSE), "\t")[[1]]
if (!"NAME" %in% header) stop("crosscheck_genomestudio_betas: Methylation_Profile.txt has no NAME column.")

our_samples   <- intersect(colnames(betas), names(group_to_name))
sample_names  <- unname(group_to_name[our_samples])
avg_beta_cols <- paste0(sample_names, ".AVG_Beta")
have_cols     <- intersect(avg_beta_cols, header)
if (!length(have_cols))
    stop("crosscheck_genomestudio_betas: none of this cohort's samples have a matching ",
         "<Sample_Name>.AVG_Beta column in Methylation_Profile.txt — wrong file/wave?")
use_cols <- utils::head(have_cols, N_SAMPLES)
use_groups <- our_samples[match(sub("\\.AVG_Beta$", "", use_cols), sample_names)]

cat("crosscheck: comparing", length(use_cols), "sample(s):", paste(use_groups, collapse = ", "), "\n")

## ---- Read only NAME + the chosen AVG_Beta columns (never the whole file) ----
dt <- fread(profile_path, select = c("NAME", use_cols), sep = "\t", header = TRUE)

## ---- Pick a probe subset present in BOTH our betas and the profile ----------
common_probes <- intersect(rownames(betas), dt$NAME)
if (!length(common_probes))
    stop("crosscheck_genomestudio_betas: zero probe overlap between our betas' rownames and ",
         "Methylation_Profile.txt's NAME column — check ARRAY_VERSION/canonicalization.")
set.seed(1)
probe_subset <- sample(common_probes, min(N_PROBES, length(common_probes)))

dt_sub <- dt[dt$NAME %in% probe_subset, ]
dt_sub <- dt_sub[!duplicated(dt_sub$NAME), ]  # NAME can repeat if a v2 replicate wasn't yet collapsed upstream
rownames(dt_sub) <- dt_sub$NAME

## ---- Compare, per sample ----------------------------------------------------
results <- lapply(seq_along(use_groups), function(i) {
    ours  <- betas[dt_sub$NAME, use_groups[i]]
    theirs <- as.numeric(dt_sub[[use_cols[i]]])
    ok <- is.finite(ours) & is.finite(theirs)
    data.frame(Sample_Group = use_groups[i], n_probes = sum(ok),
               pearson_r = suppressWarnings(cor(ours[ok], theirs[ok])),
               mean_abs_diff = mean(abs(ours[ok] - theirs[ok])))
})
out <- do.call(rbind, results)
print(out)

low_r <- out$pearson_r < 0.9
if (any(low_r, na.rm = TRUE))
    cat("\nWARNING: pearson_r < 0.9 for:", paste(out$Sample_Group[which(low_r)], collapse = ", "),
        "— GenomeStudio's beta formula differs from minfi/noob/dasen by design, so some deviation",
        "is expected, but a very low r may indicate a sample/probe identity mismatch worth checking.\n")
cat("\ncrosscheck_genomestudio_betas: done -", nrow(out), "sample(s),", length(probe_subset), "probe(s) each\n")
