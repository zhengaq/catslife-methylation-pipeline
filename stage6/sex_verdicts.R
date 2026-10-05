### stage6/sex_verdicts.R: one verdict per array of interest, combining the PCA sex call
### (pca_sex_scores.csv), the intensity sex call (sex_intensity.csv), SAMPLE_SWAPS_FILE, and
### the other arrays of the same person. Arrays of interest are every SAMPLE_SWAPS_FILE row
### plus every array any check flags: a PCA or intensity sex that differs from admin sex, an
### ambiguous or atypical PCA position, a non-XX/XY intensity pattern, or PCA and intensity
### disagreeing. Nothing is dropped; the verdicts are for review.
### Needs: pca_sex_batch.R and sex_intensity.R outputs, sex_qc.csv, SAMPLE_SWAPS_FILE.
### Outputs: REPORT_DIR/{sex_identity_verdicts.csv, sex_identity_memo.md}

source("config.R"); source("stage6/sex_helpers.R")

need_f <- file.path(REPORT_DIR, c("pca_sex_scores.csv", "sex_intensity.csv"))
if (any(!file.exists(need_f)))
    stop("stage6/sex_verdicts.R: missing ", paste(need_f[!file.exists(need_f)], collapse = ", "),
         "; run stage6/pca_sex_batch.R and stage6/sex_intensity.R first.")
sc  <- read.csv(need_f[1], colClasses = "character")
si  <- read.csv(need_f[2], colClasses = "character")
sqc <- read_sex_qc()
swaps <- read_sample_swaps(SAMPLE_SWAPS_FILE)

## one row per array known to any input ----
v <- data.frame(Sample = unique(c(sqc$Sample, sc$Sample, si$Sample)), stringsAsFactors = FALSE)
qi <- match(v$Sample, sqc$Sample); ci <- match(v$Sample, sc$Sample); ii <- match(v$Sample, si$Sample)
v$Subject_ID_sheet <- sqc$Subject_ID_sheet[qi]
v$Subject_ID       <- sqc$Subject_ID[qi]
v$IndividualID     <- sqc$IndividualID[qi]
v$identity_action  <- sqc$identity_action[qi]
v$C12_note         <- swaps$notes[match(sheet_label(v$Subject_ID_sheet), swaps$from)]
v$Sex_admin_sheet  <- sqc$Sex_admin_sheet[qi]
v$Sex_admin        <- sqc$Sex_admin[qi]
v$relabel_type     <- relabel_type(v$identity_action, v$Sex_admin_sheet, v$Sex_admin)
v$Sex_pca   <- sc$Sex_pca[ci]
v$pca_score <- as.numeric(sc$score[ci]); v$pca_position <- as.numeric(sc$position[ci])
v$P_male    <- as.numeric(sc$P_male[ci]); v$pca_robust_z <- as.numeric(sc$robust_z[ci])
v$pca_ambiguous <- sc$ambiguous[ci] %in% "TRUE"; v$pca_atypical <- sc$atypical[ci] %in% "TRUE"
v$Sex_int   <- si$Sex_int[ii]; v$int_pattern <- si$pattern[ii]
v$yMinusX   <- as.numeric(si$yMinusX[ii]); v$XCI_beta <- as.numeric(si$XCI_beta[ii])
v$Y_detected <- as.numeric(si$Y_detected[ii])
v$stage1_qc_failed <- si$stage1_qc_failed[ii] %in% "TRUE"

## The array's sex from the methods: the shared call when PCA and intensity agree, the one
## available when only one ran (an array dropped by stage-1 QC has no PCA call), else NA.
v$Sex_methyl <- ifelse(is.na(v$Sex_pca), v$Sex_int, ifelse(is.na(v$Sex_int) | v$Sex_int == v$Sex_pca, v$Sex_pca, NA))
v$methods_disagree <- !is.na(v$Sex_pca) & !is.na(v$Sex_int) & v$Sex_pca != v$Sex_int

## other arrays of the same person (post-correction labels) ----
other <- lapply(seq_len(nrow(v)), function(i) {
    if (is.na(v$IndividualID[i])) return(v[0, ])
    v[v$IndividualID %in% v$IndividualID[i] & v$Sample != v$Sample[i], ]
})
v$other_arrays <- vapply(other, function(o) if (!nrow(o)) NA_character_ else
    paste0(o$Sample, " (", o$Subject_ID, "): PCA ", o$Sex_pca, ", intensity ", o$Sex_int, collapse = "; "),
    character(1))
v$other_reads_admin <- vapply(seq_len(nrow(v)), function(i) {
    o <- other[[i]]; nrow(o) > 0 && any(o$Sex_methyl %in% v$Sex_admin[i])
}, logical(1))
v$other_all_opposite <- vapply(seq_len(nrow(v)), function(i) {
    o <- other[[i]]; nrow(o) > 0 && !is.na(v$Sex_admin[i]) && all(!is.na(o$Sex_methyl) & o$Sex_methyl != v$Sex_admin[i])
}, logical(1))

## arrays of interest ----
v$flagged <- with(v, (!is.na(Sex_admin) & !is.na(Sex_methyl) & Sex_methyl != Sex_admin) |
                     (!is.na(Sex_admin) & !is.na(Sex_pca) & Sex_pca != Sex_admin) |
                     (!is.na(Sex_admin) & !is.na(Sex_int) & Sex_int != Sex_admin) |
                     pca_ambiguous | pca_atypical | methods_disagree |
                     (!is.na(int_pattern) & !int_pattern %in% c("XX", "XY")))
v <- v[v$flagged | !is.na(v$identity_action), ]

## verdicts, in priority order ----
verdict <- function(r) {
    normal <- r$int_pattern %in% c("XX", "XY") || is.na(r$int_pattern)
    pca_ok <- !(r$pca_ambiguous || r$pca_atypical)
    if (r$relabel_type %in% "cross-sex relabel")
        return(if (!r$methods_disagree && r$Sex_methyl %in% r$Sex_admin && normal && pca_ok)
                   "C12 relabel confirmed by sex" else "C12 relabel not confirmed by sex: review")
    if (r$methods_disagree) return("unresolved: PCA and intensity sex disagree")
    if (r$int_pattern %in% c("XXY-like", "X0-like"))
        return(paste0("possible sex-chromosome anomaly (", r$int_pattern, ")"))
    if (r$int_pattern %in% "atypical" || !pca_ok)
        return(if (r$stage1_qc_failed) "atypical sex signal on an array that failed stage-1 QC"
               else "atypical sex signal: possible mixed DNA")
    if (is.na(r$Sex_admin)) return("no admin sex (array not in the phenotype file)")
    if (!is.na(r$Sex_methyl) && r$Sex_methyl != r$Sex_admin) {
        if (r$other_reads_admin)
            return(paste0("likely sample mix-up: another array of this person reads ", r$Sex_admin))
        if (r$other_all_opposite)
            return(paste0("possible admin sex error: every array of this person reads ", r$Sex_methyl))
        return("sex mismatch on this person's only array: sample mix-up or admin sex error")
    }
    if (r$relabel_type %in% "same-sex relabel (not testable by sex)") return("C12 relabel, same-sex pair: not testable by sex")
    if (r$identity_action %in% "flag") return("C12 flag: sex consistent")
    "consistent"
}
v$verdict <- vapply(seq_len(nrow(v)), function(i) verdict(v[i, ]), character(1))
v$in_C12  <- !is.na(v$identity_action)
v$evidence <- with(v, paste0(
    "PCA ", ifelse(is.na(Sex_pca), "n/a", paste0(Sex_pca, " (position ", sprintf("%.2f", pca_position),
        ", P_male ", ifelse(P_male < 0.001, "<0.001", ifelse(P_male > 0.999, ">0.999", sprintf("%.3f", P_male))), ", z ", sprintf("%.1f", pca_robust_z), ")")),
    "; intensity ", ifelse(is.na(Sex_int), "n/a", paste0(Sex_int, " (", int_pattern, ", yMed-xMed ",
        sprintf("%.2f", yMinusX), ", XCI beta ", sprintf("%.2f", XCI_beta), ", Y detected ",
        sprintf("%.0f%%", 100 * Y_detected), ")"))))
v <- v[order(!v$in_C12, v$verdict, v$Subject_ID_sheet), c(
    "Sample", "Subject_ID_sheet", "Subject_ID", "IndividualID", "in_C12", "identity_action", "relabel_type",
    "C12_note", "Sex_admin_sheet", "Sex_admin", "Sex_pca", "Sex_int", "int_pattern", "verdict", "evidence",
    "other_arrays", "pca_score", "pca_position", "P_male", "pca_robust_z", "pca_ambiguous", "pca_atypical",
    "yMinusX", "XCI_beta", "Y_detected", "stage1_qc_failed")]
write.csv(v, file.path(REPORT_DIR, "sex_identity_verdicts.csv"), row.names = FALSE)
cat("sex_verdicts:", nrow(v), "array(s) of interest\n"); print(table(v$verdict))

## memo for the lab ----
md_table <- function(d) {
    if (!nrow(d)) return("None.\n")
    esc <- function(x) gsub("\\|", "/", ifelse(is.na(x), "", as.character(x)))
    paste0("| ", paste(names(d), collapse = " | "), " |\n|", strrep("---|", ncol(d)), "\n",
           paste0("| ", apply(d, 1, function(r) paste(esc(r), collapse = " | ")), " |", collapse = "\n"), "\n")
}
pp   <- read.csv(file.path(REPORT_DIR, "pca_sex_prepost.csv"), colClasses = "character")
conc <- read.csv(file.path(REPORT_DIR, "sex_intensity_concordance.csv"), colClasses = "character")
agr  <- conc[conc$first == "agree", ]
xs   <- pp[pp$relabel_type %in% "cross-sex relabel", ]
new  <- v[!v$in_C12, c("Sample", "Subject_ID", "Sex_admin", "Sex_pca", "Sex_int", "verdict", "evidence", "other_arrays")]
c12  <- v[v$in_C12, c("Sample", "Subject_ID_sheet", "Subject_ID", "C12_note", "Sex_admin_sheet", "Sex_admin",
                      "Sex_pca", "Sex_int", "verdict")]
memo <- c(
    "# Methylation sex check: sample identity review", "",
    paste0("Generated ", format(Sys.Date()), " from the stage-6 sex checks. Two independent methods call each ",
           "array's sex: its position on the sex principal component of the methylation data, and the ",
           "chrX/chrY probe intensities. Arrays are listed when either method disagrees with the admin sex, ",
           "gives an atypical signal, or when the array is in the sample-swap file. Nothing has been removed."), "",
    "## Summary", "",
    paste0("- Arrays with a PCA sex call: ", nrow(sc), "; with an intensity sex call: ", nrow(si), "."),
    paste0("- PCA and intensity calls agree on ", agr$n[1], " of ", agr$total[1], " arrays."),
    paste0("- Cross-sex relabels in the swap file: ", nrow(xs), "; detected under the original label and ",
           "consistent after the relabel: ", sum(xs$prepost == "resolved by relabel"), "."),
    paste0("- Same-sex relabels cannot be checked by sex: ", sum(pp$relabel_type %in% "same-sex relabel (not testable by sex)"), "."),
    paste0("- Arrays not in the swap file with a sex finding: ", nrow(new), "."), "",
    "## New findings (not in the swap file)", "", md_table(new),
    "## Swap-file rows", "", md_table(c12),
    "## Limits", "",
    "- A swap between two people of the same sex does not change either array's sex and is invisible here; genotype matching is needed for those.",
    "- A sample whose array reads the other sex can be a mix-up or an admin sex error. Another array of the same person (other wave or duplicate) separates the two; without one the finding is unresolved.",
    "- An atypical signal (between the sexes on the PCA, or an unusual X/Y intensity pattern) fits mixed DNA from two people; it can also come from a sex-chromosome anomaly or a poor-quality array.")
writeLines(memo, file.path(REPORT_DIR, "sex_identity_memo.md"))
cat("stage6/sex_verdicts: wrote sex_identity_verdicts.csv and sex_identity_memo.md to", REPORT_DIR, "\n")
