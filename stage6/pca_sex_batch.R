### stage6/pca_sex_batch.R: structure check. Re-runs the stage-2 PCA on the most-variable
### CpGs twice, with and without the sex chromosomes, and colors PC1/PC2 by sex and by
### Sample_Plate. Shows how much of the leading variation is sex, and the batch structure
### underneath it. Calls each sample's sex from its position on the sex principal component
### without using labels, and lays the calls out by person so that a person whose arrays
### disagree in sex, or move atypically on the sex component, stands out.
### Needs the stage-1 dasen betas, sex_qc.csv from build_phenotype_file.R, the array
### annotation package and minfi.
### Outputs: REPORT_DIR/{PCA_sex_batch.pdf, pca_sex_batch_summary.csv, pca_sex_scores.csv,
###          pca_sex_mismatch.csv, pca_sex_by_person.csv}
### Stops (after writing every output) if a cross-sex relabel in SAMPLE_SWAPS_FILE is not
### detected under its sheet label, or is still a mismatch under its corrected label.

source("config.R")
suppressMessages(library(minfi))

if (!file.exists(F_DASENB))
    stop("stage6/pca_sex_batch.R: missing ", F_DASENB, "; run stage 1 first.")

## helpers ----
## Two-component 1-D Gaussian mixture fitted by EM, without labels. Returns the component
## means, sds and weights (component 1 = lower mean) and each value's posterior probability
## of belonging to component 2. Posteriors are computed in log space, so values far from both
## components still get a finite, well-defined probability.
fit_two_gaussians <- function(x, max_iter = 1000L, tol = 1e-10) {
    ok <- is.finite(x)
    if (sum(ok) < 4) stop("fit_two_gaussians: fewer than 4 finite values")
    xs <- x[ok]
    km <- stats::kmeans(xs, centers = range(xs))
    mu <- as.numeric(km$centers); g <- km$cluster
    floor_sd <- 1e-6 * diff(range(xs))
    sd <- pmax(c(stats::sd(xs[g == 1]), stats::sd(xs[g == 2])), floor_sd, na.rm = TRUE)
    w  <- as.numeric(table(factor(g, levels = 1:2))) / length(xs)
    post2 <- function(v, mu, sd, w) {
        l1 <- log(w[1]) + stats::dnorm(v, mu[1], sd[1], log = TRUE)
        l2 <- log(w[2]) + stats::dnorm(v, mu[2], sd[2], log = TRUE)
        1 / (1 + exp(l1 - l2))
    }
    ll_old <- -Inf
    for (it in seq_len(max_iter)) {
        r2 <- post2(xs, mu, sd, w); r1 <- 1 - r2
        w  <- c(mean(r1), mean(r2))
        mu <- c(sum(r1 * xs) / sum(r1), sum(r2 * xs) / sum(r2))
        sd <- pmax(sqrt(c(sum(r1 * (xs - mu[1])^2) / sum(r1), sum(r2 * (xs - mu[2])^2) / sum(r2))), floor_sd)
        ll <- sum(log(w[1] * stats::dnorm(xs, mu[1], sd[1]) + w[2] * stats::dnorm(xs, mu[2], sd[2])))
        if (is.finite(ll) && abs(ll - ll_old) < tol * (1 + abs(ll))) break
        ll_old <- ll
    }
    if (mu[1] > mu[2]) { mu <- rev(mu); sd <- rev(sd); w <- rev(w) }
    p <- rep(NA_real_, length(x)); p[ok] <- post2(xs, mu, sd, w)
    list(mean = mu, sd = sd, weight = w, p_upper = p, iterations = it)
}

## Robust z of each value against a reference group (median / MAD). The MAD is floored so a
## near-constant group does not turn ordinary values into extreme z.
robust_z <- function(x, ref) {
    ref <- ref[is.finite(ref)]
    s <- max(stats::mad(ref), 1e-6 * max(1, abs(stats::median(ref))))
    (x - stats::median(ref)) / s
}

## Classify each SAMPLE_SWAPS_FILE action by what a sex check can see: a relabel whose sheet
## and corrected persons differ in admin sex is testable by sex, a same-sex relabel is not.
relabel_type <- function(identity_action, sex_sheet, sex_admin) {
    ifelse(!identity_action %in% "relabel",
           ifelse(identity_action %in% "flag", "flag only", NA_character_),
    ifelse(is.na(sex_sheet) | is.na(sex_admin), "relabel, sex unknown",
    ifelse(sex_sheet != sex_admin, "cross-sex relabel", "same-sex relabel")))
}

## inputs ----
dv <- load_one(F_DASENB)
M  <- canonicalize_v2_probe_ids(dv$M)          # bare cg ids, matches the annotation
rm(dv); gc()

targets <- load_targets()[, c("Sample_Group", "DNA_Source", "Sample_Plate")]
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)

## sex_qc.csv (Sample == Sentrix id == beta colname) with the sheet-label provenance columns.
sqc_f <- file.path(REPORT_DIR, "sex_qc.csv")
if (!file.exists(sqc_f))
    stop("stage6/pca_sex_batch.R: missing ", sqc_f, "; run scripts/build/build_phenotype_file.R first.")
sqc  <- read.csv(sqc_f, stringsAsFactors = FALSE, colClasses = "character")
need <- c("Sample", "Subject_ID", "Subject_ID_sheet", "identity_action", "IndividualID",
          "Sex_admin", "Sex_admin_sheet")
if (length(miss <- setdiff(need, names(sqc))))
    stop("stage6/pca_sex_batch.R: ", sqc_f, " lacks column(s) ", paste(miss, collapse = ", "),
         "; rebuild it with scripts/build/build_phenotype_file.R.")
sqc[sqc == ""] <- NA_character_

## sex-chromosome probe set (from the array annotation) ----
anno_pkg <- if (ARRAY_VERSION == "v2") {
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
} else {
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
}
if (!requireNamespace(anno_pkg, quietly = TRUE))
    stop("stage6/pca_sex_batch.R: annotation package not installed: ", anno_pkg)
## Attach the annotation package with library(): minfi::getAnnotation() -> updateObject()
## looks it up as "package:<anno_pkg>" on the search list.
suppressMessages(library(anno_pkg, character.only = TRUE))
anno      <- minfi::getAnnotation(get(anno_pkg))
anno_bare <- sub("_[A-Za-z0-9]+$", "", rownames(anno))
chr_of    <- setNames(as.character(anno$chr), anno_bare)
probe_chr <- chr_of[rownames(M)]
is_sexchr <- probe_chr %in% c("chrX", "chrY")
cat(sprintf("probes: %d total | chrX/chrY %d | unmapped-to-anno %d\n",
            nrow(M), sum(is_sexchr, na.rm = TRUE), sum(is.na(probe_chr))))

## PCA on the top-var CpGs of a probe set (mirrors stage 2) ----
run_pca <- function(Msub) {
    vary <- matrixStats::rowVars(Msub)
    keep <- rownames(Msub)[order(vary, decreasing = TRUE)[seq_len(min(PCA_NCPG, nrow(Msub)))]]
    Xs   <- scale(t(Msub[keep, ]))
    Xs   <- Xs[, colSums(!is.finite(Xs)) == 0, drop = FALSE]
    prcomp(Xs, scale. = FALSE, rank. = 10)
}
## Restrict to blood/saliva samples, as stage 2 does (drops Cell_Line controls).
keep_s  <- targets$Sample_Group[targets$DNA_Source %in% c("Buffy_Coat", "PBMC", "Saliva")]
M       <- M[, colnames(M) %in% keep_s, drop = FALSE]
pca_all <- run_pca(M)                              # baseline: includes sex chromosomes
pca_aut <- run_pca(M[!(is_sexchr %in% TRUE), ])    # autosomes only

## per-sample sex + plate aligned to PCA row order ----
samp  <- rownames(pca_all$x)
qi    <- match(samp, sqc$Sample)
sx    <- factor(sqc$Sex_admin[qi], levels = c("F", "M"))
plate <- factor(targets$Sample_Plate[match(samp, targets$Sample_Group)])

## Variance explained by a grouping factor g on a PC vector (one-way ANOVA R^2).
r2 <- function(pc, g) { ok <- !is.na(g); if (length(unique(g[ok])) < 2) return(NA_real_)
    summary(lm(pc[ok] ~ g[ok]))$r.squared }
summ <- function(p, tag) data.frame(
    pca = tag, PC = 1:min(8, ncol(p$x)),
    var_explained = round(summary(p)$importance[2, 1:min(8, ncol(p$x))], 3),
    r2_sex   = round(sapply(1:min(8, ncol(p$x)), function(k) r2(p$x[, k], sx)),    3),
    r2_plate = round(sapply(1:min(8, ncol(p$x)), function(k) r2(p$x[, k], plate)), 3))
out <- rbind(summ(pca_all, "all_cpg_incl_sexchr"), summ(pca_aut, "autosomes_only"))
write.csv(out, file.path(REPORT_DIR, "pca_sex_batch_summary.csv"), row.names = FALSE)
cat("PC1 baseline: var =", out$var_explained[1], "| r2 sex =", out$r2_sex[1], "| r2 plate =", out$r2_plate[1], "\n")
cat("PC1 autosomes-only: var =", out$var_explained[out$pca == "autosomes_only"][1],
    "| r2 sex =", out$r2_sex[out$pca == "autosomes_only"][1], "\n")

## Sex call on the sex component ----
## The sex component is the PC (of the first 8, sex chromosomes included) with the highest R^2
## on admin sex. The call itself uses no labels: a two-component Gaussian mixture on that PC
## gives each sample a posterior probability of belonging to the male cluster (the component
## holding most admin-male samples), so the call can be scored against any label set. A sample
## is ambiguous when that probability is between SEX_AMBIG_P and 1 - SEX_AMBIG_P, and atypical
## when it lies more than SEX_ATYPICAL_Z robust SDs from the centre of its called cluster; tight
## clusters give near-certain posteriors even to samples far from both, so the atypical rule is
## what finds partial or mixed signals.
r2_all <- out$r2_sex[out$pca == "all_cpg_incl_sexchr"]
k <- if (any(is.finite(r2_all))) which.max(r2_all) else NA_integer_
if (is.na(k) || r2_all[k] < 0.5)
    stop("stage6/pca_sex_batch.R: no leading PC tracks admin sex (max R^2 = ",
         if (is.na(k)) "NA" else round(r2_all[k], 2), "); the PCA sex check cannot run. ",
         "Check that sex_qc.csv's Sample ids match the beta matrix columns.")
pc  <- pca_all$x[, k]
fit <- fit_two_gaussians(pc)
upper_is_m <- mean(pc[sx %in% "M"]) > mean(pc[sx %in% "F"])
p_male  <- if (upper_is_m) fit$p_upper else 1 - fit$p_upper
sex_pca <- ifelse(p_male >= 0.5, "M", "F")
ctr_f   <- median(pc[sex_pca == "F"]); ctr_m <- median(pc[sex_pca == "M"])
z_call  <- ifelse(sex_pca == "M", robust_z(pc, pc[sex_pca == "M"]), robust_z(pc, pc[sex_pca == "F"]))
ambig   <- p_male > SEX_AMBIG_P & p_male < 1 - SEX_AMBIG_P
atyp    <- abs(z_call) > SEX_ATYPICAL_Z
cat(sprintf("sex component: PC%d (R^2 sex %.2f) | cluster centres F %.1f, M %.1f | %d ambiguous, %d atypical\n",
            k, r2_all[k], ctr_f, ctr_m, sum(ambig), sum(atyp)))

## Per-sample scores and call ----
sexpc_cols <- setNames(as.data.frame(round(pca_all$x, 4)), paste0("all_PC", seq_len(ncol(pca_all$x))))
autpc_cols <- setNames(as.data.frame(round(pca_aut$x[samp, ], 4)), paste0("aut_PC", seq_len(ncol(pca_aut$x))))
scores <- data.frame(
    Sample = samp, Sample_Plate = as.character(plate),
    Subject_ID = sqc$Subject_ID[qi], Subject_ID_sheet = sqc$Subject_ID_sheet[qi],
    identity_action = sqc$identity_action[qi], IndividualID = sqc$IndividualID[qi],
    Sex_admin = sqc$Sex_admin[qi], Sex_admin_sheet = sqc$Sex_admin_sheet[qi],
    sex_PC = k, score = round(pc, 3),
    position = round((pc - ctr_f) / (ctr_m - ctr_f), 3),    # 0 = female centre, 1 = male centre
    P_male = signif(p_male, 4), Sex_pca = sex_pca, robust_z = round(z_call, 2),
    ambiguous = ambig, atypical = atyp, stringsAsFactors = FALSE)
scores$relabel_type <- relabel_type(scores$identity_action, scores$Sex_admin_sheet, scores$Sex_admin)
scores <- cbind(scores, sexpc_cols, autpc_cols)
rownames(scores) <- NULL
write.csv(scores, file.path(REPORT_DIR, "pca_sex_scores.csv"), row.names = FALSE)

## The call scored against the sheet labels (pre) and the corrected labels (post).
mm_pre  <- !is.na(scores$Sex_admin_sheet) & scores$Sex_pca != scores$Sex_admin_sheet
mm_post <- !is.na(scores$Sex_admin) & scores$Sex_pca != scores$Sex_admin
cat("sheet labels: ", sum(mm_pre), " mismatch(es) | corrected labels: ", sum(mm_post),
    " mismatch(es)\n", sep = "")

## Post-correction mismatches (the labels the pipeline uses) ----
mism_cols <- c("Sample", "Subject_ID", "Sex_admin", "Sex_pca", "sex_PC", "score", "position",
               "P_male", "robust_z", "ambiguous", "atypical")
mism <- scores[mm_post, mism_cols]
write.csv(mism, file.path(REPORT_DIR, "pca_sex_mismatch.csv"), row.names = FALSE)
cat(sprintf("sex mismatch: %d sample(s) whose PCA sex differs from admin sex\n", nrow(mism)))

## By person ----
## One row per array, grouped by person (the Subject_ID's random_id) and ordered by wave, with
## duplicate aliquots after the wave's main array. Listed: every person with two or more arrays,
## and any single-array person whose array mismatches its admin sex, is atypical, or is in
## SAMPLE_SWAPS_FILE. Per-person columns: the called sex of each array in that order, whether the
## calls differ, and the largest change in sex-PC position between two of the person's arrays.
## That change is scored against the within-person differences of every multi-array person (its
## pairwise position differences, centred on their median, scaled by their MAD), so a systematic
## shift between waves does not count as atypical; a person is flagged when the score exceeds
## SEX_ATYPICAL_Z.
bp <- scores[!is.na(scores$Subject_ID), ]
bp$mismatch <- (mm_post)[!is.na(scores$Subject_ID)]
bp$Person   <- subject_base_id(bp$Subject_ID)
bp$Wave     <- subject_wave(bp$Subject_ID)
bp <- bp[order(bp$Person, bp$Wave, grepl("D$", bp$Subject_ID), bp$Sample), ]
grp  <- split(seq_len(nrow(bp)), bp$Person)
pair_d <- lapply(grp, function(i) if (length(i) < 2) numeric() else {
    pr <- utils::combn(i, 2); bp$position[pr[2, ]] - bp$position[pr[1, ]] })
d_all <- unlist(pair_d, use.names = FALSE)
if (length(d_all) < 3)
    cat("by person: fewer than 3 within-person position differences; position change not scored\n")
person <- data.frame(
    Person   = as.integer(names(grp)),
    n_arrays = lengths(grp),
    sex_calls = vapply(grp, function(i) paste(bp$Sex_pca[i], collapse = "/"), character(1)),
    sex_differs = vapply(grp, function(i) if (length(i) < 2) NA else length(unique(bp$Sex_pca[i])) > 1, logical(1)),
    position_change = vapply(grp, function(i) if (length(i) < 2) NA_real_ else round(diff(range(bp$position[i])), 3),
                             numeric(1)),
    position_change_z = vapply(pair_d, function(d) if (!length(d) || length(d_all) < 3) NA_real_
                               else round(max(abs(robust_z(d, d_all))), 2), numeric(1)),
    stringsAsFactors = FALSE)
person$position_change_atypical <- person$position_change_z > SEX_ATYPICAL_Z
bp <- merge(bp, person, by = "Person", sort = FALSE)
bp <- bp[order(bp$Person, bp$Wave, grepl("D$", bp$Subject_ID), bp$Sample), ]
single_flag <- tapply(bp$mismatch | bp$atypical | !is.na(bp$identity_action), bp$Person, any)
keep_p <- person$Person[person$n_arrays >= 2 | single_flag[as.character(person$Person)]]
bp <- bp[bp$Person %in% keep_p, c(
    "Person", "IndividualID", "n_arrays", "Wave", "Subject_ID", "Sample", "Sample_Plate",
    "Subject_ID_sheet", "identity_action", "relabel_type", "Sex_admin_sheet", "Sex_admin", "Sex_pca",
    "score", "position", "P_male", "robust_z", "atypical",
    "sex_calls", "sex_differs", "position_change", "position_change_z", "position_change_atypical")]
write.csv(bp, file.path(REPORT_DIR, "pca_sex_by_person.csv"), row.names = FALSE)
pp_kept <- person[person$Person %in% keep_p, ]
cat(sprintf("by person: %d person(s) listed (%d with 2+ arrays); %d whose calls differ between arrays, %d with an atypical position change\n",
            nrow(pp_kept), sum(pp_kept$n_arrays >= 2), sum(pp_kept$sex_differs %in% TRUE),
            sum(pp_kept$position_change_atypical %in% TRUE)))

## plots ----
scatter <- function(p, g, ttl, leg) {
    plot(p$x[, 1], p$x[, 2], col = as.integer(g), pch = 19, cex = 0.6,
         xlab = "PC 1", ylab = "PC 2", main = ttl)
    if (!is.null(leg)) legend("topleft", legend = levels(g), col = seq_along(levels(g)),
                              pch = 19, bty = "n", cex = 0.8, title = leg)
}
pdf(file.path(REPORT_DIR, "PCA_sex_batch.pdf"), height = 9, width = 13)
par(mfrow = c(2, 3))
plot(summary(pca_all)$importance[2, 1:10], type = "b", ylim = c(0, 1),
     ylab = "Prop. variance", xlab = "PC", main = "Scree: all CpGs (incl. sex chr)")
scatter(pca_all, sx,    "All CpGs, colored by sex",   "Sex")
scatter(pca_all, plate, "All CpGs, colored by plate", NULL)
plot(summary(pca_aut)$importance[2, 1:10], type = "b", ylim = c(0, 1),
     ylab = "Prop. variance", xlab = "PC", main = "Scree: autosomes only")
scatter(pca_aut, sx,    "Autosomes, colored by sex",   "Sex")
scatter(pca_aut, plate, "Autosomes, colored by plate", NULL)

## Page 2: the sex component, and each flagged or relabeled array against both label sets.
par(mfrow = c(1, 2), mar = c(5, 9, 4, 2), oma = c(2.5, 0, 0, 0))
hist(pc, breaks = 120, col = "grey80", border = NA, xlab = paste0("PC", k, " score"),
     main = paste0("Sex component (PC", k, ")"))
abline(v = c(ctr_f, ctr_m), lty = 2)
rug(pc[atyp | ambig], col = "red", lwd = 2)
legend("top", c("cluster centres", "ambiguous / atypical"), lty = c(2, 1), col = c("black", "red"),
       bty = "n", cex = 0.8)
show <- scores[mm_pre | mm_post | scores$relabel_type %in% "cross-sex relabel" |
               scores$atypical | scores$ambiguous, ]
if (nrow(show)) {
    show <- show[order(show$score), ]
    y <- seq_len(nrow(show))
    sex_col <- c(F = "firebrick", M = "steelblue")
    plot(show$score, y, xlim = range(pc), yaxt = "n", pch = NA, ylab = "",
         xlab = paste0("PC", k, " score"), main = "Flagged and relabeled arrays")
    axis(2, at = y, las = 1, cex.axis = 0.7,
         labels = ifelse(show$Subject_ID_sheet %in% show$Subject_ID & show$Subject_ID_sheet == show$Subject_ID,
                         show$Subject_ID, paste0(show$Subject_ID_sheet, " -> ", show$Subject_ID)))
    abline(v = c(ctr_f, ctr_m), lty = 2, col = "grey60")
    points(show$score, y, pch = 1,  cex = 1.8, col = sex_col[show$Sex_admin_sheet])
    points(show$score, y, pch = 19, cex = 0.9, col = sex_col[show$Sex_admin])
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    legend("bottom", c("admin sex of sheet label (open ring)", "admin sex of corrected label (filled dot)",
                       "female", "male"), pch = c(1, 19, 15, 15),
           col = c("black", "black", sex_col), bty = "n", horiz = TRUE, cex = 0.85, text.width = NA)
} else plot.new()
dev.off()
cat("stage6/pca_sex_batch: wrote PCA_sex_batch.pdf, pca_sex_batch_summary.csv, pca_sex_scores.csv,",
    "pca_sex_mismatch.csv, pca_sex_by_person.csv to", REPORT_DIR, "\n")

## Known-positive check ----
## Every cross-sex relabel is a sample whose sheet label names a person of the other sex, so the
## call must contradict the sheet label and agree with the corrected one. A relabel that is not
## detected, or any relabel that leaves a mismatch the sheet label did not have, means the check,
## the beta-to-label alignment, or the correction itself is wrong.
is_xs <- scores$relabel_type %in% "cross-sex relabel"
ok_xs <- mm_pre & !mm_post
fail  <- (is_xs & !ok_xs) | (!is.na(scores$identity_action) & !mm_pre & mm_post)
cat(sprintf("known positives: %d cross-sex relabel(s), %d detected under the sheet label and resolved by the relabel\n",
            sum(is_xs), sum(is_xs & ok_xs)))
if (any(fail)) {
    f <- scores[fail, ]
    stop("stage6/pca_sex_batch.R: the PCA sex call disagrees with SAMPLE_SWAPS_FILE for ",
         paste0(f$Sample, " (sheet ", f$Subject_ID_sheet, " [", f$Sex_admin_sheet, "] -> ", f$Subject_ID,
                " [", f$Sex_admin, "], called ", f$Sex_pca, ")", collapse = ", "),
         ". A cross-sex relabel should mismatch its sheet label's admin sex and match its corrected ",
         "label's; see pca_sex_by_person.csv.")
}
