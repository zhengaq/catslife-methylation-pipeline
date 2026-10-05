### stage6/pca_sex_batch.R: structure check. Re-runs the stage-2 PCA on the most-variable
### CpGs twice, with and without the sex chromosomes, and colors PC1/PC2 by sex and by
### Sample_Plate. Shows how much of the leading variation is sex, and the batch structure
### underneath it. Calls each sample's sex from its position on the sex principal component
### without using labels, then scores that call against the sample sheet's labels before
### (pre) and after (post) the SAMPLE_SWAPS_FILE corrections.
### Needs the stage-1 dasen betas, sex_qc.csv from build_phenotype_file.R, the array
### annotation package and minfi.
### Outputs: REPORT_DIR/{PCA_sex_batch.pdf, pca_sex_batch_summary.csv, pca_sex_scores.csv,
###          pca_sex_mismatch.csv, pca_sex_prepost.csv, pca_sex_confusion.csv}
### Stops (after writing every output) if a cross-sex relabel in SAMPLE_SWAPS_FILE is not
### detected under its sheet label, or is still a mismatch under its corrected label.

source("config.R"); source("stage6/sex_helpers.R")
suppressMessages(library(minfi))

if (!file.exists(F_DASENB))
    stop("stage6/pca_sex_batch.R: missing ", F_DASENB, "; run stage 1 first.")

## inputs ----
dv <- load_one(F_DASENB)
M  <- canonicalize_v2_probe_ids(dv$M)          # bare cg ids, matches the annotation
rm(dv); gc()

targets <- load_targets()[, c("Sample_Group", "DNA_Source", "Sample_Plate")]
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)

sqc <- read_sex_qc()                           # Sample == Sentrix id == beta colname

## sex-chromosome probe set (from the array annotation) ----
anno      <- array_annotation()
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
## holding most admin-male samples). The same call is then scored against two label sets, so
## pre and post differ only in the labels. A sample is ambiguous when that probability is
## between SEX_AMBIG_P and 1 - SEX_AMBIG_P, and atypical when it lies more than SEX_ATYPICAL_Z
## robust SDs from the centre of its called cluster; tight clusters give near-certain posteriors
## even to samples far from both, so the atypical rule is what finds partial or mixed signals.
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
    ambiguous = ambig, atypical = atyp,
    mismatch_pre  = !is.na(sqc$Sex_admin_sheet[qi]) & sex_pca != sqc$Sex_admin_sheet[qi],
    mismatch_post = !is.na(sqc$Sex_admin[qi]) & sex_pca != sqc$Sex_admin[qi],
    stringsAsFactors = FALSE)
scores$relabel_type <- relabel_type(scores$identity_action, scores$Sex_admin_sheet, scores$Sex_admin)
scores$prepost <- with(scores, ifelse(
    mismatch_pre & !mismatch_post, "resolved by relabel",
    ifelse(!mismatch_pre & mismatch_post, "introduced by relabel",
    ifelse(mismatch_pre & mismatch_post, "mismatch under both labels",
    ifelse(is.na(Sex_admin) & is.na(Sex_admin_sheet), "no admin sex", "consistent under both labels")))))
scores <- cbind(scores, sexpc_cols, autpc_cols)
rownames(scores) <- NULL
write.csv(scores, file.path(REPORT_DIR, "pca_sex_scores.csv"), row.names = FALSE)

## Post-correction mismatches (the labels the pipeline uses) ----
mism_cols <- c("Sample", "Subject_ID", "Sex_admin", "Sex_pca", "sex_PC", "score", "position",
               "P_male", "robust_z", "ambiguous", "atypical")
mism <- scores[scores$mismatch_post, mism_cols]
write.csv(mism, file.path(REPORT_DIR, "pca_sex_mismatch.csv"), row.names = FALSE)
cat(sprintf("sex mismatch (post): %d sample(s) whose PCA sex differs from admin sex\n", nrow(mism)))

## Pre/post: the same call scored against the sheet labels (pre) and the corrected labels (post) ----
pp_rows <- with(scores, !is.na(identity_action) | mismatch_pre | mismatch_post | ambiguous | atypical)
pp <- scores[pp_rows, c("Sample", "Subject_ID_sheet", "Subject_ID", "identity_action", "relabel_type",
                        "Sex_admin_sheet", "Sex_admin", "Sex_pca", "score", "position", "P_male",
                        "robust_z", "ambiguous", "atypical", "mismatch_pre", "mismatch_post", "prepost")]
pp <- pp[order(is.na(pp$identity_action), pp$relabel_type, pp$Subject_ID_sheet), ]
write.csv(pp, file.path(REPORT_DIR, "pca_sex_prepost.csv"), row.names = FALSE)

conf <- do.call(rbind, lapply(c("pre", "post"), function(lab) {
    adm <- if (lab == "pre") scores$Sex_admin_sheet else scores$Sex_admin
    adm[is.na(adm)] <- "none"
    tb  <- as.data.frame(table(Sex_admin = adm, Sex_pca = scores$Sex_pca), stringsAsFactors = FALSE)
    tb$n_ambiguous_or_atypical <- mapply(function(a, s)
        sum(adm == a & scores$Sex_pca == s & (scores$ambiguous | scores$atypical)), tb$Sex_admin, tb$Sex_pca)
    cbind(labels = lab, tb)
}))
names(conf)[names(conf) == "Freq"] <- "n"
write.csv(conf, file.path(REPORT_DIR, "pca_sex_confusion.csv"), row.names = FALSE)
cat("pre/post: mismatches pre =", sum(scores$mismatch_pre), "| post =", sum(scores$mismatch_post),
    "| resolved by relabel =", sum(scores$prepost == "resolved by relabel"),
    "| introduced by relabel =", sum(scores$prepost == "introduced by relabel"), "\n")

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
show <- pp[pp$mismatch_pre | pp$mismatch_post | pp$relabel_type %in% "cross-sex relabel" |
           pp$atypical | pp$ambiguous, ]
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
    "pca_sex_mismatch.csv, pca_sex_prepost.csv, pca_sex_confusion.csv to", REPORT_DIR, "\n")

## Known-positive check ----
## Every cross-sex relabel is a sample whose sheet label names a person of the other sex, so the
## call must contradict the sheet label (pre) and agree with the corrected one (post). A relabel
## that is not detected, or that leaves a mismatch, means the check, the beta-to-label alignment,
## or the correction itself is wrong.
is_xs <- scores$relabel_type %in% "cross-sex relabel"
ok_xs <- scores$prepost %in% "resolved by relabel"
fail  <- scores[(is_xs & !ok_xs) | scores$prepost %in% "introduced by relabel", ]
cat(sprintf("known positives: %d cross-sex relabel(s), %d detected pre and resolved post\n",
            sum(is_xs), sum(is_xs & ok_xs)))
if (nrow(fail))
    stop("stage6/pca_sex_batch.R: the PCA sex call disagrees with SAMPLE_SWAPS_FILE for ",
         paste0(fail$Sample, " (sheet ", fail$Subject_ID_sheet, " -> ", fail$Subject_ID, ": ",
                fail$prepost, ")", collapse = ", "),
         ". Inspect pca_sex_prepost.csv: a cross-sex relabel should mismatch under its sheet label ",
         "and match under its corrected label.")
