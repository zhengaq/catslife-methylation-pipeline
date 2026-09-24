### stage6/pca_sex_batch.R: structure check. Re-runs the stage-2 PCA on the most-variable
### CpGs twice, with and without the sex chromosomes, and colors PC1/PC2 by sex and by
### Sample_Plate. Shows how much of the leading variation is sex, and the batch structure
### underneath it, and lists samples whose methylation sex disagrees with their admin sex.
### Needs the stage-1 dasen betas, the array annotation package and minfi.
### Outputs: REPORT_DIR/{PCA_sex_batch.pdf, pca_sex_batch_summary.csv, pca_sex_mismatch.csv}

source("config.R")
suppressMessages(library(minfi))

if (!file.exists(F_DASENB))
    stop("stage6/pca_sex_batch.R: missing ", F_DASENB, "; run stage 1 first.")

## inputs ----
dv <- load_one(F_DASENB)
M  <- canonicalize_v2_probe_ids(dv$M)          # bare cg ids, matches the annotation
rm(dv); gc()

targets <- load_targets()[, c("Sample_Group", "DNA_Source", "Sample_Plate")]
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)

sqc <- read.csv(file.path(REPORT_DIR, "sex_qc.csv"), stringsAsFactors = FALSE)
sex <- setNames(sqc$Sex_admin, sqc$Sample)     # Sample == Sentrix id == beta colname
sid <- setNames(sqc$Subject_ID, sqc$Sample)

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
sx    <- factor(sex[samp], levels = c("F", "M"))
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

## Sex mismatches: samples on the other sex's side of the sex component, i.e. a mislabeled
## (swapped) sample or a sex-chromosome anomaly. The sex component is the PC (of the first 8,
## sex chromosomes included) with the highest R^2 on admin sex; each sample is assigned to the
## nearer of the female and male means. Skipped when no PC tracks sex (R^2 < 0.5).
r2_all <- out$r2_sex[out$pca == "all_cpg_incl_sexchr"]
k      <- if (any(is.finite(r2_all))) which.max(r2_all) else NA_integer_
mism_f <- file.path(REPORT_DIR, "pca_sex_mismatch.csv")
mism   <- data.frame(Sample = character(), Subject_ID = character(), Sex_admin = character(),
                     Sex_pca = character(), PC = integer(), score = numeric())
if (!is.na(k) && r2_all[k] >= 0.5) {
    pc   <- pca_all$x[, k]
    mu   <- tapply(pc, sx, mean)
    call <- ifelse(abs(pc - mu[["F"]]) < abs(pc - mu[["M"]]), "F", "M")
    bad  <- !is.na(sx) & call != as.character(sx)
    mism <- data.frame(Sample = samp[bad], Subject_ID = unname(sid[samp[bad]]),
                       Sex_admin = as.character(sx[bad]), Sex_pca = call[bad], PC = k,
                       score = round(pc[bad], 3))
    cat(sprintf("sex mismatch: PC%d (R^2 sex %.2f); %d sample(s) whose PCA sex differs from admin sex\n",
                k, r2_all[k], nrow(mism)))
} else cat("sex mismatch: skipped, no leading PC tracks admin sex (max R^2 =",
           if (is.na(k)) "NA" else round(r2_all[k], 2), ")\n")
write.csv(mism, mism_f, row.names = FALSE)

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
dev.off()
cat("stage6/pca_sex_batch: wrote PCA_sex_batch.pdf, pca_sex_batch_summary.csv, pca_sex_mismatch.csv to",
    REPORT_DIR, "\n")
