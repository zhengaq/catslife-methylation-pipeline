### stage6/pca_sex_batch.R — Stage 6 structure check. Re-runs the stage-2 PCA two ways
### (all most-variable CpGs incl. sex chromosomes == the stage-2 baseline; then autosomes
### only) and re-colors PC1/PC2 by SEX and by Sample_Plate instead of tissue. Tests whether
### the stage-2 PC1 (~80% of variance) is sex, and surfaces the technical/batch structure
### that sex masks. Needs the stage-1 dasen betas + the array annotation + minfi.
### Outputs: REPORT_DIR/PCA_sex_batch.pdf + REPORT_DIR/pca_sex_batch_summary.csv

source("config.R")
suppressMessages(library(minfi))

if (!file.exists(F_DASENB))
    stop("stage6/pca_sex_batch.R: missing ", F_DASENB, " — run stage 1 first.")

## ---- inputs -------------------------------------------------------------
dv <- load_one(F_DASENB)
M  <- canonicalize_v2_probe_ids(dv$M)          # bare cg ids, matches the annotation
rm(dv); gc()

targets <- load_targets()[, c("Sample_Group", "DNA_Source", "Sample_Plate")]
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)

sqc <- read.csv(file.path(REPORT_DIR, "sex_qc.csv"), stringsAsFactors = FALSE)
sex <- setNames(sqc$Sex_admin, sqc$Sample)     # Sample == Sentrix id == beta colname

## ---- sex-chromosome probe set (from the array annotation) ---------------
anno_pkg <- if (ARRAY_VERSION == "v2") {
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
} else {
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
}
if (!requireNamespace(anno_pkg, quietly = TRUE))
    stop("stage6/pca_sex_batch.R: annotation package not installed: ", anno_pkg)
anno      <- minfi::getAnnotation(getExportedValue(anno_pkg, anno_pkg))
anno_bare <- sub("_[A-Za-z0-9]+$", "", rownames(anno))
chr_of    <- setNames(as.character(anno$chr), anno_bare)
probe_chr <- chr_of[rownames(M)]
is_sexchr <- probe_chr %in% c("chrX", "chrY")
cat(sprintf("probes: %d total | chrX/chrY %d | unmapped-to-anno %d\n",
            nrow(M), sum(is_sexchr, na.rm = TRUE), sum(is.na(probe_chr))))

## ---- PCA on the top-var CpGs of a probe set (mirrors stage 2) -----------
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

## ---- per-sample sex + plate aligned to PCA row order --------------------
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

## ---- plots --------------------------------------------------------------
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
scatter(pca_all, sx,    "All CpGs — colored by SEX",   "Sex")
scatter(pca_all, plate, "All CpGs — colored by PLATE", NULL)
plot(summary(pca_aut)$importance[2, 1:10], type = "b", ylim = c(0, 1),
     ylab = "Prop. variance", xlab = "PC", main = "Scree: autosomes only")
scatter(pca_aut, sx,    "Autosomes — colored by SEX",   "Sex")
scatter(pca_aut, plate, "Autosomes — colored by PLATE", NULL)
dev.off()
cat("stage6/pca_sex_batch: wrote PCA_sex_batch.pdf + pca_sex_batch_summary.csv to", REPORT_DIR, "\n")
