### stage6/sex_intensity.R: sex call from the raw IDAT intensities, independent of the
### methylation PCA. Reads the IDATs in batches of SEX_INT_CHUNK samples and records, per
### sample:
###   xMed, yMed      median log2 total intensity (Meth + Unmeth) of the chrX / chrY probes;
###                   Sex_int applies minfi::getSex()'s rule (yMed - xMed < -2 => F), which is
###                   per-sample, so batching does not change it;
###   x_rel, y_rel    xMed / yMed minus the sample's median autosomal intensity, which removes
###                   array-wide brightness so X and Y copy number can be read separately;
###   Y_detected      fraction of chrY probes detected (detection p < DETP_THRESHOLD);
###   XCI_beta        mean raw beta of chrX CpG-island probes (about 0.5 with an inactive X,
###                   near 0 without one).
### Each of x_rel, y_rel and XCI_beta is classed F-typical, M-typical or atypical against the
### samples of each called sex (SEX_ATYPICAL_Z robust SDs); together they give a pattern:
### XX, XY, XXY-like (male Y, female X dosage and X inactivation), X0-like (no Y, male X dosage,
### no X inactivation), or atypical (anything else, e.g. mixed DNA or a failed array).
### Needs the IDATs (or the test profile's injected RGChannelSet), the array annotation, minfi;
### sex_qc.csv; and, for the PCA agreement, pca_sex_scores.csv from pca_sex_batch.R.
### Outputs: REPORT_DIR/{sex_intensity.csv, sex_intensity_concordance.csv, sex_intensity.pdf}

source("config.R"); source("stage6/sex_helpers.R")
suppressMessages(library(minfi))

targets <- load_targets()
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)
targets <- targets[targets$DNA_Source %in% c("Buffy_Coat", "PBMC", "Saliva"), ]
sqc     <- read_sex_qc()

## probe sets ----
## The autosomal reference is a fixed random 50,000 autosomal probes: a median over that many
## is indistinguishable from one over all of them and keeps each batch's matrices small.
anno  <- array_annotation()
chr   <- as.character(anno$chr)
x_ids <- rownames(anno)[chr == "chrX"]
y_ids <- rownames(anno)[chr == "chrY"]
xci_ids <- rownames(anno)[chr == "chrX" & anno$Relation_to_Island %in% "Island"]
auto  <- rownames(anno)[chr %in% paste0("chr", 1:22)]
set.seed(20261005)
a_ids <- sample(auto, min(50000L, length(auto)))
cat(sprintf("probes: chrX %d (CpG island %d) | chrY %d | autosomal reference %d\n",
            length(x_ids), length(xci_ids), length(y_ids), length(a_ids)))

## per-batch intensity summaries ----
batch_stats <- function(tg) {
    rg <- read_rgSet_batch(tg)
    ms <- preprocessRaw(rg)
    rn <- rownames(ms)
    need <- unique(c(a_ids, x_ids, y_ids))
    keep <- rn %in% need
    Mth <- getMeth(ms)[keep, , drop = FALSE]; Un <- getUnmeth(ms)[keep, , drop = FALSE]
    rm(ms); gc()
    CN  <- log2(pmax(Mth + Un, 1))
    rows <- function(ids) which(rownames(CN) %in% ids)
    xci  <- rows(xci_ids)
    beta <- Mth[xci, , drop = FALSE] / (Mth[xci, , drop = FALSE] + Un[xci, , drop = FALSE] + 100)
    detP <- minfi::detectionP(rg)
    yd   <- which(rownames(detP) %in% y_ids)
    out <- data.frame(
        Sample     = colnames(rg),
        aMed       = matrixStats::colMedians(CN, rows = rows(a_ids), na.rm = TRUE),
        xMed       = matrixStats::colMedians(CN, rows = rows(x_ids), na.rm = TRUE),
        yMed       = matrixStats::colMedians(CN, rows = rows(y_ids), na.rm = TRUE),
        Y_detected = colMeans(detP[yd, , drop = FALSE] < DETP_THRESHOLD, na.rm = TRUE),
        XCI_beta   = colMeans(beta, na.rm = TRUE),
        stringsAsFactors = FALSE)
    rm(rg, Mth, Un, CN, beta, detP); gc()
    out
}
n   <- nrow(targets)
idx <- split(seq_len(n), ceiling(seq_len(n) / max(1L, SEX_INT_CHUNK)))
cat("sex_intensity:", n, "samples in", length(idx), "batch(es) of up to", SEX_INT_CHUNK, "\n")
res <- do.call(rbind, lapply(seq_along(idx), function(b) {
    cat("  batch", b, "/", length(idx), "\n"); utils::flush.console()
    batch_stats(targets[idx[[b]], , drop = FALSE])
}))

## calls and patterns ----
res$x_rel   <- res$xMed - res$aMed
res$y_rel   <- res$yMed - res$aMed
res$yMinusX <- res$yMed - res$xMed
res$Sex_int <- ifelse(res$yMinusX < -2, "F", "M")
if (length(unique(res$Sex_int)) < 2)
    stop("stage6/sex_intensity.R: every sample is called ", res$Sex_int[1],
         " from intensities; the X/Y reference groups cannot be formed.")
typ <- function(v) sex_typical(v, v[res$Sex_int == "F"], v[res$Sex_int == "M"], SEX_ATYPICAL_Z)
res$X_class   <- typ(res$x_rel)
res$Y_class   <- typ(res$y_rel)
res$XCI_class <- typ(res$XCI_beta)
res$pattern <- with(res, ifelse(is.na(X_class) | is.na(Y_class) | is.na(XCI_class), NA_character_,
    ifelse(X_class == "F" & Y_class == "F" & XCI_class == "F", "XX",
    ifelse(X_class == "M" & Y_class == "M" & XCI_class == "M", "XY",
    ifelse(Y_class == "M" & X_class == "F" & XCI_class == "F", "XXY-like",
    ifelse(Y_class == "F" & X_class == "M" & XCI_class == "M", "X0-like", "atypical"))))))

## context: admin sex, stage-1 QC, PCA call ----
qi <- match(res$Sample, sqc$Sample)
res$Subject_ID <- sqc$Subject_ID[qi]
res$Sex_admin  <- sqc$Sex_admin[qi]
failed_f <- file.path(INTERMEDIATE_DIR, "methylation_data_detP.failedsamp.txt")
failed   <- if (file.exists(failed_f)) read.delim(failed_f, colClasses = "character")$Sample_Group else character()
res$stage1_qc_failed <- res$Sample %in% failed
scores_f <- file.path(REPORT_DIR, "pca_sex_scores.csv")
if (file.exists(scores_f)) {
    sc <- read.csv(scores_f, colClasses = "character")
    res$Sex_pca <- sc$Sex_pca[match(res$Sample, sc$Sample)]
} else {
    cat("sex_intensity: ", scores_f, " not found; run pca_sex_batch.R for the PCA agreement\n")
    res$Sex_pca <- NA_character_
}
num <- c("aMed", "xMed", "yMed", "x_rel", "y_rel", "yMinusX", "XCI_beta", "Y_detected")
res[num] <- lapply(res[num], round, 4)
res <- res[, c("Sample", "Subject_ID", "Sex_admin", "Sex_int", "Sex_pca", "pattern", "stage1_qc_failed",
               "xMed", "yMed", "yMinusX", "aMed", "x_rel", "y_rel", "XCI_beta", "Y_detected",
               "X_class", "Y_class", "XCI_class")]
write.csv(res, file.path(REPORT_DIR, "sex_intensity.csv"), row.names = FALSE)

## agreement ----
## One block per comparison: the 2 x 2 cross-tabulation, then the agreement count and rate.
compare <- function(a, b, name) {
    ok <- !is.na(a) & !is.na(b)
    tb <- as.data.frame(table(first = factor(a[ok], c("F", "M")), second = factor(b[ok], c("F", "M"))),
                        stringsAsFactors = FALSE)
    data.frame(comparison = name, first = c(tb$first, "agree"), second = c(tb$second, "of"),
               n = c(tb$Freq, sum(a[ok] == b[ok])), total = sum(ok),
               rate = c(rep(NA, nrow(tb)), round(mean(a[ok] == b[ok]), 5)))
}
conc <- rbind(compare(res$Sex_pca, res$Sex_int, "PCA (first) vs intensity (second)"),
              compare(res$Sex_admin, res$Sex_int, "admin (first) vs intensity (second)"))
write.csv(conc, file.path(REPORT_DIR, "sex_intensity_concordance.csv"), row.names = FALSE)
pa <- conc[conc$first == "agree", ]
cat(sprintf("intensity sex agrees with PCA sex %d/%d (%.2f%%), with admin sex %d/%d (%.2f%%) | patterns: %s\n",
            pa$n[1], pa$total[1], 100 * pa$rate[1], pa$n[2], pa$total[2], 100 * pa$rate[2],
            paste(names(table(res$pattern)), table(res$pattern), sep = " ", collapse = ", ")))

## plots ----
flag <- with(res, (!is.na(Sex_admin) & Sex_int != Sex_admin) | !pattern %in% c("XX", "XY") |
                  (!is.na(Sex_pca) & Sex_pca != Sex_int))
col_admin <- c(F = "firebrick", M = "steelblue")[res$Sex_admin]; col_admin[is.na(col_admin)] <- "grey50"
lab <- ifelse(is.na(res$Subject_ID), res$Sample, res$Subject_ID)
pdf(file.path(REPORT_DIR, "sex_intensity.pdf"), height = 6, width = 13)
par(mfrow = c(1, 3), oma = c(2.5, 0, 0, 0))
panel <- function(x, y, xl, yl, ttl) {
    plot(x, y, col = col_admin, pch = 19, cex = 0.5, xlab = xl, ylab = yl, main = ttl)
    if (any(flag)) {
        points(x[flag], y[flag], pch = 1, cex = 1.6)
        text(x[flag], y[flag], lab[flag], pos = 4, cex = 0.6)
    }
}
panel(res$xMed, res$yMed, "xMed (log2 intensity)", "yMed (log2 intensity)", "getSex: chrX vs chrY")
abline(a = -2, b = 1, lty = 2)
panel(res$x_rel, res$y_rel, "chrX - autosomes", "chrY - autosomes", "Copy number relative to autosomes")
panel(res$XCI_beta, res$y_rel, "mean beta, chrX CpG islands", "chrY - autosomes", "X inactivation vs Y")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", c("admin F", "admin M", "no admin sex", "flagged"), pch = c(19, 19, 19, 1),
       col = c("firebrick", "steelblue", "grey50", "black"), bty = "n", horiz = TRUE, cex = 0.9)
dev.off()
cat("stage6/sex_intensity: wrote sex_intensity.csv, sex_intensity_concordance.csv, sex_intensity.pdf to",
    REPORT_DIR, "\n")
