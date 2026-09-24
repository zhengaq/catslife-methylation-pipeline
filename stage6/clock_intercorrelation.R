### stage6/clock_intercorrelation.R: sensitivity check on how the 15 clocks relate to one another.
### Builds the 15 x 15 correlation matrix of the clock values and of their age accelerations
### (dnaMethyAge's residual of each clock on chronological age), then asks whether the matrix holds
### up when
###   - Spearman replaces Pearson (outliers, non-linearity);
###   - the sample is cut to one row per person, then one per family (repeat waves, twins and
###     siblings are not independent);
###   - each wave is analysed alone (needs a Wave or Subject_ID column).
### Every subset uses the same rows for all 15 clocks: samples with all 15 values (complete cases).
### 95% CIs for the full-sample Pearson matrix come from a bootstrap that resamples whole families.
### Base R only; reads DERIVED_DIR/mAge_clocks.csv.
### Outputs: SENS_DIR/clock_intercorrelation_{value,acceleration}_matrix.csv (Pearson, full sample),
###          SENS_DIR/clock_intercorrelation_long.csv (every pair x measure x method x subset),
###          SENS_DIR/clock_intercorrelation_sensitivity.csv (shift of each subset vs the full sample),
###          REPORT_DIR/clock_intercorrelation.pdf (heatmaps + full-sample vs one-per-family scatter)

source("config.R")

CLOCKS_FILE <- file.path(DERIVED_DIR, "mAge_clocks.csv")
if (!file.exists(CLOCKS_FILE))
    stop("stage6/clock_intercorrelation.R: missing ", CLOCKS_FILE, "; run stage 5 first.")
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
BOOT_B <- as.integer(Sys.getenv("METHYL_INTERCOR_BOOT", "1000"))   # family-bootstrap draws; 0 = no CIs

m      <- read.csv(CLOCKS_FILE, check.names = FALSE, stringsAsFactors = FALSE)
idcol  <- if ("aid"    %in% names(m)) "aid"    else "IndividualID"
famcol <- if ("pfamid" %in% names(m)) "pfamid" else "FamilyID"
if ("clock_excluded" %in% names(m)) m <- m[!(m$clock_excluded %in% TRUE), , drop = FALSE]

## Clock value columns (as in validity_clocks.R); stage 5 writes each clock's acceleration column
## immediately after its value column.
clocks <- grep("(_mAge|_mitoticdivisions|Dunedin_Pace)$", names(m), value = TRUE)
accel  <- names(m)[match(clocks, names(m)) + 1L]
if (length(clocks) < 2 || anyNA(accel) || !all(grepl("Acceleration$", accel)))
    stop("stage6/clock_intercorrelation.R: expected each clock column to be followed by its ",
         "*_Acceleration column in ", basename(CLOCKS_FILE))
LABELS <- c(Horvath_mAge = "Horvath", Hannum_mAge = "Hannum", Horvath2_mAge = "Horvath2",
            ZhangQ_mAge = "ZhangQ", PhenoAge_mAge = "PhenoAge", PCGrimAge_mAge = "PCGrimAge",
            Dunedin_Pace = "DunedinPACE", ZhangY_mAge = "ZhangY", LuA_mAge = "DNAmTL",
            epiTOC_mitoticdivisions = "epiTOC", epiTOC2_mitoticdivisions = "epiTOC2",
            ShirebyG2020_mAge = "Shireby", PedBE_mAge = "PedBE", PanM2_mAge = "PanM2", PanM3_mAge = "PanM3")
lab <- function(x) ifelse(x %in% names(LABELS), LABELS[x], x)

m$.wave <- if ("Wave" %in% names(m)) m$Wave else
           if ("Subject_ID" %in% names(m)) subject_wave(m$Subject_ID) else NA
## A missing family id would drop the row from the family bootstrap; give it its own cluster.
m$.fam  <- ifelse(is.na(m[[famcol]]) | m[[famcol]] == "", paste0(".row", seq_len(nrow(m))), m[[famcol]])

MEASURES <- list(value = clocks, acceleration = accel)
num_mat <- function(d, cols) {
    x <- vapply(cols, function(cl) suppressWarnings(as.numeric(d[[cl]])), numeric(nrow(d)))
    x <- matrix(x, nrow = nrow(d)); colnames(x) <- clocks   # both measures keyed by clock name
    x
}
one_per <- function(d, key) {                 # one random row per key level (as in validity_clocks.R)
    k <- d[[key]]; d <- d[!is.na(k) & k != "", , drop = FALSE]
    s <- split(seq_len(nrow(d)), d[[key]])
    d[vapply(s, function(ix) if (length(ix) == 1) ix else sample(ix, 1), integer(1)), , drop = FALSE]
}
## Subsets are drawn from the complete-case rows, so a person whose one sample lacks a clock is
## represented by their other sample rather than lost.
subsets_of <- function(d) {
    set.seed(123)
    person <- one_per(d, idcol)
    out <- list(all = d, one_per_person = person, one_per_family = one_per(person, famcol))
    for (w in sort(unique(stats::na.omit(d$.wave)))) out[[paste0("wave", w)]] <- d[d$.wave %in% w, , drop = FALSE]
    out
}
UT <- upper.tri(diag(length(clocks)))
UT_IX <- which(UT, arr.ind = TRUE)            # column-major, same order as r[UT]
pairs_long <- function(r, n, measure, method, subset)
    data.frame(measure = measure, method = method, subset = subset,
               clock_1 = clocks[UT_IX[, 1]], clock_2 = clocks[UT_IX[, 2]], r = r[UT], n = n)
family_boot <- function(x, fam, B) {
    fams <- split(seq_len(nrow(x)), fam)
    draws <- replicate(B, {
        ix <- unlist(fams[sample.int(length(fams), replace = TRUE)], use.names = FALSE)
        suppressWarnings(cor(x[ix, , drop = FALSE]))[UT]
    })
    apply(draws, 1, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
}

long <- list(); full <- list()
cat("=== stage 6: clock inter-correlation on", basename(CLOCKS_FILE), "(", nrow(m), "rows,",
    length(clocks), "clocks ) ===\n")
for (ms in names(MEASURES)) {
    x_all <- num_mat(m, MEASURES[[ms]])
    d <- m[stats::complete.cases(x_all), , drop = FALSE]
    subs <- subsets_of(d)
    cat(sprintf("%-12s complete cases %d | %s\n", ms, nrow(d),
                paste(sprintf("%s n=%d", names(subs), vapply(subs, nrow, integer(1))), collapse = ", ")))
    for (sn in names(subs)) {
        x <- num_mat(subs[[sn]], MEASURES[[ms]])
        for (meth in c("pearson", "spearman")) {
            r <- suppressWarnings(cor(x, method = meth))
            row <- pairs_long(r, nrow(x), ms, meth, sn)
            row$ci_lo <- NA_real_; row$ci_hi <- NA_real_
            if (sn == "all" && meth == "pearson") {
                full[[ms]] <- r
                if (BOOT_B > 0) {
                    set.seed(123)
                    ci <- family_boot(x, subs$all$.fam, BOOT_B)
                    row$ci_lo <- ci[1, ]; row$ci_hi <- ci[2, ]
                }
            }
            long[[length(long) + 1]] <- row
        }
    }
}
long <- do.call(rbind, long)
long$r <- round(long$r, 4); long$ci_lo <- round(long$ci_lo, 4); long$ci_hi <- round(long$ci_hi, 4)
write.csv(long, file.path(SENS_DIR, "clock_intercorrelation_long.csv"), row.names = FALSE)
for (ms in names(full)) {
    r <- round(full[[ms]], 3); dimnames(r) <- list(lab(clocks), lab(clocks))
    write.csv(r, file.path(SENS_DIR, paste0("clock_intercorrelation_", ms, "_matrix.csv")))
}

## Sensitivity summary: how far each subset / method moves the 105 pairwise correlations away from
## the full-sample Pearson matrix of the same measure.
key  <- paste(long$measure, long$clock_1, long$clock_2)
base <- long[long$subset == "all" & long$method == "pearson", ]
long$r_full <- base$r[match(key, paste(base$measure, base$clock_1, base$clock_2))]
grp  <- unique(long[, c("measure", "method", "subset")])
grp  <- grp[!(grp$subset == "all" & grp$method == "pearson"), ]
sens <- do.call(rbind, lapply(seq_len(nrow(grp)), function(i) {
    g <- long[long$measure == grp$measure[i] & long$method == grp$method[i] & long$subset == grp$subset[i], ]
    d <- g$r - g$r_full
    data.frame(grp[i, ], n = g$n[1], mean_abs_shift = round(mean(abs(d)), 3),
               max_abs_shift = round(max(abs(d)), 3), pairs_shift_gt_0.1 = sum(abs(d) > 0.1),
               largest_shift_pair = paste(lab(g$clock_1[which.max(abs(d))]), lab(g$clock_2[which.max(abs(d))]), sep = " / "))
}))
write.csv(sens, file.path(SENS_DIR, "clock_intercorrelation_sensitivity.csv"), row.names = FALSE)

## Figures ----
heat <- function(r, ord, title) {
    r <- r[ord, ord]; k <- ncol(r)
    pal <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
    op <- par(mar = c(7, 7, 3, 1)); on.exit(par(op))
    z <- t(r)[, k:1]                           # row 1 at the top
    image(1:k, 1:k, z, zlim = c(-1, 1), col = pal, axes = FALSE, xlab = "", ylab = "", main = title)
    axis(1, at = 1:k, labels = lab(colnames(r)), las = 2, cex.axis = 0.8, tick = FALSE)
    axis(2, at = 1:k, labels = rev(lab(rownames(r))), las = 1, cex.axis = 0.8, tick = FALSE)
    for (i in 1:k) for (j in 1:k)
        text(i, j, sprintf("%.2f", z[i, j]), cex = 0.55, col = if (abs(z[i, j]) > 0.6) "white" else "black")
}
## One ordering for both heatmaps (clustered on the acceleration matrix, 1 - r distance) so the
## two pages can be compared cell by cell.
ord <- stats::hclust(stats::as.dist(1 - full$acceleration), method = "average")$order
n_of <- function(ms) base$n[base$measure == ms][1]
grDevices::pdf(file.path(REPORT_DIR, "clock_intercorrelation.pdf"), width = 8.5, height = 8)
heat(full$value,        ord, sprintf("Clock values: Pearson r (n = %d)", n_of("value")))
heat(full$acceleration, ord, sprintf("Age acceleration: Pearson r (n = %d)", n_of("acceleration")))
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (ms in names(MEASURES)) {
    g <- long[long$measure == ms & long$method == "pearson" & long$subset == "one_per_family", ]
    ttl <- c(value = "Clock values", acceleration = "Age acceleration")[[ms]]
    plot(g$r_full, g$r, xlim = c(-1, 1), ylim = c(-1, 1), pch = 19, cex = 0.6,
         xlab = "Full sample r", ylab = "One per family r", main = ttl)
    abline(0, 1, lty = 3); abline(h = 0, v = 0, col = "grey80")
}
par(op)
invisible(grDevices::dev.off())

## Console summary ----
off <- function(r) r[UT]
for (ms in names(full)) {
    r <- off(full[[ms]])
    cat(sprintf("%-12s median |r| %.2f (range %.2f to %.2f)\n", ms, median(abs(r)), min(r), max(r)))
}
acc <- base[base$measure == "acceleration", ]
acc <- acc[order(-acc$r), ]
cat("acceleration, strongest pairs:", paste(sprintf("%s/%s %.2f", lab(head(acc$clock_1, 3)), lab(head(acc$clock_2, 3)), head(acc$r, 3)), collapse = "; "), "\n")
worst <- sens[which.max(sens$max_abs_shift), ]
cat(sprintf("sensitivity: largest shift vs full-sample Pearson = %.3f (%s, %s, %s: %s)\n",
            worst$max_abs_shift, worst$measure, worst$method, worst$subset, worst$largest_shift_pair))
if ("LuA_mAge" %in% clocks)
    cat("note: DNAmTL (LuA_mAge) is telomere length; negative correlations with the age clocks are expected\n")
cat("stage6/clock_intercorrelation: wrote 4 tables to", SENS_DIR, "and clock_intercorrelation.pdf to", REPORT_DIR, "\n")
