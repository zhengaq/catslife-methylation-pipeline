### PCA on the dasen-normalized M-values, colored by tissue. Standalone diagnostic; feeds nothing downstream.
cat("Beginning analyses\n", date(), "\n\n")
source("config.R")
suppressMessages(library("minfi"))


########################################################################################################
### Sample annotations
########################################################################################################
targets <- load_targets()
targets <- targets[, c("Sample_Group", "DNA_Source")]   ### Array is a minfi pData column, not in the raw sheet
targets$DNA_Source <- canonicalize_dna_source(targets$DNA_Source)


########################################################################################################
###  LOAD IN DATA (DASEN-NORMALIZED DATA)
########################################################################################################
dasen.values <- load_one(F_DASENB)
cat("Done loading dasen values\n", date(), "\n")

dasen.IDs <- as.matrix(colnames(dasen.values$b), ncol = 1); colnames(dasen.IDs) <- "Sample_Group"
IDs <- merge(dasen.IDs, targets, sort = FALSE)
cat("Has the ID key information been added to the Sample names, without sorting the sample names?\n")
cat(all(dasen.IDs[, 1] == IDs[, 1]), "\n")


########################################################################################################
###  RUN PCA  (default: top PCA_NCPG most-variable CpGs; METHYL_PCA_SUBSET=FALSE => all CpGs)
########################################################################################################
if (PCA_SUBSET) {
    ## Most-variable CpGs: leading PCs are driven by the highest-variance probes, so this is the
    ## standard input for a structure/QC PCA (and memory-safe vs all ~800k). See README.runtime.md.
    vary     <- matrixStats::rowVars(dasen.values$M)
    keep.cpg <- rownames(dasen.values$M)[order(vary, decreasing = TRUE)[seq_len(min(PCA_NCPG, nrow(dasen.values$M)))]]
    Mpca     <- t(dasen.values$M[keep.cpg, ])
} else {
    Mpca <- t(dasen.values$M)
}
rm(dasen.values); gc()   ### loaded betas/M no longer needed once transposed into Mpca

Mpca.s <- scale(Mpca); rm(Mpca); gc()   ### column-wise center/scale in C; apply(, 2, scale) OOMs at cohort scale
## scaling turns constant CpGs (sd 0) into non-finite columns; drop them so prcomp doesn't choke
Mpca.s <- Mpca.s[, colSums(!is.finite(Mpca.s)) == 0, drop = FALSE]

sel    <- which(IDs[, 'DNA_Source'] %in% c("Buffy_Coat", "PBMC", "Saliva"))
Mpca   <- Mpca.s[sel, , drop = FALSE]
colpch <- as.factor(IDs$DNA_Source[sel])
save(Mpca, file = file.path(ANALYSIS_DIR, "testMpca.RDat"))

cat("Starting PCA\n", date(), "\n")
pca <- prcomp(Mpca, scale. = FALSE, rank. = 10)
save(pca, file = file.path(REPORT_DIR, "dasen_Mpca_pca.RDat"))

npc <- min(10, ncol(pca$x))
imp <- summary(pca)$importance
pdf(file.path(REPORT_DIR, "PCA_plots.pdf"), height = 10, width = 10)
par(mfrow = c(3, 3), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
    plot(imp[2, 1:npc], ylim = c(0, 1), type = 'b', ylab = "Proportion of variance explained")
    points(imp[3, 1:npc], col = 2, type = 'b')
    legend("right", legend = c("Per axis", "Cumulative"), title = "Variance Explained", bty = 'n', col = 1:2, pch = 1:2)
    plot(pca$x[, 1], pca$x[, 2], col = colpch, pch = as.numeric(colpch), xlab = "PC 1", ylab = "PC 2")
    legend("topleft", legend = levels(colpch), title = "Tissue", bty = 'n', col = seq_along(levels(colpch)), pch = seq_along(levels(colpch)))
    for (pr in list(c(1, 3), c(1, 4), c(2, 3), c(3, 4), c(4, 5), c(5, 6), c(6, 7), c(7, 8))) {
        if (max(pr) <= npc)
            plot(pca$x[, pr[1]], pca$x[, pr[2]], col = colpch, pch = as.numeric(colpch),
                 xlab = paste("PC", pr[1]), ylab = paste("PC", pr[2]))
    }
dev.off()

cat("Completed PCA\n", date(), "\n\n")
