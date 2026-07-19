#!/usr/bin/env Rscript
### Estimate cell-type proportions, residualize betas on them, then on plate batch.
### Blood (PBMC + Buffy Coat) and Saliva are handled separately, each with its own
### reference. Chunked by CpG range for SLURM parallelization. Example (5 chunks):
###   Rscript 3.methylation_adjust.chunked.R --nparts 5 --part $i      # i = 1..5
### then merge the chunks with 4.methylation_merge.chunked.R.

suppressMessages(library("optparse"))
cat("Beginning analyses\n", date(), "\n\n")

source("config.R")

# command-line options
option_list <- list(
  make_option("--nparts", action = "store", default = NPARTS, type = 'numeric', help = "how many parts to chunk?"),
  make_option("--part",   action = "store", default = 1,      type = 'numeric', help = "the chunk being run")
)
opt <- parse_args(OptionParser(option_list = option_list))
str(opt)

suppressMessages({
    library("minfi"); library("wateRmelon"); library("EpiDISH")
    library("FlowSorted.Blood.EPIC"); library("BeadSorted.Saliva.EPIC")
    library(ExperimentHub); library(ewastools)
})

## one-line TSV writer (no trailing tab — avoids a phantom empty column on read-back)
wl <- function(vals, file, append = TRUE) cat(paste(vals, collapse = "\t"), "\n", file = file, sep = "", append = append)


########################################################################################################
### Sample annotations
########################################################################################################
targetscsv <- read_sample_sheet(SAMPLE_SHEET)
targetscsv$DNA_Source <- canonicalize_dna_source(targetscsv$DNA_Source)
targets    <- load_targets()
targets    <- targets[, c("Sample_Group", "DNA_Source")]   ### Array is a minfi pData column, not in the raw sheet
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
if (all(dasen.IDs[, 1] != IDs[, 1])) {
    cat("IDs don't match!!!\n", file = stderr()); quit(save = 'no')
}


########################################################################################################
### ESTIMATE CELL TYPE PROPORTIONS FOR BLOOD AND SALIVA SEPARATELY
########################################################################################################
cat("\nStarting Cell Type Proportion Estimation\n", date(), "\n")

saliva.samp <- which(IDs$DNA_Source == "Saliva")
blood.samp  <- which(IDs$DNA_Source %in% c("Buffy_Coat", "PBMC"))
has_blood  <- length(blood.samp)  > 0L
has_saliva <- length(saliva.samp) > 0L
if (!has_blood && !has_saliva)
    stop("stage 3: no blood (Buffy_Coat/PBMC) or Saliva samples in DNA_Source; nothing to adjust")
cat(sprintf("Tissue tracks present: blood n=%d, saliva n=%d\n", length(blood.samp), length(saliva.samp)))

b.blood  <- if (has_blood)  dasen.values$b[, blood.samp,  drop = FALSE] else NULL
b.saliva <- if (has_saliva) dasen.values$b[, saliva.samp, drop = FALSE] else NULL
rm(dasen.values); gc()   ### only the per-tissue beta submatrices are used from here; free the full betas/M list

### estimate cell proportions (only for the tissue tracks actually present)
data(centEpiFibIC.m)
data(centDHSbloodDMC.m)
cp.list <- list()
if (has_blood) {
    cellprop.blood <- epidish(b.blood, as.matrix(centDHSbloodDMC.m), method = "RPC")
    cp.b <- cbind(as.character(colnames(b.blood)), cellprop.blood$estF); colnames(cp.b)[1] <- "Sample_Group"
    cp.list$blood <- cp.b
}
if (has_saliva) {
    cellprop.saliva <- estimateLC(b.saliva, ref = "salivaEPIC", constrain = TRUE)
    cp.s <- cbind(as.character(colnames(b.saliva)), cellprop.saliva); colnames(cp.s)[1] <- "Sample_Group"
    cp.list$saliva <- cp.s
}
cp.all <- if (length(cp.list) == 2L) merge(cp.list$blood, cp.list$saliva, all = TRUE) else cp.list[[1]]
cp.IDs <- merge(cp.all, IDs)
write.table(cp.IDs, file = CELL_PROPORTIONS_FILE, row.names = F, col.names = T, sep = "\t", quote = F)
cat("Completed Cell Type Proportion Estimation\n", date(), "\n\n")


########################################################################################################
### Residualize Beta values on cell type proportions, for saliva & blood separately
########################################################################################################
cat("Residualizing B values\n")
MBmat.blood  <- b.blood
MBmat.saliva <- b.saliva

cat("\nStarting B-value residualization\n", date(), "\n\n")
### When both tissues are present their CpG rownames must match (fail loud); test.sites = the shared CpG set
if (has_blood && has_saliva && !all(rownames(MBmat.blood) == rownames(MBmat.saliva))) {
    cat("**** ROWNAMES DO NOT MATCH BETWEEN BLOOD AND SALIVA\nQUITTING NOW", date(), "\n"); quit(save = "no")
}
test.sites <- rownames(if (has_blood) MBmat.blood else MBmat.saliva)

### Optionally cap the per-CpG loop for speed (RESID_CPG_LIMIT=0 means all CpGs).
if (RESID_CPG_LIMIT > 0 && length(test.sites) > RESID_CPG_LIMIT) {
    cat("Capping residualization to first", RESID_CPG_LIMIT, "CpGs\n")
    test.sites <- test.sites[seq_len(RESID_CPG_LIMIT)]
}

### Set up start and end CpGs for this chunk:
ntests <- length(test.sites)
tstart <- (opt$part - 1) * floor(ntests / opt$nparts) + 1
tend   <- opt$part * floor(ntests / opt$nparts)
if (opt$part == opt$nparts) tend <- ntests
cat("Starting residualization for\nchunk\tstart\tend\n")
cat(opt$part, "\t", tstart, "\t", tend, "\n")

blood.resid.file  <- file.path(REPORT_DIR, paste0("B.residualized.blood.",  opt$part, ".txt"))
saliva.resid.file <- file.path(REPORT_DIR, paste0("B.residualized.saliva.", opt$part, ".txt"))

### Adjust blood samples using monocytes as the (omitted) reference cell type
if (has_blood) {
    wl(c("CpG", colnames(MBmat.blood)), blood.resid.file, append = FALSE)
    for (i in tstart:tend) {
        fit <- lm(MBmat.blood[i, ] ~ cellprop.blood$estF[, 'B'] + cellprop.blood$estF[, 'NK'] + cellprop.blood$estF[, 'CD4T'] + cellprop.blood$estF[, 'CD8T'] + cellprop.blood$estF[, 'Eosino'] + cellprop.blood$estF[, 'Neutro'])
        wl(c(test.sites[i], fit$resid + coef(fit)[1]), blood.resid.file)
    }
    cat("Done with blood adjustment\n", date(), "\n")
}

### Adjust saliva samples by epithelial cell proportion:
if (has_saliva) {
    wl(c("CpG", colnames(MBmat.saliva)), saliva.resid.file, append = FALSE)
    for (i in tstart:tend) {
        fit <- lm(MBmat.saliva[i, ] ~ cellprop.saliva$Epithelial.cells)
        wl(c(test.sites[i], fit$resid + coef(fit)[1]), saliva.resid.file)
    }
    cat("Done with saliva adjustment\n", date(), "\n")
}
cat("Completed chunk", opt$part, "of", opt$nparts, "of B residualization\n", date(), "\n")
rm(list = intersect(c("b.blood", "b.saliva", "MBmat.blood", "MBmat.saliva", "cellprop.blood", "cellprop.saliva"), ls())); gc()   ### residualized betas are on disk now; free before the plate-batch pass


########################################################################################################
### Adjust for batches (plate)
########################################################################################################
### Drop any stale adjusted file for an absent tissue so stage 4 does not merge it
for (tissue in c("blood", "saliva")[!c(has_blood, has_saliva)]) {
    stale <- file.path(REPORT_DIR, paste0("B.adjusted.regression.", tissue, ".", opt$part, ".txt"))
    if (file.exists(stale)) file.remove(stale)
}
for (tissue in c("blood", "saliva")[c(has_blood, has_saliva)]) {
    resid.file <- file.path(REPORT_DIR, paste0("B.residualized.", tissue, ".", opt$part, ".txt"))
    adj.file   <- file.path(REPORT_DIR, paste0("B.adjusted.regression.", tissue, ".", opt$part, ".txt"))
    d <- read.table(resid.file, header = TRUE, check.names = FALSE)

    rownames(d) <- d[, 1]   ### CpG names from the first column, then drop it
    d <- d[, -1]
    cat("Completed reading in data\n", date(), "\n")

    d.IDs <- as.matrix(colnames(d), ncol = 1); colnames(d.IDs) <- "Sample_Group"
    d.IDs <- strip_x_prefix(d.IDs)   ### undo R's X-prefix on numeric-leading sample names
    IDs   <- merge(d.IDs, targetscsv, sort = FALSE)
    cat("Sample names matched to the sheet without reordering?\n")
    cat(all(d.IDs[, 1] == IDs[, 1]), "\n")
    if (any(d.IDs[, 1] != IDs[, 1])) {
        cat("IDs don't match!!!\n", file = stderr()); quit(save = 'no')
    }

    wl(c("CpG", colnames(d)), adj.file, append = FALSE)
    for (cpg in 1:nrow(d)) {
        fit <- lm(as.numeric(d[cpg, ]) ~ IDs$Sample_Plate)
        wl(c(rownames(d)[cpg], resid(fit) + coef(fit)[1]), adj.file)
    }
}
cat("Completed analyses\n", date(), "\n\n")
