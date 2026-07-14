### Script to QC and process CATSLife Methylation data
### following https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylation_tutorials.html
### also uses https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
### and https://www.rdocumentation.org/packages/minfi/versions/1.18.4
###
### Dependencies come from the R environment set up per the README; this script no
### longer installs packages. Paths, filenames and QC thresholds come from config.R.

cat("Beginning analyses\n", date(), "\n\n")

source("config.R")
library("EpiDISH")
library("minfi")
library("wateRmelon")

########################################################################################################
### Sample annotations
########################################################################################################
targets <- load_targets()


########################################################################################################
### 1. Read methylation data  (IDATs named in the sample sheet)
########################################################################################################
rgSet <- load_raw_rgSet()
save(rgSet, file = F_RAW)

### Information on the rgSet
pd <- pData(rgSet)
pd$DNA_Source <- canonicalize_dna_source(pd$DNA_Source)  ### canonicalize ONCE, here, before anything derives from it
tis <- data.frame(as.factor(pd$DNA_Source), as.numeric(as.factor(pd$DNA_Source))); colnames(tis) <- c("tissue", "tnumeric")
ctl.samp <- which(tis$tissue == "Cell_Line")
cat("Completed reading in all the array data\n", date(), "\n\n")


########################################################################################################
### 2. FILTER POOR-PERFORMING PROBES and SAMPLES. REQUIRE PROBE MISSINGNESS <0.01 using a detP>0.05 threshold
### Calculate the detection p-values for probes and individuals.
### Requires mapping probes to genome/annotation, then filtering.
########################################################################################################
detP <- minfi::detectionP(rgSet)   ### qualify: ewastools (loaded in stage 3) also exports detectionP and would mask minfi's in a shared session
save(detP, file = F_DETP)
str(detP)
cat("completed detP calculations\ndetP dimensions:\n")
cat(dim(detP), "\n")
cat("does the order of the samples match between rgSet and pd/detP:\n")
print(all(colnames(rgSet) == pd$Sample_Group)) ### Sample order still matches
print(all(colnames(detP)  == pd$Sample_Group)) ### Sample order still matches


# examine mean detection p-values across all samples to identify any failed samples
samp.missingness.thresh <- SAMPLE_MISSINGNESS
samp.mean.detP <- colMeans(detP)
m <- nrow(detP)
samp.missing <- colSums(detP < DETP_THRESHOLD)
keep.samp <- samp.missing >= (1 - samp.missingness.thresh) * m
cat("sample QC: 1% probe missingness threshold at SAMPLE LEVEL\nremove keep\n")
cat(table(keep.samp), "\n")
rm.samp <- which(!keep.samp)
if (length(rm.samp) > 0) {
    print(cbind(names(rm.samp), tis[which(!keep.samp), ]))
    write.table(pd[rm.samp, ], file.path(ANALYSIS_DIR, "methylation_data_detP.failedsamp.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
}

# examine mean detection p-values across all probes to identify any failed probes
probe.missingness.thresh <- PROBE_MISSINGNESS
n <- ncol(detP)
probe.missing <- rowSums(detP < DETP_THRESHOLD)
keep.probe <- probe.missing >= (1 - probe.missingness.thresh) * n
cat("sample QC: 1% missingness threshold at PROBE LEVEL\nremove keep\n")
cat(table(keep.probe), "\n")
rm.probe <- which(!keep.probe)
if (length(rm.probe) > 0) {
    write.table(as.matrix(names(rm.probe), ncol = 1), file.path(ANALYSIS_DIR, "methylation_data_detP.failedprobe.txt"), row.names = F, col.names = T, sep = "\t", quote = F)
}


pdf(file.path(REPORT_DIR, "minfi_QC.missingness_detP.pdf"), height = 4, width = 8)
    par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
    hist(1 - samp.missing / m, breaks = 40, main = 'Sample Missingness (detP>0.05)', xlab = "Missingness Per Sample"); abline(v = samp.missingness.thresh, col = "red")
    text(0.01, 100, paste0("Probe detP>0.05\n", samp.missingness.thresh, " Sample Missingness\nexclude\tretain\n", table(keep.samp)[1], "\t", table(keep.samp)[2]), pos = 4)
    hist(1 - probe.missing / n, breaks = 100, main = 'Probe Missingness (detP>0.05)', xlab = "Missingness Per Probe"); abline(v = samp.missingness.thresh, col = "red")
    text(0.2, 4e5, paste0("Probe detP>0.05\n", probe.missingness.thresh, " Probe Missingness\nexclude\tretain\n", table(keep.probe)[1], "\t", table(keep.probe)[2]), pos = 4)
dev.off()


### Sample-missingness table (written first, so the by-tissue comparison can use it in-memory)
samp.miss.prop <- 1 - samp.missing / m
miss.data <- cbind(samp.miss.prop, pd)
write.table(miss.data, file.path(REPORT_DIR, "sample_missingness.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

### Compare missingness across tissues:
d <- as.data.frame(miss.data)
hist.all <- hist(d$samp.miss.prop, breaks = 100, plot = FALSE)
fit <- lm(d$samp.miss.prop ~ d$DNA_Source)
print(summary(fit))
jpeg(file.path(REPORT_DIR, "Sample_missingness_tissue.jpg"), width = 6, height = 5, units = 'in', res = 400)
hist(d$samp.miss.prop[which(d$DNA_Source == "PBMC")],       breaks = hist.all$breaks, col = rgb(1, 0, 0, .2), border = rgb(1, 0, 0, .2), ylim = c(0, 60), main = "Sample Missingness by Tissue", xlab = "Missingness Per Sample")
hist(d$samp.miss.prop[which(d$DNA_Source == "Saliva")],     breaks = hist.all$breaks, col = rgb(0, 0, 1, .2), border = rgb(0, 0, 1, .2), add = T)
hist(d$samp.miss.prop[which(d$DNA_Source == "Buffy_Coat")], breaks = hist.all$breaks, col = rgb(0, 1, 0, .2), border = rgb(0, 1, 0, .2), add = T)
abline(v = 0.01, col = 2)
legend("topright", c("PBMC", "Saliva", "Buffy_Coat"), col = c(rgb(1, 0, 0, .2), rgb(0, 0, 1, .2), rgb(0, 1, 0, .2)), pch = 15, bty = 'n')
dev.off()


### Map rgSet to genome to get probeIDs, then subset: remove failed samples, failed probes, and probes on SNPs
grgSet <- mapToGenome(rgSet)        ##  Mapping to genome. REMOVES a few SNPs

### Fail loud on an array-version mismatch instead of silently mis-mapping/dropping probes
### under the wrong manifest (minfi resolves the manifest from the IDAT header, not a flag).
resolved.array <- annotation(grgSet)["array"]
cat("Array annotation resolved to:", resolved.array, annotation(grgSet)["annotation"], "\n")
if (ARRAY_VERSION == "v2" && !grepl("EPICv2", resolved.array))
    stop("ARRAY_VERSION=v2 but minfi resolved the array annotation to '", resolved.array,
         "' — check the raw IDATs / installed manifest packages.")
if (ARRAY_VERSION == "v1" && grepl("EPICv2", resolved.array))
    stop("ARRAY_VERSION=v1 but minfi resolved the array annotation to '", resolved.array,
         "' — the IDATs look like EPIC v2; set METHYL_ARRAY_VERSION=v2.")

grgSet.ns <- dropLociWithSnps(grgSet)  ### DROPPING SITES ON SNPS

rm.probes.rgSet <- which(rownames(grgSet.ns) %in% names(rm.probe))
rgSetflt <- grgSet.ns[keep_idx(nrow(grgSet.ns), rm.probes.rgSet), keep_idx(ncol(grgSet.ns), rm.samp)]

### Also remove failed samples from targets/pd/tis:
cat("does the order of the samples match between grgSet.ns and pd/tis/targets:\n")
print(all(colnames(grgSet.ns) == targets$Sample_Group))
pd      <- pd[keep_idx(nrow(pd), rm.samp), ]
targets <- targets[keep_idx(nrow(targets), rm.samp), ]
tis     <- tis[keep_idx(nrow(tis), rm.samp), ]

rgSetflt
cat("are there any probes that are in the remove list still in the filtered rgSetflt?:\n")
print(any(rownames(rgSetflt) %in% names(rm.probe))) ### Should be FALSE

save(rgSetflt, file = F_RGFLT)
cat("Completed filtering probes and samples on detP\n", date(), "\n\n")


########################################################################################################
### 3. Background correction ('noob') and normalize with 'dasen' method (from wateRmelon)
########################################################################################################
rgSet.rmsamp <- rgSet[, keep_idx(ncol(rgSet), rm.samp)] ### RGChannelSet with failed samples removed
MSetNoob <- preprocessNoob(rgSet.rmsamp, verbose = TRUE) ### preprocessNoob needs an RGChannelSet (problem probes removed later)
save(MSetNoob, file = F_NOOB)
MSetNoob
cat("Completed noob background correction\n", date(), "\n\n")

### Remove failed probes from the detP tests here:
MSetNoob.flt <- MSetNoob
if (any(rownames(MSetNoob) %in% names(rm.probe))) {
    rm.probes.MSetNoob <- which(rownames(MSetNoob) %in% names(rm.probe))
    MSetNoob.flt <- MSetNoob[-rm.probes.MSetNoob, ]
}
print(any(rownames(MSetNoob.flt) %in% names(rm.probe))) ### FALSE => failed-probe removal worked
save(MSetNoob.flt, file = F_NOOBFLT)


### Normalize using the dasen method in wateRmelon
dasen.melon <- dasen(MSetNoob.flt)
save(dasen.melon, file = F_DASEN)
dasen.melon
str(dasen.melon)
cat("Completed dasen normalization\n", date(), "\n\n")


### Plot the quality metrics:
qu <- qual(betas(MSetNoob.flt), betas(dasen.melon)) ### A couple of metrics between normalized & non-normalized betas
str(qu)
pdf(file.path(REPORT_DIR, "minfi_QC.wateRmelonQC.pdf"), height = 8, width = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
boxplot(qu[, 1] ~ pd$DNA_Source)
boxplot(qu[, 2] ~ pd$DNA_Source)
plot(qu[, 1], qu[, 2], col = tis$tnumeric, pch = tis$tnumeric); abline(0, 1, col = 1, lty = 3, lwd = 2)
legend('bottomright', legend = levels(tis$tissue), col = 1:4, pch = 1:4, bty = 'n', title = 'tissue')
dev.off()


print(all(colnames(dasen.melon) == pd$Sample_Group))

betas.m <- getBeta(dasen.melon)
Mvalues <- getM(dasen.melon)
### Canonicalize EPIC v2 replicate-probe IDs (strip "_TC21"-style suffixes, collapsing to the
### bare cg######## id clock lists / cell-type references expect) ONCE, here, so every
### downstream reader (stages 3-5) sees canonical ids unconditionally. No-op for v1 data.
betas.m <- canonicalize_v2_probe_ids(betas.m)
Mvalues <- canonicalize_v2_probe_ids(Mvalues)
dasen.values <- list()
dasen.values$b <- betas.m
dasen.values$M <- Mvalues
save(dasen.values, file = F_DASENB)

### Compare normalization by tissue (excluding cell-line controls):
rm.qual <- which(!rownames(pd) %in% rownames(qu))
if (length(rm.qual)) pd <- pd[-rm.qual, ]
print(all(rownames(pd) == rownames(qu)))
ctl.samp <- which(pd$DNA_Source == "Cell_Line")
pd.noctl  <- pd[keep_idx(nrow(pd), ctl.samp), ]
qu.noctl  <- qu[keep_idx(nrow(qu), ctl.samp), ]
tis.noctl <- as.factor(pd.noctl$DNA_Source)
tnumeric  <- as.numeric(tis.noctl)
jpeg(file.path(REPORT_DIR, "minfi_QC.wateRmelonQC.byTissue.jpg"), height = 6, width = 6, units = 'in', res = 400)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
boxplot(qu.noctl[, 1] ~ pd.noctl$DNA_Source, ylab = "wateRmelon qual Method rmsd", xlab = "Tissue")
boxplot(qu.noctl[, 2] ~ pd.noctl$DNA_Source, ylab = "wateRmelon qual Method sdd", xlab = "Tissue")
plot(qu.noctl[, 1], qu.noctl[, 2], col = tnumeric, pch = tnumeric, ylab = "wateRmelon qual Method sdd", xlab = "wateRmelon qual Method rmsd"); abline(0, 1, col = 1, lty = 3, lwd = 2)
legend('bottomright', legend = unique(pd.noctl$DNA_Source), col = 1:3, pch = 1:3, bty = 'n', title = 'tissue')
dev.off()

cat("Completed analyses\n", date(), "\n\n")
