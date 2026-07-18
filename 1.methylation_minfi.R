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
### RESUME TIERS. dasen.melon can be reached three ways, cheapest-resume first:
###   (0) F_DASEN   -- the dasen checkpoint: skip everything; only betas/M extraction + save remain.
###   (1) F_NOOBFLT -- the noob checkpoint: re-run only dasen (skips read/detP/QC/mapToGenome/noob).
###   (2) fresh     -- read IDATs -> detP QC -> map + drop SNP loci -> noob -> dasen.
### MSetNoob.flt (the pre-normalization set, needed for the wateRmelon `qual` QC plots) is available
### in tiers 1/2 but NOT tier 0, so those QC plots are skipped on an F_DASEN resume.
########################################################################################################
dasen.melon <- NULL

if (RESUME && file.exists(F_DASEN)) {
    cat("RESUME: loading dasen.melon from", F_DASEN, "- only betas/M extraction remains\n")
    dasen.melon  <- load_one(F_DASEN)
    MSetNoob.flt <- NULL
} else if (RESUME && file.exists(F_NOOBFLT)) {
    cat("RESUME: loading noob checkpoint from", F_NOOBFLT, "\n")
    ckpt <- load_one(F_NOOBFLT)
    if (!is.list(ckpt) || is.null(ckpt$mset))
        stop("F_NOOBFLT is an old-format checkpoint (a bare MethylSet, predating the resume + ",
             "SNP-drop redesign). Delete it and re-run stage 1 fresh (METHYL_RESUME=FALSE) so the ",
             "SNP-locus drop is applied to the deliverable.")
    MSetNoob.flt <- ckpt$mset
    pd  <- ckpt$pd
    tis <- ckpt$tis
    rm(ckpt)
    cat("RESUME: skipping read, detP, QC, mapToGenome, dropLociWithSnps, noob\n")
} else {

    ####################################################################################################
    ### Sample annotations
    ####################################################################################################
    targets <- load_targets()

    ####################################################################################################
    ### 1. Read methylation data  (IDATs named in the sample sheet)
    ####################################################################################################
    if (RESUME && file.exists(F_RAW)) {
        cat("RESUME: loading rgSet from", F_RAW, "\n"); rgSet <- load_one(F_RAW)
    } else {
        rgSet <- load_raw_rgSet()
        if (SAVE_INTERMEDIATES) save(rgSet, file = F_RAW)
    }

    ### Information on the rgSet
    pd <- pData(rgSet)
    pd$DNA_Source <- canonicalize_dna_source(pd$DNA_Source)  ### canonicalize ONCE, here, before anything derives from it
    tis <- data.frame(as.factor(pd$DNA_Source), as.numeric(as.factor(pd$DNA_Source))); colnames(tis) <- c("tissue", "tnumeric")
    cat("Completed reading in all the array data\n", date(), "\n\n")

    ####################################################################################################
    ### 2. FILTER POOR-PERFORMING PROBES and SAMPLES. REQUIRE PROBE MISSINGNESS <0.01 using a detP>0.05 threshold
    ### Calculate the detection p-values for probes and individuals.
    ####################################################################################################
    if (RESUME && file.exists(F_DETP)) {
        cat("RESUME: loading detP from", F_DETP, "\n"); detP <- load_one(F_DETP)
    } else {
        detP <- detectionP_chunked(rgSet)   ### chunked by sample to cap peak memory (see config.R)
        if (SAVE_INTERMEDIATES) save(detP, file = F_DETP)
    }
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
    rm(detP); gc()   ### QC drop-sets computed; free detP (~10GB) before the memory-heavy steps

    pdf(file.path(REPORT_DIR, "minfi_QC.missingness_detP.pdf"), height = 4, width = 8)
        par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
        hist(1 - samp.missing / m, breaks = 40, main = 'Sample Missingness (detP>0.05)', xlab = "Missingness Per Sample"); abline(v = samp.missingness.thresh, col = "red")
        text(0.01, 100, paste0("Probe detP>0.05\n", samp.missingness.thresh, " Sample Missingness\nexclude\tretain\n", table(keep.samp)[1], "\t", table(keep.samp)[2]), pos = 4)
        hist(1 - probe.missing / n, breaks = 100, main = 'Probe Missingness (detP>0.05)', xlab = "Missingness Per Probe"); abline(v = probe.missingness.thresh, col = "red")
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

    ####################################################################################################
    ### 3. Map to genome + drop SNP loci, then background-correct (noob) and build the QC-filtered set.
    ### mapToGenome + dropLociWithSnps define the probes to KEEP (mappable, not on a SNP). noob needs
    ### the raw RGChannelSet, so it runs separately on rgSet; the SNP-surviving + detP-passing probe set
    ### is then applied to the noob output so the SNP-locus drop actually reaches the deliverable.
    ####################################################################################################
    grgSet <- mapToGenome(rgSet)        ## Map to genome; also drops probes that do not map.

    ### Fail loud on an array-version mismatch instead of silently mis-mapping/dropping probes
    ### under the wrong manifest (minfi resolves the manifest from the IDAT header, not a flag).
    resolved.array <- annotation(grgSet)["array"]
    cat("Array annotation resolved to:", resolved.array, annotation(grgSet)["annotation"], "\n")
    if (ARRAY_VERSION == "v2" && !grepl("EPICv2", resolved.array))
        stop("ARRAY_VERSION=v2 but minfi resolved the array annotation to '", resolved.array,
             "'; check the raw IDATs / installed manifest packages.")
    if (ARRAY_VERSION == "v1" && grepl("EPICv2", resolved.array))
        stop("ARRAY_VERSION=v1 but minfi resolved the array annotation to '", resolved.array,
             "'; the IDATs look like EPIC v2, set METHYL_ARRAY_VERSION=v2.")

    grgSet.ns <- dropLociWithSnps(grgSet)  ### DROP SNP-affected loci (CpG-site + single-base-extension SNPs)
    cat("probes: mapped", nrow(grgSet), "-> after dropLociWithSnps", nrow(grgSet.ns),
        "(dropped", nrow(grgSet) - nrow(grgSet.ns), "SNP-affected)\n")
    stopifnot(all(colnames(grgSet.ns) == pd$Sample_Group))   ### sample order matches metadata (fail loud)
    rm(grgSet); gc()   ### mapToGenome result consumed

    ### The probes to keep in the deliverable: mappable + not-on-SNP (rownames(grgSet.ns)) AND
    ### detP-passing (not in rm.probe). Applied to the noob output below; this is the step that
    ### carries the SNP drop into MSetNoob.flt (and thus into dasen_betas.RDat).
    keep.probes <- setdiff(rownames(grgSet.ns), names(rm.probe))
    rm(grgSet.ns); gc()   ### genome-mapped set consumed (only the probe id list is carried forward)

    ### Filter sample metadata to the retained samples:
    pd  <- pd[keep_idx(nrow(pd), rm.samp), ]
    tis <- tis[keep_idx(nrow(tis), rm.samp), ]
    cat("Completed filtering probes and samples on detP + SNP loci\n", date(), "\n\n")

    ### noob background correction (needs the raw RGChannelSet; failed samples removed first):
    rgSet.rmsamp <- rgSet[, keep_idx(ncol(rgSet), rm.samp)] ### RGChannelSet with failed samples removed
    rm(rgSet); gc()   ### raw rgSet no longer needed; free ~15GB immediately before preprocessNoob
    MSetNoob <- preprocessNoob_chunked(rgSet.rmsamp) ### chunked by sample (config.R) to cap peak memory
    rm(rgSet.rmsamp); gc()   ### noob input consumed
    if (SAVE_INTERMEDIATES) save(MSetNoob, file = F_NOOB)
    cat("Completed noob background correction\n", date(), "\n\n")

    ### Apply the detP + SNP-locus filter to the noob output. keep.probes came from the SNP-dropped,
    ### genome-mapped set above; intersecting with rownames(MSetNoob) via a logical mask guards
    ### against any manifest id mismatch and preserves the MethylSet's row order.
    keep.mask <- rownames(MSetNoob) %in% keep.probes
    stopifnot(sum(keep.mask) > 0)
    MSetNoob.flt <- MSetNoob[keep.mask, ]
    cat("probes: noob", nrow(MSetNoob), "-> after detP + SNP filter", nrow(MSetNoob.flt),
        "(dropped", nrow(MSetNoob) - nrow(MSetNoob.flt), ")\n")
    stopifnot(!any(rownames(MSetNoob.flt) %in% names(rm.probe)))   ### no detP-failed probe survives
    rm(MSetNoob); gc()   ### MSetNoob.flt carries forward; free MSetNoob (~25GB) before dasen

    ### Self-contained checkpoint: the QC-filtered MethylSet + the retained-sample metadata, so a
    ### resume can jump straight to dasen without re-reading IDATs or recomputing detP (see the top).
    if (SAVE_INTERMEDIATES) {
        ckpt <- list(mset = MSetNoob.flt, pd = pd, tis = tis)
        save(ckpt, file = F_NOOBFLT)
        rm(ckpt)
    }
}

########################################################################################################
### Normalize with dasen (tiers 1/2 only; tier 0 already loaded dasen.melon from F_DASEN). Default:
### dasen_stream() (streaming/low-memory reimplementation in config.R), output identical to
### wateRmelon::dasen but with a ~40-60GB peak instead of the >100GB stock needs at cohort scale.
### F_DASEN is written BEFORE the memory-heavy QC/beta extraction below, so a crash there can resume.
########################################################################################################
if (is.null(dasen.melon)) {
    dasen.melon <- if (DASEN_STREAM) dasen_stream(MSetNoob.flt) else dasen(MSetNoob.flt)
    if (SAVE_INTERMEDIATES) save(dasen.melon, file = F_DASEN)
    print(dasen.melon)
    cat("Completed dasen normalization\n", date(), "\n\n")
    ### Sample identity must line up between the normalized set and the metadata (fail loud):
    stopifnot(all(colnames(dasen.melon) == pd$Sample_Group))
}

########################################################################################################
### 4. Betas / M-values (the deliverable) + wateRmelon QC. Free the big MethylSets as soon as each is
### drained so this back half stays well under the memory cap (previously MSetNoob.flt + dasen.melon +
### several beta matrices were all held at once -- the spike that killed the run after dasen).
########################################################################################################
betas.m <- getBeta(dasen.melon)   ### normalized betas: the deliverable AND the dasen side of qual

### wateRmelon `qual` compares the pre-normalization betas against the normalized betas. The
### pre-normalization betas need MSetNoob.flt, which is absent when resuming from F_DASEN (tier 0),
### so those QC plots are skipped there. getBeta() == wateRmelon::betas() for a MethylSet (offset 100).
qc.available <- !is.null(MSetNoob.flt)
if (qc.available) {
    b.noob <- getBeta(MSetNoob.flt)   ### pre-normalization betas (qual input only)
    rm(MSetNoob.flt); gc()            ### drained; free ~2 full matrices before qual
    qu <- qual(b.noob, betas.m)       ### metrics between non-normalized & normalized betas
    rm(b.noob); gc()
    str(qu)
    pdf(file.path(REPORT_DIR, "minfi_QC.wateRmelonQC.pdf"), height = 8, width = 8)
        par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(1.75, .75, 0))
        boxplot(qu[, 1] ~ pd$DNA_Source)
        boxplot(qu[, 2] ~ pd$DNA_Source)
        plot(qu[, 1], qu[, 2], col = tis$tnumeric, pch = tis$tnumeric); abline(0, 1, col = 1, lty = 3, lwd = 2)
        legend('bottomright', legend = levels(tis$tissue), col = 1:4, pch = 1:4, bty = 'n', title = 'tissue')
    dev.off()

    ### Compare normalization by tissue (excluding cell-line controls):
    rm.qual <- which(!rownames(pd) %in% rownames(qu))
    pd <- pd[keep_idx(nrow(pd), rm.qual), ]
    print(all(rownames(pd) == rownames(qu)))
    ctl.samp  <- which(pd$DNA_Source == "Cell_Line")
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
} else {
    cat("RESUME from F_DASEN: skipping wateRmelon QC plots (pre-normalization betas unavailable)\n")
}

Mvalues <- getM(dasen.melon)
rm(dasen.melon); gc()   ### betas/M extracted; free the normalized MethylSet before canonicalization

### Canonicalize EPIC v2 replicate-probe IDs (strip "_TC21"-style suffixes, collapsing to the
### bare cg######## id clock lists / cell-type references expect) ONCE, here, so every
### downstream reader (stages 3-5) sees canonical ids unconditionally. No-op for v1 data.
betas.m <- canonicalize_v2_probe_ids(betas.m)
Mvalues <- canonicalize_v2_probe_ids(Mvalues)
dasen.values <- list()
dasen.values$b <- betas.m
dasen.values$M <- Mvalues
save(dasen.values, file = F_DASENB)
cat("Completed analyses\n", date(), "\n\n")
