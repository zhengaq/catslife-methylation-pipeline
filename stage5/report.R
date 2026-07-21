### stage5/report.R — assemble the core stage-5 tables (output/tables/) and figures
### (output/reports/): descriptives, clock-vs-age validity, per-clock violin plots,
### and cross-tissue rank corrplots.
source("config.R"); source("stage5/helpers.R")
suppressMessages({ library(dplyr); library(ggplot2); library(corrplot) })

TABLES_DIR <- file.path(ANALYSIS_DIR, "tables"); dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
have <- function(f) file.exists(file.path(ANALYSIS_DIR, f))

m  <- read.csv(file.path(ANALYSIS_DIR, "mAge_clocks.csv"))
ph <- read.csv(PHENOTYPE_FILE)
if ("Zygosity" %in% names(ph)) m <- merge(m, unique(ph[, c("Sample", "Zygosity")]), by = "Sample", all.x = TRUE)
clock_cols <- intersect(LME_CLOCKS, names(m))

## ---- Descriptives: N / age by tissue x zygosity ---------------------------
grp <- if ("Zygosity" %in% names(m)) c("DNA_Source", "Zygosity") else "DNA_Source"
desc <- m %>% group_by(across(all_of(grp))) %>%
  summarise(n_samples = n(), n_individuals = n_distinct(aid),
            age_mean = round(mean(Age, na.rm = TRUE), 1), age_sd = round(sd(Age, na.rm = TRUE), 1),
            .groups = "drop")
write.csv(desc, file.path(TABLES_DIR, "descriptives.csv"), row.names = FALSE)

## ---- Per-clock descriptives: cohort (one/person) + full-sample + per-wave ---
if ("Wave" %in% names(ph)) m <- merge(m, unique(ph[, c("Sample", "Wave")]), by = "Sample", all.x = TRUE)
desc1 <- function(v) data.frame(n = sum(!is.na(v)), mean = round(mean(v, na.rm = TRUE), 3),
                                sd = round(sd(v, na.rm = TRUE), 3))
strat_rows <- function(df, stratum) do.call(rbind, lapply(clock_cols, function(cl)
  cbind(clock = cl, stratum = stratum, desc1(df[[cl]]))))
## "cohort" = one value per individual (repeated samples collapsed to a per-person mean), so
## multiply-sampled people are not over-weighted; "all" keeps every sample; "wave*" splits by wave.
cohort <- do.call(rbind, lapply(clock_cols, function(cl)
  cbind(clock = cl, stratum = "cohort", desc1(tapply(m[[cl]], m$aid, mean, na.rm = TRUE)))))
clock_desc <- rbind(cohort, strat_rows(m, "all"))
if ("Wave" %in% names(m))
  for (w in sort(unique(stats::na.omit(m$Wave))))
    clock_desc <- rbind(clock_desc, strat_rows(m[m$Wave %in% w, , drop = FALSE], paste0("wave", w)))
clock_desc <- clock_desc[order(clock_desc$clock, clock_desc$stratum), ]
write.csv(clock_desc, file.path(TABLES_DIR, "clock_descriptives.csv"), row.names = FALSE)
cat("report: wrote clock_descriptives.csv -", length(unique(clock_desc$clock)), "clocks x",
    length(unique(clock_desc$stratum)), "strata (cohort/all/wave)\n")

## ---- Clock validity: correlation of each clock with chronological age ------
## cor(..., use="complete.obs") ERRORS (not just warns) on zero complete pairs —
## e.g. a clock compute_clocks() had to NA-fill entirely. NA that one clock, not
## the whole pipeline.
validity <- do.call(rbind, lapply(clock_cols, function(cl) data.frame(
  clock = cl, r_age = tryCatch(suppressWarnings(cor(m[[cl]], m$Age, use = "complete.obs")),
                               error = function(e) NA_real_))))
write.csv(validity, file.path(TABLES_DIR, "clock_age_validity.csv"), row.names = FALSE)

## ---- Violin plots (one random member per family) --------------------------
set.seed(123)
fam1 <- m %>% group_by(pfamid) %>% slice_sample(n = 1) %>% select(pfamid, aid) %>% ungroup()
mf <- m %>% semi_join(fam1, by = c("pfamid", "aid"))
mf$DNA_Source <- factor(mf$DNA_Source, levels = DNA_SOURCES)
for (cl in clock_cols) {
  p <- ggplot(mf, aes(x = DNA_Source, y = .data[[cl]], fill = DNA_Source)) +
    geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_classic() + theme(legend.position = "none") + labs(title = cl, x = "DNA source", y = cl)
  ggsave(file.path(REPORT_DIR, paste0("violin_", cl, ".png")), p, width = 5, height = 4)
}

## ---- Rank corrplots (from population's rank_corr.rds) ----------------------
## Skip clocks with no finite cross-tissue overlap — corrplot() errors on an all-NA matrix.
if (have("rank_corr.rds")) {
  R <- readRDS(file.path(ANALYSIS_DIR, "rank_corr.rds"))
  for (cl in names(R)) {
    if (!any(is.finite(R[[cl]]$mean))) { cat("rank corrplot: skip", cl, "(no finite cross-tissue overlap)\n"); next }
    png(file.path(REPORT_DIR, paste0("rank_", cl, ".png")), width = 600, height = 600)
    corrplot(R[[cl]]$mean, method = "circle", type = "upper", tl.col = "black",
             addCoef.col = "black", number.cex = 0.8, title = paste("Rank corr -", cl), mar = c(0, 0, 1, 0))
    dev.off()
  }
}

cat("report: tables ->", TABLES_DIR, "; figures ->", REPORT_DIR, "\n")
