### stage6/sex_helpers.R: shared helpers for the stage-6 sex checks (pca_sex_batch.R,
### sex_intensity.R, sex_verdicts.R). Sourced after config.R.

## Two-component 1-D Gaussian mixture fitted by EM, without labels. Returns the component
## means, sds and weights (component 1 = lower mean) and each value's posterior probability
## of belonging to component 2. Posteriors are computed in log space, so values far from both
## components still get a finite, well-defined probability.
fit_two_gaussians <- function(x, max_iter = 1000L, tol = 1e-10) {
    ok <- is.finite(x)
    if (sum(ok) < 4) stop("fit_two_gaussians: fewer than 4 finite values")
    xs <- x[ok]
    km <- stats::kmeans(xs, centers = range(xs))
    mu <- as.numeric(km$centers); g <- km$cluster
    floor_sd <- 1e-6 * diff(range(xs))
    sd <- pmax(c(stats::sd(xs[g == 1]), stats::sd(xs[g == 2])), floor_sd, na.rm = TRUE)
    w  <- as.numeric(table(factor(g, levels = 1:2))) / length(xs)
    post2 <- function(v, mu, sd, w) {
        l1 <- log(w[1]) + stats::dnorm(v, mu[1], sd[1], log = TRUE)
        l2 <- log(w[2]) + stats::dnorm(v, mu[2], sd[2], log = TRUE)
        1 / (1 + exp(l1 - l2))
    }
    ll_old <- -Inf
    for (it in seq_len(max_iter)) {
        r2 <- post2(xs, mu, sd, w); r1 <- 1 - r2
        w  <- c(mean(r1), mean(r2))
        mu <- c(sum(r1 * xs) / sum(r1), sum(r2 * xs) / sum(r2))
        sd <- pmax(sqrt(c(sum(r1 * (xs - mu[1])^2) / sum(r1), sum(r2 * (xs - mu[2])^2) / sum(r2))), floor_sd)
        ll <- sum(log(w[1] * stats::dnorm(xs, mu[1], sd[1]) + w[2] * stats::dnorm(xs, mu[2], sd[2])))
        if (is.finite(ll) && abs(ll - ll_old) < tol * (1 + abs(ll))) break
        ll_old <- ll
    }
    if (mu[1] > mu[2]) { mu <- rev(mu); sd <- rev(sd); w <- rev(w) }
    p <- rep(NA_real_, length(x)); p[ok] <- post2(xs, mu, sd, w)
    list(mean = mu, sd = sd, weight = w, p_upper = p, iterations = it)
}

## Robust z of each value against a reference group (median / MAD). The MAD is floored so a
## near-constant group does not turn ordinary values into extreme z.
robust_z <- function(x, ref) {
    ref <- ref[is.finite(ref)]
    s <- max(stats::mad(ref), 1e-6 * max(1, abs(stats::median(ref))))
    (x - stats::median(ref)) / s
}

## Classify each value of a sex-informative measure against its female and male reference
## groups: "F" or "M" when within `z_max` robust SDs of that group's centre (the nearer one if
## both), else "atypical" (outside both sexes' typical range).
sex_typical <- function(x, ref_f, ref_m, z_max) {
    zf <- abs(robust_z(x, ref_f)); zm <- abs(robust_z(x, ref_m))
    out <- ifelse(zf <= z_max & (zf <= zm | zm > z_max), "F",
           ifelse(zm <= z_max, "M", "atypical"))
    out[!is.finite(x)] <- NA_character_
    out
}

## Probe annotation of the array version in use (rownames = the array's probe ids). The
## package is attached with library(): minfi::getAnnotation() -> updateObject() looks it up as
## "package:<pkg>" on the search list.
array_annotation <- function() {
    pkg <- if (ARRAY_VERSION == "v2") "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
           else "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    if (!requireNamespace(pkg, quietly = TRUE)) stop("annotation package not installed: ", pkg)
    suppressMessages(library(pkg, character.only = TRUE))
    minfi::getAnnotation(get(pkg))
}

## Classify each SAMPLE_SWAPS_FILE action by what a sex check can see. A relabel whose sheet and
## corrected persons differ in admin sex is testable by sex; a same-sex relabel is invisible to
## any sex check.
relabel_type <- function(identity_action, sex_sheet, sex_admin) {
    ifelse(!identity_action %in% "relabel",
           ifelse(identity_action %in% "flag", "flag only", NA_character_),
    ifelse(is.na(sex_sheet) | is.na(sex_admin), "relabel, sex unknown",
    ifelse(sex_sheet != sex_admin, "cross-sex relabel", "same-sex relabel (not testable by sex)")))
}

## Read sex_qc.csv and require the label-provenance columns the pre/post check uses.
read_sex_qc <- function(path = file.path(REPORT_DIR, "sex_qc.csv")) {
    if (!file.exists(path))
        stop("missing ", path, "; run scripts/build/build_phenotype_file.R first.")
    sqc  <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
    need <- c("Sample", "Subject_ID", "Subject_ID_sheet", "identity_action", "IndividualID",
              "Sex_admin", "Sex_admin_sheet")
    miss <- setdiff(need, names(sqc))
    if (length(miss))
        stop(path, " lacks column(s) ", paste(miss, collapse = ", "),
             "; rebuild it with scripts/build/build_phenotype_file.R.")
    sqc[sqc == ""] <- NA_character_
    sqc
}
