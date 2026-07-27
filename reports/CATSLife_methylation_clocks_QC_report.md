# Epigenetic clocks in the CATSLife methylation cohort: quality control and validity

## 1. Dataset

DNA methylation was assayed on the Illumina Infinium MethylationEPIC v2.0 array in whole blood (buffy coat) from participants of the Colorado Adoption/Twin Study of Lifespan behavioral development and cognitive aging (CATSLife). After quality control the analytic set comprises **1,650 samples from 1,009 individuals in 538 families**, with a mean age of 35.3 years (SD 5.7). The cohort is genetically informative and longitudinal: it contains monozygotic and dizygotic twins, non-twin biological siblings, and adoptive siblings, and 632 individuals were sampled at two waves roughly six years apart.

| Family type | Samples | Individuals | Mean age (SD) |
|---|--:|--:|--:|
| Monozygotic twins | 521 | 321 | 31.8 (3.3) |
| Dizygotic twins | 444 | 279 | 31.7 (3.3) |
| Biological siblings | 308 | 189 | 40.1 (4.5) |
| Adoptive siblings | 248 | 156 | 40.8 (4.6) |
| Unclassified | 129 | 64 | 39.5 (4.9) |

## 2. Processing

Raw IDATs were read with `minfi`. Samples and probes with more than 1% of calls undetected (detection *p* ≥ 0.05) were removed, dropping 44 of 1,689 samples. The remaining data were background-corrected (`noob`), probes overlapping SNPs or mapping ambiguously were dropped (`dropLociWithSnps`), and signals were normalised with `dasen` (`wateRmelon`). EPIC v2 replicate-probe identifiers were collapsed to their base identifiers, retaining the replicate with the lowest missingness.

The fifteen epigenetic clocks were then computed with `dnaMethyAge` on the **normalised, unadjusted** betas. Because the published clocks are fixed-weight predictors trained on normalised input, cell-composition and plate effects are not removed before the clock is applied; they are instead estimated separately (EpiDISH cell proportions; plate batch) and made available as covariates for age-acceleration models.

## 3. Completeness

Of the 1,650 samples, 1,605 (97.3%) carry a value for all fifteen clocks. The 45 with any missing clock have identifiable causes: 38 were dropped at the detection-*p* step and therefore have no betas; three are curated sex-discordant individuals whose clocks were set to missing; and four lack only PCGrimAge, the clock with the widest probe requirement. All excluded and flagged samples are retained in the released table, so alternative inclusion criteria can be applied.

## 4. Clock performance

For each clock, the correlation with chronological age, the technical reliability across intentional duplicate pairs (intraclass correlation, ICC(1,1)), and the within-person stability across the two waves are given below.

| Clock | Reports | Age *r* | Reliability (ICC) | Wave-to-wave *r* |
|---|---|--:|--:|--:|
| PCGrimAge | age / mortality | 0.83 | **0.98** | 0.90 |
| ZhangQ | age | **0.95** | 0.74 | 0.95 |
| DunedinPACE | pace of aging | 0.17 | **0.92** | 0.76 |
| DNAmTL (LuA2019) | telomere length (kb) | −0.63 | **0.90** | 0.84 |
| Horvath2 (Skin & Blood) | age | 0.92 | 0.55 | 0.90 |
| Hannum | age | 0.84 | 0.65 | 0.87 |
| Shireby (cortical) | age | 0.74 | 0.69 | 0.79 |
| PanMammalian 2 | age | 0.70 | 0.70 | 0.70 |
| PedBE | pediatric age | 0.67 | 0.73 | 0.65 |
| PanMammalian 3 | age | 0.54 | 0.70 | 0.59 |
| epiTOC2 | mitotic divisions | 0.17 | 0.56 | 0.61 |
| epiTOC | mitotic divisions | 0.14 | 0.45 | 0.52 |
| ZhangY | mortality score | 0.27 | 0.63 | 0.61 |
| Horvath (2013) | age | 0.80 | **0.39** | 0.75 |
| PhenoAge | phenotypic age | 0.74 | **0.34** | 0.76 |

The trained age clocks correlate between 0.54 and 0.95 with chronological age. The low values for DunedinPACE, epiTOC/epiTOC2, and ZhangY are by design, since these estimate the pace of aging, mitotic history, and mortality risk rather than age. Even for the age clocks the correlations are attenuated, because the cohort's narrow age range (SD 5.7 years) restricts the variance available to predict.

Accuracy is only one axis, and technical reliability varies widely without following it. The principal-component and second-generation clocks are the most reliable (PCGrimAge ICC 0.98, DunedinPACE 0.92, DNAmTL 0.90), whereas the first-generation Horvath (0.39) and PhenoAge (0.34) clocks are considerably noisier, in line with their known reliability limitations. These estimates come from twelve duplicate pairs and so indicate rank order rather than exact values.

DNAmTL differs in kind from the other measures. Labelled `LuA2019` in `dnaMethyAge`, it estimates telomere length in kilobases (range 6.6–8.0, median 7.4) rather than age, and its negative correlation with chronological age is the expected biology rather than a failure.

## 5. Properties relevant to analysis

**Relatedness.** The cohort is predominantly twins and siblings, so its samples are not independent. This dependence does not distort the correlations above: recomputing them on one randomly chosen sample per person, and then per family, changes every clock by at most 0.04. Analyses must nonetheless account for family clustering, through a family random effect or cluster-robust standard errors, while the same family and zygosity structure supports twin and family designs directly.

**Repeated measures.** The two waves, a mean of 6.1 years apart, allow within-person change models. Across-wave stability is high for the reliable clocks (ZhangQ 0.95, Horvath2 and PCGrimAge 0.90), and for the trained age clocks methylation age advanced close to one year per chronological year (0.8–1.2), as expected.

**Batch structure.** Principal-components analysis shows that the leading axis of variation is sex, which accounts for 99% of the first component, as expected when sex chromosomes are retained. Once they are removed, the dominant axis becomes sequencing plate (86% of the first autosomal component), a genuine batch effect. Because the clocks use unadjusted betas, analyses of age acceleration should include cell composition and plate as covariates; intrinsic and extrinsic variants, adjusting for plate and cell proportions or for plate alone, are provided for this.

**Sex.** Sex was not verified independently from chromosomal intensity, so the three sex-discordant exclusions rest on curated flags rather than a genotype- or intensity-based determination, and a small number of further discordances cannot be excluded.

## 6. Recommendations

- **Prefer reliable measures.** Where a construct allows a choice, favour PCGrimAge, DunedinPACE, ZhangQ, and DNAmTL; treat the Horvath (2013) and PhenoAge age estimates as noisy in this cohort, and interpret single-clock findings from them cautiously.
- **Model relatedness.** Use family-clustered or mixed models, and do not treat samples as independent.
- **Adjust age acceleration** for cell composition and plate, choosing the intrinsic or extrinsic variant to match the question.
- **Honour the exclusion flags.** Curated sex-discordant individuals and detection-*p* failures are retained but flagged; apply the provided `clock_excluded` indicator unless there is reason to override it.
- **Use DNAmTL as telomere length,** not as an age or age-acceleration measure.
