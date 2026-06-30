# dash

Depth-Aware Structural-zero Hellinger (DASH) test for microbiome differential
abundance.

DASH is a two-part, per-taxon test:

1. **Prevalence arm.** A structural-zero score test from a per-taxon
   zero-inflated Tobit model.
2. **Abundance arm.** A weighted least squares test on Hellinger-Riemann
   intrinsic coordinates (HRIC), with an HC3 sandwich standard error and a robust
   affine bias correction.

The two arms are combined per taxon with a Bonferroni min-P rule
(`p_bonf_minp`) and a Cauchy combination (`p_cauchy`).

## Install

```r
# install.packages("devtools")
devtools::install_local("dash")   # or devtools::install("dash")
```

## Use

```r
library(dash)

# counts: samples x taxa integer matrix; group: two-level indicator
res <- dash(counts, group)

# optional covariate adjustment
res <- dash(counts, group, covariates = z)

# multiple-testing correction is applied by the caller
res$q_bonf_minp <- p.adjust(res$p_bonf_minp, "BH")
```

`dash()` is the single entry point. With the defaults `min_positive = 3`,
`d0 = 1`, and `m_abund = 5`, the per-taxon p-values reproduce the simulation
results.

## Development

The `R/` sources carry roxygen2 comments. Regenerate `man/` and `NAMESPACE`
with:

```r
devtools::document("dash")
devtools::check("dash")
```

`MASS` is an optional dependency: when available the affine bias-correction
line is fit with `MASS::lqs` (least trimmed squares); otherwise DASH falls back
to a base-R Theil-Sen fit.
