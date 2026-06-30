# DASH

**DASH** is an R package for two-group microbiome differential abundance analysis in sparse count data.

DASH stands for **Depth-Aware Structural-zero and Hellinger**. The method is designed for settings where zeros are common and biologically ambiguous: an observed zero may reflect limited sequencing depth, or it may indicate structural absence of a taxon in a biological group. DASH combines a depth-aware prevalence component with an abundance component based on Hellinger-Riemann intrinsic coordinates.

## Overview

For each taxon, DASH forms two complementary tests:

1. **Prevalence arm**
   A structural-zero score test derived from a per-taxon zero-inflated Tobit model. Library size enters through a depth-aware censoring threshold, so a zero observed in a deep sample contributes stronger evidence for true absence than a zero observed in a shallow sample.

2. **Abundance arm**
   A weighted least squares test on Hellinger-Riemann intrinsic coordinates (HRIC). The abundance arm uses posterior-presence weights from the prevalence model, includes log library size as a technical nuisance covariate, uses an HC3 sandwich standard error, and applies a robust affine correction for compositional background bias.

The two component p-values are combined per taxon by:

* `p_bonf_minp`: Bonferroni min-P combination, used as the primary DASH omnibus p-value.
* `p_cauchy`: Cauchy combination, provided as a sensitivity analysis.

Multiple testing across taxa is left to the user. In typical differential abundance analysis, apply Benjamini-Hochberg correction to `p_bonf_minp`.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("yiqianomics/DASH")
```

Alternatively, if you have cloned the repository locally:

```r
# install.packages("devtools")
devtools::install_local("DASH")
```

Then load the package:

```r
library(dash)
```

## Quick start

```r
library(dash)

set.seed(1)

n <- 80
p <- 40

counts <- matrix(rpois(n * p, lambda = 5), nrow = n, ncol = p)
colnames(counts) <- paste0("taxon_", seq_len(p))

group <- factor(
  rep(c("Control", "Case"), each = n / 2),
  levels = c("Control", "Case")
)

res <- dash(counts, group)

res$q_bonf_minp <- p.adjust(res$p_bonf_minp, method = "BH")

head(res[order(res$q_bonf_minp), ])
```

## Input format

The main function is:

```r
dash(counts, group, covariates = NULL, d0 = 1, m_abund = 5L, min_positive = 3L)
```

### `counts`

`counts` should be a numeric count matrix or data frame with:

* rows = samples
* columns = taxa
* non-negative count values
* at least two taxa after filtering
* positive total abundance for retained samples

Column names are used as taxon identifiers. If column names are absent, DASH assigns names of the form `tax1`, `tax2`, and so on.

### `group`

`group` must contain exactly two distinct groups.

Numeric group variables are mapped to 0/1 using the larger value as the comparison group. Factor or character group variables are mapped according to their two levels, with the second level treated as the comparison group.

For clarity, it is recommended to specify factor levels explicitly:

```r
group <- factor(group, levels = c("Control", "Case"))
```

### `covariates`

Optional covariates can be supplied as a vector, matrix, or data frame with one row per sample. These covariates are used for adjustment in both the prevalence and abundance arms.

For categorical covariates, use `model.matrix()` to construct a numeric design matrix:

```r
metadata <- data.frame(
  age = rnorm(n),
  batch = factor(rep(1:4, length.out = n))
)

Z <- model.matrix(~ age + batch, data = metadata)[, -1, drop = FALSE]

res <- dash(counts, group, covariates = Z)
```

## Output

`dash()` returns a data frame with one row per retained taxon:

| Column         | Description                                                                   |
| -------------- | ----------------------------------------------------------------------------- |
| `taxon`        | Taxon identifier                                                              |
| `p_prevalence` | P-value from the prevalence arm                                               |
| `p_abundance`  | P-value from the abundance arm; set to 1 when the abundance arm is not formed |
| `p_bonf_minp`  | Primary DASH omnibus p-value using Bonferroni min-P combination               |
| `p_cauchy`     | DASH omnibus p-value using Cauchy combination                                 |

Example:

```r
res <- dash(counts, group)
res$q_dash <- p.adjust(res$p_bonf_minp, method = "BH")

sig <- subset(res, q_dash < 0.05)
sig
```

The returned object also stores three attributes:

```r
attr(res, "group_levels")
attr(res, "kept_taxa")
attr(res, "kept_samples")
```

These record the group coding, retained taxa, and retained samples after preprocessing.

## Default preprocessing

By default, DASH uses the same preprocessing choices as the simulation implementation:

```r
dash(counts, group, d0 = 1, m_abund = 5L, min_positive = 3L)
```

The defaults mean:

* `min_positive = 3`: retain taxa with at least three positive counts across all samples.
* `m_abund = 5`: form the abundance arm only for taxa with at least five positive counts in each group.
* `d0 = 1`: use 1 as the detection constant in the Tobit censoring threshold.

Samples with zero total abundance after taxon filtering are automatically dropped, because HRIC requires a positive row total.

To test all supplied taxa without the default positive-count retention filter:

```r
res <- dash(counts, group, min_positive = 0)
```

## Method summary

DASH targets two complementary signals:

* a change in structural absence or prevalence;
* a change in abundance conditional on non-structural presence.

The prevalence arm uses a depth-aware zero-inflated Tobit model to estimate posterior structural-zero probabilities. The abundance arm uses HRIC coordinates computed from the observed composition, with posterior-presence weights and robust standard errors. The final per-taxon DASH p-value combines the two arms while avoiding an unnecessary two-test penalty when the abundance arm is not formed.