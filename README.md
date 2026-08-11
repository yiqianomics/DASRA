# DORAM

DORAM performs taxon-level association testing for two complementary features
of microbiome count data:

- **Occupancy:** whether latent structural absence differs between two groups.
- **Abundance:** whether zero-inclusive read-relative abundance differs between
  two groups.

It uses each raw count together with the sample's original sequencing depth,
so an observed zero can be distinguished probabilistically from a
finite-depth sampling zero. DORAM also reports a covariance-aware joint test of
the two endpoints.

## Installation

Install DORAM from GitHub with `remotes`:

```r
install.packages("remotes") # run once if remotes is not installed
remotes::install_github("yiqianomics/DORAM")
library(DORAM)
```

Alternatively, users who already have `devtools` can install DORAM with:

```r
devtools::install_github("yiqianomics/DORAM")
```

## Quick start

`counts` is a samples-by-taxa raw count table. `metadata` is a data frame with
the same sample IDs in its row names.

```r
fit <- DORAM(
  counts = count_table,
  metadata = sample_data,
  group = "condition",
  reference = "control",
  library_size = "original_read_depth",
  covariates = c("age", "sex"),
  cores = 4L
)

fit
fit$results
```

The main function is:

```r
DORAM(
  counts,
  metadata,
  group,
  library_size,
  covariates = NULL,
  reference = NULL,
  taxa = NULL,
  cores = 1L,
  verbose = FALSE
)
```

## Input requirements

- `counts` must contain finite nonnegative integer counts, with samples in rows
  and taxa in columns.
- `counts` and `metadata` must have unique row names containing exactly the
  same sample IDs. Metadata is aligned to the count table by ID.
- `group` names a metadata column with exactly two observed groups. Specify
  `reference` for character or factor groups.
- `library_size` names the metadata column containing each sample's original,
  unfiltered sequencing depth. Do not recompute it after selecting taxa.
- `covariates` optionally names metadata columns for adjustment. Numeric and
  categorical covariates are supported; do not include the tested group again.
- `taxa` optionally selects taxa to analyze. By default, every count-table
  column is analyzed, and each analyzed taxon produces one row in `results`.

DORAM does not normalize counts, add pseudocounts, discard zeros, or
automatically filter taxa.

## Results

`fit$results` is the main output. It contains one row per analyzed taxon:

| Column | Meaning |
|---|---|
| `taxon` | Taxon name |
| `p_occupancy`, `q_occupancy` | Raw and BH-adjusted occupancy p-values |
| `p_abundance`, `q_abundance` | Raw and BH-adjusted abundance p-values |
| `p_joint`, `q_joint` | Raw and BH-adjusted joint p-values |

BH adjustment is performed separately for the occupancy, abundance, and joint
families. The set selected by `taxa` defines each multiplicity-adjustment
family.

Additional components are available when needed:

```r
fit$details       # estimates, statistics, degrees of freedom, and status
fit$descriptives  # zeros, positive counts, and mean Y/N by group
fit$diagnostics   # numerical diagnostics for unavailable tests
```

If a taxon does not contain enough information for reliable latent-model
inference, its unavailable p- and q-values are reported as `NA`; the taxon is
retained and the reason is recorded in `fit$details`.

## Interpretation

The occupancy test concerns latent structural absence, not observed
zero-versus-nonzero prevalence. Because occupancy equals one minus structural
absence, the two-sided null is equivalently equality of occupancy between
groups.

The abundance test uses `Y/N` from every sample, including zero counts. It is
therefore different from a positive-count-only abundance test. The estimand is
read-relative abundance, not absolute abundance.

The joint test evaluates whether either endpoint differs between groups while
accounting for their estimated covariance; it is not a combination of two
marginal p-values.

DORAM assumes independent samples and a scientifically reasonable
structural-absence/count-mixture model. As with other read-relative methods, it
does not by itself remove compositional closure.
