# DORAM

DORAM performs depth-aware, taxonwise inference for two ecological endpoints
in microbiome count data: latent occupancy and zero-inclusive read-relative
abundance. It models each raw count together with the original sequencing
depth, allowing an observed zero to be either structural absence or a
finite-depth sampling zero.

For every taxon, DORAM reports an occupancy test, an ecological abundance
test, and a covariance-aware joint test. It does not normalize counts, add
pseudocounts, discard zeros, or use resampling calibration.

## Installation

```r
remotes::install_github("yiqianomics/DORAM")
library(DORAM)
```

## Quick start

`counts` is a samples-by-taxa table of raw counts. `metadata` is a data frame
whose row names identify the same samples. `group`, `library_size`, and
`covariates` name columns in `metadata`.

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

The public interface is intentionally compact:

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

## Inputs

### Counts and sample IDs

- `counts` must have samples in rows, taxa in columns, and finite
  nonnegative-integer entries.
- Both `counts` and `metadata` must have unique, nonempty sample row names and
  exactly the same sample-ID set.
- Metadata rows may be in a different order. DORAM matches them to `counts` by
  ID; it never silently drops samples or binds them by position.

### Original sequencing depth

`library_size` must name a metadata column containing the original,
unfiltered, positive integer sequencing depth for every sample. Do not replace
it with the row sum of a taxon-filtered table. DORAM checks numerical
compatibility between counts and depth, but it cannot verify the column's
provenance.

### Group and reference

`group` must name a metadata column with exactly two observed levels. Supply
`reference` for character or factor groups. It may be omitted only when the
group is logical (`FALSE` is the reference) or numeric 0/1 (0 is the
reference). The fitted comparison is recorded in `fit$contrast`.

### Covariates

`covariates` is `NULL` or a vector of metadata column names. Numeric variables
are used directly. Logical, factor, and character variables are expanded into
treatment-coded indicator columns. Factor level order determines its baseline;
character levels use radix-sorted order. This coding is fixed by DORAM and does
not depend on the user's global contrast options.

Do not repeat the tested group in `covariates`. DORAM constructs one numeric
covariate design and uses it in the structural, conditional-present, and
ecological components.

### Taxa

By default, every count-table column is tested. `taxa` can select unique taxon
names or column indices. The selected set defines the family used for BH
adjustment, so any analysis filter should be specified before examining
group-association results. Taxon selection never changes the original library
sizes.

## Results

`fit$results` contains one row per taxon:

| Column | Meaning |
|---|---|
| `p_occupancy`, `q_occupancy` | Raw and BH-adjusted occupancy p-values |
| `ecological_difference` | Adjusted comparison-minus-reference difference in all-sample `Y/N` |
| `p_ecological`, `q_ecological` | Raw and BH-adjusted ecological p-values |
| `p_joint`, `q_joint` | Raw and BH-adjusted joint p-values |
| `status_*` | `ok` or an explicit reason that an endpoint was unavailable |

The detailed three-row-per-taxon table is `fit$tests`; it includes the test
statistic, degrees of freedom, log p-value, and endpoint availability.
`as.data.frame(fit)` returns the concise `fit$results` table.

BH adjustment is performed separately for occupancy, ecological, and joint
families. Every selected taxon remains in the corresponding family denominator.
An unavailable endpoint retains `NA` p- and q-values and is never counted as a
discovery.

DORAM fails closed: it does not drop a failed taxon, substitute p = 1, retry
with another calibration, or collapse the joint test to one degree of freedom.
Ecological inference can remain available when occupancy is unavailable, but
the joint test requires both components.

Additional information is available without cluttering the printed result:

```r
fit$descriptives       # zeros, positives, and mean Y/N by group
fit$diagnostics        # numerical fit and endpoint diagnostics
fit$posterior          # restricted-null structural probabilities
```

## Endpoint interpretation

- **Occupancy** is a one-degree-of-freedom test of equality of the latent
  structural-absence probability. Because occupancy is one minus structural
  absence, this is equivalently an occupancy-equality test. It is not a test of
  observed zero/nonzero prevalence.
- **Ecological abundance** is a one-degree-of-freedom adjusted group contrast
  in `Y/N` using every sample, including all zeros. It is read-relative, not
  positive-only or absolute abundance.
- **Joint** is a covariance-aware two-degree-of-freedom test that neither
  endpoint changes. It uses the cross-endpoint covariance and is not a
  combination of marginal p-values.

Occupancy and ecological abundance can capture overlapping aspects of the same
ecological change. They should not be described as independent biological
mechanisms.

## Restricted-null structural probabilities

The fitted object contains three samples-by-taxa matrices:

```r
fit$posterior$rho_null
fit$posterior$gamma_null
fit$posterior$tau_null
```

`rho_null` is the fitted structural-absence mixing probability under the
occupancy null. `gamma_null` is the corresponding posterior structural-absence
probability after observing count and depth; it is exactly zero for a positive
count. `tau_null` is `1 - gamma_null`.

These are plug-in quantities from the restricted model with the occupancy
group effect fixed at zero. They do not include parameter uncertainty, are not
posteriors from an unrestricted alternative model, and must not be interpreted
as estimated group effects. DORAM therefore does not report an occupancy odds
ratio or confidence interval.

## Scope

- The latent occupancy interpretation depends on the structural-absence and
  conditional-present count mixture being scientifically reasonable.
- Exact sequencing depth addresses finite-depth sampling but does not by
  itself remove latent depth dependence. Any depth adjustment should be
  scientifically prespecified and changes the conditional estimand.
- The ecological endpoint is read-relative and does not remove compositional
  closure.
- Empirical-sandwich calibration robustifies score variance; it cannot repair
  a misspecified latent count model.
- Inference is first-order asymptotic and assumes independent samples.
  Clustered, repeated, longitudinal, or matched observations require a method
  that represents their dependence structure.
