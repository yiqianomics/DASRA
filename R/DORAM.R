#' Test latent occupancy and zero-inclusive abundance
#'
#' @description
#' `DORAM()` tests two endpoints for every selected taxon. The
#' occupancy endpoint concerns latent structural absence rather than observed
#' zero/nonzero prevalence. The abundance endpoint compares the read-relative
#' abundance `Y/N` across all samples, including zeros. A covariance-aware
#' two-degree-of-freedom joint test evaluates whether either endpoint changes.
#'
#' The function uses raw counts and the original sequencing depth directly. It
#' does not normalize counts, add pseudocounts, remove zeros, or automatically
#' filter taxa.
#'
#' @param counts A numeric matrix or data frame with samples in rows and taxa
#'   in columns. Entries must be finite nonnegative integers. Unique sample row
#'   names are required.
#' @param metadata A data frame with samples in rows and variables in columns.
#'   Its row names must contain exactly the same sample IDs as `counts`; row
#'   order may differ because DORAM matches samples by ID.
#' @param group Name of the metadata column defining the two groups to compare.
#' @param library_size Name of the metadata column containing each sample's
#'   original, unfiltered sequencing depth. This column is required even when
#'   `counts` contains only selected taxa.
#' @param covariates Optional character vector naming metadata columns for
#'   adjustment. Numeric variables are used directly. Logical, factor, and
#'   character variables are expanded deterministically into treatment-coded
#'   indicator columns. Factor level order defines its baseline; character
#'   levels use radix-sorted order. The tested `group` cannot also be a
#'   covariate.
#' @param reference Reference level for the group comparison. It is required
#'   for character and factor groups. It may be omitted only for logical groups
#'   (where `FALSE` is the reference) or numeric 0/1 groups (where 0 is the
#'   reference).
#' @param taxa Optional unique taxon names or integer column indices to test.
#'   The default tests every taxon. This selection defines each BH family but
#'   never changes the supplied original library sizes.
#' @param cores Positive integer number of taxa to fit in parallel. The default
#'   is 1. Parallel taxon fitting uses process forking on non-Windows systems;
#'   Windows runs serially and warns when a value greater than 1 is supplied.
#' @param verbose If `TRUE`, print taxon-level progress. The default is
#'   `FALSE`.
#'
#' @return An object of class `DORAM` with the following components:
#'
#' - `results`: one row per taxon with occupancy, abundance, and joint p- and
#'   BH-adjusted q-values.
#' - `details`: long-form estimates, test statistics, degrees of freedom,
#'   p-values, q-values, and endpoint availability.
#' - `descriptives`: taxon-level zero counts, positive counts, and group means.
#' - `diagnostics`: fit, endpoint, and parameter-boundary diagnostics.
#' - `posterior`: restricted-null `rho_null`, `gamma_null`, and `tau_null`
#'   matrices.
#' - `contrast`: reference and comparison group labels.
#' - `sample_id` and `taxa`: analyzed sample and taxon identifiers.
#' - `settings`: input-column, covariate-coding, multiplicity, and parallel
#'   settings.
#' - `call`: the matched function call.
#'
#' @details
#' The occupancy statistic is a nuisance-projected score test fitted under the
#' restriction of no group effect on structural absence. Present-component
#' group location and scale effects remain unrestricted nuisance parameters.
#' Therefore, `posterior$gamma_null`, `posterior$rho_null`, and
#' `posterior$tau_null` are plug-in quantities under the restricted occupancy
#' null; they are not posterior predictions from an unrestricted group-effect
#' model. DORAM does not report an occupancy odds ratio or confidence interval.
#'
#' The abundance estimate in `details` is the adjusted comparison-minus-reference
#' coefficient for all-sample `Y/N`. It is read-relative, not absolute
#' abundance, and no compositional-closure correction is performed. The joint
#' statistic uses the estimated cross-endpoint covariance; it is not a
#' combination of two marginal p-values.
#'
#' Multiplicity adjustment is performed separately for occupancy, abundance,
#' and joint tests with the Benjamini--Hochberg procedure. Every selected taxon
#' remains in its family denominator. An unavailable test retains `NA` for its
#' p- and q-values and is not a discovery.
#'
#' @examples
#' set.seed(20260810)
#' n <- 180L
#' sample_id <- sprintf("sample_%03d", seq_len(n))
#' group <- rbinom(n, 1L, 0.5)
#' age_z <- rnorm(n)
#' reads <- pmax(1500L, round(exp(rnorm(n, log(8000), 0.35))))
#'
#' # The first taxon has a structural-absence group effect, the second has a
#' # conditional-present abundance effect, and the third is a null taxon.
#' alpha0 <- c(-1.6, -1.2, -1.4)
#' alpha_z <- c(0.25, -0.15, 0.15)
#' delta <- c(2, 0, 0)
#' beta0 <- c(-6.2, -6.5, -6.3)
#' beta_z <- c(0.15, 0.20, -0.10)
#' zeta <- c(0, 0.65, 0)
#' sigma <- c(0.65, 0.70, 0.75)
#'
#' rho <- sapply(seq_len(3), function(j) {
#'   plogis(alpha0[j] + alpha_z[j] * age_z + delta[j] * group)
#' })
#' structural <- matrix(rbinom(n * 3, 1L, as.vector(rho)), n, 3)
#' latent_abundance <- sapply(seq_len(3), function(j) {
#'   rnorm(n, beta0[j] + beta_z[j] * age_z + zeta[j] * group, sigma[j])
#' })
#' taxon_probability <- plogis(latent_abundance) * (1 - structural)
#' stopifnot(max(rowSums(taxon_probability)) < 1)
#'
#' # Draw the three focal taxa jointly and leave the remaining probability
#' # mass as an implicit background community.
#' count_table <- t(vapply(seq_len(n), function(i) {
#'   probabilities <- c(
#'     taxon_probability[i, ],
#'     1 - sum(taxon_probability[i, ])
#'   )
#'   as.integer(rmultinom(1, reads[i], probabilities)[1:3, 1])
#' }, integer(3)))
#' colnames(count_table) <- c(
#'   "joint_signal", "abundance_signal", "null_taxon"
#' )
#' rownames(count_table) <- sample_id
#' sample_data <- data.frame(
#'   diagnosis = factor(
#'     ifelse(group == 0, "control", "case"),
#'     levels = c("control", "case")
#'   ),
#'   reads = reads,
#'   age_z = age_z,
#'   row.names = sample_id
#' )
#'
#' \dontrun{
#' fit <- DORAM(
#'   counts = count_table,
#'   metadata = sample_data,
#'   group = "diagnosis",
#'   reference = "control",
#'   library_size = "reads",
#'   covariates = "age_z",
#'   cores = 3L
#' )
#' fit
#' fit$results
#' fit$details[, c(
#'   "taxon", "endpoint", "available", "p_value", "q_value", "status"
#' )]
#' }
#'
#' @importFrom stats constrOptim dbinom dnorm integrate lm.fit median pchisq
#' @importFrom stats plogis pnorm qlogis sd setNames uniroot
#' @importFrom Rcpp evalCpp
#' @useDynLib DORAM, .registration = TRUE
#' @importFrom utils tail
#' @export
DORAM <- function(
    counts,
    metadata,
    group,
    library_size,
    covariates = NULL,
    reference = NULL,
    taxa = NULL,
    cores = 1L,
    verbose = FALSE) {
  call <- match.call()
  cores <- .doram_cores(cores)
  verbose <- .doram_logical_scalar(verbose, "verbose")

  counts <- .doram_counts_matrix(counts)
  sample_id <- rownames(counts)
  aligned <- .doram_metadata(metadata, sample_id)
  metadata <- aligned$data

  group_info <- .doram_group(metadata, group, reference)
  depth_info <- .doram_library_size(metadata, library_size, counts)
  covariate_info <- .doram_covariates(
    metadata,
    covariates = covariates,
    group = group_info$column
  )
  Z <- covariate_info$matrix
  X <- matrix(
    1,
    nrow = nrow(counts),
    ncol = 1L,
    dimnames = list(sample_id, "Intercept")
  )
  if (!is.null(Z)) X <- cbind(X, Z)

  tolerance <- CN_D3_CONTRACT$precheck$rank_tolerance
  if (qr(X, tol = tolerance)$rank < ncol(X)) {
    .doram_abort("the intercept-plus-covariate design matrix is rank deficient")
  }
  group_residual <- stats::lm.fit(X, group_info$encoded)$residuals
  if (!all(is.finite(group_residual)) || sum(group_residual^2) <= tolerance) {
    .doram_abort("the group is duplicated or perfectly explained by covariates")
  }

  selected_taxa <- .doram_taxa(taxa, colnames(counts))
  selected_counts <- counts[, selected_taxa, drop = FALSE]
  original_depth <- depth_info$values

  fit_taxon <- function(taxon) {
    tryCatch(
      .doram_fit_one(
        taxon = taxon,
        y = selected_counts[, taxon],
        library_size = original_depth,
        group = group_info$encoded,
        covariates = Z,
        subject_id = sample_id,
        X = X
      ),
      error = function(e) {
        .doram_exception_fit(taxon, nrow(selected_counts), conditionMessage(e))
      }
    )
  }
  fitted_taxa <- .doram_map_taxa(
    selected_taxa,
    fit_taxon,
    cores = cores,
    verbose = verbose
  )
  names(fitted_taxa) <- selected_taxa

  details <- do.call(rbind, lapply(fitted_taxa, `[[`, "rows"))
  rownames(details) <- NULL
  details <- .doram_adjust_results(details)
  results <- .doram_primary_results(details, selected_taxa)
  descriptives <- .doram_descriptives(
    selected_counts,
    library_size = original_depth,
    group = group_info$encoded,
    taxa = selected_taxa
  )
  rownames(descriptives) <- NULL

  gamma_null <- do.call(cbind, lapply(fitted_taxa, `[[`, "gamma_null"))
  rho_null <- do.call(cbind, lapply(fitted_taxa, `[[`, "rho_null"))
  colnames(gamma_null) <- colnames(rho_null) <- selected_taxa
  rownames(gamma_null) <- rownames(rho_null) <- sample_id
  posterior <- list(
    gamma_null = gamma_null,
    rho_null = rho_null,
    tau_null = 1 - gamma_null
  )

  fit_diagnostics <- do.call(
    rbind,
    lapply(fitted_taxa, `[[`, "fit_diagnostic")
  )
  endpoint_diagnostics <- do.call(
    rbind,
    lapply(fitted_taxa, `[[`, "endpoint_diagnostics")
  )
  boundary_diagnostics <- do.call(
    rbind,
    lapply(fitted_taxa, `[[`, "boundary_diagnostics")
  )
  rownames(fit_diagnostics) <- rownames(endpoint_diagnostics) <- NULL
  rownames(boundary_diagnostics) <- NULL

  structure(
    list(
      call = call,
      results = results,
      details = details,
      descriptives = descriptives,
      diagnostics = list(
        fit = fit_diagnostics,
        endpoint = endpoint_diagnostics,
        boundary = boundary_diagnostics
      ),
      posterior = posterior,
      contrast = c(
        reference = group_info$reference,
        comparison = group_info$comparison
      ),
      sample_id = sample_id,
      taxa = selected_taxa,
      settings = list(
        group_column = group_info$column,
        library_size_column = depth_info$column,
        covariate_columns = covariate_info$columns,
        covariate_design = covariate_info$design_names,
        categorical_levels = covariate_info$levels,
        metadata_reordered = aligned$reordered,
        family_size = length(selected_taxa),
        p_adjust_method = "BH",
        cores = cores
      )
    ),
    class = c("DORAM", "list")
  )
}

#' @noRd
#' @export
print.DORAM <- function(x, ...) {
  if (!inherits(x, "DORAM")) .doram_abort("x must be a DORAM fit")
  cat(
    "DORAM fit: ", x$contrast[["comparison"]], " vs ",
    x$contrast[["reference"]], "\n",
    length(x$sample_id), " samples; ", length(x$taxa), " taxa\n",
    sep = ""
  )
  shown <- utils::head(x$results, 6L)
  display_columns <- c(
    "taxon", "q_occupancy", "q_abundance", "q_joint"
  )
  print(shown[, display_columns, drop = FALSE], row.names = FALSE, digits = 4)
  if (nrow(x$results) > nrow(shown)) {
    cat("... ", nrow(x$results) - nrow(shown),
        " more taxa; see fit$results\n", sep = "")
  }
  p_columns <- c("p_occupancy", "p_abundance", "p_joint")
  if (anyNA(x$results[, p_columns, drop = FALSE])) {
    cat("Unavailable tests are NA; see fit$details for more information.\n")
  }
  invisible(x)
}

#' @noRd
#' @export
as.data.frame.DORAM <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}
