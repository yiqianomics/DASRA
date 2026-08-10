#' Test latent occupancy and zero-inclusive ecological abundance
#'
#' @description
#' `DORAM()` tests two ecological endpoints for every selected taxon. The
#' occupancy endpoint concerns latent structural absence rather than observed
#' zero/nonzero prevalence. The ecological endpoint compares the read-relative
#' abundance `Y/N` across all samples, including zeros. A covariance-aware
#' two-degree-of-freedom joint test evaluates whether either endpoint changes.
#'
#' The function uses raw counts and the original sequencing depth directly. It
#' does not normalize counts, add pseudocounts, remove zeros, filter taxa, or
#' use resampling. Tests that cannot be calibrated safely are returned as
#' unavailable rather than being retried with a different procedure.
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
#'   is 1.
#' @param verbose If `TRUE`, print taxon-level progress. The default is
#'   `FALSE`.
#'
#' @return An object of class `DORAM`. Its main component, `results`, contains
#'   one row per taxon with raw p-values, BH-adjusted q-values, endpoint status,
#'   and the ecological difference (comparison minus reference). The `tests`
#'   component contains the corresponding statistics and degrees of freedom.
#'   `descriptives`, `diagnostics`, and restricted-null posterior matrices are
#'   retained for inspection but are not printed by default.
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
#' The ecological difference is the adjusted comparison-minus-reference
#' coefficient for all-sample `Y/N`. It is read-relative, not absolute
#' abundance, and no compositional-closure correction is performed. The joint
#' statistic uses the estimated cross-endpoint covariance; it is not a
#' combination of two marginal p-values.
#'
#' Multiplicity adjustment is performed separately for occupancy, ecological,
#' and joint tests with the Benjamini--Hochberg procedure. Every selected taxon
#' remains in its family denominator. An unavailable test retains `NA` for its
#' p- and q-values and is not a discovery.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120
#' sample_id <- sprintf("sample_%03d", seq_len(n))
#' reads <- round(exp(rnorm(n, log(5000), 0.35)))
#' diagnosis <- rep(c("control", "case"), each = n / 2)
#' age <- round(rnorm(n, 50, 12))
#' sex <- factor(rep(c("female", "male"), length.out = n))
#'
#' structural <- rbinom(n, 1, plogis(-0.8 + 0.01 * age))
#' probability <- plogis(-7 + 0.35 * (diagnosis == "case"))
#' y <- ifelse(structural == 1, 0, rbinom(n, reads, probability))
#' count_table <- cbind(focal_taxon = y)
#' rownames(count_table) <- sample_id
#' sample_data <- data.frame(
#'   diagnosis = diagnosis, reads = reads, age = age, sex = sex,
#'   row.names = sample_id
#' )
#'
#' fit <- DORAM(
#'   counts = count_table,
#'   metadata = sample_data,
#'   group = "diagnosis",
#'   reference = "control",
#'   library_size = "reads",
#'   covariates = c("age", "sex")
#' )
#' fit$results
#' }
#'
#' @importFrom stats constrOptim dbinom dnorm integrate lm.fit median pchisq
#' @importFrom stats plogis pnorm qlogis sd setNames uniroot
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

  tests <- do.call(rbind, lapply(fitted_taxa, `[[`, "rows"))
  rownames(tests) <- NULL
  tests <- .doram_adjust_results(tests)
  results <- .doram_primary_results(tests, selected_taxa)
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
  rownames(fit_diagnostics) <- rownames(endpoint_diagnostics) <- NULL

  structure(
    list(
      call = call,
      results = results,
      tests = tests,
      descriptives = descriptives,
      diagnostics = list(
        fit = fit_diagnostics,
        endpoint = endpoint_diagnostics
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
        sampling_design = "iid_random_design",
        strata = NULL,
        integration_level = "production",
        p_adjust_method = "BH",
        alpha = 0.05,
        keep_full_fits = FALSE,
        cores = cores
      )
    ),
    class = c("DORAM", "list")
  )
}

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
    "taxon", "q_occupancy", "ecological_difference", "q_ecological",
    "q_joint"
  )
  print(shown[, display_columns, drop = FALSE], row.names = FALSE, digits = 4)
  if (nrow(x$results) > nrow(shown)) {
    cat("... ", nrow(x$results) - nrow(shown),
        " more taxa; see fit$results\n", sep = "")
  }
  status_columns <- c("status_occupancy", "status_ecological", "status_joint")
  if (any(as.matrix(x$results[, status_columns, drop = FALSE]) != "ok")) {
    cat("Some tests were unavailable; see the status columns in fit$results.\n")
  }
  invisible(x)
}

#' @export
as.data.frame.DORAM <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}
