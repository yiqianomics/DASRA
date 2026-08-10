# Internal validation, alignment, fitting, and result helpers.
#
# The public interface is intentionally small. These helpers keep the user-
# facing code readable while preserving the verified taxon-level likelihood
# and analytic score construction without modification.

.doram_abort <- function(...) stop(..., call. = FALSE)

.doram_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) != 1L || !is.finite(x)) return(default)
  as.numeric(x)
}

.doram_log_scalar <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return(NA_real_)
  x <- as.numeric(x)
  if (is.finite(x) || (is.infinite(x) && x < 0)) x else NA_real_
}

.doram_integer_scalar <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) != 1L || !is.finite(x)) return(default)
  as.integer(x)
}

.doram_text_scalar <- function(x, default = "") {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- as.character(x[[1L]])
  if (nzchar(value)) value else default
}

.doram_failure_reason <- function(x, default = "unavailable") {
  reason <- .doram_text_scalar(x$failure_code)
  if (!nzchar(reason)) reason <- .doram_text_scalar(x$reason)
  if (!nzchar(reason)) reason <- default
  reason
}

.doram_counts_matrix <- function(counts) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    .doram_abort("counts must be a samples-by-taxa matrix or data frame")
  }
  if (length(dim(counts)) != 2L || !nrow(counts) || !ncol(counts)) {
    .doram_abort("counts must contain at least one sample and one taxon")
  }
  if (is.data.frame(counts) && !all(vapply(counts, is.numeric, logical(1)))) {
    .doram_abort("every counts column must be numeric")
  }
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) .doram_abort("counts must be numeric")
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != round(counts))) {
    .doram_abort("counts must contain finite nonnegative integers")
  }

  sample_id <- rownames(counts)
  if (is.null(sample_id) || anyNA(sample_id) || any(!nzchar(sample_id)) ||
      anyDuplicated(sample_id)) {
    .doram_abort("counts must have unique, nonempty sample row names")
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- sprintf("taxon_%05d", seq_len(ncol(counts)))
  }
  if (anyNA(colnames(counts)) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(colnames(counts))) {
    .doram_abort("counts must have unique, nonempty taxon column names")
  }
  counts
}

.doram_metadata <- function(metadata, sample_id) {
  if (!is.data.frame(metadata)) {
    .doram_abort("metadata must be a data frame with samples in rows")
  }
  if (!nrow(metadata) || !ncol(metadata)) {
    .doram_abort("metadata must contain at least one sample and one variable")
  }
  metadata_id <- rownames(metadata)
  if (is.null(metadata_id) || anyNA(metadata_id) || any(!nzchar(metadata_id)) ||
      anyDuplicated(metadata_id)) {
    .doram_abort("metadata must have unique, nonempty sample row names")
  }
  if (is.null(colnames(metadata)) || anyNA(colnames(metadata)) ||
      any(!nzchar(colnames(metadata))) || anyDuplicated(colnames(metadata))) {
    .doram_abort("metadata must have unique, nonempty column names")
  }

  missing_metadata <- setdiff(sample_id, metadata_id)
  extra_metadata <- setdiff(metadata_id, sample_id)
  if (length(missing_metadata) || length(extra_metadata)) {
    detail <- character()
    if (length(missing_metadata)) {
      detail <- c(detail, paste0(
        "missing from metadata: ",
        paste(utils::head(missing_metadata, 5L), collapse = ", ")
      ))
    }
    if (length(extra_metadata)) {
      detail <- c(detail, paste0(
        "missing from counts: ",
        paste(utils::head(extra_metadata, 5L), collapse = ", ")
      ))
    }
    .doram_abort(
      "counts and metadata must contain exactly the same sample IDs (",
      paste(detail, collapse = "; "), ")"
    )
  }

  reordered <- !identical(metadata_id, sample_id)
  metadata <- metadata[match(sample_id, metadata_id), , drop = FALSE]
  if (!identical(rownames(metadata), sample_id)) {
    .doram_abort("internal sample alignment failed")
  }
  list(data = metadata, reordered = reordered)
}

.doram_column_name <- function(x, metadata, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .doram_abort(label, " must be one metadata column name")
  }
  if (!x %in% colnames(metadata)) {
    .doram_abort(label, " column '", x, "' was not found in metadata")
  }
  x
}

.doram_group <- function(metadata, group, reference = NULL) {
  group <- .doram_column_name(group, metadata, "group")
  raw <- metadata[[group]]
  if (!(is.numeric(raw) || is.logical(raw) || is.factor(raw) ||
        is.character(raw))) {
    .doram_abort("the group column must be numeric, logical, factor, or character")
  }
  if (anyNA(raw)) .doram_abort("the group column cannot contain missing values")

  values <- as.character(raw)
  observed <- unique(values)
  if (length(observed) != 2L) {
    .doram_abort("the group column must contain exactly two observed levels")
  }
  if (is.null(reference)) {
    if (is.logical(raw) && setequal(unique(raw), c(FALSE, TRUE))) {
      reference <- "FALSE"
    } else if (is.numeric(raw) && setequal(unique(as.numeric(raw)), c(0, 1))) {
      reference <- "0"
    } else {
      .doram_abort(
        "reference is required when the group column is character, factor, ",
        "or is not numeric 0/1"
      )
    }
  }
  if (length(reference) != 1L || is.na(reference)) {
    .doram_abort("reference must name exactly one observed group level")
  }
  reference <- as.character(reference)
  if (!reference %in% observed) {
    .doram_abort("reference '", reference, "' is not an observed group level")
  }
  comparison <- setdiff(observed, reference)
  if (length(comparison) != 1L) .doram_abort("could not determine the comparison group")
  encoded <- as.integer(values == comparison)
  names(encoded) <- rownames(metadata)
  list(
    column = group,
    encoded = encoded,
    reference = reference,
    comparison = comparison,
    labels = values
  )
}

.doram_library_size <- function(metadata, library_size, counts) {
  library_size <- .doram_column_name(library_size, metadata, "library_size")
  depth <- metadata[[library_size]]
  if (!is.numeric(depth) || is.logical(depth) || length(depth) != nrow(counts) ||
      any(!is.finite(depth)) || any(depth <= 0) || any(depth != round(depth))) {
    .doram_abort(
      "the library_size column must contain one finite positive integer per sample"
    )
  }
  depth <- as.numeric(depth)
  names(depth) <- rownames(metadata)
  upper <- CN_D3_CONTRACT$integration$certified_library_size[["upper"]]
  if (any(depth > upper)) {
    .doram_abort("library_size exceeds the certified count-likelihood envelope")
  }
  if (any(apply(counts, 1L, max) > depth)) {
    .doram_abort("no taxon count may exceed its original library size")
  }
  if (any(rowSums(counts) > depth)) {
    .doram_abort("taxon counts cannot sum to more than the original library size")
  }
  list(column = library_size, values = depth)
}

.doram_covariates <- function(metadata, covariates, group) {
  sample_id <- rownames(metadata)
  if (is.null(covariates)) {
    return(list(columns = character(), matrix = NULL, design_names = character(),
                levels = list()))
  }
  if (!is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates)) || anyDuplicated(covariates)) {
    .doram_abort("covariates must be NULL or unique metadata column names")
  }
  missing_columns <- setdiff(covariates, colnames(metadata))
  if (length(missing_columns)) {
    .doram_abort(
      "covariate columns were not found in metadata: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (group %in% covariates) {
    .doram_abort("the tested group column cannot also be included in covariates")
  }

  raw <- metadata[, covariates, drop = FALSE]
  supported <- vapply(
    raw,
    function(x) is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x),
    logical(1)
  )
  if (!all(supported)) {
    .doram_abort(
      "covariates must be numeric, logical, factor, or character columns: ",
      paste(covariates[!supported], collapse = ", ")
    )
  }
  if (anyNA(raw)) .doram_abort("selected covariates cannot contain missing values")

  coding_levels <- list()
  model_data <- raw
  for (name in covariates) {
    x <- raw[[name]]
    if (is.character(x)) {
      levels_x <- sort(unique(x), method = "radix")
      x <- factor(x, levels = levels_x)
    } else if (is.factor(x)) {
      levels_x <- levels(droplevels(x))
      x <- factor(as.character(x), levels = levels_x, ordered = FALSE)
    } else if (is.logical(x)) {
      levels_x <- c("FALSE", "TRUE")[c(FALSE, TRUE) %in% unique(x)]
      x <- factor(as.character(x), levels = levels_x, ordered = FALSE)
    } else {
      levels_x <- NULL
    }
    if (!is.null(levels_x)) {
      if (length(levels_x) < 2L) {
        .doram_abort("categorical covariate '", name, "' has only one observed level")
      }
      coding_levels[[name]] <- levels_x
    }
    model_data[[name]] <- x
  }

  contrasts_arg <- lapply(coding_levels, function(levels_x) {
    stats::contr.treatment(levels_x, base = 1L)
  })
  design <- tryCatch(
    stats::model.matrix(
      ~ .,
      data = model_data,
      contrasts.arg = if (length(contrasts_arg)) contrasts_arg else NULL,
      na.action = stats::na.fail
    ),
    error = function(e) e
  )
  if (inherits(design, "error")) {
    .doram_abort("could not construct the covariate design: ", conditionMessage(design))
  }
  if (!identical(rownames(design), sample_id) || nrow(design) != nrow(metadata)) {
    .doram_abort("covariate construction changed sample identity or order")
  }
  intercept <- which(colnames(design) == "(Intercept)")
  if (length(intercept)) design <- design[, -intercept, drop = FALSE]
  if (!ncol(design)) .doram_abort("covariates produced an empty design matrix")
  storage.mode(design) <- "double"
  if (any(!is.finite(design))) {
    .doram_abort("covariates must produce a finite numeric design")
  }
  if (anyNA(colnames(design)) || any(!nzchar(colnames(design))) ||
      anyDuplicated(colnames(design))) {
    .doram_abort("covariates produced duplicate or invalid design column names")
  }
  rownames(design) <- sample_id
  list(
    columns = covariates,
    matrix = design,
    design_names = colnames(design),
    levels = coding_levels
  )
}

.doram_taxa <- function(taxa, taxon_names) {
  if (is.null(taxa)) return(taxon_names)
  if (is.numeric(taxa) && !is.logical(taxa)) {
    if (!length(taxa) || any(!is.finite(taxa)) || any(taxa != round(taxa)) ||
        any(taxa < 1L) || any(taxa > length(taxon_names))) {
      .doram_abort("numeric taxa must be valid integer column indices")
    }
    taxa <- taxon_names[as.integer(taxa)]
  } else {
    if (!is.character(taxa) || !length(taxa) || anyNA(taxa) ||
        any(!nzchar(taxa))) {
      .doram_abort("taxa must be NULL, taxon names, or integer column indices")
    }
    missing_taxa <- setdiff(taxa, taxon_names)
    if (length(missing_taxa)) {
      .doram_abort("unknown taxa: ", paste(missing_taxa, collapse = ", "))
    }
  }
  if (anyDuplicated(taxa)) .doram_abort("taxa must select unique columns")
  taxa
}

.doram_cores <- function(cores) {
  if (!is.numeric(cores) || is.logical(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores != round(cores) || cores < 1L) {
    .doram_abort("cores must be one positive integer")
  }
  as.integer(cores)
}

.doram_logical_scalar <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .doram_abort(label, " must be one logical value (TRUE or FALSE)")
  }
  x
}

.doram_endpoint_row <- function(result, taxon, endpoint, estimate = NA_real_) {
  available <- isTRUE(result$available) &&
    length(result$p_value) == 1L && is.finite(result$p_value) &&
    length(result$statistic) == 1L && is.finite(result$statistic) &&
    length(result$df) == 1L && is.finite(result$df)
  data.frame(
    taxon = taxon,
    endpoint = endpoint,
    available = available,
    estimate = if (available) estimate else NA_real_,
    statistic = if (available) .doram_scalar(result$statistic) else NA_real_,
    df = if (available) .doram_integer_scalar(result$df) else NA_integer_,
    p_value = if (available) .doram_scalar(result$p_value) else NA_real_,
    log_p_value = if (available) .doram_log_scalar(result$log_p_value) else NA_real_,
    q_value = NA_real_,
    status = if (available) "ok" else .doram_failure_reason(result),
    stringsAsFactors = FALSE
  )
}

.doram_endpoint_diagnostic <- function(result, taxon, endpoint) {
  available <- isTRUE(result$available)
  data.frame(
    taxon = taxon,
    endpoint = endpoint,
    available = available,
    status = if (available) "ok" else .doram_failure_reason(result),
    max_leverage = .doram_scalar(result$max_leverage),
    score_ess = .doram_scalar(result$score_ess),
    condition_ratio = .doram_scalar(result$condition_ratio),
    calibration = {
      value <- .doram_text_scalar(result$calibration)
      if (nzchar(value)) value else NA_character_
    },
    stage = .doram_text_scalar(result$stage),
    stringsAsFactors = FALSE
  )
}

.doram_exception_fit <- function(taxon, n, message) {
  rows <- do.call(rbind, lapply(c("occupancy", "ecological", "joint"), function(endpoint) {
    data.frame(
      taxon = taxon, endpoint = endpoint, available = FALSE,
      estimate = NA_real_, statistic = NA_real_, df = NA_integer_,
      p_value = NA_real_, log_p_value = NA_real_, q_value = NA_real_,
      status = "taxon_fit_exception", stringsAsFactors = FALSE
    )
  }))
  endpoint_diagnostics <- do.call(rbind, lapply(
    c("occupancy", "ecological", "joint"),
    function(endpoint) data.frame(
      taxon = taxon, endpoint = endpoint, available = FALSE,
      status = "taxon_fit_exception", max_leverage = NA_real_,
      score_ess = NA_real_, condition_ratio = NA_real_,
      calibration = NA_character_, stage = "taxon_wrapper",
      stringsAsFactors = FALSE
    )
  ))
  list(
    rows = rows,
    gamma_null = rep(NA_real_, n),
    rho_null = rep(NA_real_, n),
    fit_diagnostic = data.frame(
      taxon = taxon, available = FALSE, status = "taxon_fit_exception",
      error_message = as.character(message), stationarity = NA_real_,
      converged_starts = NA_integer_,
      integration_max_subject_loglik_difference = NA_real_,
      integration_max_subject_score_difference = NA_real_,
      integration_max_total_score_difference = NA_real_,
      sensitivity_asymmetry = NA_real_,
      nuisance_condition_number = NA_real_,
      delta_null_constraint = 0,
      stringsAsFactors = FALSE
    ),
    endpoint_diagnostics = endpoint_diagnostics
  )
}

.doram_fit_one <- function(taxon, y, library_size, group, covariates,
                           subject_id, X) {
  count_data <- cn_prepare_data(
    y = y,
    N = library_size,
    g = group,
    z = covariates,
    subject_id = subject_id
  )
  endpoints <- tryCatch(
    d3a_fit_endpoints(
      count_data = count_data,
      X = X,
      strata = NULL,
      sampling_design = "iid_random_design",
      integration_level = "production"
    ),
    error = function(e) e
  )
  if (inherits(endpoints, "error")) {
    return(.doram_exception_fit(taxon, length(y), conditionMessage(endpoints)))
  }

  ecological_estimate <- if (isTRUE(endpoints$ecological_rows$available)) {
    .doram_scalar(endpoints$ecological_rows$beta_hat)
  } else {
    NA_real_
  }
  rows <- rbind(
    .doram_endpoint_row(endpoints$occupancy, taxon, "occupancy"),
    .doram_endpoint_row(
      endpoints$ecological, taxon, "ecological",
      estimate = ecological_estimate
    ),
    .doram_endpoint_row(endpoints$joint, taxon, "joint")
  )
  endpoint_diagnostics <- rbind(
    .doram_endpoint_diagnostic(endpoints$occupancy, taxon, "occupancy"),
    .doram_endpoint_diagnostic(endpoints$ecological, taxon, "ecological"),
    .doram_endpoint_diagnostic(endpoints$joint, taxon, "joint")
  )

  fit_available <- isTRUE(endpoints$fit$available)
  gamma_null <- rho_null <- rep(NA_real_, length(y))
  if (fit_available) {
    gamma_null <- as.numeric(endpoints$fit$evaluation$gamma)
    rho_null <- as.numeric(endpoints$fit$evaluation$rho)
  }
  projected_available <- isTRUE(endpoints$occupancy_projected$available)
  fit_status <- if (fit_available) "ok" else .doram_failure_reason(endpoints$fit)
  fit_diagnostic <- data.frame(
    taxon = taxon,
    available = fit_available,
    status = fit_status,
    error_message = "",
    stationarity = if (fit_available) {
      .doram_scalar(endpoints$fit$stationarity)
    } else NA_real_,
    converged_starts = if (fit_available) {
      .doram_integer_scalar(endpoints$fit$n_converged)
    } else NA_integer_,
    integration_max_subject_loglik_difference = if (fit_available) {
      .doram_scalar(endpoints$fit$integration_audit$maximum_subject_loglik_difference)
    } else NA_real_,
    integration_max_subject_score_difference = if (fit_available) {
      .doram_scalar(endpoints$fit$integration_audit$maximum_subject_score_difference)
    } else NA_real_,
    integration_max_total_score_difference = if (fit_available) {
      .doram_scalar(endpoints$fit$integration_audit$maximum_total_score_difference)
    } else NA_real_,
    sensitivity_asymmetry = if (projected_available) {
      .doram_scalar(endpoints$occupancy_projected$sensitivity$asymmetry)
    } else NA_real_,
    nuisance_condition_number = if (projected_available) {
      .doram_scalar(endpoints$occupancy_projected$nuisance_health$condition_number)
    } else NA_real_,
    delta_null_constraint = 0,
    stringsAsFactors = FALSE
  )
  list(
    rows = rows,
    gamma_null = gamma_null,
    rho_null = rho_null,
    fit_diagnostic = fit_diagnostic,
    endpoint_diagnostics = endpoint_diagnostics
  )
}

.doram_adjust_results <- function(results) {
  for (endpoint in c("occupancy", "ecological", "joint")) {
    index <- which(results$endpoint == endpoint)
    results$q_value[index] <- stats::p.adjust(
      results$p_value[index], method = "BH", n = length(index)
    )
  }
  results
}

.doram_primary_results <- function(tests, taxa) {
  take <- function(endpoint) {
    x <- tests[tests$endpoint == endpoint, , drop = FALSE]
    x[match(taxa, x$taxon), , drop = FALSE]
  }
  occupancy <- take("occupancy")
  ecological <- take("ecological")
  joint <- take("joint")
  data.frame(
    taxon = taxa,
    p_occupancy = occupancy$p_value,
    q_occupancy = occupancy$q_value,
    status_occupancy = occupancy$status,
    ecological_difference = ecological$estimate,
    p_ecological = ecological$p_value,
    q_ecological = ecological$q_value,
    status_ecological = ecological$status,
    p_joint = joint$p_value,
    q_joint = joint$q_value,
    status_joint = joint$status,
    stringsAsFactors = FALSE
  )
}

.doram_descriptives <- function(counts, library_size, group, taxa) {
  do.call(rbind, lapply(taxa, function(taxon) {
    y <- counts[, taxon]
    reference <- group == 0L
    comparison <- group == 1L
    data.frame(
      taxon = taxon,
      n_samples = length(y),
      zeros_reference = sum(y[reference] == 0L),
      zeros_comparison = sum(y[comparison] == 0L),
      positives_reference = sum(y[reference] > 0L),
      positives_comparison = sum(y[comparison] > 0L),
      mean_y_over_n_reference = mean(y[reference] / library_size[reference]),
      mean_y_over_n_comparison = mean(y[comparison] / library_size[comparison]),
      stringsAsFactors = FALSE
    )
  }))
}

.doram_map_taxa <- function(taxa, FUN, cores, verbose) {
  serial <- cores == 1L || length(taxa) == 1L ||
    identical(.Platform$OS.type, "windows")
  if (serial) {
    if (cores > 1L && identical(.Platform$OS.type, "windows")) {
      warning("parallel taxon fitting is unavailable on Windows; using one core",
              call. = FALSE)
    }
    out <- vector("list", length(taxa))
    for (i in seq_along(taxa)) {
      if (isTRUE(verbose)) {
        message(sprintf("DORAM: fitting taxon %d/%d (%s)",
                        i, length(taxa), taxa[[i]]))
      }
      out[[i]] <- FUN(taxa[[i]])
    }
    return(out)
  }
  parallel::mclapply(
    taxa,
    FUN,
    mc.cores = min(cores, length(taxa)),
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
}
