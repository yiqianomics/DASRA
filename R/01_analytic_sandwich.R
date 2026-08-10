# DORAM-D3 analytic efficient-score/sandwich inference overlay.
#
# This file intentionally does not modify the frozen count-native 0.7 source.
# It consumes already projected subject-level score rows and supplies the next
# candidate inference layer:
#   * iid/random-design global centering;
#   * independent-sampling-strata within-stratum centering;
#   * HC0 empirical meat;
#   * analytic chi-square upper-tail calibration.
#
# There is no bootstrap, multiplier, resampling seed, penalty, ridge,
# pseudoinverse, rank reduction, or data-dependent calibration choice here.

D3A_CONTRACT <- list(
  schema = "doram-d3-asymptotic-sandwich-overlay-0.1.0",
  primary_reference = "analytic_chisq",
  calibration_label = "asymptotic_empirical_sandwich_chisq",
  meat = "empirical_centered_hc0",
  condition_ratio_min = 1e-10,
  supported_dimensions = c(1L, 2L),
  count_source_schema = "doram-d3-count-native-prototype-0.7.0"
)

d3a_fail <- function(reason, stage, details = list()) {
  c(list(
    available = FALSE,
    testable = FALSE,
    p_value = NA_real_,
    log_p_value = NA_real_,
    statistic = NA_real_,
    df = NA_integer_,
    calibration = D3A_CONTRACT$calibration_label,
    failure_code = as.character(reason),
    reason = as.character(reason),
    stage = as.character(stage)
  ), details)
}

d3a_quadratic <- function(U, B,
                          condition_ratio_min =
                            D3A_CONTRACT$condition_ratio_min) {
  U <- as.numeric(U)
  B <- as.matrix(B)
  q <- length(U)
  if (!q %in% D3A_CONTRACT$supported_dimensions ||
      nrow(B) != q || ncol(B) != q ||
      any(!is.finite(U)) || any(!is.finite(B))) {
    return(list(ok = FALSE, reason = "invalid_quadratic_inputs"))
  }
  if (length(condition_ratio_min) != 1L ||
      !is.finite(condition_ratio_min) || condition_ratio_min <= 0 ||
      condition_ratio_min >= 1) {
    return(list(ok = FALSE, reason = "invalid_condition_ratio_contract"))
  }
  B <- (B + t(B)) / 2
  diagonal <- diag(B)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    return(list(ok = FALSE, reason = "nonpositive_meat_diagonal"))
  }
  scale <- sqrt(diagonal)
  correlation <- B / tcrossprod(scale)
  correlation <- (correlation + t(correlation)) / 2
  eigenvalues <- tryCatch(
    eigen(correlation, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, q)
  )
  condition_ratio <- if (length(eigenvalues) == q &&
      all(is.finite(eigenvalues)) && max(eigenvalues) > 0) {
    min(eigenvalues) / max(eigenvalues)
  } else {
    NA_real_
  }
  if (!is.finite(condition_ratio) || min(eigenvalues) <= 0 ||
      condition_ratio < condition_ratio_min) {
    return(list(
      ok = FALSE,
      reason = "unstable_standardized_meat",
      scale = scale,
      correlation = correlation,
      eigenvalues = eigenvalues,
      condition_ratio = condition_ratio
    ))
  }
  factor <- tryCatch(chol(correlation), error = function(e) NULL)
  if (is.null(factor)) {
    return(list(
      ok = FALSE,
      reason = "standardized_meat_cholesky_failed",
      scale = scale,
      correlation = correlation,
      eigenvalues = eigenvalues,
      condition_ratio = condition_ratio
    ))
  }
  U_scaled <- U / scale
  whitened <- tryCatch(
    forwardsolve(t(factor), U_scaled),
    error = function(e) rep(NA_real_, q)
  )
  statistic <- sum(whitened^2)
  if (!is.finite(statistic) || statistic < 0) {
    return(list(
      ok = FALSE,
      reason = "invalid_quadratic_statistic",
      scale = scale,
      correlation = correlation,
      eigenvalues = eigenvalues,
      condition_ratio = condition_ratio
    ))
  }
  list(
    ok = TRUE,
    statistic = statistic,
    scale = scale,
    correlation = correlation,
    eigenvalues = eigenvalues,
    condition_ratio = condition_ratio,
    factor = factor,
    U_scaled = U_scaled
  )
}

d3a_center_rows <- function(phi, strata = NULL,
                            sampling_design = "iid_random_design") {
  phi <- as.matrix(phi)
  n <- nrow(phi)
  if (length(sampling_design) != 1L ||
      !sampling_design %in% c("iid_random_design",
                              "independent_strata_fixed_counts")) {
    return(list(ok = FALSE, reason = "invalid_sampling_design"))
  }
  if (identical(sampling_design, "iid_random_design")) {
    if (!is.null(strata)) {
      return(list(ok = FALSE,
                  reason = "iid_design_must_not_supply_sampling_strata"))
    }
    centered <- sweep(phi, 2L, colMeans(phi), FUN = "-")
    return(list(
      ok = TRUE,
      mode = "iid_random_design",
      strata = rep("__iid__", n),
      stratum_names = "__iid__",
      stratum_sizes = setNames(n, "__iid__"),
      sample_stratum_proportions = setNames(1, "__iid__"),
      target_stratum_weights = setNames(1, "__iid__"),
      centered = centered,
      stratum_means = matrix(colMeans(phi), nrow = 1L,
                             dimnames = list("__iid__", colnames(phi)))
    ))
  }
  if (is.null(strata) || length(strata) != n || anyNA(strata)) {
    return(list(ok = FALSE, reason = "invalid_sampling_strata"))
  }
  strata <- as.character(strata)
  if (any(!nzchar(strata))) {
    return(list(ok = FALSE, reason = "empty_sampling_stratum_label"))
  }
  levels <- sort(unique(strata), method = "radix")
  indices <- lapply(levels, function(level) which(strata == level))
  sizes <- vapply(indices, length, integer(1))
  names(sizes) <- levels
  if (any(sizes < 2L)) {
    return(list(
      ok = FALSE,
      reason = "sampling_stratum_has_fewer_than_two_subjects",
      stratum_sizes = sizes
    ))
  }
  centered <- phi
  means <- matrix(NA_real_, nrow = length(levels), ncol = ncol(phi),
                  dimnames = list(levels, colnames(phi)))
  for (j in seq_along(indices)) {
    index <- indices[[j]]
    means[j, ] <- colMeans(phi[index, , drop = FALSE])
    centered[index, ] <- sweep(phi[index, , drop = FALSE], 2L,
                               means[j, ], FUN = "-")
  }
  list(
    ok = TRUE,
    mode = "independent_strata_fixed_counts",
    strata = strata,
    stratum_names = levels,
    stratum_sizes = sizes,
    sample_stratum_proportions = sizes / n,
    target_stratum_weights = sizes / n,
    centered = centered,
    stratum_means = means
  )
}

d3a_score_test <- function(phi, subject_id, strata = NULL,
                           sampling_design = "iid_random_design",
                           endpoint = "unspecified",
                           condition_ratio_min =
                             D3A_CONTRACT$condition_ratio_min) {
  phi <- as.matrix(phi)
  storage.mode(phi) <- "double"
  subject_id <- as.character(subject_id)
  n <- nrow(phi)
  q <- ncol(phi)
  if (n < 3L || !q %in% D3A_CONTRACT$supported_dimensions ||
      any(!is.finite(phi))) {
    return(d3a_fail("invalid_score_matrix", "geometry"))
  }
  if (length(subject_id) != n || anyNA(subject_id) ||
      any(!nzchar(subject_id)) || anyDuplicated(subject_id)) {
    return(d3a_fail("invalid_subject_id", "geometry"))
  }
  canonical_order <- order(subject_id, method = "radix")
  phi_canonical <- phi[canonical_order, , drop = FALSE]
  subject_canonical <- subject_id[canonical_order]
  strata_canonical <- if (is.null(strata)) NULL else strata[canonical_order]
  centered_result <- d3a_center_rows(
    phi_canonical, strata_canonical, sampling_design = sampling_design
  )
  if (!isTRUE(centered_result$ok)) {
    return(d3a_fail(
      centered_result$reason,
      "geometry",
      centered_result[setdiff(names(centered_result), c("ok", "reason"))]
    ))
  }
  centered_canonical <- centered_result$centered
  U <- colSums(phi_canonical)
  B <- crossprod(centered_canonical)
  quadratic <- d3a_quadratic(U, B, condition_ratio_min)
  if (!isTRUE(quadratic$ok)) {
    return(d3a_fail(
      quadratic$reason,
      "geometry",
      quadratic[setdiff(names(quadratic), c("ok", "reason"))]
    ))
  }
  standardized <- sweep(centered_canonical, 2L, quadratic$scale, FUN = "/")
  qr_standardized <- qr(standardized, LAPACK = TRUE)
  orthonormal <- qr.Q(qr_standardized, complete = FALSE)
  if (ncol(orthonormal) < q) {
    return(d3a_fail("standardized_score_qr_rank_failure", "geometry"))
  }
  leverage_canonical <- rowSums(
    orthonormal[, seq_len(q), drop = FALSE]^2
  )
  if (any(!is.finite(leverage_canonical)) ||
      abs(sum(leverage_canonical) - q) > 1e-8 * max(1, q)) {
    return(d3a_fail("invalid_score_leverage", "geometry"))
  }
  centered <- matrix(NA_real_, nrow = n, ncol = q,
                     dimnames = dimnames(phi))
  centered[canonical_order, ] <- centered_canonical
  leverage <- numeric(n)
  leverage[canonical_order] <- leverage_canonical
  log_p_value <- pchisq(
    quadratic$statistic, df = q, lower.tail = FALSE, log.p = TRUE
  )
  p_value <- pchisq(
    quadratic$statistic, df = q, lower.tail = FALSE, log.p = FALSE
  )
  valid_log_p <- is.finite(log_p_value) ||
    (is.infinite(log_p_value) && log_p_value < 0 && identical(p_value, 0))
  if (!valid_log_p || is.na(p_value) || p_value < 0 || p_value > 1) {
    return(d3a_fail("invalid_chisq_upper_tail", "calibration"))
  }
  # Leverage and score ESS are diagnostics.  There is no arbitrary hard cutoff.
  list(
    available = TRUE,
    testable = TRUE,
    reason = NA_character_,
    failure_code = "",
    stage = "complete",
    endpoint = as.character(endpoint),
    calibration = D3A_CONTRACT$calibration_label,
    primary_reference = D3A_CONTRACT$primary_reference,
    statistic = quadratic$statistic,
    df = as.integer(q),
    p_value = p_value,
    log_p_value = log_p_value,
    n = n,
    q = q,
    subject_id = subject_id,
    canonical_subject_id = subject_canonical,
    canonical_order = canonical_order,
    phi = phi,
    centered = centered,
    U = U,
    B = B,
    meat = D3A_CONTRACT$meat,
    sampling_mode = centered_result$mode,
    sampling_strata = if (is.null(strata)) NULL else as.character(strata),
    stratum_names = centered_result$stratum_names,
    stratum_sizes = centered_result$stratum_sizes,
    sample_stratum_proportions =
      centered_result$sample_stratum_proportions,
    target_stratum_weights = centered_result$target_stratum_weights,
    stratum_means = centered_result$stratum_means,
    centered_stratum_residual = max(vapply(
      split(seq_len(n), centered_result$strata),
      function(index) max(abs(colSums(centered_canonical[index, , drop = FALSE]))),
      numeric(1)
    )),
    scale = quadratic$scale,
    standardized_meat = quadratic$correlation,
    standardized_meat_eigenvalues = quadratic$eigenvalues,
    condition_ratio = quadratic$condition_ratio,
    condition_number = 1 / quadratic$condition_ratio,
    leverage = leverage,
    max_leverage = max(leverage),
    score_ess = q^2 / sum(leverage^2)
  )
}

d3a_fwl_rows <- function(y, N, group, X, subject_id) {
  y <- as.numeric(y)
  N <- as.numeric(N)
  group <- as.numeric(group)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  subject_id <- as.character(subject_id)
  n <- length(y)
  if (length(N) != n || length(group) != n || nrow(X) != n ||
      length(subject_id) != n || anyNA(subject_id) ||
      any(!nzchar(subject_id)) || anyDuplicated(subject_id) ||
      n < ncol(X) + 3L || any(!is.finite(c(y, N, group, X))) ||
      any(N <= 0) || any(N != round(N)) || any(y != round(y)) ||
      any(y < 0 | y > N) || any(!group %in% c(0, 1)) ||
      qr(X)$rank != ncol(X)) {
    return(d3a_fail("invalid_ecological_inputs", "ecological_rows"))
  }
  has_intercept <- any(vapply(
    seq_len(ncol(X)),
    function(j) max(abs(X[, j] - 1)) <= 100 * .Machine$double.eps,
    logical(1)
  ))
  if (!has_intercept) {
    return(d3a_fail("missing_ecological_intercept", "ecological_rows"))
  }
  R <- y / N
  qx <- qr(X)
  theta_null <- qr.coef(qx, R)
  group_residual <- group - qr.fitted(qx, group)
  group_ss <- sum(group_residual^2)
  if (!is.finite(group_ss) || group_ss <= 1e-12 * max(1, sum(group^2))) {
    return(d3a_fail("no_residual_group_variation", "ecological_rows"))
  }
  residual_null <- R - drop(X %*% theta_null)
  phi <- as.numeric(group_residual * residual_null)
  if (any(!is.finite(phi))) {
    return(d3a_fail("nonfinite_ecological_score_rows", "ecological_rows"))
  }
  full_design <- cbind(X, group = group)
  qfull <- qr(full_design)
  if (qfull$rank != ncol(full_design)) {
    return(d3a_fail("rank_deficient_full_ecological_design", "ecological_rows"))
  }
  full_coef <- qr.coef(qfull, R)
  score_equation_error <- max(abs(crossprod(X, residual_null)))
  fwl_sum_error <- abs(sum(phi) - sum(group * residual_null))
  tolerance <- 1e-9 * max(1, sum(abs(phi)), sum(abs(group * residual_null)),
                          sum(abs(R)), sum(abs(residual_null)))
  if (!is.finite(score_equation_error) || !is.finite(fwl_sum_error) ||
      score_equation_error > tolerance || fwl_sum_error > tolerance) {
    return(d3a_fail("ecological_fwl_identity_failure", "ecological_rows"))
  }
  list(
    available = TRUE,
    testable = NA,
    reason = NA_character_,
    failure_code = "",
    stage = "ecological_rows",
    endpoint = "ecological",
    subject_id = subject_id,
    phi = matrix(phi, ncol = 1L,
                 dimnames = list(NULL, "ecological")),
    R = R,
    theta_null = theta_null,
    group_residual = group_residual,
    residual_null = residual_null,
    beta_hat = unname(full_coef[[ncol(full_design)]]),
    score_equation_error = score_equation_error,
    fwl_sum_error = fwl_sum_error
  )
}

d3a_projected_endpoint_test <- function(
    projected, endpoint, strata = NULL,
    sampling_design = "iid_random_design") {
  if (!isTRUE(projected$available)) return(projected)
  if (is.null(projected$phi) || is.null(projected$subject_id)) {
    return(d3a_fail("projected_endpoint_missing_rows", "projection"))
  }
  result <- d3a_score_test(
    phi = projected$phi,
    subject_id = projected$subject_id,
    strata = strata,
    sampling_design = sampling_design,
    endpoint = endpoint
  )
  if (isTRUE(result$available)) result$projected <- projected
  result
}

d3a_joint_test <- function(occupancy_projected, ecological_rows,
                           strata = NULL,
                           sampling_design = "iid_random_design") {
  if (!isTRUE(occupancy_projected$available)) {
    return(d3a_fail(
      paste0("occupancy_component_",
             if (!is.null(occupancy_projected$failure_code))
               occupancy_projected$failure_code else "unavailable"),
      "joint_formation"
    ))
  }
  if (!isTRUE(ecological_rows$available)) {
    return(d3a_fail(
      paste0("ecological_component_",
             if (!is.null(ecological_rows$failure_code))
               ecological_rows$failure_code else "unavailable"),
      "joint_formation"
    ))
  }
  occ_id <- as.character(occupancy_projected$subject_id)
  eco_id <- as.character(ecological_rows$subject_id)
  if (!identical(occ_id, eco_id)) {
    return(d3a_fail("joint_subject_alignment_mismatch", "joint_formation"))
  }
  occ_phi <- as.numeric(occupancy_projected$phi)
  eco_phi <- as.numeric(ecological_rows$phi)
  if (length(occ_phi) != length(eco_phi)) {
    return(d3a_fail("joint_score_length_mismatch", "joint_formation"))
  }
  d3a_score_test(
    phi = cbind(occupancy = occ_phi, ecological = eco_phi),
    subject_id = eco_id,
    strata = strata,
    sampling_design = sampling_design,
    endpoint = "joint"
  )
}

d3a_count_source_is_compatible <- function() {
  exists("CN_D3_CONTRACT", inherits = TRUE) &&
    identical(CN_D3_CONTRACT$schema_version,
              D3A_CONTRACT$count_source_schema) &&
    exists("cn_fit_occupancy_null", mode = "function", inherits = TRUE) &&
    exists("cn_project_occupancy_score", mode = "function", inherits = TRUE)
}

d3a_fit_endpoints <- function(count_data, X, strata = NULL,
                              sampling_design = "iid_random_design",
                              integration_level = "production") {
  if (!d3a_count_source_is_compatible()) {
    stop("frozen count-native 0.7 source is not loaded or is incompatible",
         call. = FALSE)
  }
  fit <- cn_fit_occupancy_null(
    count_data, integration_level = integration_level
  )
  occupancy_projected <- if (isTRUE(fit$available)) {
    cn_project_occupancy_score(fit, count_data)
  } else {
    fit
  }
  if (isTRUE(occupancy_projected$available)) {
    occupancy_projected$phi <- matrix(
      occupancy_projected$phi, ncol = 1L,
      dimnames = list(NULL, "occupancy")
    )
  }
  ecological_rows <- d3a_fwl_rows(
    y = count_data$y,
    N = count_data$N,
    group = count_data$g,
    X = X,
    subject_id = count_data$subject_id
  )
  occupancy <- d3a_projected_endpoint_test(
    occupancy_projected, "occupancy", strata, sampling_design
  )
  ecological <- d3a_projected_endpoint_test(
    ecological_rows, "ecological", strata, sampling_design
  )
  joint <- d3a_joint_test(
    occupancy_projected, ecological_rows, strata, sampling_design
  )
  list(
    fit = fit,
    occupancy_projected = occupancy_projected,
    ecological_rows = ecological_rows,
    occupancy = occupancy,
    ecological = ecological,
    joint = joint
  )
}
