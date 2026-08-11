# Compiled adaptive quadrature for the present-count component.
#
# The count likelihood, posterior-mode calculation, error checks, and
# inferential calculations remain in R. The compiled kernel evaluates and
# adaptively refines the vector Gauss--Kronrod panels.

cn_integrated_present_compiled <- function(
    y, N, eta, sigma, level = c("production", "audit")) {
  library_envelope <- cn_library_envelope_check(N)
  eta_envelope <- cn_eta_envelope_check(eta)
  if (length(y) != 1L || !is.finite(y) || y < 0 || y != round(y) ||
      y > N || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
      !isTRUE(library_envelope$ok) || !isTRUE(eta_envelope$ok)) {
    cn_stop("integration inputs are outside the supported numerical range")
  }

  control <- cn_integration_control(level)
  mode <- cn_posterior_mode(y, N, eta, sigma)
  h_mode <- eta + sigma * mode
  p_mode <- plogis(h_mode)
  q_mode <- plogis(-h_mode)
  mode_residual <- sigma * cn_binomial_residual(y, N, h_mode) - mode
  mode_backward_error <- abs(mode_residual) / max(
    1, abs(sigma * y), abs(sigma * N * p_mode), abs(mode)
  )
  if (!is.finite(mode_backward_error) ||
      mode_backward_error > 100 * .Machine$double.eps) {
    cn_stop("posterior mode exceeds the backward-error contract")
  }

  local_scale <- 1 / sqrt(1 + sigma^2 * N * p_mode * q_mode)
  if (!is.finite(local_scale) || local_scale <= 0) {
    cn_stop("invalid integration scale")
  }

  tail_error <- cn_tail_moment_bounds(control$tail_radius, mode_residual)
  kernel <- doram_integrate_moments_cpp(
    y, N, eta, sigma, mode, p_mode, mode_residual, local_scale,
    control$relative_tolerance, control$tail_radius,
    control$maximum_panels, tail_error
  )
  value <- kernel$value
  error <- kernel$error
  transformed_error <- setNames(
    kernel$transformed_error,
    c("log_K", "score_eta", "score_log_sigma")
  )
  denominator <- value[1L]
  mean_u <- value[2L] / denominator
  second_u <- value[3L] / denominator
  score_eta <- (mode + mean_u) / sigma
  score_log_sigma <- mode^2 + 2 * mode * mean_u + second_u - 1
  anchor <- cn_anchor_at_mode(y, N, eta, sigma, mode)
  if (any(!is.finite(c(
    denominator, mean_u, second_u, score_eta, score_log_sigma, anchor,
    transformed_error
  )))) {
    cn_stop("nonfinite vector integration result")
  }

  list(
    log_K = anchor + log(denominator),
    score_eta = score_eta,
    score_log_sigma = score_log_sigma,
    mode = mode,
    mode_residual = mode_residual,
    mode_backward_error = mode_backward_error,
    scale = local_scale,
    integration_level = control$level,
    integration_absolute_error = setNames(
      error,
      c("denominator", "first_centered_moment", "second_centered_moment")
    ),
    integration_transformed_error = transformed_error,
    integration_panels = kernel$panels,
    integration_evaluations = kernel$evaluations,
    integration_tail_bound = tail_error
  )
}
