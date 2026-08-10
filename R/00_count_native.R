# Standalone DORAM-D3 count-native prototype.
#
# This file is deliberately independent of the installed dash package and of
# the frozen D1/D2 verification code.  It implements only the taxonwise
# zero-inflated binomial-logistic-normal likelihood and the restricted
# occupancy score H0: delta = 0, with zeta and kappa left unrestricted.

if (!requireNamespace("statmod", quietly = TRUE)) {
  stop("count_native.R requires the statmod package", call. = FALSE)
}
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("count_native.R requires the numDeriv package", call. = FALSE)
}

CN_D3_CONTRACT <- list(
  schema_version = "doram-d3-count-native-prototype-0.7.0",
  integration = list(
    engine = "mode_centered_vector_gk15",
    production_relative_tolerance = 1e-8,
    audit_relative_tolerance = 1e-10,
    tail_radius = 12,
    maximum_panels = 20000L,
    diagnostic_agh_nodes = 31L,
    maximum_subject_loglik_difference = 1e-7,
    maximum_subject_score_difference = 1e-6,
    maximum_total_score_difference = 1e-5,
    certified_library_size = c(lower = 1, upper = 1e9),
    certified_eta_predictor = c(lower = -20, upper = 10),
    mode_tolerance = 1e-12,
    mode_max_iterations = 100L
  ),
  parameter_bounds = list(
    structural = c(lower = -10, upper = 10),
    location = c(lower = -20, upper = 10),
    log_scale = c(lower = log(0.1), upper = log(8))
  ),
  precheck = list(
    rank_tolerance = 1e-10,
    minimum_zero_per_group = 2L,
    minimum_positive_per_group = 4L
  ),
  optimizer = list(
    method = "constrOptim-BFGS",
    inner_method = "BFGS",
    maxit = 3000L,
    relative_tolerance = 1e-12,
    barrier_mu = 1e-6,
    barrier_outer_iterations = 100L,
    barrier_outer_tolerance = 1e-10,
    strict_feasibility_relative_margin = 1e-8,
    invalid_objective_penalty = 1e100,
    stationarity_relative_score_norm = 1e-5,
    boundary_relative_distance = 1e-6,
    near_mode_objective_relative_tolerance = 1e-8,
    distinct_mode_parameter_tolerance = 1e-3,
    distinct_mode_fitted_tolerance = 1e-4,
    interior_polish = list(
      direction = "scaled_score_only",
      score_scale_floor = 100 * .Machine$double.eps,
      fraction_to_boundary = 0.995,
      maximum_iterations = 25L,
      maximum_backtracks = 50L,
      armijo = 1e-4
    )
  ),
  sensitivity = list(
    asymmetry_tolerance = 2e-5,
    scaled_minimum_eigenvalue = 1e-8,
    maximum_condition_number = 1e8,
    solve_relative_agreement = 1e-6,
    score_identity_relative_tolerance = 1e-10
  ),
  applicability = list(
    maximum_centered_score_leverage = 1 / 3,
    minimum_score_effective_sample_size = 12,
    minimum_second_moment = 100 * .Machine$double.eps
  )
)

.cn_gh_cache <- new.env(parent = emptyenv())

cn_stop <- function(...) stop(..., call. = FALSE)

cn_failure <- function(code, ..., stage = NA_character_) {
  structure(
    c(list(
      available = FALSE,
      p_value = NA_real_,
      p_multiplier = NA_real_,
      p_chisq = NA_real_,
      statistic = NA_real_,
      failure_code = as.character(code),
      reason = as.character(code),
      stage = as.character(stage)
    ), list(...)),
    class = c("cn_d3_refusal", "list")
  )
}

cn_softplus <- function(x) {
  out <- numeric(length(x))
  hi <- x > 0
  out[hi] <- x[hi] + log1p(exp(-x[hi]))
  out[!hi] <- log1p(exp(x[!hi]))
  out
}

cn_log_sum_exp <- function(x) {
  if (!length(x) || anyNA(x)) cn_stop("log-sum-exp received an invalid vector")
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}

cn_logspace_add <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  both_negative_infinite <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
  out[both_negative_infinite] <- -Inf
  out
}

cn_log_binom_kernel <- function(z, y, N, eta, sigma) {
  h <- eta + sigma * z
  log_p <- -cn_softplus(-h)
  log_one_minus_p <- -cn_softplus(h)
  out <- rep(lchoose(N, y), length(h))
  if (y > 0) out <- out + y * log_p
  if (N - y > 0) out <- out + (N - y) * log_one_minus_p
  out
}

cn_binomial_residual <- function(y, N, h) {
  out <- numeric(length(h))
  nonnegative <- h >= 0
  out[nonnegative] <-
    (y - N) + N * plogis(-h[nonnegative])
  out[!nonnegative] <- y - N * plogis(h[!nonnegative])
  out
}

cn_gh_rule <- function(q) {
  if (length(q) != 1L || !is.finite(q) || q < 3 || q != as.integer(q)) {
    cn_stop("q must be one integer of at least three")
  }
  key <- as.character(as.integer(q))
  if (!exists(key, envir = .cn_gh_cache, inherits = FALSE)) {
    rule <- statmod::gauss.quad.prob(as.integer(q), dist = "normal")
    if (any(!is.finite(rule$nodes)) || any(!is.finite(rule$weights)) ||
        any(rule$weights <= 0) || abs(sum(rule$weights) - 1) > 1e-12) {
      cn_stop("invalid standard-normal Gauss-Hermite rule")
    }
    assign(key, rule, envir = .cn_gh_cache)
  }
  get(key, envir = .cn_gh_cache, inherits = FALSE)
}

cn_posterior_mode <- function(
    y, N, eta, sigma,
    tol = CN_D3_CONTRACT$integration$mode_tolerance,
    maxit = CN_D3_CONTRACT$integration$mode_max_iterations) {
  if (length(y) != 1L || length(N) != 1L || length(eta) != 1L ||
      length(sigma) != 1L || any(!is.finite(c(y, N, eta, sigma))) ||
      N <= 0 || y < 0 || y > N || sigma <= 0) {
    cn_stop("invalid posterior-mode inputs")
  }
  score <- function(z) {
    sigma * cn_binomial_residual(y, N, eta + sigma * z) - z
  }
  backward_error <- function(z) {
    h <- eta + sigma * z
    p <- plogis(h)
    abs(score(z)) / max(1, abs(sigma * y), abs(sigma * N * p), abs(z))
  }
  backward_tolerance <- 100 * .Machine$double.eps
  span <- 8
  while ((score(-span) <= 0 || score(span) >= 0) && span < 1024) span <- span * 2
  if (score(-span) <= 0 || score(span) >= 0) {
    cn_stop("could not bracket the unique latent posterior mode")
  }
  pseudo_p <- (y + 0.5) / (N + 1)
  z <- (qlogis(pseudo_p) - eta) / sigma
  z <- min(max(z, -span), span)
  lower <- -span
  upper <- span
  for (iter in seq_len(maxit)) {
    s <- score(z)
    if (!is.finite(s)) cn_stop("nonfinite posterior-mode score")
    if (backward_error(z) <= backward_tolerance) return(z)
    if (s > 0) lower <- z else upper <- z
    h <- eta + sigma * z
    p <- plogis(h)
    q <- plogis(-h)
    curvature <- 1 + sigma^2 * N * p * q
    proposal <- z + s / curvature
    if (!is.finite(proposal) || proposal <= lower || proposal >= upper) {
      proposal <- (lower + upper) / 2
    }
    if (abs(proposal - z) <= tol * (1 + abs(proposal))) break
    z <- proposal
  }
  root <- uniroot(
    score, interval = c(lower, upper), tol = .Machine$double.eps,
    maxiter = 2000L
  )
  if (!is.finite(root$root) || backward_error(root$root) > backward_tolerance) {
    cn_stop("posterior-mode solver did not meet its backward-error contract")
  }
  root$root
}

cn_integration_control <- function(level = c("production", "audit")) {
  level <- match.arg(level)
  if (level == "production") {
    relative_tolerance <-
      CN_D3_CONTRACT$integration$production_relative_tolerance
  } else {
    relative_tolerance <- CN_D3_CONTRACT$integration$audit_relative_tolerance
  }
  list(
    level = level,
    relative_tolerance = relative_tolerance,
    tail_radius = CN_D3_CONTRACT$integration$tail_radius,
    maximum_panels = CN_D3_CONTRACT$integration$maximum_panels
  )
}

cn_library_envelope_check <- function(N) {
  envelope <- CN_D3_CONTRACT$integration$certified_library_size
  N <- as.numeric(N)
  list(
    ok = length(N) > 0L && all(is.finite(N)) && all(N == round(N)) &&
      all(N >= envelope["lower"] & N <= envelope["upper"]),
    observed_range = if (length(N) && all(is.finite(N))) range(N) else
      c(NA_real_, NA_real_),
    envelope = envelope
  )
}

cn_eta_envelope_check <- function(eta) {
  envelope <- CN_D3_CONTRACT$integration$certified_eta_predictor
  eta <- as.numeric(eta)
  list(
    ok = length(eta) > 0L && all(is.finite(eta)) &&
      all(eta >= envelope["lower"] & eta <= envelope["upper"]),
    fitted_range = if (length(eta) && all(is.finite(eta))) range(eta) else
      c(NA_real_, NA_real_),
    envelope = envelope
  )
}

cn_eta_envelope_condition <- function(check) {
  structure(
    list(
      message = "present-location predictor is outside the certified envelope",
      call = NULL,
      fitted_eta_range = check$fitted_range,
      certified_eta_predictor = check$envelope
    ),
    class = c("cn_eta_envelope_error", "error", "condition")
  )
}

cn_eta_trial_violation_summary <- function(candidates) {
  violations <- unlist(lapply(candidates, function(candidate) {
    candidate$eta_trial_violations
  }), recursive = FALSE)
  ranges <- lapply(violations, function(violation) {
    as.numeric(violation$trial_eta_range)
  })
  finite_ranges <- ranges[vapply(ranges, function(x) {
    length(x) == 2L && all(is.finite(x))
  }, logical(1))]
  observed <- if (length(finite_ranges)) {
    range(unlist(finite_ranges, use.names = FALSE))
  } else c(NA_real_, NA_real_)
  list(
    count = length(violations),
    trial_eta_range = observed,
    certified_eta_predictor =
      CN_D3_CONTRACT$integration$certified_eta_predictor,
    violations = violations
  )
}

# The production likelihood uses a mode-centered vector Gauss--Kronrod rule.
# It integrates the denominator and the first two posterior moments together
# on a finite interval.  Strong log concavity supplies explicit omitted-tail
# bounds.  The Fisher scores then follow from the exact integration-by-parts
# identities E(y-Np | y,S=0)=E(Z | y,S=0)/sigma and
# E{sigma*Z*(y-Np) | y,S=0}=E(Z^2 | y,S=0)-1.  This avoids cancellation in
# signed score integrals and remains stable at original shotgun-scale depths.
cn_expm1_minus_x <- function(x) {
  x <- as.numeric(x)
  out <- expm1(x) - x
  small <- is.finite(x) & abs(x) <= 0.5
  out[small] <- vapply(x[small], function(value) {
    if (value == 0) return(0)
    term <- value^2 / 2
    total <- term
    for (order in 3:100) {
      term <- term * value / order
      updated <- total + term
      if (abs(term) <= .Machine$double.eps *
          max(abs(updated), .Machine$double.xmin)) {
        total <- updated
        break
      }
      total <- updated
    }
    total
  }, numeric(1))
  out
}

cn_log1p_minus_x <- function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x <= -1)) {
    cn_stop("log1p remainder received an input outside its domain")
  }
  out <- log1p(x) - x
  small <- abs(x) <= 0.5
  out[small] <- vapply(x[small], function(value) {
    if (value == 0) return(0)
    term <- -value^2 / 2
    total <- term
    for (order in 3:100) {
      term <- term * (-value) * (order - 1) / order
      updated <- total + term
      if (abs(term) <= .Machine$double.eps *
          max(abs(updated), .Machine$double.xmin)) {
        total <- updated
        break
      }
      total <- updated
    }
    total
  }, numeric(1))
  out
}

cn_softplus_bregman <- function(h_mode, displacement, p_mode) {
  displacement <- as.numeric(displacement)
  if (h_mode > 0) {
    base <- -h_mode
    stable_displacement <- -displacement
    derivative <- plogis(base)
  } else {
    base <- h_mode
    stable_displacement <- displacement
    derivative <- plogis(base)
  }
  out <- cn_softplus(base + stable_displacement) - cn_softplus(base) -
    derivative * stable_displacement
  small <- abs(stable_displacement) <= 0.5
  if (any(small)) {
    exponential_remainder <- cn_expm1_minus_x(stable_displacement[small])
    x <- derivative * expm1(stable_displacement[small])
    log_remainder <- cn_log1p_minus_x(x)
    out[small] <- log_remainder + derivative * exponential_remainder
  }
  rounding_limit <- 100 * .Machine$double.eps *
    pmax(1, abs(cn_softplus(base + stable_displacement)),
         abs(cn_softplus(base)), abs(derivative * stable_displacement))
  if (any(!is.finite(out)) || any(out < -rounding_limit)) {
    cn_stop("softplus Bregman remainder violated convexity")
  }
  pmax(out, 0)
}

cn_log_ratio_at_mode <- function(
    u, y, N, eta, sigma, mode, p_mode = NULL, mode_residual = NULL) {
  h_mode <- eta + sigma * mode
  if (is.null(p_mode)) p_mode <- plogis(h_mode)
  if (is.null(mode_residual)) {
    mode_residual <- sigma * (y - N * p_mode) - mode
  }
  displacement <- sigma * u
  divergence <- cn_softplus_bregman(h_mode, displacement, p_mode)
  upper_bound <- mode_residual * u - u^2 / 2
  out <- upper_bound - N * divergence
  if (any(!is.finite(out)) ||
      any(out > upper_bound + 1e-10 * pmax(1, abs(upper_bound)))) {
    cn_stop("anchored posterior ratio violated strong log concavity")
  }
  out
}

cn_anchor_at_mode <- function(y, N, eta, sigma, mode) {
  h_mode <- eta + sigma * mode
  stable_probability <- if (h_mode >= 0) plogis(-h_mode) else plogis(h_mode)
  stable_count <- if (h_mode >= 0) N - y else y
  log_binomial <- if (stable_probability > 0 && stable_probability < 1) {
    # Use binomial symmetry so the probability passed to dbinom is at most
    # one half; representing 1-q near one would lose q-sized anchor changes.
    dbinom(stable_count, size = N, prob = stable_probability, log = TRUE)
  } else {
    # At a rounded probability boundary, retain the finite softplus tail that
    # dbinom(..., prob=0/1) would discard.  The y=0/N guards avoid 0*(-Inf).
    cn_log_binom_kernel(mode, y, N, eta, sigma)
  }
  if (length(log_binomial) != 1L || !is.finite(log_binomial)) {
    cn_stop("nonfinite posterior-mode binomial anchor")
  }
  unname(log_binomial + dnorm(mode, log = TRUE))
}

cn_vector_integrand <- function(
    u, y, N, eta, sigma, mode, p_mode, mode_residual) {
  base <- exp(cn_log_ratio_at_mode(
    u, y, N, eta, sigma, mode, p_mode, mode_residual
  ))
  cbind(base, base * u, base * u^2)
}

cn_gk15_panel <- function(
    a, b, y, N, eta, sigma, mode, p_mode, mode_residual) {
  nodes <- c(
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
    -0.7415311855993945, -0.5860872354676911, -0.4058451513773972,
    -0.2077849550078985, 0,
     0.2077849550078985,  0.4058451513773972,  0.5860872354676911,
     0.7415311855993945,  0.8648644233597691,  0.9491079123427585,
     0.9914553711208126
  )
  kronrod_weights <- c(
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502,
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854,
    0.2044329400752989, 0.2094821410847278,
    0.2044329400752989, 0.1903505780647854, 0.1690047266392679,
    0.1406532597155259, 0.1047900103222502, 0.06309209262997855,
    0.02293532201052922
  )
  gauss_weights <- c(
    0.1294849661688697, 0.2797053914892767, 0.3818300505051189,
    0.4179591836734694, 0.3818300505051189, 0.2797053914892767,
    0.1294849661688697
  )
  midpoint <- (a + b) / 2
  halfwidth <- (b - a) / 2
  u <- midpoint + halfwidth * nodes
  values <- cn_vector_integrand(
    u, y, N, eta, sigma, mode, p_mode, mode_residual
  )
  kronrod <- halfwidth * colSums(values * kronrod_weights)
  gauss <- halfwidth * colSums(
    values[seq(2L, 14L, by = 2L), , drop = FALSE] * gauss_weights
  )
  list(a = a, b = b, value = kronrod, error = abs(kronrod - gauss))
}

cn_tail_moment_bounds <- function(radius, mode_residual = 0) {
  r <- as.numeric(mode_residual)
  exponential_anchor <- exp(r^2 / 2)
  right_probability <- sqrt(2 * pi) * pnorm(r - radius)
  left_probability <- sqrt(2 * pi) * pnorm(-radius - r)
  right_density <- exp(-(radius - r)^2 / 2)
  left_density <- exp(-(radius + r)^2 / 2)
  denominator <- exponential_anchor *
    (right_probability + left_probability)
  absolute_first <- exponential_anchor * (
    right_density + r * right_probability +
      left_density - r * left_probability
  )
  second <- exponential_anchor * (
    (radius + r) * right_density +
      (radius - r) * left_density +
      (1 + r^2) * (right_probability + left_probability)
  )
  c(denominator, absolute_first, second)
}

cn_initial_vector_breaks <- function(local_scale, radius) {
  positive <- local_scale
  while (tail(positive, 1L) < radius / 2) {
    positive <- c(positive, 2 * tail(positive, 1L))
  }
  positive <- unique(pmin(positive, radius))
  positive <- positive[positive > 0 & positive < radius]
  sort(unique(c(-radius, -rev(positive), 0, positive, radius)))
}

cn_integrated_present <- function(
    y, N, eta, sigma, level = c("production", "audit")) {
  library_envelope <- cn_library_envelope_check(N)
  eta_envelope <- cn_eta_envelope_check(eta)
  if (length(y) != 1L || !is.finite(y) || y < 0 || y != round(y) ||
      y > N || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
      !isTRUE(library_envelope$ok) || !isTRUE(eta_envelope$ok)) {
    cn_stop("integration inputs are outside the certified numerical envelope")
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
  breaks <- cn_initial_vector_breaks(local_scale, control$tail_radius)
  panels <- lapply(seq_len(length(breaks) - 1L), function(index) {
    cn_gk15_panel(
      breaks[index], breaks[index + 1L], y, N, eta, sigma, mode,
      p_mode, mode_residual
    )
  })
  tail_error <- cn_tail_moment_bounds(control$tail_radius, mode_residual)
  converged <- FALSE
  repeat {
    value <- unname(Reduce(`+`, lapply(panels, `[[`, "value")))
    error <- unname(
      Reduce(`+`, lapply(panels, `[[`, "error")) + tail_error
    )
    denominator <- value[1L]
    if (is.finite(denominator) && denominator > error[1L]) {
      mean_u_now <- value[2L] / denominator
      second_u_now <- value[3L] / denominator
      log_error_bound <- -log1p(-error[1L] / denominator)
      mean_u_error_bound <-
        (error[2L] + abs(mean_u_now) * error[1L]) /
        (denominator - error[1L])
      second_u_error_bound <-
        (error[3L] + abs(second_u_now) * error[1L]) /
        (denominator - error[1L])
      transformed_error <- c(
        log_K = log_error_bound,
        score_eta = mean_u_error_bound / sigma,
        score_log_sigma =
          2 * abs(mode) * mean_u_error_bound + second_u_error_bound
      )
      if (all(is.finite(transformed_error)) &&
          all(transformed_error <= control$relative_tolerance)) {
        converged <- TRUE
        break
      }
    } else {
      transformed_error <- rep(Inf, 3L)
    }
    if (length(panels) >= control$maximum_panels) break
    ratios <- vapply(panels, function(panel) {
      if (!is.finite(denominator) || denominator <= 0) {
        return(max(panel$error))
      }
      local_mean <-
        (panel$error[2L] + abs(value[2L] / denominator) * panel$error[1L]) /
        denominator
      local_second <-
        (panel$error[3L] + abs(value[3L] / denominator) * panel$error[1L]) /
        denominator
      max(
        panel$error[1L] / denominator,
        local_mean / sigma,
        2 * abs(mode) * local_mean + local_second
      )
    }, numeric(1))
    split_index <- which.max(ratios)
    old <- panels[[split_index]]
    midpoint <- (old$a + old$b) / 2
    children <- list(
      cn_gk15_panel(old$a, midpoint, y, N, eta, sigma, mode,
                    p_mode, mode_residual),
      cn_gk15_panel(midpoint, old$b, y, N, eta, sigma, mode,
                    p_mode, mode_residual)
    )
    panels <- append(panels[-split_index], children)
  }
  if (!converged) {
    cn_stop("vector Gauss--Kronrod integration did not meet its error contract")
  }
  denominator <- value[1L]
  mean_u <- value[2L] / denominator
  second_u <- value[3L] / denominator
  score_eta <- (mode + mean_u) / sigma
  score_log_sigma <- mode^2 + 2 * mode * mean_u + second_u - 1
  anchor <- cn_anchor_at_mode(y, N, eta, sigma, mode)
  if (any(!is.finite(c(
    denominator, mean_u, second_u, score_eta, score_log_sigma, anchor,
    transformed_error
  )))) cn_stop("nonfinite vector integration result")
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
      error, c("denominator", "first_centered_moment", "second_centered_moment")
    ),
    integration_transformed_error = transformed_error,
    integration_panels = length(panels),
    integration_evaluations = 15L * (length(breaks) - 1L) +
      30L * (length(panels) - (length(breaks) - 1L)),
    integration_tail_bound = tail_error
  )
}

cn_adaptive_present <- function(
    y, N, eta, sigma,
    q = CN_D3_CONTRACT$integration$diagnostic_agh_nodes,
    keep_nodes = FALSE) {
  mode <- cn_posterior_mode(y, N, eta, sigma)
  p_mode <- plogis(eta + sigma * mode)
  scale <- 1 / sqrt(1 + sigma^2 * N * p_mode * (1 - p_mode))
  if (!is.finite(scale) || scale <= 0) cn_stop("invalid adaptive quadrature scale")
  gh <- cn_gh_rule(q)
  zstar <- mode + scale * gh$nodes
  log_likelihood <- cn_log_binom_kernel(zstar, y, N, eta, sigma)
  log_target <- log_likelihood + dnorm(zstar, log = TRUE)
  log_proposal <- dnorm(zstar, mean = mode, sd = scale, log = TRUE)
  log_terms <- log(gh$weights) + log_target - log_proposal
  log_K <- cn_log_sum_exp(log_terms)
  posterior_weight <- exp(log_terms - log_K)
  if (!is.finite(log_K) || any(!is.finite(posterior_weight)) ||
      abs(sum(posterior_weight) - 1) > 1e-10) {
    cn_stop("nonfinite adaptive quadrature result")
  }
  p <- plogis(eta + sigma * zstar)
  d <- y - N * p
  score_eta <- sum(posterior_weight * d)
  score_log_sigma <- sum(posterior_weight * d * sigma * zstar)
  ans <- list(
    log_K = log_K,
    score_eta = score_eta,
    score_log_sigma = score_log_sigma,
    mode = mode,
    scale = scale
  )
  if (keep_nodes) {
    ans[c("zstar", "posterior_weight", "d")] <-
      list(zstar, posterior_weight, d)
  }
  ans
}

cn_reference_present <- function(
    y, N, eta, sigma, rel.tol = 2e-11, subdivisions = 3000L) {
  mode <- cn_posterior_mode(y, N, eta, sigma)
  p_mode <- plogis(eta + sigma * mode)
  scale <- 1 / sqrt(1 + sigma^2 * N * p_mode * (1 - p_mode))
  anchor <- cn_log_binom_kernel(mode, y, N, eta, sigma) + dnorm(mode, log = TRUE)
  base <- function(t) {
    z <- mode + scale * t
    exp(cn_log_binom_kernel(z, y, N, eta, sigma) +
          dnorm(z, log = TRUE) - anchor) * scale
  }
  integrate_one <- function(fun) {
    integrate(
      fun, lower = -Inf, upper = Inf, rel.tol = rel.tol,
      abs.tol = rel.tol / 100,
      subdivisions = subdivisions, stop.on.error = TRUE
    )$value
  }
  denominator <- integrate_one(base)
  if (!is.finite(denominator) || denominator <= 0) {
    cn_stop("reference integration failed")
  }
  numerator_eta <- integrate_one(function(t) {
    z <- mode + scale * t
    base(t) * (y - N * plogis(eta + sigma * z))
  })
  numerator_scale <- integrate_one(function(t) {
    z <- mode + scale * t
    d <- y - N * plogis(eta + sigma * z)
    base(t) * d * sigma * z
  })
  list(
    log_K = anchor + log(denominator),
    score_eta = numerator_eta / denominator,
    score_log_sigma = numerator_scale / denominator,
    mode = mode,
    scale = scale
  )
}

cn_observation_eval <- function(
    y, N, nu, eta, sigma,
    integration_level = c("production", "audit")) {
  integration_level <- match.arg(integration_level)
  present <- cn_integrated_present(
    y, N, eta, sigma, level = integration_level
  )
  log_rho <- -cn_softplus(-nu)
  log_one_minus_rho <- -cn_softplus(nu)
  rho <- plogis(nu)
  if (y == 0) {
    log_present_zero <- log_one_minus_rho + present$log_K
    log_f <- cn_logspace_add(log_rho, log_present_zero)
    gamma <- exp(log_rho - log_f)
    tau <- exp(log_present_zero - log_f)
    posterior_sum <- gamma + tau
    gamma <- gamma / posterior_sum
    tau <- tau / posterior_sum
  } else {
    log_f <- log_one_minus_rho + present$log_K
    gamma <- 0
    tau <- 1
  }
  r <- gamma - rho
  v <- tau * present$score_eta
  h <- tau * present$score_log_sigma
  if (any(!is.finite(c(log_f, rho, gamma, tau, r, v, h)))) {
    cn_stop("nonfinite observed-likelihood or score contribution")
  }
  list(
    loglik = log_f,
    r = r,
    v = v,
    h = h,
    rho = rho,
    gamma = gamma,
    tau = tau,
    log_K = present$log_K,
    mode = present$mode,
    adaptive_scale = present$scale
  )
}

cn_covariate_matrix <- function(z, n, prefix) {
  if (is.null(z)) return(NULL)
  Z <- as.matrix(z)
  storage.mode(Z) <- "double"
  if (nrow(Z) != n || !ncol(Z) || any(!is.finite(Z))) {
    cn_stop(prefix, " must be a finite matrix with one row per observation")
  }
  if (is.null(colnames(Z))) colnames(Z) <- paste0(prefix, seq_len(ncol(Z)))
  colnames(Z) <- make.unique(paste0(prefix, ":", colnames(Z)))
  Z
}

cn_prepare_data <- function(
    y, N, g, z = NULL, z_rho = z, z_eta = z, subject_id = NULL) {
  y <- as.numeric(y)
  N <- as.numeric(N)
  g <- as.numeric(g)
  n <- length(y)
  if (!n || length(N) != n || length(g) != n) {
    cn_stop("y, N and g must have the same positive length")
  }
  if (any(!is.finite(y)) || any(y < 0) || any(y != round(y))) {
    cn_stop("y must contain nonnegative integer counts")
  }
  if (any(!is.finite(N)) || any(N <= 0) || any(N != round(N))) {
    cn_stop("N must contain positive original integer library sizes")
  }
  if (any(y > N)) cn_stop("each count y must be no greater than N")
  if (any(!is.finite(g)) || any(!g %in% c(0, 1))) {
    cn_stop("g must be coded 0/1")
  }
  if (is.null(subject_id)) subject_id <- seq_len(n)
  if (length(subject_id) != n || anyNA(subject_id) || anyDuplicated(subject_id)) {
    cn_stop("subject_id must be nonmissing, unique, and have one value per row")
  }
  subject_id <- as.character(subject_id)
  Zrho <- cn_covariate_matrix(z_rho, n, "rho")
  Zeta <- cn_covariate_matrix(z_eta, n, "eta")
  Xrho <- matrix(1, nrow = n, ncol = 1L,
                 dimnames = list(NULL, "Intercept"))
  Xeta <- matrix(1, nrow = n, ncol = 1L,
                 dimnames = list(NULL, "Intercept"))
  if (!is.null(Zrho)) Xrho <- cbind(Xrho, Zrho)
  if (!is.null(Zeta)) Xeta <- cbind(Xeta, Zeta)
  structure(
    list(
      y = y, N = N, g = g, n = n, subject_id = subject_id,
      Xrho = Xrho, Xeta = Xeta
    ),
    class = "cn_d3_data"
  )
}

cn_layout <- function(data) {
  if (!inherits(data, "cn_d3_data")) cn_stop("data must come from cn_prepare_data")
  p_alpha <- ncol(data$Xrho)
  p_beta <- ncol(data$Xeta)
  natural_names <- c(
    "delta", "zeta",
    paste0("alpha:", colnames(data$Xrho)),
    paste0("beta:", colnames(data$Xeta)),
    "omega", "kappa"
  )
  idx <- list(
    delta = 1L,
    zeta = 2L,
    alpha = 2L + seq_len(p_alpha),
    beta = 2L + p_alpha + seq_len(p_beta),
    omega = 2L + p_alpha + p_beta + 1L,
    kappa = 2L + p_alpha + p_beta + 2L
  )
  list(
    p_alpha = p_alpha,
    p_beta = p_beta,
    natural_names = natural_names,
    idx = idx,
    p = length(natural_names),
    opt_names = c(
      "zeta",
      natural_names[idx$alpha], natural_names[idx$beta],
      "log_sigma_g0", "log_sigma_g1"
    )
  )
}

cn_named_theta <- function(theta, layout) {
  if (!is.null(names(theta))) {
    if (!setequal(names(theta), layout$natural_names)) {
      cn_stop("named theta does not match the count-native parameter layout")
    }
    theta <- theta[layout$natural_names]
  }
  theta <- as.numeric(theta)
  if (length(theta) != layout$p || any(!is.finite(theta))) {
    cn_stop("theta has the wrong length or contains a nonfinite value")
  }
  setNames(theta, layout$natural_names)
}

cn_opt_to_natural <- function(opt_par, layout) {
  if (is.null(names(opt_par)) || !setequal(names(opt_par), layout$opt_names)) {
    cn_stop("optimizer parameter names do not match the count-native layout")
  }
  opt_par <- opt_par[layout$opt_names]
  theta <- setNames(numeric(layout$p), layout$natural_names)
  theta["delta"] <- 0
  theta["zeta"] <- opt_par["zeta"]
  theta[layout$natural_names[layout$idx$alpha]] <-
    opt_par[layout$natural_names[layout$idx$alpha]]
  theta[layout$natural_names[layout$idx$beta]] <-
    opt_par[layout$natural_names[layout$idx$beta]]
  theta["omega"] <- opt_par["log_sigma_g0"]
  theta["kappa"] <- opt_par["log_sigma_g1"] - opt_par["log_sigma_g0"]
  theta
}

cn_natural_to_opt <- function(theta, layout) {
  theta <- cn_named_theta(theta, layout)
  out <- setNames(numeric(length(layout$opt_names)), layout$opt_names)
  out["zeta"] <- theta["zeta"]
  out[layout$natural_names[layout$idx$alpha]] <-
    theta[layout$natural_names[layout$idx$alpha]]
  out[layout$natural_names[layout$idx$beta]] <-
    theta[layout$natural_names[layout$idx$beta]]
  out["log_sigma_g0"] <- theta["omega"]
  out["log_sigma_g1"] <- theta["omega"] + theta["kappa"]
  out
}

cn_score_natural_to_opt <- function(total_score, layout) {
  if (!is.null(names(total_score))) {
    if (!setequal(names(total_score), layout$natural_names)) {
      cn_stop("named score does not match the count-native parameter layout")
    }
    total_score <- total_score[layout$natural_names]
  }
  total_score <- as.numeric(total_score)
  if (length(total_score) != layout$p || any(!is.finite(total_score))) {
    cn_stop("score has the wrong length or contains a nonfinite value")
  }
  total_score <- setNames(total_score, layout$natural_names)
  out <- setNames(numeric(length(layout$opt_names)), layout$opt_names)
  out["zeta"] <- total_score["zeta"]
  out[layout$natural_names[layout$idx$alpha]] <-
    total_score[layout$natural_names[layout$idx$alpha]]
  out[layout$natural_names[layout$idx$beta]] <-
    total_score[layout$natural_names[layout$idx$beta]]
  out["log_sigma_g0"] <- total_score["omega"] - total_score["kappa"]
  out["log_sigma_g1"] <- total_score["kappa"]
  out
}

cn_eval <- function(
    theta, data,
    integration_level = c("production", "audit")) {
  integration_level <- match.arg(integration_level)
  layout <- cn_layout(data)
  theta <- cn_named_theta(theta, layout)
  idx <- layout$idx
  nu <- as.numeric(
    data$Xrho %*% theta[idx$alpha] + theta["delta"] * data$g
  )
  eta <- as.numeric(
    data$Xeta %*% theta[idx$beta] + theta["zeta"] * data$g
  )
  log_sigma <- as.numeric(theta["omega"] + theta["kappa"] * data$g)
  sigma <- exp(log_sigma)
  if (any(!is.finite(c(nu, eta, sigma))) || any(sigma <= 0)) {
    cn_stop("nonfinite count-native linear predictor")
  }
  eta_envelope <- cn_eta_envelope_check(eta)
  if (!isTRUE(eta_envelope$ok)) {
    stop(cn_eta_envelope_condition(eta_envelope))
  }
  rows <- lapply(seq_len(data$n), function(i) {
    cn_observation_eval(
      data$y[i], data$N[i], nu[i], eta[i], sigma[i],
      integration_level = integration_level
    )
  })
  pull <- function(name) vapply(rows, `[[`, numeric(1), name)
  r <- pull("r")
  v <- pull("v")
  h <- pull("h")
  score <- matrix(
    0, nrow = data$n, ncol = layout$p,
    dimnames = list(NULL, layout$natural_names)
  )
  score[, "delta"] <- data$g * r
  score[, "zeta"] <- data$g * v
  score[, layout$natural_names[idx$alpha]] <- data$Xrho * r
  score[, layout$natural_names[idx$beta]] <- data$Xeta * v
  score[, "omega"] <- h
  score[, "kappa"] <- data$g * h
  loglik_i <- pull("loglik")
  if (any(!is.finite(score)) || any(!is.finite(loglik_i))) {
    cn_stop("nonfinite count-native likelihood evaluation")
  }
  list(
    loglik = sum(loglik_i),
    objective = -sum(loglik_i),
    loglik_i = loglik_i,
    score = score,
    total_score = colSums(score),
    theta = theta,
    layout = layout,
    nu = nu,
    eta = eta,
    log_sigma = log_sigma,
    sigma = sigma,
    rho = pull("rho"),
    gamma = pull("gamma"),
    tau = pull("tau"),
    log_K = pull("log_K"),
    posterior_mode = pull("mode"),
    adaptive_scale = pull("adaptive_scale"),
    integration_level = integration_level
  )
}

cn_integration_audit <- function(theta, data) {
  production <- cn_eval(theta, data, integration_level = "production")
  audit <- cn_eval(theta, data, integration_level = "audit")
  list(
    production_control = cn_integration_control("production"),
    audit_control = cn_integration_control("audit"),
    loglik_difference = production$loglik - audit$loglik,
    maximum_subject_loglik_difference =
      max(abs(production$loglik_i - audit$loglik_i)),
    maximum_subject_score_difference =
      max(abs(production$score - audit$score)),
    maximum_total_score_difference =
      max(abs(production$total_score - audit$total_score)),
    production = production,
    audit = audit
  )
}

cn_rank <- function(X, tol = CN_D3_CONTRACT$precheck$rank_tolerance) {
  qr(X, tol = tol)$rank
}

cn_precheck <- function(data) {
  if (!inherits(data, "cn_d3_data")) {
    return(cn_failure("invalid_data_object", stage = "precheck"))
  }
  failures <- character()
  add <- function(code) failures <<- unique(c(failures, code))
  library_envelope <- cn_library_envelope_check(data$N)
  if (!isTRUE(library_envelope$ok)) {
    add("library_size_outside_certified_envelope")
  }
  if (length(unique(data$g)) != 2L) add("both_groups_required")
  if (cn_rank(data$Xrho) < ncol(data$Xrho)) add("structural_design_rank")
  if (cn_rank(data$Xeta) < ncol(data$Xeta)) add("present_design_rank")
  if (cn_rank(cbind(data$Xrho, GroupTarget = data$g)) <= cn_rank(data$Xrho)) {
    add("occupancy_group_not_identified")
  }
  if (cn_rank(cbind(data$Xeta, GroupLocation = data$g)) <= cn_rank(data$Xeta)) {
    add("location_group_not_identified")
  }
  positive <- data$y > 0
  present_design <- cbind(GroupLocation = data$g, data$Xeta)
  if (sum(positive) <= ncol(present_design) ||
      cn_rank(present_design[positive, , drop = FALSE]) < ncol(present_design)) {
    add("positive_present_design_rank")
  }
  for (gg in c(0, 1)) {
    if (sum(!positive & data$g == gg) <
        CN_D3_CONTRACT$precheck$minimum_zero_per_group) {
      add(paste0("insufficient_zero_group", gg))
    }
    if (sum(positive & data$g == gg) <
        CN_D3_CONTRACT$precheck$minimum_positive_per_group) {
      add(paste0("insufficient_positive_group", gg))
    }
  }
  if (length(failures)) {
    return(cn_failure(
      paste(failures, collapse = "+"),
      failure_codes = failures,
      observed_library_size_range = library_envelope$observed_range,
      certified_library_size = library_envelope$envelope,
      stage = "precheck"
    ))
  }
  list(available = TRUE, failure_code = "", failure_codes = character())
}

cn_opt_bounds <- function(layout) {
  lower <- upper <- setNames(numeric(length(layout$opt_names)), layout$opt_names)
  for (name in layout$opt_names) {
    if (grepl("^alpha:", name)) {
      lower[name] <- CN_D3_CONTRACT$parameter_bounds$structural["lower"]
      upper[name] <- CN_D3_CONTRACT$parameter_bounds$structural["upper"]
    } else if (name == "zeta" || grepl("^beta:", name)) {
      lower[name] <- CN_D3_CONTRACT$parameter_bounds$location["lower"]
      upper[name] <- CN_D3_CONTRACT$parameter_bounds$location["upper"]
    } else if (name %in% c("log_sigma_g0", "log_sigma_g1")) {
      lower[name] <- CN_D3_CONTRACT$parameter_bounds$log_scale["lower"]
      upper[name] <- CN_D3_CONTRACT$parameter_bounds$log_scale["upper"]
    } else {
      cn_stop("unclassified optimizer parameter: ", name)
    }
  }
  list(lower = lower, upper = upper)
}

cn_optimizer_eta_map <- function(data, layout) {
  map <- matrix(
    0, nrow = data$n, ncol = length(layout$opt_names),
    dimnames = list(NULL, layout$opt_names)
  )
  map[, "zeta"] <- data$g
  beta_names <- layout$natural_names[layout$idx$beta]
  map[, beta_names] <- data$Xeta
  map
}

cn_optimizer_constraints <- function(data, layout, bounds = cn_opt_bounds(layout)) {
  # These linear inequalities are a numerical-domain guard, not a constrained
  # estimand.  A formed fit must later be separated from every constraint and
  # satisfy the original (barrier-free) nuisance score equations.
  p <- length(layout$opt_names)
  identity <- diag(p)
  dimnames(identity) <- list(NULL, layout$opt_names)
  eta_map <- cn_optimizer_eta_map(data, layout)
  envelope <- CN_D3_CONTRACT$integration$certified_eta_predictor
  ui <- rbind(identity, -identity, eta_map, -eta_map)
  rownames(ui) <- c(
    paste0("coefficient_lower:", layout$opt_names),
    paste0("coefficient_upper:", layout$opt_names),
    paste0("eta_lower:", seq_len(data$n)),
    paste0("eta_upper:", seq_len(data$n))
  )
  ci <- c(
    unname(bounds$lower), -unname(bounds$upper),
    rep(unname(envelope["lower"]), data$n),
    rep(-unname(envelope["upper"]), data$n)
  )
  names(ci) <- rownames(ui)
  list(
    ui = ui, ci = ci, eta_map = eta_map, bounds = bounds,
    envelope = envelope,
    coefficient_rows = seq_len(2L * p),
    eta_lower_rows = 2L * p + seq_len(data$n),
    eta_upper_rows = 2L * p + data$n + seq_len(data$n)
  )
}

cn_constraint_slack <- function(opt_par, constraints) {
  opt_names <- colnames(constraints$ui)
  if (is.null(names(opt_par)) || !setequal(names(opt_par), opt_names)) {
    cn_stop("optimizer parameter names do not match the constraint layout")
  }
  opt_par <- as.numeric(opt_par[opt_names])
  setNames(
    drop(constraints$ui %*% opt_par - constraints$ci),
    rownames(constraints$ui)
  )
}

cn_strictly_feasible_start <- function(start, anchor, constraints) {
  # The feasible set is convex.  Contracting along the line from a strictly
  # feasible anchor preserves the linear predictor model; it never clips an
  # individual fitted eta or changes the likelihood.
  relative_margin <-
    CN_D3_CONTRACT$optimizer$strict_feasibility_relative_margin
  margin <- relative_margin * (1 + abs(constraints$ci))
  names(margin) <- rownames(constraints$ui)
  anchor_slack <- cn_constraint_slack(anchor, constraints)
  if (any(!is.finite(anchor_slack)) || any(anchor_slack <= margin)) {
    cn_stop("the count-native optimizer anchor is not strictly feasible")
  }
  start_slack <- cn_constraint_slack(start, constraints)
  if (all(is.finite(start_slack)) && all(start_slack > margin)) return(start)

  opt_names <- colnames(constraints$ui)
  start <- as.numeric(start[opt_names])
  anchor <- as.numeric(anchor[opt_names])
  slope <- drop(constraints$ui %*% (start - anchor))
  decreasing <- slope < 0
  maximum_fraction <- 1
  if (any(decreasing)) {
    maximum_fraction <- min(
      (anchor_slack[decreasing] - margin[decreasing]) / (-slope[decreasing])
    )
  }
  if (!is.finite(maximum_fraction) || maximum_fraction <= 0) {
    cn_stop("could not contract an optimizer start into the feasible interior")
  }
  fraction <- min(1, 0.95 * maximum_fraction)
  candidate <- setNames(anchor + fraction * (start - anchor), opt_names)
  candidate_slack <- cn_constraint_slack(candidate, constraints)
  if (any(!is.finite(candidate_slack)) || any(candidate_slack <= margin)) {
    cn_stop("contracted optimizer start is not strictly feasible")
  }
  candidate
}

cn_eta_boundary_status <- function(eta, subject_id) {
  envelope <- CN_D3_CONTRACT$integration$certified_eta_predictor
  relative <- CN_D3_CONTRACT$optimizer$boundary_relative_distance
  lower_slack <- eta - envelope["lower"]
  upper_slack <- envelope["upper"] - eta
  lower_limit <- relative * (1 + abs(envelope["lower"]))
  upper_limit <- relative * (1 + abs(envelope["upper"]))
  active_lower <- which(lower_slack <= lower_limit)
  active_upper <- which(upper_slack <= upper_limit)
  list(
    active = length(active_lower) > 0L || length(active_upper) > 0L,
    active_lower = active_lower,
    active_upper = active_upper,
    active_subject_id = unique(as.character(subject_id[c(
      active_lower, active_upper
    )])),
    minimum_lower_slack = min(lower_slack),
    minimum_upper_slack = min(upper_slack),
    fitted_eta_range = range(eta),
    certified_eta_predictor = envelope,
    boundary_limit = c(lower = lower_limit, upper = upper_limit)
  )
}

cn_start_set <- function(data, layout) {
  positive <- data$y > 0
  empirical_logit <- qlogis((data$y[positive] + 0.5) / (data$N[positive] + 1))
  design <- cbind(GroupLocation = data$g, data$Xeta)
  fit <- lm.fit(design[positive, , drop = FALSE], empirical_logit)
  coefficient <- fit$coefficients
  coefficient[!is.finite(coefficient)] <- 0
  residual_scale <- sd(fit$residuals)
  if (!is.finite(residual_scale) || residual_scale < 0.1) residual_scale <- 0.5
  residual_scale <- min(max(residual_scale, 0.2), 2)
  zero_fraction <- mean(!positive)
  envelope <- CN_D3_CONTRACT$integration$certified_eta_predictor
  anchor_eta <- median(empirical_logit)
  if (!is.finite(anchor_eta)) anchor_eta <- mean(envelope)
  anchor_eta <- min(max(
    anchor_eta, envelope["lower"] + 1
  ), envelope["upper"] - 1)
  start_spec <- matrix(c(
    0.25, 0.5, 0.00,  0.00, 0.00,
    0.50, 0.5, 0.25,  0.25, 0.00,
    0.50, 1.0, 0.50, -0.25, 0.00,
    0.25, 2.0, 1.00,  0.00, 0.25,
    0.75, 0.5, 0.50,  0.50, 0.00,
    0.75, 2.0, 0.50, -0.50, 0.50
  ), ncol = 5L, byrow = TRUE,
  dimnames = list(NULL, c(
    "structural_fraction", "scale_multiplier", "location_weight",
    "zeta_offset", "intercept_offset"
  )))
  bounds <- cn_opt_bounds(layout)
  constraints <- cn_optimizer_constraints(data, layout, bounds)
  beta_names <- layout$natural_names[layout$idx$beta]
  lapply(seq_len(nrow(start_spec)), function(k) {
    theta <- setNames(numeric(layout$p), layout$natural_names)
    theta["delta"] <- 0
    weight <- start_spec[k, "location_weight"]
    anchor_beta <- setNames(numeric(length(beta_names)), beta_names)
    anchor_beta[1L] <- anchor_eta
    theta["zeta"] <- weight * coefficient[1L] +
      start_spec[k, "zeta_offset"]
    theta[beta_names] <- anchor_beta + weight *
      (setNames(coefficient[-1L], beta_names) - anchor_beta)
    theta[beta_names[1L]] <- theta[beta_names[1L]] +
      start_spec[k, "intercept_offset"]
    structural_probability <- start_spec[k, "structural_fraction"] * zero_fraction
    structural_probability <- min(max(structural_probability, 0.005), 0.95)
    theta[layout$natural_names[layout$idx$alpha[1L]]] <-
      qlogis(structural_probability)
    sigma_start <- min(max(
      residual_scale * start_spec[k, "scale_multiplier"],
      exp(CN_D3_CONTRACT$parameter_bounds$log_scale["lower"] + 0.1)
    ), exp(CN_D3_CONTRACT$parameter_bounds$log_scale["upper"] - 0.1))
    theta["omega"] <- log(sigma_start)
    theta["kappa"] <- 0
    opt <- cn_natural_to_opt(theta, layout)
    opt <- setNames(
      pmin(pmax(opt, bounds$lower + 0.01), bounds$upper - 0.01),
      layout$opt_names
    )
    anchor <- opt
    anchor["zeta"] <- 0
    anchor[beta_names] <- anchor_beta
    cn_strictly_feasible_start(opt, anchor, constraints)
  })
}

cn_relative_stationarity <- function(score, parameter_names = colnames(score)) {
  score <- as.matrix(score)
  if (is.null(colnames(score)) || !length(parameter_names) ||
      !all(parameter_names %in% colnames(score))) {
    cn_stop("stationarity parameters do not match score columns")
  }
  selected <- score[, parameter_names, drop = FALSE]
  score_norm <- sqrt(colSums(selected^2))
  if (any(!is.finite(selected)) || any(!is.finite(score_norm)) ||
      any(score_norm <= 0)) {
    return(list(
      valid = FALSE, maximum = Inf,
      components = setNames(rep(Inf, length(parameter_names)), parameter_names),
      score_norm = score_norm
    ))
  }
  components <- abs(colSums(selected)) / score_norm
  list(
    valid = all(is.finite(components)),
    maximum = max(components),
    components = components,
    score_norm = score_norm
  )
}

cn_better_nonconverged_indices <- function(
    candidates, finite, converged, best_objective) {
  unresolved_index <- which(finite & !converged)
  if (!length(unresolved_index)) return(integer())
  unresolved_objective <- vapply(
    candidates[unresolved_index], `[[`, numeric(1), "objective"
  )
  tolerance <-
    CN_D3_CONTRACT$optimizer$near_mode_objective_relative_tolerance *
    (1 + abs(best_objective))
  unresolved_index[unresolved_objective < best_objective - tolerance]
}

cn_active_boundary_parameters <- function(opt_par, bounds) {
  distance_lower <- opt_par - bounds$lower
  distance_upper <- bounds$upper - opt_par
  boundary_limit <- CN_D3_CONTRACT$optimizer$boundary_relative_distance *
    (1 + pmax(abs(bounds$lower), abs(bounds$upper)))
  names(opt_par)[distance_lower <= boundary_limit |
                   distance_upper <= boundary_limit]
}

cn_near_mode_comparison <- function(selected, other, reference_eval, other_eval) {
  parameter_difference <- max(abs(other$theta - selected$theta))
  fitted_difference <- max(abs(c(
    other_eval$rho - reference_eval$rho,
    other_eval$eta - reference_eval$eta,
    other_eval$log_sigma - reference_eval$log_sigma,
    other_eval$gamma - reference_eval$gamma
  )))
  list(
    distinct =
      parameter_difference >
        CN_D3_CONTRACT$optimizer$distinct_mode_parameter_tolerance &&
      fitted_difference >
        CN_D3_CONTRACT$optimizer$distinct_mode_fitted_tolerance,
    parameter_difference = parameter_difference,
    fitted_difference = fitted_difference
  )
}

cn_interior_polish_controls <- function() {
  controls <- CN_D3_CONTRACT$optimizer$interior_polish
  expected_names <- c(
    "direction", "score_scale_floor", "fraction_to_boundary",
    "maximum_iterations", "maximum_backtracks", "armijo"
  )
  if (!is.list(controls) || !identical(names(controls), expected_names) ||
      !identical(controls$direction, "scaled_score_only") ||
      length(controls$score_scale_floor) != 1L ||
      !is.finite(controls$score_scale_floor) ||
      controls$score_scale_floor <= 0 ||
      length(controls$fraction_to_boundary) != 1L ||
      !is.finite(controls$fraction_to_boundary) ||
      controls$fraction_to_boundary <= 0 ||
      controls$fraction_to_boundary >= 1 ||
      length(controls$maximum_iterations) != 1L ||
      !is.finite(controls$maximum_iterations) ||
      controls$maximum_iterations < 1 ||
      controls$maximum_iterations !=
        as.integer(controls$maximum_iterations) ||
      length(controls$maximum_backtracks) != 1L ||
      !is.finite(controls$maximum_backtracks) ||
      controls$maximum_backtracks < 0 ||
      controls$maximum_backtracks !=
        as.integer(controls$maximum_backtracks) ||
      length(controls$armijo) != 1L || !is.finite(controls$armijo) ||
      controls$armijo <= 0 || controls$armijo >= 1) {
    cn_stop("invalid count-native interior-polish contract")
  }
  controls$maximum_iterations <- as.integer(controls$maximum_iterations)
  controls$maximum_backtracks <- as.integer(controls$maximum_backtracks)
  controls
}

cn_interior_polish_constraint_margin <- function(constraints) {
  if (!is.list(constraints) || !is.matrix(constraints$ui) ||
      !is.matrix(constraints$eta_map) || is.null(rownames(constraints$ui))) {
    cn_stop("unexpected interior-polish constraint layout")
  }
  p <- ncol(constraints$ui)
  if (nrow(constraints$ui) != 2L * p + 2L * nrow(constraints$eta_map)) {
    cn_stop("unexpected interior-polish constraint dimensions")
  }
  relative <- CN_D3_CONTRACT$optimizer$boundary_relative_distance
  coefficient_limit <- relative *
    (1 + pmax(abs(constraints$bounds$lower),
              abs(constraints$bounds$upper)))
  eta_envelope <- CN_D3_CONTRACT$integration$certified_eta_predictor
  eta_lower_limit <- as.numeric(
    relative * (1 + abs(eta_envelope["lower"]))
  )
  eta_upper_limit <- as.numeric(
    relative * (1 + abs(eta_envelope["upper"]))
  )
  margin <- c(
    coefficient_limit,
    coefficient_limit,
    rep(eta_lower_limit, nrow(constraints$eta_map)),
    rep(eta_upper_limit, nrow(constraints$eta_map))
  )
  names(margin) <- rownames(constraints$ui)
  if (length(margin) != nrow(constraints$ui) ||
      any(!is.finite(margin)) || any(margin <= 0)) {
    cn_stop("invalid accepted-interior constraint margin")
  }
  margin
}

cn_interior_polish_state <- function(
    opt_par, data, layout, constraints) {
  if (is.null(names(opt_par)) ||
      !setequal(names(opt_par), layout$opt_names)) {
    return(list(ok = FALSE, reason = "invalid_optimizer_coordinates"))
  }
  opt_par <- setNames(
    as.numeric(opt_par[layout$opt_names]), layout$opt_names
  )
  if (any(!is.finite(opt_par))) {
    return(list(
      ok = FALSE, reason = "invalid_optimizer_coordinates",
      opt_par = opt_par
    ))
  }
  slack <- tryCatch(
    cn_constraint_slack(opt_par, constraints), error = function(e) e
  )
  if (inherits(slack, "error")) {
    return(list(
      ok = FALSE, reason = "constraint_evaluation_failed",
      error = conditionMessage(slack), opt_par = opt_par
    ))
  }
  margin <- cn_interior_polish_constraint_margin(constraints)
  if (any(!is.finite(slack)) || any(slack <= margin)) {
    return(list(
      ok = FALSE, reason = "outside_accepted_interior",
      opt_par = opt_par, slack = slack, margin = margin
    ))
  }
  theta <- tryCatch(
    cn_opt_to_natural(opt_par, layout), error = function(e) e
  )
  if (inherits(theta, "error")) {
    return(list(
      ok = FALSE, reason = "coordinate_mapping_failed",
      error = conditionMessage(theta), opt_par = opt_par,
      slack = slack, margin = margin
    ))
  }
  evaluation <- tryCatch(
    cn_eval(theta, data, integration_level = "production"),
    error = function(e) e
  )
  if (inherits(evaluation, "error")) {
    return(list(
      ok = FALSE, reason = "evaluation_failed",
      error = conditionMessage(evaluation), opt_par = opt_par,
      theta = theta, slack = slack, margin = margin
    ))
  }
  nuisance <- setdiff(layout$natural_names, "delta")
  stationarity <- cn_relative_stationarity(evaluation$score, nuisance)
  list(
    ok = isTRUE(stationarity$valid),
    reason = if (isTRUE(stationarity$valid)) "" else
      "invalid_nuisance_score_norm",
    opt_par = opt_par,
    theta = theta,
    evaluation = evaluation,
    nuisance = nuisance,
    stationarity = stationarity,
    slack = slack,
    margin = margin
  )
}

cn_interior_polish_step_limit <- function(
    state, direction_opt, constraints, controls) {
  if (length(direction_opt) != ncol(constraints$ui) ||
      any(!is.finite(direction_opt))) {
    return(NA_real_)
  }
  slope <- drop(constraints$ui %*% as.numeric(direction_opt))
  decreasing <- slope < 0
  if (!any(decreasing)) return(1)
  available <- state$slack[decreasing] - state$margin[decreasing]
  maximum <- min(available / (-slope[decreasing]))
  if (!is.finite(maximum) || maximum <= 0) return(NA_real_)
  min(1, controls$fraction_to_boundary * maximum)
}

cn_interior_polish_direction <- function(state, layout, controls) {
  score <- state$evaluation$score[, state$nuisance, drop = FALSE]
  total <- colSums(score)
  scale <- pmax(colSums(score^2), controls$score_scale_floor)
  direction_natural <- total / scale
  full <- setNames(numeric(layout$p), layout$natural_names)
  full[state$nuisance] <- direction_natural
  proposed <- state$theta + full
  proposed["delta"] <- 0
  direction_opt <- tryCatch(
    cn_natural_to_opt(proposed, layout) - state$opt_par,
    error = function(e) rep(NA_real_, length(layout$opt_names))
  )
  names(direction_opt) <- layout$opt_names
  list(
    type = "scaled_score",
    natural = direction_natural,
    opt = direction_opt,
    ascent = sum(total * direction_natural)
  )
}

cn_interior_polish <- function(
    opt_par, data, layout, constraints) {
  controls <- cn_interior_polish_controls()
  state <- cn_interior_polish_state(
    opt_par, data, layout, constraints
  )
  if (!isTRUE(state$ok)) {
    return(list(
      ok = FALSE, reason = state$reason, initial = state,
      final = state, trace = list(),
      direction_contract = controls$direction
    ))
  }
  threshold <- CN_D3_CONTRACT$optimizer$stationarity_relative_score_norm
  initial <- state
  trace <- list()
  if (state$stationarity$maximum <= threshold) {
    return(list(
      ok = TRUE, reason = "already_stationary", initial = initial,
      final = state, trace = trace, iterations = 0L,
      direction_contract = controls$direction
    ))
  }
  for (iteration in seq_len(controls$maximum_iterations)) {
    direction <- cn_interior_polish_direction(state, layout, controls)
    if (!is.finite(direction$ascent) || direction$ascent <= 0 ||
        any(!is.finite(direction$opt))) {
      return(list(
        ok = FALSE, reason = "no_ascent_direction",
        initial = initial, final = state, trace = trace,
        direction_contract = controls$direction
      ))
    }
    step <- cn_interior_polish_step_limit(
      state, direction$opt, constraints, controls
    )
    if (!is.finite(step) || step <= 0) {
      return(list(
        ok = FALSE, reason = "no_feasible_step",
        initial = initial, final = state, trace = trace,
        direction_contract = controls$direction
      ))
    }
    accepted <- FALSE
    trial <- NULL
    evaluated_step <- NA_real_
    backtrack_used <- NA_integer_
    for (backtrack in seq.int(0L, controls$maximum_backtracks)) {
      evaluated_step <- step
      backtrack_used <- as.integer(backtrack)
      candidate <- state$opt_par + evaluated_step * direction$opt
      trial <- cn_interior_polish_state(
        candidate, data, layout, constraints
      )
      target_loglik <- state$evaluation$loglik +
        controls$armijo * evaluated_step * direction$ascent
      if (isTRUE(trial$ok) &&
          trial$evaluation$loglik >= target_loglik) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    trace[[length(trace) + 1L]] <- data.frame(
      iteration = iteration,
      direction = direction$type,
      accepted = accepted,
      step = evaluated_step,
      backtracks = backtrack_used,
      objective_before = state$evaluation$objective,
      objective_after = if (isTRUE(trial$ok))
        trial$evaluation$objective else NA_real_,
      stationarity_before = state$stationarity$maximum,
      stationarity_after = if (isTRUE(trial$ok))
        trial$stationarity$maximum else NA_real_,
      minimum_slack_after_margin = if (isTRUE(trial$ok))
        min(trial$slack - trial$margin) else NA_real_,
      stringsAsFactors = FALSE
    )
    if (!accepted) {
      return(list(
        ok = FALSE, reason = "line_search_failed",
        initial = initial, final = state, trace = trace,
        direction_contract = controls$direction
      ))
    }
    state <- trial
    if (state$stationarity$maximum <= threshold) {
      return(list(
        ok = TRUE, reason = "stationary", initial = initial,
        final = state, trace = trace, iterations = as.integer(iteration),
        direction_contract = controls$direction
      ))
    }
  }
  list(
    ok = FALSE, reason = "iteration_limit", initial = initial,
    final = state, trace = trace,
    iterations = controls$maximum_iterations,
    direction_contract = controls$direction
  )
}

cn_interior_polish_audit_ok <- function(audit) {
  if (is.null(audit) || inherits(audit, "error")) return(FALSE)
  values <- c(
    audit$maximum_subject_loglik_difference,
    audit$maximum_subject_score_difference,
    audit$maximum_total_score_difference
  )
  limits <- c(
    CN_D3_CONTRACT$integration$maximum_subject_loglik_difference,
    CN_D3_CONTRACT$integration$maximum_subject_score_difference,
    CN_D3_CONTRACT$integration$maximum_total_score_difference
  )
  length(values) == 3L && length(limits) == 3L &&
    isTRUE(all(is.finite(values)) && all(is.finite(limits)) &&
             all(values <= limits))
}

cn_interior_polish_candidate_status <- function(
    candidate, data, layout, bounds, constraints) {
  finite <- is.list(candidate) && is.null(candidate$error) &&
    length(candidate$objective) == 1L &&
    is.finite(candidate$objective) &&
    length(candidate$theta) == layout$p &&
    all(is.finite(candidate$theta)) &&
    length(candidate$opt_par) == length(layout$opt_names) &&
    all(is.finite(candidate$opt_par))
  if (!finite || !isTRUE(candidate$convergence == 0L)) {
    return(list(kind = "not_converged", candidate = candidate))
  }
  active <- cn_active_boundary_parameters(candidate$opt_par, bounds)
  if (length(active)) {
    return(list(
      kind = "noninterior", reason = "active_coefficient_boundary",
      active = active, candidate = candidate,
      objective = candidate$objective
    ))
  }
  state <- cn_interior_polish_state(
    candidate$opt_par, data, layout, constraints
  )
  if (!isTRUE(state$ok)) {
    kind <- if (identical(state$reason, "outside_accepted_interior"))
      "noninterior" else "unresolved_interior"
    return(list(
      kind = kind, reason = state$reason, candidate = candidate,
      objective = candidate$objective, state = state
    ))
  }
  eta_boundary <- cn_eta_boundary_status(
    state$evaluation$eta, data$subject_id
  )
  if (isTRUE(eta_boundary$active)) {
    return(list(
      kind = "noninterior",
      reason = "active_present_location_predictor_envelope",
      candidate = candidate, objective = candidate$objective,
      state = state, eta_boundary = eta_boundary
    ))
  }
  list(
    kind = "interior", candidate = candidate, state = state,
    objective = candidate$objective, eta_boundary = eta_boundary
  )
}

cn_interior_polish_candidate <- function(
    status, data, layout, bounds, constraints) {
  if (!identical(status$kind, "interior")) return(status)
  polish <- cn_interior_polish(
    status$candidate$opt_par, data, layout, constraints
  )
  if (!isTRUE(polish$ok)) {
    return(list(
      kind = "unresolved_interior", reason = polish$reason,
      candidate = status$candidate, polish = polish
    ))
  }
  final <- polish$final
  active <- cn_active_boundary_parameters(final$opt_par, bounds)
  eta_boundary <- cn_eta_boundary_status(
    final$evaluation$eta, data$subject_id
  )
  slack <- cn_constraint_slack(final$opt_par, constraints)
  audit <- tryCatch(
    cn_integration_audit(final$theta, data), error = function(e) e
  )
  threshold <- CN_D3_CONTRACT$optimizer$stationarity_relative_score_norm
  valid <- !length(active) && !isTRUE(eta_boundary$active) &&
    length(slack) == nrow(constraints$ui) && length(slack) > 0L &&
    all(is.finite(slack)) && all(slack > 0) &&
    cn_interior_polish_audit_ok(audit) &&
    is.finite(final$stationarity$maximum) &&
    final$stationarity$maximum <= threshold &&
    final$evaluation$objective <=
      polish$initial$evaluation$objective +
        .Machine$double.eps *
          (1 + abs(polish$initial$evaluation$objective))
  if (!valid) {
    return(list(
      kind = "unresolved_interior",
      reason = "post_polish_contract_failure",
      candidate = status$candidate, polish = polish,
      active = active, eta_boundary = eta_boundary,
      slack = slack, audit = audit
    ))
  }
  candidate <- status$candidate
  candidate$opt_par <- final$opt_par
  candidate$theta <- final$theta
  candidate$objective <- final$evaluation$objective
  candidate$constraint_slack <- slack
  candidate$barrier_value <- NA_real_
  candidate$message <- "barrier-free scaled-score interior polish"
  candidate$interior_polish <- polish
  list(
    kind = "resolved_interior", candidate = candidate,
    state = final, evaluation = final$evaluation,
    eta_boundary = eta_boundary, slack = slack, audit = audit,
    polish = polish, objective = final$evaluation$objective
  )
}

cn_repair_nonstationary_restricted_fit <- function(baseline, data) {
  if (isTRUE(baseline$available) ||
      !identical(baseline$failure_code, "nonstationary_restricted_fit")) {
    return(baseline)
  }
  if (!is.list(baseline$candidates) || !length(baseline$candidates)) {
    return(cn_failure(
      "missing_restricted_candidates_for_interior_polish",
      baseline_fit = baseline, stage = "restricted_fit"
    ))
  }
  layout <- cn_layout(data)
  bounds <- cn_opt_bounds(layout)
  constraints <- cn_optimizer_constraints(data, layout, bounds)
  statuses <- lapply(
    baseline$candidates,
    cn_interior_polish_candidate_status,
    data = data, layout = layout, bounds = bounds,
    constraints = constraints
  )
  polished <- lapply(
    statuses,
    cn_interior_polish_candidate,
    data = data, layout = layout, bounds = bounds,
    constraints = constraints
  )
  unresolved <- which(vapply(
    polished,
    function(x) identical(x$kind, "unresolved_interior"),
    logical(1)
  ))
  if (length(unresolved)) {
    return(cn_failure(
      "unresolved_restricted_candidate_polish",
      baseline_fit = baseline,
      unresolved_indices = unresolved,
      unresolved_reasons = vapply(
        polished[unresolved], `[[`, character(1), "reason"
      ),
      candidate_repairs = polished,
      stage = "restricted_fit"
    ))
  }
  resolved <- which(vapply(
    polished,
    function(x) identical(x$kind, "resolved_interior"),
    logical(1)
  ))
  if (!length(resolved)) {
    return(cn_failure(
      "no_resolved_interior_restricted_candidate",
      baseline_fit = baseline, candidate_repairs = polished,
      stage = "restricted_fit"
    ))
  }
  objectives <- vapply(
    polished[resolved], `[[`, numeric(1), "objective"
  )
  selected_index <- resolved[which.min(objectives)]
  selected <- polished[[selected_index]]
  best <- selected$objective
  tolerance <-
    CN_D3_CONTRACT$optimizer$near_mode_objective_relative_tolerance *
    (1 + abs(best))

  finite <- vapply(baseline$candidates, function(candidate) {
    is.null(candidate$error) && is.finite(candidate$objective) &&
      all(is.finite(candidate$theta))
  }, logical(1))
  converged <- finite & vapply(baseline$candidates, function(candidate) {
    isTRUE(candidate$convergence == 0L)
  }, logical(1))

  noninterior <- which(vapply(
    polished,
    function(x) identical(x$kind, "noninterior"),
    logical(1)
  ))
  competitive <- noninterior[vapply(polished[noninterior], function(x) {
    length(x$objective) == 1L && isTRUE(is.finite(x$objective)) &&
      isTRUE(x$objective <= best + tolerance)
  }, logical(1))]
  if (length(competitive)) {
    return(cn_failure(
      "competitive_noninterior_restricted_mode",
      baseline_fit = baseline,
      competitive_indices = competitive,
      competitive_objectives = vapply(
        polished[competitive], `[[`, numeric(1), "objective"
      ),
      selected_repaired_index = selected_index,
      selected_repaired_objective = best,
      candidate_repairs = polished,
      stage = "restricted_fit"
    ))
  }

  near <- resolved[vapply(polished[resolved], function(x) {
    x$objective - best <= tolerance
  }, logical(1))]
  if (length(near) > 1L) {
    for (index in setdiff(near, selected_index)) {
      other <- polished[[index]]
      comparison <- cn_near_mode_comparison(
        selected$candidate, other$candidate,
        selected$evaluation, other$evaluation
      )
      if (isTRUE(comparison$distinct)) {
        return(cn_failure(
          "distinct_near_equal_restricted_modes",
          baseline_fit = baseline,
          selected_repaired_index = selected_index,
          competing_repaired_index = index,
          parameter_difference = comparison$parameter_difference,
          fitted_difference = comparison$fitted_difference,
          candidate_repairs = polished,
          stage = "restricted_fit"
        ))
      }
    }
  }

  repaired_candidates <- baseline$candidates
  for (index in resolved) {
    repaired_candidates[[index]] <- polished[[index]]$candidate
  }
  final <- selected$state
  structure(list(
    available = TRUE,
    failure_code = "",
    theta = final$theta,
    opt_par = final$opt_par,
    objective = final$evaluation$objective,
    convergence = 0L,
    evaluation = final$evaluation,
    integration_audit = selected$audit,
    layout = layout,
    bounds = bounds,
    constraints = constraints,
    final_constraint_slack = selected$slack,
    eta_boundary = selected$eta_boundary,
    candidates = repaired_candidates,
    barrier_candidates = baseline$candidates,
    selected_index = selected_index,
    n_converged = sum(converged),
    penalty_evaluations = if (
      length(baseline$penalty_evaluations) == 1L &&
      is.finite(baseline$penalty_evaluations)
    ) as.integer(baseline$penalty_evaluations) else NA_integer_,
    stationarity = final$stationarity$maximum,
    stationarity_components = final$stationarity$components,
    nuisance_score_norm = final$stationarity$score_norm,
    integration_level = "production",
    interior_polish = list(
      attempted = TRUE,
      baseline_failure = baseline$failure_code,
      candidate_repairs = polished,
      selected_index = selected_index
    )
  ), class = c("cn_d3_restricted_fit", "list"))
}

cn_fit_occupancy_null <- function(
    data, integration_level = c("production", "audit")) {
  integration_level <- match.arg(integration_level)
  precheck <- cn_precheck(data)
  if (!isTRUE(precheck$available)) return(precheck)
  layout <- cn_layout(data)
  bounds <- cn_opt_bounds(layout)
  constraints <- cn_optimizer_constraints(data, layout, bounds)
  starts <- cn_start_set(data, layout)
  penalty_evaluations <- 0L
  fit_one <- function(start) {
    eta_trial_violations <- list()
    evaluate <- function(theta) {
      tryCatch(
        cn_eval(theta, data, integration_level = integration_level),
        cn_eta_envelope_error = function(e) {
          eta_trial_violations[[length(eta_trial_violations) + 1L]] <<-
            list(
              trial_eta_range = e$fitted_eta_range,
              certified_eta_predictor = e$certified_eta_predictor
            )
          e
        },
        error = function(e) e
      )
    }
    objective <- function(par, ...) {
      theta <- cn_opt_to_natural(setNames(par, layout$opt_names), layout)
      evaluation <- evaluate(theta)
      if (inherits(evaluation, "error") || !is.finite(evaluation$objective)) {
        penalty_evaluations <<- penalty_evaluations + 1L
        return(CN_D3_CONTRACT$optimizer$invalid_objective_penalty)
      }
      evaluation$objective
    }
    gradient <- function(par, ...) {
      theta <- cn_opt_to_natural(setNames(par, layout$opt_names), layout)
      evaluation <- evaluate(theta)
      if (inherits(evaluation, "error")) stop(evaluation)
      -cn_score_natural_to_opt(evaluation$total_score, layout)
    }
    result <- tryCatch(
      constrOptim(
        theta = as.numeric(start), f = objective, grad = gradient,
        ui = constraints$ui, ci = constraints$ci,
        mu = CN_D3_CONTRACT$optimizer$barrier_mu,
        method = CN_D3_CONTRACT$optimizer$inner_method,
        outer.iterations =
          CN_D3_CONTRACT$optimizer$barrier_outer_iterations,
        outer.eps = CN_D3_CONTRACT$optimizer$barrier_outer_tolerance,
        control = list(
          maxit = CN_D3_CONTRACT$optimizer$maxit,
          reltol = CN_D3_CONTRACT$optimizer$relative_tolerance
        )
      ),
      error = function(e) structure(
        list(error = conditionMessage(e)), class = "cn_optimizer_error"
      )
    )
    if (inherits(result, "cn_optimizer_error")) {
      return(list(
        error = result$error,
        eta_trial_violations = eta_trial_violations
      ))
    }
    opt_par <- setNames(result$par, layout$opt_names)
    list(
      opt_par = opt_par,
      theta = cn_opt_to_natural(opt_par, layout),
      objective = result$value,
      convergence = result$convergence,
      message = result$message,
      counts = result$counts,
      outer_iterations = result$outer.iterations,
      barrier_value = result$barrier.value,
      constraint_slack = cn_constraint_slack(opt_par, constraints),
      eta_trial_violations = eta_trial_violations
    )
  }
  candidates <- lapply(starts, fit_one)
  finite <- vapply(candidates, function(candidate) {
    is.null(candidate$error) && is.finite(candidate$objective) &&
      all(is.finite(candidate$theta))
  }, logical(1))
  eta_trial_details <- function() {
    eta_summary <- cn_eta_trial_violation_summary(candidates)
    list(
      trial_eta_range = eta_summary$trial_eta_range,
      certified_eta_predictor = eta_summary$certified_eta_predictor,
      eta_trial_violations = eta_summary$violations
    )
  }
  if (!any(finite)) {
    return(cn_failure(
      "all_optimizer_candidates_failed", candidates = candidates,
      penalty_evaluations = penalty_evaluations,
      optimizer_trial_diagnostics = eta_trial_details(),
      stage = "restricted_fit"
    ))
  }
  converged <- finite & vapply(candidates, function(candidate) {
    isTRUE(candidate$convergence == 0L)
  }, logical(1))
  if (!any(converged)) {
    return(cn_failure(
      "nonconverged_restricted_fit", candidates = candidates,
      penalty_evaluations = penalty_evaluations,
      optimizer_trial_diagnostics = eta_trial_details(),
      stage = "restricted_fit"
    ))
  }
  converged_index <- which(converged)
  objective <- vapply(candidates[converged_index], `[[`, numeric(1), "objective")
  selected_index <- converged_index[which.min(objective)]
  selected <- candidates[[selected_index]]

  best <- selected$objective
  better_nonconverged <- cn_better_nonconverged_indices(
    candidates, finite, converged, best
  )
  if (length(better_nonconverged)) {
    return(cn_failure(
      "better_nonconverged_candidate",
      unresolved_indices = better_nonconverged,
      unresolved_objectives = vapply(
        candidates[better_nonconverged], `[[`, numeric(1), "objective"
      ),
      candidates = candidates, selected_index = selected_index,
      stage = "restricted_fit"
    ))
  }
  near_index <- converged_index[vapply(candidates[converged_index], function(candidate) {
    candidate$objective - best <=
      CN_D3_CONTRACT$optimizer$near_mode_objective_relative_tolerance * (1 + abs(best))
  }, logical(1))]
  if (length(near_index) > 1L) {
    reference_eval <- cn_eval(
      selected$theta, data, integration_level = integration_level
    )
    for (index in setdiff(near_index, selected_index)) {
      other <- candidates[[index]]
      other_eval <- cn_eval(
        other$theta, data, integration_level = integration_level
      )
      comparison <- cn_near_mode_comparison(
        selected, other, reference_eval, other_eval
      )
      if (isTRUE(comparison$distinct)) {
        return(cn_failure(
          "distinct_near_equal_restricted_modes",
          candidates = candidates, selected_index = selected_index,
          parameter_difference = comparison$parameter_difference,
          fitted_difference = comparison$fitted_difference,
          penalty_evaluations = penalty_evaluations,
          stage = "restricted_fit"
        ))
      }
    }
  }

  active <- cn_active_boundary_parameters(selected$opt_par, bounds)
  if (length(active)) {
    return(cn_failure(
      "active_restricted_boundary", active_parameters = active,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, stage = "restricted_fit"
    ))
  }

  evaluation <- tryCatch(
    cn_eval(
      selected$theta, data, integration_level = integration_level
    ),
    error = function(e) e
  )
  if (inherits(evaluation, "error")) {
    if (inherits(evaluation, "cn_eta_envelope_error")) {
      return(cn_failure(
        "present_location_predictor_outside_certified_envelope",
        fitted_eta_range = evaluation$fitted_eta_range,
        certified_eta_predictor = evaluation$certified_eta_predictor,
        candidates = candidates, selected_index = selected_index,
        theta = selected$theta, stage = "restricted_fit"
      ))
    }
    return(cn_failure(
      "final_restricted_evaluation_failed",
      evaluation_error = conditionMessage(evaluation),
      candidates = candidates, selected_index = selected_index,
      stage = "restricted_fit"
    ))
  }
  eta_envelope <- cn_eta_envelope_check(evaluation$eta)
  if (!isTRUE(eta_envelope$ok)) {
    return(cn_failure(
      "present_location_predictor_outside_certified_envelope",
      fitted_eta_range = eta_envelope$fitted_range,
      certified_eta_predictor = eta_envelope$envelope,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, evaluation = evaluation,
      stage = "restricted_fit"
    ))
  }
  final_constraint_slack <- cn_constraint_slack(selected$opt_par, constraints)
  if (any(!is.finite(final_constraint_slack)) ||
      any(final_constraint_slack <= 0)) {
    return(cn_failure(
      "restricted_fit_constraint_violation",
      constraint_slack = final_constraint_slack,
      fitted_eta_range = eta_envelope$fitted_range,
      certified_eta_predictor = eta_envelope$envelope,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, evaluation = evaluation,
      stage = "restricted_fit"
    ))
  }
  eta_boundary <- cn_eta_boundary_status(
    evaluation$eta, data$subject_id
  )
  if (isTRUE(eta_boundary$active)) {
    return(cn_failure(
      "active_present_location_predictor_envelope",
      fitted_eta_range = eta_boundary$fitted_eta_range,
      certified_eta_predictor = eta_boundary$certified_eta_predictor,
      minimum_eta_lower_slack = eta_boundary$minimum_lower_slack,
      minimum_eta_upper_slack = eta_boundary$minimum_upper_slack,
      active_eta_lower_rows = eta_boundary$active_lower,
      active_eta_upper_rows = eta_boundary$active_upper,
      active_eta_subject_id = eta_boundary$active_subject_id,
      eta_boundary_limit = eta_boundary$boundary_limit,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, evaluation = evaluation,
      stage = "restricted_fit"
    ))
  }
  integration_audit <- NULL
  if (integration_level == "production") {
    integration_audit <- tryCatch(
      cn_integration_audit(selected$theta, data), error = function(e) e
    )
    if (inherits(integration_audit, "error")) {
      return(cn_failure(
        "integration_audit_failed",
        integration_audit_error = conditionMessage(integration_audit),
        candidates = candidates, selected_index = selected_index,
        theta = selected$theta, evaluation = evaluation,
        stage = "restricted_fit"
      ))
    }
    audit_contract <- CN_D3_CONTRACT$integration
    if (integration_audit$maximum_subject_loglik_difference >
          audit_contract$maximum_subject_loglik_difference ||
        integration_audit$maximum_subject_score_difference >
          audit_contract$maximum_subject_score_difference ||
        integration_audit$maximum_total_score_difference >
          audit_contract$maximum_total_score_difference) {
      return(cn_failure(
        "integration_audit_disagreement",
        integration_audit = integration_audit,
        candidates = candidates, selected_index = selected_index,
        theta = selected$theta, evaluation = evaluation,
        stage = "restricted_fit"
      ))
    }
  }
  # This is the original likelihood score, not the constrOptim barrier
  # gradient.  Together with inactive constraints it ensures that a formed
  # solution is an interior stationary point of the unchanged D3 model.
  nuisance_names <- setdiff(layout$natural_names, "delta")
  stationarity_result <- cn_relative_stationarity(
    evaluation$score, nuisance_names
  )
  if (!isTRUE(stationarity_result$valid)) {
    return(cn_failure(
      "invalid_nuisance_score_norm",
      nuisance_score_norm = stationarity_result$score_norm,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, evaluation = evaluation,
      stage = "restricted_fit"
    ))
  }
  nuisance_score_norm <- stationarity_result$score_norm
  stationarity_components <- stationarity_result$components
  stationarity <- stationarity_result$maximum
  if (!is.finite(stationarity) ||
      stationarity >
      CN_D3_CONTRACT$optimizer$stationarity_relative_score_norm) {
    baseline_failure <- cn_failure(
      "nonstationary_restricted_fit", stationarity = stationarity,
      stationarity_components = stationarity_components,
      nuisance_score_norm = nuisance_score_norm,
      candidates = candidates, selected_index = selected_index,
      theta = selected$theta, evaluation = evaluation,
      stage = "restricted_fit"
    )
    if (integration_level != "production") return(baseline_failure)
    return(tryCatch(
      cn_repair_nonstationary_restricted_fit(
        baseline_failure, data
      ),
      error = function(e) cn_failure(
        "restricted_interior_polish_exception",
        baseline_fit = baseline_failure,
        interior_polish_error = conditionMessage(e),
        stage = "restricted_fit"
      )
    ))
  }
  structure(list(
    available = TRUE,
    failure_code = "",
    theta = selected$theta,
    opt_par = selected$opt_par,
    objective = selected$objective,
    convergence = selected$convergence,
    evaluation = evaluation,
    integration_audit = integration_audit,
    layout = layout,
    bounds = bounds,
    constraints = constraints,
    final_constraint_slack = final_constraint_slack,
    eta_boundary = eta_boundary,
    candidates = candidates,
    selected_index = selected_index,
    n_converged = sum(converged),
    penalty_evaluations = penalty_evaluations,
    stationarity = stationarity,
    stationarity_components = stationarity_components,
    nuisance_score_norm = nuisance_score_norm,
    integration_level = integration_level
  ), class = c("cn_d3_restricted_fit", "list"))
}

cn_sensitivity <- function(
    theta, data, integration_level = c("production", "audit")) {
  integration_level <- match.arg(integration_level)
  layout <- cn_layout(data)
  theta <- cn_named_theta(theta, layout)
  score_map <- function(x) {
    x <- setNames(x, layout$natural_names)
    as.numeric(cn_eval(
      x, data, integration_level = integration_level
    )$total_score)
  }
  jacobian <- numDeriv::jacobian(
    func = score_map, x = as.numeric(theta), method = "Richardson",
    method.args = list(eps = 1e-4, d = 0.1, r = 4L, v = 2L)
  )
  dimnames(jacobian) <- list(layout$natural_names, layout$natural_names)
  A_raw <- -jacobian
  asymmetry <- max(abs(A_raw - t(A_raw))) / max(1, max(abs(A_raw)))
  A <- (A_raw + t(A_raw)) / 2
  list(A = A, A_raw = A_raw, jacobian = jacobian, asymmetry = asymmetry)
}

cn_matrix_health <- function(M, label) {
  if (!is.matrix(M) || nrow(M) != ncol(M) || any(!is.finite(M))) {
    return(list(ok = FALSE, failure_code = paste0("nonfinite_", label)))
  }
  diagonal <- diag(M)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    return(list(ok = FALSE, failure_code = paste0("nonpositive_", label)))
  }
  D_inverse <- diag(1 / sqrt(diagonal), nrow = length(diagonal))
  scaled <- D_inverse %*% M %*% D_inverse
  values <- tryCatch(
    eigen((scaled + t(scaled)) / 2, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(scaled))
  )
  minimum <- if (length(values)) min(values) else NA_real_
  condition <- if (all(is.finite(values)) && minimum > 0) max(values) / minimum else Inf
  ok <- all(is.finite(values)) &&
    minimum > CN_D3_CONTRACT$sensitivity$scaled_minimum_eigenvalue &&
    condition < CN_D3_CONTRACT$sensitivity$maximum_condition_number
  list(
    ok = ok,
    failure_code = if (ok) "" else paste0("weak_or_singular_", label),
    scaled_minimum_eigenvalue = minimum,
    condition_number = condition,
    scaled_matrix = scaled
  )
}

cn_project_occupancy_score <- function(
    fit, data = NULL, sensitivity_override = NULL) {
  if (!isTRUE(fit$available)) return(fit)
  if (is.null(data)) cn_stop("data is required for occupancy projection")
  evaluation <- fit$evaluation
  layout <- fit$layout
  sensitivity <- if (is.null(sensitivity_override)) {
    tryCatch(
      cn_sensitivity(
        fit$theta, data, integration_level = fit$integration_level
      ),
      error = function(e) e
    )
  } else {
    sensitivity_override
  }
  if (inherits(sensitivity, "error")) {
    return(cn_failure(
      "sensitivity_jacobian_failed",
      sensitivity_error = conditionMessage(sensitivity), fit = fit,
      stage = "score_projection"
    ))
  }
  if (!is.finite(sensitivity$asymmetry) ||
      sensitivity$asymmetry > CN_D3_CONTRACT$sensitivity$asymmetry_tolerance) {
    return(cn_failure(
      "full_sensitivity_asymmetry",
      sensitivity_asymmetry = sensitivity$asymmetry,
      sensitivity = sensitivity, fit = fit, stage = "score_projection"
    ))
  }
  target <- "delta"
  nuisance <- setdiff(layout$natural_names, target)
  A <- sensitivity$A
  A_nn <- A[nuisance, nuisance, drop = FALSE]
  A_tn <- A[target, nuisance, drop = FALSE]
  A_nt <- A[nuisance, target, drop = FALSE]
  nuisance_health <- cn_matrix_health(A_nn, "nuisance_sensitivity")
  if (!nuisance_health$ok) {
    return(cn_failure(
      nuisance_health$failure_code,
      nuisance_health = nuisance_health, sensitivity = sensitivity,
      fit = fit, stage = "score_projection"
    ))
  }
  rhs <- as.numeric(t(A_tn))
  C_cholesky <- tryCatch({
    factor <- chol(t(A_nn))
    backsolve(factor, forwardsolve(t(factor), rhs))
  }, error = function(e) rep(NA_real_, length(rhs)))
  C_qr <- tryCatch(
    qr.solve(t(A_nn), rhs, tol = 1e-10),
    error = function(e) rep(NA_real_, length(rhs))
  )
  solve_agreement <- max(abs(C_cholesky - C_qr)) /
    max(1, max(abs(C_qr)))
  if (any(!is.finite(C_cholesky)) || any(!is.finite(C_qr)) ||
      !is.finite(solve_agreement) ||
      solve_agreement > CN_D3_CONTRACT$sensitivity$solve_relative_agreement) {
    return(cn_failure(
      "nuisance_solve_disagreement", solve_agreement = solve_agreement,
      nuisance_health = nuisance_health, sensitivity = sensitivity,
      fit = fit, stage = "score_projection"
    ))
  }
  C <- setNames(C_cholesky, nuisance)
  phi <- as.numeric(evaluation$score[, target] -
    evaluation$score[, nuisance, drop = FALSE] %*% C)
  U <- sum(phi)
  U_direct <- evaluation$total_score[target] -
    sum(C * evaluation$total_score[nuisance])
  score_identity <- abs(U - U_direct) / max(1, abs(U))
  if (any(!is.finite(phi)) || !is.finite(U)) {
    return(cn_failure(
      "nonfinite_projected_score", fit = fit, sensitivity = sensitivity,
      stage = "score_projection"
    ))
  }
  if (!is.finite(score_identity) ||
      score_identity > CN_D3_CONTRACT$sensitivity$score_identity_relative_tolerance) {
    return(cn_failure(
      "stack_score_identity_failure", score_identity = score_identity,
      fit = fit, sensitivity = sensitivity, stage = "score_projection"
    ))
  }
  A_eff <- as.numeric(A[target, target] - C %*% A_nt)
  structure(list(
    available = TRUE,
    failure_code = "",
    phi = phi,
    subject_id = data$subject_id,
    U = U,
    C = C,
    A_eff_diagnostic = A_eff,
    score_identity = score_identity,
    solve_agreement = solve_agreement,
    nuisance_health = nuisance_health,
    sensitivity = sensitivity,
    fit = fit
  ), class = c("cn_d3_projected_score", "list"))
}

