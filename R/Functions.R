# Copyright (C) 2026 Yiqian Zhang
# SPDX-License-Identifier: GPL-3.0-or-later

#' DASH: Depth-Aware Structural-absence and Hellinger Framework
#'
#' @description
#' DASH provides two-component differential-abundance inference for sparse
#' microbiome count data. For each retained taxon, it reports a depth-aware,
#' nuisance-adjusted structural-absence score; a posterior-presence-weighted
#' Hellinger-Riemann intrinsic-coordinate abundance contrast; and a primary
#' fixed two-component Bonferroni min-P omnibus p-value.
#'
#' @section Main functions:
#' - [dash()] runs the differential-abundance analysis.
#' - [hric_transform()] computes Hellinger-Riemann intrinsic coordinates.
#' - [dash_bootstrap()] performs a stratified full-pipeline bootstrap
#'   sensitivity analysis for the abundance component.
#' - [dash_label_swap_check()] checks numerical invariance to reversing the
#'   binary group coding.
#'
#' @section Statistical scope:
#' The current implementation supports independent samples and a binary group
#' comparison with optional covariate adjustment. Structural absence is a
#' model-implied latent state rather than a directly observed biological fact.
#' The positive-state model is an upper-truncated-normal working likelihood for
#' transformed integer counts. The analytic HC3 abundance variance is
#' conditional on fitted auxiliary quantities and the cross-fitted affine
#' background; [dash_bootstrap()] is provided as a full-pipeline sensitivity
#' analysis.
#'
#' @keywords internal
#' @md
"_PACKAGE"

# Internal utilities

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Validate and coerce an integer-valued scalar
#' @keywords internal
#' @noRd
validate_integer_scalar <- function(x, name, minimum = 0L) {
  if (
    length(x) != 1L ||
      !is.numeric(x) ||
      !is.finite(x) ||
      abs(x - round(x)) > sqrt(.Machine$double.eps) ||
      x < minimum ||
      x > .Machine$integer.max
  ) {
    stop(
      sprintf("`%s` must be one integer-valued number >= %s.", name, minimum),
      call. = FALSE
    )
  }

  as.integer(round(x))
}

#' Validate a non-missing logical scalar
#' @keywords internal
#' @noRd
validate_flag <- function(x, name) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }

  x
}

#' Validate zero-model optimization controls
#' @keywords internal
#' @noRd
validate_fit_control <- function(control) {
  defaults <- list(
    maxit = 1000L,
    factr = 1e7,
    nm_maxit = 1000L,
    nm_reltol = 1e-9,
    gradient_mean_tol = 1e-5,
    boundary_tol = 1e-6,
    information_rel_tol = 1e-10,
    expand_bounds = TRUE
  )

  if (!is.list(control)) {
    stop("`fit_control` must be a list.", call. = FALSE)
  }
  if (length(control) > 0L) {
    if (is.null(names(control)) || any(names(control) == "")) {
      stop("Every `fit_control` element must be named.", call. = FALSE)
    }
    if (anyDuplicated(names(control))) {
      stop("`fit_control` contains duplicated names.", call. = FALSE)
    }
    unknown <- setdiff(names(control), names(defaults))
    if (length(unknown) > 0L) {
      stop(
        "Unknown `fit_control` element(s): ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
  }

  out <- utils::modifyList(defaults, control)
  out$maxit <- validate_integer_scalar(out$maxit, "fit_control$maxit", 1L)
  out$nm_maxit <- validate_integer_scalar(
    out$nm_maxit,
    "fit_control$nm_maxit",
    1L
  )

  validate_positive <- function(x, name, upper = Inf) {
    if (
      length(x) != 1L ||
        !is.numeric(x) ||
        !is.finite(x) ||
        x <= 0 ||
        x > upper
    ) {
      stop(
        sprintf("`%s` must be one finite number in (0, %s].", name, upper),
        call. = FALSE
      )
    }
    as.numeric(x)
  }

  out$factr <- validate_positive(out$factr, "fit_control$factr")
  out$nm_reltol <- validate_positive(
    out$nm_reltol,
    "fit_control$nm_reltol",
    1
  )
  out$gradient_mean_tol <- validate_positive(
    out$gradient_mean_tol,
    "fit_control$gradient_mean_tol",
    1
  )
  out$boundary_tol <- validate_positive(
    out$boundary_tol,
    "fit_control$boundary_tol",
    1
  )
  out$information_rel_tol <- validate_positive(
    out$information_rel_tol,
    "fit_control$information_rel_tol",
    1
  )
  out$expand_bounds <- validate_flag(
    out$expand_bounds,
    "fit_control$expand_bounds"
  )
  out
}

#' Clamp a vector to a closed interval
#' @keywords internal
#' @noRd
clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

#' Numerically stable log(exp(a) + exp(b))
#' @keywords internal
#' @noRd
logspace_add <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  both_neg_inf <- is.infinite(m) & m < 0
  out[both_neg_inf] <- -Inf
  out
}

#' Stable lower-tail inverse Mills ratio phi(x) / Phi(x)
#' @keywords internal
#' @noRd
stable_mills <- function(x) {
  x <- as.numeric(x)
  out <- exp(stats::dnorm(x, log = TRUE) - stats::pnorm(x, log.p = TRUE))

  bad <- !is.finite(out)
  if (any(bad)) {
    ## For a very negative argument, phi(x)/Phi(x) = -x + O(1/|x|).
    xb <- x[bad]
    approx <- ifelse(
      xb < 0,
      -xb + 1 / pmax(-xb, 1),
      0
    )
    out[bad] <- approx
  }

  pmax(out, 0)
}

#' Run an expression with a deterministic seed without changing the caller RNG
#' @keywords internal
#' @noRd
with_preserved_seed <- function(seed, expr) {
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

#' Encode a binary group as 0/1
#' @keywords internal
#' @noRd
encode_group <- function(group, n, case_level = NULL) {
  if (length(group) != n) {
    stop("`group` length must equal the number of samples.")
  }
  if (anyNA(group)) {
    stop("`group` contains missing values.")
  }

  if (is.factor(group) || is.character(group) || is.logical(group)) {
    x <- as.character(group)
    lev <- if (is.factor(group)) levels(group) else sort(unique(x))
    lev <- lev[lev %in% unique(x)]

    if (length(lev) != 2L) {
      stop("`group` must have exactly two observed levels.")
    }

    if (!is.null(case_level)) {
      case_level <- as.character(case_level)
      if (!case_level %in% lev) {
        stop("`case_level` is not an observed group level.")
      }
      control_level <- setdiff(lev, case_level)
      lev <- c(control_level, case_level)
    }

    g <- as.integer(x == lev[2L])
    return(list(g = g, levels = lev))
  }

  gv <- as.numeric(group)
  if (any(!is.finite(gv))) {
    stop("`group` contains non-finite values.")
  }
  lev_num <- sort(unique(gv))
  if (length(lev_num) != 2L) {
    stop("`group` must take exactly two distinct values.")
  }

  if (!is.null(case_level)) {
    case_num <- as.numeric(case_level)
    if (length(case_num) != 1L || !is.finite(case_num) || !case_num %in% lev_num) {
      stop("`case_level` is not an observed numeric group value.")
    }
    lev_num <- c(setdiff(lev_num, case_num), case_num)
  }

  list(
    g = as.integer(gv == lev_num[2L]),
    levels = as.character(lev_num)
  )
}

#' Convert covariates to a numeric design block
#'
#' Continuous columns are centered and scaled. Binary indicator columns are
#' retained on their original 0/1 or two-value scale. Factors are expanded with
#' treatment contrasts by model.matrix(). Constant columns are removed.
#' @keywords internal
#' @noRd
prepare_covariates <- function(covariates, n) {
  if (is.null(covariates)) {
    return(NULL)
  }

  if (is.atomic(covariates) && is.null(dim(covariates))) {
    covariates <- data.frame(covariate = covariates)
  }

  if (is.data.frame(covariates)) {
    if (nrow(covariates) != n) {
      stop("`covariates` must have one row per sample.")
    }
    if (anyNA(covariates)) {
      stop("`covariates` contains missing values; handle missingness before DASH.")
    }

    dat <- covariates
    for (j in seq_along(dat)) {
      if (is.character(dat[[j]])) {
        dat[[j]] <- factor(dat[[j]])
      }
      if (is.ordered(dat[[j]])) {
        dat[[j]] <- factor(
          as.character(dat[[j]]),
          levels = levels(dat[[j]])
        )
      }
      if (is.factor(dat[[j]])) {
        dat[[j]] <- droplevels(dat[[j]])
      }
    }

    constant <- vapply(
      dat,
      function(x) length(unique(x)) < 2L,
      logical(1)
    )
    if (any(constant)) {
      warning(
        "Removing constant covariate column(s): ",
        paste(names(dat)[constant], collapse = ", "),
        call. = FALSE
      )
      dat <- dat[, !constant, drop = FALSE]
    }
    if (ncol(dat) == 0L) {
      return(NULL)
    }

    factor_names <- names(dat)[vapply(dat, is.factor, logical(1))]
    contrasts_arg <- if (length(factor_names) == 0L) {
      NULL
    } else {
      out <- lapply(dat[factor_names], function(x) {
        stats::contr.treatment(levels(x), base = 1L)
      })
      names(out) <- factor_names
      out
    }

    Z <- stats::model.matrix(
      ~ .,
      data = dat,
      na.action = stats::na.fail,
      contrasts.arg = contrasts_arg
    )
    if ("(Intercept)" %in% colnames(Z)) {
      Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    }
  } else {
    Z <- as.matrix(covariates)
    if (nrow(Z) != n) {
      stop("`covariates` must have one row per sample.")
    }
    if (!is.numeric(Z)) {
      stop("Matrix covariates must be numeric; use a data frame for factors.")
    }
  }

  if (ncol(Z) == 0L) {
    return(NULL)
  }
  if (anyNA(Z) || any(!is.finite(Z))) {
    stop("`covariates` contains missing or non-finite values.")
  }

  Z <- matrix(
    as.numeric(Z),
    nrow = nrow(Z),
    ncol = ncol(Z),
    dimnames = dimnames(Z)
  )

  if (is.null(colnames(Z))) {
    colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  }
  colnames(Z) <- make.unique(make.names(colnames(Z)))

  variable <- vapply(
    seq_len(ncol(Z)),
    function(j) {
      s <- stats::sd(Z[, j])
      is.finite(s) && s > 0
    },
    logical(1)
  )

  if (any(!variable)) {
    warning(
      "Removing constant covariate column(s): ",
      paste(colnames(Z)[!variable], collapse = ", "),
      call. = FALSE
    )
    Z <- Z[, variable, drop = FALSE]
  }

  if (ncol(Z) == 0L) {
    return(NULL)
  }

  centers <- numeric(ncol(Z))
  scales <- rep(1, ncol(Z))

  for (j in seq_len(ncol(Z))) {
    uj <- sort(unique(Z[, j]))
    if (length(uj) == 2L) {
      ## Put every binary numeric indicator on a common 0/1 scale.
      centers[j] <- uj[1L]
      scales[j] <- uj[2L] - uj[1L]
      Z[, j] <- as.numeric(Z[, j] == uj[2L])
    } else if (length(uj) > 2L) {
      centers[j] <- mean(Z[, j])
      scales[j] <- stats::sd(Z[, j])
      Z[, j] <- (Z[, j] - centers[j]) / scales[j]
    }
  }

  attr(Z, "centers") <- centers
  attr(Z, "scales") <- scales
  Z
}

#' Construct an intercept/group/covariate design matrix
#' @keywords internal
#' @noRd
design_matrix <- function(g = NULL, z = NULL, include_group = TRUE) {
  n <- if (!is.null(g)) {
    length(g)
  } else if (!is.null(z)) {
    nrow(z)
  } else {
    stop("At least one of `g` or `z` must determine the number of rows.")
  }

  X <- matrix(1, nrow = n, ncol = 1L)
  colnames(X) <- "Intercept"

  if (include_group) {
    if (is.null(g)) {
      stop("`g` must be supplied when include_group = TRUE.")
    }
    X <- cbind(X, Group = as.numeric(g))
  }

  if (!is.null(z)) {
    X <- cbind(X, z)
  }

  X
}

#' Check full column rank
#' @keywords internal
#' @noRd
is_full_rank <- function(X, tol = 1e-10) {
  if (nrow(X) < ncol(X)) {
    return(FALSE)
  }
  qr(X, tol = tol)$rank == ncol(X)
}

#' Invert a symmetric positive-definite matrix with an eigenvalue check
#' @keywords internal
#' @noRd
invert_spd <- function(A, rel_tol = 1e-10) {
  A <- (A + t(A)) / 2
  ee <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)

  invalid <- list(
    valid = FALSE,
    inverse = NULL,
    condition = Inf,
    min_eigen = NA_real_,
    max_eigen = NA_real_
  )

  if (is.null(ee) || any(!is.finite(ee$values))) {
    return(invalid)
  }

  max_eig <- max(ee$values)
  min_eig <- min(ee$values)

  if (
    !is.finite(max_eig) || max_eig <= 0 ||
      !is.finite(min_eig) || min_eig <= rel_tol * max_eig
  ) {
    invalid$min_eigen <- min_eig
    invalid$max_eigen <- max_eig
    return(invalid)
  }

  inv <- ee$vectors %*%
    (diag(1 / ee$values, nrow = length(ee$values)) %*% t(ee$vectors))

  list(
    valid = TRUE,
    inverse = inv,
    condition = max_eig / min_eig,
    min_eigen = min_eig,
    max_eigen = max_eig
  )
}


# P-value combinations

#' Cauchy combination of p-values
#' @keywords internal
#' @noRd
cauchy_combination <- function(ps, eps = 1e-15) {
  ps <- as.numeric(ps)
  if (length(ps) == 0L || any(!is.finite(ps))) {
    return(NA_real_)
  }
  ps <- clamp(ps, eps, 1 - eps)
  out <- 0.5 - atan(mean(tan((0.5 - ps) * pi))) / pi
  clamp(out, 0, 1)
}

#' Fixed two-component Bonferroni min-P combination
#' @keywords internal
#' @noRd
bonferroni_two_component <- function(p_absence, p_abundance) {
  pa <- if (is.finite(p_absence)) p_absence else 1
  pb <- if (is.finite(p_abundance)) p_abundance else 1
  min(1, 2 * min(pa, pb))
}


# Hellinger-Riemann intrinsic coordinates

#' Compute Hellinger-Riemann intrinsic coordinates
#'
#' Closes each row to a composition, applies the square-root map, and evaluates
#' the spherical logarithm map at the uniform composition. Observed zeros remain
#' exact boundary values; no pseudocount is added.
#'
#' @param X A non-negative numeric matrix with samples in rows and components in
#'   columns. Every row must have a positive total and at least two components
#'   must be supplied.
#'
#' @return A numeric matrix with the same dimensions and dimnames as `X`.
#'
#' @examples
#' x <- rbind(
#'   sample1 = c(taxon_a = 3, taxon_b = 1, taxon_c = 0),
#'   sample2 = c(taxon_a = 0, taxon_b = 2, taxon_c = 5)
#' )
#' hric_transform(x)
#'
#' @family DASH functions
#' @export
#' @md
hric_transform <- function(X) {
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("`X` must be numeric.")
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` contains missing or non-finite values.")
  }
  if (any(X < 0)) {
    stop("`X` must contain non-negative values only.")
  }
  if (ncol(X) < 2L) {
    stop("`X` must contain at least two components.")
  }

  rs <- rowSums(X)
  if (any(rs <= 0)) {
    stop("Every row of `X` must have positive total abundance.")
  }

  Pi <- sweep(X, 1L, rs, FUN = "/")
  sqrt_Pi <- sqrt(pmax(Pi, 0))
  k <- ncol(Pi)
  sqrt_pi0 <- rep(1 / sqrt(k), k)

  inner <- as.vector(sqrt_Pi %*% sqrt_pi0)
  inner <- clamp(inner, 0, 1)
  radius <- sqrt(pmax(0, 1 - inner^2))

  scale_factor <- numeric(length(radius))
  small <- radius < 1e-8
  r2 <- radius[small]^2
  scale_factor[small] <- 1 + r2 / 6 + 3 * r2^2 / 40
  scale_factor[!small] <- asin(radius[!small]) / radius[!small]

  centered <- sqrt_Pi - outer(inner, sqrt_pi0)
  out <- sweep(centered, 1L, scale_factor, FUN = "*")
  dimnames(out) <- dimnames(X)
  out
}


# Depth-aware structural-absence working model

#' Evaluate the upper-truncated-normal zero model
#'
#' Conditional on presence, L is normal with mean eta and standard deviation
#' sigma, truncated to (-Inf, 0]. An observed zero is structural absence or a
#' present state with L <= a_i = log(d0 / D_i). A positive transformed count is
#' ell_i = log(y_i / D_i) in \eqn{[a_i,0]}..
#'
#' @keywords internal
#' @noRd
zero_model_state <- function(
    par,
    y,
    D,
    X_rho,
    X_eta,
    d0 = 1,
    need_hessian = FALSE) {

  n <- length(y)
  p_rho <- ncol(X_rho)
  p_eta <- ncol(X_eta)
  expected_length <- p_rho + p_eta + 1L

  invalid <- function(reason) {
    list(valid = FALSE, reason = reason, nll = 1e100)
  }

  if (length(par) != expected_length || any(!is.finite(par))) {
    return(invalid("invalid_parameter_vector"))
  }
  if (length(D) != n || any(!is.finite(D)) || any(D <= d0)) {
    return(invalid("invalid_library_size_or_detection_boundary"))
  }

  pos <- y > 0
  a <- log(d0 / D)
  ell <- numeric(n)
  ell[pos] <- log(y[pos] / D[pos])

  if (any(pos & (ell < a - 1e-10 | ell > 1e-10))) {
    return(invalid("positive_observation_outside_working_support"))
  }

  alpha <- par[seq_len(p_rho)]
  eta_coef <- par[p_rho + seq_len(p_eta)]
  omega <- par[p_rho + p_eta + 1L]
  sigma <- exp(omega)

  if (!is.finite(sigma) || sigma <= 0) {
    return(invalid("invalid_scale"))
  }

  q_rho <- as.numeric(X_rho %*% alpha)
  rho <- stats::plogis(q_rho)
  eta <- as.numeric(X_eta %*% eta_coef)

  u <- -eta / sigma
  tval <- (a - eta) / sigma
  zval <- (ell - eta) / sigma

  log_Phi_u <- stats::pnorm(u, log.p = TRUE)
  log_Phi_t <- stats::pnorm(tval, log.p = TRUE)
  log_F <- pmin(0, log_Phi_t - log_Phi_u)

  log_rho <- stats::plogis(q_rho, log.p = TRUE)
  log_one_minus_rho <- stats::plogis(
    q_rho,
    lower.tail = FALSE,
    log.p = TRUE
  )

  log_L_zero <- logspace_add(
    log_rho,
    log_one_minus_rho + log_F
  )

  log_L_positive <-
    log_one_minus_rho +
    stats::dnorm(zval, log = TRUE) -
    omega -
    log_Phi_u

  loglik_i <- ifelse(pos, log_L_positive, log_L_zero)
  if (any(!is.finite(loglik_i))) {
    return(invalid("nonfinite_loglikelihood"))
  }

  gamma <- numeric(n)
  gamma[!pos] <- exp(log_rho[!pos] - log_L_zero[!pos])
  gamma <- clamp(gamma, 0, 1)
  posterior_presence <- ifelse(pos, 1, 1 - gamma)

  lambda_u <- stable_mills(u)
  lambda_t <- stable_mills(tval)
  lambda_prime_u <- -lambda_u * (u + lambda_u)
  lambda_prime_t <- -lambda_t * (tval + lambda_t)

  structural_residual <- gamma - rho

  score_eta_scalar <- numeric(n)
  score_omega <- numeric(n)

  score_eta_scalar[pos] <- (zval[pos] + lambda_u[pos]) / sigma
  score_omega[pos] <-
    -1 + zval[pos]^2 + u[pos] * lambda_u[pos]

  A <- (lambda_u - lambda_t) / sigma
  B <- u * lambda_u - tval * lambda_t

  score_eta_scalar[!pos] <- posterior_presence[!pos] * A[!pos]
  score_omega[!pos] <- posterior_presence[!pos] * B[!pos]

  score_alpha <- sweep(X_rho, 1L, structural_residual, FUN = "*")
  score_eta <- sweep(X_eta, 1L, score_eta_scalar, FUN = "*")
  score_i <- cbind(score_alpha, score_eta, logSigma = score_omega)
  colnames(score_i) <- c(
    paste0("rho_", colnames(X_rho)),
    paste0("eta_", colnames(X_eta)),
    "logSigma"
  )

  out <- list(
    valid = TRUE,
    reason = "ok",
    nll = -sum(loglik_i),
    loglik_i = loglik_i,
    gradient_nll = -colSums(score_i),
    score_i = score_i,
    rho = rho,
    eta = eta,
    sigma = sigma,
    gamma = gamma,
    posterior_presence = posterior_presence,
    a = a,
    ell = ell,
    u = u,
    t = tval,
    z = zval,
    lambda_u = lambda_u,
    lambda_t = lambda_t,
    A = A,
    B = B,
    positive = pos
  )

  if (!need_hessian) {
    return(out)
  }

  h_qq <- numeric(n)
  h_qeta <- numeric(n)
  h_qomega <- numeric(n)
  h_etaeta <- numeric(n)
  h_etaomega <- numeric(n)
  h_omegaomega <- numeric(n)

  ## Positive observations.
  h_qq[pos] <- -rho[pos] * (1 - rho[pos])
  h_qeta[pos] <- 0
  h_qomega[pos] <- 0
  h_etaeta[pos] <- -(1 + lambda_prime_u[pos]) / sigma^2
  h_etaomega[pos] <-
    (-2 * zval[pos] - lambda_u[pos] - u[pos] * lambda_prime_u[pos]) /
    sigma
  h_omegaomega[pos] <-
    -2 * zval[pos]^2 -
    u[pos] * lambda_u[pos] -
    u[pos]^2 * lambda_prime_u[pos]

  ## Observed zeros.
  zero <- !pos
  gamma_var <- gamma * (1 - gamma)

  A_eta <- (lambda_prime_t - lambda_prime_u) / sigma^2
  A_omega <-
    (-lambda_u + lambda_t - u * lambda_prime_u +
       tval * lambda_prime_t) /
    sigma
  B_omega <-
    -u * lambda_u - u^2 * lambda_prime_u +
    tval * lambda_t + tval^2 * lambda_prime_t

  h_qq[zero] <-
    gamma_var[zero] - rho[zero] * (1 - rho[zero])
  h_qeta[zero] <- -gamma_var[zero] * A[zero]
  h_qomega[zero] <- -gamma_var[zero] * B[zero]
  h_etaeta[zero] <-
    gamma_var[zero] * A[zero]^2 +
    posterior_presence[zero] * A_eta[zero]
  h_etaomega[zero] <-
    gamma_var[zero] * A[zero] * B[zero] +
    posterior_presence[zero] * A_omega[zero]
  h_omegaomega[zero] <-
    gamma_var[zero] * B[zero]^2 +
    posterior_presence[zero] * B_omega[zero]

  H_aa <- crossprod(
    X_rho,
    sweep(X_rho, 1L, h_qq, FUN = "*")
  )
  H_ae <- crossprod(
    X_rho,
    sweep(X_eta, 1L, h_qeta, FUN = "*")
  )
  H_aw <- crossprod(X_rho, h_qomega)
  H_ee <- crossprod(
    X_eta,
    sweep(X_eta, 1L, h_etaeta, FUN = "*")
  )
  H_ew <- crossprod(X_eta, h_etaomega)
  H_ww <- sum(h_omegaomega)

  H <- rbind(
    cbind(H_aa, H_ae, H_aw),
    cbind(t(H_ae), H_ee, H_ew),
    cbind(t(H_aw), t(H_ew), matrix(H_ww, 1L, 1L))
  )
  H <- (H + t(H)) / 2
  dimnames(H) <- list(colnames(score_i), colnames(score_i))

  out$hessian_loglik <- H
  out$h_qq <- h_qq
  out$h_qeta <- h_qeta
  out$h_qomega <- h_qomega
  out$h_etaeta <- h_etaeta
  out$h_etaomega <- h_etaomega
  out$h_omegaomega <- h_omegaomega
  out
}

#' Fit a taxon-specific depth-aware zero model
#' @keywords internal
#' @noRd
fit_zero_model <- function(
    y,
    D,
    z = NULL,
    g = NULL,
    d0 = 1,
    include_group_rho = FALSE,
    include_group_eta = FALSE,
    control = list()) {

  control <- validate_fit_control(control)

  y <- as.numeric(y)
  D <- as.numeric(D)
  n <- length(y)

  make_zero_design <- function(include_group) {
    if (include_group) {
      design_matrix(g = g, z = z, include_group = TRUE)
    } else {
      X <- matrix(1, nrow = n, ncol = 1L)
      colnames(X) <- "Intercept"
      if (!is.null(z)) {
        X <- cbind(X, z)
      }
      X
    }
  }

  X_rho <- make_zero_design(include_group_rho)
  X_eta <- make_zero_design(include_group_eta)
  p_rho <- ncol(X_rho)
  p_eta <- ncol(X_eta)
  npar <- p_rho + p_eta + 1L

  empty_result <- function(reason) {
    list(
      valid = FALSE,
      reason = reason,
      par = rep(NA_real_, npar),
      alpha = rep(NA_real_, p_rho),
      eta_coef = rep(NA_real_, p_eta),
      sigma = NA_real_,
      X_rho = X_rho,
      X_eta = X_eta,
      state = NULL,
      information = NULL,
      information_inverse = NULL,
      optimization = list(
        method = NA_character_,
        convergence = NA_integer_,
        objective = NA_real_,
        n_starts = 0L,
        n_converged = 0L,
        used_fallback = FALSE,
        used_expanded_bounds = FALSE,
        boundary = NA,
        gradient_mean_norm = NA_real_,
        information_condition = Inf
      )
    )
  }

  if (length(D) != n || any(!is.finite(D)) || any(D <= d0)) {
    return(empty_result("invalid_library_size_or_detection_boundary"))
  }
  if (sum(y > 0) < 2L) {
    return(empty_result("fewer_than_two_positive_observations"))
  }
  if (!is_full_rank(X_rho) || !is_full_rank(X_eta)) {
    return(empty_result("rank_deficient_zero_model_design"))
  }
  if (n <= npar + 1L) {
    return(empty_result("insufficient_samples_for_zero_model"))
  }

  zero_fraction <- mean(y == 0)
  init_eta <- mean(log(y[y > 0] / D[y > 0]))
  if (!is.finite(init_eta)) {
    init_eta <- -5
  }

  ## Coefficient-wise numerical boxes. The present-state intercept is separated
  ## from its slopes to allow large group effects without forcing the intercept
  ## to an implausibly broad range.
  lower <- c(
    rep(-15, p_rho),
    -30,
    if (p_eta > 1L) rep(-20, p_eta - 1L) else numeric(0),
    log(0.05)
  )
  upper <- c(
    rep(15, p_rho),
    10,
    if (p_eta > 1L) rep(20, p_eta - 1L) else numeric(0),
    log(10)
  )

  lower_expanded <- c(
    rep(-25, p_rho),
    -50,
    if (p_eta > 1L) rep(-30, p_eta - 1L) else numeric(0),
    log(0.02)
  )
  upper_expanded <- c(
    rep(25, p_rho),
    20,
    if (p_eta > 1L) rep(30, p_eta - 1L) else numeric(0),
    log(20)
  )

  start_spec <- rbind(
    c(structural_fraction = 0.50, sigma = 1.0),
    c(structural_fraction = 0.25, sigma = 0.5),
    c(structural_fraction = 0.25, sigma = 2.0),
    c(structural_fraction = 0.75, sigma = 0.5),
    c(structural_fraction = 0.75, sigma = 2.0)
  )

  make_start <- function(structural_fraction, sigma_start, lo, hi) {
    p_start <- clamp(structural_fraction * zero_fraction, 0.02, 0.90)
    start <- c(
      stats::qlogis(p_start),
      rep(0, p_rho - 1L),
      init_eta,
      rep(0, p_eta - 1L),
      log(sigma_start)
    )
    pmin(pmax(start, lo), hi)
  }

  starts <- lapply(seq_len(nrow(start_spec)), function(k) {
    make_start(
      structural_fraction = start_spec[k, "structural_fraction"],
      sigma_start = start_spec[k, "sigma"],
      lo = lower,
      hi = upper
    )
  })
  start_keys <- vapply(
    starts,
    function(x) paste(format(x, digits = 16), collapse = "|"),
    character(1)
  )
  starts <- starts[!duplicated(start_keys)]

  make_objective <- function(lo, hi) {
    cache_par <- NULL
    cache_state <- NULL

    evaluate <- function(par) {
      if (!is.null(cache_par) && identical(par, cache_par)) {
        return(cache_state)
      }
      cache_par <<- par
      cache_state <<- zero_model_state(
        par = par,
        y = y,
        D = D,
        X_rho = X_rho,
        X_eta = X_eta,
        d0 = d0,
        need_hessian = FALSE
      )
      cache_state
    }

    list(
      fn = function(par) {
        st <- evaluate(par)
        if (!isTRUE(st$valid)) 1e100 else st$nll
      },
      gr = function(par) {
        st <- evaluate(par)
        if (!isTRUE(st$valid)) rep(0, length(par)) else st$gradient_nll
      },
      boxed_fn = function(par) {
        if (length(par) != length(lo) || any(!is.finite(par))) {
          return(1e100)
        }
        below <- pmax(lo - par, 0)
        above <- pmax(par - hi, 0)
        if (any(below > 0) || any(above > 0)) {
          return(1e100 + 1e6 * sum(below^2 + above^2))
        }
        st <- evaluate(par)
        if (!isTRUE(st$valid)) 1e100 else st$nll
      }
    )
  }

  run_lbfgsb <- function(start, lo, hi, objective) {
    tryCatch(
      stats::optim(
        par = pmin(pmax(start, lo), hi),
        fn = objective$fn,
        gr = objective$gr,
        method = "L-BFGS-B",
        lower = lo,
        upper = hi,
        control = list(
          maxit = as.integer(control$maxit),
          factr = control$factr
        )
      ),
      error = function(e) NULL
    )
  }

  finite_fit <- function(fit) {
    !is.null(fit) && is.finite(fit$value) && all(is.finite(fit$par))
  }
  converged_fit <- function(fit) {
    finite_fit(fit) && isTRUE(fit$convergence == 0L)
  }
  choose_best <- function(fits, predicate = finite_fit) {
    ok <- vapply(fits, predicate, logical(1))
    if (!any(ok)) {
      return(NULL)
    }
    fits_ok <- fits[ok]
    values <- vapply(fits_ok, function(x) x$value, numeric(1))
    fits_ok[[which.min(values)]]
  }

  objective_primary <- make_objective(lower, upper)
  fits_primary <- lapply(
    starts,
    run_lbfgsb,
    lo = lower,
    hi = upper,
    objective = objective_primary
  )

  used_fallback <- FALSE
  best_converged <- choose_best(fits_primary, converged_fit)
  best_finite <- choose_best(fits_primary, finite_fit)

  ## If bounded L-BFGS-B does not converge, use Nelder-Mead only to obtain a
  ## candidate, then polish and validate it with bounded L-BFGS-B.
  if (is.null(best_converged)) {
    used_fallback <- TRUE
    nm_start <- if (!is.null(best_finite)) best_finite$par else starts[[1L]]

    nm_fit <- tryCatch(
      stats::optim(
        par = nm_start,
        fn = objective_primary$boxed_fn,
        method = "Nelder-Mead",
        control = list(
          maxit = as.integer(control$nm_maxit),
          reltol = control$nm_reltol
        )
      ),
      error = function(e) NULL
    )

    if (finite_fit(nm_fit)) {
      polished <- run_lbfgsb(
        nm_fit$par,
        lo = lower,
        hi = upper,
        objective = objective_primary
      )
      fits_primary <- c(fits_primary, list(polished))
      best_converged <- choose_best(fits_primary, converged_fit)
      best_finite <- choose_best(fits_primary, finite_fit)
    }
  }

  validate_candidate <- function(fit, lo, hi, method_name, expanded) {
    if (!converged_fit(fit)) {
      return(list(valid = FALSE, reason = "no_converged_fit"))
    }

    st <- zero_model_state(
      par = fit$par,
      y = y,
      D = D,
      X_rho = X_rho,
      X_eta = X_eta,
      d0 = d0,
      need_hessian = TRUE
    )
    if (!isTRUE(st$valid)) {
      return(list(valid = FALSE, reason = st$reason))
    }

    distance_to_box <- pmin(fit$par - lo, hi - fit$par)
    scale_box <- 1 + abs(fit$par)
    boundary <- any(distance_to_box <= control$boundary_tol * scale_box)
    gradient_mean_norm <- max(abs(st$gradient_nll)) / n

    info_result <- invert_spd(
      -st$hessian_loglik,
      rel_tol = control$information_rel_tol
    )

    valid <-
      !boundary &&
      is.finite(gradient_mean_norm) &&
      gradient_mean_norm <= control$gradient_mean_tol &&
      isTRUE(info_result$valid)

    reason <- if (valid) {
      "ok"
    } else if (boundary) {
      "boundary_solution"
    } else if (!is.finite(gradient_mean_norm) ||
               gradient_mean_norm > control$gradient_mean_tol) {
      "gradient_not_small"
    } else {
      "ill_conditioned_information"
    }

    list(
      valid = valid,
      reason = reason,
      fit = fit,
      state = st,
      information = -st$hessian_loglik,
      information_inverse = info_result$inverse,
      information_condition = info_result$condition,
      boundary = boundary,
      gradient_mean_norm = gradient_mean_norm,
      method = method_name,
      expanded = expanded
    )
  }

  assessed <- validate_candidate(
    best_converged,
    lower,
    upper,
    method_name = if (used_fallback) "L-BFGS-B after Nelder-Mead" else "L-BFGS-B",
    expanded = FALSE
  )

  ## A boundary or otherwise inadequate primary fit is retried under a broader
  ## box. Formal inference is still withheld unless the expanded fit is an
  ## interior, converged, well-conditioned solution.
  if (!isTRUE(assessed$valid) && isTRUE(control$expand_bounds)) {
    objective_expanded <- make_objective(lower_expanded, upper_expanded)

    expanded_starts <- list(
      if (!is.null(best_converged)) best_converged$par else
        if (!is.null(best_finite)) best_finite$par else starts[[1L]],
      make_start(0.25, 1.0, lower_expanded, upper_expanded),
      make_start(0.75, 2.0, lower_expanded, upper_expanded)
    )

    expanded_fits <- lapply(
      expanded_starts,
      run_lbfgsb,
      lo = lower_expanded,
      hi = upper_expanded,
      objective = objective_expanded
    )
    best_expanded <- choose_best(expanded_fits, converged_fit)

    assessed_expanded <- validate_candidate(
      best_expanded,
      lower_expanded,
      upper_expanded,
      method_name = "expanded-bound L-BFGS-B",
      expanded = TRUE
    )

    if (isTRUE(assessed_expanded$valid)) {
      assessed <- assessed_expanded
    }
  }

  if (!isTRUE(assessed$valid)) {
    candidate <- if (!is.null(best_converged)) best_converged else best_finite
    out <- empty_result(assessed$reason %||% "zero_model_fit_failed")
    if (!is.null(candidate)) {
      out$par <- candidate$par
      out$alpha <- candidate$par[seq_len(p_rho)]
      out$eta_coef <- candidate$par[p_rho + seq_len(p_eta)]
      out$sigma <- exp(candidate$par[p_rho + p_eta + 1L])
      out$optimization$method <- if (used_fallback) {
        "non-inferential fallback candidate"
      } else {
        "non-inferential L-BFGS-B candidate"
      }
      out$optimization$convergence <- candidate$convergence
      out$optimization$objective <- candidate$value
    }
    out$optimization$n_starts <- length(starts)
    out$optimization$n_converged <- sum(vapply(fits_primary, converged_fit, logical(1)))
    out$optimization$used_fallback <- used_fallback
    out$optimization$used_expanded_bounds <- FALSE
    out$optimization$boundary <- assessed$boundary %||% NA
    out$optimization$gradient_mean_norm <- assessed$gradient_mean_norm %||% NA_real_
    out$optimization$information_condition <- assessed$information_condition %||% Inf
    return(out)
  }

  par <- assessed$fit$par
  names(par) <- colnames(assessed$state$score_i)

  alpha <- par[seq_len(p_rho)]
  eta_coef <- par[p_rho + seq_len(p_eta)]
  names(alpha) <- colnames(X_rho)
  names(eta_coef) <- colnames(X_eta)

  list(
    valid = TRUE,
    reason = "ok",
    par = par,
    alpha = alpha,
    eta_coef = eta_coef,
    sigma = assessed$state$sigma,
    X_rho = X_rho,
    X_eta = X_eta,
    state = assessed$state,
    information = assessed$information,
    information_inverse = assessed$information_inverse,
    optimization = list(
      method = assessed$method,
      convergence = assessed$fit$convergence,
      objective = assessed$fit$value,
      n_starts = length(starts),
      n_converged = sum(vapply(fits_primary, converged_fit, logical(1))),
      used_fallback = used_fallback,
      used_expanded_bounds = assessed$expanded,
      boundary = assessed$boundary,
      gradient_mean_norm = assessed$gradient_mean_norm,
      information_condition = assessed$information_condition
    )
  )
}

#' Nuisance-adjusted structural-absence score test
#' @keywords internal
#' @noRd
structural_absence_test <- function(
    y,
    D,
    g,
    z = NULL,
    d0 = 1,
    fit_control = list()) {

  n <- length(y)
  fit <- fit_zero_model(
    y = y,
    D = D,
    z = z,
    g = g,
    d0 = d0,
    include_group_rho = FALSE,
    include_group_eta = TRUE,
    control = fit_control
  )

  invalid <- function(reason = fit$reason) {
    list(
      formed = FALSE,
      p = NA_real_,
      z = NA_real_,
      U = NA_real_,
      V = NA_real_,
      gamma = rep(NA_real_, n),
      rho = rep(NA_real_, n),
      adjusted_score = rep(NA_real_, n),
      fit = fit,
      reason = reason
    )
  }

  if (!isTRUE(fit$valid)) {
    return(invalid())
  }

  st <- fit$state
  s_delta <- as.numeric(g) * (st$gamma - st$rho)

  ## Derivative of the structural group score with respect to every nuisance
  ## parameter. The scalar derivatives h_qq, h_qeta, and h_qomega are supplied
  ## by the analytic observed-data Hessian.
  d_delta <- cbind(
    sweep(fit$X_rho, 1L, as.numeric(g) * st$h_qq, FUN = "*"),
    sweep(fit$X_eta, 1L, as.numeric(g) * st$h_qeta, FUN = "*"),
    logSigma = as.numeric(g) * st$h_qomega
  )
  colnames(d_delta) <- colnames(st$score_i)

  ## If H is the nuisance log-likelihood Hessian and I = -H, then
  ## A_delta,theta A_theta,theta^{-1} equals
  ## -sum(d s_delta / d theta) I^{-1}.
  projection <- -as.numeric(
    matrix(colSums(d_delta), nrow = 1L) %*%
      fit$information_inverse
  )

  adjusted_score <- s_delta - as.numeric(st$score_i %*% projection)
  U <- sum(adjusted_score)
  V <- sum(adjusted_score^2)

  if (!is.finite(V) || V <= 0 || !is.finite(U)) {
    return(invalid("nonpositive_adjusted_score_variance"))
  }

  z_stat <- U / sqrt(V)
  p_value <- stats::pchisq(z_stat^2, df = 1, lower.tail = FALSE)

  list(
    formed = TRUE,
    p = clamp(p_value, 0, 1),
    z = z_stat,
    U = U,
    V = V,
    gamma = st$gamma,
    rho = st$rho,
    adjusted_score = adjusted_score,
    projection = projection,
    nuisance_score_balance = max(abs(colSums(st$score_i))) / n,
    fit = fit,
    reason = "ok"
  )
}

#' Auxiliary full zero-model fit for abundance weights and censoring severity
#' @keywords internal
#' @noRd
auxiliary_zero_fit <- function(
    y,
    D,
    g,
    z = NULL,
    d0 = 1,
    fit_control = list()) {

  n <- length(y)
  fit <- fit_zero_model(
    y = y,
    D = D,
    z = z,
    g = g,
    d0 = d0,
    include_group_rho = TRUE,
    include_group_eta = TRUE,
    control = fit_control
  )

  invalid <- function(reason = fit$reason) {
    list(
      valid = FALSE,
      gamma = rep(NA_real_, n),
      censor_correction = rep(NA_real_, n),
      structural_group_log_odds = NA_real_,
      structural_group_se = NA_real_,
      latent_group_coefficient = NA_real_,
      latent_group_se = NA_real_,
      fit = fit,
      reason = reason
    )
  }

  if (!isTRUE(fit$valid)) {
    return(invalid())
  }

  st <- fit$state
  censor_correction <- ifelse(
    y > 0,
    0,
    st$sigma * st$lambda_t
  )
  censor_correction[!is.finite(censor_correction)] <- NA_real_

  rho_group_local <- match("Group", colnames(fit$X_rho))
  eta_group_local <- match("Group", colnames(fit$X_eta))

  p_rho <- ncol(fit$X_rho)
  rho_group_global <- rho_group_local
  eta_group_global <- p_rho + eta_group_local

  delta_hat <- if (is.finite(rho_group_local)) {
    fit$alpha[rho_group_local]
  } else {
    NA_real_
  }
  delta_se <- if (is.finite(rho_group_global)) {
    sqrt(fit$information_inverse[rho_group_global, rho_group_global])
  } else {
    NA_real_
  }

  eta_group_hat <- if (is.finite(eta_group_local)) {
    fit$eta_coef[eta_group_local]
  } else {
    NA_real_
  }
  eta_group_se <- if (is.finite(eta_group_global)) {
    sqrt(fit$information_inverse[eta_group_global, eta_group_global])
  } else {
    NA_real_
  }

  list(
    valid = TRUE,
    gamma = st$gamma,
    censor_correction = censor_correction,
    structural_group_log_odds = unname(delta_hat),
    structural_group_se = unname(delta_se),
    latent_group_coefficient = unname(eta_group_hat),
    latent_group_se = unname(eta_group_se),
    fit = fit,
    reason = "ok"
  )
}


# Posterior-presence-weighted HRIC abundance component

#' Fit the affine closure background on a prespecified training set
#' @keywords internal
#' @noRd
fit_affine_background <- function(
    beta_raw,
    v_ref,
    train_idx,
    lts_seed = 1907L,
    min_train = 5L,
    inlier_cutoff = 2.5,
    matrix_tol = 1e-12) {

  invalid <- function(reason) {
    list(
      valid = FALSE,
      reason = reason,
      H = integer(0),
      GH = NULL,
      BH = NULL,
      coefficients = c(NA_real_, NA_real_),
      training_n = 0L,
      inlier_n = 0L
    )
  }

  train_idx <- intersect(
    as.integer(train_idx),
    which(is.finite(beta_raw) & is.finite(v_ref))
  )

  if (length(train_idx) < min_train) {
    return(invalid("insufficient_affine_training_taxa"))
  }
  if (stats::sd(v_ref[train_idx]) <= 1e-10) {
    return(invalid("negligible_affine_reference_variation"))
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for the LTS affine selection.")
  }

  x <- v_ref[train_idx]
  y <- beta_raw[train_idx]
  lts_quantile <- floor(length(y) / 2) + 1L

  initial_fit <- with_preserved_seed(
    lts_seed,
    tryCatch(
      MASS::lqs(
        y ~ x,
        method = "lts",
        quantile = lts_quantile,
        control = list(
          psamp = NA,
          nsamp = "best",
          adjust = TRUE
        )
      ),
      error = function(e) NULL
    )
  )

  if (is.null(initial_fit)) {
    return(invalid("lts_fit_failed"))
  }

  initial_coef <- stats::coef(initial_fit)
  if (length(initial_coef) != 2L || any(!is.finite(initial_coef))) {
    return(invalid("nonfinite_lts_coefficients"))
  }

  initial_residual <- y - (initial_coef[1L] + initial_coef[2L] * x)
  robust_scale <- stats::mad(
    initial_residual,
    center = stats::median(initial_residual),
    constant = 1.4826
  )

  if (!is.finite(robust_scale) || robust_scale <= 0) {
    robust_scale <- stats::sd(initial_residual)
  }
  if (!is.finite(robust_scale) || robust_scale <= 0) {
    return(invalid("nonpositive_affine_residual_scale"))
  }

  inlier_local <- abs(initial_residual) <= inlier_cutoff * robust_scale
  H <- train_idx[inlier_local]

  if (length(H) < 4L) {
    return(invalid("too_few_affine_inliers"))
  }

  GH <- cbind(Intercept = 1, v_ref = v_ref[H])
  if (!is_full_rank(GH)) {
    return(invalid("rank_deficient_affine_refit"))
  }

  Sgg <- crossprod(GH)
  if (!is.finite(rcond(Sgg)) || rcond(Sgg) <= matrix_tol) {
    return(invalid("ill_conditioned_affine_refit"))
  }

  Sgg_inv <- tryCatch(solve(Sgg), error = function(e) NULL)
  if (is.null(Sgg_inv) || any(!is.finite(Sgg_inv))) {
    return(invalid("affine_refit_inversion_failed"))
  }

  BH <- Sgg_inv %*% t(GH)
  coefficients <- as.numeric(BH %*% beta_raw[H])

  list(
    valid = TRUE,
    reason = "ok",
    H = H,
    GH = GH,
    BH = BH,
    coefficients = coefficients,
    training_n = length(train_idx),
    inlier_n = length(H),
    robust_scale = robust_scale,
    lts_quantile = lts_quantile
  )
}

#' Deterministic cross-fitting partitions for the affine background
#' @keywords internal
#' @noRd
make_affine_partitions <- function(ok_idx, requested_folds = 5L) {
  ok_idx <- sort(unique(as.integer(ok_idx)))
  m <- length(ok_idx)

  if (m < 6L) {
    return(list())
  }

  ## For small eligible sets, leave-one-taxon-out gives at least five
  ## training taxa. Otherwise, deterministic K-fold cross-fitting is used.
  if (m < 10L) {
    return(lapply(seq_along(ok_idx), function(k) {
      list(
        id = k,
        targets = ok_idx[k],
        training = ok_idx[-k],
        scheme = "leave-one-taxon-out"
      )
    }))
  }

  K <- max(2L, min(as.integer(requested_folds), m))
  fold_id <- rep(seq_len(K), length.out = m)

  lapply(seq_len(K), function(k) {
    list(
      id = k,
      targets = ok_idx[fold_id == k],
      training = ok_idx[fold_id != k],
      scheme = paste0(K, "-fold cross-fit")
    )
  })
}

#' Posterior-presence-weighted HRIC abundance test
#' @keywords internal
#' @noRd
hric_wls_p <- function(
    Y,
    D,
    g,
    W,
    z = NULL,
    m_abund = 5L,
    abundance_eligibility = c("pooled", "per_group", "none"),
    Ccor = NULL,
    affine_folds = 5L,
    lts_seed = 1907L,
    matrix_tol = 1e-12) {

  Y <- as.matrix(Y)
  W <- as.matrix(W)
  n <- nrow(Y)
  J <- ncol(Y)
  N_retained <- rowSums(Y)
  abundance_eligibility <- match.arg(abundance_eligibility)

  if (!all(dim(W) == dim(Y))) {
    stop("`W` must have the same dimensions as `Y`.")
  }
  if (length(D) != n || any(!is.finite(D)) || any(D <= 0)) {
    stop("`D` must contain one positive finite library size per sample.")
  }

  if (!is.null(Ccor)) {
    Ccor <- as.matrix(Ccor)
    if (!all(dim(Ccor) == dim(Y))) {
      stop("`Ccor` must have the same dimensions as `Y`.")
    }
  }

  M <- hric_transform(Y)

  Z <- z
  raw_logD <- log(D)
  sd_logD <- stats::sd(raw_logD)

  if (is.finite(sd_logD) && sd_logD > 1e-12) {
    logD <- (raw_logD - mean(raw_logD)) / sd_logD
    Z_try <- if (is.null(Z)) {
      matrix(logD, ncol = 1L, dimnames = list(NULL, "logDepth"))
    } else {
      cbind(Z, logDepth = logD)
    }

    X_try <- design_matrix(g = g, z = Z_try, include_group = TRUE)
    if (is_full_rank(X_try)) {
      Z <- Z_try
    } else {
      warning(
        "The log-depth column is collinear with the supplied design and was omitted.",
        call. = FALSE
      )
    }
  }

  X_base <- design_matrix(g = g, z = Z, include_group = TRUE)
  if (!is_full_rank(X_base)) {
    stop("The abundance base design is rank deficient.")
  }

  beta_raw <- rep(NA_real_, J)
  se_raw_hc3 <- rep(NA_real_, J)
  Phi_raw <- matrix(
    NA_real_,
    nrow = n,
    ncol = J,
    dimnames = dimnames(Y)
  )
  n_pos0 <- colSums(Y[g == 0, , drop = FALSE] > 0)
  n_pos1 <- colSums(Y[g == 1, , drop = FALSE] > 0)
  design_condition <- rep(NA_real_, J)
  censor_column_used <- rep(FALSE, J)

  for (j in seq_len(J)) {
    count_eligible <- switch(
      abundance_eligibility,
      pooled = (n_pos0[j] + n_pos1[j]) >= 2L * m_abund,
      per_group = n_pos0[j] >= m_abund && n_pos1[j] >= m_abund,
      none = TRUE
    )
    if (!isTRUE(count_eligible)) {
      next
    }

    w <- as.numeric(W[, j])
    if (any(!is.finite(w))) {
      next
    }
    w <- clamp(w, 0, 1)

    X <- X_base

    if (!is.null(Ccor)) {
      c_j <- as.numeric(Ccor[, j])
      if (all(is.finite(c_j))) {
        use_w <- w > 1e-8
        if (sum(use_w) > ncol(X_base) + 2L) {
          c_bar <- stats::weighted.mean(c_j[use_w], w[use_w])
          c_wsd <- sqrt(stats::weighted.mean(
            (c_j[use_w] - c_bar)^2,
            w[use_w]
          ))

          if (is.finite(c_wsd) && c_wsd > 1e-12) {
            X_try <- cbind(X_base, censorSeverity = c_j)
            XtWX_try <- crossprod(X_try * sqrt(w))
            if (
              is_full_rank(X_try) &&
              is.finite(rcond(XtWX_try)) &&
              rcond(XtWX_try) > matrix_tol
            ) {
              X <- X_try
              censor_column_used[j] <- TRUE
            }
          }
        }
      }
    }

    group_index <- match("Group", colnames(X))
    if (!is.finite(group_index)) {
      next
    }
    if (sum(w > 1e-8) <= ncol(X) + 1L) {
      next
    }

    sqrt_w <- sqrt(w)
    Xs <- X * sqrt_w
    ys <- M[, j] * sqrt_w
    XtWX <- crossprod(Xs)

    rc <- rcond(XtWX)
    if (!is.finite(rc) || rc <= matrix_tol) {
      next
    }
    design_condition[j] <- 1 / rc

    XtWX_inv <- tryCatch(solve(XtWX), error = function(e) NULL)
    if (is.null(XtWX_inv) || any(!is.finite(XtWX_inv))) {
      next
    }

    bhat <- as.numeric(XtWX_inv %*% crossprod(Xs, ys))
    residual <- as.numeric(M[, j] - X %*% bhat)

    leverage <- rowSums((Xs %*% XtWX_inv) * Xs)
    leverage <- clamp(leverage, 0, 1 - 1e-8)

    group_loading <- as.numeric(X %*% XtWX_inv[, group_index, drop = FALSE])
    phi_j <- group_loading * w * residual / (1 - leverage)

    if (any(!is.finite(phi_j))) {
      next
    }

    var_raw <- sum(phi_j^2)
    if (!is.finite(var_raw) || var_raw <= 0) {
      next
    }

    beta_raw[j] <- bhat[group_index]
    se_raw_hc3[j] <- sqrt(var_raw)
    Phi_raw[, j] <- phi_j
  }

  ## Equal weighting of the two groups makes the reference invariant to swapping
  ## the group labels, even when group sizes differ.
  P_retained <- sweep(Y, 1L, N_retained, FUN = "/")
  mean0 <- colMeans(P_retained[g == 0, , drop = FALSE])
  mean1 <- colMeans(P_retained[g == 1, , drop = FALSE])
  v_ref <- sqrt(pmax(0, 0.5 * (mean0 + mean1)))

  ok_idx <- which(
    is.finite(beta_raw) &
      is.finite(v_ref) &
      colSums(is.finite(Phi_raw)) == n
  )
  partitions <- make_affine_partitions(ok_idx, requested_folds = affine_folds)

  p_value <- rep(NA_real_, J)
  tested <- rep(FALSE, J)
  beta_corrected <- rep(NA_real_, J)
  se_joint <- rep(NA_real_, J)
  z_stat <- rep(NA_real_, J)
  background_intercept <- rep(NA_real_, J)
  background_slope <- rep(NA_real_, J)
  affine_partition <- rep(NA_integer_, J)
  affine_training_n <- rep(NA_integer_, J)
  affine_inlier_n <- rep(NA_integer_, J)
  affine_scheme <- rep(NA_character_, J)

  for (part in partitions) {
    lf <- fit_affine_background(
      beta_raw = beta_raw,
      v_ref = v_ref,
      train_idx = part$training,
      lts_seed = as.integer(lts_seed + part$id),
      matrix_tol = matrix_tol
    )

    if (!isTRUE(lf$valid)) {
      next
    }

    phi_H <- Phi_raw[, lf$H, drop = FALSE]
    if (any(!is.finite(phi_H))) {
      next
    }

    for (j in part$targets) {
      if (
        !is.finite(beta_raw[j]) ||
          !is.finite(v_ref[j]) ||
          any(!is.finite(Phi_raw[, j]))
      ) {
        next
      }

      g_j <- c(1, v_ref[j])
      fitted_background <- sum(g_j * lf$coefficients)
      beta_corr_j <- beta_raw[j] - fitted_background

      affine_loading <- as.numeric(t(g_j) %*% lf$BH)
      phi_background <- as.numeric(phi_H %*% affine_loading)
      phi_corrected <- Phi_raw[, j] - phi_background
      var_joint <- sum(phi_corrected^2)

      if (!is.finite(var_joint) || var_joint <= 0) {
        next
      }

      se_j <- sqrt(var_joint)
      z_j <- beta_corr_j / se_j
      p_j <- 2 * stats::pnorm(-abs(z_j))

      beta_corrected[j] <- beta_corr_j
      se_joint[j] <- se_j
      z_stat[j] <- z_j
      p_value[j] <- clamp(p_j, 0, 1)
      tested[j] <- TRUE
      background_intercept[j] <- lf$coefficients[1L]
      background_slope[j] <- lf$coefficients[2L]
      affine_partition[j] <- part$id
      affine_training_n[j] <- lf$training_n
      affine_inlier_n[j] <- lf$inlier_n
      affine_scheme[j] <- part$scheme
    }
  }

  names(p_value) <- colnames(Y)
  names(tested) <- colnames(Y)
  names(beta_raw) <- colnames(Y)
  names(beta_corrected) <- colnames(Y)
  names(se_joint) <- colnames(Y)
  names(z_stat) <- colnames(Y)

  list(
    p = p_value,
    tested = tested,
    beta_raw = beta_raw,
    se_raw_hc3 = se_raw_hc3,
    beta_corrected = beta_corrected,
    se_joint = se_joint,
    z = z_stat,
    v_reference = v_ref,
    background_intercept = background_intercept,
    background_slope = background_slope,
    affine_partition = affine_partition,
    affine_training_n = affine_training_n,
    affine_inlier_n = affine_inlier_n,
    affine_scheme = affine_scheme,
    n_positive_group0 = n_pos0,
    n_positive_group1 = n_pos1,
    censor_column_used = censor_column_used,
    design_condition = design_condition,
    abundance_eligibility = abundance_eligibility
  )
}


# Main analysis

#' Run DASH differential-abundance analysis
#'
#' For each retained taxon, `dash()` fits a depth-aware structural-absence
#' working model, forms a nuisance-adjusted structural-absence score test, and
#' estimates a posterior-presence-weighted HRIC abundance contrast with a
#' cross-fitted affine closure correction. The primary omnibus p-value is a
#' fixed two-component Bonferroni min-P combination.
#'
#' @details
#' The zero model uses `library_size` as the pre-filter sequencing-depth
#' quantity. HRIC is computed from the composition obtained after DASH taxon
#' retention. A component that cannot be estimated is returned as `NA`; for the
#' primary omnibus test, that component is assigned p = 1 while the factor of
#' two is retained. The Cauchy combination is returned only when both component
#' p-values are available.
#'
#' Structural absence is a model-implied latent state. The positive-state model
#' is an upper-truncated-normal working likelihood for transformed integer
#' counts rather than an exact likelihood for the original count process. The
#' abundance HC3 standard error is conditional on fitted auxiliary quantities
#' and the cross-fitted affine background. Use [dash_bootstrap()] to assess the
#' effect of refitting the complete pipeline.
#'
#' @param counts A non-negative integer-valued matrix or data frame with samples
#'   in rows and taxa in columns.
#' @param group A binary group indicator with one value per sample.
#' @param covariates An optional vector, matrix, or data frame of adjustment
#'   covariates. Character columns are treated as factors, factors use treatment
#'   contrasts, binary numeric columns are mapped to 0/1, and other numeric
#'   columns are centered and scaled.
#' @param library_size An optional positive vector containing the sequencing
#'   depth before DASH taxon filtering. When omitted, row sums of the supplied
#'   count matrix before DASH filtering are used.
#' @param case_level An optional observed group value to code as 1. When omitted,
#'   the second observed factor level or larger numeric value is coded as 1.
#' @param d0 A positive detection constant used in the depth-specific boundary
#'   `log(d0 / library_size)`. The default is 1.
#' @param m_abund A positive integer controlling abundance-component sparsity
#'   eligibility. With `abundance_eligibility = "pooled"`, at least
#'   `2 * m_abund` observed positives are required across all samples.
#' @param abundance_eligibility One of `"pooled"`, `"per_group"`, or `"none"`.
#'   The default `"pooled"` applies a group-blind pooled positive-count filter;
#'   `"per_group"` requires at least `m_abund` positives in each group; and
#'   `"none"` relies only on downstream estimability checks.
#' @param min_positive A non-negative integer giving the minimum number of
#'   observed positives required to retain a taxon.
#' @param affine_folds An integer of at least 2 giving the requested number of
#'   deterministic cross-fitting folds when at least ten taxa have estimable raw
#'   abundance coefficients. Smaller eligible sets use leave-one-taxon-out
#'   cross-fitting.
#' @param lts_seed A non-negative integer seed used internally for deterministic
#'   least-trimmed-squares selection. The caller's random-number state is
#'   restored after each internal fit.
#' @param fit_control An optional named list overriding zero-model optimization
#'   controls. Supported elements are `maxit`, `factr`, `nm_maxit`, `nm_reltol`,
#'   `gradient_mean_tol`, `boundary_tol`, `information_rel_tol`, and
#'   `expand_bounds`.
#' @param return_posteriors Logical. When `TRUE`, fitted posterior matrices and
#'   retained analysis inputs are attached as a `posteriors` attribute.
#'
#' @return A data frame with one row per retained taxon. Primary inferential
#'   columns are `p_absence`, `p_abundance`, `p_bonf_minp`, and `p_cauchy`.
#'   Effect and uncertainty columns include `z_absence`,
#'   `beta_abundance_raw`, `beta_abundance_corrected`,
#'   `se_abundance_raw_hc3`, `se_abundance_joint_hc3`, and `z_abundance`.
#'   Additional columns report component availability, affine-background
#'   diagnostics, positive-count summaries, and zero-model fit diagnostics.
#'   Attributes record group coding, retained sample and taxon indices,
#'   sequencing depths, retained-table totals, processed covariates, and the
#'   inferential scope of the working models.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 24
#' p <- 8
#' group <- rep(0:1, each = n / 2)
#' counts <- matrix(stats::rpois(n * p, lambda = 3), nrow = n, ncol = p)
#' counts[matrix(stats::runif(n * p) < 0.35, nrow = n)] <- 0
#' counts[rowSums(counts) == 0, 1] <- 1
#' colnames(counts) <- paste0("taxon_", seq_len(p))
#' depth <- rowSums(counts) + 1000
#'
#' result <- dash(
#'   counts = counts,
#'   group = group,
#'   library_size = depth,
#'   min_positive = 2,
#'   m_abund = 2,
#'   affine_folds = 3
#' )
#' result$q_bonf_minp <- stats::p.adjust(result$p_bonf_minp, method = "BH")
#' head(result[c("taxon", "p_absence", "p_abundance", "p_bonf_minp")])
#' }
#'
#' @seealso [hric_transform()], [dash_bootstrap()],
#'   [dash_label_swap_check()], [stats::p.adjust()]
#' @family DASH functions
#' @export
#' @md
dash <- function(
    counts,
    group,
    covariates = NULL,
    library_size = NULL,
    case_level = NULL,
    d0 = 1,
    m_abund = 5L,
    abundance_eligibility = c("pooled", "per_group", "none"),
    min_positive = 3L,
    affine_folds = 5L,
    lts_seed = 1907L,
    fit_control = list(),
    return_posteriors = FALSE) {

  Y_full <- as.matrix(counts)
  if (!is.numeric(Y_full)) {
    stop("`counts` must be a numeric matrix or data frame.")
  }
  if (anyNA(Y_full) || any(!is.finite(Y_full))) {
    stop("`counts` contains missing or non-finite values.")
  }
  if (any(Y_full < 0)) {
    stop("`counts` must be non-negative.")
  }
  if (any(abs(Y_full - round(Y_full)) > 1e-8)) {
    stop("`counts` must contain integer-valued counts.")
  }
  storage.mode(Y_full) <- "double"

  n_original <- nrow(Y_full)
  J_original <- ncol(Y_full)
  if (n_original < 3L || J_original < 2L) {
    stop("DASH requires at least three samples and two taxa.")
  }

  if (is.null(colnames(Y_full))) {
    colnames(Y_full) <- paste0("tax", seq_len(J_original))
  }
  if (anyDuplicated(colnames(Y_full))) {
    warning("Duplicate taxon names were made unique.", call. = FALSE)
    colnames(Y_full) <- make.unique(colnames(Y_full))
  }

  if (!is.numeric(d0) || length(d0) != 1L || !is.finite(d0) || d0 <= 0) {
    stop("`d0` must be one positive finite number.", call. = FALSE)
  }
  m_abund <- validate_integer_scalar(m_abund, "m_abund", minimum = 1L)
  min_positive <- validate_integer_scalar(
    min_positive,
    "min_positive",
    minimum = 0L
  )
  affine_folds <- validate_integer_scalar(
    affine_folds,
    "affine_folds",
    minimum = 2L
  )
  lts_seed <- validate_integer_scalar(lts_seed, "lts_seed", minimum = 0L)
  return_posteriors <- validate_flag(return_posteriors, "return_posteriors")
  if (!is.list(fit_control)) {
    stop("`fit_control` must be a list.", call. = FALSE)
  }
  abundance_eligibility <- match.arg(abundance_eligibility)

  group_info <- encode_group(group, n_original, case_level = case_level)
  g_full <- group_info$g
  Z_full <- prepare_covariates(covariates, n_original)

  D_full <- if (is.null(library_size)) {
    rowSums(Y_full)
  } else {
    as.numeric(library_size)
  }

  if (length(D_full) != n_original || any(!is.finite(D_full)) || any(D_full <= 0)) {
    stop("`library_size` must contain one positive finite value per sample.")
  }
  if (any(D_full + 1e-8 < rowSums(Y_full))) {
    stop("Every `library_size` must be at least the supplied count-table row sum.")
  }
  if (any(D_full <= d0)) {
    stop("Every sequencing depth must exceed `d0` for the continuous working support.")
  }

  positive_counts <- Y_full[Y_full > 0]
  if (length(positive_counts) > 0L && d0 > min(positive_counts) + 1e-8) {
    stop("`d0` cannot exceed an observed positive count under the working support.")
  }

  keep_tax <- if (min_positive > 0L) {
    colSums(Y_full > 0) >= min_positive
  } else {
    rep(TRUE, J_original)
  }

  Y <- Y_full[, keep_tax, drop = FALSE]
  if (ncol(Y) < 2L) {
    stop("Fewer than two taxa remain after taxon retention.")
  }

  keep_sample <- rowSums(Y) > 0
  Y <- Y[keep_sample, , drop = FALSE]
  D <- D_full[keep_sample]
  g <- g_full[keep_sample]
  Z <- if (is.null(Z_full)) NULL else Z_full[keep_sample, , drop = FALSE]

  if (length(unique(g)) != 2L) {
    stop("Only one group remains after dropping retained-table-empty samples.")
  }

  ## Global fixed-design checks before taxon-wise fitting.
  X_rho_null <- matrix(1, nrow = nrow(Y), ncol = 1L)
  colnames(X_rho_null) <- "Intercept"
  if (!is.null(Z)) {
    X_rho_null <- cbind(X_rho_null, Z)
  }
  X_full <- design_matrix(g = g, z = Z, include_group = TRUE)

  if (!is_full_rank(X_rho_null)) {
    stop("The intercept/covariate structural-null design is rank deficient.")
  }
  if (!is_full_rank(X_full)) {
    stop("Group is collinear with the intercept/covariates after preprocessing.")
  }

  J <- ncol(Y)
  taxa <- colnames(Y)

  absence_list <- lapply(seq_len(J), function(j) {
    structural_absence_test(
      y = Y[, j],
      D = D,
      g = g,
      z = Z,
      d0 = d0,
      fit_control = fit_control
    )
  })

  auxiliary_list <- lapply(seq_len(J), function(j) {
    auxiliary_zero_fit(
      y = Y[, j],
      D = D,
      g = g,
      z = Z,
      d0 = d0,
      fit_control = fit_control
    )
  })

  p_absence <- vapply(absence_list, function(x) x$p, numeric(1))
  z_absence <- vapply(absence_list, function(x) x$z, numeric(1))
  absence_formed <- vapply(absence_list, function(x) x$formed, logical(1))

  Gamma_aux <- do.call(cbind, lapply(auxiliary_list, function(x) x$gamma))
  Ccor <- do.call(cbind, lapply(auxiliary_list, function(x) x$censor_correction))
  dimnames(Gamma_aux) <- dimnames(Y)
  dimnames(Ccor) <- dimnames(Y)
  W <- 1 - Gamma_aux

  abundance <- hric_wls_p(
    Y = Y,
    D = D,
    g = g,
    W = W,
    z = Z,
    m_abund = as.integer(m_abund),
    abundance_eligibility = abundance_eligibility,
    Ccor = Ccor,
    affine_folds = as.integer(affine_folds),
    lts_seed = as.integer(lts_seed)
  )

  p_abundance <- abundance$p
  abundance_formed <- abundance$tested

  ## Fixed two-component family: an unavailable component is assigned one, and
  ## the Bonferroni factor remains two. This avoids data-dependent removal of the
  ## multiplicity penalty.
  p_bonf <- vapply(seq_len(J), function(j) {
    bonferroni_two_component(p_absence[j], p_abundance[j])
  }, numeric(1))

  p_cauchy <- vapply(seq_len(J), function(j) {
    if (isTRUE(absence_formed[j]) && isTRUE(abundance_formed[j])) {
      cauchy_combination(c(p_absence[j], p_abundance[j]))
    } else {
      NA_real_
    }
  }, numeric(1))

  null_fit_valid <- vapply(
    absence_list,
    function(x) isTRUE(x$fit$valid),
    logical(1)
  )
  null_fit_reason <- vapply(
    absence_list,
    function(x) x$fit$reason %||% NA_character_,
    character(1)
  )
  null_fit_method <- vapply(
    absence_list,
    function(x) x$fit$optimization$method %||% NA_character_,
    character(1)
  )
  null_fit_boundary <- vapply(
    absence_list,
    function(x) x$fit$optimization$boundary %||% NA,
    logical(1)
  )
  null_fit_gradient <- vapply(
    absence_list,
    function(x) x$fit$optimization$gradient_mean_norm %||% NA_real_,
    numeric(1)
  )
  null_fit_condition <- vapply(
    absence_list,
    function(x) x$fit$optimization$information_condition %||% Inf,
    numeric(1)
  )

  aux_fit_valid <- vapply(auxiliary_list, function(x) x$valid, logical(1))
  aux_fit_reason <- vapply(
    auxiliary_list,
    function(x) x$fit$reason %||% NA_character_,
    character(1)
  )
  aux_fit_method <- vapply(
    auxiliary_list,
    function(x) x$fit$optimization$method %||% NA_character_,
    character(1)
  )
  aux_fit_boundary <- vapply(
    auxiliary_list,
    function(x) x$fit$optimization$boundary %||% NA,
    logical(1)
  )
  aux_fit_gradient <- vapply(
    auxiliary_list,
    function(x) x$fit$optimization$gradient_mean_norm %||% NA_real_,
    numeric(1)
  )
  aux_fit_condition <- vapply(
    auxiliary_list,
    function(x) x$fit$optimization$information_condition %||% Inf,
    numeric(1)
  )

  structural_group_log_odds <- vapply(
    auxiliary_list,
    function(x) x$structural_group_log_odds,
    numeric(1)
  )
  structural_group_se <- vapply(
    auxiliary_list,
    function(x) x$structural_group_se,
    numeric(1)
  )
  latent_group_coefficient <- vapply(
    auxiliary_list,
    function(x) x$latent_group_coefficient,
    numeric(1)
  )
  latent_group_se <- vapply(
    auxiliary_list,
    function(x) x$latent_group_se,
    numeric(1)
  )

  out <- data.frame(
    taxon = taxa,
    p_absence = unname(p_absence),
    z_absence = unname(z_absence),
    absence_formed = unname(absence_formed),
    structural_group_log_odds_aux = unname(structural_group_log_odds),
    structural_group_se_aux = unname(structural_group_se),
    latent_group_coefficient_aux = unname(latent_group_coefficient),
    latent_group_se_aux = unname(latent_group_se),
    p_abundance = unname(p_abundance),
    beta_abundance_raw = unname(abundance$beta_raw),
    se_abundance_raw_hc3 = unname(abundance$se_raw_hc3),
    beta_abundance_corrected = unname(abundance$beta_corrected),
    se_abundance_joint_hc3 = unname(abundance$se_joint),
    z_abundance = unname(abundance$z),
    abundance_formed = unname(abundance_formed),
    p_bonf_minp = unname(p_bonf),
    p_cauchy = unname(p_cauchy),
    n_positive_group0 = unname(abundance$n_positive_group0),
    n_positive_group1 = unname(abundance$n_positive_group1),
    censor_column_used = unname(abundance$censor_column_used),
    affine_partition = unname(abundance$affine_partition),
    affine_scheme = unname(abundance$affine_scheme),
    affine_training_taxa = unname(abundance$affine_training_n),
    affine_inlier_taxa = unname(abundance$affine_inlier_n),
    affine_reference_root_abundance = unname(abundance$v_reference),
    affine_background_intercept = unname(abundance$background_intercept),
    affine_background_slope = unname(abundance$background_slope),
    null_fit_valid = unname(null_fit_valid),
    null_fit_reason = unname(null_fit_reason),
    null_fit_method = unname(null_fit_method),
    null_fit_boundary = unname(null_fit_boundary),
    null_fit_gradient_mean_norm = unname(null_fit_gradient),
    null_fit_information_condition = unname(null_fit_condition),
    auxiliary_fit_valid = unname(aux_fit_valid),
    auxiliary_fit_reason = unname(aux_fit_reason),
    auxiliary_fit_method = unname(aux_fit_method),
    auxiliary_fit_boundary = unname(aux_fit_boundary),
    auxiliary_fit_gradient_mean_norm = unname(aux_fit_gradient),
    auxiliary_fit_information_condition = unname(aux_fit_condition),
    stringsAsFactors = FALSE
  )


  attr(out, "group_levels") <- group_info$levels
  attr(out, "kept_taxa") <- which(keep_tax)
  attr(out, "kept_samples") <- which(keep_sample)
  attr(out, "library_size") <- D
  attr(out, "retained_table_total") <- rowSums(Y)
  attr(out, "covariate_matrix") <- Z
  attr(out, "working_likelihood") <-
    "upper-truncated-normal mixed discrete-continuous working likelihood"
  attr(out, "abundance_variance_scope") <-
    "conditional HC3; use dash_bootstrap() for full-pipeline sensitivity"
  attr(out, "abundance_eligibility") <- abundance_eligibility

  if (isTRUE(return_posteriors)) {
    Gamma_abs <- do.call(cbind, lapply(absence_list, function(x) x$gamma))
    Rho_abs <- do.call(cbind, lapply(absence_list, function(x) x$rho))
    dimnames(Gamma_abs) <- dimnames(Y)
    dimnames(Rho_abs) <- dimnames(Y)

    attr(out, "posteriors") <- list(
      gamma_absence_null = Gamma_abs,
      rho_absence_null = Rho_abs,
      gamma_auxiliary = Gamma_aux,
      censor_severity = Ccor,
      counts_retained = Y,
      group_binary = g,
      library_size = D,
      retained_total = rowSums(Y)
    )
  }

  out
}


# Full-pipeline bootstrap sensitivity analysis

#' Full-pipeline bootstrap sensitivity analysis for DASH
#'
#' Resamples samples with replacement within group, conditional on the samples
#' and taxa retained by the primary analysis. Every bootstrap replicate refits
#' both zero models, posterior weights, censoring covariates, HRIC regressions,
#' cross-fitted affine backgrounds, and corrected abundance coefficients.
#'
#' @details
#' Percentile intervals are calculated from finite bootstrap corrected
#' coefficients. The bootstrap p-value is a two-sided Wald sensitivity
#' calculation using the bootstrap standard error and the primary corrected
#' coefficient. Because robust affine selection is discontinuous, this function
#' is intended as a sensitivity analysis rather than universally exact
#' post-selection inference.
#'
#' @param B An integer of at least 20 giving the number of bootstrap replicates.
#' @param seed A non-negative integer bootstrap seed.
#' @param min_success_fraction The minimum fraction of finite bootstrap
#'   corrected coefficients required for a taxon-specific bootstrap result.
#' @inheritParams dash
#'
#' @return The data frame returned by [dash()] with additional columns for the
#'   number of successful taxon-specific bootstrap replicates, bootstrap
#'   standard errors, percentile intervals, bootstrap abundance p-values, and
#'   bootstrap Bonferroni min-P omnibus p-values. The complete matrix of
#'   bootstrap corrected coefficients is stored in the
#'   `bootstrap_beta_matrix` attribute.
#'
#' @examples
#' \dontrun{
#' boot <- dash_bootstrap(
#'   counts = counts,
#'   group = group,
#'   library_size = depth,
#'   B = 500,
#'   seed = 20260710
#' )
#' }
#'
#' @family DASH functions
#' @export
#' @md
dash_bootstrap <- function(
    counts,
    group,
    covariates = NULL,
    library_size = NULL,
    case_level = NULL,
    d0 = 1,
    m_abund = 5L,
    abundance_eligibility = c("pooled", "per_group", "none"),
    min_positive = 3L,
    affine_folds = 5L,
    lts_seed = 1907L,
    fit_control = list(),
    B = 200L,
    seed = 20260710L,
    min_success_fraction = 0.8) {

  B <- validate_integer_scalar(B, "B", minimum = 20L)
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  if (!is.finite(min_success_fraction) ||
      min_success_fraction <= 0 || min_success_fraction > 1) {
    stop("`min_success_fraction` must lie in (0, 1].")
  }
  abundance_eligibility <- match.arg(abundance_eligibility)

  base <- dash(
    counts = counts,
    group = group,
    covariates = covariates,
    library_size = library_size,
    case_level = case_level,
    d0 = d0,
    m_abund = m_abund,
    abundance_eligibility = abundance_eligibility,
    min_positive = min_positive,
    affine_folds = affine_folds,
    lts_seed = lts_seed,
    fit_control = fit_control,
    return_posteriors = FALSE
  )

  Y_full <- as.matrix(counts)
  D_full <- if (is.null(library_size)) rowSums(Y_full) else as.numeric(library_size)
  group_info <- encode_group(group, nrow(Y_full), case_level = case_level)
  Z_full <- prepare_covariates(covariates, nrow(Y_full))

  keep_sample <- attr(base, "kept_samples")
  keep_tax <- attr(base, "kept_taxa")

  Y <- Y_full[keep_sample, keep_tax, drop = FALSE]
  D <- D_full[keep_sample]
  g <- group_info$g[keep_sample]
  Z <- if (is.null(Z_full)) NULL else Z_full[keep_sample, , drop = FALSE]

  idx0 <- which(g == 0)
  idx1 <- which(g == 1)
  taxa <- base$taxon
  boot_beta <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(taxa),
    dimnames = list(NULL, taxa)
  )
  replicate_ok <- rep(FALSE, B)

  with_preserved_seed(seed, {
    for (b in seq_len(B)) {
      boot_index <- c(
        sample(idx0, length(idx0), replace = TRUE),
        sample(idx1, length(idx1), replace = TRUE)
      )

      fit_b <- tryCatch(
        dash(
          counts = Y[boot_index, , drop = FALSE],
          group = g[boot_index],
          covariates = if (is.null(Z)) NULL else Z[boot_index, , drop = FALSE],
          library_size = D[boot_index],
          d0 = d0,
          m_abund = m_abund,
          abundance_eligibility = abundance_eligibility,
          min_positive = 0L,
          affine_folds = affine_folds,
          lts_seed = as.integer(lts_seed + b),
          fit_control = fit_control,
          return_posteriors = FALSE
        ),
        error = function(e) NULL
      )

      if (!is.null(fit_b)) {
        match_idx <- match(taxa, fit_b$taxon)
        boot_beta[b, ] <- fit_b$beta_abundance_corrected[match_idx]
        replicate_ok[b] <- TRUE
      }
    }
  })

  n_success <- colSums(is.finite(boot_beta))
  min_success <- max(20L, ceiling(B * min_success_fraction))

  se_boot <- rep(NA_real_, length(taxa))
  ci_low <- rep(NA_real_, length(taxa))
  ci_high <- rep(NA_real_, length(taxa))
  p_boot <- rep(NA_real_, length(taxa))

  for (j in seq_along(taxa)) {
    vals <- boot_beta[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) < min_success) {
      next
    }

    se_j <- stats::sd(vals)
    if (!is.finite(se_j) || se_j <= 0) {
      next
    }

    se_boot[j] <- se_j
    ci <- stats::quantile(vals, probs = c(0.025, 0.975), names = FALSE)
    ci_low[j] <- ci[1L]
    ci_high[j] <- ci[2L]
    p_boot[j] <- 2 * stats::pnorm(
      -abs(base$beta_abundance_corrected[j] / se_j)
    )
  }

  base$bootstrap_successes_abundance <- n_success
  base$se_abundance_full_pipeline_bootstrap <- se_boot
  base$ci95_abundance_bootstrap_low <- ci_low
  base$ci95_abundance_bootstrap_high <- ci_high
  base$p_abundance_full_pipeline_bootstrap <- p_boot

  base$p_bonf_minp_full_pipeline_bootstrap <- vapply(
    seq_len(nrow(base)),
    function(j) {
      bonferroni_two_component(
        base$p_absence[j],
        base$p_abundance_full_pipeline_bootstrap[j]
      )
    },
    numeric(1)
  )

  attr(base, "bootstrap_replicates_requested") <- B
  attr(base, "bootstrap_replicates_completed") <- sum(replicate_ok)
  attr(base, "bootstrap_seed") <- seed
  attr(base, "bootstrap_beta_matrix") <- boot_beta
  base
}


# Group-label-swap diagnostic

#' Check numerical invariance to reversing group coding
#'
#' Runs [dash()] with the original and reversed binary group coding. Corrected
#' abundance coefficients should reverse sign, whereas component p-values,
#' omnibus p-values, and abundance standard errors should remain unchanged up
#' to numerical tolerance.
#'
#' @inheritParams dash
#'
#' @return A data frame containing taxon-specific absolute discrepancies. The
#'   `maximum_discrepancies` attribute contains the largest finite discrepancy
#'   for each checked quantity.
#'
#' @examples
#' \dontrun{
#' check <- dash_label_swap_check(
#'   counts = counts,
#'   group = group,
#'   library_size = depth
#' )
#' attr(check, "maximum_discrepancies")
#' }
#'
#' @family DASH functions
#' @export
#' @md
dash_label_swap_check <- function(
    counts,
    group,
    covariates = NULL,
    library_size = NULL,
    case_level = NULL,
    d0 = 1,
    m_abund = 5L,
    abundance_eligibility = c("pooled", "per_group", "none"),
    min_positive = 3L,
    affine_folds = 5L,
    lts_seed = 1907L,
    fit_control = list()) {

  n <- nrow(as.matrix(counts))
  abundance_eligibility <- match.arg(abundance_eligibility)
  group_info <- encode_group(group, n, case_level = case_level)
  g <- group_info$g

  fit_forward <- dash(
    counts = counts,
    group = g,
    covariates = covariates,
    library_size = library_size,
    d0 = d0,
    m_abund = m_abund,
    abundance_eligibility = abundance_eligibility,
    min_positive = min_positive,
    affine_folds = affine_folds,
    lts_seed = lts_seed,
    fit_control = fit_control
  )

  fit_swapped <- dash(
    counts = counts,
    group = 1L - g,
    covariates = covariates,
    library_size = library_size,
    d0 = d0,
    m_abund = m_abund,
    abundance_eligibility = abundance_eligibility,
    min_positive = min_positive,
    affine_folds = affine_folds,
    lts_seed = lts_seed,
    fit_control = fit_control
  )

  j2 <- match(fit_forward$taxon, fit_swapped$taxon)

  out <- data.frame(
    taxon = fit_forward$taxon,
    absence_p_absolute_difference = abs(
      fit_forward$p_absence - fit_swapped$p_absence[j2]
    ),
    abundance_p_absolute_difference = abs(
      fit_forward$p_abundance - fit_swapped$p_abundance[j2]
    ),
    omnibus_p_absolute_difference = abs(
      fit_forward$p_bonf_minp - fit_swapped$p_bonf_minp[j2]
    ),
    corrected_beta_sign_reversal_error = abs(
      fit_forward$beta_abundance_corrected +
        fit_swapped$beta_abundance_corrected[j2]
    ),
    abundance_se_absolute_difference = abs(
      fit_forward$se_abundance_joint_hc3 -
        fit_swapped$se_abundance_joint_hc3[j2]
    ),
    stringsAsFactors = FALSE
  )

  finite_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else max(x)
  }

  attr(out, "maximum_discrepancies") <- c(
    p_absence = finite_max(out$absence_p_absolute_difference),
    p_abundance = finite_max(out$abundance_p_absolute_difference),
    p_omnibus = finite_max(out$omnibus_p_absolute_difference),
    beta_sign_reversal = finite_max(out$corrected_beta_sign_reversal_error),
    abundance_se = finite_max(out$abundance_se_absolute_difference)
  )

  out
}
