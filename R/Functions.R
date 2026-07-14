## =============================================================================
##  DASH: Depth-Aware Structural-zero Hellinger test
##  Reference implementation -- single-file companion to the manuscript
## =============================================================================
##
##  Author : Yiqian Zhang
##  License: GPL (>= 3)
##
##  DASH is a two-part, per-taxon differential abundance test for microbiome
##  count data. For each taxon it forms two complementary tests and combines
##  them into a single omnibus test:
##
##    * Prevalence arm. A structural-zero score test derived from a per-taxon
##      structural-zero mixture model with a left-censored normal observation
##      component. Library size enters through a depth-aware censoring threshold,
##      so an observed zero in a deep sample gives stronger evidence of true
##      absence than a zero in a shallow sample.
##    * Abundance arm. A weighted least squares test on Hellinger-Riemann
##      intrinsic coordinates (HRIC), with auxiliary posterior-presence weights,
##      a censoring-correction nuisance covariate, a robust affine (LTS)
##      correction for the compositional background, and a joint sample-level
##      HC3 variance for the affine-corrected coefficient.
##
##  The two component p-values are combined by a Bonferroni min-P rule (primary)
##  and a Cauchy combination (sensitivity). Multiple testing across taxa is left
##  to the caller; apply Benjamini-Hochberg, stats::p.adjust(., "BH"), to the
##  reported column.
##
##  The zero model can be estimated either by direct observed-data likelihood
##  optimization or by EM. Direct optimization reproduces the original DASH
##  implementation and is the default. The abundance component uses the revised
##  joint-HC3 variance. The default values are d0 = 1, m_abund = 5, and
##  min_positive = 3.
##
##  Entry points
##    dash()           the test: input coercion, feature retention, both arms,
##                     and the per-taxon combination, returned as a tidy table
##    hric_transform() the HRIC transform, usable on its own
##
##  Dependencies: base R, 'stats', and 'MASS'. The affine correction uses
##  MASS::lqs for least trimmed squares inlier selection.
##
##  Contents
##    1. Internal utilities
##    2. P-value combination rules
##    3. HRIC transform
##    4. Prevalence arm: structural-zero mixture and the structural-zero score test
##    5. Abundance arm: HRIC weighted least squares with auxiliary weights,
##       censoring correction, and affine bias correction
##    6. DASH: the test
## =============================================================================

#' @importFrom stats optim qlogis plogis pnorm dnorm pchisq sd mad median qnorm coef p.adjust
NULL


## -----------------------------------------------------------------------------
## 1. Internal utilities
## -----------------------------------------------------------------------------

#' Clamp a vector to a range
#' @keywords internal
#' @noRd
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' Standardize a vector; return zeros if the s.d. is not positive
#' @keywords internal
#' @noRd
safe_scale <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }
  
  as.numeric((x - m) / s)
}

#' Solve a matrix, falling back to a QR solve
#' @keywords internal
#' @noRd
safe_solve <- function(A) tryCatch(solve(A), error = function(e) qr.solve(A, diag(ncol(A))))

#' Stable lower-tail inverse Mills ratio phi(t) / Phi(t)
#' @keywords internal
#' @noRd
lower_mills_ratio <- function(t) {
  t <- as.numeric(t)
  log_ratio <- dnorm(t, log = TRUE) - pnorm(t, log.p = TRUE)
  ratio <- exp(log_ratio)
  
  bad <- !is.finite(ratio)
  use_asymptotic <- bad & is.finite(t) & t < 0
  ratio[use_asymptotic] <- -t[use_asymptotic]
  ratio[!is.finite(ratio) | ratio < 0] <- 0
  ratio
}

#' Coerce covariates to a named matrix (or NULL)
#' @keywords internal
#' @noRd
make_cov_matrix <- function(z = NULL) { if (is.null(z)) return(NULL); Z <- as.matrix(z); if (ncol(Z) == 0) return(NULL); colnames(Z) <- paste0("z", seq_len(ncol(Z))); Z }

#' Design matrix: intercept, optional group, optional covariates
#' @keywords internal
#' @noRd
design_matrix <- function(g = NULL, z = NULL, include_group = TRUE) {
  nn <- if (!is.null(g)) length(g) else nrow(as.matrix(z))
  X <- matrix(1, nn, 1); colnames(X) <- "Intercept"
  if (include_group && !is.null(g)) X <- cbind(X, Group = as.numeric(g))
  Z <- make_cov_matrix(z)
  if (!is.null(Z)) X <- cbind(X, Z)
  X
}

#' Residualize the group indicator on the measured covariates
#' @keywords internal
#' @noRd
residualize_group <- function(g, z = NULL) {
  if (is.null(z)) return(as.numeric(g) - mean(g))
  Z <- make_cov_matrix(z); X <- cbind(Intercept = 1, Z)
  bh <- safe_solve(crossprod(X)) %*% crossprod(X, as.numeric(g))
  as.numeric(g) - as.numeric(X %*% bh)
}

#' Robust empirical-null scale (defined for completeness; not used in the
#' DASH testing path, matching the simulation).
#' @keywords internal
#' @noRd
empirical_null_scale <- function(z, min_taxa = 30) {
  zz <- z[is.finite(z)]
  
  if (length(zz) < min_taxa) return(1)
  
  s <- median(abs(zz - median(zz, na.rm = TRUE)), na.rm = TRUE) / qnorm(0.75)
  
  if (!is.finite(s) || s <= 0) s <- 1
  
  ## Only correct anti-conservative overdispersion.
  ## Do not make the method more liberal when s < 1.
  max(1, min(s, 3))
}

## -----------------------------------------------------------------------------
## 2. P-value combination rules
## -----------------------------------------------------------------------------

#' Cauchy combination of p-values
#' @keywords internal
#' @noRd
cauchy_combination <- function(ps) { ps <- ps[is.finite(ps)]; if (!length(ps)) return(NA_real_); ps <- pmin(pmax(ps, 1e-15), 1 - 1e-15); 0.5 - atan(mean(tan((0.5 - ps) * pi))) / pi }

#' Bonferroni min-P combination of p-values
#' @keywords internal
#' @noRd
bonf_minp <- function(ps) { ps <- ps[is.finite(ps)]; if (!length(ps)) return(NA_real_); min(1, length(ps) * min(ps)) }

## -----------------------------------------------------------------------------
## 3. HRIC transform
## -----------------------------------------------------------------------------

#' Hellinger-Riemann intrinsic coordinates (HRIC) transform
#'
#' Map a non-negative abundance matrix to Hellinger-Riemann intrinsic
#' coordinates (HRIC). Each row is closed to a composition, mapped to the square-root
#' simplex, and contrasted against the uniform composition with an arcsine
#' (geodesic) scaling.
#'
#' @param X A numeric matrix of non-negative abundances (samples x taxa), with
#'   at least two columns and a positive total in every row.
#' @return A matrix of HRIC coordinates with the same dimensions and names as
#'   `X`.
#' @export
hric_transform <- function(X) {
  X <- as.matrix(X); if (!is.numeric(X)) stop("X must be numeric."); if (any(X < 0, na.rm = TRUE)) stop("X must contain non-negative values only."); if (anyNA(X)) stop("X contains NA values.")
  k <- ncol(X); if (k < 2) stop("X must contain at least two components.")
  rs <- rowSums(X); if (any(rs <= 0)) stop("Each row must have positive total abundance.")
  Pi <- sweep(X, 1, rs, FUN = "/"); Pi <- pmax(Pi, 0)
  sqrt_Pi <- sqrt(Pi); sqrt_pi0 <- rep(1 / sqrt(k), k)
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0); c_pi <- pmin(pmax(c_pi, 0), 1)
  s_pi <- sqrt(pmax(1 - c_pi^2, 0)); scale_factor <- asin(s_pi) / s_pi; scale_factor[s_pi < 1e-12] <- 1
  centered <- sqrt_Pi - outer(c_pi, sqrt_pi0); T <- sweep(centered, 1, scale_factor, "*")
  rownames(T) <- rownames(X); colnames(T) <- colnames(X); T
}

## -----------------------------------------------------------------------------
## 4. Prevalence arm: structural-zero mixture
## -----------------------------------------------------------------------------

#' Fit the per-taxon structural-zero mixture model by direct optimization
#' of the observed-data likelihood
#' @keywords internal
#' @noRd
fit_zero_optim <- function(y, N, z = NULL, d0 = 1.0,
                           g = NULL,
                           include_group_rho = FALSE,
                           include_group_eta = FALSE) {
  y <- as.numeric(y)
  N <- as.numeric(N)
  pos <- y > 0
  
  ell <- ifelse(pos, log(pmax(y, 1) / N), 0)
  a_i <- log(d0 / N)
  
  make_zero_design <- function(include_group) {
    if (include_group && is.null(g)) {
      stop("g must be supplied when include_group=TRUE.")
    }
    
    if (include_group) {
      design_matrix(g = g, z = z, include_group = TRUE)
    } else if (!is.null(z)) {
      design_matrix(g = NULL, z = z, include_group = FALSE)
    } else {
      X <- matrix(1, length(y), 1)
      colnames(X) <- "Intercept"
      X
    }
  }
  
  X_rho <- make_zero_design(include_group_rho)
  X_eta <- make_zero_design(include_group_eta)
  
  npar_rho <- ncol(X_rho)
  npar_eta <- ncol(X_eta)
  
  p0 <- min(max(mean(!pos) * 0.5, 0.02), 0.9)
  init_eta <- if (any(pos)) {
    mean(log(pmax(y[pos] / N[pos], 1e-8)))
  } else {
    -5
  }
  
  init <- c(
    c(qlogis(p0), rep(0, npar_rho - 1)),
    c(init_eta, rep(0, npar_eta - 1)),
    log(1.0)
  )
  
  lower <- c(rep(-10, npar_rho), rep(-30, npar_eta), log(0.1))
  upper <- c(rep(10, npar_rho), rep(5, npar_eta), log(8.0))
  
  ## This is the original direct observed-data likelihood objective.
  nll <- function(par) {
    alpha <- par[seq_len(npar_rho)]
    eta_coef <- par[npar_rho + seq_len(npar_eta)]
    sig <- exp(par[npar_rho + npar_eta + 1L])
    
    if (!is.finite(sig) || sig <= 0) return(1e10)
    
    rho <- clamp(
      plogis(as.numeric(X_rho %*% alpha)),
      1e-12,
      1 - 1e-12
    )
    eta <- as.numeric(X_eta %*% eta_coef)
    
    tval <- (a_i - eta) / sig
    L_zero <- rho + (1 - rho) * pnorm(tval)
    L_positive <- (1 - rho) * dnorm((ell - eta) / sig) / sig
    L <- ifelse(pos, L_positive, L_zero)
    
    -sum(log(pmax(L, 1e-300)))
  }
  
  ## Preserve the original bounded L-BFGS-B fit and Nelder-Mead fallback.
  primary_fit <- tryCatch(
    optim(
      init,
      nll,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 500L, factr = 1e7)
    ),
    error = function(e) NULL
  )
  
  fallback_fit <- NULL
  if (!is.null(primary_fit) && is.finite(primary_fit$value)) {
    fit_used <- primary_fit
    par <- primary_fit$par
    optimizer_used <- "L-BFGS-B"
  } else {
    fallback_fit <- tryCatch(
      optim(
        init,
        nll,
        method = "Nelder-Mead",
        control = list(maxit = 500L, reltol = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fallback_fit)) {
      fit_used <- fallback_fit
      par <- fallback_fit$par
      optimizer_used <- "Nelder-Mead"
    } else {
      fit_used <- NULL
      par <- init
      optimizer_used <- "initialization"
    }
  }
  
  alpha <- par[seq_len(npar_rho)]
  eta_coef <- par[npar_rho + seq_len(npar_eta)]
  sig <- max(exp(par[npar_rho + npar_eta + 1L]), 0.1)
  
  rho <- clamp(
    plogis(as.numeric(X_rho %*% alpha)),
    1e-12,
    1 - 1e-12
  )
  eta <- as.numeric(X_eta %*% eta_coef)
  tval <- (a_i - eta) / sig
  Phi <- pnorm(tval)
  denom <- pmax(rho + (1 - rho) * Phi, 1e-300)
  gamma <- ifelse(pos, 0, rho / denom)
  
  final_nll <- tryCatch(nll(par), error = function(e) NA_real_)
  valid <-
    is.finite(final_nll) &&
    is.finite(sig) &&
    sig > 0 &&
    all(is.finite(rho)) &&
    all(is.finite(eta)) &&
    all(is.finite(gamma))
  
  optimization_succeeded <- !is.null(fit_used)
  converged <-
    optimization_succeeded &&
    isTRUE(fit_used$convergence == 0)
  
  iterations <- if (
    optimization_succeeded &&
    !is.null(fit_used$counts) &&
    !is.null(fit_used$counts[["function"]])
  ) {
    as.integer(fit_used$counts[["function"]])
  } else {
    NA_integer_
  }
  
  list(
    alpha = alpha,
    eta_coef = eta_coef,
    sig = sig,
    X_rho = X_rho,
    X_eta = X_eta,
    rho = rho,
    eta = eta,
    gamma = gamma,
    logLik = if (is.finite(final_nll)) -final_nll else NA_real_,
    valid = valid,
    converged = converged,
    optimization_succeeded = optimization_succeeded,
    iterations = iterations,
    optimizer = optimizer_used,
    method = "optim"
  )
}

#' Fit the per-taxon structural-zero mixture model by EM
#' @keywords internal
#' @noRd
fit_zero_em <- function(y, N, z = NULL, d0 = 1.0,
                        g = NULL,
                        include_group_rho = FALSE,
                        include_group_eta = FALSE,
                        maxit = 1000L,
                        tol = 1e-8) {
  y <- as.numeric(y)
  N <- as.numeric(N)
  
  pos <- y > 0
  zero <- !pos
  
  ell <- ifelse(pos, log(pmax(y, 1) / N), 0)
  a_i <- log(d0 / N)
  
  make_zero_design <- function(include_group) {
    if (include_group && is.null(g)) {
      stop("g must be supplied when include_group=TRUE.")
    }
    
    if (include_group) {
      design_matrix(g = g, z = z, include_group = TRUE)
    } else if (!is.null(z)) {
      design_matrix(g = NULL, z = z, include_group = FALSE)
    } else {
      X <- matrix(1, length(y), 1)
      colnames(X) <- "Intercept"
      X
    }
  }
  
  X_rho <- make_zero_design(include_group_rho)
  X_eta <- make_zero_design(include_group_eta)
  
  npar_rho <- ncol(X_rho)
  npar_eta <- ncol(X_eta)
  
  p0 <- min(max(mean(zero) * 0.5, 0.02), 0.9)
  
  init_eta <- if (any(pos)) {
    mean(log(pmax(y[pos] / N[pos], 1e-8)))
  } else {
    -5
  }
  
  alpha <- c(qlogis(p0), rep(0, npar_rho - 1))
  eta_coef <- c(init_eta, rep(0, npar_eta - 1))
  omega <- log(1.0)
  
  lower_rho <- rep(-10, npar_rho)
  upper_rho <- rep(10, npar_rho)
  
  lower_eta <- c(rep(-30, npar_eta), log(0.1))
  upper_eta <- c(rep(5, npar_eta), log(8.0))
  
  log_expit <- function(x) {
    -log1p(exp(-abs(x))) + pmin(x, 0)
  }
  
  log1m_expit <- function(x) {
    -log1p(exp(-abs(x))) - pmax(x, 0)
  }
  
  e_step <- function(alpha, eta_coef, omega) {
    sig <- exp(omega)
    
    rho <- clamp(
      plogis(as.numeric(X_rho %*% alpha)),
      1e-12,
      1 - 1e-12
    )
    
    eta <- as.numeric(X_eta %*% eta_coef)
    tval <- (a_i - eta) / sig
    Phi <- pnorm(tval)
    
    gamma <- numeric(length(y))
    
    if (any(zero)) {
      denom <- pmax(
        rho[zero] + (1 - rho[zero]) * Phi[zero],
        1e-300
      )
      
      gamma[zero] <- rho[zero] / denom
    }
    
    gamma[pos] <- 0
    
    list(
      gamma = gamma,
      rho = rho,
      eta = eta,
      sig = sig,
      tval = tval
    )
  }
  
  observed_loglik <- function(alpha, eta_coef, omega) {
    sig <- exp(omega)
    
    lp_rho <- as.numeric(X_rho %*% alpha)
    rho <- clamp(plogis(lp_rho), 1e-12, 1 - 1e-12)
    eta <- as.numeric(X_eta %*% eta_coef)
    
    log_contribution <- numeric(length(y))
    
    if (any(zero)) {
      Phi0 <- pnorm(
        (a_i[zero] - eta[zero]) / sig
      )
      
      log_contribution[zero] <- log(
        pmax(
          rho[zero] +
            (1 - rho[zero]) * Phi0,
          1e-300
        )
      )
    }
    
    if (any(pos)) {
      rpos <- (ell[pos] - eta[pos]) / sig
      
      log_contribution[pos] <-
        log1m_expit(lp_rho[pos]) +
        dnorm(rpos, log = TRUE) -
        omega
    }
    
    sum(log_contribution)
  }
  
  ## Negative structural-absence Q function.
  nQ_rho <- function(alpha_candidate, gamma) {
    lp <- as.numeric(X_rho %*% alpha_candidate)
    
    value <- -sum(
      gamma * log_expit(lp) +
        (1 - gamma) * log1m_expit(lp)
    )
    
    if (!is.finite(value)) 1e100 else value
  }
  
  gr_nQ_rho <- function(alpha_candidate, gamma) {
    lp <- as.numeric(X_rho %*% alpha_candidate)
    rho <- plogis(lp)
    
    -as.numeric(
      crossprod(X_rho, gamma - rho)
    )
  }
  
  ## Negative present-conditional Q function.
  nQ_eta <- function(par_candidate, gamma) {
    eta_candidate <- as.numeric(
      X_eta %*% par_candidate[seq_len(npar_eta)]
    )
    
    omega_candidate <- par_candidate[npar_eta + 1L]
    sig_candidate <- exp(omega_candidate)
    
    value <- 0
    
    if (any(zero)) {
      t0 <- (
        a_i[zero] - eta_candidate[zero]
      ) / sig_candidate
      
      value <- value + sum(
        (1 - gamma[zero]) *
          pnorm(t0, log.p = TRUE)
      )
    }
    
    if (any(pos)) {
      rpos <- (
        ell[pos] - eta_candidate[pos]
      ) / sig_candidate
      
      value <- value + sum(
        dnorm(rpos, log = TRUE) -
          omega_candidate
      )
    }
    
    value <- -value
    
    if (!is.finite(value)) 1e100 else value
  }
  
  gr_nQ_eta <- function(par_candidate, gamma) {
    eta_candidate <- as.numeric(
      X_eta %*% par_candidate[seq_len(npar_eta)]
    )
    
    omega_candidate <- par_candidate[npar_eta + 1L]
    sig_candidate <- exp(omega_candidate)
    
    grad_eta <- rep(0, npar_eta)
    grad_omega <- 0
    
    if (any(pos)) {
      rpos <- (
        ell[pos] - eta_candidate[pos]
      ) / sig_candidate
      
      grad_eta <- grad_eta +
        as.numeric(
          crossprod(
            X_eta[pos, , drop = FALSE],
            (
              ell[pos] - eta_candidate[pos]
            ) / sig_candidate^2
          )
        )
      
      grad_omega <- grad_omega +
        sum(rpos^2 - 1)
    }
    
    if (any(zero)) {
      t0 <- (
        a_i[zero] - eta_candidate[zero]
      ) / sig_candidate
      
      mills <- lower_mills_ratio(t0)
      
      w0 <- 1 - gamma[zero]
      
      grad_eta <- grad_eta -
        as.numeric(
          crossprod(
            X_eta[zero, , drop = FALSE],
            w0 * mills / sig_candidate
          )
        )
      
      grad_omega <- grad_omega -
        sum(w0 * t0 * mills)
    }
    
    -c(grad_eta, grad_omega)
  }
  
  ll_old <- observed_loglik(
    alpha,
    eta_coef,
    omega
  )
  
  converged <- FALSE
  iter <- 0L
  
  for (m in seq_len(maxit)) {
    iter <- m
    
    ## E-step.
    current_e <- e_step(
      alpha,
      eta_coef,
      omega
    )
    
    gamma <- current_e$gamma
    
    ## Structural-absence M-step.
    rho_fit <- tryCatch(
      optim(
        par = alpha,
        fn = nQ_rho,
        gr = gr_nQ_rho,
        gamma = gamma,
        method = "L-BFGS-B",
        lower = lower_rho,
        upper = upper_rho,
        control = list(
          maxit = 500L,
          factr = 1e7,
          pgtol = 1e-8
        )
      ),
      error = function(e) NULL
    )
    
    rho_ok <-
      !is.null(rho_fit) &&
      rho_fit$convergence == 0L &&
      is.finite(rho_fit$value) &&
      all(is.finite(rho_fit$par))
    
    if (!rho_ok) break
    
    alpha_new <- rho_fit$par
    
    ## Present-conditional M-step.
    eta_start <- c(eta_coef, omega)
    
    eta_fit <- tryCatch(
      optim(
        par = eta_start,
        fn = nQ_eta,
        gr = gr_nQ_eta,
        gamma = gamma,
        method = "L-BFGS-B",
        lower = lower_eta,
        upper = upper_eta,
        control = list(
          maxit = 500L,
          factr = 1e7,
          pgtol = 1e-8
        )
      ),
      error = function(e) NULL
    )
    
    eta_ok <-
      !is.null(eta_fit) &&
      eta_fit$convergence == 0L &&
      is.finite(eta_fit$value) &&
      all(is.finite(eta_fit$par))
    
    if (!eta_ok) break
    
    eta_coef_new <-
      eta_fit$par[seq_len(npar_eta)]
    
    omega_new <-
      eta_fit$par[npar_eta + 1L]
    
    ll_new <- observed_loglik(
      alpha_new,
      eta_coef_new,
      omega_new
    )
    
    ## The observed likelihood should not decrease under an exact EM update.
    if (
      !is.finite(ll_new) ||
      ll_new < ll_old - 1e-7
    ) {
      break
    }
    
    rel_change <-
      abs(ll_new - ll_old) /
      (1 + abs(ll_old))
    
    alpha <- alpha_new
    eta_coef <- eta_coef_new
    omega <- omega_new
    ll_old <- ll_new
    
    if (rel_change < tol) {
      converged <- TRUE
      break
    }
  }
  
  ## Final E-step at the final parameter estimate.
  final_e <- e_step(
    alpha,
    eta_coef,
    omega
  )
  
  em_output_valid <-
    isTRUE(converged) &&
    is.finite(ll_old) &&
    is.finite(final_e$sig) &&
    final_e$sig > 0 &&
    all(is.finite(final_e$rho)) &&
    all(is.finite(final_e$eta)) &&
    all(is.finite(final_e$gamma))
  
  list(
    alpha = alpha,
    eta_coef = eta_coef,
    sig = final_e$sig,
    X_rho = X_rho,
    X_eta = X_eta,
    rho = final_e$rho,
    eta = final_e$eta,
    gamma = final_e$gamma,
    logLik = ll_old,
    valid = em_output_valid,
    converged = converged,
    optimization_succeeded = converged,
    iterations = iter,
    optimizer = "EM",
    method = "em"
  )
}

#' Fit the per-taxon zero model using the requested numerical method
#' @keywords internal
#' @noRd
fit_zero_model <- function(y, N, z = NULL, d0 = 1.0,
                           g = NULL,
                           include_group_rho = FALSE,
                           include_group_eta = FALSE,
                           method = c("optim", "em")) {
  method <- match.arg(method)
  
  if (identical(method, "optim")) {
    return(
      fit_zero_optim(
        y = y,
        N = N,
        z = z,
        d0 = d0,
        g = g,
        include_group_rho = include_group_rho,
        include_group_eta = include_group_eta
      )
    )
  }
  
  fit_zero_em(
    y = y,
    N = N,
    z = z,
    d0 = d0,
    g = g,
    include_group_rho = include_group_rho,
    include_group_eta = include_group_eta
  )
}

#' Prevalence score p-value and null-model posterior structural-zero probabilities
#' @keywords internal
#' @noRd
zero_prev_and_gamma <- function(y, N, g, z = NULL, d0 = 1.0,
                                zero_method = c("optim", "em")) {
  zero_method <- match.arg(zero_method)
  
  y <- as.numeric(y)
  N <- as.numeric(N)
  
  ## Null zero-model fit for the structural-prevalence score.
  ## The group is not included in either the structural-absence logit or
  ## the present-taxon latent abundance mean.
  fit <- fit_zero_model(
    y = y,
    N = N,
    z = z,
    d0 = d0,
    g = g,
    include_group_rho = FALSE,
    include_group_eta = FALSE,
    method = zero_method
  )
  
  rho <- fit$rho
  gamma <- fit$gamma
  
  ## Do not form the prevalence test when the selected fitting method fails
  ## to produce finite fitted zero-model quantities.
  if (!isTRUE(fit$valid)) {
    return(
      list(
        p = 1,
        gamma = gamma,
        rho = rho,
        valid = FALSE,
        converged = isTRUE(fit$converged),
        optimization_succeeded = isTRUE(fit$optimization_succeeded),
        iterations = fit$iterations,
        method = fit$method,
        optimizer = fit$optimizer
      )
    )
  }
  
  ## Depth enters through a_i = log(d0 / N_i), not by residualizing G on logN.
  gp <- residualize_group(g, z = z)
  u <- gp * (gamma - rho)
  U <- sum(u)
  V <- sum(u^2)
  
  p <- if (
    V <= 0 ||
    !is.finite(V) ||
    !is.finite(U)
  ) {
    1
  } else {
    pchisq(U^2 / V, df = 1, lower.tail = FALSE)
  }
  
  list(
    p = p,
    gamma = gamma,
    rho = rho,
    valid = TRUE,
    converged = isTRUE(fit$converged),
    optimization_succeeded = isTRUE(fit$optimization_succeeded),
    iterations = fit$iterations,
    method = fit$method,
    optimizer = fit$optimizer
  )
}

#' Auxiliary posterior structural-zero probabilities and censoring correction
#' @keywords internal
#' @noRd
zero_aux_gamma_censor <- function(y, N, g, z = NULL, d0 = 1.0,
                                  zero_method = c("optim", "em")) {
  zero_method <- match.arg(zero_method)
  
  y <- as.numeric(y)
  N <- as.numeric(N)
  
  ## Auxiliary zero-model fit for abundance weights.
  ## The group is included in the structural-absence logit so that
  ## group-specific structural prevalence is not forced into the abundance arm.
  ## The present-taxon latent abundance mean remains free of the group.
  fit <- fit_zero_model(
    y = y,
    N = N,
    z = z,
    d0 = d0,
    g = g,
    include_group_rho = TRUE,
    include_group_eta = FALSE,
    method = zero_method
  )
  
  rho <- fit$rho
  gamma <- fit$gamma
  eta <- fit$eta
  sig <- fit$sig
  
  pos <- y > 0
  a <- log(d0 / N)
  tval <- (a - eta) / sig
  
  ## Inverse-Mills left-censoring covariate:
  ## E[L | L <= a] = eta - sigma * phi(t) / Phi(t).
  ratio <- lower_mills_ratio(tval)
  censor_corr <- ifelse(pos, 0, sig * ratio)
  censor_corr[!is.finite(censor_corr)] <- 0
  
  list(
    gamma = gamma,
    censor_corr = censor_corr,
    rho = rho,
    valid = isTRUE(fit$valid),
    converged = isTRUE(fit$converged),
    optimization_succeeded = isTRUE(fit$optimization_succeeded),
    iterations = fit$iterations,
    method = fit$method,
    optimizer = fit$optimizer
  )
}

## -----------------------------------------------------------------------------
## 5. Abundance arm: HRIC weighted least squares
## -----------------------------------------------------------------------------

#' Robust affine (LTS) line fit of HRIC group coefficients on root abundance
#' @keywords internal
#' @noRd
## Abundance test: weighted least squares on HRIC coordinates with auxiliary
## posterior-presence weights w_ij = 1 - gamma_ij^aux. Positive counts get weight 1,
## while zeros likely to be structural after allowing group-specific prevalence
## receive small weights. The abundance regression also includes a taxon-specific
## left-censoring correction covariate when it is usable.
## Group effects are corrected by a robust affine fit across taxa, and the final
## uncertainty is computed by a joint sample-level HC3 variance. W is an n x J
## matrix of posterior-presence weights.
##
## Robust affine correction of the raw HRIC group coefficients b_j using
## baseline root abundance v0_j.
##
## An initial LTS fit is used to identify the affine-inlier set H. After H is
## selected, the final
## affine intercept and slope are obtained by ordinary least squares on H:
##
##   a_hat = (G_H^T G_H)^{-1} G_H^T beta_hat_H = B_H beta_hat_H.
##
## The same final refit is used for both the corrected coefficient
##
##   beta_tilde_j = beta_hat_j - g_j^T a_hat
##
## and the joint sample-level HC3 variance. The function therefore returns
## the corrected coefficients, the inlier indices H, and the matrix B_H.
hric_line_fit <- function(b, v0) {
  J <- length(b)
  bc <- rep(NA_real_, J)
  
  ok_idx <- which(
    is.finite(b) &
      is.finite(v0)
  )
  
  invalid_result <- function() {
    list(
      valid = FALSE,
      bc = bc,
      H = integer(0),
      GH = NULL,
      BH = NULL,
      coef = c(NA_real_, NA_real_)
    )
  }
  
  ## The joint affine variance requires a nonsingular
  ## intercept-plus-slope refit.
  if (
    length(ok_idx) < 5L ||
    stats::sd(v0[ok_idx]) < 1e-10
  ) {
    return(invalid_result())
  }
  
  x <- v0[ok_idx]
  y <- b[ok_idx]
  
  initial_fit <- NULL
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package 'MASS' is required for the LTS affine correction."
    )
  }
  
  initial_fit <- tryCatch(
    MASS::lqs(
      y ~ x,
      method = "lts"
    ),
    error = function(e) NULL
  )
  
  if (is.null(initial_fit)) {
    return(invalid_result())
  }
  
  initial_coef <- stats::coef(initial_fit)
  
  chat_init <- unname(initial_coef[1])
  lhat_init <- unname(initial_coef[2])
  
  if (!all(is.finite(c(chat_init, lhat_init)))) {
    return(invalid_result())
  }
  
  initial_residual <- y - (
    chat_init +
      lhat_init * x
  )
  
  ## Robust reweighting used to identify the affine-inlier set H.
  robust_scale <- stats::mad(initial_residual)
  
  if (
    !is.finite(robust_scale) ||
    robust_scale <= 0
  ) {
    robust_scale <- stats::sd(initial_residual)
  }
  
  if (
    !is.finite(robust_scale) ||
    robust_scale <= 0
  ) {
    robust_scale <- 1
  }
  
  H_local <-
    abs(initial_residual) <=
    2.5 * robust_scale
  
  if (sum(H_local) < 4L) {
    return(invalid_result())
  }
  
  H_idx <- ok_idx[H_local]
  
  GH <- cbind(
    Intercept = 1,
    v0 = v0[H_idx]
  )
  
  if (
    nrow(GH) < 3L ||
    qr(GH)$rank < 2L
  ) {
    return(invalid_result())
  }
  
  Sgg <- crossprod(GH)
  
  Sgg_inv <- tryCatch(
    solve(Sgg),
    error = function(e) NULL
  )
  
  if (
    is.null(Sgg_inv) ||
    any(!is.finite(Sgg_inv))
  ) {
    return(invalid_result())
  }
  
  ## Final OLS refit on H. The same BH must be used for both
  ## the corrected point estimate and the joint variance.
  BH <- Sgg_inv %*% t(GH)
  
  coef_final <- as.numeric(
    BH %*% b[H_idx]
  )
  
  G_ok <- cbind(
    1,
    v0[ok_idx]
  )
  
  bc[ok_idx] <-
    b[ok_idx] -
    as.numeric(G_ok %*% coef_final)
  
  list(
    valid = TRUE,
    bc = bc,
    H = H_idx,
    GH = GH,
    BH = BH,
    coef = coef_final
  )
}

#' DASH abundance arm: weighted least squares on HRIC coordinates
#' @keywords internal
#' @noRd
hric_wls_p <- function(Y, g, W, z = NULL, m_abund = 5L, Ccor = NULL) {
  Y <- as.matrix(Y)
  N <- rowSums(Y)
  
  if (!is.null(Ccor)) {
    Ccor <- as.matrix(Ccor)
    
    if (!all(dim(Ccor) == dim(Y))) {
      stop("Ccor must have the same dimensions as Y.")
    }
  }
  
  ## Manuscript version:
  ## HRIC is computed from the observed composition. Exact zeros remain boundary
  ## values of the HRIC coordinate and are not posterior-smoothed.
  M <- hric_transform(Y)
  
  ## Add log library size as a technical nuisance in the abundance arm.
  Z <- make_cov_matrix(z)
  
  raw_logN <- log(pmax(N, 1))
  sd_logN <- stats::sd(raw_logN)
  
  ## Include log library size only when it has non-negligible variation.
  if (
    is.finite(sd_logN) &&
    sd_logN > 1e-12
  ) {
    logN <- as.numeric(
      (raw_logN - mean(raw_logN)) /
        sd_logN
    )
    
    if (is.null(Z)) {
      Z <- matrix(
        logN,
        ncol = 1,
        dimnames = list(NULL, "logN")
      )
    } else {
      Z <- cbind(
        Z,
        logN = logN
      )
    }
  }
  
  X_base <- design_matrix(
    g = g,
    z = Z,
    include_group = TRUE
  )
  
  n <- nrow(M)
  J <- ncol(M)
  
  b <- rep(NA_real_, J)
  
  Phi_raw <- matrix(
    NA_real_,
    nrow = n,
    ncol = J,
    dimnames = list(
      rownames(Y),
      colnames(Y)
    )
  )
  
  for (j in seq_len(J)) {
    ## Abundance eligibility: at least m_abund observed positives in each group.
    ## If this fails, the abundance component is not formed for taxon j.
    n_pos0 <- sum(Y[, j] > 0 & g == 0)
    n_pos1 <- sum(Y[, j] > 0 & g == 1)
    
    if (n_pos0 < m_abund || n_pos1 < m_abund) next
    
    w <- W[, j]
    w[!is.finite(w) | w < 0] <- 0
    w <- pmin(w, 1)
    
    ## Add taxon-specific left-censoring correction as a nuisance covariate.
    ## The correction is used only when it has non-negligible weighted variation
    ## and does not create rank deficiency. Otherwise, fall back to the base
    ## abundance regression rather than dropping the taxon.
    X <- X_base
    
    if (!is.null(Ccor)) {
      c_j <- as.numeric(Ccor[, j])
      c_j[!is.finite(c_j)] <- 0
      
      use_w <- w > 1e-8
      
      if (sum(use_w) > ncol(X_base) + 2) {
        c_bar <- stats::weighted.mean(c_j[use_w], w[use_w])
        c_wsd <- sqrt(stats::weighted.mean((c_j[use_w] - c_bar)^2, w[use_w]))
        
        if (is.finite(c_wsd) && c_wsd > 1e-12) {
          X_try <- cbind(X_base, censorCorr = c_j)
          Xs_try <- X_try * sqrt(w)
          XtWX_try <- crossprod(Xs_try)
          
          if (qr(XtWX_try)$rank == ncol(XtWX_try)) {
            X <- X_try
          }
        }
      }
    }
    
    gidx <- match("Group", colnames(X))
    if (!is.finite(gidx)) next
    
    if (sum(w > 1e-8) <= ncol(X) + 1) next
    
    sqrtw <- sqrt(w)
    
    Xs <- X * sqrtw
    ys <- M[, j] * sqrtw
    
    XtWX <- crossprod(Xs)
    
    rnk <- qr(XtWX)$rank
    if (rnk < ncol(XtWX)) next
    
    XtWXinv <- tryCatch(solve(XtWX), error = function(e) NULL)
    if (is.null(XtWXinv)) next
    
    bhat <- as.numeric(XtWXinv %*% crossprod(Xs, ys))
    e <- as.numeric(M[, j] - X %*% bhat)
    
    ## Weighted HC3 leverage.
    h <- rowSums(
      (Xs %*% XtWXinv) * Xs
    )
    
    h <- pmin(
      pmax(h, 0),
      1 - 1e-8
    )
    
    ## Sample-level HC3 contribution to the raw group coefficient:
    ##
    ## phi_ij =
    ## e_G^T A_j^{-1} x_ij
    ## times
    ## w_ij e_ij / (1 - h_ij).
    group_loading <- as.numeric(
      X %*%
        XtWXinv[, gidx, drop = FALSE]
    )
    
    phi_j <-
      group_loading *
      w *
      e /
      (1 - h)
    
    vv_raw <- sum(phi_j^2)
    
    if (
      !is.finite(vv_raw) ||
      vv_raw <= 0 ||
      any(!is.finite(phi_j))
    ) {
      next
    }
    
    b[j] <- bhat[gidx]
    Phi_raw[, j] <- phi_j
  }
  
  ## Baseline root abundance for the affine HRIC correction.
  P0 <- sweep(Y[g == 0, , drop = FALSE], 1, pmax(N[g == 0], 1), "/")
  v0 <- sqrt(pmax(colMeans(P0), 0))
  
  ## Robust inlier selection followed by the final OLS affine refit.
  lf <- hric_line_fit(
    b = b,
    v0 = v0
  )
  
  p <- rep(1, J)
  tested <- rep(FALSE, J)
  se_joint <- rep(NA_real_, J)
  z_stat <- rep(NA_real_, J)
  
  if (isTRUE(lf$valid)) {
    ## n x |H| matrix:
    ## row i contains sample i's HC3 contributions
    ## to all raw coefficients used in the affine refit.
    phi_H <- Phi_raw[
      ,
      lf$H,
      drop = FALSE
    ]
    
    if (all(is.finite(phi_H))) {
      candidate <- which(
        is.finite(b) &
          is.finite(lf$bc) &
          is.finite(v0)
      )
      
      for (j in candidate) {
        if (any(!is.finite(Phi_raw[, j]))) {
          next
        }
        
        g_j <- c(
          1,
          v0[j]
        )
        
        ## g_j^T B_H: loading of each affine-inlier
        ## raw coefficient on taxon j's fitted background.
        affine_loading_j <- as.numeric(
          t(g_j) %*% lf$BH
        )
        
        ## Sample-level contribution to the fitted affine
        ## background evaluated at v0_j.
        phi_background_j <- as.numeric(
          phi_H %*% affine_loading_j
        )
        
        ## Sample-level contribution to the corrected coefficient:
        ## phi_corr_ij = phi_ij - g_j^T B_H phi_iH.
        phi_corr_j <-
          Phi_raw[, j] -
          phi_background_j
        
        var_joint_j <- sum(
          phi_corr_j^2
        )
        
        if (
          !is.finite(var_joint_j) ||
          var_joint_j <= 0
        ) {
          next
        }
        
        se_joint[j] <- sqrt(var_joint_j)
        
        z_stat[j] <-
          lf$bc[j] /
          se_joint[j]
        
        p[j] <-
          2 * pnorm(
            -abs(z_stat[j])
          )
        
        tested[j] <- is.finite(p[j])
      }
    }
  }
  
  bad_p <- !is.finite(p)
  p[bad_p] <- 1
  tested[bad_p] <- FALSE
  
  names(p) <- colnames(Y)
  names(tested) <- colnames(Y)
  names(se_joint) <- colnames(Y)
  names(z_stat) <- colnames(Y)
  
  list(
    p = p,
    tested = tested,
    beta_raw = b,
    beta_corrected = lf$bc,
    se_joint = se_joint,
    z = z_stat
  )
}

## -----------------------------------------------------------------------------
## 6. DASH: the test
## -----------------------------------------------------------------------------

#' Depth-Aware Structural-zero Hellinger (DASH) differential abundance test
#'
#' DASH is a two-part, per-taxon differential abundance test. The prevalence arm
#' is a structural-zero score test from a per-taxon structural-zero mixture model
#' with a left-censored normal observation component. The zero model may be
#' estimated either by direct maximization of the observed-data likelihood or
#' by EM, with direct optimization used by default. The abundance arm is a
#' weighted least squares test on Hellinger-Riemann intrinsic coordinates (HRIC)
#' with auxiliary posterior-presence weights, a censoring-correction nuisance
#' covariate, a robust affine bias correction, and a joint sample-level HC3
#' standard error for the affine-corrected coefficient.
#' The two arms are combined per taxon by a Bonferroni min-P rule (primary) and by
#' a Cauchy combination (sensitivity). Multiple testing across taxa is the
#' caller's responsibility; apply [stats::p.adjust()] with method `"BH"` to
#' whichever column you report, as in the simulation study.
#'
#' The function applies the simulation preprocessing (retain taxa with at least
#' `min_positive` positive counts, drop empty samples) and then runs the test.
#' The default analysis uses `min_positive = 3`, `d0 = 1`, and
#' `m_abund = 5`. Results also depend on the retained count table, group
#' coding and covariates. MASS is required for the initial LTS fit.
#'
#' @param counts Integer count matrix or data frame, samples (rows) by taxa
#'   (columns). Column names are used as taxon identifiers; if absent they are
#'   set to `tax1`, `tax2`, ...
#' @param group Group indicator with exactly two distinct values. Numeric 0/1 is
#'   used as is. A factor or character vector is mapped to 0/1 by its two sorted
#'   levels, so the second level is treated as the comparison ("case") group.
#' @param covariates Optional measured covariates (vector, matrix, or data
#'   frame), one row per sample, supplied to both arms for adjustment.
#' @param zero_method Numerical method used to estimate each per-taxon
#'   structural-zero mixture model. `"optim"` directly maximizes the
#'   observed-data likelihood using the original bounded L-BFGS-B procedure
#'   with a Nelder-Mead fallback and is the default. `"em"` uses the
#'   expectation-maximization algorithm. Both methods target the same
#'   observed-data likelihood.
#' @param d0 Detection constant for the sample-specific censoring threshold. Default 1.
#' @param m_abund Abundance-eligibility threshold: the abundance arm is formed
#'   for a taxon only if it has at least `m_abund` positive counts in each
#'   group. Default 5.
#' @param min_positive Retain a taxon only if it has at least this many positive
#'   counts across all samples, matching the simulation. Set to 0 to disable
#'   retention and test the supplied taxa as given. Default 3. Samples with no
#'   positive counts are always dropped, since the HRIC transform requires a
#'   positive total in every sample.
#'
#' @return A data frame with one row per retained taxon and columns:
#'   \describe{
#'     \item{taxon}{taxon identifier}
#'     \item{p_prevalence}{prevalence-arm p-value}
#'     \item{prevalence_formed}{whether the prevalence arm was successfully formed}
#'     \item{p_abundance}{abundance-arm p-value (1 when the arm is not formed)}
#'     \item{abundance_formed}{whether the abundance arm was successfully formed}
#'     \item{p_bonf_minp}{DASH (Bonf-minP) omnibus p-value}
#'     \item{p_cauchy}{DASH (Cauchy) omnibus p-value}
#'   }
#'   The retained sample and taxon indices are attached as the attributes
#'   `kept_samples` and `kept_taxa`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80; p <- 40
#' Y <- matrix(rpois(n * p, lambda = 5), n, p)
#' colnames(Y) <- paste0("tax", seq_len(p))
#' grp <- rep(0:1, each = n / 2)
#' res <- dash(Y, grp)
#' res$q_bonf_minp <- p.adjust(res$p_bonf_minp, "BH")
#' head(res)
#' }
#'
#' @export
dash <- function(counts, group, covariates = NULL,
                 d0 = 1, m_abund = 5L, min_positive = 3L,
                 zero_method = c("optim", "em")) {
  zero_method <- match.arg(zero_method)
  
  Y <- as.matrix(counts)
  if (!is.numeric(Y)) {
    stop("`counts` must be a numeric count matrix (samples x taxa).")
  }
  n <- nrow(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("tax", seq_len(ncol(Y)))
  
  ## Group indicator -> 0/1.
  if (is.factor(group) || is.character(group)) {
    gf <- factor(group)
    if (nlevels(gf) != 2L) {
      stop("`group` must have exactly two distinct values.")
    }
    g <- as.integer(gf) - 1L
    group_levels <- levels(gf)
  } else {
    gv <- as.numeric(group)
    u <- sort(unique(gv[is.finite(gv)]))
    if (length(u) != 2L) {
      stop("`group` must take exactly two distinct values.")
    }
    g <- as.integer(gv == u[2L])
    group_levels <- as.character(u)
  }
  if (length(g) != n) {
    stop("`group` length must equal the number of samples (rows of `counts`).")
  }
  
  ## Covariates passed through unchanged; the arms rename the columns internally.
  z <- if (is.null(covariates)) NULL else as.matrix(covariates)
  if (!is.null(z) && nrow(z) != n) {
    stop("`covariates` must have one row per sample.")
  }
  
  ## Feature retention then empty-sample drop.
  keep_tax <- if (min_positive > 0L) {
    colSums(Y > 0) >= min_positive
  } else {
    rep(TRUE, ncol(Y))
  }
  Y <- Y[, keep_tax, drop = FALSE]
  if (ncol(Y) < 2L) {
    stop("Fewer than two taxa remain after retention; cannot run DASH.")
  }
  
  keep_samp <- rowSums(Y) > 0
  Y <- Y[keep_samp, , drop = FALSE]
  g <- g[keep_samp]
  if (!is.null(z)) z <- z[keep_samp, , drop = FALSE]
  if (length(unique(g)) < 2L) {
    stop("Only one group remains after dropping empty samples.")
  }
  
  ## --- Two-part test -------------------------------------------------------
  J <- ncol(Y)
  N <- rowSums(Y)
  taxa_names <- colnames(Y)
  
  ## Null zero-model fit for the structural-prevalence score.
  prev_list <- lapply(seq_len(J), function(j) {
    zero_prev_and_gamma(
      y = as.numeric(Y[, j]),
      N = N,
      g = g,
      z = z,
      d0 = d0,
      zero_method = zero_method
    )
  })
  
  p_pr <- vapply(
    prev_list,
    function(o) o$p,
    numeric(1)
  )
  
  prev_ok <- vapply(
    prev_list,
    function(o) isTRUE(o$valid),
    logical(1)
  )
  
  names(p_pr) <- taxa_names
  names(prev_ok) <- taxa_names
  
  ## Auxiliary zero-model fit for abundance weights and censoring correction.
  ## This fit allows the structural-absence probability to differ by group.
  aux_list <- lapply(seq_len(J), function(j) {
    zero_aux_gamma_censor(
      y = as.numeric(Y[, j]),
      N = N,
      g = g,
      z = z,
      d0 = d0,
      zero_method = zero_method
    )
  })
  
  Gamma_aux <- vapply(
    aux_list,
    function(o) o$gamma,
    numeric(nrow(Y))
  )
  
  Ccor <- vapply(
    aux_list,
    function(o) o$censor_corr,
    numeric(nrow(Y))
  )
  
  aux_ok <- vapply(
    aux_list,
    function(o) isTRUE(o$valid),
    logical(1)
  )
  
  colnames(Gamma_aux) <- taxa_names
  rownames(Gamma_aux) <- rownames(Y)
  
  colnames(Ccor) <- taxa_names
  rownames(Ccor) <- rownames(Y)
  
  ## Presence weights for abundance regression use the auxiliary posterior
  ## structural-zero probabilities, not the null-model probabilities.
  W <- 1 - Gamma_aux
  
  ## A taxon with an invalid auxiliary zero-model fit is not eligible
  ## for the abundance component.
  if (any(!aux_ok)) {
    W[, !aux_ok] <- 0
    Ccor[, !aux_ok] <- 0
  }
  
  ## Abundance component:
  ## observed HRIC + log-depth nuisance + auxiliary presence weights +
  ## censoring-correction covariate + abundance eligibility + joint HC3/Z
  ## inference. Returns the p-values and a logical flag marking the taxa
  ## for which the abundance component was actually formed (eligible + fitted).
  ab <- hric_wls_p(
    Y = Y,
    g = g,
    W = W,
    z = z,
    m_abund = m_abund,
    Ccor = Ccor
  )
  p_ab <- ab$p
  tested_ab <- ab$tested
  names(p_ab) <- taxa_names
  names(tested_ab) <- taxa_names
  
  ## Conservative NA handling.
  p_pr[!is.finite(p_pr)] <- 1
  p_ab[!is.finite(p_ab)] <- 1
  
  ## Omnibus combination:
  ## combine the two component p-values when both components are formed;
  ## otherwise use the single successfully formed component.
  p_cauchy <- vapply(
    seq_len(J),
    function(j) {
      has_prev <- isTRUE(prev_ok[j])
      has_abund <- isTRUE(tested_ab[j])
      
      if (has_prev && has_abund) {
        cauchy_combination(
          c(p_pr[j], p_ab[j])
        )
      } else if (has_prev) {
        unname(p_pr[j])
      } else if (has_abund) {
        unname(p_ab[j])
      } else {
        1
      }
    },
    numeric(1)
  )
  
  p_bonf <- vapply(
    seq_len(J),
    function(j) {
      has_prev <- isTRUE(prev_ok[j])
      has_abund <- isTRUE(tested_ab[j])
      
      if (has_prev && has_abund) {
        min(
          1,
          2 * min(p_pr[j], p_ab[j])
        )
      } else if (has_prev) {
        unname(p_pr[j])
      } else if (has_abund) {
        unname(p_ab[j])
      } else {
        1
      }
    },
    numeric(1)
  )
  
  ## --- Tidy output ---------------------------------------------------------
  out <- data.frame(
    taxon              = taxa_names,
    p_prevalence       = unname(p_pr),
    prevalence_formed  = unname(prev_ok),
    p_abundance        = unname(p_ab),
    abundance_formed   = unname(tested_ab),
    p_bonf_minp        = unname(p_bonf),
    p_cauchy           = unname(p_cauchy),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "zero_method") <- zero_method
  attr(out, "group_levels") <- group_levels
  attr(out, "kept_taxa") <- which(keep_tax)
  attr(out, "kept_samples") <- which(keep_samp)
  out
}
