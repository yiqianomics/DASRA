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
##      zero-inflated Tobit model. Library size enters through a depth-aware
##      censoring threshold, so an observed zero in a deep sample gives stronger
##      evidence of true absence than a zero in a shallow sample.
##    * Abundance arm. A weighted least squares test on Hellinger-Riemann
##      intrinsic coordinates (HRIC), with posterior-presence weights, an HC3
##      sandwich standard error, and a robust affine (LTS) correction for the
##      compositional background whose line-fit variance is propagated into the
##      standard error.
##
##  The two component p-values are combined by a Bonferroni min-P rule (primary)
##  and a Cauchy combination (sensitivity). Multiple testing across taxa is left
##  to the caller; apply Benjamini-Hochberg, stats::p.adjust(., "BH"), to the
##  reported column.
##
##  The numerics are identical to the implementation used in the simulation
##  studies. The defaults d0 = 1, m_abund = 5, and min_positive = 3 reproduce
##  the simulation results exactly.
##
##  Entry points
##    dash()           the test: input coercion, feature retention, both arms,
##                     and the per-taxon combination, returned as a tidy table
##    hric_transform() the HRIC transform, usable on its own
##
##  Dependencies: base R and 'stats'. 'MASS' is optional: when available the
##  affine correction uses MASS::lqs (least trimmed squares); otherwise a
##  base-R Theil-Sen fit is used.
##
##  Contents
##    1. Internal utilities
##    2. P-value combination rules
##    3. HRIC transform
##    4. Prevalence arm: zero-inflated Tobit and the structural-zero score test
##    5. Abundance arm: HRIC weighted least squares with affine bias correction
##    6. DASH: the test
## =============================================================================

#' @importFrom stats optim qlogis plogis pnorm dnorm pchisq pt sd mad median qnorm coef p.adjust
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
## 4. Prevalence arm: zero-inflated Tobit
## -----------------------------------------------------------------------------

#' Fit the per-taxon zero-inflated Tobit null model
#' @keywords internal
#' @noRd
fit_tobit_null <- function(y, N, z = NULL, d0 = 1.0) {
  y <- as.numeric(y); N <- as.numeric(N); pos <- y > 0
  ell <- ifelse(pos, log(pmax(y, 1) / N), 0); cc <- log(d0 / N)
  if (is.null(z)) {
    X <- matrix(1, length(y), 1)
    colnames(X) <- "Intercept"
  } else {
    X <- design_matrix(g = NULL, z = z, include_group = FALSE)
  }
  npar_x <- ncol(X)
  p0 <- min(max(mean(!pos) * 0.5, 0.02), 0.9); init_eta <- if (any(pos)) mean(log(pmax(y[pos] / N[pos], 1e-8))) else -5
  init <- c(c(qlogis(p0), rep(0, npar_x - 1)), c(init_eta, rep(0, npar_x - 1)), log(1.0))
  lower <- c(rep(-10, npar_x), rep(-30, npar_x), log(0.1)); upper <- c(rep(10, npar_x), rep(5, npar_x), log(8.0))
  nll <- function(par) {
    alpha <- par[seq_len(npar_x)]; eta_coef <- par[npar_x + seq_len(npar_x)]; sig <- exp(par[2 * npar_x + 1])
    if (!is.finite(sig) || sig <= 0) return(1e10)
    rho <- clamp(plogis(as.numeric(X %*% alpha)), 1e-12, 1 - 1e-12); eta <- as.numeric(X %*% eta_coef)
    a <- (cc - eta) / sig; Lz <- rho + (1 - rho) * pnorm(a); Lp <- (1 - rho) * dnorm((ell - eta) / sig) / sig; L <- ifelse(pos, Lp, Lz)
    -sum(log(pmax(L, 1e-300)))
  }
  out <- tryCatch(optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 500, factr = 1e7)), error = function(e) NULL)
  par <- if (!is.null(out) && is.finite(out$value)) out$par else tryCatch(optim(init, nll, method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-8))$par, error = function(e) init)
  list(alpha = par[seq_len(npar_x)], eta_coef = par[npar_x + seq_len(npar_x)], sig = max(exp(par[2 * npar_x + 1]), 0.1), X = X)
}

#' Prevalence score p-value and posterior structural-zero probabilities
#' @keywords internal
#' @noRd
## Per-taxon Tobit fit: returns the prevalence score p-value and the posterior
## structural-zero probabilities gamma_ij (0 for positive counts). gamma is reused as
## 1 - gamma for the abundance weights, so the Tobit is fit only once per taxon.


tobit_prev_and_gamma <- function(y, N, g, z = NULL, d0 = 1.0) {
  y <- as.numeric(y)
  N <- as.numeric(N)
  
  fit <- fit_tobit_null(y, N, z = z, d0 = d0)
  X <- fit$X
  
  rho <- clamp(plogis(as.numeric(X %*% fit$alpha)), 1e-10, 1 - 1e-10)
  eta <- as.numeric(X %*% fit$eta_coef)
  sig <- fit$sig
  
  pos <- y > 0
  a <- log(d0 / N)
  
  Phi <- pnorm((a - eta) / sig)
  denom <- pmax(rho + (1 - rho) * Phi, 1e-300)
  
  gamma <- ifelse(pos, 0, rho / denom)
  
  ## Manuscript version:
  ## residualize group only on measured covariates z.
  ## Depth enters the prevalence model through the Tobit censoring threshold a_i.
  gp <- residualize_group(g, z = z)
  
  u <- gp * (gamma - rho)
  
  U <- sum(u)
  V <- sum(u^2)
  
  p <- if (V <= 0 || !is.finite(V)) {
    1
  } else {
    pchisq(U^2 / V, df = 1, lower.tail = FALSE)
  }
  
  list(
    p = p,
    gamma = gamma
  )
}

## -----------------------------------------------------------------------------
## 5. Abundance arm: HRIC weighted least squares
## -----------------------------------------------------------------------------

#' Robust affine (LTS) line fit of HRIC group coefficients on root abundance
#' @keywords internal
#' @noRd
## Abundance test: weighted least squares on HRIC coordinates with posterior-presence
## weights w_ij = 1 - gamma_ij (positives get weight 1, confident structural zeros ~0).
## Group effect uses an HC3 heteroscedasticity-robust sandwich SE, then an affine
## (intercept + slope) bias correction across taxa. W is an n x J matrix of weights.
##
## Robust LTS line fit of the HRIC group coefficients b_j on baseline root abundance
## v0_j. Returns the bias-corrected residual b_j - (chat + lhat v0_j) together with
## the 2x2 covariance Sigma of the line estimate (chat, lhat), computed on the LTS
## inliers H as Sigma = sigma2H * (sum_{k in H} g_k g_k^T)^{-1} with g_k = (1, v0_k)
## and sigma2H = sum_{k in H} r_k^2 / (|H| - 2). This lets the abundance arm add the
## line-fit variance g_j^T Sigma g_j to the HC3 variance of beta-hat when forming the
## standard error of beta-tilde. Uses MASS::lqs when available; falls back to a
## Theil-Sen fit (base R), and to an intercept-only median if the regressor has no
## usable spread. The line is fit on abundance-eligible taxa only (b is NA otherwise);
## differential taxa act as sparse high-residual outliers.
hric_line_fit <- function(b, v0) {
  J <- length(b)
  bc <- rep(NA_real_, J)
  Sigma <- matrix(0, 2, 2)
  ok <- is.finite(b) & is.finite(v0)
  
  ## Intercept-only fallback (too few taxa, or no spread in v0).
  if (sum(ok) < 5 || stats::sd(v0[ok]) < 1e-10) {
    bb <- b[ok]
    med <- median(bb, na.rm = TRUE)
    if (is.finite(med)) {
      bc[ok] <- bb - med
      r <- bb - med
      scl <- stats::mad(r)
      if (!is.finite(scl) || scl <= 0) scl <- stats::sd(r)
      if (!is.finite(scl) || scl <= 0) scl <- 1
      H <- abs(r) <= 2.5 * scl
      if (sum(H) < 3) H <- rep(TRUE, length(r))
      s2 <- sum(r[H]^2) / max(sum(H) - 1, 1)
      Sigma[1, 1] <- s2 / sum(H)   # variance of the common location only
    }
    return(list(bc = bc, Sigma = Sigma))
  }
  
  x <- v0[ok]; y <- b[ok]
  
  fit <- NULL
  if (requireNamespace("MASS", quietly = TRUE)) {
    fit <- tryCatch(MASS::lqs(y ~ x, method = "lts"), error = function(e) NULL)
  }
  
  if (!is.null(fit)) {
    co <- stats::coef(fit)
    chat <- co[1]; lhat <- co[2]
  } else {
    ## Theil-Sen fallback (no dependency).
    dx <- outer(x, x, "-"); dy <- outer(y, y, "-"); ut <- upper.tri(dx)
    sx <- dx[ut]; good <- abs(sx) > 1e-12
    if (!any(good)) {
      med <- median(y)
      bc[ok] <- y - med
      r <- y - med
      s2 <- sum(r^2) / max(length(r) - 1, 1)
      Sigma[1, 1] <- s2 / length(r)
      return(list(bc = bc, Sigma = Sigma))
    }
    lhat <- median(dy[ut][good] / sx[good])
    chat <- median(y - lhat * x)
  }
  
  r_all <- y - (chat + lhat * x)
  bc[ok] <- r_all
  
  ## LTS reweighting: inliers within 2.5 robust scale units.
  scl <- stats::mad(r_all)
  if (!is.finite(scl) || scl <= 0) scl <- stats::sd(r_all)
  if (!is.finite(scl) || scl <= 0) scl <- 1
  H <- abs(r_all) <= 2.5 * scl
  if (sum(H) < 4) H <- rep(TRUE, length(r_all))
  
  ## Sigma = sigma2H * (sum_{k in H} g_k g_k^T)^{-1}, g_k = (1, v0_k).
  gH <- cbind(1, x[H])
  sigma2H <- sum(r_all[H]^2) / max(sum(H) - 2, 1)
  Sgg <- crossprod(gH)
  Sigma <- tryCatch(
    sigma2H * solve(Sgg),
    error = function(e) {
      d <- diag(Sgg); d[d <= 0] <- 1
      sigma2H * diag(1 / d, 2)
    }
  )
  
  list(bc = bc, Sigma = Sigma)
}

#' DASH abundance arm: weighted least squares on HRIC coordinates
#' @keywords internal
#' @noRd
hric_wls_p <- function(Y, g, W, z = NULL, m_abund = 5L) {
  Y <- as.matrix(Y)
  N <- rowSums(Y)
  
  ## Manuscript version:
  ## HRIC is computed from the observed composition. Exact zeros remain boundary
  ## values of the HRIC coordinate and are not posterior-smoothed.
  M <- hric_transform(Y)
  
  ## Add log library size as a technical nuisance in the abundance arm.
  Z <- make_cov_matrix(z)
  logN <- safe_scale(log(pmax(N, 1)))
  
  if (is.null(Z)) {
    Z <- matrix(logN, ncol = 1)
    colnames(Z) <- "logN"
  } else {
    Z <- cbind(Z, logN = logN)
  }
  
  X <- design_matrix(g = g, z = Z, include_group = TRUE)
  
  J <- ncol(M)
  gidx <- match("Group", colnames(X))
  
  b <- rep(NA_real_, J)
  se <- rep(NA_real_, J)
  df <- rep(NA_real_, J)
  
  for (j in seq_len(J)) {
    ## Abundance eligibility: at least m_abund observed positives in each group.
    ## If this fails, the abundance component is not formed for taxon j.
    n_pos0 <- sum(Y[, j] > 0 & g == 0)
    n_pos1 <- sum(Y[, j] > 0 & g == 1)
    
    if (n_pos0 < m_abund || n_pos1 < m_abund) next
    
    w <- W[, j]
    w[!is.finite(w) | w < 0] <- 0
    w <- pmin(w, 1)
    
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
    
    ## HC3 leverage-adjusted sandwich.
    h <- rowSums((Xs %*% XtWXinv) * Xs)
    h <- pmin(pmax(h, 0), 0.99)
    
    meat <- crossprod(X, (w^2 * e^2 / (1 - h)^2) * X)
    V <- XtWXinv %*% meat %*% XtWXinv
    
    vv <- V[gidx, gidx]
    if (!is.finite(vv) || vv <= 0) next
    
    b[j] <- bhat[gidx]
    se[j] <- sqrt(vv)
    df[j] <- max(sum(w > 1e-8) - rnk, 1)
  }
  
  ## Baseline root abundance for the affine HRIC correction.
  P0 <- sweep(Y[g == 0, , drop = FALSE], 1, pmax(N[g == 0], 1), "/")
  v0 <- sqrt(pmax(colMeans(P0), 0))
  
  ## Affine bias correction (LTS line) plus the covariance of the line estimate.
  lf <- hric_line_fit(b, v0)
  bc <- lf$bc
  
  ## Propagate line-fit uncertainty into the SE of beta-tilde:
  ## Var(beta-tilde_j) = se_HC3(beta-hat_j)^2 + g_j^T Sigma g_j, g_j = (1, v0_j).
  G2 <- cbind(1, v0)
  var_line <- rowSums((G2 %*% lf$Sigma) * G2)
  var_line[!is.finite(var_line) | var_line < 0] <- 0
  se_tilde <- sqrt(se^2 + var_line)
  
  p <- rep(1, J)
  
  ## Taxa for which the abundance component was actually formed (eligible + fitted).
  tested <- is.finite(b) & is.finite(se) & se > 0 & is.finite(df) & is.finite(bc)
  
  ## Finite-sample t calibration on the bias-corrected coefficient.
  p[tested] <- 2 * pt(abs(bc[tested] / se_tilde[tested]),
                      df = df[tested], lower.tail = FALSE)
  
  p[!is.finite(p)] <- 1
  names(p) <- colnames(Y)
  names(tested) <- colnames(Y)
  
  list(p = p, tested = tested)
}

## -----------------------------------------------------------------------------
## 6. DASH: the test
## -----------------------------------------------------------------------------

#' Depth-Aware Structural-zero Hellinger (DASH) differential abundance test
#'
#' DASH is a two-part, per-taxon differential abundance test. The prevalence arm
#' is a structural-zero score test from a per-taxon zero-inflated Tobit model.
#' The abundance arm is a weighted least squares test on Hellinger-Riemann
#' intrinsic coordinates (HRIC) with an HC3 sandwich standard error and a robust
#' affine bias correction. The two arms are combined per taxon by a Bonferroni
#' min-P rule (primary) and by a Cauchy combination (sensitivity). Multiple
#' testing across taxa is the caller's responsibility; apply
#' [stats::p.adjust()] with method `"BH"` to whichever column you report, as in
#' the simulation study.
#'
#' The function applies the simulation preprocessing (retain taxa with at least
#' `min_positive` positive counts, drop empty samples) and then runs the test.
#' With `min_positive = 3`, `d0 = 1`, and `m_abund = 5` the per-taxon p-values
#' reproduce the simulation results exactly.
#'
#' @param counts Integer count matrix or data frame, samples (rows) by taxa
#'   (columns). Column names are used as taxon identifiers; if absent they are
#'   set to `tax1`, `tax2`, ...
#' @param group Group indicator with exactly two distinct values. Numeric 0/1 is
#'   used as is. A factor or character vector is mapped to 0/1 by its two sorted
#'   levels, so the second level is treated as the comparison ("case") group.
#' @param covariates Optional measured covariates (vector, matrix, or data
#'   frame), one row per sample, supplied to both arms for adjustment.
#' @param d0 Detection constant for the Tobit censoring threshold. Default 1.
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
#'     \item{p_abundance}{abundance-arm p-value (1 when the arm is not formed)}
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
                 d0 = 1, m_abund = 5L, min_positive = 3L) {
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
  
  ## Fit Tobit once per taxon:
  ## prevalence p-value + posterior structural-zero gamma.
  prev_list <- lapply(seq_len(J), function(j) {
    tobit_prev_and_gamma(
      y = as.numeric(Y[, j]),
      N = N,
      g = g,
      z = z,
      d0 = d0
    )
  })
  
  p_pr <- vapply(prev_list, function(o) o$p, numeric(1))
  names(p_pr) <- taxa_names
  
  Gamma <- vapply(prev_list, function(o) o$gamma, numeric(nrow(Y)))
  
  colnames(Gamma) <- taxa_names
  rownames(Gamma) <- rownames(Y)
  
  ## Presence weights for abundance regression.
  W <- 1 - Gamma
  
  ## Abundance component:
  ## observed HRIC + log-depth nuisance + abundance eligibility + HC3/t with
  ## line-fit variance. Returns the p-values and a logical flag marking the taxa
  ## for which the abundance component was actually formed (eligible + fitted).
  ab <- hric_wls_p(
    Y = Y,
    g = g,
    W = W,
    z = z,
    m_abund = m_abund
  )
  p_ab <- ab$p
  tested_ab <- ab$tested
  names(p_ab) <- taxa_names
  names(tested_ab) <- taxa_names
  
  ## Conservative NA handling.
  p_pr[!is.finite(p_pr)] <- 1
  p_ab[!is.finite(p_ab)] <- 1
  
  ## Omnibus combination:
  ## when both components are formed, combine the two;
  ## when only the prevalence component is formed there is no second test and no
  ## multiplicity, so the omnibus p-value equals p_prev (no factor of two).
  p_cauchy <- vapply(seq_len(J), function(j) {
    if (isTRUE(unname(tested_ab[j]))) cauchy_combination(c(p_pr[j], p_ab[j])) else unname(p_pr[j])
  }, numeric(1))
  
  p_bonf <- vapply(seq_len(J), function(j) {
    if (isTRUE(unname(tested_ab[j]))) min(1, 2 * min(p_pr[j], p_ab[j])) else unname(p_pr[j])
  }, numeric(1))
  
  ## --- Tidy output ---------------------------------------------------------
  out <- data.frame(
    taxon        = taxa_names,
    p_prevalence = unname(p_pr),
    p_abundance  = unname(p_ab),
    p_bonf_minp  = unname(p_bonf),
    p_cauchy     = unname(p_cauchy),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "group_levels") <- group_levels
  attr(out, "kept_taxa") <- which(keep_tax)
  attr(out, "kept_samples") <- which(keep_samp)
  out
}
