# DASRA: Depth-Aware Structural-Absence and Relative-Abundance Analysis

#' DASRA: Depth-Aware Structural-Absence and Relative-Abundance Analysis
#'
#' DASRA provides depth-aware inference for structural-absence and
#' relative-abundance associations in microbiome count data. It reports the
#' two components separately and combines them in an omnibus analysis.
#'
#' @keywords internal
#' @import stats
#' @importFrom Rcpp evalCpp
#' @useDynLib DASRA, .registration = TRUE
#' @noRd
"_PACKAGE"

# --------------------------------------------------------------------------
# Internal numerical engine
# --------------------------------------------------------------------------

count_clamp <- function(x, lo, hi) {
    pmin(pmax(x, lo), hi)
}

count_softplus <- function(x) {
    pmax(x, 0) + log1p(exp(-abs(x)))
}

count_logspace_add <- function(a, b) {
    m <- pmax(a, b)
    out <- m + log(exp(a - m) + exp(b - m))
    both_inf <- !is.finite(m)
    out[both_inf] <- m[both_inf]
    out
}

count_row_log_sum_exp <- function(x) {
    x <- as.matrix(x)
    m <- apply(x, 1L, max)
    out <- m + log(rowSums(exp(x - m)))
    bad <- !is.finite(m)
    out[bad] <- m[bad]
    out
}

count_numeric_gradient <- function(par, fn,
                                   lower = rep(-Inf, length(par)),
                                   upper = rep(Inf, length(par))) {
    grad <- numeric(length(par))
    f0 <- fn(par)
    for (k in seq_along(par)) {
        h <- 1e-05 * (1 + abs(par[k]))
        can_left <- par[k] - h >= lower[k]
        can_right <- par[k] + h <= upper[k]
        if (can_left && can_right) {
            p_left <- p_right <- par
            p_left[k] <- par[k] - h
            p_right[k] <- par[k] + h
            grad[k] <- (fn(p_right) - fn(p_left)) / (2 * h)
        } else if (can_right) {
            p_right <- par
            p_right[k] <- par[k] + h
            grad[k] <- (fn(p_right) - f0) / h
        } else if (can_left) {
            p_left <- par
            p_left[k] <- par[k] - h
            grad[k] <- (f0 - fn(p_left)) / h
        } else {
            grad[k] <- NA
        }
    }
    grad
}

.dasra_count_gh_cache <- new.env(parent = emptyenv())

make_count_gh_rule <- function(Q = 21L) {
    Q <- as.integer(Q)
    if (length(Q) != 1L || !is.finite(Q) || Q < 3L) {
        stop("Q must be one integer at least 3.")
    }
    cache_key <- as.character(Q)
    if (exists(cache_key, envir = .dasra_count_gh_cache,
               inherits = FALSE)) {
        return(get(cache_key, envir = .dasra_count_gh_cache,
                   inherits = FALSE))
    }
    J <- matrix(0, Q, Q)
    if (Q > 1L) {
        off <- sqrt(seq_len(Q - 1L) / 2)
        J[cbind(seq_len(Q - 1L), 2:Q)] <- off
        J[cbind(2:Q, seq_len(Q - 1L))] <- off
    }
    ee <- eigen(J, symmetric = TRUE)
    ord <- order(ee$values)
    node <- ee$values[ord]
    weight <- ee$vectors[1L, ord]^2
    weight <- weight / sum(weight)
    rule <- list(
        Q = Q,
        node = node,
        weight = weight,
        log_weight = log(weight),
        log_raw_weight = log(weight) + 0.5 * log(pi)
    )
    assign(cache_key, rule, envir = .dasra_count_gh_cache)
    rule
}

.dasra_structural_gh_cache <- new.env(parent = emptyenv())

make_structural_gh_rule <- function(Q = 1001L) {
    Q <- as.integer(Q)
    if (length(Q) != 1L || !is.finite(Q) || Q < 3L) {
        stop("Q must be one integer at least 3.")
    }
    cache_key <- as.character(Q)
    if (exists(cache_key, envir = .dasra_structural_gh_cache,
               inherits = FALSE)) {
        return(get(cache_key, envir = .dasra_structural_gh_cache,
                   inherits = FALSE))
    }
    if (!requireNamespace("statmod", quietly = TRUE)) {
        stop("Package 'statmod' is required for structural quadrature.")
    }
    package_rule <- statmod::gauss.quad(Q, kind = "hermite")
    node <- as.numeric(package_rule$nodes)
    if (length(node) != Q || any(!is.finite(node)) ||
        is.unsorted(node, strictly = TRUE)) {
        stop("The Gauss-Hermite node rule is invalid.")
    }
    log_raw_weight <- dasra_gh_log_weights_cpp(node, Q)
    if (length(log_raw_weight) != Q || any(!is.finite(log_raw_weight))) {
        stop("The Gauss-Hermite log weights are invalid.")
    }
    log_weight <- log_raw_weight - 0.5 * log(pi)
    rule <- list(
        Q = Q,
        node = node,
        weight = exp(log_weight),
        log_weight = log_weight,
        log_raw_weight = log_raw_weight,
        provider = "statmod::gauss.quad",
        weight_scale = "log-domain Hermite recurrence"
    )
    assign(cache_key, rule, envir = .dasra_structural_gh_cache)
    rule
}

count_latent_mode <- function(y, N, eta, sigma, iterations = 80L,
                              score_tolerance = 1e-12,
                              bracket_tolerance = 1e-12) {
    y <- as.numeric(y)
    N <- as.numeric(N)
    eta <- as.numeric(eta)
    if (length(N) != length(y) || length(eta) != length(y) ||
        length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("Invalid inputs to count_latent_mode.")
    }
    prior_precision <- 1 / sigma^2
    pseudo_p <- count_clamp((y + 0.5) / (N + 1), 1e-12, 1 - 1e-12)
    x <- (qlogis(pseudo_p) + prior_precision * eta) / (1 + prior_precision)
    lower <- eta - sigma^2 * (N - y)
    upper <- eta + sigma^2 * y
    x <- pmin(pmax(x, lower), upper)
    for (iter in seq_len(as.integer(iterations))) {
        p <- plogis(x)
        score <- y - N * p - (x - eta) * prior_precision
        curvature <- N * p * (1 - p) + prior_precision
        positive <- score > 0
        lower[positive] <- x[positive]
        upper[!positive] <- x[!positive]
        newton <- x + score / curvature
        midpoint <- lower + 0.5 * (upper - lower)
        bracket_width <- upper - lower
        use_newton <- is.finite(newton) & newton > lower & newton < upper &
            abs(newton - x) <= 0.5 * bracket_width
        next_x <- ifelse(use_newton, newton, midpoint)
        newton_step_size <- abs(score / curvature)
        converged <- newton_step_size <= score_tolerance * (1 + abs(x)) |
            abs(upper - lower) <= bracket_tolerance * (1 + abs(x))
        next_x[converged] <- x[converged]
        x <- next_x
        if (all(converged)) break
    }
    p <- plogis(x)
    score <- y - N * p - (x - eta) * prior_precision
    curvature <- N * p * (1 - p) + prior_precision
    newton_step_size <- abs(score / curvature)
    if (any(!is.finite(x)) || any(!is.finite(curvature)) ||
        any(curvature <= 0) ||
        any(newton_step_size > 1e-08 * (1 + abs(x)))) {
        stop("Safeguarded latent-mode solver failed to converge.")
    }
    list(mode = x, curvature = curvature, score = score,
         bracket_width = upper - lower)
}

count_log_hy_adaptive <- function(y, N, eta, sigma, gh,
                                  return_nodes = FALSE) {
    y <- as.numeric(y)
    N <- as.numeric(N)
    eta <- as.numeric(eta)
    n <- length(y)
    if (length(N) != n || length(eta) != n) {
        stop("y, N, and eta must have the same length.")
    }
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be finite and positive.")
    }
    if (any(!is.finite(y)) || any(!is.finite(N)) || any(y < 0) ||
        any(N < 0) || any(y > N) ||
        any(abs(y - round(y)) > 1e-08) ||
        any(abs(N - round(N)) > 1e-08)) {
        stop("y and N must be finite integer counts satisfying 0 <= y <= N.")
    }
    mode_info <- count_latent_mode(y, N, eta, sigma)
    mode <- mode_info$mode
    curvature <- mode_info$curvature
    Q <- length(gh$node)
    scale <- sqrt(2 / curvature)
    x_node <- matrix(mode, nrow = n, ncol = Q) +
        matrix(scale, nrow = n, ncol = Q) *
        matrix(gh$node, nrow = n, ncol = Q, byrow = TRUE)
    log_kernel <- sweep(x_node, 1L, y, "*") -
        sweep(count_softplus(x_node), 1L, N, "*") -
        sweep((x_node - matrix(eta, nrow = n, ncol = Q))^2, 1L,
              rep(2 * sigma^2, n), "/") -
        log(sigma) - 0.5 * log(2 * pi)
    log_choose <- lgamma(N + 1) - lgamma(y + 1) - lgamma(N - y + 1)
    log_terms <- log_kernel +
        matrix(gh$node^2, nrow = n, ncol = Q, byrow = TRUE) +
        matrix(gh$log_raw_weight, nrow = n, ncol = Q, byrow = TRUE)
    log_hy <- log_choose + 0.5 * log(2 / curvature) +
        count_row_log_sum_exp(log_terms)
    if (!return_nodes) {
        return(log_hy)
    }
    list(log_hy = log_hy, log_terms = log_terms,
         log_term_normalizer = count_row_log_sum_exp(log_terms),
         x_node = x_node, mode = mode, curvature = curvature)
}

count_log_hy_adaptive_R <- count_log_hy_adaptive

count_log_hy_adaptive <- function(y, N, eta, sigma, gh,
                                  return_nodes = FALSE) {
    if (isTRUE(return_nodes)) {
        return(count_log_hy_adaptive_R(
            y = y, N = N, eta = eta, sigma = sigma, gh = gh,
            return_nodes = TRUE
        ))
    }
    dasra_count_log_hy_adaptive_cpp(
        y = as.numeric(y),
        N = as.numeric(N),
        eta = as.numeric(eta),
        sigma = as.numeric(sigma),
        node = as.numeric(gh$node),
        log_raw_weight = as.numeric(gh$log_raw_weight),
        iterations = 80L,
        score_tolerance = 1e-12,
        bracket_tolerance = 1e-12
    )
}

count_log_hy <- function(y, N, eta, sigma, gh, return_nodes = FALSE) {
    y <- as.numeric(y)
    N <- as.numeric(N)
    eta <- as.numeric(eta)
    n <- length(y)
    if (length(N) != n || length(eta) != n) {
        stop("y, N, and eta must have the same length.")
    }
    if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
        stop("sigma must be finite and positive.")
    }
    if (any(!is.finite(y)) || any(!is.finite(N)) || any(y < 0) ||
        any(N < 0) || any(y > N) ||
        any(abs(y - round(y)) > 1e-08) ||
        any(abs(N - round(N)) > 1e-08)) {
        stop("y and N must be finite integer counts satisfying 0 <= y <= N.")
    }
    Q <- length(gh$node)
    x_node <- matrix(eta, nrow = n, ncol = Q) +
        sqrt(2) * sigma * matrix(gh$node, nrow = n, ncol = Q, byrow = TRUE)
    log_p <- -count_softplus(-x_node)
    log_one_minus_p <- -count_softplus(x_node)
    log_kernel <- sweep(log_p, 1L, y, "*") +
        sweep(log_one_minus_p, 1L, N - y, "*")
    log_choose <- lgamma(N + 1) - lgamma(y + 1) - lgamma(N - y + 1)
    log_terms <- sweep(log_kernel, 1L, log_choose, "+") +
        matrix(gh$log_weight, nrow = n, ncol = Q, byrow = TRUE)
    log_hy <- count_row_log_sum_exp(log_terms)
    if (!return_nodes) {
        return(log_hy)
    }
    list(log_hy = log_hy, log_terms = log_terms,
         log_term_normalizer = log_hy, x_node = x_node)
}

count_design_matrix <- function(g = NULL, z = NULL, include_group = FALSE,
                                n = NULL) {
    if (is.null(n)) {
        if (!is.null(g)) {
            n <- length(g)
        } else if (!is.null(z)) {
            n <- nrow(as.matrix(z))
        } else {
            stop("n is required when both g and z are NULL.")
        }
    }
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "Intercept"
    if (include_group) {
        if (is.null(g) || length(g) != n) {
            stop("g must be supplied with one entry per sample when group is included.")
        }
        X <- cbind(X, Group = as.numeric(g))
    }
    if (!is.null(z)) {
        Z <- as.matrix(z)
        if (nrow(Z) != n) {
            stop("z must have one row per sample.")
        }
        if (ncol(Z) > 0L) {
            colnames(Z) <- paste0("z", seq_len(ncol(Z)))
            X <- cbind(X, Z)
        }
    }
    X
}

# Equal-weight Cauchy combination used for the sensitivity omnibus.
cauchy_combination <- function(ps) {
    ps <- ps[is.finite(ps)]
    if (!length(ps)) return(NA_real_)
    ps <- pmin(pmax(ps, 1e-15), 1 - 1e-15)
    0.5 - atan(mean(tan((0.5 - ps) * pi))) / pi
}

# --------------------------------------------------------------------------
# Relative-abundance engine
# --------------------------------------------------------------------------

.dasra_abundance_control <- function() {
    list(
        quadrature_Q = 41L,
        effect_quadrature_Q = 41L,
        root_tolerance = 1e-5,
        parameter_boundary_tolerance = 1e-5,
        jacobian_condition_limit = 1e10
    )
}

.dasra_abundance_log1mexp <- function(log_x) {
    log_x <- as.numeric(log_x)
    if (any(log_x > 1e-12, na.rm = TRUE)) {
        stop("The log probability must not exceed zero.")
    }
    log_x <- pmin(log_x, 0)
    cutoff <- -log(2)
    output <- numeric(length(log_x))
    near_zero <- log_x > cutoff
    output[near_zero] <- log(-expm1(log_x[near_zero]))
    output[!near_zero] <- log1p(-exp(log_x[!near_zero]))
    output
}

.dasra_abundance_numeric_gradient <- function(
        par, fn, lower = rep(-Inf, length(par)),
        upper = rep(Inf, length(par))) {
    par <- as.numeric(par)
    gradient <- numeric(length(par))
    f0 <- fn(par)
    for (k in seq_along(par)) {
        step <- 1e-4 * (1 + abs(par[k]))
        left_ok <- par[k] - step > lower[k]
        right_ok <- par[k] + step < upper[k]
        if (left_ok && right_ok) {
            left <- right <- par
            left[k] <- par[k] - step
            right[k] <- par[k] + step
            gradient[k] <- (fn(right) - fn(left)) / (2 * step)
        } else if (right_ok) {
            right <- par
            right[k] <- par[k] + step
            gradient[k] <- (fn(right) - f0) / step
        } else if (left_ok) {
            left <- par
            left[k] <- par[k] - step
            gradient[k] <- (f0 - fn(left)) / step
        } else {
            gradient[k] <- NA_real_
        }
    }
    gradient
}

.dasra_abundance_numeric_jacobian <- function(
        par, fn, lower = rep(-Inf, length(par)),
        upper = rep(Inf, length(par))) {
    par <- as.numeric(par)
    f0 <- as.numeric(fn(par))
    jacobian <- matrix(NA_real_, nrow = length(f0), ncol = length(par))
    for (k in seq_along(par)) {
        step <- 1e-4 * (1 + abs(par[k]))
        left_ok <- par[k] - step > lower[k]
        right_ok <- par[k] + step < upper[k]
        if (left_ok && right_ok) {
            left <- right <- par
            left[k] <- par[k] - step
            right[k] <- par[k] + step
            jacobian[, k] <- (fn(right) - fn(left)) / (2 * step)
        } else if (right_ok) {
            right <- par
            right[k] <- par[k] + step
            jacobian[, k] <- (fn(right) - f0) / step
        } else if (left_ok) {
            left <- par
            left[k] <- par[k] - step
            jacobian[, k] <- (f0 - fn(left)) / step
        }
    }
    jacobian
}

.dasra_abundance_condition_number <- function(x) {
    singular_values <- tryCatch(
        svd(x, nu = 0, nv = 0)$d,
        error = function(e) numeric()
    )
    if (!length(singular_values) || any(!is.finite(singular_values)) ||
        min(singular_values) <= 0) {
        return(Inf)
    }
    max(singular_values) / min(singular_values)
}

.dasra_abundance_marginal <- function(
        y, N, eta, sigma, gh, need_moments = TRUE) {
    dasra_count_moments_adaptive_cpp(
        y = as.numeric(y),
        N = as.numeric(N),
        eta = as.numeric(eta),
        sigma = as.numeric(sigma),
        node = as.numeric(gh$node),
        log_raw_weight = as.numeric(gh$log_raw_weight),
        need_moments = isTRUE(need_moments),
        iterations = 80L,
        score_tolerance = 1e-12,
        bracket_tolerance = 1e-12
    )
}

.dasra_abundance_mean_log_relative <- function(location, sigma, gh) {
    dasra_mean_log_relative_cpp(
        location = as.numeric(location),
        sigma = as.numeric(sigma),
        z = sqrt(2) * as.numeric(gh$node),
        weight = as.numeric(gh$weight)
    )
}

.dasra_abundance_designs <- function(group, z) {
    group <- as.numeric(group)
    n <- length(group)
    if (is.null(z)) {
        z <- matrix(numeric(), nrow = n, ncol = 0L)
    } else {
        z <- as.matrix(z)
    }
    if (nrow(z) != n) {
        stop("The abundance covariate matrix has the wrong number of rows.")
    }
    if (!ncol(z)) {
        X_b <- matrix(1, nrow = n, ncol = 1L,
                      dimnames = list(NULL, "Intercept"))
        X_rho <- cbind(Intercept = 1, Group = group)
    } else {
        X_b <- cbind(Intercept = 1, z)
        colnames(X_b) <- c("Intercept", paste0("z", seq_len(ncol(z))))
        X_rho <- cbind(Intercept = 1, Group = group, z)
        colnames(X_rho) <- c(
            "Intercept", "Group", paste0("z", seq_len(ncol(z)))
        )
    }
    list(X_b = X_b, X_rho = X_rho, X_full = X_rho)
}

.dasra_abundance_layout <- function(X_b, X_rho) {
    nb <- ncol(X_b)
    na <- ncol(X_rho)
    list(
        nb = nb,
        na = na,
        b = seq_len(nb),
        omega = nb + 1L,
        a = nb + 1L + seq_len(na),
        zeta = nb + na + 2L,
        dimension = nb + na + 2L
    )
}

.dasra_abundance_bounds <- function(layout) {
    lower <- upper <- rep(NA_real_, layout$dimension)
    lower[layout$b] <- c(-20, rep(-5, layout$nb - 1L))
    upper[layout$b] <- c(5, rep(5, layout$nb - 1L))
    lower[layout$omega] <- log(0.15)
    upper[layout$omega] <- log(4)
    lower[layout$a] <- -10
    upper[layout$a] <- 10
    lower[layout$zeta] <- -5
    upper[layout$zeta] <- 5
    list(lower = lower, upper = upper)
}

.dasra_abundance_state <- function(theta, data, gh, return_psi = TRUE) {
    layout <- data$layout
    b <- theta[layout$b]
    omega <- theta[layout$omega]
    a <- theta[layout$a]
    zeta <- theta[layout$zeta]

    sigma <- exp(omega)
    eta <- as.numeric(data$X_b %*% b + zeta * data$group)
    lp_rho <- as.numeric(data$X_rho %*% a)
    rho <- count_clamp(plogis(lp_rho), 1e-10, 1 - 1e-10)

    marginal_y <- .dasra_abundance_marginal(
        y = data$y, N = data$N, eta = eta, sigma = sigma,
        gh = gh, need_moments = TRUE
    )
    marginal_0 <- .dasra_abundance_marginal(
        y = rep(0, length(data$y)), N = data$N, eta = eta,
        sigma = sigma, gh = gh, need_moments = TRUE
    )

    log_h0 <- pmin(marginal_0$log_h, -1e-12)
    h0 <- count_clamp(exp(log_h0), 0, 1 - 1e-12)
    r <- 1 - h0
    if (any(!is.finite(r)) || any(r <= 1e-12)) {
        stop("The conditional detection probability is numerically degenerate.")
    }

    positive <- data$y > 0
    gamma <- numeric(length(data$y))
    zero <- !positive
    if (any(zero)) {
        gamma[zero] <- plogis(qlogis(rho[zero]) - log_h0[zero])
    }
    presence_weight <- 1 - gamma

    base <- list(
        eta = eta,
        sigma = sigma,
        rho = rho,
        gamma = gamma,
        presence_weight = presence_weight,
        h0 = h0,
        r = r,
        marginal_y = marginal_y,
        marginal_0 = marginal_0
    )
    if (!isTRUE(return_psi)) return(base)

    n <- length(data$y)
    mark_eta <- numeric(n)
    mark_omega <- numeric(n)
    if (any(positive)) {
        ratio <- h0[positive] / r[positive]
        mark_eta[positive] <- marginal_y$score_eta[positive] +
            ratio * marginal_0$score_eta[positive]
        mark_omega[positive] <- marginal_y$score_omega[positive] +
            ratio * marginal_0$score_omega[positive]
    }

    psi_b <- data$X_b * mark_eta
    psi_omega <- mark_omega
    psi_a <- data$X_rho * (gamma - rho)
    psi_zeta <- data$group * presence_weight * marginal_y$score_eta
    base$psi <- cbind(psi_b, omega = psi_omega, psi_a, zeta = psi_zeta)
    base
}

.dasra_abundance_initial_mark <- function(data, gh) {
    positive <- data$y > 0
    if (sum(positive) <= ncol(data$X_full) + 1L ||
        qr(data$X_full[positive, , drop = FALSE])$rank <
            ncol(data$X_full)) {
        return(NULL)
    }

    X <- data$X_full
    y <- data$y
    N <- data$N
    objective <- function(par) {
        coefficient <- par[seq_len(ncol(X))]
        sigma <- exp(par[ncol(X) + 1L])
        eta <- as.numeric(X %*% coefficient)
        marginal_y <- .dasra_abundance_marginal(
            y[positive], N[positive], eta[positive], sigma, gh,
            need_moments = FALSE
        )
        marginal_0 <- .dasra_abundance_marginal(
            rep(0, sum(positive)), N[positive], eta[positive], sigma, gh,
            need_moments = FALSE
        )
        if (any(marginal_0$log_h >= -1e-12)) return(1e8)
        log_r <- .dasra_abundance_log1mexp(marginal_0$log_h)
        value <- -mean(marginal_y$log_h - log_r)
        if (is.finite(value)) value else 1e8
    }

    observed_logit <- qlogis(count_clamp(
        (y[positive] + 0.5) / (N[positive] + 1),
        1e-8, 1 - 1e-8
    ))
    coefficient_start <- tryCatch(
        as.numeric(qr.solve(
            data$X_full[positive, , drop = FALSE], observed_logit
        )),
        error = function(e) rep(0, ncol(data$X_full))
    )
    coefficient_start[!is.finite(coefficient_start)] <- 0
    coefficient_start[1L] <- count_clamp(coefficient_start[1L], -20, 5)
    if (length(coefficient_start) > 1L) {
        coefficient_start[-1L] <- count_clamp(
            coefficient_start[-1L], -5, 5
        )
    }

    residual <- observed_logit -
        as.numeric(data$X_full[positive, , drop = FALSE] %*%
                       coefficient_start)
    sigma_start <- count_clamp(stats::sd(residual), 0.3, 2)
    if (!is.finite(sigma_start)) sigma_start <- 1

    starts <- lapply(c(0.5, sigma_start, 1, 1.5), function(s) {
        c(coefficient_start, log(s))
    })
    starts[[length(starts) + 1L]] <- c(
        c(mean(observed_logit), rep(0, ncol(X) - 1L)), log(1)
    )
    lower <- c(-20, rep(-5, ncol(X) - 1L), log(0.15))
    upper <- c(5, rep(5, ncol(X) - 1L), log(4))

    fits <- lapply(starts, function(start) {
        start <- count_clamp(start, lower, upper)
        tryCatch(
            optim(
                start, objective, method = "L-BFGS-B",
                lower = lower, upper = upper,
                control = list(maxit = 300L, factr = 1e7)
            ),
            error = function(e) NULL
        )
    })
    valid <- vapply(fits, function(fit) {
        !is.null(fit) && is.finite(fit$value) && all(is.finite(fit$par))
    }, logical(1))
    if (!any(valid)) return(NULL)
    candidates <- fits[valid]
    candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
}

.dasra_abundance_initial_detection <- function(
        data, b, omega, zeta, gh) {
    sigma <- exp(omega)
    eta <- as.numeric(data$X_b %*% b + zeta * data$group)
    marginal_0 <- .dasra_abundance_marginal(
        rep(0, length(data$y)), data$N, eta, sigma, gh,
        need_moments = FALSE
    )
    if (any(marginal_0$log_h >= -1e-12)) return(NULL)
    r <- 1 - exp(marginal_0$log_h)
    detected <- as.numeric(data$y > 0)

    objective <- function(a) {
        rho <- count_clamp(
            plogis(as.numeric(data$X_rho %*% a)), 1e-10, 1 - 1e-10
        )
        q <- count_clamp((1 - rho) * r, 1e-12, 1 - 1e-12)
        value <- -mean(detected * log(q) + (1 - detected) * log1p(-q))
        if (is.finite(value)) value else 1e8
    }

    observed_zero <- mean(data$y == 0)
    starts <- lapply(c(0.10, 0.30, 0.50, 0.70), function(multiplier) {
        probability <- count_clamp(
            multiplier * max(observed_zero, 0.05), 0.01, 0.90
        )
        c(qlogis(probability), rep(0, ncol(data$X_rho) - 1L))
    })
    fits <- lapply(starts, function(start) {
        tryCatch(
            optim(
                start, objective, method = "L-BFGS-B",
                lower = rep(-10, length(start)),
                upper = rep(10, length(start)),
                control = list(maxit = 300L, factr = 1e7)
            ),
            error = function(e) NULL
        )
    })
    valid <- vapply(fits, function(fit) {
        !is.null(fit) && is.finite(fit$value) && all(is.finite(fit$par))
    }, logical(1))
    if (!any(valid)) return(NULL)
    candidates <- fits[valid]
    candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
}

.dasra_abundance_effect <- function(theta, data, gh_effect) {
    layout <- data$layout
    b <- theta[layout$b]
    sigma <- exp(theta[layout$omega])
    zeta <- theta[layout$zeta]
    baseline <- as.numeric(data$X_b %*% b)
    mean(
        .dasra_abundance_mean_log_relative(
            baseline + zeta, sigma, gh_effect
        ) - .dasra_abundance_mean_log_relative(
            baseline, sigma, gh_effect
        )
    )
}

.dasra_abundance_empty_fit <- function(status, n = 0L) {
    list(
        available = FALSE,
        status = as.character(status)[1L],
        theta = NULL,
        raw_delta = NA_real_,
        raw_se = NA_real_,
        raw_p = NA_real_,
        phi = if (n > 0L) rep(NA_real_, n) else numeric(),
        score_residue = NA_real_,
        jacobian_condition = NA_real_,
        mean_presence_weight = NA_real_,
        mean_zero_presence_weight = NA_real_
    )
}

.dasra_abundance_fit_taxon <- function(
        y, N, group, z, gh_fit, gh_effect, control) {
    y <- as.numeric(y)
    N <- as.numeric(N)
    group <- as.numeric(group)
    n <- length(y)
    if (is.null(z)) {
        z <- matrix(numeric(), nrow = n, ncol = 0L)
    } else {
        z <- as.matrix(z)
    }

    if (length(N) != n || length(group) != n || nrow(z) != n ||
        any(!is.finite(y)) || any(!is.finite(N)) ||
        any(y < 0) || any(N <= 0) || any(y > N)) {
        return(.dasra_abundance_empty_fit("invalid_input", n))
    }

    design <- .dasra_abundance_designs(group, z)
    positive <- y > 0
    if (sum(positive) < 3L) {
        return(.dasra_abundance_empty_fit(
            "fewer_than_three_positive_counts", n
        ))
    }
    if (length(unique(group[positive])) < 2L) {
        return(.dasra_abundance_empty_fit(
            "positive_counts_in_one_group_only", n
        ))
    }
    if (qr(design$X_full)$rank < ncol(design$X_full) ||
        qr(design$X_rho)$rank < ncol(design$X_rho) ||
        qr(design$X_full[positive, , drop = FALSE])$rank <
            ncol(design$X_full)) {
        return(.dasra_abundance_empty_fit("rank_deficient_design", n))
    }

    layout <- .dasra_abundance_layout(design$X_b, design$X_rho)
    bounds <- .dasra_abundance_bounds(layout)
    data <- list(
        y = y,
        N = N,
        group = group,
        X_b = design$X_b,
        X_rho = design$X_rho,
        X_full = design$X_full,
        layout = layout
    )

    mark_fit <- .dasra_abundance_initial_mark(data, gh_fit)
    if (is.null(mark_fit)) {
        return(.dasra_abundance_empty_fit(
            "positive_mark_initialization_failed", n
        ))
    }

    full_coefficient <- mark_fit$par[seq_len(ncol(design$X_full))]
    b_start <- c(full_coefficient[1L], full_coefficient[-c(1L, 2L)])
    if (ncol(design$X_full) == 2L) b_start <- full_coefficient[1L]
    zeta_start <- full_coefficient[2L]
    omega_start <- mark_fit$par[ncol(design$X_full) + 1L]

    detection_fit <- .dasra_abundance_initial_detection(
        data = data, b = b_start, omega = omega_start,
        zeta = zeta_start, gh = gh_fit
    )
    if (is.null(detection_fit)) {
        return(.dasra_abundance_empty_fit(
            "detection_initialization_failed", n
        ))
    }

    theta_start <- numeric(layout$dimension)
    theta_start[layout$b] <- b_start
    theta_start[layout$omega] <- omega_start
    theta_start[layout$a] <- detection_fit$par
    theta_start[layout$zeta] <- zeta_start
    theta_start <- count_clamp(theta_start, bounds$lower, bounds$upper)

    initial_state <- tryCatch(
        .dasra_abundance_state(
            theta_start, data, gh_fit, return_psi = TRUE
        ),
        error = function(e) NULL
    )
    if (is.null(initial_state)) {
        return(.dasra_abundance_empty_fit(
            "initial_estimating_equations_failed", n
        ))
    }
    equation_scale <- sqrt(colMeans(initial_state$psi^2))
    equation_scale[
        !is.finite(equation_scale) | equation_scale < 1e-3
    ] <- 1

    equation_mean <- function(theta) {
        state <- tryCatch(
            .dasra_abundance_state(
                theta, data, gh_fit, return_psi = TRUE
            ),
            error = function(e) NULL
        )
        if (is.null(state)) return(rep(1e3, layout$dimension))
        colMeans(state$psi)
    }
    equation_scaled <- function(theta) equation_mean(theta) / equation_scale
    candidates <- list(theta_start)

    root_1 <- tryCatch(
        nleqslv::nleqslv(
            theta_start, equation_scaled,
            method = "Broyden", global = "dbldog",
            control = list(
                ftol = 1e-9, xtol = 1e-9, maxit = 150L,
                allowSingular = FALSE
            )
        ),
        error = function(e) NULL
    )
    if (!is.null(root_1) && all(is.finite(root_1$x))) {
        candidates[[length(candidates) + 1L]] <- root_1$x
    }

    best_before_optim <- candidates[[which.min(vapply(
        candidates,
        function(x) {
            value <- equation_scaled(x)
            if (any(!is.finite(value))) Inf else sum(value^2)
        },
        numeric(1)
    ))]]
    best_before_optim <- count_clamp(
        best_before_optim, bounds$lower, bounds$upper
    )

    score_objective <- function(theta) {
        value <- equation_scaled(theta)
        if (any(!is.finite(value))) return(1e12)
        sum(value^2) / 2
    }
    squared_fit <- tryCatch(
        optim(
            best_before_optim, score_objective, method = "L-BFGS-B",
            lower = bounds$lower, upper = bounds$upper,
            control = list(maxit = 500L, factr = 1e7)
        ),
        error = function(e) NULL
    )
    if (!is.null(squared_fit) && is.finite(squared_fit$value) &&
        all(is.finite(squared_fit$par))) {
        candidates[[length(candidates) + 1L]] <- squared_fit$par
        root_2 <- tryCatch(
            nleqslv::nleqslv(
                squared_fit$par, equation_scaled,
                method = "Newton", global = "dbldog",
                control = list(
                    ftol = 1e-10, xtol = 1e-10, maxit = 100L,
                    allowSingular = FALSE
                )
            ),
            error = function(e) NULL
        )
        if (!is.null(root_2) && all(is.finite(root_2$x))) {
            candidates[[length(candidates) + 1L]] <- root_2$x
        }
    }

    candidate_metric <- vapply(candidates, function(theta) {
        if (any(theta < bounds$lower) || any(theta > bounds$upper)) {
            return(Inf)
        }
        value <- equation_mean(theta)
        if (any(!is.finite(value))) Inf else max(abs(value))
    }, numeric(1))
    if (!any(is.finite(candidate_metric))) {
        return(.dasra_abundance_empty_fit(
            "estimating_equation_solver_failed", n
        ))
    }

    theta_hat <- candidates[[which.min(candidate_metric)]]
    residue <- min(candidate_metric)
    if (!is.finite(residue) || residue > control$root_tolerance) {
        return(.dasra_abundance_empty_fit(
            "estimating_equation_residue", n
        ))
    }

    boundary_distance <- pmin(
        theta_hat - bounds$lower, bounds$upper - theta_hat
    )
    if (any(boundary_distance <= control$parameter_boundary_tolerance)) {
        return(.dasra_abundance_empty_fit("parameter_on_boundary", n))
    }

    fitted_state <- tryCatch(
        .dasra_abundance_state(
            theta_hat, data, gh_fit, return_psi = TRUE
        ),
        error = function(e) NULL
    )
    if (is.null(fitted_state) || any(!is.finite(fitted_state$psi))) {
        return(.dasra_abundance_empty_fit("fitted_state_failed", n))
    }

    jacobian_mean <- .dasra_abundance_numeric_jacobian(
        theta_hat, equation_mean,
        lower = bounds$lower, upper = bounds$upper
    )
    condition <- .dasra_abundance_condition_number(jacobian_mean)
    if (any(!is.finite(jacobian_mean)) || !is.finite(condition) ||
        condition > control$jacobian_condition_limit) {
        return(.dasra_abundance_empty_fit(
            "singular_or_ill_conditioned_jacobian", n
        ))
    }

    delta_hat <- .dasra_abundance_effect(theta_hat, data, gh_effect)
    delta_gradient <- .dasra_abundance_numeric_gradient(
        theta_hat,
        function(theta) .dasra_abundance_effect(
            theta, data, gh_effect
        ),
        lower = bounds$lower, upper = bounds$upper
    )
    if (!is.finite(delta_hat) || any(!is.finite(delta_gradient))) {
        return(.dasra_abundance_empty_fit("effect_gradient_failed", n))
    }

    A_sum <- n * jacobian_mean
    parameter_contribution <- tryCatch(
        solve(A_sum, t(fitted_state$psi)),
        error = function(e) NULL
    )
    if (is.null(parameter_contribution) ||
        any(!is.finite(parameter_contribution))) {
        return(.dasra_abundance_empty_fit("influence_solve_failed", n))
    }

    phi <- -as.numeric(delta_gradient %*% parameter_contribution)
    raw_variance <- sum(phi^2)
    if (!is.finite(raw_variance) || raw_variance <= 0) {
        return(.dasra_abundance_empty_fit(
            "nonpositive_raw_variance", n
        ))
    }

    raw_se <- sqrt(raw_variance)
    raw_z <- delta_hat / raw_se
    raw_p <- 2 * pnorm(-abs(raw_z))
    list(
        available = TRUE,
        status = "ok",
        theta = theta_hat,
        raw_delta = delta_hat,
        raw_se = raw_se,
        raw_p = raw_p,
        phi = phi,
        score_residue = residue,
        jacobian_condition = condition,
        mean_presence_weight = mean(fitted_state$presence_weight),
        mean_zero_presence_weight = if (any(y == 0)) {
            mean(fitted_state$presence_weight[y == 0])
        } else {
            NA_real_
        }
    )
}

.dasra_abundance_lts_pilot <- function(values) {
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    m <- length(values)
    if (m < 3L) return(NA_real_)
    h <- floor(m / 2L) + 1L
    ordered <- sort(values)
    cumulative <- c(0, cumsum(ordered))
    cumulative_square <- c(0, cumsum(ordered^2))
    starts <- seq_len(m - h + 1L)
    sums <- cumulative[starts + h] - cumulative[starts]
    sums_square <- cumulative_square[starts + h] -
        cumulative_square[starts]
    sse <- sums_square - sums^2 / h
    best <- starts[which.min(sse)]
    mean(ordered[best:(best + h - 1L)])
}

.dasra_abundance_correct <- function(
        fits, taxa, n_samples, keep_diagnostics) {
    p_taxa <- length(fits)
    raw_delta <- vapply(fits, function(x) x$raw_delta, numeric(1))
    raw_se <- vapply(fits, function(x) x$raw_se, numeric(1))
    raw_p <- vapply(fits, function(x) x$raw_p, numeric(1))
    available <- vapply(
        fits, function(x) isTRUE(x$available), logical(1)
    )
    fit_status <- vapply(
        fits, function(x) as.character(x$status)[1L], character(1)
    )

    corrected_p <- rep(1, p_taxa)
    corrected_estimate <- corrected_se <- corrected_z <-
        rep(NA_real_, p_taxa)
    formed <- rep(FALSE, p_taxa)
    reason <- fit_status
    background_size <- rep(NA_integer_, p_taxa)
    background_pilot <- background_estimate <- rep(NA_real_, p_taxa)
    reference_taxa <- vector("list", p_taxa)

    eligible <- which(
        available & is.finite(raw_delta) & is.finite(raw_se) & raw_se > 0
    )
    if (length(eligible) < 5L) {
        reason[eligible] <- "fewer_than_five_eligible_taxa"
    } else {
        phi_matrix <- matrix(NA_real_, nrow = n_samples, ncol = p_taxa)
        for (j in eligible) phi_matrix[, j] <- fits[[j]]$phi
        threshold <- sqrt(2 * log(max(n_samples * length(eligible), 3)))

        for (j in eligible) {
            others <- setdiff(eligible, j)
            if (length(others) < 3L) {
                reason[j] <- "fewer_than_three_candidate_reference_taxa"
                next
            }
            pilot <- .dasra_abundance_lts_pilot(raw_delta[others])
            if (!is.finite(pilot)) {
                reason[j] <- "background_pilot_unavailable"
                next
            }
            studentized <- abs(raw_delta[others] - pilot) /
                pmax(raw_se[others], 1e-10)
            reference <- others[studentized <= threshold]
            if (length(reference) < 3L) {
                reason[j] <- "fewer_than_three_reference_taxa"
                next
            }

            background <- mean(raw_delta[reference])
            phi_background <- rowMeans(
                phi_matrix[, reference, drop = FALSE]
            )
            corrected_phi <- phi_matrix[, j] - phi_background
            variance <- sum(corrected_phi^2)
            if (!is.finite(variance) || variance <= 0) {
                reason[j] <- "nonpositive_corrected_variance"
                next
            }

            corrected <- raw_delta[j] - background
            se <- sqrt(variance)
            z_value <- corrected / se
            p_value <- 2 * pnorm(-abs(z_value))
            if (!is.finite(p_value) || p_value < 0 || p_value > 1) {
                reason[j] <- "nonfinite_p_value"
                next
            }

            corrected_p[j] <- p_value
            corrected_estimate[j] <- corrected
            corrected_se[j] <- se
            corrected_z[j] <- z_value
            formed[j] <- TRUE
            reason[j] <- "ok"
            background_size[j] <- length(reference)
            background_pilot[j] <- pilot
            background_estimate[j] <- background
            reference_taxa[[j]] <- taxa[reference]
        }
    }

    diagnostics <- NULL
    if (keep_diagnostics) {
        diagnostics <- data.frame(
            taxon = taxa,
            raw_estimate = raw_delta,
            raw_standard_error = raw_se,
            raw_p_value = raw_p,
            corrected_estimate = corrected_estimate,
            corrected_standard_error = corrected_se,
            corrected_p_value = corrected_p,
            formed = formed,
            reason = reason,
            background_size = background_size,
            background_pilot = background_pilot,
            background_estimate = background_estimate,
            score_residue = vapply(
                fits, function(x) x$score_residue, numeric(1)
            ),
            jacobian_condition = vapply(
                fits, function(x) x$jacobian_condition, numeric(1)
            ),
            mean_presence_weight = vapply(
                fits, function(x) x$mean_presence_weight, numeric(1)
            ),
            mean_zero_presence_weight = vapply(
                fits, function(x) x$mean_zero_presence_weight, numeric(1)
            ),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        diagnostics$reference_taxa <- I(reference_taxa)
        rownames(diagnostics) <- taxa
        diagnostics <- list(taxon = diagnostics, raw_fits = fits)
    }

    list(
        p = corrected_p,
        formed = formed,
        reason = reason,
        estimate = corrected_estimate,
        se = corrected_se,
        z = corrected_z,
        diagnostics = diagnostics
    )
}

.dasra_abundance_arm <- function(Y, N, g, z, keep_diagnostics) {
    taxa <- colnames(Y)
    n_samples <- nrow(Y)
    control <- .dasra_abundance_control()
    gh_fit <- make_count_gh_rule(control$quadrature_Q)
    gh_effect <- make_count_gh_rule(control$effect_quadrature_Q)
    fits <- vector("list", ncol(Y))
    for (j in seq_len(ncol(Y))) {
        fits[[j]] <- tryCatch(
            .dasra_abundance_fit_taxon(
                y = as.numeric(Y[, j]),
                N = N,
                group = g,
                z = z,
                gh_fit = gh_fit,
                gh_effect = gh_effect,
                control = control
            ),
            error = function(e) {
                .dasra_abundance_empty_fit(
                    paste0(
                        "relative_abundance_fit_error: ",
                        conditionMessage(e)
                    ),
                    n_samples
                )
            }
        )
    }
    names(fits) <- taxa
    .dasra_abundance_correct(
        fits = fits,
        taxa = taxa,
        n_samples = n_samples,
        keep_diagnostics = keep_diagnostics
    )
}

# --------------------------------------------------------------------------
# Structural-absence score engine
# --------------------------------------------------------------------------

zt_log1mexp <- function(log_x) {
    log_x <- as.numeric(log_x)
    log_x[log_x > 0 & log_x < 1e-10] <- 0
    answer <- rep(NaN, length(log_x))
    valid <- is.finite(log_x) & log_x <= 0
    low <- valid & log_x <= log(0.5)
    high <- valid & !low
    answer[low] <- log1p(-exp(log_x[low]))
    answer[high] <- log(-expm1(log_x[high]))
    answer
}

zt_clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

zt_central_derivative_matrix <- function(fn, par, lower = NULL,
                                         upper = NULL, rel_step = 1e-04) {
    par <- as.numeric(par)
    f0 <- as.numeric(fn(par))
    answer <- matrix(NA, nrow = length(f0), ncol = length(par))
    for (k in seq_along(par)) {
        h <- rel_step * max(1, abs(par[k]))
        if (!is.null(lower)) h <- min(h, 0.45 * (par[k] - lower[k]))
        if (!is.null(upper)) h <- min(h, 0.45 * (upper[k] - par[k]))
        if (!is.finite(h) || h <= .Machine$double.eps^0.4) next
        plus <- minus <- par
        plus[k] <- par[k] + h
        minus[k] <- par[k] - h
        answer[, k] <- (as.numeric(fn(plus)) - as.numeric(fn(minus))) /
            (2 * h)
    }
    answer
}

zt_numeric_gradient <- function(fn, par, lower = NULL, upper = NULL,
                                rel_step = 1e-05) {
    as.numeric(zt_central_derivative_matrix(function(x) fn(x), par, lower,
                                            upper, rel_step))
}

zt_inference_steps <- function(par, n, lower, upper, base = 1e-04,
                               reference_n = 120L) {
    par <- as.numeric(par)
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    n <- as.integer(n)
    reference_n <- as.integer(reference_n)
    if (length(par) != length(lower) || length(par) != length(upper) ||
        length(n) != 1L || !is.finite(n) || n < 1L ||
        length(reference_n) != 1L || !is.finite(reference_n) ||
        reference_n < 1L || length(base) != 1L || !is.finite(base) ||
        base <= 0 || any(!is.finite(par)) || any(par <= lower) ||
        any(par >= upper)) {
        return(NULL)
    }
    rate <- base * min(1, (reference_n / n)^(1 / 3))
    step <- rate * pmax(1, abs(par))
    step <- pmin(step, 0.45 * (par - lower), 0.45 * (upper - par))
    if (any(!is.finite(step)) ||
        any(step <= .Machine$double.eps^0.4)) {
        return(NULL)
    }
    list(base = base, reference_n = reference_n, rate = rate, step = step)
}

zt_central_derivative_matrix_fixed <- function(fn, par, step, f0 = NULL) {
    par <- as.numeric(par)
    step <- as.numeric(step)
    if (length(step) != length(par) || any(!is.finite(step)) ||
        any(step <= 0)) {
        stop("Invalid fixed finite-difference steps.")
    }
    if (is.null(f0)) f0 <- fn(par)
    f0 <- as.numeric(f0)
    answer <- matrix(NA_real_, nrow = length(f0), ncol = length(par))
    for (k in seq_along(par)) {
        plus <- minus <- par
        plus[k] <- par[k] + step[k]
        minus[k] <- par[k] - step[k]
        f_plus <- as.numeric(fn(plus))
        f_minus <- as.numeric(fn(minus))
        if (length(f_plus) != length(f0) ||
            length(f_minus) != length(f0)) {
            stop("A finite-difference function changed output length.")
        }
        answer[, k] <- (f_plus - f_minus) / (2 * step[k])
    }
    answer
}

zt_central_hessian_fixed <- function(fn, par, step, f0 = NULL) {
    par <- as.numeric(par)
    step <- as.numeric(step)
    p <- length(par)
    if (length(step) != p || any(!is.finite(step)) || any(step <= 0)) {
        stop("Invalid fixed finite-difference steps.")
    }
    if (is.null(f0)) f0 <- fn(par)
    f0 <- as.numeric(f0)
    if (length(f0) != 1L || !is.finite(f0)) {
        stop("The Hessian criterion must be finite and scalar.")
    }
    answer <- matrix(NA_real_, p, p)
    for (k in seq_len(p)) {
        plus <- minus <- par
        plus[k] <- par[k] + 2 * step[k]
        minus[k] <- par[k] - 2 * step[k]
        answer[k, k] <-
            (fn(plus) - 2 * f0 + fn(minus)) / (4 * step[k]^2)
    }
    if (p > 1L) {
        for (k in seq_len(p - 1L)) {
            for (l in (k + 1L):p) {
                pp <- pm <- mp <- mm <- par
                pp[c(k, l)] <- par[c(k, l)] + step[c(k, l)]
                pm[k] <- par[k] + step[k]
                pm[l] <- par[l] - step[l]
                mp[k] <- par[k] - step[k]
                mp[l] <- par[l] + step[l]
                mm[c(k, l)] <- par[c(k, l)] - step[c(k, l)]
                value <- (fn(pp) - fn(pm) - fn(mp) + fn(mm)) /
                    (4 * step[k] * step[l])
                answer[k, l] <- answer[l, k] <- value
            }
        }
    }
    (answer + t(answer)) / 2
}

zt_beta_components <- function(beta, y, N, X_eta, gh) {
    p_eta <- ncol(X_eta)
    if (length(beta) != p_eta + 1L || any(!is.finite(beta))) return(NULL)
    sigma <- exp(beta[p_eta + 1L])
    if (!is.finite(sigma) || sigma <= 0) return(NULL)
    eta <- as.numeric(X_eta %*% beta[seq_len(p_eta)])
    log_h0 <- count_log_hy_adaptive(
        y = rep(0, length(y)), N = N, eta = eta, sigma = sigma, gh = gh
    )
    log_hy <- log_h0
    pos <- y > 0
    if (any(pos)) {
        log_hy[pos] <- count_log_hy_adaptive(
            y = y[pos], N = N[pos], eta = eta[pos], sigma = sigma, gh = gh
        )
    }
    if (any(!is.finite(log_hy)) || any(!is.finite(log_h0)) ||
        any(log_h0 >= -1e-14)) {
        return(NULL)
    }
    log_r <- zt_log1mexp(log_h0)
    if (any(!is.finite(log_r))) return(NULL)
    list(eta = eta, sigma = sigma, log_hy = log_hy, log_h0 = log_h0,
         log_r = log_r, r = exp(log_r))
}

zt_beta_detection_components <- function(beta, N, X_eta, gh) {
    p_eta <- ncol(X_eta)
    if (length(beta) != p_eta + 1L || any(!is.finite(beta))) return(NULL)
    sigma <- exp(beta[p_eta + 1L])
    if (!is.finite(sigma) || sigma <= 0) return(NULL)
    eta <- as.numeric(X_eta %*% beta[seq_len(p_eta)])
    log_h0 <- count_log_hy_adaptive(
        y = rep(0, length(N)), N = N, eta = eta, sigma = sigma, gh = gh
    )
    if (any(!is.finite(log_h0)) || any(log_h0 >= -1e-14)) return(NULL)
    log_r <- zt_log1mexp(log_h0)
    if (any(!is.finite(log_r))) return(NULL)
    list(eta = eta, sigma = sigma, log_h0 = log_h0,
         log_r = log_r, r = exp(log_r))
}

zt_beta_loglik_by_sample <- function(beta, y, N, X_eta, gh) {
    out <- numeric(length(y))
    comp <- zt_beta_components(beta, y, N, X_eta, gh)
    if (is.null(comp)) return(rep(-1e+10, length(y)))
    pos <- y > 0
    out[pos] <- comp$log_hy[pos] - comp$log_r[pos]
    out
}

zt_beta_loglik_by_sample_inference <- function(beta, y, N, X_eta, gh) {
    out <- numeric(length(y))
    comp <- zt_beta_components(beta, y, N, X_eta, gh)
    if (is.null(comp)) return(rep(NA_real_, length(y)))
    pos <- y > 0
    out[pos] <- comp$log_hy[pos] - comp$log_r[pos]
    out
}

zt_beta_nll <- function(beta, y, N, X_eta, gh) {
    value <- -sum(zt_beta_loglik_by_sample(beta, y, N, X_eta, gh))
    if (is.finite(value)) value else 1e+12
}

zt_select_optim <- function(starts, fn, lower, upper, maxit = 500L,
                            gr = NULL) {
    fits <- lapply(starts, function(start) {
        control <- list(maxit = maxit, factr = 1e+07, pgtol = 1e-06)
        if (is.null(gr)) {
            control$ndeps <- rep(1e-05, length(start))
        }
        tryCatch(
            optim(
                par = zt_clamp(start, lower, upper),
                fn = fn,
                gr = gr,
                method = "L-BFGS-B",
                lower = lower,
                upper = upper,
                control = control
            ),
            error = function(e) NULL
        )
    })
    finite <- vapply(fits, function(x) {
        !is.null(x) && is.finite(x$value) && all(is.finite(x$par))
    }, logical(1))
    converged <- finite & vapply(fits, function(x) {
        !is.null(x) && identical(x$convergence, 0L)
    }, logical(1))
    pool <- if (any(converged)) which(converged) else which(finite)
    if (!length(pool)) return(NULL)
    values <- vapply(pool, function(k) fits[[k]]$value, numeric(1))
    best <- fits[[pool[which.min(values)]]]
    best$n_starts <- length(starts)
    best$n_finite <- sum(finite)
    best$n_converged <- sum(converged)
    best
}

zt_fit_beta <- function(y, N, X_eta, gh, maxit = 500L) {
    p_eta <- ncol(X_eta)
    pos <- y > 0
    n_positive <- sum(pos)
    pseudo_prop <- zt_clamp((y[pos] + 0.5) / (N[pos] + 1), 1e-10,
                            1 - 1e-10)
    regression <- tryCatch(
        lm.fit(X_eta[pos, , drop = FALSE], qlogis(pseudo_prop))$coefficients,
        error = function(e) rep(NA_real_, p_eta)
    )
    if (length(regression) != p_eta) regression <- rep(NA_real_, p_eta)
    if (!is.finite(regression[1L])) {
        regression[1L] <- mean(qlogis(pseudo_prop))
    }
    bad_slope <- !is.finite(regression)
    bad_slope[1L] <- FALSE
    regression[bad_slope] <- 0
    intercept_only <- c(mean(qlogis(pseudo_prop)), rep(0, p_eta - 1L))
    lower <- c(-30, rep(-10, p_eta - 1L), log(0.1))
    upper <- c(5, rep(10, p_eta - 1L), log(8))
    starts <- list(
        c(regression, log(0.5)),
        c(regression, log(1)),
        c(regression, log(2)),
        c(intercept_only, log(0.5)),
        c(intercept_only, log(2))
    )
    starts <- lapply(starts, zt_clamp, lower = lower, upper = upper)
    keys <- vapply(starts, function(x) paste(signif(x, 12), collapse = "|"),
                   character(1))
    starts <- starts[!duplicated(keys)]

    fn_unscaled <- function(beta) zt_beta_nll(beta, y, N, X_eta, gh)
    fn_optimization <- function(beta) fn_unscaled(beta) / n_positive
    fit <- zt_select_optim(starts, fn_optimization, lower, upper, maxit)
    if (is.null(fit)) {
        return(list(
            ok = FALSE,
            reason = "conditional_present_no_finite_fit"
        ))
    }

    gradient <- zt_numeric_gradient(fn_optimization, fit$par, lower, upper)
    hessian <- tryCatch(
        optimHess(fit$par, fn_optimization),
        error = function(e) NULL
    )
    eigenvalues <- if (!is.null(hessian) && all(is.finite(hessian))) {
        tryCatch(
            eigen((hessian + t(hessian)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values,
            error = function(e) NULL
        )
    } else {
        NULL
    }
    bound_distance <- pmin(fit$par - lower, upper - fit$par)
    gradient_ok <- all(is.finite(gradient)) && max(abs(gradient)) <= 0.05
    information_ok <- length(eigenvalues) == length(fit$par) &&
        min(eigenvalues) > 0
    condition_ok <- information_ok &&
        max(eigenvalues) / min(eigenvalues) <= 1e+10
    ok <- identical(fit$convergence, 0L) && gradient_ok &&
        all(bound_distance > 1e-05) && information_ok && condition_ok
    reason <- if (!identical(fit$convergence, 0L)) {
        "conditional_present_nonconvergence"
    } else if (!gradient_ok) {
        "conditional_present_gradient"
    } else if (any(bound_distance <= 1e-05)) {
        "conditional_present_boundary"
    } else if (!information_ok) {
        "conditional_present_information_nonpositive"
    } else if (!condition_ok) {
        "conditional_present_information_ill_conditioned"
    } else {
        "ok"
    }
    list(
        ok = ok,
        reason = reason,
        par = fit$par,
        nll = fn_unscaled,
        optimization_objective = fn_optimization,
        lower = lower,
        upper = upper,
        gradient = gradient,
        hessian = hessian,
        eigenvalues = eigenvalues,
        optimization = fit
    )
}

zt_stable_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

zt_detection_state_from_components <- function(alpha, y, X_rho, comp) {
    if (is.null(comp)) return(NULL)
    lp_rho <- as.numeric(X_rho %*% alpha)
    rho <- plogis(lp_rho)
    log_rho <- -zt_stable_softplus(-lp_rho)
    log_q <- comp$log_r - zt_stable_softplus(lp_rho)
    if (any(!is.finite(log_q)) || any(log_q >= 0)) return(NULL)
    log_zero_prob <- zt_log1mexp(log_q)
    if (any(!is.finite(log_zero_prob))) return(NULL)
    detected <- y > 0
    structural_residual <- numeric(length(y))
    structural_residual[detected] <- -rho[detected]
    structural_residual[!detected] <- exp(
        log_rho[!detected] + log_q[!detected] - log_zero_prob[!detected]
    )
    if (any(!is.finite(structural_residual))) return(NULL)
    gamma <- structural_residual + rho
    list(rho = rho, q = exp(log_q), gamma = gamma,
         structural_residual = structural_residual, component = comp,
         log_q = log_q, log_zero_prob = log_zero_prob)
}

zt_detection_state <- function(alpha, beta, y, N, X_rho, X_eta, gh) {
    comp <- zt_beta_components(beta, y, N, X_eta, gh)
    zt_detection_state_from_components(alpha, y, X_rho, comp)
}

zt_detection_nll_from_components <- function(alpha, y, X_rho, comp) {
    state <- zt_detection_state_from_components(alpha, y, X_rho, comp)
    if (is.null(state)) return(1e+300)
    detected <- y > 0
    value <- -sum(ifelse(detected, state$log_q, state$log_zero_prob))
    if (is.finite(value)) value else 1e+300
}

zt_detection_nll <- function(alpha, beta, y, N, X_rho, X_eta, gh) {
    comp <- zt_beta_detection_components(beta, N, X_eta, gh)
    zt_detection_nll_from_components(alpha, y, X_rho, comp)
}

zt_alpha_scores_from_components <- function(alpha, y, X_rho, comp) {
    state <- zt_detection_state_from_components(alpha, y, X_rho, comp)
    if (is.null(state)) return(matrix(NA, length(y), ncol(X_rho)))
    X_rho * state$structural_residual
}

zt_alpha_scores <- function(alpha, beta, y, N, X_rho, X_eta, gh) {
    comp <- zt_beta_detection_components(beta, N, X_eta, gh)
    zt_alpha_scores_from_components(alpha, y, X_rho, comp)
}

zt_target_score_from_components <- function(alpha, y, g, X_rho, comp) {
    state <- zt_detection_state_from_components(alpha, y, X_rho, comp)
    if (is.null(state)) return(rep(NA, length(y)))
    as.numeric(g) * state$structural_residual
}

zt_target_score <- function(alpha, beta, y, N, g, X_rho, X_eta, gh) {
    comp <- zt_beta_detection_components(beta, N, X_eta, gh)
    zt_target_score_from_components(alpha, y, g, X_rho, comp)
}

zt_fit_alpha <- function(beta, y, N, X_rho, X_eta, gh, maxit = 500L) {
    p_alpha <- ncol(X_rho)
    zero_fraction <- mean(y == 0)
    lower <- rep(-10, p_alpha)
    upper <- rep(10, p_alpha)
    structural_starts <- c(0.05, 0.25, 0.5, 0.75) *
        max(zero_fraction, 0.05)
    starts <- lapply(structural_starts, function(p) {
        c(qlogis(zt_clamp(p, 0.01, 0.9)), rep(0, p_alpha - 1L))
    })
    conditional_present <- zt_beta_detection_components(beta, N, X_eta, gh)
    if (is.null(conditional_present)) {
        return(list(
            ok = FALSE,
            reason = "structural_absence_no_finite_fit"
        ))
    }
    n_observations <- length(y)
    fn_unscaled <- function(alpha) {
        zt_detection_nll_from_components(
            alpha, y, X_rho, conditional_present
        )
    }
    fn_optimization <- function(alpha) fn_unscaled(alpha) / n_observations
    gr_optimization <- function(alpha) {
        -colSums(zt_alpha_scores_from_components(
            alpha, y, X_rho, conditional_present
        )) / n_observations
    }
    fit <- zt_select_optim(
        starts, fn_optimization, lower, upper, maxit,
        gr = gr_optimization
    )
    if (is.null(fit)) {
        return(list(
            ok = FALSE,
            reason = "structural_absence_no_finite_fit"
        ))
    }

    gradient <- zt_numeric_gradient(fn_optimization, fit$par, lower, upper)
    hessian <- tryCatch(
        optimHess(fit$par, fn_optimization),
        error = function(e) NULL
    )
    eigenvalues <- if (!is.null(hessian) && all(is.finite(hessian))) {
        tryCatch(
            eigen((hessian + t(hessian)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values,
            error = function(e) NULL
        )
    } else {
        NULL
    }
    bound_distance <- pmin(fit$par - lower, upper - fit$par)
    gradient_ok <- all(is.finite(gradient)) && max(abs(gradient)) <= 0.05
    information_ok <- length(eigenvalues) == length(fit$par) &&
        min(eigenvalues) > 0
    condition_ok <- information_ok &&
        max(eigenvalues) / min(eigenvalues) <= 1e+10
    ok <- identical(fit$convergence, 0L) && gradient_ok &&
        all(bound_distance > 1e-05) && information_ok && condition_ok
    reason <- if (!identical(fit$convergence, 0L)) {
        "structural_absence_nonconvergence"
    } else if (!gradient_ok) {
        "structural_absence_gradient"
    } else if (any(bound_distance <= 1e-05)) {
        "structural_absence_boundary"
    } else if (!information_ok) {
        "structural_absence_information_nonpositive"
    } else if (!condition_ok) {
        "structural_absence_information_ill_conditioned"
    } else {
        "ok"
    }
    list(
        ok = ok,
        reason = reason,
        par = fit$par,
        nll = fn_unscaled,
        optimization_objective = fn_optimization,
        lower = lower,
        upper = upper,
        gradient = gradient,
        hessian = hessian,
        eigenvalues = eigenvalues,
        optimization = fit
    )
}

zt_unpack_theta <- function(theta, p_beta) {
    list(beta = theta[seq_len(p_beta)], alpha = theta[-seq_len(p_beta)])
}

zt_unavailable <- function(reason, n, diagnostics = list()) {
    list(p = 1, tested = FALSE, reason = reason, gamma = rep(NA, n),
         rho = rep(NA, n), diagnostics = diagnostics)
}

zt_count_structural_test <- function(y, N, g, z = NULL, Q = 1001L,
                                     min_positive_samples = 3L,
                                     derivative_base = 1e-04,
                                     derivative_reference_n = 120L,
                                     keep_fit = FALSE) {
    y <- as.numeric(y)
    N <- as.numeric(N)
    g <- as.numeric(g)
    n <- length(y)
    if (length(N) != n || length(g) != n || any(!is.finite(c(y, N, g))) ||
        any(y < 0) || any(N <= 0) || any(y > N) ||
        any(abs(y - round(y)) > 1e-08) ||
        any(abs(N - round(N)) > 1e-08) ||
        !identical(sort(unique(g)), c(0, 1))) {
        stop("Invalid inputs to zt_count_structural_test.")
    }
    if (length(min_positive_samples) != 1L ||
        !is.finite(min_positive_samples) || min_positive_samples < 1L ||
        abs(min_positive_samples - round(min_positive_samples)) > 1e-08) {
        stop("min_positive_samples must be one positive integer.")
    }
    min_positive_samples <- as.integer(round(min_positive_samples))
    if (!is.null(z)) {
        z_check <- as.matrix(z)
        if (!is.numeric(z_check) || nrow(z_check) != n ||
            any(!is.finite(z_check))) {
            stop("z must be a finite numeric covariate matrix with one row per sample.")
        }
    }

    X_rho <- count_design_matrix(g = g, z = z, include_group = FALSE, n = n)
    X_eta <- count_design_matrix(g = g, z = z, include_group = TRUE, n = n)
    n_positive <- sum(y > 0)
    n_pos0 <- sum(y > 0 & g == 0)
    n_pos1 <- sum(y > 0 & g == 1)
    n_zero <- sum(y == 0)
    base_diag <- list(
        n_positive = n_positive,
        n_pos0 = n_pos0,
        n_pos1 = n_pos1,
        n_zero = n_zero
    )
    if (n_positive < min_positive_samples) {
        return(zt_unavailable("insufficient_positive_support", n, base_diag))
    }
    if (n_zero < 1L) {
        return(zt_unavailable("no_observed_zeros", n, base_diag))
    }
    if (qr(X_rho)$rank < ncol(X_rho) || qr(X_eta)$rank < ncol(X_eta)) {
        return(zt_unavailable("model_design_rank_deficient", n, base_diag))
    }
    if (qr(X_eta[y > 0, , drop = FALSE])$rank < ncol(X_eta)) {
        return(zt_unavailable("positive_part_design_rank_deficient", n,
                              base_diag))
    }

    gh <- make_structural_gh_rule(Q)
    beta_fit <- zt_fit_beta(y, N, X_eta, gh)
    if (!isTRUE(beta_fit$ok)) {
        details <- if (keep_fit) c(base_diag, list(beta_fit = beta_fit)) else base_diag
        return(zt_unavailable(beta_fit$reason, n, details))
    }
    beta <- beta_fit$par
    alpha_fit <- zt_fit_alpha(beta, y, N, X_rho, X_eta, gh)
    if (!isTRUE(alpha_fit$ok)) {
        details <- if (keep_fit) {
            c(base_diag, list(beta_fit = beta_fit, alpha_fit = alpha_fit))
        } else {
            base_diag
        }
        return(zt_unavailable(alpha_fit$reason, n, details))
    }
    alpha <- alpha_fit$par
    detection_component <- zt_beta_detection_components(
        beta, N, X_eta, gh
    )
    if (is.null(detection_component)) {
        return(zt_unavailable("detection_state_evaluation_failed", n,
                              base_diag))
    }

    p_beta <- length(beta)
    p_alpha <- length(alpha)
    p_theta <- p_beta + p_alpha
    theta <- c(beta, alpha)
    lower <- c(beta_fit$lower, alpha_fit$lower)
    upper <- c(beta_fit$upper, alpha_fit$upper)
    step_info <- zt_inference_steps(
        theta, n, lower, upper, base = derivative_base,
        reference_n = derivative_reference_n
    )
    if (is.null(step_info)) {
        return(zt_unavailable("inference_derivative_step_unavailable", n,
                              base_diag))
    }
    theta_step <- step_info$step
    beta_step <- theta_step[seq_len(p_beta)]
    alpha_step <- theta_step[p_beta + seq_len(p_alpha)]

    beta_loglik_base <- zt_beta_loglik_by_sample_inference(
        beta, y, N, X_eta, gh
    )
    beta_score <- zt_central_derivative_matrix_fixed(
        function(b) zt_beta_loglik_by_sample_inference(
            b, y, N, X_eta, gh
        ),
        beta, beta_step, f0 = beta_loglik_base
    )
    alpha_score <- zt_alpha_scores_from_components(
        alpha, y, X_rho, detection_component
    )
    if (any(!is.finite(beta_score)) || any(!is.finite(alpha_score))) {
        return(zt_unavailable("score_evaluation_failed", n, base_diag))
    }
    psi <- cbind(beta_score, alpha_score)

    A <- matrix(0, p_theta, p_theta)
    beta_loglik_sum <- function(b) {
        sum(zt_beta_loglik_by_sample_inference(b, y, N, X_eta, gh))
    }
    beta_jacobian <- tryCatch(
        zt_central_hessian_fixed(
            beta_loglik_sum, beta, beta_step,
            f0 = sum(beta_loglik_base)
        ),
        error = function(e) NULL
    )
    if (is.null(beta_jacobian) || any(!is.finite(beta_jacobian))) {
        return(zt_unavailable("conditional_present_score_jacobian_failed", n,
                              base_diag))
    }
    A[seq_len(p_beta), seq_len(p_beta)] <-
        (beta_jacobian + t(beta_jacobian)) / 2
    alpha_beta_sum_fn <- function(b) {
        perturbed_component <- zt_beta_detection_components(
            b, N, X_eta, gh
        )
        colSums(zt_alpha_scores_from_components(
            alpha, y, X_rho, perturbed_component
        ))
    }
    alpha_alpha_sum_fn <- function(a) {
        colSums(zt_alpha_scores_from_components(
            a, y, X_rho, detection_component
        ))
    }
    alpha_score_sum <- colSums(alpha_score)
    alpha_beta_jacobian <- zt_central_derivative_matrix_fixed(
        alpha_beta_sum_fn, beta, beta_step, f0 = alpha_score_sum
    )
    alpha_alpha_jacobian <- zt_central_derivative_matrix_fixed(
        alpha_alpha_sum_fn, alpha, alpha_step, f0 = alpha_score_sum
    )
    A[p_beta + seq_len(p_alpha), seq_len(p_beta)] <-
        alpha_beta_jacobian
    A[p_beta + seq_len(p_alpha), p_beta + seq_len(p_alpha)] <-
        alpha_alpha_jacobian
    if (any(!is.finite(A))) {
        return(zt_unavailable("stacked_jacobian_failed", n, base_diag))
    }
    singular_values <- tryCatch(
        svd(A, nu = 0, nv = 0)$d,
        error = function(e) NULL
    )
    if (is.null(singular_values) || min(singular_values) <= 0 ||
        max(singular_values) / min(singular_values) > 1e+10) {
        details <- if (keep_fit) {
            c(base_diag, list(jacobian_singular_values = singular_values))
        } else {
            base_diag
        }
        return(zt_unavailable("stacked_jacobian_singular", n, details))
    }

    target <- zt_target_score_from_components(
        alpha, y, g, X_rho, detection_component
    )
    target_beta_fn <- function(b) {
        perturbed_component <- zt_beta_detection_components(
            b, N, X_eta, gh
        )
        zt_target_score_from_components(
            alpha, y, g, X_rho, perturbed_component
        )
    }
    target_alpha_fn <- function(a) {
        zt_target_score_from_components(
            a, y, g, X_rho, detection_component
        )
    }
    M <- c(
        colSums(zt_central_derivative_matrix_fixed(
            target_beta_fn, beta, beta_step, f0 = target
        )),
        colSums(zt_central_derivative_matrix_fixed(
            target_alpha_fn, alpha, alpha_step, f0 = target
        ))
    )
    if (any(!is.finite(target)) || any(!is.finite(M))) {
        return(zt_unavailable("target_derivative_failed", n, base_diag))
    }
    adjustment <- tryCatch(
        as.numeric(M %*% solve(A)),
        error = function(e) NULL
    )
    if (is.null(adjustment) || any(!is.finite(adjustment))) {
        return(zt_unavailable("stacked_projection_failed", n, base_diag))
    }
    influence <- target - as.numeric(psi %*% adjustment)
    if (any(!is.finite(influence))) {
        return(zt_unavailable("nonfinite_adjusted_score", n, base_diag))
    }

    U <- sum(target)
    V <- sum(influence^2)
    if (!is.finite(U) || !is.finite(V) || V <= 0) {
        return(zt_unavailable("nonpositive_robust_variance", n, base_diag))
    }
    statistic <- U^2 / V
    p <- pchisq(statistic, 1, lower.tail = FALSE)
    if (!is.finite(p)) {
        return(zt_unavailable("nonfinite_p_value", n, base_diag))
    }

    state <- zt_detection_state(alpha, beta, y, N, X_rho, X_eta, gh)
    if (is.null(state)) {
        return(zt_unavailable("detection_state_evaluation_failed", n,
                              base_diag))
    }

    fitted_sigma <- exp(beta[length(beta)])
    diagnostics <- c(base_diag, list(
        U = U,
        V = V,
        statistic = statistic,
        jacobian_condition = max(singular_values) / min(singular_values)
    ))

    if (keep_fit && is.finite(fitted_sigma) && fitted_sigma > 2) {
        comparison_Q <- max(2001L, as.integer(Q))
        quadrature_diagnostic <- tryCatch(
            {
                high_gh <- make_structural_gh_rule(comparison_Q)
                base_conditional <- zt_beta_loglik_by_sample(
                    beta, y, N, X_eta, gh
                )
                high_conditional <- zt_beta_loglik_by_sample(
                    beta, y, N, X_eta, high_gh
                )
                high_state <- zt_detection_state(
                    alpha, beta, y, N, X_rho, X_eta, high_gh
                )
                if (is.null(high_state)) {
                    stop("The higher-order detection state could not be evaluated.")
                }
                differences <- c(
                    conditional = max(abs(base_conditional - high_conditional)),
                    detection = max(
                        abs(state$component$r - high_state$component$r)
                    ),
                    posterior = max(abs(state$gamma - high_state$gamma))
                )
                if (any(!is.finite(differences))) {
                    stop("The higher-order quadrature differences are non-finite.")
                }
                list(ok = TRUE, differences = differences,
                     error = NA_character_)
            },
            error = function(e) {
                list(ok = FALSE, differences = rep(NA_real_, 3L),
                     error = conditionMessage(e))
            }
        )
        diagnostics$quadrature_checked <- TRUE
        diagnostics$quadrature_check_succeeded <- quadrature_diagnostic$ok
        diagnostics$quadrature_check_error <- quadrature_diagnostic$error
        diagnostics$quadrature_comparison_Q <- comparison_Q
        diagnostics$quadrature_conditional_max_abs <-
            unname(quadrature_diagnostic$differences[1L])
        diagnostics$quadrature_detection_max_abs <-
            unname(quadrature_diagnostic$differences[2L])
        diagnostics$quadrature_posterior_max_abs <-
            unname(quadrature_diagnostic$differences[3L])
    } else {
        diagnostics$quadrature_checked <- FALSE
        diagnostics$quadrature_check_succeeded <- NA
        diagnostics$quadrature_check_error <- NA_character_
        diagnostics$quadrature_comparison_Q <- NA_integer_
        diagnostics$quadrature_conditional_max_abs <- NA_real_
        diagnostics$quadrature_detection_max_abs <- NA_real_
        diagnostics$quadrature_posterior_max_abs <- NA_real_
    }

    if (keep_fit) {
        diagnostics$fit <- list(
            conditional_present = beta_fit,
            structural_absence = alpha_fit,
            stacked_jacobian = A,
            target_derivative = M,
            nuisance_adjustment = adjustment,
            estimating_functions = psi,
            adjusted_score = influence,
            derivative_step = theta_step,
            X_rho = X_rho,
            X_eta = X_eta,
            quadrature = gh
        )
    }
    list(
        p = p,
        tested = TRUE,
        reason = "ok",
        gamma = state$gamma,
        rho = state$rho,
        diagnostics = diagnostics
    )
}

# --------------------------------------------------------------------------
# Public interface
# --------------------------------------------------------------------------

.dasra_encode_group <- function(group, reference, n) {
    if (length(group) != n || anyNA(group)) {
        stop("`group` must have one non-missing value per sample.",
             call. = FALSE)
    }
    if (is.factor(group)) {
        value <- droplevels(group)
        level_values <- levels(value)
        if (length(level_values) != 2L) {
            stop("`group` must have exactly two observed levels.",
                 call. = FALSE)
        }
        key <- as.character(value)
        if (is.null(reference)) reference_key <- level_values[1L]
    } else if (is.logical(group)) {
        value <- as.logical(group)
        if (!identical(sort(unique(value)), c(FALSE, TRUE))) {
            stop("`group` must contain both groups.", call. = FALSE)
        }
        key <- as.character(value)
        level_values <- c("FALSE", "TRUE")
        if (is.null(reference)) reference_key <- "FALSE"
    } else if (is.numeric(group)) {
        value <- as.numeric(group)
        if (any(!is.finite(value))) {
            stop("`group` must be finite.", call. = FALSE)
        }
        observed <- sort(unique(value))
        if (length(observed) != 2L) {
            stop("`group` must have exactly two observed values.",
                 call. = FALSE)
        }
        key <- as.character(value)
        level_values <- as.character(observed)
        if (is.null(reference)) {
            if (!identical(observed, c(0, 1))) {
                stop("Supply `reference` when numeric `group` is not coded 0/1.",
                     call. = FALSE)
            }
            reference_key <- "0"
        }
    } else if (is.character(group)) {
        key <- as.character(group)
        level_values <- unique(key)
        if (length(level_values) != 2L) {
            stop("`group` must have exactly two observed values.",
                 call. = FALSE)
        }
        if (is.null(reference)) {
            stop("Supply `reference` for a character `group`.",
                 call. = FALSE)
        }
    } else {
        stop("`group` must be a factor, character, logical, or numeric vector.",
             call. = FALSE)
    }
    if (!is.null(reference)) {
        if (length(reference) != 1L || is.na(reference)) {
            stop("`reference` must be one non-missing group value.",
                 call. = FALSE)
        }
        reference_key <- as.character(reference)
    }
    if (!(reference_key %in% level_values)) {
        stop("`reference` is not an observed group value.", call. = FALSE)
    }
    comparison_key <- setdiff(level_values, reference_key)
    if (length(comparison_key) != 1L) {
        stop("Could not determine the comparison group.", call. = FALSE)
    }
    list(
        g = as.integer(key != reference_key),
        reference = reference_key,
        comparison = comparison_key
    )
}

.dasra_validate_design <- function(g, z) {
    n <- length(g)
    restricted <- count_design_matrix(g = g, z = z, include_group = FALSE,
                                      n = n)
    full <- count_design_matrix(g = g, z = z, include_group = TRUE, n = n)
    if (qr(restricted)$rank < ncol(restricted)) {
        stop("The restricted covariate design is rank deficient.",
             call. = FALSE)
    }
    if (qr(full)$rank < ncol(full)) {
        stop("The group-adjusted covariate design is rank deficient.",
             call. = FALSE)
    }
    invisible(NULL)
}

# Numerical failures are recorded at the taxon level so that other taxa can
# still be analyzed.
.dasra_structural_arm <- function(Y, N, g, z, keep_diagnostics) {
    J <- ncol(Y)
    n_samples <- nrow(Y)
    fits <- vector("list", J)
    for (j in seq_len(J)) {
        fits[[j]] <- tryCatch(
            zt_count_structural_test(
                y = as.numeric(Y[, j]),
                N = N,
                g = g,
                z = z,
                Q = 1001L,
                min_positive_samples = 3L,
                derivative_base = 1e-4,
                derivative_reference_n = 120L,
                keep_fit = keep_diagnostics
            ),
            error = function(e) {
                zt_unavailable(
                    paste0("structural_fit_error: ", conditionMessage(e)),
                    n_samples
                )
            }
        )
    }
    names(fits) <- colnames(Y)
    p <- vapply(fits, function(x) as.numeric(x$p)[1L], numeric(1))
    formed <- vapply(fits, function(x) {
        isTRUE(x$tested) && length(x$p) == 1L && is.finite(x$p) &&
            x$p >= 0 && x$p <= 1
    }, logical(1))
    reason <- vapply(fits, function(x) {
        if (isTRUE(x$tested) && length(x$p) == 1L && is.finite(x$p) &&
            x$p >= 0 && x$p <= 1) {
            "ok"
        } else {
            as.character(x$reason)[1L]
        }
    }, character(1))
    p[!formed] <- 1
    score_z <- vapply(fits, function(x) {
        if (!isTRUE(x$tested)) return(NA_real_)
        U <- as.numeric(x$diagnostics$U)
        V <- as.numeric(x$diagnostics$V)
        if (length(U) != 1L || length(V) != 1L || !is.finite(U) ||
            !is.finite(V) || V <= 0) {
            return(NA_real_)
        }
        U / sqrt(V)
    }, numeric(1))
    list(
        p = p,
        formed = formed,
        reason = reason,
        score_z = score_z,
        diagnostics = if (keep_diagnostics) fits else NULL
    )
}

#' Depth-Aware Structural-Absence and Relative-Abundance Analysis
#'
#' Fits DASRA to microbiome count data and tests a binary group contrast in
#' structural-absence probability, relative abundance conditional on presence,
#' or both. Sequencing depth enters the latent-state count model so that an
#' observed zero can receive different posterior support for structural absence
#' and nondetection across samples.
#'
#' The structural-absence component estimates the conditional-present count
#' distribution from positive counts, estimates structural-absence nuisance
#' parameters from the detection indicators, and evaluates a nuisance-adjusted
#' score for the group effect. Positive `z_structural_absence` values indicate
#' greater structural absence in the comparison group.
#'
#' The relative-abundance component estimates a covariate-standardized
#' comparison-minus-reference difference in mean log relative abundance
#' conditional on taxon presence. Zero-count samples enter its estimating
#' equations through their posterior probability of presence and their
#' posterior latent-abundance score. The taxon-specific effects are centered
#' against a target-excluded robust cross-taxon reference background, with
#' sample-aligned influence contributions used for the corrected standard error.
#' A positive `estimate_relative_abundance` means that the target taxon's
#' present-conditional log-relative-abundance contrast exceeds this shared
#' compositional background. This correction is designed for analyses in which
#' a strict majority of eligible taxa share a common compositional background;
#' the reported abundance effect is a reference-centered relative contrast.
#'
#' Taxa with positive counts in fewer than three samples are excluded before
#' component fitting and omitted from the multiple-testing families. Every
#' retained taxon remains in each requested family; an unavailable component is
#' assigned a conservative p-value of one.
#'
#' @param counts Raw non-negative integer counts in a matrix or data frame.
#' @param metadata Sample metadata in a data frame. Row names must contain all
#'   sample names in `counts`.
#' @param formula A one-sided formula containing `group` as an additive main
#'   effect and any adjustment terms, for example `~ disease + age + sex`.
#'   Interactions involving `group`, offsets, and random-effect terms are not
#'   supported.
#' @param group Name of the binary metadata variable to test.
#' @param library_size Either the name of a numeric column in `metadata` or a
#'   numeric vector named by sample containing the original sequencing depth
#'   for each sample.
#' @param taxa_are_rows Logical. If `TRUE`, taxa are rows and samples are
#'   columns in `counts`. If `FALSE`, samples are rows.
#' @param reference Optional reference value for `group`. A factor uses its
#'   first level by default; logical and numeric 0/1 groups use `FALSE` or 0.
#'   Character groups and other numeric codings require an explicit reference.
#' @param p_adjust_method Method passed to [stats::p.adjust()] separately for
#'   each requested component and omnibus family.
#' @param component Analysis to run. `"all"` fits both components and reports
#'   the omnibus analyses. `"structural_absence"` fits the structural-absence
#'   component. `"relative_abundance"` fits the relative-abundance component.
#' @param full_output Logical. If `TRUE`, the returned object includes detailed
#'   fits for the requested components.
#'
#' @return An object of class `dasra` with elements:
#'   \describe{
#'     \item{results}{A taxon-level table for the requested analyses.}
#'     \item{diagnostics}{Taxon retention, support counts, component-formation
#'       status, and machine-readable formation reasons.}
#'     \item{settings}{The fitted contrast and analysis settings.}
#'     \item{call}{The matched function call.}
#'     \item{fits}{Detailed component fits when `full_output = TRUE`.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(2026)
#' n <- 80
#' taxa <- paste0("Taxon_", seq_len(8))
#' samples <- paste0("Sample_", seq_len(n))
#' group <- rep(c(0, 1), each = n / 2)
#' library_size <- sample(seq(8000L, 12000L, by = 500L), n, replace = TRUE)
#'
#' baseline <- seq(-6.0, -5.2, length.out = length(taxa))
#' abundance_shift <- c(0.25, -0.20, rep(0, length(taxa) - 2L))
#' absence_reference <- c(0.15, 0.20, 0.18, 0.22, 0.16, 0.24, 0.19, 0.21)
#' absence_comparison <- c(0.30, 0.20, 0.10, 0.22, 0.16, 0.24, 0.19, 0.21)
#'
#' probability <- matrix(0, nrow = n, ncol = length(taxa))
#' for (j in seq_along(taxa)) {
#'   absent_probability <- ifelse(
#'     group == 0, absence_reference[j], absence_comparison[j]
#'   )
#'   present <- runif(n) > absent_probability
#'   latent_abundance <- baseline[j] + abundance_shift[j] * group +
#'     rnorm(n, sd = 0.35)
#'   probability[, j] <- present * plogis(latent_abundance)
#' }
#'
#' count_by_sample <- t(vapply(seq_len(n), function(i) {
#'   draw <- rmultinom(
#'     1,
#'     size = library_size[i],
#'     prob = c(probability[i, ], 1 - sum(probability[i, ]))
#'   )
#'   draw[seq_along(taxa), 1]
#' }, numeric(length(taxa))))
#' counts <- t(count_by_sample)
#' rownames(counts) <- taxa
#' colnames(counts) <- samples
#'
#' metadata <- data.frame(
#'   group = factor(group, levels = c(0, 1), labels = c("control", "case")),
#'   reads = library_size,
#'   row.names = samples
#' )
#'
#' fit <- dasra(
#'   counts = counts,
#'   metadata = metadata,
#'   formula = ~ group,
#'   group = "group",
#'   library_size = "reads",
#'   component = "all"
#' )
#' fit
#' head(fit$results)
#' }
#'
#' @export
dasra <- function(counts, metadata, formula, group, library_size,
                  taxa_are_rows = TRUE, reference = NULL,
                  p_adjust_method = "BH",
                  component = c("all", "structural_absence",
                                "relative_abundance"),
                  full_output = FALSE) {
    call <- match.call()
    component <- match.arg(component)
    raw_counts <- as.matrix(counts)
    if (length(dim(raw_counts)) != 2L || !is.numeric(raw_counts) ||
        nrow(raw_counts) < 2L || ncol(raw_counts) < 2L) {
        stop("`counts` must be a numeric matrix or data frame with at least two rows and columns.",
             call. = FALSE)
    }
    if (any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
        any(abs(raw_counts - round(raw_counts)) > 1e-8)) {
        stop("`counts` must contain finite, non-negative integer counts.",
             call. = FALSE)
    }
    if (length(taxa_are_rows) != 1L || is.na(taxa_are_rows) ||
        !is.logical(taxa_are_rows)) {
        stop("`taxa_are_rows` must be TRUE or FALSE.", call. = FALSE)
    }
    if (isTRUE(taxa_are_rows)) {
        taxa <- rownames(raw_counts)
        samples <- colnames(raw_counts)
        Y <- t(raw_counts)
    } else {
        taxa <- colnames(raw_counts)
        samples <- rownames(raw_counts)
        Y <- raw_counts
    }
    if (is.null(taxa) || anyNA(taxa) || any(!nzchar(taxa)) ||
        anyDuplicated(taxa)) {
        stop("Taxon names must be unique and non-empty.", call. = FALSE)
    }
    if (is.null(samples) || anyNA(samples) || any(!nzchar(samples)) ||
        anyDuplicated(samples)) {
        stop("Sample names in `counts` must be unique and non-empty.",
             call. = FALSE)
    }
    colnames(Y) <- taxa
    rownames(Y) <- samples
    if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata)
    metadata_columns <- colnames(metadata)
    if (is.null(metadata_columns) || anyNA(metadata_columns) ||
        any(!nzchar(metadata_columns)) || anyDuplicated(metadata_columns)) {
        stop("`metadata` column names must be unique and non-empty.",
             call. = FALSE)
    }
    metadata_names <- rownames(metadata)
    if (is.null(metadata_names) || anyNA(metadata_names) ||
        any(!nzchar(metadata_names)) || anyDuplicated(metadata_names)) {
        stop("`metadata` must have unique, non-empty sample row names.",
             call. = FALSE)
    }
    missing_metadata <- setdiff(samples, metadata_names)
    if (length(missing_metadata)) {
        shown <- paste(
            missing_metadata[seq_len(min(5L, length(missing_metadata)))],
            collapse = ", "
        )
        stop(sprintf("Metadata are missing for count samples: %s%s",
                     shown,
                     if (length(missing_metadata) > 5L) ", ..." else ""),
             call. = FALSE)
    }
    extra_metadata <- setdiff(metadata_names, samples)
    if (length(extra_metadata)) {
        warning(
            sprintf("Ignoring %d metadata row(s) not present in `counts`.",
                    length(extra_metadata)),
            call. = FALSE
        )
    }
    metadata <- metadata[samples, , drop = FALSE]
    if (!inherits(formula, "formula") || length(formula) != 2L) {
        stop("`formula` must be a one-sided formula.", call. = FALSE)
    }
    if (length(group) != 1L || is.na(group) || !is.character(group) ||
        !(group %in% colnames(metadata))) {
        stop("`group` must name one column of `metadata`.", call. = FALSE)
    }
    if (grepl("\\|", paste(deparse(formula), collapse = ""))) {
        stop("Random-effect terms are not supported.", call. = FALSE)
    }
    terms_object <- stats::terms(formula, data = metadata)
    if (length(attr(terms_object, "offset"))) {
        stop("Offset terms are not supported in `formula`.", call. = FALSE)
    }
    missing_formula_variables <- setdiff(
        all.vars(terms_object), colnames(metadata)
    )
    if (length(missing_formula_variables)) {
        stop(
            sprintf(
                "Every variable in `formula` must be a column of `metadata`; missing: %s.",
                paste(missing_formula_variables, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    if (!identical(as.integer(attr(terms_object, "intercept")), 1L)) {
        stop("`formula` must retain its intercept.", call. = FALSE)
    }
    term_labels <- attr(terms_object, "term.labels")
    factor_map <- attr(terms_object, "factors")
    unquote_name <- function(x) sub("^`(.*)`$", "\\1", x)
    group_rows <- which(unquote_name(rownames(factor_map)) == group)
    group_terms <- if (length(group_rows) == 1L) {
        which(factor_map[group_rows, ] != 0)
    } else {
        integer()
    }
    main_group_terms <- group_terms[
        colSums(factor_map[, group_terms, drop = FALSE] != 0) == 1L &
            unquote_name(term_labels[group_terms]) == group
    ]
    if (length(main_group_terms) != 1L) {
        stop("`formula` must contain `group` as one additive main effect.",
             call. = FALSE)
    }
    if (length(setdiff(group_terms, main_group_terms))) {
        stop("Interactions involving `group` are not supported.",
             call. = FALSE)
    }
    model_frame <- tryCatch(
        stats::model.frame(
            terms_object, data = metadata, na.action = stats::na.fail,
            drop.unused.levels = TRUE
        ),
        error = function(e) {
            stop(sprintf("Could not construct the model from `metadata`: %s",
                         conditionMessage(e)), call. = FALSE)
        }
    )
    model_matrix <- stats::model.matrix(terms_object, data = model_frame)
    assignment <- attr(model_matrix, "assign")
    group_term <- main_group_terms
    z <- model_matrix[
        , !(assignment %in% c(0L, group_term)), drop = FALSE
    ]
    if (!ncol(z)) z <- NULL
    if (!is.null(z) && any(!is.finite(z))) {
        stop("The covariate model matrix contains non-finite values.",
             call. = FALSE)
    }
    group_info <- .dasra_encode_group(metadata[[group]], reference, nrow(Y))
    g <- group_info$g
    .dasra_validate_design(g, z)
    if (is.character(library_size)) {
        if (length(library_size) != 1L ||
            !(library_size %in% colnames(metadata))) {
            stop("Character `library_size` must name one metadata column.",
                 call. = FALSE)
        }
        depth_column <- metadata[[library_size]]
        if (!is.numeric(depth_column) || is.factor(depth_column) ||
            is.logical(depth_column) || inherits(depth_column, "Date") ||
            !is.null(dim(depth_column))) {
            stop("The metadata column named by `library_size` must be a numeric vector.",
                 call. = FALSE)
        }
        N <- as.numeric(depth_column)
        library_size_source <- paste0("metadata$", library_size)
    } else {
        if (!is.numeric(library_size) || is.factor(library_size) ||
            is.logical(library_size) || inherits(library_size, "Date") ||
            !is.null(dim(library_size))) {
            stop("Numeric `library_size` must be a vector named by sample.",
                 call. = FALSE)
        }
        depth_names <- names(library_size)
        if (is.null(depth_names) ||
            length(depth_names) != length(library_size) ||
            anyNA(depth_names) || any(!nzchar(depth_names)) ||
            anyDuplicated(depth_names)) {
            stop("Numeric `library_size` must have unique, non-empty sample names.",
                 call. = FALSE)
        }
        if (!all(samples %in% depth_names)) {
            stop("Named `library_size` is missing count samples.",
                 call. = FALSE)
        }
        N <- as.numeric(library_size[samples])
        library_size_source <- "numeric vector"
    }
    if (length(N) != nrow(Y) || any(!is.finite(N)) || any(N <= 0) ||
        any(abs(N - round(N)) > 1e-8)) {
        stop("`library_size` must provide one finite, positive integer per sample.",
             call. = FALSE)
    }
    if (any(N < rowSums(Y))) {
        stop("`library_size` cannot be smaller than the corresponding count total.",
             call. = FALSE)
    }
    if (length(p_adjust_method) != 1L || is.na(p_adjust_method) ||
        !(p_adjust_method %in% stats::p.adjust.methods)) {
        stop("`p_adjust_method` must be one of `stats::p.adjust.methods`.",
             call. = FALSE)
    }
    if (length(full_output) != 1L || is.na(full_output) ||
        !is.logical(full_output)) {
        stop("`full_output` must be TRUE or FALSE.", call. = FALSE)
    }
    full_output <- isTRUE(full_output)
    run_structural <- component %in% c("all", "structural_absence")
    run_abundance <- component %in% c("all", "relative_abundance")
    run_omnibus <- identical(component, "all")

    J <- ncol(Y)
    retained <- colSums(Y > 0) >= 3L

    structural <- NULL
    if (run_structural) {
        structural <- list(
            p = rep(NA_real_, J),
            formed = rep(FALSE, J),
            reason = rep("fewer_than_three_positive_samples", J),
            score_z = rep(NA_real_, J),
            diagnostics = NULL
        )
    }

    abundance <- NULL
    if (run_abundance) {
        abundance <- list(
            p = rep(NA_real_, J),
            formed = rep(FALSE, J),
            reason = rep("fewer_than_three_positive_samples", J),
            estimate = rep(NA_real_, J),
            se = rep(NA_real_, J),
            z = rep(NA_real_, J),
            diagnostics = NULL
        )
    }

    if (any(retained)) {
        Y_retained <- Y[, retained, drop = FALSE]

        if (run_structural) {
            structural_retained <- .dasra_structural_arm(
                Y_retained, N, g, z, full_output
            )
            structural$p[retained] <- structural_retained$p
            structural$formed[retained] <- structural_retained$formed
            structural$reason[retained] <- structural_retained$reason
            structural$score_z[retained] <- structural_retained$score_z
            structural$diagnostics <- structural_retained$diagnostics
        }

        if (run_abundance) {
            abundance_retained <- .dasra_abundance_arm(
                Y_retained, N, g, z, full_output
            )
            abundance$p[retained] <- abundance_retained$p
            abundance$formed[retained] <- abundance_retained$formed
            abundance$reason[retained] <- abundance_retained$reason
            abundance$estimate[retained] <- abundance_retained$estimate
            abundance$se[retained] <- abundance_retained$se
            abundance$z[retained] <- abundance_retained$z
            abundance$diagnostics <- abundance_retained$diagnostics
        }
    }

    adjust_family <- function(p_values) {
        adjusted <- rep(NA_real_, J)
        family <- retained & is.finite(p_values)
        if (any(family)) {
            adjusted[family] <- stats::p.adjust(
                p_values[family], method = p_adjust_method
            )
        }
        adjusted
    }

    results <- data.frame(
        taxon = taxa,
        retained = retained,
        row.names = taxa,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (run_omnibus) {
        p_omnibus <- rep(NA_real_, J)
        p_omnibus_cauchy <- rep(NA_real_, J)
        retained_index <- which(retained)
        if (length(retained_index)) {
            p_omnibus[retained_index] <- pmin(
                1,
                2 * pmin(
                    structural$p[retained_index],
                    abundance$p[retained_index]
                )
            )
            p_omnibus_cauchy[retained_index] <- vapply(
                retained_index,
                function(j) {
                    if (!abundance$formed[j]) {
                        structural$p[j]
                    } else {
                        cauchy_combination(
                            c(structural$p[j], abundance$p[j])
                        )
                    }
                },
                numeric(1)
            )
        }
        components_used <- rep("not_retained", J)
        components_used[retained] <- ifelse(
            structural$formed[retained] & abundance$formed[retained],
            "both",
            ifelse(
                structural$formed[retained],
                "structural_absence",
                ifelse(
                    abundance$formed[retained],
                    "relative_abundance",
                    "none"
                )
            )
        )
        results$p_omnibus <- p_omnibus
        results$q_omnibus <- adjust_family(p_omnibus)
        results$p_omnibus_cauchy <- p_omnibus_cauchy
        results$q_omnibus_cauchy <- adjust_family(p_omnibus_cauchy)
        results$components_used <- components_used
    }

    if (run_structural) {
        results$p_structural_absence <- structural$p
        results$q_structural_absence <- adjust_family(structural$p)
        results$z_structural_absence <- structural$score_z
    }

    if (run_abundance) {
        results$p_relative_abundance <- abundance$p
        results$q_relative_abundance <- adjust_family(abundance$p)
        results$estimate_relative_abundance <- abundance$estimate
        results$se_relative_abundance <- abundance$se
        results$z_relative_abundance <- abundance$z
    }

    diagnostics <- data.frame(
        taxon = taxa,
        retained = retained,
        n_positive_reference = colSums(Y[g == 0, , drop = FALSE] > 0),
        n_positive_comparison = colSums(Y[g == 1, , drop = FALSE] > 0),
        n_zero = colSums(Y == 0),
        row.names = taxa,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (run_structural) {
        diagnostics$formed_structural_absence <- structural$formed
        diagnostics$reason_structural_absence <- structural$reason
    }
    if (run_abundance) {
        diagnostics$formed_relative_abundance <- abundance$formed
        diagnostics$reason_relative_abundance <- abundance$reason
    }
    if (run_omnibus) {
        diagnostics$formed_omnibus <-
            structural$formed | abundance$formed
    }

    settings <- list(
        formula = formula,
        group = group,
        contrast = c(
            reference = group_info$reference,
            comparison = group_info$comparison
        ),
        n_samples = nrow(Y),
        n_taxa = ncol(Y),
        n_taxa_retained = sum(retained),
        taxa_are_rows = isTRUE(taxa_are_rows),
        library_size_source = library_size_source,
        p_adjust_method = p_adjust_method,
        component = component,
        min_positive_samples_retained = 3L
    )
    if (run_structural) {
        settings$structural_method <- paste(
            "two-stage zero-truncated/detection stacked score",
            "with adaptive Gauss-Hermite quadrature"
        )
        settings$structural_quadrature_Q <- 1001L
        settings$structural_derivative_base <- 1e-4
        settings$structural_derivative_reference_n <- 120L
    }
    if (run_abundance) {
        settings$relative_abundance_method <- paste(
            "depth-aware present-conditional mean-log-relative-abundance estimation",
            "with target-excluded robust compositional correction"
        )
        settings$abundance_quadrature_Q <- 41L
        settings$abundance_effect_quadrature_Q <- 41L
        settings$abundance_root_tolerance <- 1e-5
        settings$abundance_reference_method <- paste(
            "target-excluded intercept-only least-trimmed-squares pilot",
            "with studentized reference expansion"
        )
    }
    if (run_omnibus) {
        settings$omnibus <- "Bonferroni minimum-p"
        settings$sensitivity_omnibus <- "equal-weight Cauchy combination"
    }

    object <- list(
        results = results,
        diagnostics = diagnostics,
        settings = settings,
        call = call
    )
    if (full_output) {
        object$fits <- list()
        if (run_structural) {
            object$fits$structural_absence <- structural$diagnostics
        }
        if (run_abundance) {
            object$fits$relative_abundance <- abundance$diagnostics
        }
    }
    class(object) <- "dasra"
    object
}

#' @export
#' @noRd
print.dasra <- function(x, ...) {
    contrast <- x$settings$contrast
    retained_n <- sum(x$diagnostics$retained)
    adjustment_label <- if (identical(x$settings$p_adjust_method, "BH")) {
        "BH q"
    } else {
        paste0(x$settings$p_adjust_method, " adjusted p")
    }

    cat(
        "DASRA fit\n",
        sprintf(
            "  Samples: %d; taxa: %d; retained taxa: %d\n",
            x$settings$n_samples,
            x$settings$n_taxa,
            retained_n
        ),
        sprintf(
            "  Contrast: %s - %s\n",
            unname(contrast["comparison"]),
            unname(contrast["reference"])
        ),
        sep = ""
    )

    if ("formed_structural_absence" %in% names(x$diagnostics)) {
        formed <- sum(x$diagnostics$formed_structural_absence)
        discoveries <- sum(
            x$results$q_structural_absence <= 0.05,
            na.rm = TRUE
        )
        cat(
            sprintf("  Structural-absence tests formed: %d/%d\n",
                    formed, retained_n),
            sprintf("  Structural-absence discoveries (%s <= 0.05): %d\n",
                    adjustment_label, discoveries),
            sep = ""
        )
    }

    if ("formed_relative_abundance" %in% names(x$diagnostics)) {
        formed <- sum(x$diagnostics$formed_relative_abundance)
        discoveries <- sum(
            x$results$q_relative_abundance <= 0.05,
            na.rm = TRUE
        )
        cat(
            sprintf("  Relative-abundance tests formed: %d/%d\n",
                    formed, retained_n),
            sprintf("  Relative-abundance discoveries (%s <= 0.05): %d\n",
                    adjustment_label, discoveries),
            sep = ""
        )
    }

    if ("formed_omnibus" %in% names(x$diagnostics)) {
        discoveries <- sum(x$results$q_omnibus <= 0.05, na.rm = TRUE)
        cat(sprintf("  Omnibus discoveries (%s <= 0.05): %d\n",
                    adjustment_label, discoveries))
    }

    cat("  Results: x$results; formation details: x$diagnostics\n")
    invisible(x)
}

#' @export
#' @noRd
as.data.frame.dasra <- function(x, row.names = NULL, optional = FALSE,
                                ...) {
    as.data.frame(
        x$results,
        row.names = row.names,
        optional = optional,
        ...
    )
}
