test_that("Rcpp adaptive abundance marginals agree with the R implementation", {
    gh <- DASRA:::make_count_gh_rule(21L)
    y <- c(0, 1, 3, 8)
    N <- c(100, 150, 200, 300)
    eta <- c(-5.0, -4.7, -4.2, -3.8)
    sigma <- 0.8

    cpp <- DASRA:::dasra_count_moments_adaptive_cpp(
        y, N, eta, sigma, gh$node, gh$log_raw_weight, TRUE
    )
    r <- DASRA:::count_log_hy_adaptive_R(
        y, N, eta, sigma, gh, return_nodes = TRUE
    )
    weight <- exp(r$log_terms - r$log_term_normalizer)
    centered <- r$x_node - matrix(
        eta, nrow = length(eta), ncol = ncol(r$x_node)
    )
    m1 <- rowSums(weight * centered)
    m2 <- rowSums(weight * centered^2)

    expect_equal(cpp$log_h, r$log_hy, tolerance = 1e-10)
    expect_equal(cpp$posterior_m1, m1, tolerance = 1e-10)
    expect_equal(cpp$posterior_m2, m2, tolerance = 1e-10)
    expect_equal(cpp$score_eta, m1 / sigma^2, tolerance = 1e-10)
    expect_equal(cpp$score_omega, m2 / sigma^2 - 1, tolerance = 1e-10)
})

test_that("Rcpp mean log relative abundance agrees with direct quadrature", {
    gh <- DASRA:::make_count_gh_rule(31L)
    location <- c(-6, -4, -2)
    sigma <- 0.7
    direct <- vapply(location, function(mu) {
        sum(gh$weight * (-DASRA:::count_softplus(-(mu + sigma * sqrt(2) * gh$node))))
    }, numeric(1))
    cpp <- DASRA:::dasra_mean_log_relative_cpp(
        location, sigma, sqrt(2) * gh$node, gh$weight
    )
    expect_equal(cpp, direct, tolerance = 1e-12)
})

test_that("intercept-only LTS uses a strict majority subset", {
    values <- c(-0.52, -0.50, -0.48, -0.46, 0.01, 0.04, 0.06)
    pilot <- DASRA:::.dasra_abundance_lts_pilot(values)
    expect_lt(abs(pilot + 0.49), 0.04)
})

test_that("public input validation rejects unsupported group interactions", {
    counts <- matrix(
        c(3, 2, 4, 5, 6, 7, 8, 9, 2, 3, 4, 5),
        nrow = 3,
        dimnames = list(
            paste0("Taxon_", 1:3), paste0("Sample_", 1:4)
        )
    )
    metadata <- data.frame(
        group = factor(c("a", "a", "b", "b")),
        age = c(20, 30, 40, 50),
        reads = colSums(counts) + 100,
        row.names = colnames(counts)
    )
    expect_error(
        dasra(
            counts, metadata, ~ group * age, "group", "reads",
            component = "relative_abundance"
        ),
        "Interactions involving `group`"
    )
})

test_that("relative-abundance component forms on a regular example", {
    skip_on_cran()
    set.seed(123)
    n <- 72
    p <- 6
    group <- rep(c(0, 1), each = n / 2)
    depth <- rep(10000L, n)
    latent <- sapply(seq_len(p), function(j) {
        -5.7 + (j - 1) * 0.12 + c(0.25, -0.20, rep(0, p - 2))[j] * group +
            rnorm(n, sd = 0.25)
    })
    present <- matrix(runif(n * p) > 0.15, nrow = n)
    probability <- present * plogis(latent)
    count_by_sample <- t(vapply(seq_len(n), function(i) {
        draw <- rmultinom(
            1, depth[i], c(probability[i, ], 1 - sum(probability[i, ]))
        )
        draw[seq_len(p), 1]
    }, numeric(p)))
    counts <- t(count_by_sample)
    rownames(counts) <- paste0("Taxon_", seq_len(p))
    colnames(counts) <- paste0("Sample_", seq_len(n))
    metadata <- data.frame(
        group = factor(group), reads = depth,
        row.names = colnames(counts)
    )

    fit <- dasra(
        counts, metadata, ~ group, "group", "reads",
        component = "relative_abundance"
    )
    expect_s3_class(fit, "dasra")
    expect_true(all(fit$diagnostics$retained))
    expect_true(all(fit$diagnostics$formed_relative_abundance))
    expect_true(all(is.finite(fit$results$p_relative_abundance)))
})
