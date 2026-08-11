test_that("the publication API has one compact fitting entry point", {
  expected_formals <- as.pairlist(alist(
    counts = ,
    metadata = ,
    group = ,
    library_size = ,
    covariates = NULL,
    reference = NULL,
    taxa = NULL,
    cores = 1L,
    verbose = FALSE
  ))
  expect_identical(formals(DORAM), expected_formals)
  exports <- getNamespaceExports("DORAM")
  expect_true(
    length(exports) == 0L || identical(exports, "DORAM"),
    info = "a documented installation must export only DORAM()"
  )

  namespace <- asNamespace("DORAM")
  expect_false(exists("doram", envir = namespace, inherits = FALSE))
  expect_false(exists("DORAM_contract", envir = namespace, inherits = FALSE))
  expect_false(exists("coef.DORAM", envir = namespace, inherits = FALSE))
  expect_false(exists("predict.DORAM", envir = namespace, inherits = FALSE))
  expect_false(exists("summary.DORAM", envir = namespace, inherits = FALSE))
  expect_true(exists("print.DORAM", envir = namespace, inherits = FALSE))
  expect_true(exists("as.data.frame.DORAM", envir = namespace, inherits = FALSE))
})

test_that("documentation collapses to one user-facing man page after document", {
  package_root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
  man_directory <- file.path(package_root, "man")
  skip_if_not(
    dir.exists(man_directory),
    "run devtools::document() before checking the generated man-page surface"
  )
  pages <- sort(list.files(man_directory, pattern = "[.]Rd$"))
  expect_identical(pages, "DORAM.Rd")
})

test_that("the count-likelihood and analytic-inference contracts are compatible", {
  expect_identical(
    DORAM:::CN_D3_CONTRACT$schema_version,
    "doram-count-likelihood-1.0.0"
  )
  expect_identical(
    DORAM:::D3A_CONTRACT$schema,
    "doram-analytic-sandwich-1.0.0"
  )
  expect_identical(
    DORAM:::D3A_CONTRACT$calibration_label,
    "asymptotic_empirical_sandwich_chisq"
  )
  expect_identical(DORAM:::D3A_CONTRACT$primary_reference, "analytic_chisq")
  expect_true(DORAM:::d3a_count_source_is_compatible())
})

test_that("inference contracts retain the required numerical safeguards", {
  bounds <- DORAM:::CN_D3_CONTRACT$parameter_bounds
  expect_identical(names(bounds), c("structural", "location", "log_scale"))
  expect_true(all(vapply(bounds, function(x) {
    is.numeric(x) && identical(names(x), c("lower", "upper")) &&
      length(x) == 2L && all(is.finite(x)) && x[[1L]] < x[[2L]]
  }, logical(1L))))
  expect_true(
    DORAM:::CN_D3_CONTRACT$optimizer$boundary_relative_distance > 0
  )
  expect_identical(
    DORAM:::D3A_CONTRACT$count_source_schema,
    DORAM:::CN_D3_CONTRACT$schema_version
  )
})

test_that("the production endpoint call remains count-native and analytic", {
  x <- doram_toy_data(taxa = "target")
  original <- DORAM:::d3a_fit_endpoints
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    d3a_fit_endpoints = function(count_data, X, strata,
                                 sampling_design, integration_level) {
      captured$count_data <- count_data
      captured$X <- X
      captured$strata <- strata
      captured$sampling_design <- sampling_design
      captured$integration_level <- integration_level
      original(
        count_data = count_data, X = X, strata = strata,
        sampling_design = sampling_design,
        integration_level = integration_level
      )
    },
    .package = "DORAM"
  )

  doram_fit_toy(x, taxa = "target", covariates = "adjustment")

  expect_s3_class(captured$count_data, "cn_d3_data")
  expect_identical(captured$count_data$subject_id, x$ids)
  expect_equal(
    captured$count_data$N,
    x$metadata$original_library_size,
    tolerance = 0
  )
  expect_equal(captured$count_data$g, x$metadata$group, tolerance = 0)
  expect_equal(
    unname(captured$count_data$Xrho[, -1L, drop = FALSE]),
    unname(as.matrix(x$metadata[, "adjustment", drop = FALSE])),
    tolerance = 0
  )
  expect_equal(
    unname(captured$X[, -1L, drop = FALSE]),
    unname(as.matrix(x$metadata[, "adjustment", drop = FALSE])),
    tolerance = 0
  )
  expect_null(captured$strata)
  expect_identical(captured$sampling_design, "iid_random_design")
  expect_identical(captured$integration_level, "production")
})

test_that("the analytic inference call graph contains no resampling path", {
  prohibited <- c(
    "cn_occupancy_score_test", "cn_score_geometry_1d", "cn_multiplier_1d",
    "d3_ecological_score", "d3_joint_score", "d3_multiplier_test",
    "sample", "sample.int", "rbinom", "set.seed"
  )
  inference_functions <- list(
    DORAM:::d3a_quadratic,
    DORAM:::d3a_center_rows,
    DORAM:::d3a_score_test,
    DORAM:::d3a_fwl_rows,
    DORAM:::d3a_projected_endpoint_test,
    DORAM:::d3a_joint_test,
    DORAM:::d3a_fit_endpoints
  )
  used <- unique(unlist(lapply(inference_functions, function(fun) {
    intersect(all.names(body(fun), functions = TRUE), prohibited)
  })))

  expect_length(used, 0L)
  expect_false(exists("cn_occupancy_score_test", envir = asNamespace("DORAM"),
                      inherits = FALSE))
  expect_false(exists("cn_multiplier_1d", envir = asNamespace("DORAM"),
                      inherits = FALSE))
})

test_that("one-dimensional iid score geometry agrees with hand calculations", {
  phi <- matrix(c(1, 2, 4, 8, 16), ncol = 1L)
  ids <- sprintf("S%02d", seq_len(nrow(phi)))
  centered <- phi[, 1L] - mean(phi[, 1L])
  expected_statistic <- sum(phi)^2 / sum(centered^2)

  fit <- DORAM:::d3a_score_test(
    phi, ids, sampling_design = "iid_random_design",
    endpoint = "occupancy"
  )

  expect_true(fit$available)
  expect_identical(fit$df, 1L)
  expect_equal(fit$statistic, expected_statistic, tolerance = 1e-12)
  expect_equal(
    fit$p_value,
    stats::pchisq(expected_statistic, 1, lower.tail = FALSE),
    tolerance = 1e-12
  )

  permutation <- c(5, 2, 4, 1, 3)
  permuted <- DORAM:::d3a_score_test(
    phi[permutation, , drop = FALSE], ids[permutation],
    sampling_design = "iid_random_design", endpoint = "occupancy"
  )
  expect_equal(permuted$statistic, fit$statistic, tolerance = 1e-12)
  expect_equal(permuted$p_value, fit$p_value, tolerance = 1e-12)
})

test_that("two-dimensional iid score geometry uses globally centered HC0 meat", {
  phi <- cbind(
    occupancy = c(1, 2, 7, 10, 14, 23),
    abundance = c(3, 8, 4, 17, 11, 29)
  )
  ids <- sprintf("J%02d", seq_len(nrow(phi)))
  centered <- sweep(phi, 2L, colMeans(phi), FUN = "-")
  expected_U <- colSums(phi)
  expected_B <- crossprod(centered)
  expected_statistic <- drop(
    t(expected_U) %*% solve(expected_B, expected_U)
  )

  fit <- DORAM:::d3a_score_test(
    phi, ids, sampling_design = "iid_random_design", endpoint = "joint"
  )

  expect_true(fit$available)
  expect_identical(fit$df, 2L)
  expect_equal(fit$U, expected_U, tolerance = 1e-12)
  expect_equal(fit$B, expected_B, tolerance = 1e-12)
  expect_equal(fit$statistic, expected_statistic, tolerance = 1e-11)
})
