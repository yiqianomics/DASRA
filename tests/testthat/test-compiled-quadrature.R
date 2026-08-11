compare_present_integrators <- function(case, level = "production",
                                        tolerance = 2e-12) {
  reference <- DORAM:::cn_integrated_present(
    case[["y"]], case[["N"]], case[["eta"]], case[["sigma"]],
    level = level
  )
  compiled <- DORAM:::cn_integrated_present_compiled(
    case[["y"]], case[["N"]], case[["eta"]], case[["sigma"]],
    level = level
  )

  scalar_fields <- c(
    "log_K", "score_eta", "score_log_sigma", "mode", "mode_residual",
    "mode_backward_error", "scale"
  )
  vector_fields <- c(
    "integration_absolute_error", "integration_transformed_error",
    "integration_tail_bound"
  )
  for (field in c(scalar_fields, vector_fields)) {
    reference_value <- unname(reference[[field]])
    compiled_value <- unname(compiled[[field]])
    scaled_difference <- abs(reference_value - compiled_value) /
      pmax(1, abs(reference_value), abs(compiled_value))
    expect_true(
      all(is.finite(scaled_difference)) &&
        max(scaled_difference) <= tolerance,
      info = paste("quadrature field", field, "at", paste(case, collapse = ","))
    )
  }
  expect_identical(compiled$integration_level, reference$integration_level)
  expect_identical(compiled$integration_panels, reference$integration_panels)
  expect_identical(
    compiled$integration_evaluations,
    reference$integration_evaluations
  )
  invisible(compiled)
}

test_that("compiled quadrature agrees with the R reference on regular and extreme inputs", {
  fixed <- rbind(
    c(y = 0, N = 1, eta = -20, sigma = 0.1),
    c(y = 1, N = 1, eta = 10, sigma = 0.1),
    c(y = 9, N = 10, eta = -20, sigma = 8),
    c(y = 10, N = 10, eta = 10, sigma = 8),
    c(y = 0, N = 1000, eta = -12, sigma = 0.1),
    c(y = 1000, N = 1000, eta = 6, sigma = 8),
    c(y = 999, N = 1000, eta = -8, sigma = 1),
    c(y = 500, N = 1000, eta = 0, sigma = 0.75),
    c(y = 1, N = 1e6, eta = 10, sigma = 8),
    c(y = 999999, N = 1e6, eta = -20, sigma = 0.1)
  )

  set.seed(20260810)
  N <- sample(c(1, 2, 10, 100, 10000, 1e6), 12L, replace = TRUE)
  random <- cbind(
    y = vapply(N, function(size) sample.int(size + 1, 1L) - 1, numeric(1L)),
    N = N,
    eta = runif(length(N), -19.5, 9.5),
    sigma = exp(runif(length(N), log(0.101), log(7.9)))
  )

  cases <- rbind(fixed, random)
  for (index in seq_len(nrow(cases))) {
    compare_present_integrators(cases[index, ], level = "production")
  }
  for (index in seq_len(nrow(fixed))) {
    compare_present_integrators(fixed[index, ], level = "audit")
  }
})

test_that("the larger compiled-quadrature parity panel is opt-in", {
  skip_if_not(
    identical(tolower(Sys.getenv("DORAM_RUN_RCPP_PARITY")), "true"),
    "set DORAM_RUN_RCPP_PARITY=true to run the larger parity panel"
  )

  set.seed(20260811)
  N <- sample(c(1, 2, 5, 10, 100, 1000, 10000, 1e6), 1000L,
              replace = TRUE)
  cases <- cbind(
    y = vapply(N, function(size) sample.int(size + 1, 1L) - 1, numeric(1L)),
    N = N,
    eta = runif(length(N), -20, 10),
    sigma = exp(runif(length(N), log(0.1), log(8)))
  )

  for (index in seq_len(nrow(cases))) {
    compare_present_integrators(cases[index, ], level = "production")
  }
})

test_that("compiled quadrature does not create or advance R's random-number state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    seed_on_entry <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed_on_entry, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  DORAM:::cn_integrated_present_compiled(
    y = 7, N = 4000, eta = -6.1, sigma = 0.8,
    level = "production"
  )
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(20260812)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  DORAM:::cn_integrated_present_compiled(
    y = 7, N = 4000, eta = -6.1, sigma = 0.8,
    level = "production"
  )
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})

test_that("a full DORAM fit does not create or advance R's random-number state", {
  n <- 16L
  ids <- sprintf("rng_sample_%02d", seq_len(n))
  y <- rep(c(0L, 0L, 3L, 5L, 8L, 12L, 18L, 25L), 2L)
  counts <- matrix(y, ncol = 1L, dimnames = list(ids, "target"))
  metadata <- data.frame(
    group = rep(0:1, each = n / 2L),
    original_library_size = 3000L + 17L * seq_len(n),
    row.names = ids
  )

  original_integrator <- DORAM:::cn_integrated_present_compiled
  compiled_calls <- 0L
  local_mocked_bindings(
    cn_integrated_present_compiled = function(...) {
      compiled_calls <<- compiled_calls + 1L
      original_integrator(...)
    },
    .package = "DORAM"
  )

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    seed_on_entry <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", seed_on_entry, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  fit_once <- function() {
    DORAM(
      counts = counts,
      metadata = metadata,
      group = "group",
      library_size = "original_library_size",
      taxa = "target",
      cores = 1L
    )
  }

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  fit_without_seed <- fit_once()
  expect_s3_class(fit_without_seed, "DORAM")
  expect_gt(compiled_calls, 0L)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  set.seed(20260813)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  fit_with_seed <- fit_once()
  expect_s3_class(fit_with_seed, "DORAM")
  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})
