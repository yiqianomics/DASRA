test_that("the metadata API reproduces the reference AL_n240 inference fixture", {
  skip_if_not(
    identical(tolower(Sys.getenv("DORAM_RUN_SLOW_TESTS")), "true"),
    "set DORAM_RUN_SLOW_TESTS=true to run the production-integration parity fit"
  )
  fixture <- readRDS(test_path("fixtures", "al_n240_rep1_input.rds"))
  metadata <- doram_fixture_metadata(fixture)

  expect_identical(
    sort(rownames(metadata)),
    sort(as.character(fixture$subject_id))
  )
  metadata <- metadata[rev(fixture$subject_id), , drop = FALSE]

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) seed_before <- get(".Random.seed", envir = .GlobalEnv)

  fit <- DORAM(
    counts = fixture$counts,
    metadata = metadata,
    group = "group",
    library_size = "original_library_size",
    covariates = colnames(fixture$covariates),
    cores = 1L
  )

  if (had_seed) {
    expect_identical(get(".Random.seed", envir = .GlobalEnv), seed_before)
  } else {
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  }

  expect_identical(fit$sample_id, as.character(fixture$subject_id))
  expect_identical(fit$contrast, c(reference = "0", comparison = "1"))
  expect_identical(fit$taxa, "target")
  expect_identical(fit$settings$p_adjust_method, "BH")
  expect_identical(fit$settings$family_size, 1L)

  long <- fit$details[
    match(c("occupancy", "abundance", "joint"), fit$details$endpoint),
  ]
  expected_statistic <- c(
    0.0684497828571301,
    6.42728676421847,
    9.09887252710389
  )
  expected_p <- c(
    0.79360741127691903,
    0.01123801228478,
    0.0105731631815558
  )

  expect_true(all(long$available))
  expect_identical(long$status, rep("ok", 3L))
  expect_identical(long$df, c(1L, 1L, 2L))
  expect_equal(long$statistic, expected_statistic, tolerance = 1e-10)
  expect_equal(long$p_value, expected_p, tolerance = 1e-10)
  expect_equal(long$q_value, expected_p, tolerance = 1e-10)
  expect_equal(
    long$log_p_value,
    stats::pchisq(
      expected_statistic, c(1L, 1L, 2L),
      lower.tail = FALSE, log.p = TRUE
    ),
    tolerance = 1e-10
  )

  expect_equal(fit$results$p_occupancy, expected_p[[1L]], tolerance = 1e-10)
  expect_equal(fit$results$p_abundance, expected_p[[2L]], tolerance = 1e-10)
  expect_equal(fit$results$p_joint, expected_p[[3L]], tolerance = 1e-10)
  expect_identical(names(fit$results), doram_public_result_columns)

  gamma <- fit$posterior$gamma_null[, "target"]
  expect_true(all(is.finite(gamma)))
  expect_true(all(gamma >= 0 & gamma <= 1))
  expect_identical(
    unname(gamma[fixture$counts[, "target"] > 0]),
    rep(0, sum(fixture$counts[, "target"] > 0))
  )
  expect_equal(
    fit$posterior$gamma_null + fit$posterior$tau_null,
    matrix(
      1, nrow = nrow(fixture$counts), ncol = 1L,
      dimnames = dimnames(fit$posterior$gamma_null)
    ),
    tolerance = 1e-14
  )
  expect_true(all(is.finite(fit$posterior$rho_null)))
})
