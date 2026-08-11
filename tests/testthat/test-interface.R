test_that("metadata IDs are matched as an exact set and reordered safely", {
  x <- doram_toy_data()
  ordered <- doram_fit_toy(x)

  shuffled <- x
  shuffled$metadata <- x$metadata[rev(x$ids), , drop = FALSE]
  reordered <- doram_fit_toy(shuffled)

  expect_identical(reordered$sample_id, x$ids)
  expect_identical(rownames(reordered$posterior$gamma_null), x$ids)
  expect_equal(reordered$results, ordered$results, tolerance = 0)
  expect_equal(reordered$details, ordered$details, tolerance = 0)
  expect_equal(reordered$descriptives, ordered$descriptives, tolerance = 0)
  expect_equal(reordered$diagnostics, ordered$diagnostics, tolerance = 0)
  expect_equal(reordered$posterior, ordered$posterior, tolerance = 0)

  missing_metadata <- x$metadata[-1L, , drop = FALSE]
  expect_error(
    DORAM(x$counts, missing_metadata, "group", "original_library_size"),
    "same sample IDs|identical sample IDs|ID set"
  )

  extra_row <- x$metadata[1L, , drop = FALSE]
  rownames(extra_row) <- "extra_sample"
  extra_metadata <- rbind(x$metadata, extra_row)
  expect_error(
    DORAM(x$counts, extra_metadata, "group", "original_library_size"),
    "same sample IDs|identical sample IDs|ID set"
  )

  duplicate_metadata <- x$metadata
  attr(duplicate_metadata, "row.names") <- c(x$ids[-1L], x$ids[[2L]])
  expect_error(
    DORAM(x$counts, duplicate_metadata, "group", "original_library_size"),
    "unique.*row names|duplicate.*sample"
  )

  duplicate_counts <- x$counts
  rownames(duplicate_counts)[2L] <- rownames(duplicate_counts)[1L]
  expect_error(
    DORAM(duplicate_counts, x$metadata, "group", "original_library_size"),
    "unique.*row names|duplicate.*sample"
  )

  empty_counts <- x$counts
  rownames(empty_counts)[1L] <- ""
  expect_error(
    DORAM(empty_counts, x$metadata, "group", "original_library_size"),
    "nonempty.*row names|sample IDs"
  )

  expect_error(
    DORAM(x$counts, as.matrix(x$metadata), "group", "original_library_size"),
    "metadata must be a data.frame"
  )
})

test_that("metadata selectors and public controls are exact", {
  x <- doram_toy_data()

  expect_error(
    DORAM(x$counts, x$metadata, "missing_group", "original_library_size"),
    "group.*column"
  )
  expect_error(
    DORAM(x$counts, x$metadata, "group", "missing_depth"),
    "library_size.*column"
  )
  expect_error(
    DORAM(
      x$counts, x$metadata, "group", "original_library_size",
      covariates = "missing_covariate"
    ),
    "covariate.*column"
  )
  expect_error(
    DORAM(
      x$counts, x$metadata, c("group", "adjustment"),
      "original_library_size"
    ),
    "group.*one.*column|single.*group"
  )
  expect_error(
    DORAM(
      x$counts, x$metadata, "group",
      c("original_library_size", "adjustment")
    ),
    "library_size.*one.*column|single.*library_size"
  )
  expect_error(
    DORAM(
      x$counts, x$metadata, "group", "original_library_size",
      covariates = character()
    ),
    "covariates.*NULL|at least one"
  )
  expect_error(
    DORAM(
      x$counts, x$metadata, "group", "original_library_size",
      covariates = c("adjustment", "adjustment")
    ),
    "covariate.*unique|duplicate"
  )
  expect_error(doram_fit_toy(x, cores = TRUE), "cores.*positive integer")
  expect_error(doram_fit_toy(x, cores = 1.5), "cores.*positive integer")
  expect_error(
    doram_fit_toy(x, verbose = 1),
    "verbose.*logical|verbose.*TRUE or FALSE"
  )

  expect_error(
    doram_fit_toy(x, sampling_design = "iid_random_design"),
    "unused argument"
  )
  expect_error(doram_fit_toy(x, p_adjust_method = "BH"), "unused argument")
  expect_error(doram_fit_toy(x, alpha = 0.05), "unused argument")
  expect_error(doram_fit_toy(x, keep_full_fits = FALSE), "unused argument")
})

test_that("reference rules never depend silently on factor or character order", {
  x <- doram_toy_data()
  numeric_fit <- doram_fit_toy(x)
  expect_identical(
    numeric_fit$contrast,
    c(reference = "0", comparison = "1")
  )

  logical_x <- x
  logical_x$metadata$group <- as.logical(logical_x$metadata$group)
  logical_fit <- doram_fit_toy(logical_x)
  expect_identical(
    logical_fit$contrast,
    c(reference = "FALSE", comparison = "TRUE")
  )

  character_x <- x
  character_x$metadata$group <- ifelse(
    character_x$metadata$group == 0, "control", "treated"
  )
  expect_error(doram_fit_toy(character_x), "reference.*(required|supplied)")
  character_fit <- doram_fit_toy(character_x, reference = "control")
  expect_identical(
    character_fit$contrast,
    c(reference = "control", comparison = "treated")
  )

  factor_x <- character_x
  factor_x$metadata$group <- factor(
    factor_x$metadata$group, levels = c("treated", "control")
  )
  expect_error(doram_fit_toy(factor_x), "reference.*(required|supplied)")
  factor_fit <- doram_fit_toy(factor_x, reference = "control")
  expect_equal(factor_fit$results, character_fit$results, tolerance = 0)
  expect_equal(factor_fit$details, character_fit$details, tolerance = 0)
  expect_identical(factor_fit$contrast, character_fit$contrast)

  expect_error(
    doram_fit_toy(character_x, reference = "not_observed"),
    "reference.*observed"
  )
})

test_that("reversing reference changes direction, not the two-sided abundance test", {
  n <- 40L
  x <- doram_toy_data(n = n, taxa = "rare_zeros")
  depth <- x$metadata$original_library_size
  group <- x$metadata$group
  y <- pmax(1L, as.integer(round(depth * (0.002 + 0.0002 * group))))
  y[c(2L, n / 2L + 2L)] <- 0L
  x$counts[, "rare_zeros"] <- y

  reference_zero <- doram_fit_toy(x, reference = 0)
  reference_one <- doram_fit_toy(x, reference = 1)
  abundance_zero <- reference_zero$details[
    reference_zero$details$endpoint == "abundance", , drop = FALSE
  ]
  abundance_one <- reference_one$details[
    reference_one$details$endpoint == "abundance", , drop = FALSE
  ]

  expect_true(abundance_zero$available)
  expect_true(abundance_one$available)
  expect_equal(abundance_one$estimate, -abundance_zero$estimate,
               tolerance = 1e-12)
  expect_equal(abundance_one$statistic, abundance_zero$statistic,
               tolerance = 1e-12)
  expect_equal(abundance_one$p_value, abundance_zero$p_value,
               tolerance = 1e-12)
})

test_that("mixed covariates are expanded deterministically after ID alignment", {
  x <- doram_toy_data(n = 24L, taxa = "all_zero")
  index <- seq_along(x$ids)
  x$metadata$batch <- factor(
    c("B", "A", "C")[(index - 1L) %% 3L + 1L],
    levels = c("C", "B", "A", "unused")
  )
  x$metadata$logical_adjustment <- index %% 4L == 0L
  x$metadata$character_adjustment <- c("low", "high")[(index - 1L) %% 2L + 1L]

  fitted <- doram_fit_toy(
    x, covariates = c("adjustment", "batch", "logical_adjustment")
  )
  shuffled <- x
  shuffled$metadata <- x$metadata[rev(x$ids), , drop = FALSE]
  fitted_shuffled <- doram_fit_toy(
    shuffled, covariates = c("adjustment", "batch", "logical_adjustment")
  )
  expect_equal(fitted_shuffled$results, fitted$results, tolerance = 0)
  expect_equal(fitted_shuffled$details, fitted$details, tolerance = 0)
  expect_equal(fitted_shuffled$diagnostics, fitted$diagnostics, tolerance = 0)
  expect_identical(
    fitted_shuffled$settings$covariate_design,
    fitted$settings$covariate_design
  )

  dropped <- x
  dropped$metadata$batch <- droplevels(dropped$metadata$batch)
  fitted_dropped <- doram_fit_toy(
    dropped, covariates = c("adjustment", "batch", "logical_adjustment")
  )
  expect_equal(fitted_dropped$results, fitted$results, tolerance = 0)
  expect_equal(fitted_dropped$details, fitted$details, tolerance = 0)
  expect_identical(
    fitted_dropped$settings$covariate_design,
    fitted$settings$covariate_design
  )

  character_fit <- doram_fit_toy(
    x, covariates = c("adjustment", "character_adjustment")
  )
  factor_character <- x
  factor_character$metadata$character_adjustment <- factor(
    factor_character$metadata$character_adjustment,
    levels = sort(unique(factor_character$metadata$character_adjustment))
  )
  factor_character_fit <- doram_fit_toy(
    factor_character,
    covariates = c("adjustment", "character_adjustment")
  )
  expect_equal(factor_character_fit$results, character_fit$results,
               tolerance = 0)
  expect_equal(factor_character_fit$details, character_fit$details,
               tolerance = 0)
})

test_that("covariate coding is invariant to global contrast options", {
  x <- doram_toy_data(n = 24L, taxa = "all_zero")
  x$metadata$site <- factor(
    c("B", "A", "C")[(seq_along(x$ids) - 1L) %% 3L + 1L],
    levels = c("C", "B", "A")
  )

  treatment_fit <- local({
    old <- options(contrasts = c("contr.treatment", "contr.poly"))
    on.exit(options(old), add = TRUE)
    doram_fit_toy(x, covariates = c("adjustment", "site"))
  })
  sum_fit <- local({
    old <- options(contrasts = c("contr.sum", "contr.helmert"))
    on.exit(options(old), add = TRUE)
    doram_fit_toy(x, covariates = c("adjustment", "site"))
  })

  expect_identical(
    treatment_fit$settings$covariate_design,
    sum_fit$settings$covariate_design
  )
  expect_equal(treatment_fit$results, sum_fit$results, tolerance = 0)
  expect_equal(treatment_fit$details, sum_fit$details, tolerance = 0)
  expect_equal(treatment_fit$diagnostics, sum_fit$diagnostics, tolerance = 0)
})

test_that("invalid or rank-deficient covariate expansions fail before fitting", {
  x <- doram_toy_data(n = 20L, taxa = "all_zero")

  missing <- x
  missing$metadata$adjustment[1L] <- NA_real_
  expect_error(
    doram_fit_toy(missing, covariates = "adjustment"),
    "missing|nonmissing|NA"
  )

  nonfinite <- x
  nonfinite$metadata$adjustment[1L] <- Inf
  expect_error(
    doram_fit_toy(nonfinite, covariates = "adjustment"),
    "finite|invalid numeric design"
  )

  one_level <- x
  one_level$metadata$constant_factor <- factor(rep("only", nrow(one_level$metadata)))
  expect_error(
    doram_fit_toy(one_level, covariates = "constant_factor"),
    "one observed level|at least two|2 or more levels|constant|rank deficient"
  )

  constant_character <- x
  constant_character$metadata$constant_character <- "only"
  expect_error(
    doram_fit_toy(constant_character, covariates = "constant_character"),
    "one observed level|at least two|2 or more levels|constant|rank deficient"
  )

  duplicate_numeric <- x
  duplicate_numeric$metadata$adjustment_copy <- duplicate_numeric$metadata$adjustment
  expect_error(
    doram_fit_toy(
      duplicate_numeric,
      covariates = c("adjustment", "adjustment_copy")
    ),
    "rank deficient|linearly dependent"
  )

  group_copy <- x
  group_copy$metadata$group_copy <- group_copy$metadata$group
  expect_error(
    doram_fit_toy(group_copy, covariates = "group_copy"),
    "group.*explained|group.*covariate"
  )
  expect_error(
    doram_fit_toy(x, covariates = "group"),
    "group.*covariate|tested group"
  )

  duplicate_names <- x
  duplicate_names$metadata$another <- duplicate_names$metadata$adjustment
  names(duplicate_names$metadata)[
    names(duplicate_names$metadata) == "another"
  ] <- "adjustment"
  expect_error(
    doram_fit_toy(duplicate_names, covariates = "adjustment"),
    "metadata.*column names.*unique|metadata.*unique.*column names|duplicate.*column"
  )

  expanded_duplicate <- x
  expanded_duplicate$metadata$site <- factor(
    rep(c("A", "B"), length.out = nrow(expanded_duplicate$metadata))
  )
  expanded_duplicate$metadata$siteB <-
    as.numeric(expanded_duplicate$metadata$site == "B")
  expect_error(
    doram_fit_toy(expanded_duplicate, covariates = c("site", "siteB")),
    paste0(
      "expanded.*unique|duplicate.*design|rank deficient|",
      "invalid numeric design|numeric design.*unique column names"
    )
  )
})

test_that("original library size is required, retained after taxon selection, and bounded", {
  x <- doram_toy_data(n = 20L)
  x$counts[, "taxon_a"] <- 1L
  x$counts[, "taxon_b"] <- 2L

  fit <- doram_fit_toy(x, taxa = "taxon_a")
  descriptive <- fit$descriptives[fit$descriptives$taxon == "taxon_a", ]
  reference <- x$metadata$group == 0
  exposed <- x$metadata$group == 1
  expect_equal(
    descriptive$mean_y_over_n_reference,
    mean(x$counts[reference, "taxon_a"] /
           x$metadata$original_library_size[reference]),
    tolerance = 0
  )
  expect_equal(
    descriptive$mean_y_over_n_comparison,
    mean(x$counts[exposed, "taxon_a"] /
           x$metadata$original_library_size[exposed]),
    tolerance = 0
  )
  expect_identical(fit$taxa, "taxon_a")
  expect_identical(colnames(fit$posterior$gamma_null), "taxon_a")

  expect_error(
    DORAM(x$counts, x$metadata, group = "group"),
    "library_size.*missing|argument.*library_size"
  )
  expect_error(
    DORAM(x$counts, x$metadata, "group", NULL),
    "library_size.*column|required"
  )

  noninteger <- x
  noninteger$metadata$original_library_size[1L] <- 2000.5
  expect_error(doram_fit_toy(noninteger), "positive integer")

  too_large <- x
  too_large$metadata$original_library_size[1L] <- 1e9 + 1
  expect_error(doram_fit_toy(too_large), "supported numerical range")

  count_exceeds <- x
  count_exceeds$counts[1L, "taxon_a"] <-
    count_exceeds$metadata$original_library_size[1L] + 1L
  expect_error(doram_fit_toy(count_exceeds), "count.*library[_ ]size")

  sum_exceeds <- x
  sum_exceeds$metadata$original_library_size[] <- 1000L
  sum_exceeds$counts[,] <- 600L
  expect_error(doram_fit_toy(sum_exceeds), "sum.*library[_ ]size")

  expect_error(doram_fit_toy(x, taxa = "unknown_taxon"), "unknown taxa")
  selected <- doram_fit_toy(x, taxa = c("taxon_b", "taxon_a"))
  expect_identical(selected$taxa, c("taxon_b", "taxon_a"))
  expect_identical(selected$results$taxon, c("taxon_b", "taxon_a"))
  expect_identical(
    colnames(selected$posterior$gamma_null),
    c("taxon_b", "taxon_a")
  )
})

test_that("the publication result and object schemas are stable", {
  x <- doram_toy_data()
  fit <- expect_silent(doram_fit_toy(x))

  expect_s3_class(fit, "DORAM")
  expect_length(fit, length(doram_top_level_components))
  expect_setequal(names(fit), doram_top_level_components)
  expect_identical(names(fit$results), doram_public_result_columns)
  expect_identical(names(fit$details), doram_detail_columns)
  expect_identical(
    names(fit$descriptives),
    c(
      "taxon", "n_samples", "zeros_reference", "zeros_comparison",
      "positives_reference", "positives_comparison",
      "mean_y_over_n_reference", "mean_y_over_n_comparison"
    )
  )
  expect_identical(names(fit$diagnostics), c("fit", "endpoint", "boundary"))
  expect_true(all(c(
    "taxon", "available", "status", "error_message", "stationarity",
    "converged_starts", "nuisance_condition_number",
    "delta_null_constraint"
  ) %in% names(fit$diagnostics$fit)))
  expect_true(all(c(
    "taxon", "endpoint", "max_leverage",
    "score_ess", "condition_ratio", "calibration", "stage"
  ) %in% names(fit$diagnostics$endpoint)))
  expect_identical(names(fit$diagnostics$boundary), doram_boundary_columns)
  expect_identical(nrow(fit$diagnostics$boundary), 0L)
  expect_identical(
    names(fit$posterior),
    c("gamma_null", "rho_null", "tau_null")
  )
  expect_identical(dim(fit$posterior$gamma_null), c(20L, 2L))
  expect_identical(dimnames(fit$posterior$gamma_null),
                   list(x$ids, c("taxon_a", "taxon_b")))
  expect_equal(
    fit$posterior$tau_null,
    1 - fit$posterior$gamma_null,
    tolerance = 0
  )
  expect_identical(nrow(fit$results), 2L)
  expect_identical(nrow(fit$details), 6L)
  expect_identical(
    fit$details$endpoint,
    rep(c("occupancy", "abundance", "joint"), 2L)
  )
  expect_identical(as.data.frame(fit), fit$results)

  expect_identical(fit$settings$p_adjust_method, "BH")
  expect_identical(fit$settings$family_size, 2L)
})

test_that("unavailable taxa are retained and fail closed without noise", {
  x <- doram_toy_data()
  fit <- expect_silent(doram_fit_toy(x))

  expect_false(any(fit$details$available))
  expect_true(all(is.na(fit$details$statistic)))
  expect_true(all(is.na(fit$details$p_value)))
  expect_true(all(is.na(fit$details$q_value)))
  expect_true(all(nzchar(fit$details$status)))
  expect_true(all(is.na(fit$results$p_occupancy)))
  expect_true(all(is.na(fit$results$p_abundance)))
  expect_true(all(is.na(fit$results$p_joint)))
  expect_true(all(is.na(fit$posterior$gamma_null)))
  expect_true(all(is.na(fit$posterior$rho_null)))

  expect_output(print(fit), "DORAM")
})

test_that("abundance inference survives an occupancy refusal and joint fails closed", {
  n <- 40L
  x <- doram_toy_data(n = n, taxa = "rare_zeros")
  depth <- x$metadata$original_library_size
  group <- x$metadata$group
  y <- pmax(1L, as.integer(round(depth * (0.002 + 0.0002 * group))))
  y[c(2L, n / 2L + 2L)] <- 0L
  x$counts[, "rare_zeros"] <- y

  fit <- expect_silent(doram_fit_toy(x))
  occupancy <- fit$details[fit$details$endpoint == "occupancy", ]
  abundance <- fit$details[fit$details$endpoint == "abundance", ]
  joint <- fit$details[fit$details$endpoint == "joint", ]

  expect_false(occupancy$available)
  expect_true(abundance$available)
  expect_false(joint$available)
  expect_identical(abundance$df, 1L)
  expect_match(joint$status, "^occupancy_component_")
})

test_that("active optimizer boundaries are reported with interpretable units", {
  x <- doram_toy_data(n = 20L, taxa = "target")
  count_data <- DORAM:::cn_prepare_data(
    y = rep(c(0L, 1L), length.out = 20L),
    N = x$metadata$original_library_size,
    g = x$metadata$group,
    z = NULL,
    subject_id = x$ids
  )
  layout <- DORAM:::cn_layout(count_data)
  bounds <- DORAM:::cn_opt_bounds(layout)
  values <- (bounds$lower + bounds$upper) / 2
  parameter <- "log_sigma_g0"
  expect_true(parameter %in% names(values))
  values[[parameter]] <- bounds$lower[[parameter]]
  fit <- list(
    active_parameters = parameter,
    selected_index = 1L,
    candidates = list(list(opt_par = values)),
    stage = "restricted_fit"
  )

  observed <- DORAM:::.doram_boundary_rows("target", fit, count_data)

  expect_identical(names(observed), doram_boundary_columns)
  expect_identical(observed$taxon, "target")
  expect_identical(observed$endpoint, "occupancy")
  expect_identical(observed$parameter_block, "present_scale")
  expect_identical(observed$side, "lower")
  expect_equal(observed$value, bounds$lower[[parameter]], tolerance = 0)
  expect_equal(observed$bound, bounds$lower[[parameter]], tolerance = 0)
  expect_equal(observed$distance, 0, tolerance = 0)
  expect_equal(observed$transformed_value, 0.1, tolerance = 1e-14)
  expect_true(observed$distance <= observed$active_limit)
})

test_that("BH uses the full selected-taxon denominator including unavailable rows", {
  table <- data.frame(
    taxon = c("a", "b", "c"),
    endpoint = rep("joint", 3L),
    available = c(TRUE, TRUE, FALSE),
    estimate = NA_real_, statistic = c(1, 2, NA),
    df = c(2L, 2L, NA_integer_),
    p_value = c(0.01, 0.20, NA),
    log_p_value = c(log(0.01), log(0.20), NA),
    q_value = NA_real_, status = c("ok", "ok", "unavailable"),
    stringsAsFactors = FALSE
  )

  adjusted <- DORAM:::.doram_adjust_results(table)

  expect_equal(
    adjusted$q_value[1:2],
    stats::p.adjust(c(0.01, 0.20), method = "BH", n = 3L)
  )
  expect_true(is.na(adjusted$q_value[3L]))
})

test_that("default fitting is quiet, deterministic, and core invariant", {
  x <- doram_toy_data(n = 20L, taxa = c("a", "b", "c"))
  set.seed(51943)
  seed_before <- .Random.seed
  serial <- expect_silent(doram_fit_toy(x, cores = 1L))
  expect_identical(.Random.seed, seed_before)

  expect_message(
    doram_fit_toy(x, taxa = "a", verbose = TRUE),
    "fitting taxon"
  )

  if (!identical(.Platform$OS.type, "windows")) {
    parallel <- expect_silent(doram_fit_toy(x, cores = 2L))
    expect_equal(parallel$results, serial$results, tolerance = 0)
    expect_equal(parallel$details, serial$details, tolerance = 0)
    expect_equal(parallel$diagnostics, serial$diagnostics, tolerance = 0)
    expect_equal(parallel$posterior, serial$posterior, tolerance = 0)
    expect_identical(.Random.seed, seed_before)
  }
})
