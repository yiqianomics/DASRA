.doram_example_source_root <- function(required_file) {
  test_root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
  candidates <- unique(c(
    test_root,
    file.path(test_root, "00_pkg_src", "DORAM")
  ))
  present <- vapply(candidates, function(path) {
    file.exists(file.path(path, required_file))
  }, logical(1L))
  if (!any(present)) return(NULL)
  candidates[which(present)[[1L]]]
}

.doram_with_restored_rng <- function(code) {
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
  force(code)
}

.doram_roxygen_example <- function(source_file) {
  source <- readLines(source_file, warn = FALSE)
  start <- grep("^#' @examples[[:space:]]*$", source)
  stop <- grep("^#' @importFrom[[:space:]]", source)
  stop <- stop[stop > start]
  if (length(start) != 1L || !length(stop)) {
    stop("could not identify the DORAM roxygen example block", call. = FALSE)
  }
  block <- source[seq.int(start + 1L, stop[[1L]] - 1L)]
  sub("^#' ?", "", block)
}

.doram_expect_example_fit <- function(fit, count_table, sample_data) {
  taxa <- c("joint_signal", "abundance_signal", "null_taxon")
  endpoints <- c("occupancy", "abundance", "joint")

  expect_s3_class(fit, "DORAM")
  expect_identical(fit$taxa, taxa)
  expect_identical(fit$results$taxon, taxa)
  expect_identical(dim(count_table), c(180L, 3L))
  expect_identical(colnames(count_table), taxa)
  expect_identical(rownames(count_table), fit$sample_id)
  expect_identical(rownames(sample_data), fit$sample_id)
  expect_identical(nrow(fit$results), 3L)
  expect_identical(nrow(fit$details), 9L)

  expected_keys <- expand.grid(
    endpoint = endpoints,
    taxon = taxa,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  observed_keys <- fit$details[c("taxon", "endpoint")]
  expect_setequal(
    paste(observed_keys$taxon, observed_keys$endpoint, sep = "::"),
    paste(expected_keys$taxon, expected_keys$endpoint, sep = "::")
  )
  expect_true(all(fit$details$available))
  expect_identical(fit$details$status, rep("ok", 9L))
  expect_identical(
    fit$details$df,
    rep(c(1L, 1L, 2L), times = 3L)
  )
  for (column in c(
    "estimate", "statistic", "p_value", "log_p_value", "q_value"
  )) {
    if (identical(column, "estimate")) {
      values <- fit$details[[column]][fit$details$endpoint == "abundance"]
    } else {
      values <- fit$details[[column]]
    }
    expect_true(all(is.finite(values)), info = column)
  }
  expect_true(all(fit$details$p_value >= 0 & fit$details$p_value <= 1))
  expect_true(all(fit$details$q_value >= 0 & fit$details$q_value <= 1))
  expect_identical(names(fit$results), doram_public_result_columns)

  expect_true(all(fit$diagnostics$fit$available))
  expect_identical(fit$diagnostics$fit$status, rep("ok", 3L))
  expect_true(all(is.finite(fit$diagnostics$fit$stationarity)))
  expect_true(all(
    fit$diagnostics$fit$stationarity <=
      DORAM:::CN_D3_CONTRACT$optimizer$stationarity_relative_score_norm
  ))
  expect_true(all(fit$diagnostics$fit$converged_starts >= 1L))
  expect_true(all(is.finite(
    fit$diagnostics$fit$nuisance_condition_number
  )))
  expect_identical(
    fit$diagnostics$fit$delta_null_constraint,
    rep(0, 3L)
  )
  expect_true(all(fit$diagnostics$endpoint$available))
  expect_identical(fit$diagnostics$endpoint$status, rep("ok", 9L))

  expect_true(all(is.finite(fit$posterior$gamma_null)))
  expect_true(all(is.finite(fit$posterior$rho_null)))
  expect_true(all(is.finite(fit$posterior$tau_null)))
  expect_true(all(fit$posterior$gamma_null >= 0 &
                    fit$posterior$gamma_null <= 1))
  expect_true(all(fit$posterior$rho_null >= 0 &
                    fit$posterior$rho_null <= 1))
  expect_equal(
    fit$posterior$gamma_null + fit$posterior$tau_null,
    matrix(
      1, nrow = 180L, ncol = 3L,
      dimnames = dimnames(fit$posterior$gamma_null)
    ),
    tolerance = 1e-14
  )

  joint_signal <- fit$results[fit$results$taxon == "joint_signal", ]
  expect_lt(joint_signal$p_occupancy, 1e-8)
  expect_lt(joint_signal$p_abundance, 1e-3)
  expect_lt(joint_signal$p_joint, 1e-8)

  abundance_signal <- fit$results[
    fit$results$taxon == "abundance_signal",
  ]
  expect_lt(abundance_signal$p_abundance, 0.01)

  null_taxon <- fit$results[fit$results$taxon == "null_taxon", ]
  expect_true(all(unlist(null_taxon[c(
    "p_occupancy", "p_abundance", "p_joint"
  )]) > 0.05))
}

test_that("the published example has a coherent three-taxon dual-endpoint DGP", {
  source_root <- .doram_example_source_root(file.path("R", "DORAM.R"))
  skip_if(is.null(source_root), "raw DORAM roxygen source is unavailable")
  block <- .doram_roxygen_example(file.path(source_root, "R", "DORAM.R"))
  dontrun <- which(block == "\\dontrun{")
  expect_length(dontrun, 1L)

  data_code <- paste(block[seq_len(dontrun - 1L)], collapse = "\n")
  example <- new.env(parent = .GlobalEnv)
  .doram_with_restored_rng(eval(parse(text = data_code), envir = example))

  taxa <- c("joint_signal", "abundance_signal", "null_taxon")
  expect_identical(example$n, 180L)
  expect_identical(dim(example$count_table), c(180L, 3L))
  expect_identical(colnames(example$count_table), taxa)
  expect_identical(rownames(example$count_table), example$sample_id)
  expect_true(is.integer(example$count_table))
  expect_true(all(example$count_table >= 0L))
  expect_true(all(rowSums(example$count_table) <= example$reads))
  expect_lt(max(rowSums(example$taxon_probability)), 1)
  expect_identical(rownames(example$sample_data), example$sample_id)
  expect_identical(levels(example$sample_data$diagnosis), c("control", "case"))
  expect_identical(example$delta, c(2, 0, 0))
  expect_identical(example$zeta, c(0, 0.65, 0))

  group_index <- split(seq_len(example$n), example$group)
  for (j in seq_len(3L)) {
    expect_gt(sum(example$count_table[group_index[["0"]], j] == 0L), 0L)
    expect_gt(sum(example$count_table[group_index[["1"]], j] == 0L), 0L)
    expect_gt(sum(example$count_table[group_index[["0"]], j] > 0L), 0L)
    expect_gt(sum(example$count_table[group_index[["1"]], j] > 0L), 0L)
  }

  joint_index <- match("joint_signal", taxa)
  rho_reference <- plogis(
    example$alpha0[joint_index] +
      example$alpha_z[joint_index] * example$age_z
  )
  rho_comparison <- plogis(
    example$alpha0[joint_index] +
      example$alpha_z[joint_index] * example$age_z +
      example$delta[joint_index]
  )
  present_probability <- plogis(example$latent_abundance[, joint_index])
  expect_true(all(rho_comparison > rho_reference))
  expect_true(all(
    (1 - rho_comparison) * present_probability <
      (1 - rho_reference) * present_probability
  ))

  example_text <- paste(block, collapse = "\n")
  expect_match(example_text, "set.seed(20260810)", fixed = TRUE)
  expect_match(example_text, "counts = count_table", fixed = TRUE)
  expect_match(example_text, "metadata = sample_data", fixed = TRUE)
  expect_match(example_text, 'reference = "control"', fixed = TRUE)
  expect_match(example_text, 'covariates = "age_z"', fixed = TRUE)
  expect_match(example_text, "cores = 3L", fixed = TRUE)
})

test_that("the actual generated help example remains a successful three-taxon fit", {
  skip_if_not(
    identical(tolower(Sys.getenv("DORAM_RUN_EXAMPLE_TESTS")), "true"),
    paste(
      "set DORAM_RUN_EXAMPLE_TESTS=true to execute the generated",
      "production example"
    )
  )

  source_root <- .doram_example_source_root(file.path("man", "DORAM.Rd"))
  skip_if(
    is.null(source_root),
    "run devtools::document() before the generated-example regression"
  )
  rd_file <- file.path(source_root, "man", "DORAM.Rd")
  example_file <- tempfile("DORAM-example-", fileext = ".R")
  on.exit(unlink(example_file), add = TRUE)
  tools::Rd2ex(rd_file, example_file, commentDontrun = FALSE)

  example <- new.env(parent = .GlobalEnv)
  run_example <- function() {
    invisible(capture.output(
      .doram_with_restored_rng(
        sys.source(example_file, envir = example, keep.source = FALSE)
      )
    ))
  }
  if (identical(.Platform$OS.type, "windows")) {
    expect_warning(run_example(), "parallel taxon fitting is unavailable")
  } else {
    expect_warning(run_example(), NA)
  }

  expect_true(exists("fit", envir = example, inherits = FALSE))
  .doram_expect_example_fit(
    example$fit,
    example$count_table,
    example$sample_data
  )
})
