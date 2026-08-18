#!/usr/bin/env Rscript

# =============================================================================
# DASRA structural-absence FDR stress simulation
#
# Purpose
# -------
# Evaluate whether the DASRA structural-absence component requires any
# cross-taxon "compositional-bias correction." The simulation keeps the true
# structural-absence null for all or most taxa while imposing sparse, dense,
# or dispersion-only changes in latent positive abundances. Counts are then
# generated jointly by multinomial sequencing, so closure and cross-taxon
# dependence are present by construction.
#
# Primary quantities
# ------------------
# 1. Marginal structural-null rejection rate.
# 2. Empirical FDR after BH and BY adjustment.
# 3. Structural-signal power.
# 4. Taxon retention and structural-component formation rates.
# 5. Rejection rates for directly abundance-shifted structural-null taxa and
#    for closure-only structural-null taxa.
#
# The script never adds a pseudocount. Unretained taxa are excluded from the
# multiplicity family, matching DASRA. Retained but unformed structural tests
# enter the family with p = 1, matching the package implementation.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("DASRA", "ggplot2", "dplyr", "tidyr")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install the following packages before running this script: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(DASRA)
    library(ggplot2)
    library(dplyr)
})

# Prevent nested BLAS/OpenMP parallelism when outer simulation workers are used.
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

env_flag <- function(name, default = FALSE) {
    value <- Sys.getenv(name, unset = if (default) "true" else "false")
    tolower(trimws(value)) %in% c("1", "true", "t", "yes", "y")
}

env_integer <- function(name, default, minimum = 1L) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || !is.finite(value) || value < minimum) {
        stop(sprintf("%s must be an integer at least %d.", name, minimum),
             call. = FALSE)
    }
    value
}

env_numeric <- function(name, default, lower = -Inf, upper = Inf) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || !is.finite(value) ||
        value < lower || value > upper) {
        stop(
            sprintf("%s must be finite and lie in [%s, %s].",
                    name, format(lower), format(upper)),
            call. = FALSE
        )
    }
    value
}

requested_cores <- env_integer("DASRA_SA_FDR_CORES", 8L)
physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(physical_cores) || physical_cores < 1L) physical_cores <- 1L

CONFIG <- list(
    script_version = "2026-08-17-sa-fdr-v1",
    n_rep = env_integer("DASRA_SA_FDR_REPS", 200L),
    n_per_group = env_integer("DASRA_SA_FDR_N_PER_GROUP", 120L),
    n_taxa = env_integer("DASRA_SA_FDR_TAXA", 50L, minimum = 20L),
    n_cores = max(1L, min(8L, requested_cores, physical_cores)),
    alpha = env_numeric("DASRA_SA_FDR_ALPHA", 0.05, 0, 1),
    base_seed = env_integer("DASRA_SA_FDR_SEED", 20260817L),
    output_dir = Sys.getenv(
        "DASRA_SA_FDR_OUTDIR",
        file.path(getwd(), "dasra_structural_absence_fdr_simulation")
    ),
    resume = env_flag("DASRA_SA_FDR_RESUME", TRUE),
    include_heavy_tail = env_flag("DASRA_SA_FDR_HEAVY_TAIL", FALSE),

    # Sequencing-depth distribution
    depth_median = 8000,
    depth_sdlog = 0.45,
    depth_min = 1500L,
    depth_max = 40000L,

    # Structural-absence truth
    structural_signal_fraction = 0.10,
    structural_log_odds_effect = 1.20,

    # Positive-abundance stress settings
    sparse_abundance_fraction = 0.10,
    sparse_abundance_log_effect = 1.00,
    dense_abundance_fraction = 0.40,
    dense_abundance_log_effect = 0.60,
    dispersion_fraction = 0.25,
    dispersion_multiplier = 1.80,

    # Latent abundance and dependence
    abundance_sd_min = 0.45,
    abundance_sd_max = 0.75,
    sample_factor_sd = 0.65,
    other_mass_log_mean = log(5),
    other_mass_log_sd = 0.35
)

validate_config <- function(config) {
    stopifnot(
        config$n_rep >= 1L,
        config$n_per_group >= 20L,
        config$n_taxa >= 20L,
        config$n_cores >= 1L,
        config$n_cores <= 8L,
        config$alpha > 0,
        config$alpha < 1,
        config$structural_signal_fraction > 0,
        config$structural_signal_fraction < 0.5,
        config$sparse_abundance_fraction > 0,
        config$sparse_abundance_fraction < 0.5,
        config$dense_abundance_fraction >
            config$sparse_abundance_fraction,
        config$dense_abundance_fraction < 0.8,
        config$dispersion_fraction > 0,
        config$dispersion_fraction < 0.8,
        config$dispersion_multiplier > 1
    )
    invisible(TRUE)
}
validate_config(CONFIG)

SCENARIOS <- data.frame(
    scenario = c(
        "global_null",
        "sparse_abundance_only",
        "dense_abundance_only",
        "dispersion_only",
        "structural_only",
        "mixed_structural_abundance"
    ),
    structural_fraction = c(
        0,
        0,
        0,
        0,
        CONFIG$structural_signal_fraction,
        CONFIG$structural_signal_fraction
    ),
    structural_effect = c(
        0,
        0,
        0,
        0,
        CONFIG$structural_log_odds_effect,
        CONFIG$structural_log_odds_effect
    ),
    abundance_fraction = c(
        0,
        CONFIG$sparse_abundance_fraction,
        CONFIG$dense_abundance_fraction,
        0,
        0,
        CONFIG$sparse_abundance_fraction
    ),
    abundance_effect = c(
        0,
        CONFIG$sparse_abundance_log_effect,
        CONFIG$dense_abundance_log_effect,
        0,
        0,
        CONFIG$sparse_abundance_log_effect
    ),
    dispersion_fraction = c(
        0,
        0,
        0,
        CONFIG$dispersion_fraction,
        0,
        0
    ),
    dispersion_multiplier = c(
        1,
        1,
        1,
        CONFIG$dispersion_multiplier,
        1,
        1
    ),
    stringsAsFactors = FALSE
)

DGP_VARIANTS <- if (isTRUE(CONFIG$include_heavy_tail)) {
    c("lognormal", "heavy_tail")
} else {
    "lognormal"
}

dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(CONFIG$output_dir, "replicate_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

safe_namespace_body <- function(name) {
    namespace <- asNamespace("DASRA")
    if (!exists(name, envir = namespace, inherits = FALSE)) return("")
    object <- get(name, envir = namespace, inherits = FALSE)
    if (!is.function(object)) return("")
    paste(deparse(body(object)), collapse = "\n")
}

text_md5 <- function(text) {
    path <- tempfile(fileext = ".txt")
    on.exit(unlink(path), add = TRUE)
    writeLines(text, path, useBytes = TRUE)
    unname(tools::md5sum(path))
}

DASRA_FINGERPRINT <- text_md5(paste(
    as.character(utils::packageVersion("DASRA")),
    safe_namespace_body("dasra"),
    safe_namespace_body("zt_count_structural_test"),
    sep = "\n---\n"
))

design_signature <- CONFIG
design_signature$n_rep <- NULL
design_signature$n_cores <- NULL
design_signature$resume <- NULL
design_signature$output_dir <- NULL
design_signature$dasra_version <- as.character(utils::packageVersion("DASRA"))
design_signature$dasra_fingerprint <- DASRA_FINGERPRINT
design_signature$scenarios <- SCENARIOS
design_signature$dgp_variants <- DGP_VARIANTS

signature_path <- file.path(CONFIG$output_dir, "simulation_design_config.rds")
if (file.exists(signature_path)) {
    previous_signature <- readRDS(signature_path)
    if (!identical(previous_signature, design_signature)) {
        stop(
            "The output directory contains checkpoints from a different ",
            "simulation design or DASRA implementation. Use a new ",
            "DASRA_SA_FDR_OUTDIR or remove the existing directory.",
            call. = FALSE
        )
    }
} else {
    saveRDS(design_signature, signature_path)
}

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

matched_group_vector <- function(values) {
    c(values, sample(values, length(values), replace = FALSE))
}

make_base_replication <- function(rep_id, dgp, config) {
    set.seed(config$base_seed + 10000L * rep_id +
                 if (identical(dgp, "heavy_tail")) 700000L else 0L)

    n0 <- config$n_per_group
    n <- 2L * n0
    p <- config$n_taxa

    group <- c(rep(0, n0), rep(1, n0))

    z_reference <- rnorm(n0)
    z <- matched_group_vector(z_reference)

    depth_reference <- round(exp(
        rnorm(n0, log(config$depth_median), config$depth_sdlog)
    ))
    depth_reference <- as.integer(clamp(
        depth_reference, config$depth_min, config$depth_max
    ))
    depth <- matched_group_vector(depth_reference)

    # Taxon-specific baseline structural-absence probabilities.
    rho_grid <- seq(qlogis(0.10), qlogis(0.42), length.out = p)
    rho_intercept <- rho_grid + rnorm(p, 0, 0.12)
    rho_z_effect <- rnorm(p, 0, 0.12)

    # Baseline absolute abundance spans roughly two orders of magnitude.
    base_weight <- exp(
        seq(log(0.25), log(0.0025), length.out = p) +
            rnorm(p, 0, 0.12)
    )
    abundance_intercept <- log(base_weight)
    abundance_z_effect <- rnorm(p, 0, 0.10)
    factor_loading <- rnorm(p, 0, 0.25)
    abundance_sd <- runif(
        p, config$abundance_sd_min, config$abundance_sd_max
    )

    sample_factor <- rnorm(n, 0, config$sample_factor_sd)
    presence_uniform <- matrix(runif(n * p), nrow = n, ncol = p)

    if (identical(dgp, "heavy_tail")) {
        # Standardized t_4 errors have variance one.
        error_matrix <- matrix(rt(n * p, df = 4) / sqrt(2),
                               nrow = n, ncol = p)
    } else {
        error_matrix <- matrix(rnorm(n * p), nrow = n, ncol = p)
    }

    other_mass <- exp(
        config$other_mass_log_mean +
            0.10 * z +
            0.15 * sample_factor +
            rnorm(n, 0, config$other_mass_log_sd)
    )

    list(
        group = group,
        z = z,
        depth = depth,
        rho_intercept = rho_intercept,
        rho_z_effect = rho_z_effect,
        abundance_intercept = abundance_intercept,
        abundance_z_effect = abundance_z_effect,
        factor_loading = factor_loading,
        abundance_sd = abundance_sd,
        sample_factor = sample_factor,
        presence_uniform = presence_uniform,
        error_matrix = error_matrix,
        other_mass = other_mass,
        signal_order = sample.int(p),
        taxa = sprintf("Taxon_%03d", seq_len(p)),
        samples = sprintf("Sample_%03d", seq_len(n))
    )
}

select_signal_sets <- function(base, scenario_row, config) {
    p <- config$n_taxa
    order <- base$signal_order

    n_structural <- if (scenario_row$structural_fraction > 0) {
        max(1L, as.integer(round(
            scenario_row$structural_fraction * p
        )))
    } else {
        0L
    }

    n_abundance <- if (scenario_row$abundance_fraction > 0) {
        max(1L, as.integer(round(
            scenario_row$abundance_fraction * p
        )))
    } else {
        0L
    }

    n_dispersion <- if (scenario_row$dispersion_fraction > 0) {
        max(1L, as.integer(round(
            scenario_row$dispersion_fraction * p
        )))
    } else {
        0L
    }

    cursor <- 1L
    structural_set <- integer()
    abundance_set <- integer()
    dispersion_set <- integer()

    if (n_structural > 0L) {
        structural_set <- order[cursor:(cursor + n_structural - 1L)]
        cursor <- cursor + n_structural
    }

    if (n_abundance > 0L) {
        # In the mixed scenario, keep the two signal sets disjoint.
        if (cursor + n_abundance - 1L <= p) {
            abundance_set <- order[cursor:(cursor + n_abundance - 1L)]
            cursor <- cursor + n_abundance
        } else {
            abundance_set <- setdiff(order, structural_set)[
                seq_len(n_abundance)
            ]
        }
    }

    if (n_dispersion > 0L) {
        available <- setdiff(order, union(structural_set, abundance_set))
        if (length(available) < n_dispersion) available <- order
        dispersion_set <- available[seq_len(n_dispersion)]
    }

    list(
        structural = structural_set,
        abundance = abundance_set,
        dispersion = dispersion_set
    )
}

simulate_scenario <- function(base, scenario_row, rep_id, dgp, config) {
    n <- length(base$group)
    p <- config$n_taxa
    group <- base$group

    signal_sets <- select_signal_sets(base, scenario_row, config)

    structural_effect <- numeric(p)
    if (length(signal_sets$structural)) {
        directions <- rep(c(1, -1), length.out = length(
            signal_sets$structural
        ))
        structural_effect[signal_sets$structural] <-
            scenario_row$structural_effect * directions
    }

    abundance_effect <- numeric(p)
    if (length(signal_sets$abundance)) {
        abundance_effect[signal_sets$abundance] <-
            scenario_row$abundance_effect
    }

    dispersion_multiplier <- rep(1, p)
    if (length(signal_sets$dispersion)) {
        dispersion_multiplier[signal_sets$dispersion] <-
            scenario_row$dispersion_multiplier
    }

    lp_rho <- matrix(
        base$rho_intercept,
        nrow = n,
        ncol = p,
        byrow = TRUE
    ) +
        outer(base$z, base$rho_z_effect) +
        outer(group, structural_effect)

    structural_absence_probability <- plogis(lp_rho)
    structurally_absent <-
        base$presence_uniform < structural_absence_probability

    log_abundance_mean <- matrix(
        base$abundance_intercept,
        nrow = n,
        ncol = p,
        byrow = TRUE
    ) +
        outer(base$z, base$abundance_z_effect) +
        outer(base$sample_factor, base$factor_loading) +
        outer(group, abundance_effect)

    sd_matrix <- matrix(
        base$abundance_sd,
        nrow = n,
        ncol = p,
        byrow = TRUE
    )
    case_rows <- which(group == 1)
    if (length(signal_sets$dispersion)) {
        sd_matrix[case_rows, signal_sets$dispersion] <-
            sweep(
                sd_matrix[case_rows, signal_sets$dispersion, drop = FALSE],
                2L,
                dispersion_multiplier[signal_sets$dispersion],
                "*"
            )
    }

    absolute_abundance <- exp(
        log_abundance_mean + sd_matrix * base$error_matrix
    )
    absolute_abundance[structurally_absent] <- 0

    total_mass <- rowSums(absolute_abundance) + base$other_mass
    probability <- absolute_abundance / total_mass
    other_probability <- base$other_mass / total_mass

    # Joint multinomial sequencing preserves the fixed-sum constraint and
    # induces cross-taxon count dependence. The final "other" category is not
    # passed to DASRA but remains part of the library size.
    set.seed(
        config$base_seed +
            1000000L +
            10000L * rep_id +
            100L * match(scenario_row$scenario, SCENARIOS$scenario) +
            if (identical(dgp, "heavy_tail")) 50L else 0L
    )
    count_by_sample <- t(vapply(seq_len(n), function(i) {
        draw <- rmultinom(
            1L,
            size = base$depth[i],
            prob = c(probability[i, ], other_probability[i])
        )
        as.numeric(draw[seq_len(p), 1L])
    }, numeric(p)))

    counts <- t(count_by_sample)
    rownames(counts) <- base$taxa
    colnames(counts) <- base$samples

    metadata <- data.frame(
        group = factor(
            ifelse(group == 0, "control", "case"),
            levels = c("control", "case")
        ),
        z = base$z,
        reads = as.integer(base$depth),
        row.names = base$samples,
        check.names = FALSE
    )

    list(
        counts = counts,
        metadata = metadata,
        structural_alt = seq_len(p) %in% signal_sets$structural,
        abundance_alt = seq_len(p) %in% signal_sets$abundance,
        dispersion_alt = seq_len(p) %in% signal_sets$dispersion,
        structural_effect = structural_effect,
        abundance_effect = abundance_effect,
        baseline_structural_probability = plogis(base$rho_intercept),
        mean_observed_prevalence_reference =
            colMeans(count_by_sample[group == 0, , drop = FALSE] > 0),
        mean_observed_prevalence_comparison =
            colMeans(count_by_sample[group == 1, , drop = FALSE] > 0)
    )
}

fit_dasra_structural <- function(simulated) {
    fit_arguments <- list(
        counts = simulated$counts,
        metadata = simulated$metadata,
        formula = ~ group + z,
        group = "group",
        library_size = "reads",
        taxa_are_rows = TRUE,
        reference = "control",
        p_adjust_method = "BH",
        full_output = FALSE
    )

    supported_arguments <- names(formals(DASRA::dasra))
    if ("component" %in% supported_arguments) {
        fit_arguments$component <- "structural_absence"
    }

    do.call(DASRA::dasra, fit_arguments)
}

extract_structural_results <- function(fit, taxa) {
    p <- length(taxa)
    results <- fit$results
    diagnostics <- fit$diagnostics

    result_taxa <- if ("taxon" %in% names(results)) {
        as.character(results$taxon)
    } else {
        rownames(results)
    }
    result_index <- match(taxa, result_taxa)

    raw_p <- q_bh <- rep(NA_real_, p)
    retained <- rep(FALSE, p)

    if ("p_structural_absence" %in% names(results)) {
        raw_p <- as.numeric(results$p_structural_absence[result_index])
    }
    if ("q_structural_absence" %in% names(results)) {
        q_bh <- as.numeric(results$q_structural_absence[result_index])
    }
    if ("retained" %in% names(results)) {
        retained <- as.logical(results$retained[result_index])
    }

    diagnostic_taxa <- if ("taxon" %in% names(diagnostics)) {
        as.character(diagnostics$taxon)
    } else {
        rownames(diagnostics)
    }
    diagnostic_index <- match(taxa, diagnostic_taxa)

    formed <- rep(FALSE, p)
    reason <- rep("not_retained", p)

    if ("retained" %in% names(diagnostics)) {
        retained <- as.logical(diagnostics$retained[diagnostic_index])
    }
    if ("formed_structural_absence" %in% names(diagnostics)) {
        formed <- as.logical(
            diagnostics$formed_structural_absence[diagnostic_index]
        )
    } else {
        formed <- retained & is.finite(raw_p)
    }
    if ("reason_structural_absence" %in% names(diagnostics)) {
        reason <- as.character(
            diagnostics$reason_structural_absence[diagnostic_index]
        )
    }

    retained[is.na(retained)] <- FALSE
    formed[is.na(formed)] <- FALSE
    reason[is.na(reason) | !nzchar(reason)] <- ifelse(
        retained[is.na(reason) | !nzchar(reason)],
        "unknown_unformed_reason",
        "not_retained"
    )

    # Retained but unformed tests are p = 1 in the implemented family.
    raw_p[retained & !formed] <- 1
    raw_p[retained & formed & !is.finite(raw_p)] <- 1

    q_by <- rep(NA_real_, p)
    family <- retained & is.finite(raw_p)
    if (any(family)) {
        q_by[family] <- p.adjust(raw_p[family], method = "BY")
    }

    list(
        raw_p = raw_p,
        q_bh = q_bh,
        q_by = q_by,
        retained = retained,
        formed = formed,
        reason = reason
    )
}

make_failed_result <- function(taxa, message) {
    p <- length(taxa)
    list(
        raw_p = rep(NA_real_, p),
        q_bh = rep(NA_real_, p),
        q_by = rep(NA_real_, p),
        retained = rep(FALSE, p),
        formed = rep(FALSE, p),
        reason = rep(paste0("dasra_fit_error: ", message), p)
    )
}

safe_rate <- function(x) {
    if (!length(x)) return(NA_real_)
    mean(x)
}

make_metric_rows <- function(
    extracted,
    simulated,
    rep_id,
    dgp,
    scenario,
    config
) {
    structural_alt <- simulated$structural_alt
    structural_null <- !structural_alt
    p <- length(structural_alt)

    adjustments <- list(BH = extracted$q_bh, BY = extracted$q_by)

    do.call(rbind, lapply(names(adjustments), function(adjustment) {
        q_value <- adjustments[[adjustment]]
        rejected <- rep(FALSE, p)
        valid <- extracted$retained & is.finite(q_value)
        rejected[valid] <- q_value[valid] <= config$alpha

        n_rejected <- sum(rejected)
        n_false <- sum(rejected & structural_null)
        n_true <- sum(rejected & structural_alt)
        n_alt <- sum(structural_alt)
        fdp <- if (n_rejected > 0L) n_false / n_rejected else 0

        available_null <- structural_null & extracted$formed &
            is.finite(extracted$raw_p)
        raw_null_rejection_rate <- if (any(available_null)) {
            mean(extracted$raw_p[available_null] <= config$alpha)
        } else {
            NA_real_
        }

        data.frame(
            rep = rep_id,
            dgp = dgp,
            scenario = scenario,
            adjustment = adjustment,
            n_taxa = p,
            n_retained = sum(extracted$retained),
            n_formed = sum(extracted$formed),
            n_rejected = n_rejected,
            n_false = n_false,
            n_true = n_true,
            n_structural_alt = n_alt,
            fdp = fdp,
            any_false = as.integer(n_false > 0L),
            power = if (n_alt > 0L) n_true / n_alt else NA_real_,
            retention_rate = mean(extracted$retained),
            formation_rate_among_retained = if (any(extracted$retained)) {
                mean(extracted$formed[extracted$retained])
            } else {
                NA_real_
            },
            overall_formation_rate = mean(extracted$formed),
            raw_null_rejection_rate = raw_null_rejection_rate,
            stringsAsFactors = FALSE
        )
    }))
}

make_taxon_rows <- function(
    extracted,
    simulated,
    rep_id,
    dgp,
    scenario,
    config
) {
    taxa <- rownames(simulated$counts)
    data.frame(
        rep = rep_id,
        dgp = dgp,
        scenario = scenario,
        taxon = taxa,
        structural_alt = simulated$structural_alt,
        abundance_alt = simulated$abundance_alt,
        dispersion_alt = simulated$dispersion_alt,
        structural_effect = simulated$structural_effect,
        abundance_effect = simulated$abundance_effect,
        baseline_structural_probability =
            simulated$baseline_structural_probability,
        observed_prevalence_reference =
            simulated$mean_observed_prevalence_reference,
        observed_prevalence_comparison =
            simulated$mean_observed_prevalence_comparison,
        retained = extracted$retained,
        formed = extracted$formed,
        reason = extracted$reason,
        p_structural_absence = extracted$raw_p,
        q_bh = extracted$q_bh,
        q_by = extracted$q_by,
        reject_raw = extracted$formed &
            is.finite(extracted$raw_p) &
            extracted$raw_p <= config$alpha,
        reject_bh = extracted$retained &
            is.finite(extracted$q_bh) &
            extracted$q_bh <= config$alpha,
        reject_by = extracted$retained &
            is.finite(extracted$q_by) &
            extracted$q_by <= config$alpha,
        stringsAsFactors = FALSE
    )
}

run_replicate <- function(rep_id, config, scenarios, dgp_variants) {
    checkpoint_path <- file.path(
        checkpoint_dir,
        sprintf("replicate_%04d.rds", rep_id)
    )

    if (isTRUE(config$resume) && file.exists(checkpoint_path)) {
        saved <- readRDS(checkpoint_path)
        if (identical(saved$design_signature, design_signature)) {
            return(saved)
        }
    }

    metric_rows <- list()
    taxon_rows <- list()
    row_id <- 0L

    for (dgp in dgp_variants) {
        base <- make_base_replication(rep_id, dgp, config)

        for (s in seq_len(nrow(scenarios))) {
            scenario_row <- scenarios[s, , drop = FALSE]
            scenario <- scenario_row$scenario

            simulated <- simulate_scenario(
                base = base,
                scenario_row = scenario_row,
                rep_id = rep_id,
                dgp = dgp,
                config = config
            )

            extracted <- tryCatch({
                fit <- fit_dasra_structural(simulated)
                extract_structural_results(fit, rownames(simulated$counts))
            }, error = function(e) {
                make_failed_result(
                    rownames(simulated$counts),
                    conditionMessage(e)
                )
            })

            row_id <- row_id + 1L
            metric_rows[[row_id]] <- make_metric_rows(
                extracted = extracted,
                simulated = simulated,
                rep_id = rep_id,
                dgp = dgp,
                scenario = scenario,
                config = config
            )
            taxon_rows[[row_id]] <- make_taxon_rows(
                extracted = extracted,
                simulated = simulated,
                rep_id = rep_id,
                dgp = dgp,
                scenario = scenario,
                config = config
            )
        }
    }

    result <- list(
        design_signature = design_signature,
        metrics = dplyr::bind_rows(metric_rows),
        taxa = dplyr::bind_rows(taxon_rows)
    )
    saveRDS(result, checkpoint_path)
    result
}

message("DASRA structural-absence FDR stress simulation")
message("  Replications: ", CONFIG$n_rep)
message("  Samples per group: ", CONFIG$n_per_group)
message("  Taxa: ", CONFIG$n_taxa)
message("  Scenarios: ", paste(SCENARIOS$scenario, collapse = ", "))
message("  DGP variants: ", paste(DGP_VARIANTS, collapse = ", "))
message("  Workers: ", CONFIG$n_cores)
message("  Output: ", normalizePath(
    CONFIG$output_dir, winslash = "/", mustWork = FALSE
))
message("  DASRA version: ", utils::packageVersion("DASRA"))
message("  DASRA fingerprint: ", DASRA_FINGERPRINT)

replicate_ids <- seq_len(CONFIG$n_rep)

if (.Platform$OS.type == "unix" && CONFIG$n_cores > 1L) {
    replicate_results <- parallel::mclapply(
        replicate_ids,
        run_replicate,
        config = CONFIG,
        scenarios = SCENARIOS,
        dgp_variants = DGP_VARIANTS,
        mc.cores = CONFIG$n_cores,
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
    )
} else {
    if (CONFIG$n_cores > 1L) {
        warning(
            "Forked parallelism is unavailable on this platform; ",
            "running sequentially."
        )
    }
    replicate_results <- lapply(
        replicate_ids,
        run_replicate,
        config = CONFIG,
        scenarios = SCENARIOS,
        dgp_variants = DGP_VARIANTS
    )
}

replicate_metrics <- dplyr::bind_rows(lapply(
    replicate_results, `[[`, "metrics"
))
taxon_results <- dplyr::bind_rows(lapply(
    replicate_results, `[[`, "taxa"
))

safe_mean <- function(x) {
    if (!length(x) || all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    sd(x)
}

summary_results <- replicate_metrics %>%
    group_by(dgp, scenario, adjustment) %>%
    summarise(
        n_rep = n(),
        empirical_fdr = mean(fdp),
        mcse_fdr = safe_sd(fdp) / sqrt(n()),
        fdr_lower = pmax(0, empirical_fdr - 1.96 * mcse_fdr),
        fdr_upper = pmin(1, empirical_fdr + 1.96 * mcse_fdr),
        marginal_fwer = mean(any_false),
        marginal_fwer_mcse =
            sqrt(marginal_fwer * (1 - marginal_fwer) / n()),
        marginal_fwer_lower =
            pmax(0, marginal_fwer - 1.96 * marginal_fwer_mcse),
        marginal_fwer_upper =
            pmin(1, marginal_fwer + 1.96 * marginal_fwer_mcse),
        mFDR = sum(n_false) / max(sum(n_rejected), 1),
        mean_discoveries = mean(n_rejected),
        mean_false_discoveries = mean(n_false),
        structural_power = safe_mean(power),
        structural_power_mcse =
            safe_sd(power) / sqrt(sum(is.finite(power))),
        mean_retention_rate = mean(retention_rate),
        mean_formation_rate_among_retained =
            safe_mean(formation_rate_among_retained),
        mean_overall_formation_rate = mean(overall_formation_rate),
        mean_raw_null_rejection_rate =
            safe_mean(raw_null_rejection_rate),
        .groups = "drop"
    )

taxon_results <- taxon_results %>%
    mutate(
        taxon_class = case_when(
            structural_alt ~ "structural signal",
            abundance_alt ~ "abundance-shifted structural null",
            dispersion_alt ~ "dispersion-shifted structural null",
            TRUE ~ "other structural null"
        )
    )

taxon_class_summary <- taxon_results %>%
    group_by(dgp, scenario, taxon_class) %>%
    summarise(
        n_taxon_replications = n(),
        retention_rate = mean(retained),
        formation_rate = mean(formed),
        raw_rejection_rate = safe_mean(
            ifelse(formed, reject_raw, NA)
        ),
        bh_rejection_rate = mean(reject_bh),
        by_rejection_rate = mean(reject_by),
        mean_observed_prevalence_difference = mean(
            observed_prevalence_comparison -
                observed_prevalence_reference
        ),
        .groups = "drop"
    )

failure_summary <- taxon_results %>%
    filter(retained, !formed) %>%
    count(dgp, scenario, reason, name = "n") %>%
    group_by(dgp, scenario) %>%
    mutate(proportion_among_unformed = n / sum(n)) %>%
    ungroup()

write.csv(
    replicate_metrics,
    file.path(CONFIG$output_dir, "combined_replicate_metrics.csv"),
    row.names = FALSE
)
write.csv(
    taxon_results,
    file.path(CONFIG$output_dir, "combined_taxon_results.csv"),
    row.names = FALSE
)
write.csv(
    summary_results,
    file.path(CONFIG$output_dir, "combined_fdr_power_summary.csv"),
    row.names = FALSE
)
write.csv(
    taxon_class_summary,
    file.path(CONFIG$output_dir, "taxon_class_summary.csv"),
    row.names = FALSE
)
write.csv(
    failure_summary,
    file.path(CONFIG$output_dir, "structural_failure_summary.csv"),
    row.names = FALSE
)
write.csv(
    SCENARIOS,
    file.path(CONFIG$output_dir, "simulation_scenarios.csv"),
    row.names = FALSE
)

scenario_labels <- c(
    global_null = "Global null",
    sparse_abundance_only = "Sparse abundance only",
    dense_abundance_only = "Dense abundance only",
    dispersion_only = "Abundance dispersion only",
    structural_only = "Structural absence only",
    mixed_structural_abundance = "Mixed structural + abundance"
)

summary_results <- summary_results %>%
    mutate(
        scenario_label = factor(
            scenario_labels[scenario],
            levels = unname(scenario_labels)
        )
    )

base_theme <- theme_classic(base_size = 12) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 28, hjust = 1),
        plot.margin = margin(8, 12, 8, 8)
    )

fdr_plot <- ggplot(
    summary_results,
    aes(
        x = scenario_label,
        y = empirical_fdr,
        group = adjustment,
        linetype = adjustment,
        shape = adjustment
    )
) +
    geom_hline(yintercept = CONFIG$alpha, linetype = 3) +
    geom_errorbar(
        aes(ymin = fdr_lower, ymax = fdr_upper),
        width = 0.10,
        position = position_dodge(width = 0.18)
    ) +
    geom_line(position = position_dodge(width = 0.18)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.18)) +
    facet_wrap(~ dgp) +
    coord_cartesian(ylim = c(0, max(
        CONFIG$alpha * 2,
        summary_results$fdr_upper,
        na.rm = TRUE
    ))) +
    labs(
        x = NULL,
        y = "Empirical structural-absence FDR",
        linetype = "Adjustment",
        shape = "Adjustment"
    ) +
    base_theme

power_data <- summary_results %>%
    filter(scenario %in% c(
        "structural_only",
        "mixed_structural_abundance"
    ))

power_plot <- ggplot(
    power_data,
    aes(
        x = scenario_label,
        y = structural_power,
        group = adjustment,
        linetype = adjustment,
        shape = adjustment
    )
) +
    geom_errorbar(
        aes(
            ymin = pmax(0, structural_power -
                           1.96 * structural_power_mcse),
            ymax = pmin(1, structural_power +
                           1.96 * structural_power_mcse)
        ),
        width = 0.10,
        position = position_dodge(width = 0.18)
    ) +
    geom_line(position = position_dodge(width = 0.18)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.18)) +
    facet_wrap(~ dgp) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = NULL,
        y = "Structural-absence power",
        linetype = "Adjustment",
        shape = "Adjustment"
    ) +
    base_theme

availability_data <- summary_results %>%
    filter(adjustment == "BH") %>%
    select(
        dgp,
        scenario,
        scenario_label,
        mean_retention_rate,
        mean_formation_rate_among_retained,
        mean_overall_formation_rate
    ) %>%
    tidyr::pivot_longer(
        cols = c(
            mean_retention_rate,
            mean_formation_rate_among_retained,
            mean_overall_formation_rate
        ),
        names_to = "quantity",
        values_to = "rate"
    ) %>%
    mutate(
        quantity = factor(
            quantity,
            levels = c(
                "mean_retention_rate",
                "mean_formation_rate_among_retained",
                "mean_overall_formation_rate"
            ),
            labels = c(
                "Taxa retained",
                "Tests formed among retained taxa",
                "Tests formed among all taxa"
            )
        )
    )

availability_plot <- ggplot(
    availability_data,
    aes(
        x = scenario_label,
        y = rate,
        group = quantity,
        linetype = quantity,
        shape = quantity
    )
) +
    geom_line(position = position_dodge(width = 0.20)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.20)) +
    facet_wrap(~ dgp) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = NULL,
        y = "Proportion",
        linetype = NULL,
        shape = NULL
    ) +
    base_theme

class_plot_data <- taxon_class_summary %>%
    filter(taxon_class != "structural signal") %>%
    mutate(
        scenario_label = factor(
            scenario_labels[scenario],
            levels = unname(scenario_labels)
        )
    )

class_plot <- ggplot(
    class_plot_data,
    aes(
        x = scenario_label,
        y = raw_rejection_rate,
        group = taxon_class,
        linetype = taxon_class,
        shape = taxon_class
    )
) +
    geom_hline(yintercept = CONFIG$alpha, linetype = 3) +
    geom_line(position = position_dodge(width = 0.18)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.18)) +
    facet_wrap(~ dgp) +
    labs(
        x = NULL,
        y = "Raw structural-null rejection rate",
        linetype = NULL,
        shape = NULL
    ) +
    base_theme

global_null_p <- taxon_results %>%
    filter(
        scenario == "global_null",
        formed,
        is.finite(p_structural_absence)
    ) %>%
    group_by(dgp) %>%
    arrange(p_structural_absence, .by_group = TRUE) %>%
    mutate(
        expected = (row_number() - 0.5) / n()
    ) %>%
    ungroup()

qq_plot <- ggplot(
    global_null_p,
    aes(x = expected, y = p_structural_absence)
) +
    geom_abline(intercept = 0, slope = 1, linetype = 3) +
    geom_point(alpha = 0.35, size = 0.8) +
    facet_wrap(~ dgp) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
        x = "Expected uniform quantile",
        y = "Observed structural-absence p-value"
    ) +
    theme_classic(base_size = 12)

save_plot <- function(plot, stem, width, height) {
    pdf_device <- if (isTRUE(capabilities("cairo"))) {
        grDevices::cairo_pdf
    } else {
        "pdf"
    }
    ggsave(
        filename = file.path(CONFIG$output_dir, paste0(stem, ".pdf")),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = pdf_device
    )
    ggsave(
        filename = file.path(CONFIG$output_dir, paste0(stem, ".png")),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = 600
    )
}

save_plot(fdr_plot, "fig_structural_fdr", 10.0, 5.8)
save_plot(power_plot, "fig_structural_power", 7.5, 5.2)
save_plot(availability_plot, "fig_structural_availability", 10.0, 5.8)
save_plot(class_plot, "fig_structural_null_by_taxon_class", 10.5, 6.0)
save_plot(qq_plot, "fig_structural_global_null_qq", 6.5, 5.8)

session_path <- file.path(CONFIG$output_dir, "session_information.txt")
sink(session_path)
cat("DASRA structural-absence FDR stress simulation\n")
cat("Run date:", format(Sys.time()), "\n\n")
cat("Configuration:\n")
print(CONFIG)
cat("\nScenarios:\n")
print(SCENARIOS)
cat("\nDGP variants:\n")
print(DGP_VARIANTS)
cat("\nDASRA version:\n")
print(utils::packageVersion("DASRA"))
cat("\nDASRA implementation fingerprint:\n")
print(DASRA_FINGERPRINT)
cat("\nSession information:\n")
print(sessionInfo())
sink()

message("Simulation complete.")
message("Key outputs:")
message("  combined_fdr_power_summary.csv")
message("  combined_replicate_metrics.csv")
message("  combined_taxon_results.csv")
message("  taxon_class_summary.csv")
message("  structural_failure_summary.csv")
message("  fig_structural_fdr.pdf")
message("  fig_structural_power.pdf")
message("  fig_structural_availability.pdf")
message("  fig_structural_null_by_taxon_class.pdf")
message("  fig_structural_global_null_qq.pdf")
