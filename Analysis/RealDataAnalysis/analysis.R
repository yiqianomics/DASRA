options(stringsAsFactors = FALSE, warn = 1)

locate_analysis_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    candidates <- character()
    if (length(file_argument) == 1L) {
        script_path <- sub("^--file=", "", file_argument)
        candidates <- c(candidates, dirname(normalizePath(script_path)))
    }
    candidates <- unique(c(
        candidates,
        getwd(),
        file.path(getwd(), "Analysis", "RealDataAnalysis")
    ))
    expected <- file.path(
        "crc_baxter", "processed", "crc_baxter_dasra_input.rds"
    )
    matches <- candidates[file.exists(file.path(candidates, expected))]
    if (length(matches) != 1L) {
        stop("Could not identify the RealDataAnalysis directory.", call. = FALSE)
    }
    normalizePath(matches[[1L]])
}

analysis_directory <- locate_analysis_directory()
local_library <- file.path(analysis_directory, "R_lib")
if (dir.exists(local_library)) {
    .libPaths(c(normalizePath(local_library), .libPaths()))
}

required_packages <- c(
    "DASRA", "maaslin3", "ZINQ", "ANCOMBC", "MicrobiomeStat",
    "corncob", "edgeR", "DESeq2", "metagenomeSeq", "Biobase",
    "ggplot2", "dplyr", "tidyr", "patchwork"
)
missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
)]
if (length(missing_packages)) {
    stop(
        sprintf("Required R packages are unavailable: %s",
                paste(missing_packages, collapse = ", ")),
        call. = FALSE
    )
}

method_order <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "DASRA combined",
    "MaAsLin3 prevalence",
    "MaAsLin3 abundance",
    "MaAsLin3 combined",
    "ZINQ prevalence",
    "ZINQ abundance",
    "ZINQ combined",
    "ANCOM-BC2",
    "LinDA",
    "corncob",
    "edgeR",
    "DESeq2",
    "metagenomeSeq"
)

method_information <- data.frame(
    method = method_order,
    family = c(
        rep("DASRA", 3L), rep("MaAsLin3", 3L), rep("ZINQ", 3L),
        "ANCOM-BC2", "LinDA", "corncob", "edgeR", "DESeq2",
        "metagenomeSeq"
    ),
    component = c(
        "structural absence", "present-conditional abundance", "combined",
        "prevalence", "abundance", "combined",
        "prevalence", "abundance", "combined",
        rep("differential abundance", 6L)
    ),
    stringsAsFactors = FALSE
)

dataset_configurations <- list(
    crc_baxter = list(
        dataset_id = "crc_baxter",
        output_prefix = "baxter_crc",
        input_file = "crc_baxter_dasra_input.rds",
        reference = "H",
        comparison = "CRC",
        covariates = c("age_z", "sex")
    ),
    cdi_schubert = list(
        dataset_id = "cdi_schubert",
        output_prefix = "schubert_cdi",
        input_file = "cdi_schubert_dasra_input.rds",
        reference = "H",
        comparison = "CDI",
        covariates = c("age_z", "sex", "antibiotics_3mo")
    )
)

align_numeric <- function(values, taxa) {
    output <- stats::setNames(rep(NA_real_, length(taxa)), taxa)
    if (is.null(values)) return(output)
    value_names <- names(values)
    values <- suppressWarnings(as.numeric(values))
    if (is.null(value_names)) {
        if (length(values) == length(taxa)) output[] <- values
        return(output)
    }
    shared <- intersect(taxa, value_names)
    output[shared] <- values[match(shared, value_names)]
    output
}

align_logical <- function(values, taxa, default = FALSE) {
    output <- stats::setNames(rep(default, length(taxa)), taxa)
    if (is.null(values)) return(output)
    value_names <- names(values)
    values <- as.logical(values)
    if (is.null(value_names)) {
        if (length(values) == length(taxa)) output[] <- values
        return(output)
    }
    shared <- intersect(taxa, value_names)
    output[shared] <- values[match(shared, value_names)]
    output
}

align_character <- function(values, taxa, default = "not available") {
    output <- stats::setNames(rep(default, length(taxa)), taxa)
    if (is.null(values)) return(output)
    value_names <- names(values)
    values <- as.character(values)
    if (is.null(value_names)) {
        if (length(values) == length(taxa)) output[] <- values
        return(output)
    }
    shared <- intersect(taxa, value_names)
    output[shared] <- values[match(shared, value_names)]
    output
}

make_result <- function(taxa, method, family, component, p_value,
                        available, reason = NULL, estimate = NULL,
                        statistic = NULL, native_q_value = NULL,
                        effect_definition = "not reported",
                        components_used = NULL,
                        sensitivity_passed = NULL,
                        p_value_unfiltered = NULL) {
    p_value <- align_numeric(p_value, taxa)
    if (is.null(p_value_unfiltered)) p_value_unfiltered <- p_value
    p_value_unfiltered <- align_numeric(p_value_unfiltered, taxa)
    available <- align_logical(available, taxa)
    reason <- align_character(reason, taxa)
    estimate <- align_numeric(estimate, taxa)
    statistic <- align_numeric(statistic, taxa)
    native_q_value <- align_numeric(native_q_value, taxa)
    components_used <- align_character(components_used, taxa, NA_character_)
    sensitivity_passed <- align_logical(
        sensitivity_passed, taxa, default = NA
    )

    invalid_p <- available & (
        !is.finite(p_value) | p_value < 0 | p_value > 1
    )
    available[invalid_p] <- FALSE
    reason[invalid_p] <- "invalid p-value"
    missing_reason <- available & (
        is.na(reason) | !nzchar(reason) | reason == "not available"
    )
    reason[missing_reason] <- "available"

    p_for_adjustment <- ifelse(available, p_value, 1)
    q_value <- stats::p.adjust(p_for_adjustment, method = "BH")
    significant <- available & q_value <= 0.05

    data.frame(
        taxon = taxa,
        method = method,
        family = family,
        component = component,
        available = unname(available),
        reason = unname(reason),
        p_value = unname(p_value),
        p_value_unfiltered = unname(p_value_unfiltered),
        p_for_bh = unname(p_for_adjustment),
        q_value = unname(q_value),
        native_q_value = unname(native_q_value),
        significant = unname(significant),
        estimate = unname(estimate),
        statistic = unname(statistic),
        effect_definition = effect_definition,
        components_used = unname(components_used),
        sensitivity_passed = unname(sensitivity_passed),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

make_failed_results <- function(taxa, family, message_text) {
    rows <- method_information[method_information$family == family, , drop = FALSE]
    do.call(rbind, lapply(seq_len(nrow(rows)), function(index) {
        make_result(
            taxa = taxa,
            method = rows$method[[index]],
            family = family,
            component = rows$component[[index]],
            p_value = NULL,
            available = FALSE,
            reason = paste("method failed:", message_text)
        )
    }))
}

package_version_text <- function(package) {
    as.character(utils::packageVersion(package))
}

execute_family <- function(family, taxa, package, runner) {
    message("Running ", family)
    started <- proc.time()[["elapsed"]]
    outcome <- tryCatch(
        runner(),
        error = function(error) list(
            rows = make_failed_results(taxa, family, conditionMessage(error)),
            note = conditionMessage(error),
            failed = TRUE
        )
    )
    elapsed <- proc.time()[["elapsed"]] - started
    rows <- outcome$rows
    failed <- isTRUE(outcome$failed)
    if (!failed && !identical(sort(unique(rows$method)),
                              sort(method_information$method[
                                  method_information$family == family
                              ]))) {
        stop(sprintf("%s returned an unexpected result set.", family),
             call. = FALSE)
    }
    status <- if (failed) {
        "failed"
    } else if (any(!rows$available)) {
        "completed with unavailable taxa"
    } else {
        "completed"
    }
    list(
        rows = rows,
        status = data.frame(
            family = family,
            package = package,
            version = package_version_text(package),
            status = status,
            available_results = sum(rows$available),
            total_results = nrow(rows),
            discoveries = sum(rows$significant),
            elapsed_seconds = round(elapsed, 2),
            note = if (is.null(outcome$note)) "" else outcome$note,
            stringsAsFactors = FALSE
        )
    )
}

formulas_for_dataset <- function(configuration) {
    full_terms <- c(configuration$covariates, "group")
    null_terms <- configuration$covariates
    maaslin_terms <- c(
        configuration$covariates, "log_library_size_z", "group"
    )
    list(
        full = stats::as.formula(
            paste("~", paste(full_terms, collapse = " + "))
        ),
        null = stats::as.formula(
            paste("~", paste(null_terms, collapse = " + "))
        ),
        full_text = paste(full_terms, collapse = " + "),
        null_text = paste(null_terms, collapse = " + "),
        maaslin = stats::as.formula(
            paste("~", paste(maaslin_terms, collapse = " + "))
        ),
        maaslin_text = paste(maaslin_terms, collapse = " + ")
    )
}

run_dasra <- function(counts, metadata, configuration, formulas,
                      evaluation_taxa) {
    taxa <- evaluation_taxa
    fit <- DASRA::dasra(
        counts = counts,
        metadata = metadata,
        formula = formulas$full,
        group = "group",
        library_size = "library_size",
        taxa_are_rows = TRUE,
        reference = configuration$reference,
        p_adjust_method = "BH",
        component = "all",
        full_output = FALSE,
        conditional_present_starts = 1L
    )

    results <- fit$results[match(taxa, fit$results$taxon), , drop = FALSE]
    diagnostics <- fit$diagnostics[
        match(taxa, fit$diagnostics$taxon), , drop = FALSE
    ]
    if (anyNA(results$taxon) || anyNA(diagnostics$taxon)) {
        stop("DASRA did not return every taxon.", call. = FALSE)
    }

    names_by_taxon <- function(values) stats::setNames(values, taxa)
    structural_available <- names_by_taxon(
        diagnostics$retained & diagnostics$formed_structural_absence
    )
    abundance_available <- names_by_taxon(
        diagnostics$retained & diagnostics$formed_relative_abundance
    )
    combined_available <- structural_available & abundance_available

    structural_reason <- names_by_taxon(
        ifelse(
            !diagnostics$retained,
            "taxon not retained",
            diagnostics$reason_structural_absence
        )
    )
    abundance_reason <- names_by_taxon(
        ifelse(
            !diagnostics$retained,
            "taxon not retained",
            diagnostics$reason_relative_abundance
        )
    )
    combined_reason <- names_by_taxon(ifelse(
        combined_available,
        "available",
        "one or both components were unavailable"
    ))

    rows <- rbind(
        make_result(
            taxa, "DASRA structural absence", "DASRA",
            "structural absence",
            names_by_taxon(results$p_structural_absence),
            structural_available,
            structural_reason,
            statistic = names_by_taxon(results$z_structural_absence),
            native_q_value = names_by_taxon(results$q_structural_absence),
            effect_definition = paste(
                "score z; positive values indicate greater structural",
                "absence in the comparison group"
            )
        ),
        make_result(
            taxa, "DASRA present-conditional abundance", "DASRA",
            "present-conditional abundance",
            names_by_taxon(results$p_relative_abundance),
            abundance_available,
            abundance_reason,
            estimate = names_by_taxon(results$estimate_relative_abundance),
            statistic = names_by_taxon(results$z_relative_abundance),
            native_q_value = names_by_taxon(results$q_relative_abundance),
            effect_definition = paste(
                "comparison-minus-reference present-conditional log-relative",
                "abundance contrast against the shared compositional background"
            )
        ),
        make_result(
            taxa, "DASRA combined", "DASRA", "combined",
            names_by_taxon(results$p_omnibus),
            combined_available,
            combined_reason,
            native_q_value = names_by_taxon(results$q_omnibus),
            effect_definition = "Bonferroni minimum-p combination",
            components_used = names_by_taxon(results$components_used)
        )
    )
    list(
        rows = rows,
        note = paste(
            "Both components are required for the combined result;",
            "unavailable p-values remain in the BH family as one."
        )
    )
}

extract_maaslin_rows <- function(table, taxa, comparison) {
    table <- table[
        table$metadata == "group" & table$value == comparison,
        ,
        drop = FALSE
    ]
    p_value <- estimate <- joint_p <- native_q <- joint_q <-
        stats::setNames(rep(NA_real_, length(taxa)), taxa)
    available <- stats::setNames(rep(FALSE, length(taxa)), taxa)
    reason <- stats::setNames(rep("taxon not returned", length(taxa)), taxa)

    for (taxon in taxa) {
        index <- which(table$feature == taxon)
        if (length(index) != 1L) {
            if (length(index) > 1L) reason[[taxon]] <- "ambiguous contrast rows"
            next
        }
        row <- table[index, , drop = FALSE]
        row_error <- as.character(row$error[[1L]])
        p_value[[taxon]] <- row$pval_individual[[1L]]
        estimate[[taxon]] <- row$coef[[1L]]
        joint_p[[taxon]] <- row$pval_joint[[1L]]
        native_q[[taxon]] <- row$qval_individual[[1L]]
        joint_q[[taxon]] <- row$qval_joint[[1L]]
        available[[taxon]] <- is.finite(p_value[[taxon]])
        reason[[taxon]] <- if (available[[taxon]]) {
            "available"
        } else if (!is.na(row_error) && nzchar(row_error)) {
            paste("model error:", row_error)
        } else {
            "invalid p-value"
        }
    }
    list(
        p_value = p_value,
        estimate = estimate,
        joint_p = joint_p,
        native_q = native_q,
        joint_q = joint_q,
        available = available,
        reason = reason
    )
}

run_maaslin3 <- function(counts, metadata, configuration, formulas,
                         work_directory, evaluation_taxa) {
    taxa <- evaluation_taxa
    output_directory <- file.path(work_directory, "maaslin3")
    if (dir.exists(output_directory)) unlink(output_directory, recursive = TRUE)
    on.exit(try(maaslin3::maaslin_log_reset(), silent = TRUE), add = TRUE)

    fit <- maaslin3::maaslin3(
        input_data = as.data.frame(t(counts), check.names = FALSE),
        input_metadata = metadata,
        output = output_directory,
        formula = formulas$maaslin,
        min_abundance = 0,
        min_prevalence = 0,
        max_prevalence = 1.01,
        zero_threshold = 0,
        min_variance = 0,
        normalization = "TSS",
        transform = "LOG",
        correction = "BH",
        standardize = TRUE,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        warn_prevalence = TRUE,
        augment = TRUE,
        cores = 1,
        max_significance = 1,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        save_models = FALSE,
        save_plots_rds = FALSE,
        verbosity = "WARN",
        reference = paste0("group,", configuration$reference)
    )
    maaslin3::maaslin_log_reset()

    prevalence <- extract_maaslin_rows(
        fit$fit_data_prevalence$results, taxa, configuration$comparison
    )
    abundance <- extract_maaslin_rows(
        fit$fit_data_abundance$results, taxa, configuration$comparison
    )
    joint_agreement <- is.finite(prevalence$joint_p) &
        is.finite(abundance$joint_p) &
        abs(prevalence$joint_p - abundance$joint_p) <= 1e-12
    combined_available <- prevalence$available & abundance$available &
        joint_agreement
    combined_reason <- stats::setNames(ifelse(
        combined_available,
        "available",
        "both component models and a consistent joint p-value were required"
    ), taxa)

    rows <- rbind(
        make_result(
            taxa, "MaAsLin3 prevalence", "MaAsLin3", "prevalence",
            prevalence$p_value, prevalence$available, prevalence$reason,
            estimate = prevalence$estimate,
            native_q_value = prevalence$native_q,
            effect_definition = paste(
                "change in log odds of observed presence for the comparison",
                "group"
            )
        ),
        make_result(
            taxa, "MaAsLin3 abundance", "MaAsLin3", "abundance",
            abundance$p_value, abundance$available, abundance$reason,
            estimate = abundance$estimate,
            native_q_value = abundance$native_q,
            effect_definition = paste(
                "log2 relative-abundance coefficient tested against the",
                "feature-median compositional null"
            )
        ),
        make_result(
            taxa, "MaAsLin3 combined", "MaAsLin3", "combined",
            prevalence$joint_p, combined_available, combined_reason,
            native_q_value = prevalence$joint_q,
            effect_definition = "MaAsLin3 Beta(1,2) calibrated joint test",
            components_used = stats::setNames(
                ifelse(combined_available, "both", "not both"), taxa
            )
        )
    )
    list(
        rows = rows,
        note = paste(
            "TSS normalization and log2 transformation; abundance uses the",
            "default median comparison; both arms are required for the joint set."
        )
    )
}

quantile_cauchy_p <- function(p_values, taus) {
    if (length(p_values) != length(taus) ||
        any(!is.finite(p_values)) || any(p_values < 0 | p_values > 1)) {
        return(NA_real_)
    }
    p_values <- pmin(pmax(p_values, 1e-15), 1 - 1e-15)
    weights <- ifelse(taus <= 0.5, taus, 1 - taus)
    weights <- weights / sum(weights)
    statistic <- sum(weights * tan((0.5 - p_values) * pi))
    min(max(1 - stats::pcauchy(statistic), 0), 1)
}

run_zinq <- function(counts, metadata, configuration, evaluation_taxa) {
    taxa <- evaluation_taxa
    count_by_sample <- t(counts)
    sample_totals <- metadata$library_size
    relative_abundance <- sweep(count_by_sample, 1L, sample_totals, "/")
    group_binary <- as.integer(metadata$group == configuration$comparison)
    covariate_formula <- stats::as.formula(
        paste(
            "~",
            paste(
                c(configuration$covariates, "log_library_size_z"),
                collapse = " + "
            )
        )
    )
    covariate_matrix <- stats::model.matrix(
        covariate_formula, metadata
    )[, -1L, drop = FALSE]
    if (ncol(covariate_matrix)) {
        colnames(covariate_matrix) <- paste0("covariate", seq_len(
            ncol(covariate_matrix)
        ))
    }
    model_terms <- c("group", colnames(covariate_matrix))
    model_formula <- stats::as.formula(
        paste("y ~", paste(model_terms, collapse = " + "))
    )
    taus <- c(0.25, 0.50, 0.75)

    prevalence_p <- abundance_p <- combined_p <-
        stats::setNames(rep(NA_real_, length(taxa)), taxa)
    prevalence_available <- abundance_available <- combined_available <-
        stats::setNames(rep(FALSE, length(taxa)), taxa)
    prevalence_reason <- abundance_reason <- combined_reason <-
        stats::setNames(rep("model did not return a valid p-value", length(taxa)), taxa)

    for (index in seq_along(taxa)) {
        taxon <- taxa[[index]]
        model_data <- data.frame(
            y = as.numeric(relative_abundance[, taxon]),
            group = group_binary,
            stringsAsFactors = FALSE
        )
        if (ncol(covariate_matrix)) {
            model_data <- cbind(
                model_data,
                as.data.frame(covariate_matrix, check.names = FALSE)
            )
        }
        fit <- tryCatch(
            ZINQ::ZINQ_tests(
                formula.logistic = model_formula,
                formula.quantile = model_formula,
                C = "group",
                y_CorD = "C",
                data = model_data,
                taus = taus,
                seed = 2026L + index
            ),
            error = function(error) error
        )
        if (inherits(fit, "error")) {
            failure <- paste("model error:", conditionMessage(fit))
            prevalence_reason[[taxon]] <- failure
            abundance_reason[[taxon]] <- failure
            combined_reason[[taxon]] <- failure
            next
        }

        prevalence_p[[taxon]] <- as.numeric(fit$pvalue.logistic)[[1L]]
        abundance_p[[taxon]] <- quantile_cauchy_p(
            as.numeric(fit$pvalue.quantile), taus
        )
        combined_p[[taxon]] <- tryCatch(
            as.numeric(ZINQ::ZINQ_combination(
                fit, method = "Cauchy", taus = taus
            )),
            error = function(error) NA_real_
        )
        prevalence_available[[taxon]] <- is.finite(prevalence_p[[taxon]])
        abundance_available[[taxon]] <- is.finite(abundance_p[[taxon]])
        combined_available[[taxon]] <- prevalence_available[[taxon]] &
            abundance_available[[taxon]] & is.finite(combined_p[[taxon]])
        prevalence_reason[[taxon]] <- if (prevalence_available[[taxon]]) {
            "available"
        } else {
            "invalid logistic p-value"
        }
        abundance_reason[[taxon]] <- if (abundance_available[[taxon]]) {
            "available"
        } else {
            "invalid quantile component p-value"
        }
        combined_reason[[taxon]] <- if (combined_available[[taxon]]) {
            "available"
        } else {
            "both component models were not available"
        }
    }

    rows <- rbind(
        make_result(
            taxa, "ZINQ prevalence", "ZINQ", "prevalence",
            prevalence_p, prevalence_available, prevalence_reason,
            effect_definition = "Firth logistic test of observed presence"
        ),
        make_result(
            taxa, "ZINQ abundance", "ZINQ", "abundance",
            abundance_p, abundance_available, abundance_reason,
            effect_definition = paste(
                "weighted Cauchy combination of nonzero quantile tests at",
                "0.25, 0.50, and 0.75"
            )
        ),
        make_result(
            taxa, "ZINQ combined", "ZINQ", "combined",
            combined_p, combined_available, combined_reason,
            effect_definition = "native ZINQ Cauchy combination",
            components_used = stats::setNames(
                ifelse(combined_available, "both", "not both"), taxa
            )
        )
    )
    list(
        rows = rows,
        note = paste(
            "Relative abundance uses the original library size; models",
            "adjust for standardized log library size; taus 0.25, 0.50,",
            "and 0.75; native Cauchy joint test."
        )
    )
}

run_ancombc2 <- function(counts, metadata, configuration, formulas,
                         evaluation_taxa) {
    taxa <- evaluation_taxa
    fit <- ANCOMBC::ancombc2(
        data = counts,
        taxa_are_rows = TRUE,
        aggregate_data = counts,
        meta_data = metadata,
        fix_formula = formulas$full_text,
        rand_formula = NULL,
        p_adj_method = "BH",
        pseudo = 0,
        pseudo_sens = TRUE,
        prv_cut = 0,
        lib_cut = 0,
        s0_perc = 0.05,
        group = "group",
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        n_cl = 1,
        verbose = FALSE,
        global = FALSE,
        pairwise = FALSE,
        dunnet = FALSE,
        trend = FALSE,
        iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
        em_control = list(tol = 1e-05, max_iter = 100)
    )
    result <- fit$res
    contrast <- paste0("group", configuration$comparison)
    required_columns <- c(
        "taxon", paste0("lfc_", contrast), paste0("W_", contrast),
        paste0("p_", contrast), paste0("q_", contrast),
        paste0("passed_ss_", contrast)
    )
    missing_columns <- setdiff(required_columns, names(result))
    if (length(missing_columns)) {
        stop(
            sprintf("ANCOM-BC2 output is missing: %s",
                    paste(missing_columns, collapse = ", ")),
            call. = FALSE
        )
    }
    if (anyDuplicated(result$taxon)) {
        stop("ANCOM-BC2 returned duplicated taxon identifiers.", call. = FALSE)
    }
    row_index <- match(taxa, result$taxon)
    raw_p_value <- stats::setNames(
        result[[paste0("p_", contrast)]][row_index], taxa
    )
    native_q <- stats::setNames(
        result[[paste0("q_", contrast)]][row_index], taxa
    )
    estimate <- stats::setNames(
        result[[paste0("lfc_", contrast)]][row_index], taxa
    )
    statistic <- stats::setNames(
        result[[paste0("W_", contrast)]][row_index], taxa
    )
    sensitivity_passed <- stats::setNames(
        result[[paste0("passed_ss_", contrast)]][row_index], taxa
    )
    available <- stats::setNames(
        !is.na(row_index) & is.finite(raw_p_value), taxa
    )
    p_value <- raw_p_value
    sensitivity_failed <- available & !is.na(sensitivity_passed) &
        !sensitivity_passed
    p_value[sensitivity_failed] <- 1
    reason <- stats::setNames(ifelse(
        is.na(row_index),
        "taxon not returned",
        ifelse(
            sensitivity_failed,
            "pseudocount sensitivity failed; p-value set to one",
            ifelse(is.finite(raw_p_value), "available", "invalid p-value")
        )
    ), taxa)

    rows <- make_result(
        taxa, "ANCOM-BC2", "ANCOM-BC2", "differential abundance",
        p_value, available, reason,
        estimate = estimate,
        statistic = statistic,
        native_q_value = native_q,
        effect_definition = "bias-corrected natural-log fold change",
        sensitivity_passed = sensitivity_passed,
        p_value_unfiltered = raw_p_value
    )
    list(
        rows = rows,
        note = paste(
            "Structural-zero detection and negative lower-bound correction;",
            "a raw p-value is set to one when pseudocount sensitivity fails."
        )
    )
}

run_linda <- function(counts, metadata, configuration, formulas,
                      evaluation_taxa) {
    taxa <- evaluation_taxa
    fit <- MicrobiomeStat::linda(
        feature.dat = counts,
        meta.dat = metadata,
        formula = paste("~", formulas$full_text),
        feature.dat.type = "count",
        prev.filter = 0,
        mean.abund.filter = 0,
        max.abund.filter = 0,
        is.winsor = TRUE,
        outlier.pct = 0.03,
        adaptive = TRUE,
        zero.handling = "pseudo-count",
        pseudo.cnt = 0.5,
        p.adj.method = "BH",
        alpha = 0.05,
        n.cores = 1,
        verbose = FALSE
    )
    contrast <- paste0("group", configuration$comparison)
    if (!contrast %in% names(fit$output)) {
        stop("LinDA did not return the requested group contrast.", call. = FALSE)
    }
    result <- fit$output[[contrast]]
    if (anyDuplicated(rownames(result))) {
        stop("LinDA returned duplicated taxon identifiers.", call. = FALSE)
    }
    row_index <- match(taxa, rownames(result))
    p_value <- stats::setNames(result$pvalue[row_index], taxa)
    available <- stats::setNames(
        !is.na(row_index) & is.finite(p_value), taxa
    )
    reason <- stats::setNames(ifelse(
        is.na(row_index),
        "taxon not returned",
        ifelse(is.finite(p_value), "available", "invalid p-value")
    ), taxa)
    rows <- make_result(
        taxa, "LinDA", "LinDA", "differential abundance",
        p_value, available, reason,
        estimate = stats::setNames(result$log2FoldChange[row_index], taxa),
        statistic = stats::setNames(result$stat[row_index], taxa),
        native_q_value = stats::setNames(result$padj[row_index], taxa),
        effect_definition = "bias-corrected log2 fold change"
    )
    list(
        rows = rows,
        note = paste(
            "Adaptive bias correction, 3% winsorization, and a 0.5",
            "pseudo-count; package-level filters disabled."
        )
    )
}

run_corncob <- function(counts, metadata, configuration, formulas,
                        evaluation_taxa) {
    taxa <- evaluation_taxa
    p_value <- estimate <- statistic <-
        stats::setNames(rep(NA_real_, length(taxa)), taxa)
    available <- stats::setNames(rep(FALSE, length(taxa)), taxa)
    reason <- stats::setNames(rep("model unavailable", length(taxa)), taxa)
    mean_formula <- stats::as.formula(
        paste("cbind(W, M - W) ~", formulas$full_text)
    )
    dispersion_formula <- stats::as.formula(
        paste("~", formulas$full_text)
    )
    contrast <- paste0("mu.group", configuration$comparison)

    for (taxon in taxa) {
        model_data <- metadata[, c(configuration$covariates, "group"),
                               drop = FALSE]
        model_data$W <- as.numeric(counts[taxon, ])
        model_data$M <- as.numeric(metadata$library_size)
        fit <- tryCatch(
            corncob::bbdml(
                formula = mean_formula,
                phi.formula = dispersion_formula,
                data = model_data,
                method = "trust",
                robust = TRUE
            ),
            error = function(error) error
        )
        if (inherits(fit, "error")) {
            reason[[taxon]] <- paste("model error:", conditionMessage(fit))
            next
        }
        coefficient_table <- tryCatch(
            summary(fit)$coefficients,
            error = function(error) NULL
        )
        if (is.null(coefficient_table) ||
            !contrast %in% rownames(coefficient_table)) {
            reason[[taxon]] <- "group coefficient unavailable"
            next
        }
        p_column <- grep(
            "Pr\\(|p.value|p value", colnames(coefficient_table),
            ignore.case = TRUE, value = TRUE
        )
        statistic_column <- grep(
            "t value|z value|stat", colnames(coefficient_table),
            ignore.case = TRUE, value = TRUE
        )
        if (!length(p_column)) {
            reason[[taxon]] <- "group p-value unavailable"
            next
        }
        p_value[[taxon]] <- coefficient_table[contrast, p_column[[1L]]]
        estimate[[taxon]] <- coefficient_table[contrast, "Estimate"]
        if (length(statistic_column)) {
            statistic[[taxon]] <- coefficient_table[
                contrast, statistic_column[[1L]]
            ]
        }
        available[[taxon]] <- is.finite(p_value[[taxon]])
        reason[[taxon]] <- if (available[[taxon]]) {
            "available"
        } else {
            "invalid p-value"
        }
    }
    native_q <- stats::setNames(
        stats::p.adjust(ifelse(available, p_value, 1), method = "BH"),
        taxa
    )

    rows <- make_result(
        taxa, "corncob", "corncob", "differential abundance",
        p_value, available, reason,
        estimate = estimate,
        statistic = statistic,
        native_q_value = native_q,
        effect_definition = "beta-binomial mean-model logit coefficient"
    )
    list(
        rows = rows,
        note = paste(
            "Robust beta-binomial Wald inference; mean and dispersion models",
            "share the adjusted formula; original library size is the",
            "binomial denominator."
        )
    )
}

run_edger <- function(counts, metadata, configuration, formulas,
                      evaluation_taxa) {
    taxa <- evaluation_taxa
    design <- stats::model.matrix(formulas$full, metadata)
    contrast <- paste0("group", configuration$comparison)
    coefficient <- match(contrast, colnames(design))
    if (is.na(coefficient)) {
        stop("edgeR design does not contain the requested group contrast.",
             call. = FALSE)
    }
    fit_data <- edgeR::DGEList(counts = counts, group = metadata$group)
    keep <- edgeR::filterByExpr(fit_data, design = design)
    if (!any(keep)) {
        stop("edgeR filterByExpr retained no taxa.", call. = FALSE)
    }
    fit_data <- fit_data[keep, , keep.lib.sizes = FALSE]
    fit_data <- edgeR::calcNormFactors(fit_data)
    fit_data <- edgeR::estimateDisp(fit_data, design, robust = TRUE)
    fit <- edgeR::glmQLFit(fit_data, design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, coef = coefficient)
    result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
    row_index <- match(taxa, rownames(result))
    p_value <- stats::setNames(result$PValue[row_index], taxa)
    available <- stats::setNames(
        !is.na(row_index) & is.finite(p_value), taxa
    )
    reason <- stats::setNames(ifelse(
        is.na(row_index),
        "filterByExpr excluded taxon",
        ifelse(is.finite(p_value), "available", "invalid p-value")
    ), taxa)
    rows <- make_result(
        taxa, "edgeR", "edgeR", "differential abundance",
        p_value, available, reason,
        estimate = stats::setNames(result$logFC[row_index], taxa),
        statistic = stats::setNames(result$F[row_index], taxa),
        native_q_value = stats::setNames(result$FDR[row_index], taxa),
        effect_definition = "TMM-normalized log2 fold change"
    )
    list(
        rows = rows,
        note = paste(
            "TMM normalization with robust quasi-likelihood fitting;",
            "filterByExpr exclusions remain in the common BH family as one."
        )
    )
}

run_deseq2 <- function(counts, metadata, configuration, formulas,
                       evaluation_taxa) {
    taxa <- evaluation_taxa
    fit_data <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = metadata,
        design = formulas$full
    )
    fit <- DESeq2::DESeq(
        fit_data,
        test = "Wald",
        quiet = TRUE,
        fitType = "parametric",
        sfType = "poscounts",
        betaPrior = FALSE
    )
    result <- DESeq2::results(
        fit,
        contrast = c(
            "group", configuration$comparison, configuration$reference
        ),
        independentFiltering = TRUE,
        cooksCutoff = TRUE,
        alpha = 0.05
    )
    row_index <- match(taxa, rownames(result))
    p_value <- stats::setNames(result$pvalue[row_index], taxa)
    available <- stats::setNames(
        !is.na(row_index) & is.finite(p_value), taxa
    )
    reason <- stats::setNames(ifelse(
        is.na(row_index),
        "taxon not returned",
        ifelse(is.finite(p_value), "available", "invalid p-value")
    ), taxa)
    rows <- make_result(
        taxa, "DESeq2", "DESeq2", "differential abundance",
        p_value, available, reason,
        estimate = stats::setNames(result$log2FoldChange[row_index], taxa),
        statistic = stats::setNames(result$stat[row_index], taxa),
        native_q_value = stats::setNames(result$padj[row_index], taxa),
        effect_definition = "DESeq2 log2 fold change"
    )
    list(
        rows = rows,
        note = paste(
            "Positive-count size factors, parametric dispersion fit, and",
            "Wald test; Cook's handling and independent filtering retained."
        )
    )
}

run_metagenomeseq <- function(counts, metadata, configuration, formulas,
                              evaluation_taxa) {
    taxa <- evaluation_taxa
    phenotype <- metadata[, c(configuration$covariates, "group"), drop = FALSE]
    experiment <- metagenomeSeq::newMRexperiment(
        counts = counts,
        phenoData = Biobase::AnnotatedDataFrame(phenotype)
    )
    normalization_fallback <- FALSE
    normalization_quantile <- tryCatch(
        metagenomeSeq::cumNormStatFast(experiment, pFlag = FALSE),
        error = function(error) {
            normalization_fallback <<- TRUE
            0.5
        }
    )
    if (!is.finite(normalization_quantile)) {
        normalization_quantile <- 0.5
        normalization_fallback <- TRUE
    }
    experiment <- metagenomeSeq::cumNorm(
        experiment, p = normalization_quantile
    )
    design <- stats::model.matrix(
        formulas$full, data = Biobase::pData(experiment)
    )
    contrast <- paste0("group", configuration$comparison)
    coefficient <- match(contrast, colnames(design))
    if (is.na(coefficient)) {
        stop(
            "metagenomeSeq design does not contain the requested contrast.",
            call. = FALSE
        )
    }
    fit <- metagenomeSeq::fitZig(
        experiment,
        design,
        control = metagenomeSeq::zigControl(maxit = 10, verbose = FALSE)
    )
    result <- as.data.frame(metagenomeSeq::MRcoefs(
        fit,
        by = coefficient,
        number = nrow(counts),
        coef = coefficient,
        group = 4,
        adjustMethod = "BH",
        eff = 0,
        counts = 0
    ))
    p_column <- grep(
        "^pvalues$|^pvalue$|^p\\.value$|^PValue$",
        names(result), value = TRUE
    )
    if (length(p_column) != 1L) {
        stop("metagenomeSeq output has no unique p-value column.",
             call. = FALSE)
    }
    effect_column <- if (contrast %in% names(result)) {
        contrast
    } else {
        setdiff(names(result), c(p_column, "adjPvalues"))[[1L]]
    }
    row_index <- match(taxa, rownames(result))
    p_value <- stats::setNames(result[[p_column]][row_index], taxa)
    available <- stats::setNames(
        !is.na(row_index) & is.finite(p_value), taxa
    )
    reason <- stats::setNames(ifelse(
        is.na(row_index),
        "taxon not returned",
        ifelse(is.finite(p_value), "available", "invalid p-value")
    ), taxa)
    rows <- make_result(
        taxa, "metagenomeSeq", "metagenomeSeq",
        "differential abundance",
        p_value, available, reason,
        estimate = stats::setNames(result[[effect_column]][row_index], taxa),
        native_q_value = stats::setNames(
            if ("adjPvalues" %in% names(result)) {
                result$adjPvalues[row_index]
            } else {
                rep(NA_real_, length(taxa))
            },
            taxa
        ),
        effect_definition = "CSS-normalized zero-inflated Gaussian coefficient"
    )
    list(
        rows = rows,
        note = if (normalization_fallback) {
            paste(
                "cumNormStatFast was unavailable; CSS normalization used",
                "quantile 0.5000 followed by fitZig."
            )
        } else {
            sprintf(
                "CSS normalization at quantile %.4f followed by fitZig.",
                normalization_quantile
            )
        }
    )
}

validate_analysis_input <- function(input, configuration) {
    if (!is.list(input) || !all(c("counts", "metadata") %in% names(input))) {
        stop("The processed input must contain counts and metadata.",
             call. = FALSE)
    }
    counts <- input$counts
    metadata <- input$metadata
    if (!is.matrix(counts) || !is.numeric(counts) || anyNA(counts) ||
        any(!is.finite(counts)) || any(counts < 0) ||
        any(abs(counts - round(counts)) > 1e-8)) {
        stop("Counts must be a nonnegative integer-valued matrix.",
             call. = FALSE)
    }
    if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
        anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
        stop("Count row and column names must be present and unique.",
             call. = FALSE)
    }
    if (!is.data.frame(metadata) || is.null(rownames(metadata)) ||
        anyDuplicated(rownames(metadata))) {
        stop("Metadata row names must be present and unique.", call. = FALSE)
    }
    if (!setequal(colnames(counts), rownames(metadata))) {
        stop("Count samples and metadata rows do not agree.", call. = FALSE)
    }
    metadata <- metadata[colnames(counts), , drop = FALSE]
    required_columns <- c(
        "group", "library_size", configuration$covariates
    )
    missing_columns <- setdiff(required_columns, names(metadata))
    if (length(missing_columns)) {
        stop(
            sprintf("Metadata is missing: %s",
                    paste(missing_columns, collapse = ", ")),
            call. = FALSE
        )
    }
    if (any(!stats::complete.cases(metadata[, required_columns,
                                               drop = FALSE]))) {
        stop("Analysis covariates must be complete.", call. = FALSE)
    }
    metadata$group <- factor(
        as.character(metadata$group),
        levels = c(configuration$reference, configuration$comparison)
    )
    if (anyNA(metadata$group) || any(table(metadata$group) == 0L)) {
        stop("Both prespecified comparison groups must be present.",
             call. = FALSE)
    }
    if (any(!is.finite(metadata$library_size)) ||
        any(metadata$library_size <= 0)) {
        stop("Original library sizes must be positive and finite.",
             call. = FALSE)
    }
    retained_depth <- colSums(counts)
    if (any(metadata$library_size + 1e-8 < retained_depth)) {
        stop("A retained count total exceeds its original library size.",
             call. = FALSE)
    }
    metadata$log_library_size_z <- as.numeric(scale(
        log(metadata$library_size)
    ))
    counts <- round(counts)
    storage.mode(counts) <- "integer"
    remainder <- round(metadata$library_size - retained_depth)
    fitting_counts <- rbind(counts, Other_unmodeled = remainder)
    storage.mode(fitting_counts) <- "integer"
    list(
        counts = counts,
        fitting_counts = fitting_counts,
        metadata = metadata,
        retained_fraction = retained_depth / metadata$library_size,
        remainder = remainder
    )
}

make_upset_figure <- function(results, input, configuration,
                              table_directory, figure_directory) {
    expected_rows <- length(method_order) * nrow(input$counts)
    if (nrow(results) != expected_rows ||
        anyDuplicated(results[, c("taxon", "method")])) {
        stop("The UpSet input does not contain one row per taxon and method.",
             call. = FALSE)
    }
    results$method <- factor(results$method, levels = method_order)
    results <- results[order(results$method, results$taxon), , drop = FALSE]
    signature <- tidyr::pivot_wider(
        results[, c("taxon", "method", "significant")],
        names_from = "method",
        values_from = "significant",
        values_fill = FALSE
    )
    signature <- signature[, c("taxon", method_order), drop = FALSE]
    significant_matrix <- as.matrix(signature[, method_order, drop = FALSE])
    storage.mode(significant_matrix) <- "logical"
    has_discovery <- rowSums(significant_matrix) > 0L

    if (any(has_discovery)) {
        pattern <- apply(
            significant_matrix[has_discovery, , drop = FALSE],
            1L,
            function(row) paste(method_order[row], collapse = " | ")
        )
        membership <- data.frame(
            taxon = signature$taxon[has_discovery],
            methods = unname(pattern),
            stringsAsFactors = FALSE
        )
        pattern_table <- as.data.frame(table(membership$methods),
                                       stringsAsFactors = FALSE)
        names(pattern_table) <- c("methods", "n")
        pattern_table <- pattern_table[
            order(-pattern_table$n, pattern_table$methods), , drop = FALSE
        ]
        pattern_table$intersection <- paste0("I", seq_len(nrow(pattern_table)))
        membership$intersection <- pattern_table$intersection[
            match(membership$methods, pattern_table$methods)
        ]
    } else {
        pattern_table <- data.frame(
            methods = "No discoveries", n = 0L, intersection = "I1",
            stringsAsFactors = FALSE
        )
        membership <- data.frame(
            taxon = character(), methods = character(),
            intersection = character(), stringsAsFactors = FALSE
        )
    }
    pattern_table <- pattern_table[, c("intersection", "n", "methods")]
    utils::write.csv(
        pattern_table,
        file.path(
            table_directory,
            paste0(configuration$output_prefix, "_upset_intersections.csv")
        ),
        row.names = FALSE,
        na = ""
    )
    utils::write.csv(
        membership[, c("intersection", "taxon", "methods"), drop = FALSE],
        file.path(
            table_directory,
            paste0(configuration$output_prefix,
                   "_upset_intersection_membership.csv")
        ),
        row.names = FALSE,
        na = ""
    )

    positions <- data.frame(
        method = method_order,
        y = rev(seq_along(method_order)),
        stringsAsFactors = FALSE
    )
    intersections <- pattern_table$intersection
    pattern_table$intersection <- factor(
        pattern_table$intersection, levels = intersections
    )
    all_points <- merge(
        expand.grid(
            intersection = intersections,
            method = method_order,
            stringsAsFactors = FALSE
        ),
        positions,
        by = "method",
        sort = FALSE
    )
    all_points$intersection <- factor(
        all_points$intersection, levels = intersections
    )

    if (nrow(membership)) {
        active_points <- do.call(rbind, lapply(
            seq_len(nrow(pattern_table)),
            function(index) {
                methods <- strsplit(
                    pattern_table$methods[[index]], " | ", fixed = TRUE
                )[[1L]]
                data.frame(
                    intersection = pattern_table$intersection[[index]],
                    method = methods,
                    stringsAsFactors = FALSE
                )
            }
        ))
        active_points <- merge(
            active_points, positions, by = "method", sort = FALSE
        )
        active_points$intersection <- factor(
            active_points$intersection, levels = intersections
        )
        segment_data <- stats::aggregate(
            y ~ intersection, active_points,
            function(values) c(minimum = min(values), maximum = max(values))
        )
        segment_data <- data.frame(
            intersection = segment_data$intersection,
            minimum = segment_data$y[, "minimum"],
            maximum = segment_data$y[, "maximum"]
        )
    } else {
        active_points <- data.frame(
            intersection = factor(levels = intersections),
            method = character(), y = numeric()
        )
        segment_data <- data.frame(
            intersection = factor(levels = intersections),
            minimum = numeric(), maximum = numeric()
        )
    }

    method_counts <- stats::aggregate(
        significant ~ method, results, sum
    )
    method_counts <- merge(
        positions, method_counts, by = "method", all.x = TRUE, sort = FALSE
    )
    method_counts <- method_counts[match(method_order, method_counts$method), ]
    names(method_counts)[names(method_counts) == "significant"] <- "discoveries"
    maximum_method_count <- max(1, method_counts$discoveries)
    maximum_intersection <- max(1, pattern_table$n)
    method_counts$label_inside <- method_counts$discoveries >=
        0.08 * maximum_method_count
    method_counts$label_x <- ifelse(
        method_counts$label_inside,
        method_counts$discoveries / 2,
        method_counts$discoveries + 0.02 * maximum_method_count
    )

    method_colors <- c(
        "DASRA structural absence" = "#B2182B",
        "DASRA present-conditional abundance" = "#EF8A62",
        "DASRA combined" = "#D6604D",
        "MaAsLin3 prevalence" = "#2166AC",
        "MaAsLin3 abundance" = "#67A9CF",
        "MaAsLin3 combined" = "#4393C3",
        "ZINQ prevalence" = "#1B7837",
        "ZINQ abundance" = "#A6DBA0",
        "ZINQ combined" = "#5AAE61",
        "ANCOM-BC2" = "#762A83",
        "LinDA" = "#9970AB",
        "corncob" = "#5E3C99",
        "edgeR" = "#008837",
        "DESeq2" = "#E08214",
        "metagenomeSeq" = "#6B6B6B"
    )
    row_background <- data.frame(
        y = positions$y,
        fill = rep(c("grey98", "grey94"), length.out = nrow(positions))
    )
    common_theme <- ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(face = "bold"),
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            plot.margin = ggplot2::margin(6, 7, 6, 7)
        )
    method_bar <- ggplot2::ggplot(method_counts) +
        ggplot2::geom_rect(
            ggplot2::aes(
                xmin = 0, xmax = discoveries,
                ymin = y - 0.33, ymax = y + 0.33, fill = method
            )
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                x = label_x,
                y = y,
                label = discoveries,
                color = label_inside
            ),
            size = 3
        ) +
        ggplot2::scale_fill_manual(values = method_colors) +
        ggplot2::scale_color_manual(values = c("TRUE" = "white",
                                               "FALSE" = "grey20")) +
        ggplot2::scale_x_reverse(
            limits = c(maximum_method_count * 1.08, 0),
            expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0.5, length(method_order) + 0.5),
            breaks = positions$y,
            labels = positions$method,
            position = "right",
            expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::labs(x = "BH-significant taxa per method", y = NULL) +
        common_theme +
        ggplot2::theme(
            legend.position = "none",
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(
                size = 8.8, hjust = 0,
                margin = ggplot2::margin(l = 4, r = 2)
            ),
            panel.grid.major.x = ggplot2::element_line(
                color = "grey90", linewidth = 0.3
            ),
            plot.margin = ggplot2::margin(4, 3, 6, 6)
        )

    intersection_bar <- ggplot2::ggplot(
        pattern_table,
        ggplot2::aes(x = intersection, y = n)
    ) +
        ggplot2::geom_col(width = 0.66, fill = "grey22") +
        ggplot2::geom_text(
            ggplot2::aes(label = n), vjust = -0.25,
            size = 2.8, color = "grey15"
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, maximum_intersection * 1.16),
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
            title = "Significant-set intersections",
            x = NULL,
            y = "Intersection size"
        ) +
        common_theme +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                color = "grey90", linewidth = 0.35
            ),
            plot.margin = ggplot2::margin(6, 6, 4, 3)
        )

    intersection_matrix <- ggplot2::ggplot() +
        ggplot2::geom_rect(
            data = row_background,
            ggplot2::aes(
                xmin = -Inf, xmax = Inf,
                ymin = y - 0.5, ymax = y + 0.5, fill = fill
            ),
            color = NA
        ) +
        ggplot2::scale_fill_identity() +
        ggplot2::geom_point(
            data = all_points,
            ggplot2::aes(x = intersection, y = y),
            color = "grey82", size = 2.25
        ) +
        ggplot2::geom_segment(
            data = segment_data,
            ggplot2::aes(
                x = intersection, xend = intersection,
                y = minimum, yend = maximum
            ),
            color = "grey12", linewidth = 0.42
        ) +
        ggplot2::geom_point(
            data = active_points,
            ggplot2::aes(x = intersection, y = y),
            color = "grey12", size = 2.65
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0.5, length(method_order) + 0.5),
            breaks = NULL,
            expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::labs(
            x = "Intersections (sorted by size)",
            y = NULL
        ) +
        common_theme +
        ggplot2::theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 7),
            plot.margin = ggplot2::margin(4, 6, 6, 3)
        )

    figure_width <- min(
        22,
        max(17.5, 11.2 + 0.21 * nrow(pattern_table))
    )
    left_column_width <- 5.1
    figure <- patchwork::wrap_plots(
        A = patchwork::plot_spacer(),
        B = intersection_bar,
        C = method_bar,
        D = intersection_matrix,
        design = "
            AB
            CD
        ",
        widths = c(
            left_column_width,
            figure_width - left_column_width
        ),
        heights = c(0.55, 1.45)
    )
    figure_file <- file.path(
        figure_directory,
        paste0(configuration$output_prefix, "_upset.pdf")
    )
    ggplot2::ggsave(
        filename = figure_file,
        plot = figure,
        width = figure_width,
        height = 8.4,
        units = "in",
        device = grDevices::cairo_pdf,
        bg = "white"
    )
    list(
        file = figure_file,
        width = figure_width,
        intersections = nrow(pattern_table),
        taxa_in_union = nrow(membership)
    )
}

run_dataset <- function(configuration, plot_only = FALSE) {
    dataset_directory <- file.path(
        analysis_directory, configuration$dataset_id
    )
    input_path <- file.path(
        dataset_directory, "processed", configuration$input_file
    )
    table_directory <- file.path(dataset_directory, "table")
    figure_directory <- file.path(dataset_directory, "figs")
    work_directory <- file.path(dataset_directory, "work")
    dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
    dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)

    validated <- validate_analysis_input(
        readRDS(input_path), configuration
    )
    evaluation_taxa <- rownames(validated$counts)
    results_file <- file.path(
        table_directory,
        paste0(configuration$output_prefix, "_method_results_all_taxa.csv")
    )
    if (plot_only) {
        if (!file.exists(results_file)) {
            stop("Plot-only mode requires an existing method-results table.",
                 call. = FALSE)
        }
        results <- utils::read.csv(
            results_file, stringsAsFactors = FALSE, check.names = FALSE
        )
        results$significant <- as.logical(results$significant)
        figure <- make_upset_figure(
            results, validated, configuration,
            table_directory, figure_directory
        )
        message("Updated ", figure$file)
        return(invisible(figure))
    }

    formulas <- formulas_for_dataset(configuration)
    counts <- validated$fitting_counts
    metadata <- validated$metadata
    set.seed(20260820L + match(
        configuration$dataset_id, names(dataset_configurations)
    ))

    families <- list(
        execute_family(
            "DASRA", evaluation_taxa, "DASRA",
            function() run_dasra(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "MaAsLin3", evaluation_taxa, "maaslin3",
            function() run_maaslin3(
                counts, metadata, configuration, formulas,
                work_directory, evaluation_taxa
            )
        ),
        execute_family(
            "ZINQ", evaluation_taxa, "ZINQ",
            function() run_zinq(
                counts, metadata, configuration, evaluation_taxa
            )
        ),
        execute_family(
            "ANCOM-BC2", evaluation_taxa, "ANCOMBC",
            function() run_ancombc2(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "LinDA", evaluation_taxa, "MicrobiomeStat",
            function() run_linda(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "corncob", evaluation_taxa, "corncob",
            function() run_corncob(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "edgeR", evaluation_taxa, "edgeR",
            function() run_edger(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "DESeq2", evaluation_taxa, "DESeq2",
            function() run_deseq2(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        ),
        execute_family(
            "metagenomeSeq", evaluation_taxa, "metagenomeSeq",
            function() run_metagenomeseq(
                counts, metadata, configuration, formulas, evaluation_taxa
            )
        )
    )
    results <- do.call(rbind, lapply(families, `[[`, "rows"))
    status <- do.call(rbind, lapply(families, `[[`, "status"))
    results$method <- factor(results$method, levels = method_order)
    results <- results[order(results$method, results$taxon), , drop = FALSE]
    results$method <- as.character(results$method)
    results$dataset <- configuration$dataset_id
    results$reference_group <- configuration$reference
    results$comparison_group <- configuration$comparison
    results$model_formula <- ifelse(
        results$family %in% c("MaAsLin3", "ZINQ"),
        formulas$maaslin_text,
        formulas$full_text
    )
    results <- results[, c(
        "dataset", "taxon", "method", "family", "component",
        "reference_group", "comparison_group", "model_formula",
        "available", "reason", "p_value", "p_value_unfiltered",
        "p_for_bh", "q_value", "native_q_value", "significant",
        "estimate", "statistic", "effect_definition",
        "components_used", "sensitivity_passed"
    )]
    status$dataset <- configuration$dataset_id
    status <- status[, c(
        "dataset", "family", "package", "version", "status",
        "available_results", "total_results", "discoveries",
        "elapsed_seconds", "note"
    )]
    availability_summary <- stats::aggregate(
        cbind(available, significant) ~ method + family + component,
        results,
        sum
    )
    availability_summary <- availability_summary[
        match(method_order, availability_summary$method), , drop = FALSE
    ]
    discovery_summary <- availability_summary[
        , c("method", "family", "component", "significant"), drop = FALSE
    ]
    names(discovery_summary)[
        names(discovery_summary) == "significant"
    ] <- "discoveries"
    input_summary <- data.frame(
        item = c(
            "dataset", "reference samples", "comparison samples",
            "tested taxa", "normalization taxa",
            "median retained fraction of original reads",
            "minimum retained fraction of original reads",
            "multiplicity family", "significance threshold"
        ),
        value = c(
            configuration$dataset_id,
            sum(metadata$group == configuration$reference),
            sum(metadata$group == configuration$comparison),
            nrow(validated$counts),
            nrow(validated$fitting_counts),
            format(median(validated$retained_fraction), digits = 5),
            format(min(validated$retained_fraction), digits = 5),
            paste0(nrow(validated$counts),
                   " prespecified taxa per result set"),
            "BH-adjusted q <= 0.05"
        ),
        stringsAsFactors = FALSE
    )

    utils::write.csv(results, results_file, row.names = FALSE, na = "")
    status_for_export <- status[, setdiff(
        names(status), c("available_results", "total_results")
    ), drop = FALSE]
    utils::write.csv(
        status_for_export,
        file.path(
            table_directory,
            paste0(configuration$output_prefix, "_method_status.csv")
        ),
        row.names = FALSE,
        na = ""
    )
    utils::write.csv(
        discovery_summary,
        file.path(
            table_directory,
            paste0(configuration$output_prefix, "_discovery_summary.csv")
        ),
        row.names = FALSE,
        na = ""
    )
    utils::write.csv(
        input_summary,
        file.path(
            table_directory,
            paste0(configuration$output_prefix, "_analysis_input_summary.csv")
        ),
        row.names = FALSE,
        na = ""
    )

    if (any(status$status == "failed") ||
        any(availability_summary$available == 0L)) {
        stop(
            sprintf(
                "%s did not complete every required result set; inspect %s.",
                configuration$dataset_id,
                file.path(
                    table_directory,
                    paste0(configuration$output_prefix, "_method_status.csv")
                )
            ),
            call. = FALSE
        )
    }
    figure <- make_upset_figure(
        results, validated, configuration,
        table_directory, figure_directory
    )
    message(
        configuration$dataset_id, ": ", figure$taxa_in_union,
        " taxa in ", figure$intersections, " intersections; ",
        figure$file
    )
    invisible(list(
        results = results,
        status = status,
        discovery_summary = discovery_summary,
        figure = figure
    ))
}

arguments <- commandArgs(trailingOnly = TRUE)
plot_only <- "--plot-only" %in% arguments
dataset_arguments <- grep("^--dataset=", arguments, value = TRUE)
unknown_arguments <- setdiff(
    arguments, c("--plot-only", dataset_arguments)
)
if (length(unknown_arguments)) {
    stop(
        sprintf("Unknown argument: %s", unknown_arguments[[1L]]),
        call. = FALSE
    )
}
selected_datasets <- if (length(dataset_arguments)) {
    unique(sub("^--dataset=", "", dataset_arguments))
} else {
    names(dataset_configurations)
}
if (length(dataset_arguments) > 1L ||
    any(!selected_datasets %in% names(dataset_configurations))) {
    stop(
        sprintf(
            "Dataset must be one of: %s.",
            paste(names(dataset_configurations), collapse = ", ")
        ),
        call. = FALSE
    )
}
for (dataset_id in selected_datasets) {
    run_dataset(dataset_configurations[[dataset_id]], plot_only = plot_only)
}
