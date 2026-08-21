options(stringsAsFactors = FALSE, warn = 1)

# Template RDS contract:
# - probability: sample-by-taxon nonnegative matrix containing the full empirical
#   composition used for multinomial sampling. Rows are normalized here.
# - evaluation_taxa: exactly 30 column names defining the fixed test family.
# - dataset_id and label: stable machine-readable and display identifiers.
# - dataset_index: optional positive integer used in the deterministic seed base.

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x)) y else x
}

locate_analysis_directory <- function() {
    file_argument <- grep(
        "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
    )
    candidates <- character()
    if (length(file_argument) == 1L) {
        script_path <- sub("^--file=", "", file_argument)
        candidates <- c(candidates, dirname(normalizePath(script_path)))
    }
    candidates <- unique(c(
        candidates,
        file.path(getwd(), "Analysis", "NegativeControlAnalysis"),
        getwd()
    ))
    matches <- candidates[
        dir.exists(file.path(candidates, "data")) &
            dir.exists(file.path(candidates, "results"))
    ]
    if (!length(matches)) {
        stop("Could not identify the NegativeControlAnalysis directory.",
             call. = FALSE)
    }
    normalizePath(matches[[1L]])
}

analysis_directory <- locate_analysis_directory()
local_library <- file.path(analysis_directory, "R_lib")
if (dir.exists(local_library)) {
    .libPaths(c(normalizePath(local_library), .libPaths()))
}

scenario_table <- data.frame(
    scenario = c("balanced", "fourfold", "sixfold"),
    median_H = c(4000, 1500, 1000),
    median_Case = c(4000, 6000, 6000),
    scenario_index = seq_len(3L),
    stringsAsFactors = FALSE
)

method_names <- c("DASRA", "ZINQ Cauchy", "MaAsLin3 joint")
reference_group <- "H"
comparison_group <- "Case"
depth_sdlog <- 0.45
minimum_depth <- 300L
maximum_depth <- 30000L

parse_arguments <- function(arguments) {
    defaults <- list(
        dataset = NULL,
        reps = 100L,
        rep_start = 1L,
        workers = 1L,
        force = FALSE,
        keep_work = FALSE,
        help = FALSE
    )
    aliases <- c(
        "rep-start" = "rep_start",
        "keep-work" = "keep_work"
    )
    index <- 1L
    while (index <= length(arguments)) {
        argument <- arguments[[index]]
        if (!startsWith(argument, "--")) {
            stop(sprintf("Unexpected positional argument: %s", argument),
                 call. = FALSE)
        }
        item <- substring(argument, 3L)
        if (grepl("=", item, fixed = TRUE)) {
            pieces <- strsplit(item, "=", fixed = TRUE)[[1L]]
            key <- pieces[[1L]]
            value <- paste(pieces[-1L], collapse = "=")
        } else {
            key <- item
            if (key %in% c("force", "keep-work", "help")) {
                value <- "true"
            } else {
                index <- index + 1L
                if (index > length(arguments)) {
                    stop(sprintf("Missing value for --%s", key),
                         call. = FALSE)
                }
                value <- arguments[[index]]
            }
        }
        if (key %in% names(aliases)) key <- unname(aliases[[key]])
        if (!key %in% names(defaults)) {
            stop(sprintf("Unknown option: --%s", key), call. = FALSE)
        }
        if (key %in% c("force", "keep_work", "help")) {
            defaults[[key]] <- tolower(value) %in% c("1", "true", "yes")
        } else if (key %in% c("reps", "rep_start", "workers")) {
            defaults[[key]] <- suppressWarnings(as.integer(value))
        } else {
            defaults[[key]] <- value
        }
        index <- index + 1L
    }
    if (defaults$help) return(defaults)
    if (is.null(defaults$dataset) || !nzchar(defaults$dataset)) {
        stop("--dataset is required.", call. = FALSE)
    }
    for (key in c("reps", "rep_start", "workers")) {
        if (!is.finite(defaults[[key]]) || defaults[[key]] < 1L) {
            stop(sprintf("--%s must be a positive integer.",
                         gsub("_", "-", key)), call. = FALSE)
        }
    }
    defaults
}

print_usage <- function() {
    cat(paste(
        "Usage:",
        "  Rscript run_negative_control.R --dataset=<id-or-template.rds>",
        "      [--reps=100] [--rep-start=1] [--workers=1]",
        "      [--force] [--keep-work]",
        "",
        "Each replicate runs balanced, fourfold, and sixfold depth scenarios.",
        "Existing unit checkpoints are skipped unless --force is supplied.",
        sep = "\n"
    ))
}

resolve_template_path <- function(dataset) {
    if (file.exists(dataset)) return(normalizePath(dataset))
    candidates <- c(
        file.path(analysis_directory, "data", paste0(dataset, ".rds")),
        file.path(
            analysis_directory, "data", paste0(dataset, "_template.rds")
        )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing)) return(normalizePath(existing[[1L]]))

    rds_files <- list.files(
        file.path(analysis_directory, "data"), pattern = "\\.rds$",
        recursive = TRUE, full.names = TRUE
    )
    matches <- vapply(rds_files, function(path) {
        object <- tryCatch(readRDS(path), error = function(error) NULL)
        is.list(object) && identical(as.character(object$dataset_id), dataset)
    }, logical(1))
    if (sum(matches) == 1L) return(normalizePath(rds_files[matches]))
    stop(sprintf("Could not resolve one template for dataset '%s'.", dataset),
         call. = FALSE)
}

validate_template <- function(object, template_path) {
    if (!is.list(object)) {
        stop("Template RDS must contain a named list.", call. = FALSE)
    }
    probability <- object$probability %||% object$probabilities
    evaluation_taxa <- object$evaluation_taxa %||% object$selected_taxa
    dataset_id <- as.character(object$dataset_id %||%
                                   tools::file_path_sans_ext(
                                       basename(template_path)
                                   ))
    dataset_label <- as.character(
        object$label %||% object$dataset_label %||% dataset_id
    )
    if (length(dataset_id) != 1L || !nzchar(dataset_id) ||
        !grepl("^[A-Za-z0-9._-]+$", dataset_id)) {
        stop("dataset_id must be one safe, nonempty identifier.",
             call. = FALSE)
    }
    if (length(dataset_label) != 1L || !nzchar(dataset_label)) {
        stop("Template label must be one nonempty string.", call. = FALSE)
    }
    if (is.data.frame(probability)) probability <- as.matrix(probability)
    if (!is.matrix(probability) || !is.numeric(probability) ||
        is.null(rownames(probability)) || is.null(colnames(probability))) {
        stop("probability must be a named numeric matrix.", call. = FALSE)
    }
    evaluation_taxa <- as.character(evaluation_taxa)
    if (length(evaluation_taxa) != 30L || anyNA(evaluation_taxa) ||
        any(!nzchar(evaluation_taxa)) || anyDuplicated(evaluation_taxa)) {
        stop("evaluation_taxa must contain exactly 30 unique names.",
             call. = FALSE)
    }
    if (!all(evaluation_taxa %in% colnames(probability)) &&
        all(evaluation_taxa %in% rownames(probability))) {
        probability <- t(probability)
    }
    if (!all(evaluation_taxa %in% colnames(probability))) {
        stop("Every evaluation taxon must be a probability column.",
             call. = FALSE)
    }
    if ("Other_unmodeled" %in% evaluation_taxa) {
        stop("Other_unmodeled cannot be an evaluation taxon.", call. = FALSE)
    }
    if (anyNA(probability) || any(!is.finite(probability)) ||
        any(probability < 0)) {
        stop("probability entries must be finite and nonnegative.",
             call. = FALSE)
    }
    if (anyDuplicated(rownames(probability)) ||
        anyDuplicated(colnames(probability))) {
        stop("Probability row and column names must be unique.",
             call. = FALSE)
    }
    row_total <- rowSums(probability)
    if (any(!is.finite(row_total)) || any(row_total <= 0)) {
        stop("Every probability row must have positive mass.", call. = FALSE)
    }
    probability <- sweep(probability, 1L, row_total, "/")
    if (nrow(probability) < 2L || nrow(probability) %% 2L != 0L) {
        stop("The template must contain an even number of at least two samples.",
             call. = FALSE)
    }
    dataset_index <- suppressWarnings(as.integer(object$dataset_index %||% NA))
    if (length(dataset_index) != 1L || !is.finite(dataset_index) ||
        dataset_index < 1L) {
        dataset_index <- NA_integer_
    }
    list(
        probability = probability,
        evaluation_taxa = evaluation_taxa,
        dataset_id = dataset_id,
        label = dataset_label,
        dataset_index = dataset_index,
        template_path = normalizePath(template_path),
        template_md5 = unname(tools::md5sum(template_path))
    )
}

stable_text_hash <- function(text) {
    bytes <- utf8ToInt(enc2utf8(text))
    if (!length(bytes)) return(1L)
    value <- 0
    for (index in seq_along(bytes)) {
        value <- (value * 131 + bytes[[index]] + index) %% 100000L
    }
    as.integer(value)
}

dataset_seed_index <- function(template) {
    if (is.finite(template$dataset_index)) {
        return(as.integer(template$dataset_index %% 100000L))
    }
    stable_text_hash(template$dataset_id)
}

make_seed <- function(template, replicate, stream, scenario_index = 0L) {
    modulus <- 2147483646
    value <- 900000 + dataset_seed_index(template) * 1000003 +
        as.integer(replicate) * 1009 + as.integer(stream) * 97 +
        as.integer(scenario_index) * 100003
    as.integer((value %% modulus) + 1)
}

draw_balanced_group <- function(template, replicate) {
    set.seed(make_seed(template, replicate, stream = 1L))
    labels <- rep(
        c(reference_group, comparison_group),
        each = nrow(template$probability) / 2L
    )
    factor(sample(labels, replace = FALSE),
           levels = c(reference_group, comparison_group))
}

draw_scenario_counts <- function(template, group, scenario, replicate) {
    scenario_row <- scenario_table[
        scenario_table$scenario == scenario, , drop = FALSE
    ]
    if (nrow(scenario_row) != 1L) {
        stop(sprintf("Unknown scenario: %s", scenario), call. = FALSE)
    }
    set.seed(make_seed(
        template, replicate, stream = 2L,
        scenario_index = scenario_row$scenario_index
    ))
    target_median <- ifelse(
        group == reference_group,
        scenario_row$median_H,
        scenario_row$median_Case
    )
    library_size <- round(stats::rlnorm(
        length(group), meanlog = log(target_median), sdlog = depth_sdlog
    ))
    library_size <- as.integer(pmin(
        pmax(library_size, minimum_depth), maximum_depth
    ))

    evaluation_probability <- template$probability[,
        template$evaluation_taxa, drop = FALSE]
    other_probability <- pmax(0, 1 - rowSums(evaluation_probability))
    sampling_probability <- cbind(
        evaluation_probability,
        Other_unmodeled = other_probability
    )
    sampling_probability <- sweep(
        sampling_probability, 1L, rowSums(sampling_probability), "/"
    )
    counts <- vapply(seq_len(nrow(sampling_probability)), function(index) {
        as.integer(stats::rmultinom(
            1L, size = library_size[[index]],
            prob = sampling_probability[index, ]
        ))
    }, integer(ncol(sampling_probability)))
    rownames(counts) <- colnames(sampling_probability)
    colnames(counts) <- rownames(template$probability)
    storage.mode(counts) <- "integer"
    if (!identical(as.integer(colSums(counts)), library_size)) {
        stop("Generated count totals do not equal the drawn library sizes.",
             call. = FALSE)
    }
    metadata <- data.frame(
        group = factor(group, levels = c(reference_group, comparison_group)),
        library_size = library_size,
        stringsAsFactors = FALSE,
        row.names = colnames(counts)
    )
    scaled_log_depth <- as.numeric(scale(log(metadata$library_size)))
    if (any(!is.finite(scaled_log_depth))) {
        stop("Standardized log library size is not finite.", call. = FALSE)
    }
    metadata$log_library_size_z <- scaled_log_depth
    list(counts = counts, metadata = metadata)
}

failure_result <- function(taxa, reason) {
    data.frame(
        taxon = taxa,
        native_p_value = NA_real_,
        available = FALSE,
        reason = reason,
        stringsAsFactors = FALSE
    )
}

sanitize_result <- function(result, taxa, method) {
    if (!is.data.frame(result) || !all(
        c("taxon", "native_p_value", "available", "reason") %in%
            names(result)
    )) {
        result <- failure_result(taxa, "method returned an invalid result table")
    }
    aligned <- result[match(taxa, result$taxon), , drop = FALSE]
    aligned$taxon <- taxa
    missing_row <- is.na(match(taxa, result$taxon))
    aligned$native_p_value <- suppressWarnings(as.numeric(
        aligned$native_p_value
    ))
    aligned$available <- as.logical(aligned$available)
    aligned$available[is.na(aligned$available)] <- FALSE
    invalid <- missing_row | !is.finite(aligned$native_p_value) |
        aligned$native_p_value < 0 | aligned$native_p_value > 1
    aligned$available[invalid] <- FALSE
    aligned$reason <- as.character(aligned$reason)
    aligned$reason[missing_row] <- "taxon not returned"
    aligned$reason[invalid & !missing_row] <- "invalid p-value"
    aligned$reason[aligned$available &
                       (is.na(aligned$reason) | !nzchar(aligned$reason))] <-
        "available"
    aligned$p_value <- ifelse(
        aligned$available, aligned$native_p_value, 1
    )
    aligned$method <- method
    aligned[, c(
        "taxon", "method", "p_value", "native_p_value", "available",
        "reason"
    )]
}

run_with_warnings <- function(expression) {
    warnings <- character()
    value <- tryCatch(
        withCallingHandlers(
            expression,
            warning = function(warning) {
                warnings <<- c(warnings, conditionMessage(warning))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(error) error
    )
    list(value = value, warnings = unique(warnings))
}

run_dasra <- function(counts, metadata, taxa) {
    outcome <- run_with_warnings(DASRA::dasra(
        counts = counts,
        metadata = metadata,
        formula = ~ group,
        group = "group",
        library_size = "library_size",
        taxa_are_rows = TRUE,
        reference = reference_group,
        p_adjust_method = "BH",
        component = "all",
        full_output = FALSE,
        conditional_present_starts = 1L
    ))
    if (inherits(outcome$value, "error")) {
        result <- failure_result(
            taxa, paste("method error:", conditionMessage(outcome$value))
        )
        return(list(result = result, warnings = outcome$warnings))
    }
    fit <- outcome$value
    if (!is.data.frame(fit$results) || !is.data.frame(fit$diagnostics)) {
        result <- failure_result(taxa, "DASRA returned incomplete output")
        return(list(result = result, warnings = outcome$warnings))
    }
    result_index <- match(taxa, fit$results$taxon)
    diagnostic_index <- match(taxa, fit$diagnostics$taxon)
    native_p <- rep(NA_real_, length(taxa))
    available <- rep(FALSE, length(taxa))
    reason <- rep("taxon not returned", length(taxa))
    has_result <- !is.na(result_index)
    has_diagnostic <- !is.na(diagnostic_index)
    if ("p_omnibus" %in% names(fit$results)) {
        native_p[has_result] <- suppressWarnings(as.numeric(
            fit$results$p_omnibus[result_index[has_result]]
        ))
    }
    required_diagnostics <- c(
        "retained", "formed_structural_absence", "formed_relative_abundance"
    )
    if (all(required_diagnostics %in% names(fit$diagnostics))) {
        diagnostic_ok <- rep(FALSE, length(taxa))
        diagnostic_ok[has_diagnostic] <- with(
            fit$diagnostics[diagnostic_index[has_diagnostic], , drop = FALSE],
            as.logical(retained) & as.logical(formed_structural_absence) &
                as.logical(formed_relative_abundance)
        )
        available <- has_result & has_diagnostic & diagnostic_ok &
            is.finite(native_p)
        reason[has_result & has_diagnostic & diagnostic_ok] <- "available"
        failed_diagnostic <- has_diagnostic & !diagnostic_ok
        reason[failed_diagnostic] <-
            "one or both DASRA components were unavailable"
    } else {
        available <- has_result & is.finite(native_p)
        reason[has_result] <- ifelse(
            available[has_result], "available", "invalid omnibus p-value"
        )
    }
    list(
        result = data.frame(
            taxon = taxa, native_p_value = native_p,
            available = available, reason = reason,
            stringsAsFactors = FALSE
        ),
        warnings = outcome$warnings
    )
}

run_zinq <- function(counts, metadata, taxa) {
    count_by_sample <- t(counts[taxa, , drop = FALSE])
    relative_abundance <- sweep(
        count_by_sample, 1L, metadata$library_size, "/"
    )
    group_binary <- as.integer(metadata$group == comparison_group)
    covariate_matrix <- stats::model.matrix(
        ~ log_library_size_z, metadata
    )[, -1L, drop = FALSE]
    colnames(covariate_matrix) <- "covariate1"
    model_formula <- y ~ group + covariate1
    taus <- c(0.25, 0.50, 0.75)
    native_p <- rep(NA_real_, length(taxa))
    available <- rep(FALSE, length(taxa))
    reason <- rep("model did not return a valid p-value", length(taxa))
    warnings <- character()

    for (index in seq_along(taxa)) {
        taxon <- taxa[[index]]
        model_data <- data.frame(
            y = as.numeric(relative_abundance[, taxon]),
            group = group_binary,
            covariate1 = covariate_matrix[, 1L],
            stringsAsFactors = FALSE
        )
        fit_outcome <- run_with_warnings(ZINQ::ZINQ_tests(
            formula.logistic = model_formula,
            formula.quantile = model_formula,
            C = "group",
            y_CorD = "C",
            data = model_data,
            taus = taus,
            seed = 2026L + index
        ))
        warnings <- c(warnings, fit_outcome$warnings)
        if (inherits(fit_outcome$value, "error")) {
            reason[[index]] <- paste(
                "model error:", conditionMessage(fit_outcome$value)
            )
            next
        }
        combination <- run_with_warnings(ZINQ::ZINQ_combination(
            fit_outcome$value, method = "Cauchy", taus = taus
        ))
        warnings <- c(warnings, combination$warnings)
        if (inherits(combination$value, "error")) {
            reason[[index]] <- paste(
                "combination error:", conditionMessage(combination$value)
            )
            next
        }
        value <- suppressWarnings(as.numeric(combination$value))
        if (length(value) == 1L && is.finite(value) &&
            value >= 0 && value <= 1) {
            native_p[[index]] <- value
            available[[index]] <- TRUE
            reason[[index]] <- "available"
        } else {
            reason[[index]] <- "invalid native Cauchy p-value"
        }
    }
    list(
        result = data.frame(
            taxon = taxa, native_p_value = native_p,
            available = available, reason = reason,
            stringsAsFactors = FALSE
        ),
        warnings = unique(warnings)
    )
}

extract_maaslin_joint <- function(fit, taxa) {
    prevalence_table <- fit$fit_data_prevalence$results
    abundance_table <- fit$fit_data_abundance$results
    extract_table <- function(table) {
        if (!is.data.frame(table) || !all(c(
            "feature", "metadata", "value", "pval_individual", "pval_joint"
        ) %in% names(table))) {
            return(list(
                individual = rep(NA_real_, length(taxa)),
                joint = rep(NA_real_, length(taxa)),
                unique = rep(FALSE, length(taxa))
            ))
        }
        table <- table[
            table$metadata == "group" & table$value == comparison_group,
            , drop = FALSE
        ]
        individual <- joint <- rep(NA_real_, length(taxa))
        unique <- rep(FALSE, length(taxa))
        for (index in seq_along(taxa)) {
            rows <- which(table$feature == taxa[[index]])
            if (length(rows) == 1L) {
                individual[[index]] <- suppressWarnings(as.numeric(
                    table$pval_individual[rows]
                ))
                joint[[index]] <- suppressWarnings(as.numeric(
                    table$pval_joint[rows]
                ))
                unique[[index]] <- TRUE
            }
        }
        list(individual = individual, joint = joint, unique = unique)
    }
    prevalence <- extract_table(prevalence_table)
    abundance <- extract_table(abundance_table)
    agreement <- is.finite(prevalence$joint) & is.finite(abundance$joint) &
        abs(prevalence$joint - abundance$joint) <= 1e-12
    available <- prevalence$unique & abundance$unique &
        is.finite(prevalence$individual) & is.finite(abundance$individual) &
        agreement
    reason <- ifelse(
        available,
        "available",
        "both component models and a consistent joint p-value were required"
    )
    data.frame(
        taxon = taxa,
        native_p_value = prevalence$joint,
        available = available,
        reason = reason,
        stringsAsFactors = FALSE
    )
}

run_maaslin3 <- function(counts, metadata, taxa, work_directory,
                         keep_work = FALSE) {
    dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)
    output_directory <- file.path(work_directory, "maaslin3")
    if (dir.exists(output_directory)) {
        unlink(output_directory, recursive = TRUE, force = TRUE)
    }
    on.exit(try(maaslin3::maaslin_log_reset(), silent = TRUE), add = TRUE)
    if (!keep_work) {
        on.exit(unlink(work_directory, recursive = TRUE, force = TRUE),
                add = TRUE)
    }
    outcome <- run_with_warnings(maaslin3::maaslin3(
        input_data = as.data.frame(t(counts), check.names = FALSE),
        input_metadata = metadata,
        output = output_directory,
        formula = ~ log_library_size_z + group,
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
        reference = paste0("group,", reference_group)
    ))
    try(maaslin3::maaslin_log_reset(), silent = TRUE)
    if (inherits(outcome$value, "error")) {
        result <- failure_result(
            taxa, paste("method error:", conditionMessage(outcome$value))
        )
    } else {
        result <- tryCatch(
            extract_maaslin_joint(outcome$value, taxa),
            error = function(error) failure_result(
                taxa, paste("result extraction error:", conditionMessage(error))
            )
        )
    }
    list(result = result, warnings = outcome$warnings)
}

run_timed_method <- function(method, runner, taxa) {
    started <- proc.time()[["elapsed"]]
    outcome <- tryCatch(
        runner(),
        error = function(error) list(
            result = failure_result(
                taxa, paste("method error:", conditionMessage(error))
            ),
            warnings = character()
        )
    )
    elapsed <- proc.time()[["elapsed"]] - started
    list(
        result = sanitize_result(outcome$result, taxa, method),
        elapsed_seconds = elapsed,
        warnings = paste(unique(outcome$warnings), collapse = " | ")
    )
}

make_taxon_diagnostics <- function(counts, metadata, taxa) {
    tested <- counts[taxa, , drop = FALSE]
    is_H <- metadata$group == reference_group
    is_Case <- metadata$group == comparison_group
    data.frame(
        taxon = taxa,
        observed_prevalence = rowMeans(tested > 0),
        prevalence_H = rowMeans(tested[, is_H, drop = FALSE] > 0),
        prevalence_Case = rowMeans(tested[, is_Case, drop = FALSE] > 0),
        positive_samples_H = rowSums(tested[, is_H, drop = FALSE] > 0),
        positive_samples_Case = rowSums(
            tested[, is_Case, drop = FALSE] > 0
        ),
        total_count = rowSums(tested),
        stringsAsFactors = FALSE
    )
}

unit_checkpoint_path <- function(result_directory, scenario, replicate) {
    file.path(
        result_directory, "checkpoints", scenario,
        sprintf("rep_%04d.rds", replicate)
    )
}

atomic_save_rds <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
        pattern = paste0(basename(path), "."), tmpdir = dirname(path)
    )
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(object, temporary, version = 3)
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) {
            stop(sprintf("Could not write checkpoint: %s", path),
                 call. = FALSE)
        }
        unlink(temporary)
    }
    invisible(path)
}

run_unit <- function(template, group, scenario, replicate, result_directory,
                     work_root, keep_work = FALSE) {
    unit_started <- proc.time()[["elapsed"]]
    generated <- draw_scenario_counts(template, group, scenario, replicate)
    counts <- generated$counts
    metadata <- generated$metadata
    taxa <- template$evaluation_taxa
    work_directory <- file.path(
        work_root, template$dataset_id, scenario,
        sprintf("rep_%04d_pid_%d", replicate, Sys.getpid())
    )

    dasra <- run_timed_method(
        "DASRA",
        function() run_dasra(counts, metadata, taxa),
        taxa
    )
    zinq <- run_timed_method(
        "ZINQ Cauchy",
        function() run_zinq(counts, metadata, taxa),
        taxa
    )
    maaslin <- run_timed_method(
        "MaAsLin3 joint",
        function() run_maaslin3(
            counts, metadata, taxa, work_directory, keep_work
        ),
        taxa
    )
    method_results <- list(dasra, zinq, maaslin)
    taxon_pvalues <- do.call(rbind, lapply(method_results, `[[`, "result"))
    taxon_pvalues$dataset <- template$dataset_id
    taxon_pvalues$dataset_label <- template$label
    taxon_pvalues$replicate <- as.integer(replicate)
    taxon_pvalues$scenario <- scenario
    taxon_pvalues$status <- ifelse(
        taxon_pvalues$available, "available", "unavailable"
    )
    taxon_pvalues <- taxon_pvalues[, c(
        "dataset", "dataset_label", "replicate", "scenario", "method",
        "taxon", "p_value", "native_p_value", "available", "status",
        "reason"
    )]

    replicate_metrics <- do.call(rbind, lapply(method_results, function(item) {
        result <- item$result
        data.frame(
            dataset = template$dataset_id,
            dataset_label = template$label,
            replicate = as.integer(replicate),
            scenario = scenario,
            method = unique(result$method),
            n_taxa = length(taxa),
            n_available = sum(result$available),
            availability_rate = mean(result$available),
            n_rejections = sum(result$p_value <= 0.05),
            type1_error = mean(result$p_value <= 0.05),
            stringsAsFactors = FALSE
        )
    }))
    timing <- do.call(rbind, Map(function(item, method) {
        data.frame(
            dataset = template$dataset_id,
            dataset_label = template$label,
            replicate = as.integer(replicate),
            scenario = scenario,
            method = method,
            elapsed_seconds = item$elapsed_seconds,
            warnings = item$warnings,
            stringsAsFactors = FALSE
        )
    }, method_results, method_names))

    diagnostics <- make_taxon_diagnostics(counts, metadata, taxa)
    diagnostics$dataset <- template$dataset_id
    diagnostics$dataset_label <- template$label
    diagnostics$replicate <- as.integer(replicate)
    diagnostics$scenario <- scenario
    diagnostics <- diagnostics[, c(
        "dataset", "dataset_label", "replicate", "scenario", "taxon",
        "observed_prevalence", "prevalence_H", "prevalence_Case",
        "positive_samples_H", "positive_samples_Case", "total_count"
    )]
    depth_diagnostics <- data.frame(
        dataset = template$dataset_id,
        dataset_label = template$label,
        replicate = as.integer(replicate),
        scenario = scenario,
        n_H = sum(metadata$group == reference_group),
        n_Case = sum(metadata$group == comparison_group),
        median_depth_H = stats::median(
            metadata$library_size[metadata$group == reference_group]
        ),
        median_depth_Case = stats::median(
            metadata$library_size[metadata$group == comparison_group]
        ),
        depth_ratio_Case_to_H = stats::median(
            metadata$library_size[metadata$group == comparison_group]
        ) / stats::median(
            metadata$library_size[metadata$group == reference_group]
        ),
        minimum_depth = min(metadata$library_size),
        maximum_depth = max(metadata$library_size),
        group_seed = make_seed(template, replicate, stream = 1L),
        count_seed = make_seed(
            template, replicate, stream = 2L,
            scenario_index = scenario_table$scenario_index[
                match(scenario, scenario_table$scenario)
            ]
        ),
        elapsed_seconds = proc.time()[["elapsed"]] - unit_started,
        stringsAsFactors = FALSE
    )
    checkpoint <- list(
        schema_version = 1L,
        template = list(
            dataset_id = template$dataset_id,
            label = template$label,
            template_path = template$template_path,
            template_md5 = template$template_md5,
            n_samples = nrow(template$probability),
            evaluation_taxa = taxa
        ),
        taxon_pvalues = taxon_pvalues,
        replicate_metrics = replicate_metrics,
        timing = timing,
        diagnostics = diagnostics,
        depth_diagnostics = depth_diagnostics
    )
    checkpoint_path <- unit_checkpoint_path(
        result_directory, scenario, replicate
    )
    atomic_save_rds(checkpoint, checkpoint_path)
    message(sprintf(
        "%s replicate %d: %s complete (%.1f s)",
        template$label, replicate, scenario,
        depth_diagnostics$elapsed_seconds
    ))
    checkpoint_path
}

read_valid_checkpoints <- function(result_directory, template) {
    paths <- list.files(
        file.path(result_directory, "checkpoints"),
        pattern = "^rep_[0-9]+\\.rds$", recursive = TRUE,
        full.names = TRUE
    )
    objects <- lapply(paths, function(path) {
        object <- tryCatch(readRDS(path), error = function(error) NULL)
        if (is.null(object) || !identical(
            object$template$template_md5, template$template_md5
        )) {
            return(NULL)
        }
        object
    })
    Filter(Negate(is.null), objects)
}

checkpoint_is_current <- function(path, template) {
    if (!file.exists(path)) return(FALSE)
    object <- tryCatch(readRDS(path), error = function(error) NULL)
    !is.null(object) && identical(
        object$template$template_md5, template$template_md5
    )
}

atomic_write_csv <- function(table, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
        pattern = paste0(basename(path), "."), tmpdir = dirname(path)
    )
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    utils::write.csv(table, temporary, row.names = FALSE, na = "")
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) {
            stop(sprintf("Could not write output: %s", path), call. = FALSE)
        }
        unlink(temporary)
    }
    invisible(path)
}

aggregate_checkpoints <- function(result_directory, template) {
    checkpoints <- read_valid_checkpoints(result_directory, template)
    if (!length(checkpoints)) return(invisible(NULL))
    checkpoint_depth <- do.call(rbind, lapply(
        checkpoints, `[[`, "depth_diagnostics"
    ))
    units_by_replicate <- table(
        checkpoint_depth$replicate, checkpoint_depth$scenario
    )
    completed_replicates <- sum(
        rowSums(units_by_replicate > 0L) == nrow(scenario_table)
    )
    components <- c(
        "taxon_pvalues", "replicate_metrics", "timing", "diagnostics",
        "depth_diagnostics"
    )
    output_names <- c(
        taxon_pvalues = "taxon_pvalues.csv",
        replicate_metrics = "replicate_metrics.csv",
        timing = "timing.csv",
        diagnostics = "diagnostics.csv",
        depth_diagnostics = "depth_diagnostics.csv"
    )
    for (component in components) {
        table <- do.call(rbind, lapply(checkpoints, `[[`, component))
        ordering <- intersect(
            c("replicate", "scenario", "method", "taxon"), names(table)
        )
        if (length(ordering)) {
            order_arguments <- Map(function(values, variable) {
                if (identical(variable, "scenario")) {
                    match(values, scenario_table$scenario)
                } else if (identical(variable, "method")) {
                    match(values, method_names)
                } else {
                    values
                }
            }, table[ordering], ordering)
            table <- table[
                do.call(order, unname(order_arguments)), , drop = FALSE
            ]
        }
        atomic_write_csv(table, file.path(
            result_directory, output_names[[component]]
        ))
    }
    manifest <- data.frame(
        dataset = template$dataset_id,
        dataset_label = template$label,
        template_path = template$template_path,
        template_md5 = template$template_md5,
        n_samples = nrow(template$probability),
        group_size = nrow(template$probability) / 2L,
        n_taxa = length(template$evaluation_taxa),
        completed_units = length(checkpoints),
        replicates_with_any_checkpoint = nrow(units_by_replicate),
        completed_replicates = completed_replicates,
        depth_sdlog = depth_sdlog,
        minimum_depth = minimum_depth,
        maximum_depth = maximum_depth,
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        DASRA_version = as.character(utils::packageVersion("DASRA")),
        ZINQ_version = as.character(utils::packageVersion("ZINQ")),
        maaslin3_version = as.character(utils::packageVersion("maaslin3")),
        stringsAsFactors = FALSE
    )
    atomic_write_csv(manifest, file.path(result_directory, "run_manifest.csv"))
    invisible(manifest)
}

main <- function() {
    arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
    if (arguments$help) {
        print_usage()
        return(invisible(NULL))
    }
    required_packages <- c("DASRA", "ZINQ", "maaslin3")
    missing_packages <- required_packages[!vapply(
        required_packages, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing_packages)) {
        stop(sprintf(
            "Required R packages are unavailable: %s",
            paste(missing_packages, collapse = ", ")
        ), call. = FALSE)
    }

    template_path <- resolve_template_path(arguments$dataset)
    template <- validate_template(readRDS(template_path), template_path)
    result_directory <- file.path(
        analysis_directory, "results", template$dataset_id
    )
    work_root <- file.path(analysis_directory, "work")
    dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

    replicate_ids <- seq.int(
        arguments$rep_start,
        length.out = arguments$reps
    )
    groups <- stats::setNames(lapply(
        replicate_ids,
        function(replicate) draw_balanced_group(template, replicate)
    ), as.character(replicate_ids))
    units <- expand.grid(
        replicate = replicate_ids,
        scenario = scenario_table$scenario,
        stringsAsFactors = FALSE
    )
    units$scenario <- factor(
        units$scenario, levels = scenario_table$scenario
    )
    units <- units[order(units$replicate, units$scenario), , drop = FALSE]
    units$scenario <- as.character(units$scenario)
    checkpoint_paths <- mapply(
        function(scenario, replicate) unit_checkpoint_path(
            result_directory, scenario, replicate
        ),
        units$scenario, units$replicate,
        USE.NAMES = FALSE
    )
    if (!arguments$force) {
        current <- vapply(
            checkpoint_paths, checkpoint_is_current, logical(1),
            template = template
        )
        units <- units[!current, , drop = FALSE]
    }

    message(sprintf(
        "Dataset %s: %d samples (%d per arm), fixed family of %d taxa.",
        template$label, nrow(template$probability),
        nrow(template$probability) / 2L,
        length(template$evaluation_taxa)
    ))
    message(sprintf(
        "Requested replicates %d-%d; %d scenario units remain.",
        min(replicate_ids), max(replicate_ids), nrow(units)
    ))

    run_one <- function(index) {
        unit <- units[index, , drop = FALSE]
        run_unit(
            template = template,
            group = groups[[as.character(unit$replicate)]],
            scenario = unit$scenario,
            replicate = unit$replicate,
            result_directory = result_directory,
            work_root = work_root,
            keep_work = arguments$keep_work
        )
    }
    if (nrow(units)) {
        if (arguments$workers > 1L && .Platform$OS.type != "windows") {
            outcomes <- parallel::mclapply(
                seq_len(nrow(units)), run_one,
                mc.cores = min(arguments$workers, nrow(units)),
                mc.preschedule = FALSE, mc.set.seed = FALSE
            )
            failed <- vapply(outcomes, inherits, logical(1), "try-error")
            if (any(failed)) {
                stop(sprintf("%d scenario units failed before checkpointing.",
                             sum(failed)), call. = FALSE)
            }
        } else {
            for (index in seq_len(nrow(units))) {
                run_one(index)
                aggregate_checkpoints(result_directory, template)
            }
        }
    }
    manifest <- aggregate_checkpoints(result_directory, template)
    if (!is.null(manifest)) {
        message(sprintf(
            "%s: %d scenario units are checkpointed in %s",
            template$label, manifest$completed_units,
            normalizePath(result_directory)
        ))
    }
    invisible(manifest)
}

if (sys.nframe() == 0L) main()
