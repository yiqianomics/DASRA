#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

scenario_levels <- c(
  "Balanced depth",
  "Fourfold depth difference",
  "Sixfold depth difference"
)
method_levels <- c("DASRA", "ZINQ Cauchy", "MaAsLin3 joint")

usage <- function() {
  cat(paste0(
    "Summarize negative-control randomizations and draw the Type I error figure.\n\n",
    "Usage:\n",
    "  Rscript summarize_negative_control.R [options]\n\n",
    "Options:\n",
    "  --input PATH              Result file or directory; may be repeated.\n",
    "  --results-dir PATH        Default input directory.\n",
    "  --output-dir PATH         Directory for summary tables and figure.\n",
    "  --manifest PATH           Optional dataset manifest with labels/order.\n",
    "  --pdf PATH                Output PDF path.\n",
    "  --expected-reps N         Expected randomizations per cell [100].\n",
    "  --expected-datasets N     Expected number of datasets [10].\n",
    "  --expected-taxa N         Fixed taxon denominator [30].\n",
    "  --alpha X                 Raw-P rejection threshold [0.05].\n",
    "  --allow-partial           Summarize incomplete runs without failing.\n",
    "  --validate-only           Validate and write tables, but skip the PDF.\n",
    "  --help                    Show this help.\n"
  ))
}

script_directory <- function() {
  command <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command, value = TRUE)
  if (!length(file_arg)) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  normalizePath(dirname(sub("^--file=", "", file_arg[[1L]])), mustWork = TRUE)
}

parse_arguments <- function(args, script_dir) {
  defaults <- list(
    input = character(),
    results_dir = file.path(script_dir, "results"),
    output_dir = file.path(script_dir, "results", "summary"),
    manifest = NA_character_,
    pdf = NA_character_,
    expected_reps = 100L,
    expected_datasets = 10L,
    expected_taxa = 30L,
    alpha = 0.05,
    allow_partial = FALSE,
    validate_only = FALSE,
    help = FALSE
  )
  value_options <- c(
    input = "input",
    `results-dir` = "results_dir",
    `output-dir` = "output_dir",
    manifest = "manifest",
    pdf = "pdf",
    `expected-reps` = "expected_reps",
    `expected-datasets` = "expected_datasets",
    `expected-taxa` = "expected_taxa",
    alpha = "alpha"
  )
  flag_options <- c(
    `allow-partial` = "allow_partial",
    `validate-only` = "validate_only",
    help = "help"
  )

  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
    stripped <- substring(token, 3L)
    if (grepl("=", stripped, fixed = TRUE)) {
      pieces <- strsplit(stripped, "=", fixed = TRUE)[[1L]]
      name <- pieces[[1L]]
      value <- paste(pieces[-1L], collapse = "=")
    } else {
      name <- stripped
      value <- NULL
    }
    if (name %in% names(flag_options)) {
      if (!is.null(value)) {
        stop("Flag --", name, " does not take a value.", call. = FALSE)
      }
      defaults[[flag_options[[name]]]] <- TRUE
    } else if (name %in% names(value_options)) {
      if (is.null(value)) {
        i <- i + 1L
        if (i > length(args)) {
          stop("Missing value for --", name, ".", call. = FALSE)
        }
        value <- args[[i]]
      }
      key <- value_options[[name]]
      if (identical(key, "input")) {
        defaults$input <- c(defaults$input, value)
      } else {
        defaults[[key]] <- value
      }
    } else {
      stop("Unknown option: --", name, call. = FALSE)
    }
    i <- i + 1L
  }

  defaults$expected_reps <- suppressWarnings(as.integer(defaults$expected_reps))
  defaults$expected_datasets <- suppressWarnings(as.integer(defaults$expected_datasets))
  defaults$expected_taxa <- suppressWarnings(as.integer(defaults$expected_taxa))
  defaults$alpha <- suppressWarnings(as.numeric(defaults$alpha))
  if (any(is.na(c(
    defaults$expected_reps,
    defaults$expected_datasets,
    defaults$expected_taxa,
    defaults$alpha
  )))) {
    stop("Expected counts and alpha must be numeric.", call. = FALSE)
  }
  if (any(c(
    defaults$expected_reps,
    defaults$expected_datasets,
    defaults$expected_taxa
  ) <= 0L)) {
    stop("Expected counts must be positive integers.", call. = FALSE)
  }
  if (defaults$alpha <= 0 || defaults$alpha >= 1) {
    stop("--alpha must be strictly between zero and one.", call. = FALSE)
  }
  if (!length(defaults$input)) {
    defaults$input <- defaults$results_dir
  }
  defaults
}

clean_name <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

pick_column <- function(data, aliases) {
  cleaned <- clean_name(names(data))
  hit <- match(aliases, cleaned, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) NA_integer_ else hit[[1L]]
}

column_value <- function(data, aliases, default = NA) {
  index <- pick_column(data, aliases)
  if (!is.na(index)) return(data[[index]])
  if (length(default) == nrow(data)) default else rep_len(default, nrow(data))
}

canonical_scenario <- function(x) {
  key <- clean_name(as.character(x))
  answer <- rep(NA_character_, length(key))
  answer[key %in% c(
    "balanced", "balanced_depth", "equal_depth", "equal", "1x", "onefold"
  )] <- scenario_levels[[1L]]
  answer[key %in% c(
    "fourfold", "four_fold", "fourfold_depth", "fourfold_depth_difference",
    "4fold", "4_fold", "4x", "1500_6000"
  )] <- scenario_levels[[2L]]
  answer[key %in% c(
    "sixfold", "six_fold", "sixfold_depth", "sixfold_depth_difference",
    "6fold", "6_fold", "6x", "1000_6000"
  )] <- scenario_levels[[3L]]
  answer
}

canonical_method <- function(x) {
  key <- clean_name(as.character(x))
  answer <- rep(NA_character_, length(key))
  answer[key %in% c(
    "dasra", "dadra", "dasra_combined", "dasra_omnibus", "dasra_all"
  )] <- method_levels[[1L]]
  answer[key %in% c(
    "zinq", "zinq_cauchy", "cauchy", "zinq_combined"
  )] <- method_levels[[2L]]
  answer[key %in% c(
    "maaslin3", "maaslin_3", "maaslin3_joint", "maaslin_3_joint"
  )] <- method_levels[[3L]]
  answer
}

infer_method <- function(source) {
  key <- clean_name(basename(source))
  if (grepl("dasra|dadra", key)) return(method_levels[[1L]])
  if (grepl("zinq", key)) return(method_levels[[2L]])
  if (grepl("maaslin", key)) return(method_levels[[3L]])
  NA_character_
}

logical_value <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
  key <- clean_name(as.character(x))
  answer <- rep(NA, length(key))
  answer[key %in% c("true", "t", "yes", "y", "1", "available", "success", "ok")] <- TRUE
  answer[key %in% c("false", "f", "no", "n", "0", "unavailable", "failed", "error")] <- FALSE
  answer
}

method_component_is_primary <- function(method, component) {
  key <- clean_name(as.character(component))
  empty <- is.na(component) | key == ""
  keep <- empty
  keep[method == method_levels[[1L]]] <- empty[method == method_levels[[1L]]] |
    key[method == method_levels[[1L]]] %in% c("all", "combined", "omnibus", "p_omnibus")
  keep[method == method_levels[[2L]]] <- empty[method == method_levels[[2L]]] |
    key[method == method_levels[[2L]]] %in% c("cauchy", "combined", "omnibus", "p_cauchy")
  keep[method == method_levels[[3L]]] <- empty[method == method_levels[[3L]]] |
    key[method == method_levels[[3L]]] %in% c("joint", "combined", "omnibus", "p_joint")
  keep
}

normalize_long_frame <- function(data, source, object_name = "") {
  if (!is.data.frame(data) || !nrow(data)) return(NULL)
  dataset_index <- pick_column(data, c(
    "dataset", "dataset_id", "template", "template_id", "study", "study_id", "study_name"
  ))
  replicate_index <- pick_column(data, c(
    "replicate", "rep", "replicate_id", "randomization", "randomisation", "iteration", "sim"
  ))
  scenario_index <- pick_column(data, c(
    "scenario", "setting", "depth_setting", "depth_scenario", "design"
  ))
  taxon_index <- pick_column(data, c(
    "taxon", "taxon_id", "feature", "feature_id", "species"
  ))
  required_indices <- c(dataset_index, replicate_index, scenario_index, taxon_index)
  if (any(is.na(required_indices))) return(NULL)

  method_index <- pick_column(data, c("method", "method_name", "test_method", "procedure"))
  p_index <- pick_column(data, c(
    "p_value", "pvalue", "raw_p", "p_raw", "primary_p", "combined_p", "p", "pval"
  ))
  component_index <- pick_column(data, c("component", "test", "test_name", "contrast_type"))

  if (!is.na(method_index) && is.na(p_index)) {
    method_temp <- canonical_method(data[[method_index]])
    p_candidates <- list(
      DASRA = column_value(data, c("p_omnibus", "omnibus_p", "dasra_p")),
      `ZINQ Cauchy` = column_value(data, c("p_cauchy", "cauchy_p", "zinq_p")),
      `MaAsLin3 joint` = column_value(data, c("p_joint", "joint_p", "maaslin3_p"))
    )
    p_value <- rep(NA_real_, nrow(data))
    for (name in names(p_candidates)) {
      use <- method_temp == name
      p_value[use] <- suppressWarnings(as.numeric(p_candidates[[name]][use]))
    }
  } else if (!is.na(p_index)) {
    p_value <- suppressWarnings(as.numeric(data[[p_index]]))
  } else {
    p_value <- NULL
  }

  if (!is.na(method_index) && !is.null(p_value)) {
    method <- canonical_method(data[[method_index]])
    component <- if (is.na(component_index)) rep(NA_character_, nrow(data)) else data[[component_index]]
    keep_primary <- method_component_is_primary(method, component)
    long <- data.frame(
      dataset = as.character(data[[dataset_index]]),
      dataset_label = as.character(column_value(data, c(
        "dataset_label", "data_label", "display_name", "study_label"
      ), default = NA_character_)),
      dataset_order = suppressWarnings(as.integer(column_value(data, c(
        "dataset_order", "data_order", "display_order", "study_order"
      ), default = NA_integer_))),
      replicate = as.character(data[[replicate_index]]),
      scenario = canonical_scenario(data[[scenario_index]]),
      method = method,
      taxon = as.character(data[[taxon_index]]),
      p_value = p_value,
      native_p_value = suppressWarnings(as.numeric(column_value(data, c(
        "native_p_value", "native_p", "unmodified_p_value"
      ), default = p_value))),
      available_explicit = logical_value(column_value(data, c(
        "available", "is_available", "estimable", "is_estimable", "success"
      ), default = NA)),
      status = as.character(column_value(data, c(
        "status", "fit_status"
      ), default = NA_character_)),
      reason = as.character(column_value(data, c(
        "reason", "error", "message"
      ), default = NA_character_)),
      source = source,
      object = object_name,
      stringsAsFactors = FALSE
    )
    return(long[keep_primary & !is.na(long$method) & !is.na(long$scenario), , drop = FALSE])
  }

  inferred <- infer_method(source)
  if (!is.na(inferred) && !is.null(p_value)) {
    long <- data.frame(
      dataset = as.character(data[[dataset_index]]),
      dataset_label = as.character(column_value(data, c(
        "dataset_label", "data_label", "display_name", "study_label"
      ), default = NA_character_)),
      dataset_order = suppressWarnings(as.integer(column_value(data, c(
        "dataset_order", "data_order", "display_order", "study_order"
      ), default = NA_integer_))),
      replicate = as.character(data[[replicate_index]]),
      scenario = canonical_scenario(data[[scenario_index]]),
      method = inferred,
      taxon = as.character(data[[taxon_index]]),
      p_value = p_value,
      native_p_value = suppressWarnings(as.numeric(column_value(data, c(
        "native_p_value", "native_p", "unmodified_p_value"
      ), default = p_value))),
      available_explicit = logical_value(column_value(data, c(
        "available", "is_available", "estimable", "is_estimable", "success"
      ), default = NA)),
      status = as.character(column_value(data, c(
        "status", "fit_status"
      ), default = NA_character_)),
      reason = as.character(column_value(data, c(
        "reason", "error", "message"
      ), default = NA_character_)),
      source = source,
      object = object_name,
      stringsAsFactors = FALSE
    )
    return(long[!is.na(long$scenario), , drop = FALSE])
  }

  wide_p <- list(
    DASRA = pick_column(data, c("dasra_p", "p_dasra", "dasra_p_value", "p_omnibus", "omnibus_p")),
    `ZINQ Cauchy` = pick_column(data, c("zinq_p", "p_zinq", "zinq_p_value", "p_cauchy", "cauchy_p")),
    `MaAsLin3 joint` = pick_column(data, c(
      "maaslin3_p", "p_maaslin3", "maaslin3_p_value", "p_joint", "joint_p"
    ))
  )
  wide_p <- wide_p[!vapply(wide_p, is.na, logical(1L))]
  if (!length(wide_p)) return(NULL)
  pieces <- lapply(names(wide_p), function(method_name) {
    p_col <- wide_p[[method_name]]
    availability_alias <- switch(
      method_name,
      DASRA = c("dasra_available", "available_dasra"),
      `ZINQ Cauchy` = c("zinq_available", "available_zinq"),
      `MaAsLin3 joint` = c("maaslin3_available", "available_maaslin3")
    )
    data.frame(
      dataset = as.character(data[[dataset_index]]),
      dataset_label = as.character(column_value(data, c(
        "dataset_label", "data_label", "display_name", "study_label"
      ), default = NA_character_)),
      dataset_order = suppressWarnings(as.integer(column_value(data, c(
        "dataset_order", "data_order", "display_order", "study_order"
      ), default = NA_integer_))),
      replicate = as.character(data[[replicate_index]]),
      scenario = canonical_scenario(data[[scenario_index]]),
      method = method_name,
      taxon = as.character(data[[taxon_index]]),
      p_value = suppressWarnings(as.numeric(data[[p_col]])),
      native_p_value = suppressWarnings(as.numeric(data[[p_col]])),
      available_explicit = logical_value(column_value(data, availability_alias, default = NA)),
      status = NA_character_,
      reason = NA_character_,
      source = source,
      object = object_name,
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, pieces)
  long[!is.na(long$scenario), , drop = FALSE]
}

frames_from_object <- function(object, source, prefix = "root") {
  if (is.data.frame(object)) {
    normalized <- normalize_long_frame(object, source, prefix)
    if (is.null(normalized) || !nrow(normalized)) list() else list(normalized)
  } else if (is.list(object)) {
    names_object <- names(object)
    if (is.null(names_object)) names_object <- as.character(seq_along(object))
    output <- list()
    for (i in seq_along(object)) {
      child_prefix <- paste(prefix, names_object[[i]], sep = "/")
      output <- c(output, frames_from_object(object[[i]], source, child_prefix))
    }
    output
  } else {
    list()
  }
}

list_input_files <- function(paths, output_dir) {
  output_norm <- normalizePath(output_dir, mustWork = FALSE)
  generated_names <- c(
    "negative_control_taxon_pvalues.csv",
    "negative_control_replicate_type1.csv",
    "negative_control_dataset_summary.csv",
    "negative_control_availability_summary.csv",
    "negative_control_completeness.csv",
    "negative_control_type1.pdf"
  )
  files <- character()
  for (path in paths) {
    if (dir.exists(path)) {
      found <- list.files(
        path,
        pattern = "\\.(csv|csv\\.gz|tsv|tsv\\.gz|rds)$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
      files <- c(files, found)
    } else if (file.exists(path)) {
      files <- c(files, path)
    } else {
      warning("Input does not exist and was skipped: ", path, call. = FALSE)
    }
  }
  files <- unique(normalizePath(files, mustWork = TRUE))
  in_output <- startsWith(files, paste0(output_norm, .Platform$file.sep))
  files <- files[!in_output & !basename(files) %in% generated_names]
  checkpoint_marker <- paste0(.Platform$file.sep, "checkpoints", .Platform$file.sep)
  checkpoint_files <- files[grepl(checkpoint_marker, files, fixed = TRUE)]
  if (length(checkpoint_files)) {
    remove_checkpoint <- vapply(checkpoint_files, function(checkpoint) {
      dataset_dir <- sub(
        paste0(checkpoint_marker, ".*$"),
        "",
        checkpoint
      )
      aggregate_path <- file.path(dataset_dir, "taxon_pvalues.csv")
      file.exists(aggregate_path) &&
        file.info(aggregate_path)$mtime >= file.info(checkpoint)$mtime
    }, logical(1L))
    files <- setdiff(files, checkpoint_files[remove_checkpoint])
  }
  sort(files)
}

read_result_file <- function(path) {
  lower <- tolower(path)
  object <- tryCatch(
    {
      if (grepl("\\.rds$", lower)) {
        readRDS(path)
      } else if (grepl("\\.tsv(\\.gz)?$", lower)) {
        read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
      } else {
        read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      }
    },
    error = function(e) {
      warning("Could not read ", path, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(object)) list() else frames_from_object(object, path)
}

read_all_results <- function(files) {
  frames <- unlist(lapply(files, read_result_file), recursive = FALSE)
  if (!length(frames)) {
    stop(
      "No taxon-level P-value table was found. Required identifiers are dataset, replicate, scenario, and taxon, plus method and P value.",
      call. = FALSE
    )
  }
  data <- do.call(rbind, frames)
  data$dataset <- trimws(data$dataset)
  data$replicate <- trimws(data$replicate)
  data$taxon <- trimws(data$taxon)
  valid_id <- !is.na(data$dataset) & !is.na(data$replicate) & !is.na(data$taxon) &
    nzchar(data$dataset) & nzchar(data$replicate) & nzchar(data$taxon)
  data <- data[valid_id, , drop = FALSE]
  if (!nrow(data)) stop("All result rows had missing identifiers.", call. = FALSE)
  data
}

deduplicate_results <- function(data) {
  key <- paste(
    data$dataset,
    data$replicate,
    data$scenario,
    data$method,
    data$taxon,
    sep = "\r"
  )
  duplicate <- duplicated(key)
  if (any(duplicate)) {
    first <- match(key, key)
    p_first <- data$p_value[first]
    p_same <- (is.na(data$p_value) & is.na(p_first)) |
      (is.finite(data$p_value) & is.finite(p_first) & abs(data$p_value - p_first) <= 1e-12)
    available_first <- data$available_explicit[first]
    available_same <- (is.na(data$available_explicit) & is.na(available_first)) |
      (!is.na(data$available_explicit) & !is.na(available_first) &
        data$available_explicit == available_first)
    conflict <- duplicate & (!p_same | !available_same)
    if (any(conflict)) {
      index <- which(conflict)[[1L]]
      stop(
        "Conflicting duplicate records for dataset/replicate/scenario/method/taxon key: ",
        paste(data[index, c("dataset", "replicate", "scenario", "method", "taxon")], collapse = " | "),
        call. = FALSE
      )
    }
    data <- data[!duplicate, , drop = FALSE]
  }
  rownames(data) <- NULL
  data
}

find_manifest <- function(path, script_dir) {
  if (!is.na(path)) return(path)
  candidates <- c(
    file.path(script_dir, "data", "public_template_manifest.csv"),
    file.path(script_dir, "data", "dataset_manifest.csv"),
    file.path(script_dir, "data", "template_manifest.csv"),
    file.path(script_dir, "dataset_manifest.csv"),
    file.path(script_dir, "template_manifest.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1L]] else NA_character_
}

apply_dataset_metadata <- function(data, manifest_path) {
  metadata <- unique(data[c("dataset", "dataset_label", "dataset_order")])
  metadata$dataset_label[is.na(metadata$dataset_label) | !nzchar(metadata$dataset_label)] <-
    metadata$dataset[is.na(metadata$dataset_label) | !nzchar(metadata$dataset_label)]
  metadata <- metadata[order(metadata$dataset, is.na(metadata$dataset_order)), , drop = FALSE]
  metadata <- metadata[!duplicated(metadata$dataset), , drop = FALSE]

  if (!is.na(manifest_path)) {
    if (!file.exists(manifest_path)) {
      stop("Manifest does not exist: ", manifest_path, call. = FALSE)
    }
    manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
    id_index <- pick_column(manifest, c(
      "dataset", "dataset_id", "template", "template_id", "study", "study_id", "study_name"
    ))
    if (is.na(id_index)) stop("Manifest has no dataset identifier column.", call. = FALSE)
    manifest_dataset <- as.character(manifest[[id_index]])
    manifest_label <- as.character(column_value(manifest, c(
      "dataset_label", "data_label", "display_name", "study_label", "label"
    ), default = manifest_dataset))
    manifest_order <- suppressWarnings(as.integer(column_value(manifest, c(
      "dataset_order", "data_order", "display_order", "study_order", "order"
    ), default = seq_len(nrow(manifest)))))
    lookup_label <- setNames(manifest_label, manifest_dataset)
    lookup_order <- setNames(manifest_order, manifest_dataset)
    matched <- metadata$dataset %in% names(lookup_label)
    metadata$dataset_label[matched] <- lookup_label[metadata$dataset[matched]]
    metadata$dataset_order[matched] <- lookup_order[metadata$dataset[matched]]
  }

  missing_order <- is.na(metadata$dataset_order)
  start <- if (all(missing_order)) 0L else max(metadata$dataset_order[!missing_order])
  metadata$dataset_order[missing_order] <- start + seq_len(sum(missing_order))
  default_label <- metadata$dataset_label == metadata$dataset
  cleaned_label <- gsub("_([0-9]{4})_([A-Za-z])$", " \\1\\2", metadata$dataset_label)
  cleaned_label <- gsub("_", " ", cleaned_label, fixed = TRUE)
  metadata$dataset_label[default_label] <- cleaned_label[default_label]
  metadata <- metadata[order(metadata$dataset_order, metadata$dataset_label), , drop = FALSE]
  label_lookup <- setNames(metadata$dataset_label, metadata$dataset)
  order_lookup <- setNames(metadata$dataset_order, metadata$dataset)
  data$dataset_label <- unname(label_lookup[data$dataset])
  data$dataset_order <- unname(order_lookup[data$dataset])
  list(data = data, metadata = metadata)
}

build_completeness <- function(data, expected_reps, expected_taxa) {
  datasets <- unique(data$dataset)
  grid <- expand.grid(
    dataset = datasets,
    scenario = scenario_levels,
    method = method_levels,
    stringsAsFactors = FALSE
  )
  key <- interaction(
    data$dataset, data$scenario, data$method,
    drop = TRUE, lex.order = TRUE
  )
  pieces <- lapply(split(seq_len(nrow(data)), key), function(index) {
    subset <- data[index, , drop = FALSE]
    taxon_counts <- vapply(
      split(subset$taxon, subset$replicate),
      function(x) length(unique(x)),
      integer(1L)
    )
    numeric_reps <- suppressWarnings(as.integer(names(taxon_counts)))
    expected_ids <- as.character(seq_len(expected_reps))
    observed_ids <- unique(as.character(subset$replicate))
    data.frame(
      dataset = subset$dataset[[1L]],
      scenario = subset$scenario[[1L]],
      method = subset$method[[1L]],
      n_replicates = length(observed_ids),
      n_complete_taxon_sets = sum(taxon_counts == expected_taxa),
      min_taxa_per_replicate = min(taxon_counts),
      max_taxa_per_replicate = max(taxon_counts),
      replicate_ids_are_1_to_n = !anyNA(numeric_reps) &&
        setequal(as.character(numeric_reps), expected_ids),
      missing_replicate_ids = paste(setdiff(expected_ids, observed_ids), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  observed <- if (length(pieces)) do.call(rbind, pieces) else data.frame()
  completeness <- merge(
    grid, observed,
    by = c("dataset", "scenario", "method"),
    all.x = TRUE,
    sort = FALSE
  )
  for (name in c(
    "n_replicates", "n_complete_taxon_sets", "min_taxa_per_replicate", "max_taxa_per_replicate"
  )) {
    completeness[[name]][is.na(completeness[[name]])] <- 0L
  }
  completeness$replicate_ids_are_1_to_n[is.na(completeness$replicate_ids_are_1_to_n)] <- FALSE
  completeness$missing_replicate_ids[is.na(completeness$missing_replicate_ids)] <-
    paste(seq_len(expected_reps), collapse = ";")
  completeness$expected_replicates <- expected_reps
  completeness$expected_taxa_per_replicate <- expected_taxa
  completeness$complete <- completeness$n_replicates == expected_reps &
    completeness$n_complete_taxon_sets == expected_reps &
    completeness$min_taxa_per_replicate == expected_taxa &
    completeness$max_taxa_per_replicate == expected_taxa &
    completeness$replicate_ids_are_1_to_n
  completeness
}

validate_results <- function(data, metadata, completeness, options) {
  problems <- character()
  if (nrow(metadata) != options$expected_datasets) {
    problems <- c(problems, paste0(
      "found ", nrow(metadata), " datasets; expected ", options$expected_datasets
    ))
  }
  observed_scenarios <- unique(data$scenario)
  observed_methods <- unique(data$method)
  if (!setequal(observed_scenarios, scenario_levels)) {
    problems <- c(problems, "the three required depth scenarios are not all present")
  }
  if (!setequal(observed_methods, method_levels)) {
    problems <- c(problems, "the three required primary methods are not all present")
  }
  incomplete <- completeness[!completeness$complete, , drop = FALSE]
  if (nrow(incomplete)) {
    examples <- apply(
      head(incomplete[c("dataset", "scenario", "method", "n_replicates", "min_taxa_per_replicate", "max_taxa_per_replicate")], 6L),
      1L,
      paste,
      collapse = " | "
    )
    problems <- c(problems, paste0(
      nrow(incomplete), " dataset/scenario/method cells are incomplete; examples: ",
      paste(examples, collapse = "; ")
    ))
  }
  invalid_p <- is.finite(data$p_value) & (data$p_value < 0 | data$p_value > 1)
  if (any(invalid_p)) {
    problems <- c(problems, paste0(sum(invalid_p), " finite P values lie outside [0, 1]"))
  }
  if (length(problems)) {
    text <- paste(problems, collapse = "\n- ")
    if (options$allow_partial) {
      warning("Partial-results validation report:\n- ", text, call. = FALSE)
    } else {
      stop(
        "Results are not a complete ", options$expected_reps,
        " x 3 x 3 negative-control run:\n- ", text,
        "\nUse --allow-partial only for interim summaries.",
        call. = FALSE
      )
    }
  }
  invisible(problems)
}

make_replicate_summary <- function(data, expected_taxa, alpha) {
  data$valid_p <- is.finite(data$p_value) & data$p_value >= 0 & data$p_value <= 1
  data$available <- ifelse(
    is.na(data$available_explicit),
    data$valid_p,
    data$available_explicit & data$valid_p
  )
  data$analysis_p_value <- ifelse(data$available, data$p_value, 1)
  data$rejected <- data$analysis_p_value <= alpha

  key <- interaction(
    data$dataset,
    data$replicate,
    data$scenario,
    data$method,
    drop = TRUE,
    lex.order = TRUE
  )
  replicate_rows <- lapply(split(seq_len(nrow(data)), key), function(index) {
    subset <- data[index, , drop = FALSE]
    n_taxa <- length(unique(subset$taxon))
    data.frame(
      dataset = subset$dataset[[1L]],
      dataset_label = subset$dataset_label[[1L]],
      dataset_order = subset$dataset_order[[1L]],
      replicate = subset$replicate[[1L]],
      scenario = subset$scenario[[1L]],
      method = subset$method[[1L]],
      n_taxon_records = n_taxa,
      n_missing_taxon_records = max(0L, expected_taxa - n_taxa),
      n_available = sum(subset$available),
      n_unavailable = expected_taxa - sum(subset$available),
      availability_rate = sum(subset$available) / expected_taxa,
      n_rejected = sum(subset$rejected),
      type1_error = sum(subset$rejected) / expected_taxa,
      stringsAsFactors = FALSE
    )
  })
  replicate_summary <- do.call(rbind, replicate_rows)
  rownames(replicate_summary) <- NULL
  list(taxon = data, replicate = replicate_summary)
}

make_dataset_summary <- function(replicate_summary) {
  key <- interaction(
    replicate_summary$dataset,
    replicate_summary$scenario,
    replicate_summary$method,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(replicate_summary)), key), function(index) {
    subset <- replicate_summary[index, , drop = FALSE]
    n <- nrow(subset)
    mean_value <- mean(subset$type1_error)
    se_value <- if (n > 1L) stats::sd(subset$type1_error) / sqrt(n) else NA_real_
    data.frame(
      dataset = subset$dataset[[1L]],
      dataset_label = subset$dataset_label[[1L]],
      dataset_order = subset$dataset_order[[1L]],
      scenario = subset$scenario[[1L]],
      method = subset$method[[1L]],
      n_replicates = n,
      mean_type1_error = mean_value,
      monte_carlo_se = se_value,
      ci95_lower = mean_value - 1.96 * se_value,
      ci95_upper = mean_value + 1.96 * se_value,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary$scenario <- factor(summary$scenario, levels = scenario_levels)
  summary$method <- factor(summary$method, levels = method_levels)
  summary <- summary[order(
    summary$dataset_order,
    summary$method,
    summary$scenario
  ), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

make_availability_summary <- function(replicate_summary, expected_taxa) {
  key <- interaction(
    replicate_summary$dataset,
    replicate_summary$scenario,
    replicate_summary$method,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(replicate_summary)), key), function(index) {
    subset <- replicate_summary[index, , drop = FALSE]
    expected_tests <- nrow(subset) * expected_taxa
    data.frame(
      dataset = subset$dataset[[1L]],
      dataset_label = subset$dataset_label[[1L]],
      dataset_order = subset$dataset_order[[1L]],
      scenario = subset$scenario[[1L]],
      method = subset$method[[1L]],
      n_replicates_observed = nrow(subset),
      n_expected_taxon_tests_in_observed_replicates = expected_tests,
      n_taxon_records = sum(subset$n_taxon_records),
      n_missing_taxon_records = sum(subset$n_missing_taxon_records),
      n_available = sum(subset$n_available),
      n_unavailable = sum(subset$n_unavailable),
      availability_rate = sum(subset$n_available) / expected_tests,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary$scenario <- factor(summary$scenario, levels = scenario_levels)
  summary$method <- factor(summary$method, levels = method_levels)
  summary <- summary[order(
    summary$dataset_order,
    summary$method,
    summary$scenario
  ), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

write_outputs <- function(taxon, replicate, dataset, availability, completeness, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  taxon_output <- taxon[c(
    "dataset", "dataset_label", "dataset_order", "replicate", "scenario", "method",
    "taxon", "p_value", "native_p_value", "available", "analysis_p_value", "rejected",
    "status", "reason", "source", "object"
  )]
  write.csv(
    taxon_output,
    file.path(output_dir, "negative_control_taxon_pvalues.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    replicate,
    file.path(output_dir, "negative_control_replicate_type1.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    dataset,
    file.path(output_dir, "negative_control_dataset_summary.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    availability,
    file.path(output_dir, "negative_control_availability_summary.csv"),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    completeness,
    file.path(output_dir, "negative_control_completeness.csv"),
    row.names = FALSE,
    na = ""
  )
}

draw_figure <- function(dataset_summary, pdf_path, alpha) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to draw the figure.", call. = FALSE)
  }
  plot_data <- dataset_summary
  dataset_order <- unique(plot_data[order(plot_data$dataset_order), c(
    "dataset", "dataset_label", "dataset_order"
  )])
  label_levels <- rev(dataset_order$dataset_label)
  plot_data$dataset_label <- factor(plot_data$dataset_label, levels = label_levels)
  plot_data$scenario <- factor(plot_data$scenario, levels = scenario_levels)
  plot_data$method <- factor(plot_data$method, levels = method_levels)
  dodge <- ggplot2::position_dodge(width = 0.54)
  plot_upper <- max(
    alpha * 1.25,
    plot_data$mean_type1_error,
    plot_data$ci95_upper,
    na.rm = TRUE
  ) * 1.04

  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = mean_type1_error,
      y = dataset_label,
      color = scenario,
      shape = scenario
    )
  ) +
    ggplot2::geom_vline(
      xintercept = alpha,
      linetype = "dashed",
      linewidth = 0.45,
      color = "#777777"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci95_lower, xmax = ci95_upper),
      orientation = "y",
      width = 0.18,
      linewidth = 0.55,
      position = dodge
    ) +
    ggplot2::geom_point(
      size = 2.45,
      stroke = 0.9,
      position = dodge
    ) +
    ggplot2::facet_grid(. ~ method) +
    ggplot2::scale_color_manual(
      values = c(
        "Balanced depth" = "#667786",
        "Fourfold depth difference" = "#7A4E8E",
        "Sixfold depth difference" = "#D8782D"
      ),
      breaks = scenario_levels,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Balanced depth" = 1,
        "Fourfold depth difference" = 16,
        "Sixfold depth difference" = 17
      ),
      breaks = scenario_levels,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = function(x) {
        labels <- vapply(x, function(value) {
          format(round(value, 3), trim = TRUE, scientific = FALSE)
        }, character(1L))
        sub("^0\\.", ".", labels)
      },
      expand = ggplot2::expansion(mult = c(0.04, 0.08))
    ) +
    ggplot2::coord_cartesian(xlim = c(0, plot_upper), clip = "on") +
    ggplot2::labs(
      x = "Per-taxon Type I error (raw P \u2264 0.05)",
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    ggplot2::theme_classic(base_size = 10.5, base_family = "sans") +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box.spacing = grid::unit(0.15, "lines"),
      legend.spacing.x = grid::unit(0.35, "lines"),
      legend.key.width = grid::unit(1.25, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      panel.spacing.x = grid::unit(1.15, "lines"),
      axis.title.x = ggplot2::element_text(size = 11.5, margin = ggplot2::margin(t = 8)),
      axis.text.x = ggplot2::element_text(size = 9.5, color = "#222222"),
      axis.text.y = ggplot2::element_text(size = 9.7, color = "#222222"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, color = "#333333"),
      axis.line = ggplot2::element_line(linewidth = 0.5, color = "#333333"),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = pdf_path,
    plot = figure,
    device = grDevices::cairo_pdf,
    width = 10.4,
    height = 5.9,
    units = "in",
    bg = "white"
  )
  invisible(figure)
}

main <- function() {
  script_dir <- script_directory()
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE), script_dir)
  if (arguments$help) {
    usage()
    return(invisible(NULL))
  }
  arguments$output_dir <- normalizePath(arguments$output_dir, mustWork = FALSE)
  if (is.na(arguments$pdf)) {
    arguments$pdf <- file.path(arguments$output_dir, "negative_control_type1.pdf")
  }
  manifest_path <- find_manifest(arguments$manifest, script_dir)
  files <- list_input_files(arguments$input, arguments$output_dir)
  if (!length(files)) {
    stop("No CSV, TSV, or RDS result files were found.", call. = FALSE)
  }
  message("Reading ", length(files), " result file(s).")
  data <- deduplicate_results(read_all_results(files))
  metadata_result <- apply_dataset_metadata(data, manifest_path)
  data <- metadata_result$data
  metadata <- metadata_result$metadata
  completeness <- build_completeness(
    data,
    expected_reps = arguments$expected_reps,
    expected_taxa = arguments$expected_taxa
  )
  validate_results(data, metadata, completeness, arguments)
  summaries <- make_replicate_summary(
    data,
    expected_taxa = arguments$expected_taxa,
    alpha = arguments$alpha
  )
  dataset_summary <- make_dataset_summary(summaries$replicate)
  availability_summary <- make_availability_summary(
    summaries$replicate,
    expected_taxa = arguments$expected_taxa
  )
  write_outputs(
    summaries$taxon,
    summaries$replicate,
    dataset_summary,
    availability_summary,
    completeness,
    arguments$output_dir
  )
  if (!arguments$validate_only) {
    draw_figure(dataset_summary, arguments$pdf, arguments$alpha)
  }
  message(
    "Summarized ", nrow(metadata), " dataset(s), ",
    length(unique(data$replicate)), " replicate identifier(s), and ",
    format(nrow(data), big.mark = ","), " taxon-level tests."
  )
  message("Outputs: ", arguments$output_dir)
  if (!arguments$validate_only) message("Figure: ", arguments$pdf)
}

tryCatch(
  main(),
  error = function(e) {
    message("Error: ", conditionMessage(e))
    quit(status = 1L, save = "no")
  }
)
