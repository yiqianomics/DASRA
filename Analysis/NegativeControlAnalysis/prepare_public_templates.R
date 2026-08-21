#!/usr/bin/env Rscript

study_names <- c(
  "MetaCardis_2020_a",
  "LifeLinesDeep_2016",
  "AsnicarF_2021",
  "ZeeviD_2015",
  "ShaoY_2019",
  "YachidaS_2019",
  "SchirmerM_2016",
  "JieZ_2017",
  "QinJ_2012",
  "VilaAV_2018"
)

study_seeds <- stats::setNames(202608210L + seq_along(study_names), study_names)
target_sample_count <- 200L
target_taxon_count <- 30L
minimum_prevalence <- 0.10
minimum_mean_abundance <- 1e-5

usage <- function() {
  cat(
    paste0(
      "Prepare genus-level public-data templates for the negative-control analysis.\n\n",
      "Usage:\n",
      "  Rscript prepare_public_templates.R [--study all|STUDY_NAME] [--force]\n\n",
      "Options:\n",
      "  --study   Prepare all ten studies or one named study. Default: all.\n",
      "  --force   Rebuild an existing template after validating the replacement.\n",
      "  --help    Show this message.\n"
    )
  )
}

parse_cli <- function(args) {
  result <- list(study = "all", force = FALSE, help = FALSE)
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("--help", "-h")) {
      result$help <- TRUE
    } else if (argument == "--force") {
      result$force <- TRUE
    } else if (startsWith(argument, "--study=")) {
      result$study <- sub("^--study=", "", argument)
    } else if (argument == "--study") {
      if (index == length(args)) {
        stop("--study requires a value.", call. = FALSE)
      }
      index <- index + 1L
      result$study <- args[[index]]
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
    index <- index + 1L
  }
  result
}

script_directory <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) != 1L) {
    stop("Run this file with Rscript.", call. = FALSE)
  }
  script_path <- sub("^--file=", "", file_argument)
  dirname(normalizePath(script_path, mustWork = TRUE))
}

is_present_identifier <- function(x) {
  normalized <- tolower(trimws(as.character(x)))
  !is.na(x) & nzchar(normalized) &
    !normalized %in% c("na", "n/a", "nan", "none", "null", "unknown",
                       "not available", "not provided")
}

extract_long_rank <- function(lineage, prefix) {
  vapply(
    strsplit(as.character(lineage), "|", fixed = TRUE),
    function(parts) {
      match_index <- which(startsWith(parts, prefix))
      if (length(match_index) == 0L) {
        return(NA_character_)
      }
      value <- substring(parts[[match_index[[1L]]]], nchar(prefix) + 1L)
      if (is_present_identifier(value)) value else NA_character_
    },
    character(1L)
  )
}

taxonomy_column <- function(row_data, lineage, rank, prefix) {
  value <- rep(NA_character_, length(lineage))
  if (rank %in% colnames(row_data)) {
    value <- as.character(row_data[[rank]])
  }
  parsed <- extract_long_rank(lineage, prefix)
  replace_index <- !is_present_identifier(value) & is_present_identifier(parsed)
  value[replace_index] <- parsed[replace_index]
  value
}

genus_labels <- function(row_data, lineage) {
  genus <- taxonomy_column(row_data, lineage, "genus", "g__")
  unresolved <- !is_present_identifier(genus)

  for (rank_specification in list(
    c("family", "f__"),
    c("order", "o__"),
    c("class", "c__"),
    c("phylum", "p__"),
    c("superkingdom", "k__")
  )) {
    if (!any(unresolved)) {
      break
    }
    rank <- rank_specification[[1L]]
    higher_label <- taxonomy_column(
      row_data,
      lineage,
      rank,
      rank_specification[[2L]]
    )
    fill <- unresolved & is_present_identifier(higher_label)
    genus[fill] <- paste0("Unclassified_genus__within_", rank, "__", higher_label[fill])
    unresolved <- !is_present_identifier(genus)
  }

  genus[unresolved] <- "Unclassified_genus__unresolved_lineage"
  trimws(genus)
}

normalize_rows <- function(matrix) {
  totals <- rowSums(matrix)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Every selected sample must have positive finite community mass.", call. = FALSE)
  }
  matrix / totals
}

aggregate_species_to_genus <- function(experiment) {
  abundance <- as.matrix(SummarizedExperiment::assay(experiment, "relative_abundance"))
  storage.mode(abundance) <- "double"
  if (length(abundance) == 0L || any(!is.finite(abundance))) {
    stop("The relative-abundance assay is empty or non-finite.", call. = FALSE)
  }
  if (any(abundance < -1e-12)) {
    stop("The relative-abundance assay contains negative values.", call. = FALSE)
  }
  abundance[abundance < 0] <- 0

  row_data <- as.data.frame(SummarizedExperiment::rowData(experiment))
  labels <- genus_labels(row_data, rownames(abundance))
  genus_by_sample <- rowsum(abundance, group = labels, reorder = FALSE)
  sample_by_genus <- t(genus_by_sample)
  sample_by_genus <- sample_by_genus[, colSums(sample_by_genus) > 0, drop = FALSE]

  original_mass <- colSums(abundance)
  aggregated_mass <- rowSums(sample_by_genus)
  mass_tolerance <- 1e-10 * max(1, max(abs(original_mass)))
  if (max(abs(original_mass - aggregated_mass)) > mass_tolerance) {
    stop("Genus aggregation did not preserve the full community mass.", call. = FALSE)
  }

  probability <- normalize_rows(sample_by_genus)
  list(
    probability = probability,
    source_feature_count = nrow(abundance),
    genus_count = ncol(probability),
    unresolved_source_feature_count = sum(
      labels == "Unclassified_genus__unresolved_lineage"
    ),
    pre_normalization_mass = original_mass
  )
}

select_test_taxa <- function(probability) {
  prevalence <- colMeans(probability > 0)
  mean_abundance <- colMeans(probability)
  eligible <- which(
    prevalence >= minimum_prevalence &
      mean_abundance >= minimum_mean_abundance &
      is.finite(mean_abundance)
  )
  if (length(eligible) < target_taxon_count) {
    stop(
      "Only ", length(eligible), " genera passed the template filter; ",
      target_taxon_count, " are required.",
      call. = FALSE
    )
  }

  eligible <- eligible[order(mean_abundance[eligible])]
  if (length(eligible) > target_taxon_count) {
    take <- unique(round(seq(1, length(eligible), length.out = target_taxon_count)))
    eligible <- eligible[take]
  }

  summary <- data.frame(
    genus = colnames(probability),
    prevalence = unname(prevalence),
    mean_abundance = unname(mean_abundance),
    eligible = colnames(probability) %in% colnames(probability)[which(
      prevalence >= minimum_prevalence &
        mean_abundance >= minimum_mean_abundance &
        is.finite(mean_abundance)
    )],
    selected = seq_len(ncol(probability)) %in% eligible,
    stringsAsFactors = FALSE
  )

  list(names = colnames(probability)[eligible], summary = summary)
}

select_study_samples <- function(metadata, available_sample_ids, study, seed) {
  metadata <- metadata[
    metadata$study_name == study &
      tolower(trimws(as.character(metadata$body_site))) == "stool" &
      metadata$sample_id %in% available_sample_ids,
    ,
    drop = FALSE
  ]
  metadata <- metadata[order(as.character(metadata$sample_id)), , drop = FALSE]

  if (anyDuplicated(metadata$sample_id)) {
    stop("Sample identifiers are duplicated within ", study, ".", call. = FALSE)
  }

  subject_is_present <- is_present_identifier(metadata$subject_id)
  subject_key <- ifelse(
    subject_is_present,
    paste0("subject:", as.character(metadata$subject_id)),
    paste0("sample:", as.character(metadata$sample_id))
  )

  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(seed)
  randomized_rows <- sample.int(nrow(metadata), nrow(metadata), replace = FALSE)
  representative_rows <- randomized_rows[!duplicated(subject_key[randomized_rows])]
  unique_subject_count <- length(representative_rows)
  if (unique_subject_count < target_sample_count) {
    stop(
      study, " has only ", unique_subject_count,
      " eligible independent stool samples; ", target_sample_count, " are required.",
      call. = FALSE
    )
  }

  selected_rows <- sample(representative_rows, target_sample_count, replace = FALSE)
  selected <- metadata[selected_rows, , drop = FALSE]
  selected$selection_subject_key <- subject_key[selected_rows]
  selected <- selected[order(as.character(selected$sample_id)), , drop = FALSE]

  list(
    metadata = selected,
    stool_assay_sample_count = nrow(metadata),
    independent_stool_sample_count = unique_subject_count,
    subject_id_sample_count = sum(startsWith(selected$selection_subject_key, "subject:")),
    sample_id_fallback_count = sum(startsWith(selected$selection_subject_key, "sample:"))
  )
}

compact_sample_metadata <- function(metadata) {
  columns <- intersect(
    c(
      "study_name", "sample_id", "subject_id", "selection_subject_key",
      "body_site", "study_condition", "disease", "age", "age_category",
      "gender", "country", "sequencing_platform", "number_reads", "PMID",
      "NCBI_accession"
    ),
    colnames(metadata)
  )
  metadata[, columns, drop = FALSE]
}

validate_template <- function(template, expected_study = NULL) {
  required <- c(
    "dataset_id", "label", "dataset_index", "source", "sample_selection",
    "sample_metadata", "probability", "evaluation_taxa", "taxon_summary",
    "aggregation"
  )
  if (!is.list(template) || !all(required %in% names(template))) {
    return("required components are missing")
  }
  if (!is.null(expected_study) && !identical(template$dataset_id, expected_study)) {
    return("study name does not match the requested study")
  }

  probability <- template$probability
  if (!is.matrix(probability) || nrow(probability) != target_sample_count ||
      ncol(probability) < target_taxon_count) {
    return("probability matrix has the wrong dimensions")
  }
  if (any(!is.finite(probability)) || any(probability < 0)) {
    return("probability matrix is non-finite or negative")
  }
  if (max(abs(rowSums(probability) - 1)) > 1e-10) {
    return("probability rows do not sum to one")
  }
  if (anyDuplicated(rownames(probability)) || anyDuplicated(colnames(probability))) {
    return("sample or genus names are duplicated")
  }
  if (length(template$evaluation_taxa) != target_taxon_count ||
      anyDuplicated(template$evaluation_taxa) ||
      !all(template$evaluation_taxa %in% colnames(probability))) {
    return("selected genera are incomplete or invalid")
  }
  if (nrow(template$sample_metadata) != target_sample_count ||
      !identical(as.character(template$sample_metadata$sample_id), rownames(probability))) {
    return("sample metadata is not aligned to the probability matrix")
  }
  if (anyDuplicated(template$sample_metadata$selection_subject_key)) {
    return("more than one sample was retained for a subject")
  }
  NULL
}

save_template <- function(template, path) {
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  saveRDS(template, temporary_path, compress = "xz")
  temporary_template <- readRDS(temporary_path)
  problem <- validate_template(temporary_template, template$dataset_id)
  if (!is.null(problem)) {
    stop("Refusing to save an invalid template: ", problem, ".", call. = FALSE)
  }

  if (file.exists(path)) {
    saved <- file.copy(temporary_path, path, overwrite = TRUE)
    unlink(temporary_path)
  } else {
    saved <- file.rename(temporary_path, path)
  }
  if (!saved) {
    stop("Could not save template: ", basename(path), call. = FALSE)
  }

  saved_template <- readRDS(path)
  problem <- validate_template(saved_template, template$dataset_id)
  if (!is.null(problem)) {
    stop("Saved template failed validation: ", problem, ".", call. = FALSE)
  }
  invisible(path)
}

manifest_row <- function(study, template_directory) {
  filename <- paste0(study, "_genus_template.rds")
  path <- file.path(template_directory, filename)
  pending <- data.frame(
    study_order = match(study, study_names),
    study_name = study,
    status = "pending",
    template_file = file.path("data", "public_templates", filename),
    source_resource = NA_character_,
    sample_seed = unname(study_seeds[[study]]),
    n_stool_assay_samples = NA_integer_,
    n_independent_stool_samples = NA_integer_,
    n_samples = NA_integer_,
    n_genera = NA_integer_,
    n_eligible_genera = NA_integer_,
    n_test_genera = NA_integer_,
    package_version = NA_character_,
    PMID = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    return(pending)
  }

  template <- tryCatch(readRDS(path), error = identity)
  if (inherits(template, "error")) {
    pending$status <- "invalid"
    return(pending)
  }
  problem <- validate_template(template, study)
  if (!is.null(problem)) {
    pending$status <- "invalid"
    return(pending)
  }

  pending$status <- "ready"
  pending$source_resource <- template$source$resource_title
  pending$n_stool_assay_samples <- template$sample_selection$stool_assay_sample_count
  pending$n_independent_stool_samples <-
    template$sample_selection$independent_stool_sample_count
  pending$n_samples <- nrow(template$probability)
  pending$n_genera <- ncol(template$probability)
  pending$n_eligible_genera <- sum(template$taxon_summary$eligible)
  pending$n_test_genera <- length(template$evaluation_taxa)
  pending$package_version <- template$source$package_version
  pending$PMID <- paste(template$source$PMID, collapse = ";")
  pending
}

write_manifest <- function(template_directory, manifest_path) {
  manifest <- do.call(
    rbind,
    lapply(study_names, manifest_row, template_directory = template_directory)
  )
  temporary_path <- tempfile(
    pattern = ".public_template_manifest.",
    tmpdir = dirname(manifest_path)
  )
  on.exit(if (file.exists(temporary_path)) unlink(temporary_path), add = TRUE)
  utils::write.csv(manifest, temporary_path, row.names = FALSE, na = "")
  if (!file.copy(temporary_path, manifest_path, overwrite = TRUE)) {
    stop("Could not update the template manifest.", call. = FALSE)
  }
  invisible(manifest)
}

prepare_study <- function(study, metadata, template_directory, force = FALSE) {
  output_path <- file.path(
    template_directory,
    paste0(study, "_genus_template.rds")
  )
  if (file.exists(output_path) && !force) {
    existing <- tryCatch(readRDS(output_path), error = identity)
    problem <- if (inherits(existing, "error")) {
      "file cannot be read"
    } else {
      validate_template(existing, study)
    }
    if (is.null(problem)) {
      message("[", study, "] Existing validated template retained.")
      return(invisible(output_path))
    }
    stop(
      "Existing template for ", study, " is invalid (", problem,
      "). Use --force to rebuild it.",
      call. = FALSE
    )
  }

  message("[", study, "] Fetching relative-abundance data.")
  resource_pattern <- paste0("\\.", study, "\\.relative_abundance$")
  resources <- suppressMessages(
    curatedMetagenomicData::curatedMetagenomicData(
      resource_pattern,
      dryrun = FALSE,
      counts = FALSE,
      rownames = "long"
    )
  )
  if (length(resources) != 1L) {
    stop("Expected exactly one relative-abundance resource for ", study, ".",
         call. = FALSE)
  }
  resource_title <- names(resources)[[1L]]
  experiment <- resources[[1L]]

  selected <- select_study_samples(
    metadata = metadata,
    available_sample_ids = colnames(experiment),
    study = study,
    seed = unname(study_seeds[[study]])
  )
  selected_metadata <- selected$metadata
  selected_ids <- as.character(selected_metadata$sample_id)
  experiment <- experiment[, selected_ids]

  aggregated <- aggregate_species_to_genus(experiment)
  probability <- aggregated$probability
  if (!identical(rownames(probability), selected_ids)) {
    stop("Assay samples are not aligned to the deterministic selection.", call. = FALSE)
  }
  taxa <- select_test_taxa(probability)
  selected_metadata <- compact_sample_metadata(selected_metadata)

  pmid <- unique(as.character(selected_metadata$PMID))
  pmid <- pmid[is_present_identifier(pmid)]
  template <- list(
    dataset_id = study,
    label = study,
    dataset_index = match(study, study_names),
    source = list(
      package = "curatedMetagenomicData",
      package_version = as.character(
        utils::packageVersion("curatedMetagenomicData")
      ),
      resource_title = resource_title,
      data_type = "relative_abundance",
      body_site = "stool",
      PMID = sort(pmid)
    ),
    sample_selection = list(
      seed = unname(study_seeds[[study]]),
      target_sample_count = target_sample_count,
      stool_assay_sample_count = selected$stool_assay_sample_count,
      independent_stool_sample_count = selected$independent_stool_sample_count,
      subject_id_sample_count = selected$subject_id_sample_count,
      sample_id_fallback_count = selected$sample_id_fallback_count,
      rule = paste0(
        "stool samples only; one randomized sample per nonmissing subject_id, ",
        "otherwise sample_id; then a deterministic random sample of 200"
      )
    ),
    sample_metadata = selected_metadata,
    probability = probability,
    evaluation_taxa = taxa$names,
    taxon_summary = taxa$summary,
    aggregation = list(
      source_rank = "species",
      target_rank = "genus",
      source_feature_count = aggregated$source_feature_count,
      genus_count = aggregated$genus_count,
      unresolved_source_feature_count =
        aggregated$unresolved_source_feature_count,
      full_community_retained = TRUE,
      normalized_to_unit_sample_mass = TRUE,
      test_taxon_filter = list(
        minimum_prevalence = minimum_prevalence,
        minimum_mean_abundance = minimum_mean_abundance,
        selection_rule = paste0(
          "sort eligible genera by mean abundance and take round(seq(1, n, ",
          "length.out = 30))"
        )
      )
    )
  )

  save_template(template, output_path)
  message(
    "[", study, "] Ready: ", nrow(probability), " samples, ",
    ncol(probability), " genera, ", sum(taxa$summary$eligible),
    " eligible, 30 selected."
  )
  invisible(output_path)
}

arguments <- parse_cli(commandArgs(trailingOnly = TRUE))
if (arguments$help) {
  usage()
  quit(save = "no", status = 0L)
}
if (!identical(arguments$study, "all") &&
    !arguments$study %in% study_names) {
  stop(
    "Unknown study: ", arguments$study, ". Choose one of: ",
    paste(study_names, collapse = ", "), ".",
    call. = FALSE
  )
}

analysis_directory <- script_directory()
local_library <- file.path(analysis_directory, "R_lib")
if (!dir.exists(local_library)) {
  stop("The analysis-local R library is missing: R_lib", call. = FALSE)
}
.libPaths(c(local_library, .libPaths()))
if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
  stop("curatedMetagenomicData is not installed in the analysis R library.",
       call. = FALSE)
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("SummarizedExperiment is required to read the public data.",
       call. = FALSE)
}

template_directory <- file.path(analysis_directory, "data", "public_templates")
dir.create(template_directory, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(
  analysis_directory,
  "data",
  "public_template_manifest.csv"
)
metadata_environment <- new.env(parent = emptyenv())
utils::data(
  "sampleMetadata",
  package = "curatedMetagenomicData",
  envir = metadata_environment
)
if (!exists("sampleMetadata", envir = metadata_environment, inherits = FALSE)) {
  stop("curatedMetagenomicData sample metadata could not be loaded.",
       call. = FALSE)
}
metadata <- get(
  "sampleMetadata",
  envir = metadata_environment,
  inherits = FALSE
)

requested_studies <- if (identical(arguments$study, "all")) {
  study_names
} else {
  arguments$study
}

for (study in requested_studies) {
  prepare_study(
    study = study,
    metadata = metadata,
    template_directory = template_directory,
    force = arguments$force
  )
  write_manifest(template_directory, manifest_path)
}

manifest <- write_manifest(template_directory, manifest_path)
message(
  "Manifest updated: ", file.path("data", "public_template_manifest.csv"),
  " (", sum(manifest$status == "ready"), "/", nrow(manifest), " ready)."
)
