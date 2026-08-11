doram_toy_data <- function(n = 20L, taxa = c("taxon_a", "taxon_b"),
                           counts = NULL) {
  stopifnot(n %% 2L == 0L, length(taxa) >= 1L)
  ids <- sprintf("sample_%03d", seq_len(n))
  if (is.null(counts)) {
    counts <- matrix(
      0, nrow = n, ncol = length(taxa),
      dimnames = list(ids, taxa)
    )
  } else {
    counts <- as.matrix(counts)
    rownames(counts) <- ids
    if (is.null(colnames(counts))) colnames(counts) <- taxa
  }
  metadata <- data.frame(
    group = rep(0:1, each = n / 2L),
    original_library_size = 2000L + 11L * seq_len(n),
    adjustment = (seq_len(n) - mean(seq_len(n))) / n,
    row.names = ids,
    check.names = FALSE
  )
  list(counts = counts, metadata = metadata, ids = ids)
}

doram_fit_toy <- function(x, ...) {
  DORAM(
    counts = x$counts,
    metadata = x$metadata,
    group = "group",
    library_size = "original_library_size",
    ...
  )
}

doram_fixture_metadata <- function(fixture) {
  covariates <- as.data.frame(fixture$covariates, check.names = FALSE)
  metadata <- data.frame(
    group = unname(fixture$group),
    original_library_size = unname(fixture$library_size),
    covariates,
    row.names = fixture$subject_id,
    check.names = FALSE
  )
  metadata
}

doram_public_result_columns <- c(
  "taxon",
  "p_occupancy", "q_occupancy",
  "p_abundance", "q_abundance",
  "p_joint", "q_joint"
)

doram_detail_columns <- c(
  "taxon", "endpoint", "available", "estimate", "statistic", "df",
  "p_value", "log_p_value", "q_value", "status"
)

doram_boundary_columns <- c(
  "taxon", "endpoint", "stage", "parameter", "parameter_block",
  "value", "side", "bound", "distance", "active_limit",
  "transformed_value"
)

doram_top_level_components <- c(
  "results", "details", "descriptives", "diagnostics", "posterior",
  "contrast", "sample_id", "taxa", "settings", "call"
)
