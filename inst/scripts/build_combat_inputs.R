#!/usr/bin/env Rscript

# Reconstruct the two saved COMBAT example inputs from the public raw H5AD.
# This script builds inputs and validates identities; it does not fit an
# expression model or execute a breakpoint search.

options(stringsAsFactors = FALSE)

.usage <- function(status = 0L) {
  cat(paste(
    "Usage:",
    "  build_combat_inputs.R --h5ad PATH --output-dir PATH [--skip-source-md5]",
    "",
    "Fresh downloads are checked against the registered 6.4 GB source MD5.",
    "Use --skip-source-md5 only for a locally retained copy whose checksum was",
    "already verified; file size and all reconstructed semantic hashes remain checked.",
    sep = "\n"
  ))
  quit(save = "no", status = status)
}

.parse_options <- function(args) {
  out <- list(h5ad = NULL, output_dir = NULL, skip_source_md5 = FALSE)
  at <- 1L
  while (at <= length(args)) {
    key <- args[[at]]
    if (key == "--help") .usage(0L)
    if (key == "--skip-source-md5") {
      out$skip_source_md5 <- TRUE
      at <- at + 1L
      next
    }
    if (!key %in% c("--h5ad", "--output-dir") || at == length(args)) {
      stop("unknown or incomplete option: ", key, call. = FALSE)
    }
    value <- args[[at + 1L]]
    if (key == "--h5ad") out$h5ad <- value else out$output_dir <- value
    at <- at + 2L
  }
  if (is.null(out$h5ad) || is.null(out$output_dir)) .usage(2L)
  out
}

.script_path <- function() {
  raw <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", raw[grepl("^--file=", raw)])
  if (!length(file)) return(NA_character_)
  normalizePath(file[[1L]], mustWork = TRUE)
}

.resource_dir <- function() {
  installed <- system.file("extdata", "combat-reproduction",
                           package = "AnnotationBreakpoint")
  if (nzchar(installed)) return(installed)
  script <- .script_path()
  if (is.na(script)) stop("cannot locate packaged reconstruction resources", call. = FALSE)
  candidate <- file.path(dirname(dirname(script)), "extdata", "combat-reproduction")
  normalizePath(candidate, mustWork = TRUE)
}

.require_packages <- function() {
  required <- c("AnnotationBreakpoint", "Matrix", "digest", "jsonlite", "rhdf5", "yaml")
  absent <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(absent)) {
    stop("missing required package(s): ", paste(absent, collapse = ", "), call. = FALSE)
  }
  suppressPackageStartupMessages(library(Matrix))
}

.sha256 <- function(path) digest::digest(
  file = path, algo = "sha256", serialize = FALSE
)

.write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
    na = "NA", eol = "\n", fileEncoding = "UTF-8"
  )
}

.read_h5ad_column <- function(path, frame, column) {
  base <- paste0(frame, "/", column)
  listing <- rhdf5::h5ls(path, recursive = TRUE, datasetinfo = FALSE)
  direct <- listing$group == frame & listing$name == column
  if (!any(direct)) stop("H5AD column is absent: ", base, call. = FALSE)
  available <- paste0(listing$group, "/", listing$name)
  legacy_category_path <- paste0(frame, "/__categories/", column)
  if (legacy_category_path %in% available) {
    categories <- rhdf5::h5read(path, legacy_category_path)
    codes <- as.integer(rhdf5::h5read(path, base))
    out <- rep(NA_character_, length(codes))
    keep <- codes >= 0L
    out[keep] <- as.character(categories[codes[keep] + 1L])
    return(out)
  }
  if (listing$otype[which(direct)[[1L]]] == "H5I_DATASET") {
    return(rhdf5::h5read(path, base))
  }
  category_path <- paste0(base, "/categories")
  code_path <- paste0(base, "/codes")
  if (!all(c(category_path, code_path) %in% available)) {
    stop("unsupported H5AD column encoding: ", base, call. = FALSE)
  }
  categories <- rhdf5::h5read(path, category_path)
  codes <- as.integer(rhdf5::h5read(path, code_path))
  out <- rep(NA_character_, length(codes))
  keep <- codes >= 0L
  out[keep] <- as.character(categories[codes[keep] + 1L])
  out
}

.stream_csr_rows <- function(path, matrix_group, selected_rows, n_columns,
                             block_rows = 10000L, keep_columns = NULL) {
  indptr <- as.numeric(rhdf5::h5read(path, paste0(matrix_group, "/indptr")))
  n_rows <- length(indptr) - 1L
  selected_rows <- as.integer(selected_rows)
  if (!length(selected_rows) || any(selected_rows < 1L) ||
      any(selected_rows > n_rows) || anyDuplicated(selected_rows)) {
    stop("selected raw rows must be unique and within the H5AD matrix", call. = FALSE)
  }
  selected_lookup <- integer(n_rows)
  selected_lookup[selected_rows] <- seq_along(selected_rows)
  p <- integer(length(selected_rows) + 1L)
  j_parts <- vector("list", length(selected_rows))
  x_parts <- vector("list", length(selected_rows))
  starts <- seq.int(1L, n_rows, by = as.integer(block_rows))
  for (block_start in starts) {
    block_end <- min(n_rows, block_start + as.integer(block_rows) - 1L)
    chosen <- selected_rows[selected_rows >= block_start & selected_rows <= block_end]
    if (!length(chosen)) next
    data_start <- indptr[[block_start]] + 1
    data_end <- indptr[[block_end + 1L]]
    if (data_end < data_start) next
    index <- as.integer(rhdf5::h5read(
      path, paste0(matrix_group, "/indices"), start = data_start,
      count = data_end - data_start + 1
    ))
    value <- as.numeric(rhdf5::h5read(
      path, paste0(matrix_group, "/data"), start = data_start,
      count = data_end - data_start + 1
    ))
    offset <- indptr[[block_start]]
    for (row in chosen) {
      destination <- selected_lookup[[row]]
      lo <- indptr[[row]] - offset + 1
      hi <- indptr[[row + 1L]] - offset
      if (hi < lo) next
      row_index <- index[seq.int(lo, hi)] + 1L
      row_value <- value[seq.int(lo, hi)]
      if (!is.null(keep_columns)) {
        keep <- row_index %in% keep_columns
        row_index <- row_index[keep]
        row_value <- row_value[keep]
      }
      j_parts[[destination]] <- row_index - 1L
      x_parts[[destination]] <- row_value
    }
    rm(index, value)
  }
  p[-1L] <- cumsum(lengths(j_parts))
  j <- as.integer(unlist(j_parts, use.names = FALSE))
  x <- as.numeric(unlist(x_parts, use.names = FALSE))
  if (!is.null(keep_columns)) {
    mapped <- match(j + 1L, keep_columns)
    if (anyNA(mapped)) stop("selected-column mapping failed", call. = FALSE)
    j <- as.integer(mapped - 1L)
    n_columns <- length(keep_columns)
  }
  methods::new(
    "dgRMatrix", p = as.integer(p), j = j, x = x,
    Dim = as.integer(c(length(selected_rows), n_columns))
  )
}

.extract_counts <- function(path, selected_rows, feature_ids, rna_columns) {
  csr <- .stream_csr_rows(
    path, "/layers/raw", selected_rows, length(feature_ids),
    keep_columns = rna_columns
  )
  counts <- Matrix::t(methods::as(csr, "CsparseMatrix"))
  counts <- methods::as(counts, "dgCMatrix")
  rownames(counts) <- feature_ids[rna_columns]
  counts
}

.extract_columns <- function(path, selected_rows, n_features, columns) {
  as.matrix(.stream_csr_rows(
    path, "/X", selected_rows, n_features, keep_columns = columns
  ))
}

.map_json <- function(x) jsonlite::toJSON(
  as.list(x), auto_unbox = TRUE, digits = NA
)

.write_case <- function(case_id, gene, counts, cells, specimens, allowed_labels,
                        config, analysis_spec, output_root, expected) {
  input_dir <- file.path(output_root, "inputs", gene)
  spec_dir <- file.path(output_root, "analysis_specs")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(spec_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(counts, file.path(input_dir, "counts.rds"), version = 3L, compress = TRUE)
  .write_tsv(cells, file.path(input_dir, "cells.tsv"))
  .write_tsv(specimens, file.path(input_dir, "specimens.tsv"))
  .write_tsv(allowed_labels, file.path(input_dir, "allowed_labels.tsv"))
  claim <- data.frame(
    focal_gene = gene,
    focal_label = config$claim$focal_label,
    alternative_label = config$claim$alternative_label,
    case_condition = config$claim$case_condition,
    control_condition = config$claim$control_condition,
    q = as.integer(config$claim$q),
    delta = as.numeric(config$claim$delta),
    tolerance = as.numeric(config$claim$tolerance),
    feature_universe = paste(rownames(counts), collapse = "|"),
    estimator_id = config$claim$estimator_id,
    stringsAsFactors = FALSE
  )
  .write_tsv(claim, file.path(input_dir, "claim.tsv"))
  constraints <- config$constraints
  constraint_table <- data.frame(
    min_cells_CD4.TCM = constraints$min_cells_CD4.TCM,
    min_cells_CD4.NAIVE = constraints$min_cells_CD4.NAIVE,
    min_library_CD4.TCM = constraints$min_library_CD4.TCM,
    min_library_CD4.NAIVE = constraints$min_library_CD4.NAIVE,
    max_changes_global = constraints$max_changes_global,
    max_changes_per_sample = .map_json(unlist(constraints$max_changes_per_sample)),
    max_changes_per_donor = .map_json(unlist(constraints$max_changes_per_donor)),
    min_donors_per_condition = constraints$min_donors_per_condition,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  .write_tsv(constraint_table, file.path(input_dir, "constraints.tsv"))
  spec_path <- file.path(spec_dir, paste0(gene, ".yaml"))
  yaml::write_yaml(analysis_spec, spec_path)

  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  loaded <- reader(input_dir, analysis_spec = yaml::read_yaml(spec_path))
  validation <- do.call(AnnotationBreakpoint::validate_inputs, c(
    loaded[c("counts", "cells", "specimens", "allowed_labels", "claim",
             "constraints", "analysis_spec")],
    list(runtime_policy = "strict")
  ))
  if (!isTRUE(validation$ok)) {
    stop(case_id, " reconstructed input is invalid: ",
         paste(validation$errors, collapse = "; "), call. = FALSE)
  }
  hashes <- validation$input_hashes
  expected_case <- expected[expected$case_id == case_id, , drop = FALSE]
  old_components <- setdiff(expected_case$component, "old_specification_hash")
  matches <- vapply(old_components, function(component) {
    identical(
      unname(hashes[[component]]),
      expected_case$sha256[expected_case$component == component][[1L]]
    )
  }, logical(1L))
  if (!all(matches)) {
    stop(
      case_id, " reconstructed semantic hashes differ: ",
      paste(old_components[!matches], collapse = ", "), call. = FALSE
    )
  }
  hash_value <- getFromNamespace(".ab_hash", "AnnotationBreakpoint")
  historical_specification_hash <- hash_value(list(
    claim = validation$data$claim,
    analysis_spec = validation$data$analysis_spec
  ))
  expected_historical <- expected_case$sha256[
    expected_case$component == "old_specification_hash"
  ][[1L]]
  if (!identical(historical_specification_hash, expected_historical)) {
    stop(case_id, " historical specification hash differs", call. = FALSE)
  }
  data.frame(
    case_id = case_id,
    feature_id = gene,
    component = c(old_components, "old_specification_hash", "strict_specification_hash"),
    expected_sha256 = c(
      expected_case$sha256[match(old_components, expected_case$component)],
      expected_historical,
      NA_character_
    ),
    observed_sha256 = c(
      unname(hashes[old_components]),
      historical_specification_hash,
      unname(hashes[["specification_hash"]])
    ),
    matches = c(rep(TRUE, length(old_components) + 1L), NA),
    translation = c(
      rep("unchanged_component", length(old_components)),
      "historical_without_runtime_policy",
      "0.3.0_strict_runtime_policy_included"
    ),
    stringsAsFactors = FALSE
  )
}

opts <- .parse_options(commandArgs(trailingOnly = TRUE))
.require_packages()
h5ad <- normalizePath(opts$h5ad, mustWork = TRUE)
output_root <- normalizePath(opts$output_dir, mustWork = FALSE)
if (dir.exists(output_root) && length(list.files(output_root, all.files = TRUE,
                                                 no.. = TRUE))) {
  stop("output directory must be absent or empty", call. = FALSE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
resources <- .resource_dir()
config <- jsonlite::read_json(file.path(resources, "frozen_config.json"),
                              simplifyVector = FALSE)
expected <- utils::read.delim(
  file.path(resources, "expected_semantic_hashes.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
selection_connection <- gzfile(
  file.path(resources, "evaluation_cells.tsv.gz"), open = "rt"
)
selection <- tryCatch(
  utils::read.delim(
    selection_connection, stringsAsFactors = FALSE, check.names = FALSE
  ),
  finally = close(selection_connection)
)
specimens <- utils::read.delim(
  file.path(resources, "specimens.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character"
)
analysis_spec <- yaml::read_yaml(file.path(resources, "analysis_spec.yaml"))

source_bytes <- as.numeric(file.info(h5ad)$size)
if (!identical(source_bytes, as.numeric(config$source$bytes))) {
  stop("raw H5AD byte count differs from the registered source", call. = FALSE)
}
source_md5 <- if (isTRUE(opts$skip_source_md5)) {
  NA_character_
} else {
  unname(tools::md5sum(h5ad))
}
if (!isTRUE(opts$skip_source_md5) &&
    !identical(source_md5, as.character(config$source$checksum))) {
  stop("raw H5AD MD5 differs from the registered source", call. = FALSE)
}

source_cell_id <- as.character(.read_h5ad_column(h5ad, "/obs", "_index"))
selected_rows <- as.integer(selection$source_row)
if (anyNA(selected_rows) || anyDuplicated(selected_rows) ||
    !identical(source_cell_id[selected_rows], selection$cell_id)) {
  stop("frozen cell-to-source-row mapping does not match the H5AD", call. = FALSE)
}
source_fields <- list(
  sample_id = "scRNASeq_sample_ID", donor_id = "COMBAT_ID",
  condition = "Source", original_label = "Annotation_minor_subset",
  processing_pool = "Pool_ID"
)
for (field in names(source_fields)) {
  value <- as.character(.read_h5ad_column(h5ad, "/obs", source_fields[[field]]))
  if (!identical(value[selected_rows], as.character(selection[[field]]))) {
    stop("frozen selected metadata differs from raw H5AD field ", field, call. = FALSE)
  }
}

feature_ids <- as.character(.read_h5ad_column(h5ad, "/var", "_index"))
feature_types <- as.character(.read_h5ad_column(h5ad, "/var", "feature_types"))
rna_columns <- which(feature_types == config$feature_axis$rna_type)
if (length(rna_columns) != as.integer(config$feature_axis$expected_rna_features) ||
    anyDuplicated(feature_ids)) {
  stop("raw feature axis differs from the registered RNA universe", call. = FALSE)
}
marker_names <- c(config$ambiguity$marker_naive, config$ambiguity$marker_memory)
marker_columns <- match(marker_names, feature_ids)
if (anyNA(marker_columns) || any(feature_types[marker_columns] != "Antibody Capture")) {
  stop("registered ADT marker columns are absent", call. = FALSE)
}
marker_matrix <- .extract_columns(
  h5ad, selected_rows, length(feature_ids), marker_columns
)
adt_score <- asinh(marker_matrix[, 1L]) - asinh(marker_matrix[, 2L])
boundary <- as.numeric(config$ambiguity$boundary)
near <- abs(adt_score - boundary) <=
  as.numeric(config$ambiguity$near_boundary_multiplier) *
  as.numeric(config$ambiguity$scale)
label <- as.character(selection$original_label)
discordant <-
  (label == config$ambiguity$alternative_label & adt_score < boundary) |
  (label == config$ambiguity$focal_label & adt_score > boundary)
ambiguous <- near | discordant
expected_ambiguity <- as.logical(selection$adt_ambiguity)
if (!identical(ambiguous, expected_ambiguity) ||
    sum(ambiguous) != as.integer(config$dimensions$ambiguities)) {
  stop("raw ADT ambiguity reconstruction differs from the frozen panel", call. = FALSE)
}

counts <- .extract_counts(h5ad, selected_rows, feature_ids, rna_columns)
colnames(counts) <- selection$cell_id
if (!identical(dim(counts), c(
      as.integer(config$dimensions$features), as.integer(config$dimensions$cells)
    )) || length(counts@x) != as.integer(config$dimensions$nonzero_counts) ||
    any(!is.finite(counts@x)) || any(counts@x < 0) ||
    any(counts@x != floor(counts@x))) {
  stop("reconstructed sparse RNA matrix differs structurally", call. = FALSE)
}

cells <- selection[c(
  "cell_id", "sample_id", "donor_id", "condition", "original_label",
  "processing_pool", "adt_ambiguity"
)]
cells$adt_ambiguity <- as.logical(cells$adt_ambiguity)
allowed_labels <- rbind(
  data.frame(
    cell_id = cells$cell_id,
    allowed_label = cells$original_label,
    evidence_source = config$ambiguity$evidence_source_original,
    evidence_version = config$ambiguity$evidence_version_original,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cell_id = cells$cell_id[ambiguous],
    allowed_label = ifelse(
      cells$original_label[ambiguous] == config$ambiguity$alternative_label,
      config$ambiguity$focal_label, config$ambiguity$alternative_label
    ),
    evidence_source = config$ambiguity$evidence_source_alternative,
    evidence_version = config$ambiguity$evidence_version_alternative,
    stringsAsFactors = FALSE
  )
)
allowed_labels <- allowed_labels[order(
  allowed_labels$cell_id, allowed_labels$allowed_label, method = "radix"
), , drop = FALSE]

semantic <- rbind(
  .write_case(
    "REAL_FAM118A", "FAM118A", counts, cells, specimens, allowed_labels,
    config, analysis_spec, output_root, expected
  ),
  .write_case(
    "REAL_PDXK", "PDXK", counts, cells, specimens, allowed_labels,
    config, analysis_spec, output_root, expected
  )
)
.write_tsv(semantic, file.path(output_root, "semantic_hashes.tsv"))
receipt <- list(
  schema_version = "1.0.0",
  registration_id = config$registration_id,
  source = list(
    path = h5ad,
    bytes = source_bytes,
    expected_md5 = config$source$checksum,
    observed_md5 = source_md5,
    md5_check = if (isTRUE(opts$skip_source_md5)) {
      "skipped_for_previously_verified_local_source"
    } else "executed_and_matched"
  ),
  dimensions = list(
    cells = nrow(cells), features = nrow(counts), donors = nrow(specimens),
    ambiguities = sum(ambiguous), nonzero_counts = length(counts@x)
  ),
  cases = c("REAL_FAM118A", "REAL_PDXK"),
  endpoint_fits_or_searches = 0L,
  semantic_hashes = "semantic_hashes.tsv",
  file_manifest = "file_manifest.tsv"
)
jsonlite::write_json(
  receipt, file.path(output_root, "reconstruction_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA
)
tracked <- list.files(output_root, recursive = TRUE, full.names = TRUE)
tracked <- tracked[file.info(tracked)$isdir %in% FALSE]
tracked <- tracked[basename(tracked) != "file_manifest.tsv"]
file_manifest <- data.frame(
  path = substring(tracked, nchar(output_root) + 2L),
  bytes = as.numeric(file.info(tracked)$size),
  sha256 = vapply(tracked, .sha256, character(1L)),
  stringsAsFactors = FALSE
)
.write_tsv(file_manifest, file.path(output_root, "file_manifest.tsv"))
cat(
  "reconstruction_complete cells=", nrow(cells),
  " features=", nrow(counts), " donors=", nrow(specimens),
  " ambiguities=", sum(ambiguous), " cases=2\n", sep = ""
)
