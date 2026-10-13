# Internal M1 fixture helper. Normal analyses do not require expected truth tables.
read_fixture_directory <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one nonempty directory path", call. = FALSE)
  }
  fixture_dir <- normalizePath(path, mustWork = TRUE)
  if (!dir.exists(fixture_dir)) {
    stop("fixture path is not a directory", call. = FALSE)
  }
  table_names <- c(
    "counts", "cells", "specimens", "allowed_labels",
    "claim", "constraints", "expected_results"
  )
  paths <- setNames(file.path(fixture_dir, paste0(table_names, ".tsv")), table_names)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      paste0("fixture directory is missing: ", paste(basename(missing), collapse = ", ")),
      call. = FALSE
    )
  }
  read_table <- function(file) {
    utils::read.delim(
      file,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      colClasses = "character"
    )
  }
  raw_counts <- read_table(paths[["counts"]])
  if (!"gene_id" %in% names(raw_counts)) {
    stop("counts.tsv must begin with a gene_id column", call. = FALSE)
  }
  gene_ids <- raw_counts$gene_id
  count_tokens <- raw_counts[, setdiff(names(raw_counts), "gene_id"), drop = FALSE]
  numeric_counts <- suppressWarnings(matrix(
    as.numeric(as.matrix(count_tokens)),
    nrow = nrow(count_tokens),
    ncol = ncol(count_tokens),
    dimnames = list(gene_ids, names(count_tokens))
  ))
  counts <- Matrix::Matrix(numeric_counts, sparse = TRUE)
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")

  result <- list(
    counts = counts,
    cells = read_table(paths[["cells"]]),
    specimens = read_table(paths[["specimens"]]),
    allowed_labels = read_table(paths[["allowed_labels"]]),
    claim = read_table(paths[["claim"]]),
    constraints = read_table(paths[["constraints"]]),
    expected_results = read_table(paths[["expected_results"]])
  )
  attr(result, "fixture_dir") <- fixture_dir
  attr(result, "source_paths") <- paths
  attr(result, "source_hashes") <- vapply(paths, .ab_file_hash, character(1L))
  class(result) <- c("ab_fixture_inputs", "list")
  result
}

# Internal normal-input helper. Unlike fixture mode, no answer table is required.
.ab_read_input_table <- function(file) {
  utils::read.delim(
    file, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
  )
}

.ab_read_sparse_counts <- function(input_dir) {
  tsv_path <- file.path(input_dir, "counts.tsv")
  rds_path <- file.path(input_dir, "counts.rds")
  present <- c(tsv = file.exists(tsv_path), rds = file.exists(rds_path))
  if (sum(present) != 1L) {
    stop(
      if (all(present)) {
        "input directory must contain exactly one of counts.tsv or counts.rds, not both"
      } else {
        "input directory is missing exactly one required count file: counts.tsv or counts.rds"
      },
      call. = FALSE
    )
  }
  if (present[["tsv"]]) {
    raw_counts <- .ab_read_input_table(tsv_path)
    if (!"gene_id" %in% names(raw_counts)) {
      stop("counts.tsv must begin with a gene_id column", call. = FALSE)
    }
    count_tokens <- raw_counts[, setdiff(names(raw_counts), "gene_id"), drop = FALSE]
    numeric_counts <- suppressWarnings(matrix(
      as.numeric(as.matrix(count_tokens)), nrow = nrow(count_tokens),
      ncol = ncol(count_tokens),
      dimnames = list(raw_counts$gene_id, names(count_tokens))
    ))
    counts <- methods::as(Matrix::Matrix(numeric_counts, sparse = TRUE), "dgCMatrix")
    return(list(counts = counts, path = tsv_path))
  }
  counts <- readRDS(rds_path)
  if (!inherits(counts, "sparseMatrix")) {
    stop("counts.rds must contain a sparse Matrix object", call. = FALSE)
  }
  counts <- tryCatch(
    methods::as(counts, "dgCMatrix"),
    error = function(e) stop(
      paste0("counts.rds cannot be converted to dgCMatrix: ", conditionMessage(e)),
      call. = FALSE
    )
  )
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts)))) {
    stop("counts.rds must have nonempty row and column names", call. = FALSE)
  }
  list(counts = counts, path = rds_path)
}

.ab_decode_logical_column <- function(table, field) {
  if (!field %in% names(table)) return(table)
  value <- table[[field]]
  normalized <- toupper(trimws(as.character(value)))
  if (anyNA(value) || any(!normalized %in% c("TRUE", "FALSE"))) {
    stop(paste0(field, " must contain only TRUE or FALSE"), call. = FALSE)
  }
  table[[field]] <- normalized == "TRUE"
  table
}

.ab_decode_named_constraint_map <- function(token, field) {
  if (length(token) != 1L || is.na(token) || !nzchar(trimws(token))) {
    stop(paste0(field, " must be a nonempty JSON object"), call. = FALSE)
  }
  value <- tryCatch(
    jsonlite::fromJSON(token, simplifyVector = FALSE),
    error = function(e) stop(
      paste0(field, " must be valid JSON: ", conditionMessage(e)), call. = FALSE
    )
  )
  value_names <- names(value)
  if (!is.list(value) || is.null(value_names) || !length(value) ||
      anyNA(value_names) || any(!nzchar(value_names)) || anyDuplicated(value_names)) {
    stop(
      paste0(field, " must be a JSON object with unique nonempty ID keys"),
      call. = FALSE
    )
  }
  scalar_numeric <- vapply(
    value,
    function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x),
    logical(1L)
  )
  if (!all(scalar_numeric)) {
    stop(paste0(field, " JSON values must be finite numbers"), call. = FALSE)
  }
  stats::setNames(as.numeric(unlist(value, use.names = FALSE)), value_names)
}

.ab_decode_input_constraints <- function(constraints) {
  if (!is.data.frame(constraints) || nrow(constraints) != 1L) return(constraints)
  fields <- c("max_changes_per_sample", "max_changes_per_donor")
  for (field in intersect(fields, names(constraints))) {
    token <- constraints[[field]][[1L]]
    if (is.character(token) && grepl("^\\s*\\{", token)) {
      constraints[[field]] <- list(.ab_decode_named_constraint_map(token, field))
    }
  }
  constraints
}

read_input_directory <- function(path, analysis_spec = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one nonempty directory path", call. = FALSE)
  }
  input_dir <- normalizePath(path, mustWork = TRUE)
  if (!dir.exists(input_dir)) stop("input path is not a directory", call. = FALSE)
  table_names <- c("cells", "specimens", "allowed_labels", "claim", "constraints")
  paths <- setNames(file.path(input_dir, paste0(table_names, ".tsv")), table_names)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(paste0("input directory is missing: ",
                paste(basename(missing), collapse = ", ")), call. = FALSE)
  }
  count_input <- .ab_read_sparse_counts(input_dir)
  paths <- c(counts = count_input$path, paths)
  cells <- .ab_decode_logical_column(
    .ab_read_input_table(paths[["cells"]]), "adt_ambiguity"
  )
  constraints <- .ab_decode_input_constraints(
    .ab_read_input_table(paths[["constraints"]])
  )
  result <- list(
    counts = count_input$counts,
    cells = cells,
    specimens = .ab_read_input_table(paths[["specimens"]]),
    allowed_labels = .ab_read_input_table(paths[["allowed_labels"]]),
    claim = .ab_read_input_table(paths[["claim"]]),
    constraints = constraints,
    analysis_spec = analysis_spec
  )
  attr(result, "input_dir") <- input_dir
  attr(result, "source_paths") <- paths
  attr(result, "source_hashes") <- vapply(paths, .ab_file_hash, character(1L))
  class(result) <- c("ab_file_inputs", "list")
  result
}
