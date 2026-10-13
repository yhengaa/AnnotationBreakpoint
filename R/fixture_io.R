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
read_input_directory <- function(path, analysis_spec = NULL) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be one nonempty directory path", call. = FALSE)
  }
  input_dir <- normalizePath(path, mustWork = TRUE)
  if (!dir.exists(input_dir)) stop("input path is not a directory", call. = FALSE)
  table_names <- c(
    "counts", "cells", "specimens", "allowed_labels", "claim", "constraints"
  )
  paths <- setNames(file.path(input_dir, paste0(table_names, ".tsv")), table_names)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(paste0("input directory is missing: ",
                paste(basename(missing), collapse = ", ")), call. = FALSE)
  }
  read_table <- function(file) utils::read.delim(
    file, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
  )
  raw_counts <- read_table(paths[["counts"]])
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
  result <- list(
    counts = counts,
    cells = read_table(paths[["cells"]]),
    specimens = read_table(paths[["specimens"]]),
    allowed_labels = read_table(paths[["allowed_labels"]]),
    claim = read_table(paths[["claim"]]),
    constraints = read_table(paths[["constraints"]]),
    analysis_spec = analysis_spec
  )
  attr(result, "input_dir") <- input_dir
  attr(result, "source_paths") <- paths
  attr(result, "source_hashes") <- vapply(paths, .ab_file_hash, character(1L))
  class(result) <- c("ab_file_inputs", "list")
  result
}
