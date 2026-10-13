test_fixture_root <- function() {
  system.file(
    "extdata", "test-fixtures",
    package = "AnnotationBreakpoint",
    mustWork = TRUE
  )
}

read_fixture_table <- function(case_id, filename) {
  path <- file.path(test_fixture_root(), case_id, filename)
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
}

load_raw_fixture <- function(case_id) {
  count_table <- read_fixture_table(case_id, "counts.tsv")
  count_values <- as.matrix(count_table[, -1L, drop = FALSE])
  storage.mode(count_values) <- "double"
  rownames(count_values) <- count_table$gene_id
  counts <- methods::as(Matrix::Matrix(count_values, sparse = TRUE), "dgCMatrix")

  claim <- as.list(read_fixture_table(case_id, "claim.tsv")[1L, , drop = FALSE])
  claim$q <- as.integer(claim$q)
  claim$delta <- as.numeric(claim$delta)
  claim$tolerance <- as.numeric(claim$tolerance)

  constraints <- as.list(read_fixture_table(case_id, "constraints.tsv")[1L, , drop = FALSE])
  constraints <- lapply(constraints, as.integer)

  expected <- read_fixture_table(case_id, "expected_results.tsv")
  expected$rank <- as.integer(expected$rank)
  expected$cost <- as.integer(expected$cost)
  expected$feasible <- expected$feasible == "true"
  expected$beta <- as.numeric(expected$beta)
  expected$overturn <- expected$overturn == "true"

  list(
    case_id = case_id,
    counts = counts,
    cells = read_fixture_table(case_id, "cells.tsv"),
    specimens = read_fixture_table(case_id, "specimens.tsv"),
    allowed_labels = read_fixture_table(case_id, "allowed_labels.tsv"),
    claim = claim,
    constraints = constraints,
    expected = expected
  )
}

validate_raw_fixture <- function(raw) {
  validate_inputs(
    raw$counts,
    raw$cells,
    raw$specimens,
    raw$allowed_labels,
    raw$claim,
    raw$constraints
  )
}

original_assignment <- function(raw) {
  stats::setNames(raw$cells$original_label, raw$cells$cell_id)
}

assignment_from_changes <- function(raw, changes) {
  assignment <- original_assignment(raw)
  if (!identical(changes, "original")) {
    tokens <- strsplit(changes, ";", fixed = TRUE)[[1L]]
    for (token in tokens) {
      parts <- strsplit(token, "=", fixed = TRUE)[[1L]]
      assignment[[parts[[1L]]]] <- parts[[2L]]
    }
  }
  assignment
}

independent_changes <- function(raw, assignment) {
  original <- original_assignment(raw)
  assignment <- assignment[names(original)]
  changed <- sort(names(original)[assignment != original])
  if (!length(changed)) return("original")
  paste(paste0(changed, "=", assignment[changed]), collapse = ";")
}

independent_beta <- function(raw, assignment = original_assignment(raw)) {
  expressions <- numeric(nrow(raw$specimens))
  names(expressions) <- raw$specimens$sample_id
  for (sample_id in raw$specimens$sample_id) {
    sample_cells <- raw$cells$cell_id[
      raw$cells$sample_id == sample_id & assignment[raw$cells$cell_id] == "A"
    ]
    totals <- Matrix::rowSums(raw$counts[, sample_cells, drop = FALSE])
    n <- unname(totals[["G_target"]])
    library <- sum(totals[c("G_target", "G_ref1", "G_ref2")])
    expressions[[sample_id]] <- log2(1e6 * (n + 0.5) / (library + 1.5))
  }
  case_samples <- raw$specimens$sample_id[raw$specimens$condition == "case"]
  control_samples <- raw$specimens$sample_id[raw$specimens$condition == "control"]
  mean(expressions[case_samples]) - mean(expressions[control_samples])
}

candidate_row <- function(result, changes) {
  rows <- result$assignments[result$assignments$changes == changes, , drop = FALSE]
  if (nrow(rows) != 1L) stop("candidate key is not unique: ", changes)
  rows
}

expect_invalid <- function(raw, pattern = NULL) {
  validation <- validate_raw_fixture(raw)
  testthat::expect_s3_class(validation, "ab_validation")
  testthat::expect_false(validation$ok)
  testthat::expect_identical(validation$status, "invalid_input")
  if (!is.null(pattern)) {
    testthat::expect_match(paste(validation$errors, collapse = " "), pattern, ignore.case = TRUE)
  }
  invisible(validation)
}

add_fixed_cell <- function(raw, cell_id, sample_id, label, feature_counts) {
  stopifnot(!cell_id %in% raw$cells$cell_id, identical(names(feature_counts), rownames(raw$counts)))
  specimen <- raw$specimens[raw$specimens$sample_id == sample_id, , drop = FALSE]
  stopifnot(nrow(specimen) == 1L)
  new_column <- matrix(
    as.numeric(feature_counts),
    ncol = 1L,
    dimnames = list(names(feature_counts), cell_id)
  )
  raw$counts <- methods::as(
    cbind(raw$counts, Matrix::Matrix(new_column, sparse = TRUE)),
    "dgCMatrix"
  )
  raw$cells <- rbind(raw$cells, data.frame(
    cell_id = cell_id,
    sample_id = sample_id,
    donor_id = specimen$donor_id,
    condition = specimen$condition,
    original_label = label,
    stringsAsFactors = FALSE
  ))
  raw$allowed_labels <- rbind(raw$allowed_labels, data.frame(
    cell_id = cell_id,
    allowed_label = label,
    evidence_source = "synthetic_fixture_design",
    evidence_version = "v1",
    stringsAsFactors = FALSE
  ))
  raw
}

rename_cell <- function(raw, old_id, new_id) {
  stopifnot(old_id %in% raw$cells$cell_id, !new_id %in% raw$cells$cell_id)
  colnames(raw$counts)[colnames(raw$counts) == old_id] <- new_id
  raw$cells$cell_id[raw$cells$cell_id == old_id] <- new_id
  raw$allowed_labels$cell_id[raw$allowed_labels$cell_id == old_id] <- new_id
  raw
}

drop_sample <- function(raw, sample_id) {
  retained_cells <- raw$cells$cell_id[raw$cells$sample_id != sample_id]
  raw$counts <- raw$counts[, retained_cells, drop = FALSE]
  raw$cells <- raw$cells[raw$cells$sample_id != sample_id, , drop = FALSE]
  raw$allowed_labels <- raw$allowed_labels[
    raw$allowed_labels$cell_id %in% retained_cells, , drop = FALSE
  ]
  raw$specimens <- raw$specimens[raw$specimens$sample_id != sample_id, , drop = FALSE]
  raw
}

failed_fit <- function(message = "forced test callback failure") {
  structure(list(
    status = "not_assessable",
    beta = NA_real_,
    margin = NA_real_,
    claim_holds = FALSE,
    overturn = FALSE,
    error = message
  ), class = "ab_fit")
}
