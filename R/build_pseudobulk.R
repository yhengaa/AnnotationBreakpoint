.ab_pseudobulk_columns <- function(data) {
  labels <- .ab_unique_sorted(unique(c(
    data$cells$original_label,
    data$allowed_labels$allowed_label,
    data$registered_labels
  )))
  rows <- vector("list", length(data$sample_ids) * length(labels))
  at <- 1L
  for (sample_id in data$sample_ids) {
    specimen <- data$specimens[data$specimens$sample_id == sample_id, , drop = FALSE]
    for (label in labels) {
      rows[[at]] <- data.frame(
        pseudobulk_id = .ab_pair_id(sample_id, label),
        display_id = paste0(sample_id, "::", label),
        sample_id = sample_id,
        donor_id = specimen$donor_id[[1L]],
        condition = specimen$condition[[1L]],
        label = label,
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  do.call(rbind, rows)
}

.ab_build_full_counts <- function(data, assignment, columns) {
  groups <- .ab_pair_id(
    data$cells$sample_id,
    unname(assignment[data$cells$cell_id])
  )
  group_index <- match(groups, columns$pseudobulk_id)
  if (anyNA(group_index)) {
    stop("internal error: assignment refers to an unknown pseudobulk column", call. = FALSE)
  }
  incidence <- Matrix::sparseMatrix(
    i = seq_along(group_index),
    j = group_index,
    x = rep.int(1, length(group_index)),
    dims = c(length(group_index), nrow(columns)),
    dimnames = list(data$cells$cell_id, columns$pseudobulk_id)
  )
  out <- data$counts[, data$cells$cell_id, drop = FALSE] %*% incidence
  rownames(out) <- rownames(data$counts)
  colnames(out) <- columns$pseudobulk_id
  if (!inherits(out, "dgCMatrix")) out <- methods::as(out, "dgCMatrix")
  Matrix::drop0(out)
}

.ab_conservation <- function(data, pseudobulk_counts, columns) {
  rows <- vector("list", length(data$sample_ids))
  for (idx in seq_along(data$sample_ids)) {
    sample_id <- data$sample_ids[[idx]]
    input_ids <- data$cells$cell_id[data$cells$sample_id == sample_id]
    output_cols <- columns$pseudobulk_id[columns$sample_id == sample_id]
    input_total <- Matrix::rowSums(data$counts[, input_ids, drop = FALSE])
    pseudobulk_total <- Matrix::rowSums(
      pseudobulk_counts[, output_cols, drop = FALSE]
    )
    rows[[idx]] <- data.frame(
      feature_id = rownames(data$counts),
      sample_id = sample_id,
      input_total = as.numeric(input_total),
      pseudobulk_total = as.numeric(pseudobulk_total),
      conserved = as.numeric(input_total) == as.numeric(pseudobulk_total),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.ab_normalize_moves <- function(data, base_assignment, moves) {
  if (is.null(moves)) {
    stop("incremental mode requires explicit moves or a target assignment", call. = FALSE)
  }
  if (is.data.frame(moves)) {
    if (!"cell_id" %in% names(moves)) {
      stop("moves data frame must contain cell_id", call. = FALSE)
    }
    target_field <- intersect(c("new_label", "assigned_label", "to_label"), names(moves))
    if (length(target_field) != 1L) {
      stop(
        "moves must contain exactly one of new_label, assigned_label, or to_label",
        call. = FALSE
      )
    }
    cell_id <- as.character(moves$cell_id)
    new_label <- as.character(moves[[target_field]])
    if ("from_label" %in% names(moves)) {
      from_label <- as.character(moves$from_label)
      current <- unname(base_assignment[cell_id])
      if (anyNA(current) || !identical(from_label, current)) {
        stop("moves$from_label does not match the base assignment", call. = FALSE)
      }
    }
  } else if (is.atomic(moves) && !is.null(names(moves))) {
    cell_id <- names(moves)
    new_label <- as.character(moves)
  } else if (is.list(moves) &&
             all(c("cell_id", "new_label") %in% names(moves))) {
    cell_id <- as.character(moves$cell_id)
    new_label <- as.character(moves$new_label)
  } else {
    stop(
      "moves must be a named label vector or a data frame with cell_id and new_label",
      call. = FALSE
    )
  }
  if (length(cell_id) != length(new_label) || anyNA(cell_id) || anyNA(new_label) ||
      any(!nzchar(cell_id)) || any(!nzchar(new_label))) {
    stop("moves contains missing or empty cell IDs or labels", call. = FALSE)
  }
  if (anyDuplicated(cell_id)) {
    stop("moves contains conflicting duplicate cell IDs", call. = FALSE)
  }
  if (any(!cell_id %in% data$cells$cell_id)) {
    stop("moves contains unknown cell IDs", call. = FALSE)
  }
  data.frame(
    cell_id = cell_id,
    from_label = unname(base_assignment[cell_id]),
    new_label = new_label,
    stringsAsFactors = FALSE
  )
}

#' Construct sparse sample-by-label pseudobulks
#'
#' @param x A successful `ab_validation` object.
#' @param assignment A named complete assignment, or `NULL` for original labels.
#' @param method Either full reconstruction or an incremental update.
#' @param base An `ab_pseudobulk` to update in incremental mode.
#' @param moves Explicit cell-to-label moves in incremental mode.
#' @return An object of class `ab_pseudobulk`.
#' @export
build_pseudobulk <- function(x, assignment = NULL,
                             method = c("full", "incremental"),
                             base = NULL, moves = NULL) {
  method <- match.arg(method)
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    stop("x must be a successful ab_validation object", call. = FALSE)
  }
  data <- x$data
  current_hashes <- .ab_recompute_input_hashes(x)
  columns <- .ab_pseudobulk_columns(data)

  if (method == "full") {
    assignment <- .ab_normalize_assignment(data, assignment)
    aggregate <- .ab_build_full_counts(data, assignment, columns)
    applied_moves <- .ab_change_rows(data, assignment)
  } else {
    if (!inherits(base, "ab_pseudobulk")) {
      stop("base must be an ab_pseudobulk in incremental mode", call. = FALSE)
    }
    if (!identical(unname(base$input_hashes), unname(current_hashes)) ||
        !identical(names(base$input_hashes), names(current_hashes))) {
      stop("base pseudobulk was built from different immutable inputs", call. = FALSE)
    }
    base_assignment <- .ab_normalize_assignment(data, base$assignment)
    if (!is.null(assignment) && !is.null(moves)) {
      stop("provide either assignment or moves in incremental mode, not both", call. = FALSE)
    }
    if (!is.null(assignment)) {
      target_assignment <- .ab_normalize_assignment(data, assignment)
      changed_ids <- names(target_assignment)[target_assignment != base_assignment]
      moves <- setNames(unname(target_assignment[changed_ids]), changed_ids)
    }
    normalized_moves <- .ab_normalize_moves(data, base_assignment, moves)
    assignment <- base_assignment
    assignment[normalized_moves$cell_id] <- normalized_moves$new_label
    assignment <- .ab_normalize_assignment(data, assignment)

    base_columns <- base$column_metadata[, names(columns), drop = FALSE]
    if (!identical(base_columns, columns)) {
      stop("base pseudobulk column metadata does not match current inputs", call. = FALSE)
    }
    aggregate <- base$counts
    for (idx in seq_len(nrow(normalized_moves))) {
      move <- normalized_moves[idx, , drop = FALSE]
      if (move$from_label[[1L]] == move$new_label[[1L]]) next
      cell_id <- move$cell_id[[1L]]
      sample_id <- data$cells$sample_id[match(cell_id, data$cells$cell_id)]
      from_column <- .ab_pair_id(sample_id, move$from_label[[1L]])
      to_column <- .ab_pair_id(sample_id, move$new_label[[1L]])
      cell_counts <- data$counts[, cell_id, drop = FALSE]
      aggregate[, from_column] <- aggregate[, from_column, drop = FALSE] - cell_counts
      aggregate[, to_column] <- aggregate[, to_column, drop = FALSE] + cell_counts
    }
    aggregate <- Matrix::drop0(methods::as(aggregate, "dgCMatrix"))
    if (any(aggregate@x < 0)) {
      stop("incremental update produced negative pseudobulk counts", call. = FALSE)
    }
    independent_full <- .ab_build_full_counts(data, assignment, columns)
    if (!isTRUE(all.equal(aggregate, independent_full, tolerance = 0,
                          check.attributes = TRUE))) {
      stop("incremental update does not equal full reconstruction", call. = FALSE)
    }
    applied_moves <- normalized_moves
  }

  conservation <- .ab_conservation(data, aggregate, columns)
  all_conserved <- nrow(conservation) > 0L && all(conservation$conserved)
  if (!all_conserved) {
    stop("pseudobulk construction failed sample/feature conservation", call. = FALSE)
  }
  column_cells <- vapply(
    seq_len(nrow(columns)),
    function(idx) {
      sum(data$cells$sample_id == columns$sample_id[[idx]] &
            assignment == columns$label[[idx]])
    },
    integer(1L)
  )
  columns$cell_count <- unname(column_cells)
  rownames(columns) <- NULL

  structure(list(
    counts = aggregate,
    column_metadata = columns,
    columns = columns,
    assignment = assignment,
    method = method,
    moves = applied_moves,
    conservation = conservation,
    all_conserved = all_conserved,
    input_hashes = current_hashes,
    assignment_id = .ab_assignment_id(data, assignment)
  ), class = "ab_pseudobulk")
}
