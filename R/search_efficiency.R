.ab_compact_search_fit <- function(fit) {
  keep <- c(
    "status", "beta", "margin", "claim_holds", "overturn", "error",
    "failure_reason", "warnings", "assignment_id", "analysis_id", "backend",
    "backend_version", "evaluation_identity", "fit_elapsed_seconds",
    "authoritative_rebuild", "search_elapsed_seconds", "historical_fit_index",
    "historical_proposal_index"
  )
  out <- fit[intersect(keep, names(fit))]
  if (is.null(out$error)) out$error <- NA_character_
  if (is.null(out$failure_reason)) out$failure_reason <- NA_character_
  if (is.null(out$warnings)) out$warnings <- character()
  structure(out, class = "ab_fit")
}

.ab_resume_fit <- function(value, assignment_id, evaluation_identity) {
  if (!is.list(value)) {
    stop("each resume evaluation must be a compact fit list", call. = FALSE)
  }
  if (!is.null(value$assignment_id) &&
      !identical(as.character(value$assignment_id), assignment_id)) {
    stop("resume evaluation assignment_id does not match its key", call. = FALSE)
  }
  status <- as.character(value$status)
  if (length(status) != 1L || is.na(status) ||
      !status %in% c("evaluable", "not_assessable")) {
    stop("resume evaluation has an unsupported status", call. = FALSE)
  }
  if (!is.null(value$evaluation_identity) &&
      !identical(value$evaluation_identity, evaluation_identity)) {
    stop("resume evaluation identity does not match the active fit", call. = FALSE)
  }
  if (identical(status, "evaluable")) {
    scalars <- c(value$beta, value$margin)
    if (length(scalars) != 2L || any(!is.finite(as.numeric(scalars))) ||
        length(value$claim_holds) != 1L || is.na(value$claim_holds) ||
        length(value$overturn) != 1L || is.na(value$overturn)) {
      stop("evaluable resume fit has invalid scalar evidence", call. = FALSE)
    }
  }
  value$assignment_id <- assignment_id
  value$evaluation_identity <- evaluation_identity
  .ab_compact_search_fit(structure(value, class = "ab_fit"))
}

.ab_single_move_source <- function(data, current, direction) {
  candidates <- list()
  cell_ids <- data$cells$cell_id
  for (cell_id in cell_ids) {
    alternatives <- setdiff(data$allowed_by_cell[[cell_id]], current[[cell_id]])
    if (identical(direction, "outgoing_only")) {
      alternatives <- alternatives[
        current[[cell_id]] == data$original_assignment[[cell_id]] &
          data$original_assignment[[cell_id]] == data$claim$focal_label &
          alternatives == data$claim$alternative_label
      ]
    }
    for (label in alternatives) {
      candidates[[length(candidates) + 1L]] <- list(
        cell_id = cell_id,
        from_label = unname(current[[cell_id]]),
        to_label = unname(label)
      )
    }
  }
  if (!length(candidates)) return(candidates)
  cell <- vapply(candidates, `[[`, character(1L), "cell_id")
  label <- vapply(candidates, `[[`, character(1L), "to_label")
  candidates <- candidates[order(cell, label, method = "radix")]
  for (index in seq_along(candidates)) {
    candidates[[index]]$source_rank <- as.integer(index)
  }
  candidates
}

.ab_randomized_single_source <- function(data, current, direction) {
  candidates <- .ab_single_move_source(data, current, direction)
  if (length(candidates) > 1L) {
    candidates <- candidates[sample.int(length(candidates))]
  }
  candidates
}

.ab_proxy_cell_cache <- function(data) {
  cell_ids <- data$cells$cell_id
  target_by_cell <- as.numeric(data$counts[
    data$claim$focal_gene, cell_ids, drop = TRUE
  ])
  library_by_cell <- as.numeric(Matrix::colSums(data$counts[
    data$claim$feature_universe, cell_ids, drop = FALSE
  ]))
  names(target_by_cell) <- cell_ids
  names(library_by_cell) <- cell_ids
  list(
    target_by_cell = target_by_cell,
    library_by_cell = library_by_cell,
    donor_by_cell = stats::setNames(data$cells$donor_id, cell_ids),
    donor_condition = stats::setNames(
      data$specimens$condition, data$specimens$donor_id
    ),
    donors_per_condition = table(data$specimens$condition)
  )
}

.ab_proxy_ranked_source <- function(data, current, direction, proxy_cache,
                                    target_pseudocount = 0.5,
                                    library_pseudocount = 1,
                                    scale = 1e6) {
  candidates <- .ab_single_move_source(data, current, direction)
  if (!length(candidates)) return(candidates)

  cell_ids <- data$cells$cell_id
  target_by_cell <- proxy_cache$target_by_cell
  library_by_cell <- proxy_cache$library_by_cell
  focal <- data$claim$focal_label
  donor_ids <- unname(proxy_cache$donor_by_cell[cell_ids])
  donor_levels <- data$donor_ids
  focal_cells <- unname(current[cell_ids]) == focal
  target_by_donor <- vapply(donor_levels, function(donor_id) {
    sum(target_by_cell[donor_ids == donor_id & focal_cells])
  }, numeric(1L))
  library_by_donor <- vapply(donor_levels, function(donor_id) {
    sum(library_by_cell[donor_ids == donor_id & focal_cells])
  }, numeric(1L))
  names(target_by_donor) <- donor_levels
  names(library_by_donor) <- donor_levels
  donor_condition <- proxy_cache$donor_condition
  donors_per_condition <- proxy_cache$donors_per_condition

  proxy_value <- function(candidate) {
    cell_id <- candidate$cell_id
    donor_id <- unname(proxy_cache$donor_by_cell[[cell_id]])
    condition <- unname(donor_condition[[donor_id]])
    condition_sign <- if (identical(condition, data$claim$case_condition)) {
      1
    } else if (identical(condition, data$claim$control_condition)) {
      -1
    } else {
      0
    }
    target_delta <- 0
    library_delta <- 0
    if (!identical(candidate$from_label, focal) &&
        identical(candidate$to_label, focal)) {
      target_delta <- target_by_cell[[cell_id]]
      library_delta <- library_by_cell[[cell_id]]
    } else if (identical(candidate$from_label, focal) &&
               !identical(candidate$to_label, focal)) {
      target_delta <- -target_by_cell[[cell_id]]
      library_delta <- -library_by_cell[[cell_id]]
    }
    before <- log2(
      (target_by_donor[[donor_id]] + target_pseudocount) /
        (library_by_donor[[donor_id]] + library_pseudocount) * scale
    )
    after <- log2(
      (target_by_donor[[donor_id]] + target_delta + target_pseudocount) /
        (library_by_donor[[donor_id]] + library_delta +
           library_pseudocount) * scale
    )
    as.numeric(
      data$claim$q * condition_sign * (after - before) /
        unname(donors_per_condition[[condition]])
    )
  }
  scores <- vapply(candidates, proxy_value, numeric(1L))
  candidate_cells <- vapply(candidates, `[[`, character(1L), "cell_id")
  candidate_labels <- vapply(candidates, `[[`, character(1L), "to_label")
  ranked <- order(scores, candidate_cells, candidate_labels, method = "radix")
  candidates <- candidates[ranked]
  scores <- scores[ranked]
  for (index in seq_along(candidates)) {
    candidates[[index]]$proxy_score <- as.numeric(scores[[index]])
    candidates[[index]]$ranking_rank <- as.integer(index)
  }
  candidates
}
