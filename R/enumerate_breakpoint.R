.ab_call_fit_fun <- function(fit_fun, x, assignment) {
  formal_names <- tryCatch(names(formals(fit_fun)), error = function(e) NULL)
  arguments <- list(x = x, assignment = assignment)
  if (!is.null(formal_names) && !"..." %in% formal_names) {
    arguments <- arguments[intersect(names(arguments), formal_names)]
  }
  if (!length(arguments) && !is.null(formal_names) && length(formal_names) == 1L) {
    arguments <- list(assignment)
  }
  do.call(fit_fun, arguments)
}

.ab_failed_evaluation <- function(message, assignment, assignment_id,
                                  failure_reason = "fit_evaluation_failure",
                                  warnings = character()) {
  structure(list(
    status = "not_assessable",
    beta = NA_real_,
    margin = NA_real_,
    claim_holds = NA,
    overturn = FALSE,
    sample_expression = data.frame(),
    per_sample = data.frame(),
    pseudobulk_summary = data.frame(),
    error = as.character(message),
    failure_reason = as.character(failure_reason),
    warnings = unique(as.character(warnings)),
    assignment = assignment,
    assignment_id = assignment_id,
    pseudobulk = NULL
  ), class = "ab_fit")
}

.ab_standardize_fit <- function(value, data, assignment, assignment_id) {
  if (inherits(value, "error")) {
    return(.ab_failed_evaluation(
      conditionMessage(value), assignment, assignment_id,
      failure_reason = "fit_callback_failure"
    ))
  }
  if (inherits(value, "ab_fit") || is.list(value)) {
    if (!is.null(value$status) && !identical(value$status, "evaluable")) {
      if (inherits(value, "ab_fit")) {
        value$beta <- NA_real_
        value$margin <- NA_real_
        value$claim_holds <- NA
        value$overturn <- FALSE
        value$assignment <- assignment
        value$assignment_id <- assignment_id
        return(value)
      }
      message <- if (!is.null(value$error) && !is.na(value$error[[1L]])) {
        value$error[[1L]]
      } else {
        paste0("fit returned status ", value$status[[1L]])
      }
      reason <- if (!is.null(value$failure_reason) &&
                    length(value$failure_reason) == 1L &&
                    !is.na(value$failure_reason)) {
        as.character(value$failure_reason)
      } else "fit_reported_not_assessable"
      fit_warnings <- if (!is.null(value$warnings)) value$warnings else character()
      return(.ab_failed_evaluation(
        message, assignment, assignment_id, reason, fit_warnings
      ))
    }
    beta <- value$beta
  } else {
    beta <- value
  }
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta)) {
    return(.ab_failed_evaluation(
      "fit did not return one finite numeric beta", assignment, assignment_id,
      failure_reason = "invalid_fit_result"
    ))
  }
  margin <- data$claim$q * as.numeric(beta) - data$claim$delta
  if (!is.finite(margin)) {
    return(.ab_failed_evaluation(
      "fit produced a nonfinite registered margin", assignment, assignment_id,
      failure_reason = "nonfinite_registered_margin"
    ))
  }
  if (inherits(value, "ab_fit")) {
    value$beta <- as.numeric(beta)
    value$margin <- as.numeric(margin)
    value$claim_holds <- isTRUE(margin > data$claim$tolerance)
    value$overturn <- isTRUE(margin <= data$claim$tolerance)
    value$error <- NA_character_
    value$assignment <- assignment
    value$assignment_id <- assignment_id
    value$status <- "evaluable"
    return(value)
  }
  structure(list(
    status = "evaluable",
    beta = as.numeric(beta),
    margin = as.numeric(margin),
    claim_holds = isTRUE(margin > data$claim$tolerance),
    overturn = isTRUE(margin <= data$claim$tolerance),
    sample_expression = if (is.list(value) && !is.null(value$sample_expression)) {
      value$sample_expression
    } else {
      data.frame()
    },
    per_sample = if (is.list(value) && !is.null(value$per_sample)) {
      value$per_sample
    } else {
      data.frame()
    },
    pseudobulk_summary = if (is.list(value) && !is.null(value$pseudobulk_summary)) {
      value$pseudobulk_summary
    } else {
      data.frame()
    },
    error = NA_character_,
    failure_reason = NA_character_,
    warnings = if (is.list(value) && !is.null(value$warnings)) {
      unique(as.character(value$warnings))
    } else character(),
    assignment = assignment,
    assignment_id = assignment_id,
    pseudobulk = if (is.list(value)) value$pseudobulk else NULL
  ), class = "ab_fit")
}

.ab_candidate_domain <- function(data) {
  cell_ids <- data$cells$cell_id
  choices <- data$allowed_by_cell[cell_ids]
  domain_size <- prod(vapply(choices, length, numeric(1L)))
  ambiguous_count <- sum(vapply(choices, length, integer(1L)) > 1L)
  if (!is.finite(domain_size) || domain_size > 4096 || ambiguous_count > 12L) {
    stop("exact_reference_capacity_exceeded: maximum 12 ambiguous cells and 4096 assignments",
         call. = FALSE)
  }
  grid <- expand.grid(
    unname(choices),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  names(grid) <- cell_ids
  assignments <- vector("list", nrow(grid))
  change_text <- character(nrow(grid))
  cost <- integer(nrow(grid))
  assignment_id <- character(nrow(grid))
  failures <- vector("list", nrow(grid))
  for (idx in seq_len(nrow(grid))) {
    assignment <- setNames(
      vapply(grid[idx, , drop = FALSE], as.character, character(1L)),
      cell_ids
    )
    assignment <- .ab_normalize_assignment(data, assignment)
    changes <- .ab_change_rows(data, assignment)
    assignments[[idx]] <- assignment
    change_text[[idx]] <- .ab_changes_text(changes)
    cost[[idx]] <- nrow(changes)
    assignment_id[[idx]] <- .ab_assignment_id(data, assignment)
    failures[[idx]] <- .ab_constraint_failures(data, assignment)
  }
  candidate_order <- order(cost, change_text, assignment_id, method = "radix")
  assignments <- assignments[candidate_order]
  failures <- failures[candidate_order]
  table <- data.frame(
    rank = seq_along(candidate_order) - 1L,
    assignment_id = assignment_id[candidate_order],
    changes = change_text[candidate_order],
    cost = cost[candidate_order],
    feasible = !lengths(failures),
    constraint_failures = vapply(
      failures,
      function(value) paste(value, collapse = ";"),
      character(1L)
    ),
    evaluated = FALSE,
    evaluable = FALSE,
    beta = NA_real_,
    margin = NA_real_,
    claim_holds = NA,
    overturn = FALSE,
    fit_elapsed_seconds = NA_real_,
    warnings = NA_character_,
    failure_reason = NA_character_,
    evaluation_error = NA_character_,
    stringsAsFactors = FALSE
  )
  table$assignment <- I(assignments)
  list(table = table, assignments = assignments)
}

.ab_effects_by_cost <- function(candidates) {
  if (!nrow(candidates)) return(data.frame())
  rows <- lapply(split(seq_len(nrow(candidates)), candidates$cost), function(index) {
    evaluable <- candidates$evaluable[index]
    betas <- candidates$beta[index][evaluable]
    data.frame(
      cost = candidates$cost[index[[1L]]],
      assignment_count = length(index),
      feasible_count = sum(candidates$feasible[index]),
      infeasible_count = sum(!candidates$feasible[index]),
      evaluated_count = sum(candidates$evaluated[index]),
      evaluable_count = sum(evaluable),
      unevaluable_count = sum(
        candidates$evaluated[index] & !candidates$evaluable[index]
      ),
      overturn_count = sum(candidates$overturn[index]),
      minimum_beta = if (length(betas)) min(betas) else NA_real_,
      maximum_beta = if (length(betas)) max(betas) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$cost), , drop = FALSE]
}

.ab_fill_result_accounting <- function(result, candidates) {
  result$assignments <- candidates
  result$candidates <- candidates
  result$effects_by_cost <- .ab_effects_by_cost(candidates)
  result$total_assignment_count <- nrow(candidates)
  result$feasible_assignment_count <- sum(candidates$feasible)
  result$infeasible_assignment_count <- sum(!candidates$feasible)
  result$evaluable_assignment_count <- sum(candidates$evaluable)
  result$unevaluable_assignment_count <- sum(
    candidates$feasible & candidates$evaluated & !candidates$evaluable
  )
  result$unevaluated_assignment_count <- sum(
    candidates$feasible & !candidates$evaluated
  )
  result$overturn_assignment_count <- sum(candidates$overturn)
  result$search_complete <- !any(candidates$feasible & !candidates$evaluated)
  result
}

.ab_record_fit <- function(candidates, index, fit) {
  candidates$evaluated[[index]] <- TRUE
  candidates$evaluable[[index]] <- identical(fit$status, "evaluable")
  candidates$fit_elapsed_seconds[[index]] <- if (
    !is.null(fit$fit_elapsed_seconds) && length(fit$fit_elapsed_seconds) == 1L
  ) as.numeric(fit$fit_elapsed_seconds) else NA_real_
  candidates$warnings[[index]] <- if (!is.null(fit$warnings) &&
                                      length(fit$warnings)) {
    paste(unique(as.character(fit$warnings)), collapse = " | ")
  } else NA_character_
  candidates$failure_reason[[index]] <- if (!is.null(fit$failure_reason) &&
                                             length(fit$failure_reason) == 1L &&
                                             !is.na(fit$failure_reason)) {
    as.character(fit$failure_reason)
  } else NA_character_
  if (isTRUE(candidates$evaluable[[index]])) {
    candidates$beta[[index]] <- fit$beta
    candidates$margin[[index]] <- fit$margin
    candidates$claim_holds[[index]] <- fit$claim_holds
    candidates$overturn[[index]] <- fit$overturn
    candidates$evaluation_error[[index]] <- NA_character_
  } else {
    candidates$beta[[index]] <- NA_real_
    candidates$margin[[index]] <- NA_real_
    candidates$claim_holds[[index]] <- NA
    candidates$overturn[[index]] <- FALSE
    candidates$evaluation_error[[index]] <- fit$error
  }
  candidates
}

#' Exhaustively enumerate the permitted finite assignment domain
#'
#' @param x An `ab_validation` object.
#' @param fit_fun Fit function called once per attempted feasible assignment.
#' @param max_evaluations Maximum feasible fits to attempt, including original.
#' @return An object of class `ab_result`.
#' @export
enumerate_breakpoint <- function(x, fit_fun = NULL, max_evaluations = Inf,
                                 fit_identity = NULL) {
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    result <- .ab_result_skeleton(
      "invalid_input", "input_validation_failed",
      if (inherits(x, "ab_validation")) x$input_hashes else character()
    )
    result$errors <- if (inherits(x, "ab_validation")) {
      x$errors
    } else {
      "x must be an ab_validation object"
    }
    result$max_evaluations <- max_evaluations
    return(result)
  }
  if (length(max_evaluations) != 1L || is.na(max_evaluations) ||
      (is.infinite(max_evaluations) && max_evaluations < 0) ||
      (!is.infinite(max_evaluations) &&
       (!.ab_is_whole(max_evaluations) || max_evaluations < 1))) {
    result <- .ab_result_skeleton(
      "invalid_input", "invalid_max_evaluations", .ab_recompute_input_hashes(x)
    )
    result$errors <- "max_evaluations must be Inf or a positive integer"
    result$max_evaluations <- max_evaluations
    return(.ab_attach_contract(result, x$data))
  }

  data <- x$data
  current_hashes <- .ab_recompute_input_hashes(x)
  data$input_hashes <- current_hashes
  resolved_fit <- tryCatch(
    .ab_resolve_fit(x, fit_fun = fit_fun, fit_identity = fit_identity),
    error = function(error) error
  )
  if (inherits(resolved_fit, "error")) {
    result <- .ab_result_skeleton(
      "invalid_input", "invalid_fit_function", current_hashes
    )
    result$errors <- conditionMessage(resolved_fit)
    result$max_evaluations <- max_evaluations
    return(.ab_attach_contract(result, data))
  }
  fit_fun <- resolved_fit$fun
  finish <- function(result) {
    .ab_attach_contract(result, data, resolved_fit$identity)
  }
  domain <- tryCatch(.ab_candidate_domain(data), error = function(error) error)
  if (inherits(domain, "error")) {
    result <- .ab_result_skeleton(
      "not_assessable", "enumeration_capacity_exceeded", current_hashes
    )
    result$errors <- conditionMessage(domain)
    result$max_evaluations <- max_evaluations
    return(finish(result))
  }
  candidates <- domain$table
  original_index <- which(candidates$cost == 0L)
  if (length(original_index) != 1L || !candidates$feasible[[original_index]]) {
    result <- .ab_result_skeleton(
      "invalid_input", "original_assignment_infeasible", current_hashes
    )
    result$errors <- "the required original assignment is not uniquely feasible"
    result$max_evaluations <- max_evaluations
    result <- .ab_fill_result_accounting(result, candidates)
    return(finish(result))
  }

  evaluate_at <- function(index) {
    value <- tryCatch(
      .ab_call_fit_fun(fit_fun, x, candidates$assignment[[index]]),
      error = function(error) error
    )
    .ab_standardize_fit(
      value,
      data,
      candidates$assignment[[index]],
      candidates$assignment_id[[index]]
    )
  }

  original_fit <- evaluate_at(original_index)
  candidates <- .ab_record_fit(candidates, original_index, original_fit)
  if (!identical(original_fit$status, "evaluable")) {
    result <- .ab_result_skeleton(
      "not_assessable", "original_fit_not_assessable", current_hashes
    )
    result$original_fit <- original_fit
    result$errors <- original_fit$error
    result$max_evaluations <- max_evaluations
    result <- .ab_fill_result_accounting(result, candidates)
    return(finish(result))
  }
  if (!isTRUE(original_fit$claim_holds)) {
    result <- .ab_result_skeleton(
      "original_claim_not_met", "original_claim_not_met", current_hashes
    )
    result$original_fit <- original_fit
    result$max_evaluations <- max_evaluations
    result <- .ab_fill_result_accounting(result, candidates)
    return(finish(result))
  }

  evaluations_used <- 1L
  for (index in seq_len(nrow(candidates))) {
    if (index == original_index || !candidates$feasible[[index]]) next
    if (evaluations_used >= max_evaluations) break
    fit <- evaluate_at(index)
    candidates <- .ab_record_fit(candidates, index, fit)
    evaluations_used <- evaluations_used + 1L
  }

  result <- .ab_result_skeleton("no_counterexample_found", "", current_hashes)
  result$original_fit <- original_fit
  result$max_evaluations <- max_evaluations
  result <- .ab_fill_result_accounting(result, candidates)
  counterexample_index <- which(
    candidates$feasible & candidates$evaluable & candidates$overturn
  )
  failed_index <- which(
    candidates$feasible & candidates$evaluated & !candidates$evaluable
  )

  if (length(counterexample_index)) {
    best_cost <- min(candidates$cost[counterexample_index])
    minimum_index <- counterexample_index[
      candidates$cost[counterexample_index] == best_cost
    ]
    witness_index <- minimum_index[[1L]]
    witness_assignment <- candidates$assignment[[witness_index]]
    witness <- list(
      witness_id = candidates$assignment_id[[witness_index]],
      assignment_id = candidates$assignment_id[[witness_index]],
      assignment = witness_assignment,
      changes = candidates$changes[[witness_index]],
      changed_cells = .ab_changed_with_provenance(data, witness_assignment),
      cost = as.integer(best_cost),
      beta = candidates$beta[[witness_index]],
      margin = candidates$margin[[witness_index]],
      overturn = TRUE,
      input_hashes = current_hashes
    )
    result$best_verified_cost <- as.integer(best_cost)
    result$verified_upper_bound <- as.integer(best_cost)
    result$witness <- witness
    result$witness_id <- witness$witness_id
    result$minimum_cost_overturns <- candidates[minimum_index, , drop = FALSE]

    if (!result$search_complete) {
      result$status <- "verified_counterexample"
      result$termination_reason <-
        "evaluation_budget_exhausted_with_verified_counterexample"
    } else if (any(candidates$cost[failed_index] < best_cost)) {
      result$status <- "exhaustive_search_incomplete"
      result$termination_reason <-
        "cheaper_feasible_evaluation_failure_blocks_exact_minimum"
    } else {
      result$status <- "exact_minimum"
      result$exact_minimum_cost <- as.integer(best_cost)
      result$termination_reason <- "finite_domain_exhausted_exact_minimum"
    }
  } else if (!result$search_complete) {
    result$status <- "no_counterexample_found"
    result$termination_reason <-
      "evaluation_budget_exhausted_without_counterexample"
  } else if (length(failed_index)) {
    result$status <- "exhaustive_search_incomplete"
    result$termination_reason <-
      "feasible_evaluation_failures_block_exhaustive_stability"
  } else {
    result$status <- "exhaustively_stable"
    result$termination_reason <- "finite_domain_exhausted_without_counterexample"
  }
  finish(result)
}
