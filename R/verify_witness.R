.ab_compare_hashes <- function(expected, observed) {
  keys <- .ab_sort(unique(c(names(expected), names(observed))))
  expected_values <- as.character(expected[keys])
  observed_values <- as.character(observed[keys])
  data.frame(
    hash_name = keys,
    expected = unname(expected_values),
    observed = unname(observed_values),
    matches = !is.na(expected_values) & !is.na(observed_values) &
      expected_values == observed_values,
    stringsAsFactors = FALSE
  )
}

.ab_empty_changed_cells <- function(data) {
  .ab_changed_with_provenance(data, data$original_assignment)
}

.ab_compare_changed_cells <- function(provided, expected) {
  if (!is.data.frame(provided) || !all(names(expected) %in% names(provided))) {
    return(FALSE)
  }
  provided <- provided[, names(expected), drop = FALSE]
  if (nrow(provided)) {
    provided <- provided[order(provided$cell_id, method = "radix"), , drop = FALSE]
  }
  rownames(provided) <- NULL
  rownames(expected) <- NULL
  identical(provided, expected)
}

.ab_result_domain_matches <- function(result, domain) {
  candidates <- result$candidates
  if (!is.data.frame(candidates) || nrow(candidates) != nrow(domain$table)) {
    return(FALSE)
  }
  required <- c("assignment_id", "cost", "feasible", "constraint_failures")
  if (!all(required %in% names(candidates))) return(FALSE)
  identical(as.character(candidates$assignment_id), domain$table$assignment_id) &&
    identical(as.integer(candidates$cost), domain$table$cost) &&
    identical(as.logical(candidates$feasible), domain$table$feasible) &&
    identical(as.character(candidates$constraint_failures),
              domain$table$constraint_failures)
}

.ab_verify_certified_status <- function(x, result, domain, replay_fit,
                                        witness_cost = NULL) {
  if (!result$status %in% c("exact_minimum", "exhaustively_stable")) {
    return(character())
  }
  errors <- character()
  if (!isTRUE(result$search_complete)) {
    errors <- c(errors, "certified result is marked as an incomplete search")
  }
  if (!.ab_result_domain_matches(result, domain)) {
    return(c(errors, "result candidate domain does not match immutable inputs"))
  }
  replay_indices <- which(domain$table$feasible)
  if (identical(result$status, "exact_minimum")) {
    witness_index <- match(result$witness_id, domain$table$assignment_id)
    replay_indices <- unique(c(
      replay_indices[domain$table$cost[replay_indices] < witness_cost],
      witness_index
    ))
    replay_indices <- replay_indices[!is.na(replay_indices)]
  }
  replay_fits <- lapply(replay_indices, function(index) {
    replay_fit(domain$table$assignment[[index]])
  })
  replay_evaluable <- vapply(
    replay_fits, function(fit) identical(fit$status, "evaluable"), logical(1L)
  )
  if (!all(replay_evaluable)) {
    return(c(errors,
             "a feasible assignment required for certification failed replay"))
  }
  replay_overturn <- vapply(replay_fits, function(fit) fit$overturn, logical(1L))
  replay_beta <- vapply(replay_fits, function(fit) fit$beta, numeric(1L))
  recorded_index <- match(
    domain$table$assignment_id[replay_indices], result$candidates$assignment_id
  )
  recorded_beta <- result$candidates$beta[recorded_index]
  if (any(!is.finite(recorded_beta)) ||
      any(abs(recorded_beta - replay_beta) > x$data$claim$tolerance)) {
    errors <- c(errors, "recorded candidate effects do not match immutable replay")
  }
  if (identical(result$status, "exhaustively_stable")) {
    if (any(replay_overturn)) {
      errors <- c(errors, "stable result contains an overturn on full replay")
    }
    if (!is.null(result$witness) || !is.na(result$best_verified_cost) ||
        !is.na(result$exact_minimum_cost)) {
      errors <- c(errors, "stable result must not contain a witness or breakpoint cost")
    }
  } else {
    replay_costs <- domain$table$cost[replay_indices][replay_overturn]
    if (!length(replay_costs) || min(replay_costs) != witness_cost) {
      errors <- c(errors, "claimed exact minimum does not match domain replay")
    }
    if (length(result$exact_minimum_cost) != 1L ||
        is.na(result$exact_minimum_cost) ||
        result$exact_minimum_cost != witness_cost) {
      errors <- c(errors, "result exact_minimum_cost does not match witness cost")
    }
  }
  errors
}

.ab_verification_result <- function(ok, errors, status, before, after,
                                    changed_cells, hashes, hash_matches,
                                    witness_id, result_identity,
                                    fit_count, beta = NA_real_, cost = NA_integer_,
                                    evaluation_identity = NULL, started = NULL) {
  structure(list(
    ok = ok,
    errors = unique(errors),
    recomputed_cost = cost,
    recomputed_effect = beta,
    recomputed_beta = beta,
    recomputed_status = if (ok) status else "verification_failed",
    before_pseudobulks = before,
    after_pseudobulks = after,
    before_donor_pseudobulks = before,
    after_donor_pseudobulks = after,
    changed_cells = changed_cells,
    hash_matches = hash_matches,
    input_hashes = hashes,
    witness_id = witness_id,
    result_identity = result_identity,
    verification_fit_count = as.integer(fit_count),
    evaluation_identity = evaluation_identity,
    verification_elapsed_seconds = if (is.null(started)) NA_real_ else
      unname(proc.time()[["elapsed"]] - started)
  ), class = "ab_verification")
}

#' Replay and verify a reported annotation witness
#'
#' @param x A successful `ab_validation` object.
#' @param result An exact `ab_result` or practical `ab_search_result`.
#' @param witness Witness list; defaults to `result$witness`.
#' @param expected_input_hashes Expected immutable input/specification hashes.
#' @param fit_fun Optional explicit replay fit function; authoritative by default.
#' @param fit_identity Explicit versioned identity for a custom replay function.
#' @return An object of class `ab_verification`.
#' @export
verify_witness <- function(x, result, witness = result$witness,
                           expected_input_hashes = result$input_hashes,
                           fit_fun = NULL, fit_identity = NULL) {
  verification_started <- unname(proc.time()[["elapsed"]])
  errors <- character()
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    errors <- c(errors, "x must be a successful ab_validation object")
  }
  if (!inherits(result, "ab_result")) {
    errors <- c(errors, "result must be an ab_result object")
  }
  if (length(errors)) {
    return(.ab_verification_result(
      FALSE, errors, "verification_failed", data.frame(), data.frame(),
      data.frame(), character(), data.frame(), NULL, NA_character_, 0L,
      started = verification_started
    ))
  }
  data <- x$data
  result_identity <- .ab_result_identity(result)
  is_search <- inherits(result, "ab_search_result")
  if (is_search && result$status %in% c("exact_minimum", "exhaustively_stable")) {
    errors <- c(errors, "practical search results cannot carry exact or stable status")
  }
  if (is_search &&
      (!is.na(result$exact_minimum_cost) ||
       !is.na(result$certified_lower_bound))) {
    errors <- c(errors,
                "practical search results cannot carry an exact cost or certified lower bound")
  }

  observed_hashes <- tryCatch(.ab_recompute_input_hashes(x), error = function(e) e)
  if (inherits(observed_hashes, "error")) {
    errors <- c(errors, paste0("input hash recomputation failed: ",
                               conditionMessage(observed_hashes)))
    observed_hashes <- character()
  }
  hash_matches <- rbind(
    transform(.ab_compare_hashes(x$input_hashes, observed_hashes),
              source = "validation"),
    transform(.ab_compare_hashes(result$input_hashes, observed_hashes),
              source = "result"),
    transform(.ab_compare_hashes(expected_input_hashes, observed_hashes),
              source = "expected")
  )
  hash_matches <- hash_matches[, c(
    "source", "hash_name", "expected", "observed", "matches"
  ), drop = FALSE]
  if (!nrow(hash_matches) || any(!hash_matches$matches)) {
    errors <- c(errors, "immutable input/specification/constraint hashes do not match")
  }
  if (!identical(result$claim, data$claim) ||
      !identical(result$constraints, data$constraints) ||
      !identical(result$analysis_spec, data$analysis_spec)) {
    errors <- c(errors, "result analysis contract does not match validated inputs")
  }
  data$input_hashes <- observed_hashes

  authoritative <- .ab_resolve_fit(x)
  identity_type <- result$evaluation_identity$type
  if (!is.character(identity_type) || length(identity_type) != 1L ||
      is.na(identity_type) || !identity_type %in% c("authoritative", "custom")) {
    errors <- c(errors, "result evaluation identity type is missing or unsupported")
    replay_resolved <- authoritative
  } else if (identical(identity_type, "authoritative")) {
    replay_resolved <- authoritative
    if (!identical(result$evaluation_identity, authoritative$identity)) {
      errors <- c(errors, "result authoritative evaluation identity does not match x")
    }
    if (!is.null(fit_fun) || !is.null(fit_identity)) {
      errors <- c(errors,
                  "authoritative result verification cannot use a custom replay identity")
    }
  } else {
    replay_resolved <- authoritative
    if (is.null(fit_fun)) {
      errors <- c(errors,
                  "custom result verification requires the matching explicit replay fit_fun")
    } else if (!is.function(fit_fun)) {
      errors <- c(errors, "fit_fun must be a function")
    } else {
      requires_label <- !is.null(data$analysis_spec) ||
        !identical(result$evaluation_identity$label, "anonymous-custom-fit")
      if (requires_label &&
          (is.null(fit_identity) || !is.character(fit_identity) ||
           length(fit_identity) != 1L || is.na(fit_identity) ||
           !nzchar(fit_identity))) {
        errors <- c(errors,
                    "custom result verification requires its explicit fit_identity")
      }
      supplied <- tryCatch(
        .ab_resolve_fit(x, fit_fun, fit_identity), error = function(e) e
      )
      if (inherits(supplied, "error")) {
        errors <- c(errors, paste0("invalid custom replay identity: ",
                                   conditionMessage(supplied)))
      } else {
        if (!identical(supplied$identity, result$evaluation_identity)) {
          errors <- c(errors,
                      "verification fit_fun/fit_identity does not match result")
        }
        replay_resolved <- supplied
      }
    }
  }
  fit_count <- 0L
  replay_fit <- function(assignment) {
    fit_count <<- fit_count + 1L
    assignment_id <- .ab_assignment_id(data, assignment)
    value <- tryCatch(
      .ab_call_fit_fun(replay_resolved$fun, x, assignment), error = function(e) e
    )
    .ab_standardize_fit(value, data, assignment, assignment_id)
  }

  before <- tryCatch(build_pseudobulk(x, data$original_assignment, method = "full"),
                     error = function(e) e)
  before_table <- if (inherits(before, "error")) {
    errors <- c(errors, paste0("original pseudobulk replay failed: ",
                               conditionMessage(before)))
    data.frame()
  } else .ab_long_pseudobulk(data, before)
  original_replay <- replay_fit(data$original_assignment)
  reported_original <- result$original_fit
  if (!inherits(reported_original, "ab_fit") ||
      !identical(reported_original$status, original_replay$status)) {
    errors <- c(errors, "reported original fit status does not match fresh replay")
  } else if (identical(original_replay$status, "evaluable")) {
    if (!is.finite(reported_original$beta) ||
        abs(reported_original$beta - original_replay$beta) > data$claim$tolerance ||
        !is.finite(reported_original$margin) ||
        abs(reported_original$margin - original_replay$margin) >
          data$claim$tolerance ||
        !identical(reported_original$claim_holds, original_replay$claim_holds) ||
        !identical(reported_original$overturn, original_replay$overturn)) {
      errors <- c(errors, "reported original fit does not match fresh immutable replay")
    }
  }

  changed_cells <- .ab_empty_changed_cells(data)
  if (is.null(witness)) {
    if (!is.null(result$witness) || !is.null(result$witness_id) ||
        !is.na(result$best_verified_cost) ||
        !is.na(result$verified_upper_bound) ||
        !is.na(result$exact_minimum_cost) ||
        !is.na(result$certified_lower_bound)) {
      errors <- c(errors,
                  "witnessless result must not carry witness or breakpoint evidence")
    }
    if (identical(result$status, "exhaustively_stable") && !is_search) {
      domain <- tryCatch(.ab_candidate_domain(data), error = function(e) e)
      if (inherits(domain, "error")) {
        errors <- c(errors, paste0("assignment-domain replay failed: ",
                                   conditionMessage(domain)))
      } else {
        errors <- c(errors, .ab_verify_certified_status(
          x, result, domain, replay_fit
        ))
      }
      if (!identical(original_replay$status, "evaluable") ||
          !isTRUE(original_replay$claim_holds)) {
        errors <- c(errors, "original claim failed during stable-result replay")
      }
    } else if (identical(result$status, "original_claim_not_met")) {
      if (!identical(original_replay$status, "evaluable") ||
          isTRUE(original_replay$claim_holds)) {
        errors <- c(errors, "original_claim_not_met does not match fresh replay")
      }
    } else if (identical(result$status, "not_assessable")) {
      if (identical(original_replay$status, "evaluable")) {
        errors <- c(errors, "not_assessable does not match fresh replay")
      }
    } else if (identical(result$status, "no_counterexample_found")) {
      if (!identical(original_replay$status, "evaluable") ||
          !isTRUE(original_replay$claim_holds)) {
        errors <- c(errors,
                    "no_counterexample_found requires an evaluable original claim that holds")
      }
      if (any(result$candidates$overturn, na.rm = TRUE)) {
        errors <- c(errors, "no_counterexample_found contains a recorded overturn")
      }
      if (!is_search && isTRUE(result$search_complete)) {
        errors <- c(errors, "partial no-hit result is incorrectly marked complete")
      }
    } else if (identical(result$status, "exhaustive_search_incomplete") &&
               !is_search) {
      if (!all(c("feasible", "evaluated", "evaluable") %in%
               names(result$candidates)) ||
          !any(result$candidates$feasible & result$candidates$evaluated &
               !result$candidates$evaluable)) {
        errors <- c(errors, "incomplete-search accounting has no fit failure")
      }
    } else {
      errors <- c(errors, "a witness is required for this result status")
    }
    return(.ab_verification_result(
      !length(errors), errors, result$status, before_table, before_table,
      changed_cells, observed_hashes, hash_matches, NULL, result_identity,
      fit_count, original_replay$beta, evaluation_identity = replay_resolved$identity,
      started = verification_started
    ))
  }

    if (!is.list(witness)) {
    errors <- c(errors, "witness must be a named list")
    assignment <- NULL
  } else assignment <- witness$assignment
  assignment <- tryCatch(
    .ab_normalize_assignment(data, assignment),
    error = function(e) {
      errors <<- c(errors, paste0("witness assignment is invalid: ",
                                  conditionMessage(e)))
      NULL
    }
  )
  after_table <- data.frame()
  recomputed_cost <- NA_integer_
  recomputed_beta <- NA_real_
  witness_id <- NULL
  if (!is.null(assignment)) {
    if (!identical(original_replay$status, "evaluable") ||
        !isTRUE(original_replay$claim_holds)) {
      errors <- c(errors, "a witness requires an evaluable original claim that holds")
    }
    failures <- .ab_constraint_failures(data, assignment)
    if (!.ab_direction_ok(data, assignment,
                          if (is_search) result$direction else "all")) {
      failures <- c(failures, "direction restriction")
    }
    if (length(failures)) {
      errors <- c(errors, paste0("witness is infeasible: ",
                                 paste(unique(failures), collapse = ";")))
    }
    changed_cells <- .ab_changed_with_provenance(data, assignment)
    recomputed_cost <- nrow(changed_cells)
    witness_id <- .ab_assignment_id(data, assignment)
    if (is.null(witness$witness_id) || !identical(witness$witness_id, witness_id) ||
        (!is.null(witness$assignment_id) &&
         !identical(witness$assignment_id, witness_id))) {
      errors <- c(errors, "witness ID does not match its canonical assignment JSON hash")
    }
    cost_valid <- is.numeric(witness$cost) && length(witness$cost) == 1L &&
      is.finite(witness$cost) && witness$cost >= 0 &&
      witness$cost == floor(witness$cost) && witness$cost == recomputed_cost
    if (!cost_valid) errors <- c(errors, "witness cost does not match immutable replay")
    if (!is.null(witness$changed_cells) &&
        !.ab_compare_changed_cells(witness$changed_cells, changed_cells)) {
      errors <- c(errors, "witness changed-cell identity or provenance was altered")
    }
    if (!is.null(witness$input_hashes) &&
        any(!.ab_compare_hashes(witness$input_hashes, observed_hashes)$matches)) {
      errors <- c(errors, "witness immutable hashes do not match input replay")
    }
    after <- tryCatch(build_pseudobulk(x, assignment, method = "full"),
                      error = function(e) e)
    if (inherits(after, "error")) {
      errors <- c(errors, paste0("witness pseudobulk replay failed: ",
                                 conditionMessage(after)))
    } else {
      after_table <- .ab_long_pseudobulk(data, after)
      witness_replay <- replay_fit(assignment)
      if (!identical(witness_replay$status, "evaluable")) {
        errors <- c(errors, paste0("witness fit replay failed: ",
                                   witness_replay$error))
      } else {
        recomputed_beta <- witness_replay$beta
        if (!isTRUE(witness_replay$overturn)) {
          errors <- c(errors, "replayed witness does not overturn the registered claim")
        }
        if (!is.numeric(witness$beta) || length(witness$beta) != 1L ||
            !is.finite(witness$beta) ||
            abs(witness$beta - recomputed_beta) > data$claim$tolerance) {
          errors <- c(errors, "witness effect does not match immutable replay")
        }
        if (!is.numeric(witness$margin) || length(witness$margin) != 1L ||
            !is.finite(witness$margin) ||
            abs(witness$margin - witness_replay$margin) > data$claim$tolerance) {
          errors <- c(errors, "witness margin does not match immutable replay")
        }
      }
    }
    if (is.na(result$best_verified_cost) ||
        result$best_verified_cost != recomputed_cost) {
      errors <- c(errors, "result best_verified_cost does not match witness replay")
    }
    if (!identical(result$witness_id, witness_id)) {
      errors <- c(errors, "result witness_id does not match witness replay")
    }
    if (!result$status %in% c(
      "exact_minimum", "verified_counterexample", "exhaustive_search_incomplete"
    )) {
      errors <- c(errors, "result status is incompatible with a verified witness")
    }
    if (!is_search && result$status %in% c("exact_minimum", "exhaustively_stable")) {
      domain <- tryCatch(.ab_candidate_domain(data), error = function(e) e)
      if (inherits(domain, "error")) {
        errors <- c(errors, paste0("assignment-domain replay failed: ",
                                   conditionMessage(domain)))
      } else {
        errors <- c(errors, .ab_verify_certified_status(
          x, result, domain, replay_fit, witness_cost = recomputed_cost
        ))
      }
    }
  }
  .ab_verification_result(
    !length(errors), errors, result$status, before_table, after_table,
    changed_cells, observed_hashes, hash_matches, witness_id, result_identity,
    fit_count, recomputed_beta, recomputed_cost, replay_resolved$identity,
    started = verification_started
  )
}
