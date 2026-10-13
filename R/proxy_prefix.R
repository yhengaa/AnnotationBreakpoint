.ab_proxy_prefix_sizes <- function(prefix_sizes = NULL) {
  if (is.null(prefix_sizes)) {
    return(c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L))
  }
  if (!is.numeric(prefix_sizes) || !length(prefix_sizes) ||
      anyNA(prefix_sizes) || any(!is.finite(prefix_sizes)) ||
      any(prefix_sizes < 1) || any(prefix_sizes != floor(prefix_sizes)) ||
      anyDuplicated(prefix_sizes) || is.unsorted(prefix_sizes, strictly = TRUE)) {
    stop(
      "prefix_sizes must be a strictly increasing list of unique positive integers",
      call. = FALSE
    )
  }
  as.integer(prefix_sizes)
}

.ab_proxy_prefix_plan <- function(data, prefix_sizes = NULL) {
  cache <- .ab_proxy_cell_cache(data)
  ranked <- .ab_proxy_ranked_source(
    data,
    data$original_assignment,
    direction = "all",
    proxy_cache = cache,
    target_pseudocount = 0.5,
    library_pseudocount = 1,
    scale = 1e6
  )
  if (!length(ranked)) {
    ranking <- data.frame(
      proxy_rank = integer(),
      cell_id = character(),
      from_label = character(),
      to_label = character(),
      proxy_score = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    ranking <- data.frame(
      proxy_rank = seq_along(ranked),
      cell_id = vapply(ranked, `[[`, character(1L), "cell_id"),
      from_label = vapply(ranked, `[[`, character(1L), "from_label"),
      to_label = vapply(ranked, `[[`, character(1L), "to_label"),
      proxy_score = vapply(ranked, `[[`, numeric(1L), "proxy_score"),
      stringsAsFactors = FALSE
    )
  }
  if (anyDuplicated(ranking$cell_id)) {
    stop(
      "proxy_prefix requires exactly one permitted alternative per ambiguous cell",
      call. = FALSE
    )
  }
  registered <- .ab_proxy_prefix_sizes(prefix_sizes)
  available <- registered
  available <- available[available <= nrow(ranking)]
  assignments <- if (length(available)) {
    stats::setNames(lapply(available, function(k) {
      assignment <- data$original_assignment
      selected <- ranking[seq_len(k), , drop = FALSE]
      assignment[selected$cell_id] <- selected$to_label
      .ab_normalize_assignment(data, assignment)
    }), paste0("k", available))
  } else {
    list()
  }
  list(
    route = "observed_centered_proxy_ranked_cumulative_prefixes",
    ranking = ranking,
    registered_prefix_sizes = registered,
    available_prefix_sizes = available,
    assignments = assignments,
    ranking_frozen_before_fits = TRUE,
    diagnostics_included = FALSE
  )
}

.ab_proxy_prefix_candidate <- function(data, assignment, rank, k,
                                       proxy_score, failures = character(),
                                       fit = NULL, evaluated = FALSE,
                                       reason = NA_character_) {
  changes <- .ab_change_rows(data, assignment)
  evaluable <- inherits(fit, "ab_fit") && identical(fit$status, "evaluable")
  warnings <- if (!is.null(fit) && length(fit$warnings)) {
    paste(unique(as.character(fit$warnings)), collapse = " | ")
  } else NA_character_
  failure_reason <- if (!is.null(fit) && length(fit$failure_reason) == 1L) {
    as.character(fit$failure_reason)
  } else reason
  evaluation_error <- if (!is.null(fit) && length(fit$error) == 1L) {
    as.character(fit$error)
  } else NA_character_
  row <- data.frame(
    rank = as.integer(rank),
    assignment_id = .ab_assignment_id(data, assignment),
    changes = .ab_changes_text(changes),
    cost = as.integer(nrow(changes)),
    feasible = !length(failures),
    constraint_failures = paste(failures, collapse = ";"),
    evaluated = isTRUE(evaluated),
    evaluable = evaluable,
    beta = if (evaluable) as.numeric(fit$beta) else NA_real_,
    margin = if (evaluable) as.numeric(fit$margin) else NA_real_,
    claim_holds = if (evaluable) isTRUE(fit$claim_holds) else NA,
    overturn = evaluable && isTRUE(fit$overturn),
    fit_elapsed_seconds = if (!is.null(fit) &&
                              length(fit$fit_elapsed_seconds) == 1L) {
      as.numeric(fit$fit_elapsed_seconds)
    } else NA_real_,
    warnings = warnings,
    failure_reason = failure_reason,
    evaluation_error = evaluation_error,
    prefix_size = as.integer(k),
    terminal_proxy_score = as.numeric(proxy_score),
    stringsAsFactors = FALSE
  )
  row$assignment <- I(list(assignment))
  row
}

.ab_proxy_prefix_result <- function(data, hashes, resolved_identity, original_fit,
                                    candidates, witness_index = NULL,
                                    method_seed, max_fits, max_proposals,
                                    max_seconds, attempted_fits, elapsed_seconds,
                                    termination_reason, ranking,
                                    registered_sizes, available_sizes,
                                    search_censored = FALSE,
                                    verification_fits = 0L) {
  original_evaluable <- inherits(original_fit, "ab_fit") &&
    identical(original_fit$status, "evaluable")
  has_witness <- !is.null(witness_index) && length(witness_index) == 1L &&
    !is.na(witness_index)
  status <- if (!original_evaluable) {
    "not_assessable"
  } else if (!isTRUE(original_fit$claim_holds)) {
    "original_claim_not_met"
  } else if (has_witness) {
    "verified_counterexample"
  } else {
    "no_counterexample_found"
  }
  result <- .ab_search_skeleton(status, termination_reason, hashes)
  result$original_fit <- original_fit
  result$assignments <- candidates
  result$candidates <- candidates
  result$effects_by_cost <- .ab_effects_by_cost(candidates)
  result$total_assignment_count <- nrow(candidates)
  result$feasible_assignment_count <- sum(candidates$feasible, na.rm = TRUE)
  result$infeasible_assignment_count <- sum(!candidates$feasible, na.rm = TRUE)
  result$evaluable_assignment_count <- sum(candidates$evaluable, na.rm = TRUE)
  result$unevaluable_assignment_count <- sum(
    candidates$evaluated & !candidates$evaluable, na.rm = TRUE
  )
  result$unevaluated_assignment_count <- sum(!candidates$evaluated, na.rm = TRUE)
  result$overturn_assignment_count <- sum(candidates$overturn, na.rm = TRUE)
  result$method <- "proxy_prefix"
  result$seed <- as.integer(method_seed)
  result$direction <- "all"
  result$max_fits <- if (is.infinite(max_fits)) Inf else as.integer(max_fits)
  result$max_proposals <- max_proposals
  result$max_seconds <- max_seconds
  result$proposal_count <- as.integer(sum(
    is.na(candidates$failure_reason) |
      !grepl("^not_proposed:", candidates$failure_reason)
  ))
  result$unique_assignment_count <- nrow(candidates)
  result$unique_proposal_count <- result$proposal_count
  result$infeasible_proposal_count <- sum(!candidates$feasible, na.rm = TRUE)
  result$attempted_fit_count <- as.integer(attempted_fits)
  result$new_fit_count <- as.integer(attempted_fits)
  result$cumulative_fit_count <- as.integer(attempted_fits)
  result$successful_fit_count <- as.integer(original_evaluable) +
    sum(candidates$evaluable, na.rm = TRUE)
  result$failed_fit_count <- as.integer(!original_evaluable) +
    sum(candidates$evaluated & !candidates$evaluable, na.rm = TRUE)
  result$elapsed_seconds <- as.numeric(elapsed_seconds)
  result$search_censored <- isTRUE(search_censored)
  result$search_complete <- FALSE
  result$proxy_prefix <- list(
    route = "observed_centered_proxy_ranked_cumulative_prefixes",
    registered_prefix_sizes = as.integer(registered_sizes),
    available_prefix_sizes = as.integer(available_sizes),
    ranking = ranking,
    ranking_frozen_before_fits = TRUE,
    diagnostics_included = FALSE,
    no_hit_semantics = "unresolved_registered_prefix_no_hit"
  )
  result$fit_accounting <- data.frame(
    proposals = result$proposal_count,
    unique_assignments = result$unique_assignment_count,
    infeasible_proposals = result$infeasible_proposal_count,
    attempted_fits = result$attempted_fit_count,
    successful_fits = result$successful_fit_count,
    failed_fits = result$failed_fit_count,
    original_fits = as.integer(attempted_fits >= 1L),
    prefix_fits = sum(candidates$evaluated, na.rm = TRUE),
    verification_fits = as.integer(verification_fits),
    stringsAsFactors = FALSE
  )
  result$search_path <- data.frame()
  if (has_witness) {
    assignment <- candidates$assignment[[witness_index]]
    changes <- .ab_changed_with_provenance(data, assignment)
    witness_id <- candidates$assignment_id[[witness_index]]
    result$best_verified_cost <- as.integer(nrow(changes))
    result$verified_upper_bound <- result$best_verified_cost
    result$witness_id <- witness_id
    result$witness <- list(
      witness_id = witness_id,
      assignment_id = witness_id,
      assignment = assignment,
      changes = .ab_changes_text(changes),
      changed_cells = changes,
      cost = result$best_verified_cost,
      beta = candidates$beta[[witness_index]],
      margin = candidates$margin[[witness_index]],
      overturn = TRUE,
      input_hashes = hashes,
      analysis_id = data$claim$estimator_id
    )
    result$minimum_cost_overturns <- candidates[witness_index, , drop = FALSE]
    result$search_path <- data.frame(
      step = 1L,
      assignment_id = witness_id,
      cost = result$best_verified_cost,
      beta = candidates$beta[[witness_index]],
      margin = candidates$margin[[witness_index]],
      stringsAsFactors = FALSE
    )
  }
  class(result) <- c("ab_search_result", "ab_result")
  .ab_attach_contract(result, data, resolved_identity)
}

.ab_search_proxy_prefix <- function(x, data, hashes, resolved, max_fits,
                                    seed, max_proposals, max_seconds,
                                    prefix_sizes = NULL) {
  started <- unname(proc.time()[["elapsed"]])
  plan <- tryCatch(
    .ab_proxy_prefix_plan(data, prefix_sizes = prefix_sizes),
    error = function(e) e
  )
  if (inherits(plan, "error")) {
    result <- .ab_search_skeleton("invalid_input", "invalid_proxy_prefix_domain", hashes)
    result$errors <- conditionMessage(plan)
    result$method <- "proxy_prefix"
    result$seed <- as.integer(seed)
    result$direction <- "all"
    return(.ab_attach_contract(result, data, resolved$identity))
  }
  elapsed <- function() unname(proc.time()[["elapsed"]] - started)
  fit_one <- function(assignment) {
    assignment_id <- .ab_assignment_id(data, assignment)
    value <- tryCatch(
      .ab_call_fit_fun(resolved$fun, x, assignment), error = function(e) e
    )
    .ab_standardize_fit(value, data, assignment, assignment_id)
  }
  attempted <- 0L
  original_fit <- fit_one(data$original_assignment)
  attempted <- attempted + 1L
  prefix_rows <- vector("list", length(plan$available_prefix_sizes))
  hit_indices <- integer()
  stopped <- FALSE
  termination <- "registered_prefixes_exhausted_no_verified_counterexample"
  censored <- FALSE

  if (identical(original_fit$status, "evaluable") &&
      isTRUE(original_fit$claim_holds)) {
    for (index in seq_along(plan$available_prefix_sizes)) {
      k <- plan$available_prefix_sizes[[index]]
      assignment <- plan$assignments[[paste0("k", k)]]
      proxy_score <- plan$ranking$proxy_score[[k]]
      reason <- NA_character_
      failures <- .ab_constraint_failures(data, assignment)
      fit <- NULL
      evaluated <- FALSE
      if (stopped) {
        reason <- paste0("not_proposed:", termination)
      } else if (index > max_proposals) {
        termination <- "max_proposals_reached"
        censored <- TRUE
        stopped <- TRUE
        reason <- "not_proposed:max_proposals_reached"
      } else if (elapsed() >= max_seconds) {
        termination <- "max_seconds_reached"
        censored <- TRUE
        stopped <- TRUE
        reason <- "not_proposed:max_seconds_reached"
      } else {
        if (!length(failures)) {
          if (attempted >= max_fits) {
            termination <- "max_fits_reached"
            censored <- TRUE
            stopped <- TRUE
            reason <- "not_evaluated:max_fits_reached"
          } else {
            fit <- fit_one(assignment)
            attempted <- attempted + 1L
            evaluated <- TRUE
            if (identical(fit$status, "evaluable") && isTRUE(fit$overturn)) {
              hit_indices <- c(hit_indices, index)
            }
          }
        }
      }
      prefix_rows[[index]] <- .ab_proxy_prefix_candidate(
        data, assignment, index, k, proxy_score, failures, fit, evaluated, reason
      )
    }
  } else {
    termination <- if (!identical(original_fit$status, "evaluable")) {
      "original_fit_not_assessable"
    } else {
      "original_claim_not_met"
    }
    for (index in seq_along(plan$available_prefix_sizes)) {
      k <- plan$available_prefix_sizes[[index]]
      assignment <- plan$assignments[[paste0("k", k)]]
      prefix_rows[[index]] <- .ab_proxy_prefix_candidate(
        data, assignment, index, k, plan$ranking$proxy_score[[k]],
        failures = .ab_constraint_failures(data, assignment),
        reason = paste0("not_proposed:", termination)
      )
    }
  }
  candidates <- if (length(prefix_rows)) {
    do.call(rbind, prefix_rows)
  } else {
    empty <- .ab_empty_candidate_table()
    empty$prefix_size <- integer()
    empty$terminal_proxy_score <- numeric()
    empty$assignment <- I(list())
    empty
  }

  build_result <- function(witness_index = NULL, reason = termination,
                           verification_fits = 0L) {
    .ab_proxy_prefix_result(
      data, hashes, resolved$identity, original_fit, candidates,
      witness_index = witness_index, method_seed = seed,
      max_fits = max_fits, max_proposals = max_proposals,
      max_seconds = max_seconds, attempted_fits = attempted,
      elapsed_seconds = elapsed(), termination_reason = reason,
      ranking = plan$ranking,
      registered_sizes = plan$registered_prefix_sizes,
      available_sizes = plan$available_prefix_sizes,
      search_censored = censored,
      verification_fits = verification_fits
    )
  }

  search_elapsed <- elapsed()
  if (length(hit_indices)) {
    replay_errors <- character()
    replay_attempts <- list()
    replay_fit_total <- 0L
    for (index in hit_indices) {
      result <- build_result(
        index, "first_replay_verified_registered_prefix",
        verification_fits = replay_fit_total + 2L
      )
      replay_started <- unname(proc.time()[["elapsed"]])
      verification <- verify_witness(
        x, result,
        fit_fun = if (identical(resolved$identity$type, "custom")) resolved$fun else NULL,
        fit_identity = if (identical(resolved$identity$type, "custom")) {
          resolved$identity$label
        } else NULL
      )
      replay_elapsed <- unname(proc.time()[["elapsed"]] - replay_started)
      replay_fit_total <- replay_fit_total + verification$verification_fit_count
      replay_attempts[[length(replay_attempts) + 1L]] <- list(
        assignment_id = candidates$assignment_id[[index]],
        ok = isTRUE(verification$ok),
        errors = verification$errors,
        verification_fit_count = verification$verification_fit_count,
        verification_elapsed_seconds = replay_elapsed
      )
      if (isTRUE(verification$ok)) {
        result$independent_verification <- verification
        result$verification_fit_count <- as.integer(replay_fit_total)
        result$verification_embedded <- TRUE
        result$replay_attempts <- replay_attempts
        result$search_elapsed_seconds <- search_elapsed
        result$replay_elapsed_seconds <- sum(vapply(
          replay_attempts, `[[`, numeric(1L), "verification_elapsed_seconds"
        ))
        result$all_in_elapsed_seconds <- elapsed()
        return(result)
      }
      replay_errors <- c(
        replay_errors,
        paste0(candidates$assignment_id[[index]], ": ",
               paste(verification$errors, collapse = "; "))
      )
    }
    stop(
      "fitted prefix counterexamples failed independent replay: ",
      paste(replay_errors, collapse = " | "),
      call. = FALSE
    )
  }

  result <- build_result(verification_fits = 1L)
  replay_started <- unname(proc.time()[["elapsed"]])
  verification <- verify_witness(
    x, result,
    fit_fun = if (identical(resolved$identity$type, "custom")) resolved$fun else NULL,
    fit_identity = if (identical(resolved$identity$type, "custom")) {
      resolved$identity$label
    } else NULL
  )
  replay_elapsed <- unname(proc.time()[["elapsed"]] - replay_started)
  if (!isTRUE(verification$ok)) {
    stop(
      "proxy_prefix result failed independent replay: ",
      paste(verification$errors, collapse = "; "),
      call. = FALSE
    )
  }
  result$independent_verification <- verification
  result$verification_fit_count <- verification$verification_fit_count
  result$verification_embedded <- TRUE
  result$replay_attempts <- list(list(
    assignment_id = NA_character_,
    ok = TRUE,
    errors = character(),
    verification_fit_count = verification$verification_fit_count,
    verification_elapsed_seconds = replay_elapsed
  ))
  result$search_elapsed_seconds <- search_elapsed
  result$replay_elapsed_seconds <- replay_elapsed
  result$all_in_elapsed_seconds <- elapsed()
  result
}
