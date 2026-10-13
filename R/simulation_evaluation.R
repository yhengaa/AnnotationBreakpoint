.ab_sim_read_table <- function(path, numeric_columns = character(),
                               logical_columns = character()) {
  out <- utils::read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = "NA"
  )
  for (name in intersect(numeric_columns, names(out))) {
    out[[name]] <- suppressWarnings(as.numeric(out[[name]]))
  }
  for (name in intersect(logical_columns, names(out))) {
    out[[name]] <- as.logical(out[[name]])
  }
  out
}

.ab_sim_elapsed <- function(fun) {
  started <- unname(proc.time()[["elapsed"]])
  value <- fun()
  list(value = value,
       elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
}

.ab_sim_fit_count <- function(result) {
  if (inherits(result, "ab_search_result")) {
    return(as.integer(result$attempted_fit_count))
  }
  if (is.data.frame(result$candidates) && "evaluated" %in% names(result$candidates)) {
    return(as.integer(sum(result$candidates$evaluated, na.rm = TRUE)))
  }
  as.integer(!is.null(result$original_fit))
}

.ab_sim_failed_fit_count <- function(result) {
  if (inherits(result, "ab_search_result")) {
    return(as.integer(result$failed_fit_count))
  }
  if (is.data.frame(result$candidates) &&
      all(c("evaluated", "evaluable") %in% names(result$candidates))) {
    return(as.integer(sum(result$candidates$evaluated &
                            !result$candidates$evaluable, na.rm = TRUE)))
  }
  as.integer(!is.null(result$original_fit) &&
               !identical(result$original_fit$status, "evaluable"))
}

.ab_sim_run_search <- function(validation, method, common, seed, direction) {
  searched <- .ab_sim_elapsed(function() search_breakpoint(
    validation,
    method = method,
    max_fits = as.integer(common$practical_max_fits),
    seed = as.integer(seed),
    max_proposals = as.numeric(common$practical_max_proposals),
    max_seconds = as.numeric(common$practical_max_seconds),
    direction = direction
  ))
  checked <- .ab_sim_elapsed(function() verify_witness(validation, searched$value))
  list(
    result = searched$value,
    verification = checked$value,
    timing = list(
      search_seconds = searched$elapsed_seconds,
      verification_seconds = checked$elapsed_seconds,
      total_seconds = searched$elapsed_seconds + checked$elapsed_seconds
    )
  )
}

.ab_sim_relaxed_utility_validation <- function(input, probabilities,
                                               analysis_spec) {
  uncertain <- names(which(table(probabilities$cell_id) > 1L))
  cells <- input$cells
  specimens <- input$specimens
  uncertain_rows <- cells$cell_id %in% uncertain
  sample_counts <- table(factor(cells$sample_id[uncertain_rows],
                                levels = specimens$sample_id))
  donor_counts <- table(factor(cells$donor_id[uncertain_rows],
                               levels = specimens$donor_id))
  raw_constraints <- as.list(input$constraints[1L, , drop = FALSE])
  raw_constraints <- lapply(raw_constraints, function(value) value[[1L]])
  raw_constraints$max_changes_global <- as.integer(length(uncertain))
  raw_constraints$max_changes_per_sample <- setNames(
    as.numeric(sample_counts), specimens$sample_id
  )
  raw_constraints$max_changes_per_donor <- setNames(
    as.numeric(donor_counts), specimens$donor_id
  )
  validate_inputs(
    input$counts, input$cells, input$specimens, input$allowed_labels,
    input$claim, raw_constraints, analysis_spec
  )
}

.ab_sim_map_assignment <- function(probabilities, cell_order) {
  rows <- split(probabilities, factor(probabilities$cell_id, levels = cell_order))
  assignment <- vapply(rows, function(block) {
    block <- block[order(-block$probability, block$label, method = "radix"),
                   , drop = FALSE]
    block$label[[1L]]
  }, character(1L))
  setNames(assignment, cell_order)
}

.ab_sim_draw_assignments <- function(probabilities, cell_order, draws, seed) {
  if (length(draws) != 1L || is.na(draws) || draws < 0 || draws != floor(draws)) {
    stop("probability_draws must be a nonnegative integer", call. = FALSE)
  }
  rows <- split(probabilities, factor(probabilities$cell_id, levels = cell_order))
  .ab_sim_with_local_seed(seed, function() {
    lapply(seq_len(as.integer(draws)), function(draw_index) {
      assignment <- vapply(rows, function(block) {
        sample(block$label, size = 1L, prob = block$probability)
      }, character(1L))
      setNames(assignment, cell_order)
    })
  })
}

.ab_sim_utility_row <- function(family_id, realization_id, seed, method,
                                fit = NULL, validation = NULL,
                                assignment = NULL, original_assignment = NULL,
                                true_fit = NULL, draw_index = NA_integer_,
                                removed_cell_count = 0L,
                                elapsed_seconds = NA_real_,
                                forced_status = NULL,
                                forced_reason = NA_character_,
                                feasible_override = NULL) {
  validation_ok <- inherits(validation, "ab_validation") && isTRUE(validation$ok)
  feasible <- if (is.null(feasible_override)) validation_ok else {
    isTRUE(feasible_override)
  }
  fit_attempted <- !is.null(fit)
  evaluable <- fit_attempted && identical(fit$status, "evaluable")
  status <- if (!is.null(forced_status)) forced_status else if (!validation_ok) {
    "invalid_input"
  } else if (!fit_attempted) {
    "infeasible"
  } else fit$status
  beta <- if (evaluable) as.numeric(fit$beta) else NA_real_
  true_beta <- if (!is.null(true_fit) && identical(true_fit$status, "evaluable")) {
    as.numeric(true_fit$beta)
  } else NA_real_
  agreement <- if (evaluable && !is.null(true_fit) &&
                   identical(true_fit$status, "evaluable")) {
    identical(isTRUE(fit$claim_holds), isTRUE(true_fit$claim_holds))
  } else NA
  failure_reason <- if (!is.na(forced_reason)) {
    forced_reason
  } else if (!validation_ok) {
    paste(validation$errors, collapse = "; ")
  } else if (!is.null(fit) && !is.null(fit$failure_reason)) {
    paste(fit$failure_reason, collapse = "; ")
  } else NA_character_
  hamming <- if (!is.null(assignment) && !is.null(original_assignment)) {
    as.integer(sum(assignment[names(original_assignment)] != original_assignment))
  } else NA_integer_
  data.frame(
    family_id = family_id,
    realization_id = realization_id,
    seed = as.integer(seed),
    method = method,
    draw_index = as.integer(draw_index),
    status = status,
    feasible = feasible,
    fit_attempted = fit_attempted,
    fit_status = if (fit_attempted) fit$status else NA_character_,
    beta = beta,
    true_label_beta = true_beta,
    beta_minus_true = if (is.finite(beta) && is.finite(true_beta)) {
      beta - true_beta
    } else NA_real_,
    absolute_beta_error = if (is.finite(beta) && is.finite(true_beta)) {
      abs(beta - true_beta)
    } else NA_real_,
    claim_agreement = agreement,
    hamming_cost = hamming,
    removed_cell_count = as.integer(removed_cell_count),
    fit_elapsed_seconds = if (!is.null(fit) &&
                              !is.null(fit$fit_elapsed_seconds)) {
      as.numeric(fit$fit_elapsed_seconds)
    } else NA_real_,
    total_elapsed_seconds = as.numeric(elapsed_seconds),
    failure_reason = failure_reason,
    warnings = if (!is.null(fit) && length(fit$warnings)) {
      paste(unique(fit$warnings), collapse = " | ")
    } else NA_character_,
    stringsAsFactors = FALSE
  )
}

.ab_sim_fit_utility_assignment <- function(validation, assignment) {
  if (!inherits(validation, "ab_validation") || !isTRUE(validation$ok)) {
    return(list(fit = NULL, elapsed_seconds = 0, feasible = FALSE,
                failures = if (inherits(validation, "ab_validation")) {
                  validation$errors
                } else "invalid utility validation"))
  }
  normalized <- tryCatch(.ab_normalize_assignment(validation$data, assignment),
                         error = function(error) error)
  if (inherits(normalized, "error")) {
    return(list(fit = NULL, elapsed_seconds = 0, feasible = FALSE,
                failures = conditionMessage(normalized)))
  }
  failures <- .ab_constraint_failures(validation$data, normalized)
  if (length(failures)) {
    return(list(fit = NULL, elapsed_seconds = 0, feasible = FALSE,
                failures = failures))
  }
  timed <- .ab_sim_elapsed(function() fit_edger_claim(validation, normalized))
  list(fit = timed$value, elapsed_seconds = timed$elapsed_seconds,
       feasible = TRUE, failures = character())
}

.ab_sim_read_case_input <- function(directory, analysis_spec) {
  allowed_files <- c(
    "allowed_labels.tsv", "analysis_spec.yaml", "cells.tsv", "claim.tsv",
    "constraints.tsv", "counts.tsv", "specimens.tsv"
  )
  observed <- sort(list.files(directory, all.files = FALSE, no.. = TRUE),
                   method = "radix")
  if (!identical(observed, allowed_files)) {
    stop("search_input directory is not the strict six-table plus analysis-spec allowlist",
         call. = FALSE)
  }
  disk_spec <- yaml::read_yaml(file.path(directory, "analysis_spec.yaml"))
  disk_normalized <- .ab_normalize_analysis_spec(disk_spec)
  supplied_normalized <- .ab_normalize_analysis_spec(analysis_spec)
  if (!isTRUE(disk_normalized$ok) || !isTRUE(supplied_normalized$ok) ||
      !identical(.ab_hash(disk_normalized$spec),
                 .ab_hash(supplied_normalized$spec))) {
    stop("case analysis specification does not match evaluator analysis_spec",
         call. = FALSE)
  }
  read_input_directory(directory, supplied_normalized$spec)
}

.ab_evaluate_simulation_case <- function(case_dir, benchmark_spec,
                                         analysis_spec) {
  if (!is.character(case_dir) || length(case_dir) != 1L ||
      is.na(case_dir) || !dir.exists(case_dir)) {
    stop("case_dir must be one existing simulation-case directory", call. = FALSE)
  }
  if (!is.list(benchmark_spec) || is.null(benchmark_spec$common)) {
    stop("benchmark_spec must contain the registered common mapping", call. = FALSE)
  }
  root <- normalizePath(case_dir, mustWork = TRUE)
  search_dir <- file.path(root, "search_input")
  method_dir <- file.path(root, "method_input")
  truth_dir <- file.path(root, "evaluation_truth")
  input <- .ab_sim_read_case_input(search_dir, analysis_spec)
  validation <- validate_inputs(
    input$counts, input$cells, input$specimens, input$allowed_labels,
    input$claim, input$constraints, input$analysis_spec
  )
  if (!isTRUE(validation$ok)) {
    stop("written simulation search input is invalid: ",
         paste(validation$errors, collapse = "; "), call. = FALSE)
  }
  probabilities <- .ab_sim_read_table(
    file.path(method_dir, "annotation_probabilities.tsv"),
    numeric_columns = "probability"
  )
  required_probability <- c(
    "cell_id", "label", "probability", "probability_source", "probability_version"
  )
  if (!identical(names(probabilities), required_probability) ||
      any(!is.finite(probabilities$probability)) ||
      any(probabilities$probability < 0 | probabilities$probability > 1) ||
      !setequal(probabilities$cell_id, input$cells$cell_id)) {
    stop("method annotation probabilities are malformed", call. = FALSE)
  }
  probability_pairs <- paste(probabilities$cell_id, probabilities$label, sep = "\r")
  allowed_pairs <- paste(input$allowed_labels$cell_id,
                         input$allowed_labels$allowed_label, sep = "\r")
  totals <- tapply(probabilities$probability, probabilities$cell_id, sum)
  if (anyDuplicated(probability_pairs) || !setequal(probability_pairs, allowed_pairs) ||
      any(abs(totals - 1) > 1e-12)) {
    stop("method probabilities do not exactly match normalized allowed-label support",
         call. = FALSE)
  }
  metadata <- jsonlite::read_json(
    file.path(method_dir, "case_metadata.json"), simplifyVector = TRUE
  )
  family_id <- as.character(metadata$family_id)
  seed <- as.integer(metadata$seed)
  realization_id <- paste0(family_id, "__seed_", seed)
  common <- benchmark_spec$common
  methods <- as.character(unlist(common$practical_methods))
  if (!identical(methods, c("random", "greedy", "pair"))) {
    stop("registered practical methods must be random, greedy, and pair",
         call. = FALSE)
  }
  search_seed <- as.integer(seed + as.integer(common$search_seed_offset))

  exact_cap <- suppressWarnings(as.numeric(common$exact_max_evaluations))
  if (length(exact_cap) != 1L || is.na(exact_cap)) exact_cap <- Inf
  exact_timed <- .ab_sim_elapsed(function() enumerate_breakpoint(
    validation, max_evaluations = exact_cap
  ))
  exact_checked <- .ab_sim_elapsed(function() verify_witness(
    validation, exact_timed$value
  ))
  exact_timing <- list(
    search_seconds = exact_timed$elapsed_seconds,
    verification_seconds = exact_checked$elapsed_seconds,
    total_seconds = exact_timed$elapsed_seconds + exact_checked$elapsed_seconds
  )
  practical <- setNames(lapply(methods, function(method) {
    .ab_sim_run_search(validation, method, common, search_seed, "all")
  }), methods)
  outgoing <- if (isTRUE(common$outgoing_only_comparison)) {
    setNames(lapply(methods, function(method) {
      .ab_sim_run_search(validation, method, common, search_seed,
                         "outgoing_only")
    }), methods)
  } else list()

  # Truth is opened only after all optimizer runs have completed. The true-label
  # input has singleton truth support and therefore cannot share a search cache.
  generating_truth <- jsonlite::read_json(
    file.path(truth_dir, "generating_truth.json"), simplifyVector = TRUE
  )
  generation_record <- jsonlite::read_json(
    file.path(truth_dir, "generation_record.json"), simplifyVector = TRUE
  )
  generating_effects <- .ab_sim_read_table(
    file.path(truth_dir, "generating_effects.tsv"),
    numeric_columns = c("target_control_mean", "target_log2_condition_effect")
  )
  cell_truth <- .ab_sim_read_table(
    file.path(truth_dir, "cell_truth.tsv"),
    numeric_columns = "true_state",
    logical_columns = c("registered_uncertain", "locked_wrong", "is_mislabeled",
                        "correct_alternative_available")
  )
  true_input <- .ab_sim_read_case_input(file.path(truth_dir, "true_label_input"),
                                        analysis_spec)
  true_validation <- validate_inputs(
    true_input$counts, true_input$cells, true_input$specimens,
    true_input$allowed_labels, true_input$claim, true_input$constraints,
    true_input$analysis_spec
  )
  true_timed <- if (isTRUE(true_validation$ok)) {
    .ab_sim_elapsed(function() fit_edger_claim(true_validation))
  } else list(value = NULL, elapsed_seconds = 0)
  true_fit <- true_timed$value

  original_timed <- .ab_sim_elapsed(function() fit_edger_claim(validation))
  original_assignment <- validation$data$original_assignment
  utility_rows <- list(.ab_sim_utility_row(
    family_id, realization_id, seed, "original_observed_hard_labels",
    fit = original_timed$value, validation = validation,
    assignment = original_assignment, original_assignment = original_assignment,
    true_fit = true_fit, elapsed_seconds = original_timed$elapsed_seconds
  ))

  uncertain <- names(which(table(probabilities$cell_id) > 1L))
  retained <- !input$cells$cell_id %in% uncertain
  removal_input <- input
  removal_input$counts <- input$counts[, input$cells$cell_id[retained], drop = FALSE]
  removal_input$cells <- input$cells[retained, , drop = FALSE]
  removal_input$allowed_labels <- input$allowed_labels[
    input$allowed_labels$cell_id %in% removal_input$cells$cell_id, , drop = FALSE
  ]
  removal_validation <- validate_inputs(
    removal_input$counts, removal_input$cells, removal_input$specimens,
    removal_input$allowed_labels, removal_input$claim, removal_input$constraints,
    analysis_spec
  )
  removal_timed <- if (isTRUE(removal_validation$ok)) {
    .ab_sim_elapsed(function() fit_edger_claim(removal_validation))
  } else list(value = NULL, elapsed_seconds = 0)
  utility_rows[[length(utility_rows) + 1L]] <- .ab_sim_utility_row(
    family_id, realization_id, seed, "remove_all_registered_uncertain_cells",
    fit = removal_timed$value, validation = removal_validation,
    true_fit = true_fit, removed_cell_count = length(uncertain),
    elapsed_seconds = removal_timed$elapsed_seconds
  )

  utility_validation <- .ab_sim_relaxed_utility_validation(
    input, probabilities, analysis_spec
  )
  map_assignment <- .ab_sim_map_assignment(
    probabilities, validation$data$cells$cell_id
  )
  map_outcome <- .ab_sim_fit_utility_assignment(utility_validation, map_assignment)
  utility_rows[[length(utility_rows) + 1L]] <- .ab_sim_utility_row(
    family_id, realization_id, seed, "stipulated_probability_MAP_labels",
    fit = map_outcome$fit, validation = utility_validation,
    assignment = map_assignment, original_assignment = original_assignment,
    true_fit = true_fit, elapsed_seconds = map_outcome$elapsed_seconds,
    forced_status = if (map_outcome$feasible) NULL else "infeasible",
    forced_reason = if (length(map_outcome$failures)) {
      paste(map_outcome$failures, collapse = ";")
    } else NA_character_,
    feasible_override = map_outcome$feasible
  )

  draw_count <- as.integer(common$probability_draws)
  draw_seed <- as.integer(seed + as.integer(common$probability_draw_seed_offset))
  assignments <- .ab_sim_draw_assignments(
    probabilities, validation$data$cells$cell_id, draw_count, draw_seed
  )
  draw_rows <- lapply(seq_along(assignments), function(index) {
    outcome <- .ab_sim_fit_utility_assignment(utility_validation,
                                              assignments[[index]])
    .ab_sim_utility_row(
      family_id, realization_id, seed,
      "stipulated_probability_integer_assignment_draws",
      fit = outcome$fit, validation = utility_validation,
      assignment = assignments[[index]], original_assignment = original_assignment,
      true_fit = true_fit, draw_index = index,
      elapsed_seconds = outcome$elapsed_seconds,
      forced_status = if (outcome$feasible) NULL else "infeasible",
      forced_reason = if (length(outcome$failures)) {
        paste(outcome$failures, collapse = ";")
      } else NA_character_,
      feasible_override = outcome$feasible
    )
  })
  probability_draw_rows <- if (length(draw_rows)) {
    do.call(rbind, draw_rows)
  } else .ab_sim_utility_row(
    family_id, realization_id, seed,
    "stipulated_probability_integer_assignment_draws",
    validation = utility_validation, true_fit = true_fit
  )[FALSE, , drop = FALSE]
  utility_rows <- do.call(rbind, c(utility_rows, draw_rows))
  rownames(utility_rows) <- NULL
  rownames(probability_draw_rows) <- NULL

  target_effect <- generating_effects$target_log2_condition_effect[
    generating_effects$label == "A"
  ][[1L]]
  original_fit <- exact_timed$value$original_fit
  realization_summary <- data.frame(
    family_id = family_id,
    realization_id = realization_id,
    seed = seed,
    generating_effect_category = as.character(
      generating_truth$generating_effect_category
    ),
    generating_log2_effect_A = as.numeric(target_effect),
    probability_profile = as.character(generation_record$probability_profile),
    registered_ambiguous_cells = length(uncertain),
    true_label_status = if (is.null(true_fit)) "not_assessable" else true_fit$status,
    true_label_beta = if (!is.null(true_fit) &&
                         identical(true_fit$status, "evaluable")) {
      true_fit$beta
    } else NA_real_,
    original_status = if (is.null(original_fit)) "not_assessable" else original_fit$status,
    original_beta = if (!is.null(original_fit) &&
                       identical(original_fit$status, "evaluable")) {
      original_fit$beta
    } else NA_real_,
    exact_status = exact_timed$value$status,
    exact_minimum_cost = if (is.null(exact_timed$value$exact_minimum_cost) ||
                            !length(exact_timed$value$exact_minimum_cost)) {
      NA_integer_
    } else as.integer(exact_timed$value$exact_minimum_cost[[1L]]),
    exact_best_verified_cost = if (is.null(exact_timed$value$best_verified_cost) ||
                                  !length(exact_timed$value$best_verified_cost)) {
      NA_integer_
    } else as.integer(exact_timed$value$best_verified_cost[[1L]]),
    exact_total_assignments = exact_timed$value$total_assignment_count,
    exact_feasible_assignments = exact_timed$value$feasible_assignment_count,
    exact_evaluable_assignments = exact_timed$value$evaluable_assignment_count,
    exact_unevaluable_assignments = exact_timed$value$unevaluable_assignment_count,
    exact_attempted_fits = .ab_sim_fit_count(exact_timed$value),
    exact_failed_fits = .ab_sim_failed_fit_count(exact_timed$value),
    exact_search_seconds = exact_timing$search_seconds,
    exact_verification_fits = exact_checked$value$verification_fit_count,
    exact_verification_seconds = exact_timing$verification_seconds,
    exact_total_seconds = exact_timing$total_seconds,
    truth_mislabeled_cells = sum(cell_truth$is_mislabeled),
    truth_missing_correct_alternatives = sum(
      cell_truth$is_mislabeled & !cell_truth$correct_alternative_available
    ),
    stringsAsFactors = FALSE
  )
  structure(list(
    realization_summary = realization_summary,
    exact = exact_timed$value,
    exact_verification = exact_checked$value,
    exact_timing = exact_timing,
    practical = practical,
    outgoing = outgoing,
    utility_rows = utility_rows,
    probability_draw_rows = probability_draw_rows,
    true_label_fit = true_fit,
    true_label_validation = true_validation,
    generation_record = generation_record
  ), class = c("ab_simulation_evaluation", "list"))
}

.ab_sim_wilson <- function(successes, total, confidence) {
  if (!is.finite(total) || total <= 0) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt(p * (1 - p) / total + z^2 / (4 * total^2)) / denominator
  c(estimate = p, lower = max(0, center - half), upper = min(1, center + half))
}

.ab_sim_bootstrap_mean <- function(values, replicates, confidence, seed) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
  }
  if (length(values) == 1L || replicates < 1L) {
    return(c(estimate = mean(values), lower = mean(values),
             upper = mean(values), n = length(values)))
  }
  draws <- .ab_sim_with_local_seed(seed, function() {
    replicate(replicates, mean(sample(values, length(values), replace = TRUE)))
  })
  alpha <- (1 - confidence) / 2
  bounds <- stats::quantile(draws, probs = c(alpha, 1 - alpha),
                            names = FALSE, type = 8, na.rm = TRUE)
  c(estimate = mean(values), lower = bounds[[1L]], upper = bounds[[2L]],
    n = length(values))
}

.ab_sim_optimization_rows <- function(results) {
  rows <- list()
  at <- 1L
  for (evaluation in results) {
    summary <- evaluation$realization_summary[1L, , drop = FALSE]
    known_attainable <- identical(summary$exact_status[[1L]], "exact_minimum") &&
      is.finite(summary$exact_minimum_cost[[1L]])
    for (direction_name in c("all", "outgoing_only")) {
      runs <- if (direction_name == "all") evaluation$practical else evaluation$outgoing
      for (method in names(runs)) {
        run <- runs[[method]]
        result <- run$result
        verification <- run$verification
        hit <- identical(result$status, "verified_counterexample") &&
          isTRUE(verification$ok)
        rows[[at]] <- data.frame(
          family_id = summary$family_id,
          realization_id = summary$realization_id,
          seed = summary$seed,
          direction = direction_name,
          method = method,
          comparison_interpretation = if (direction_name == "all") {
            "same_objective_exact_comparison"
          } else {
            "coverage_of_bidirectional_exact_vulnerabilities"
          },
          exact_status = summary$exact_status,
          known_attainable = if (direction_name == "all") {
            known_attainable
          } else FALSE,
          bidirectional_reference_attainable = known_attainable,
          verified_hit = hit,
          exact_minimum_cost = if (direction_name == "all") {
            summary$exact_minimum_cost
          } else NA_integer_,
          verified_cost = if (hit) result$best_verified_cost else NA_integer_,
          exact_cost_gap = if (direction_name == "all" && hit && known_attainable) {
            result$best_verified_cost - summary$exact_minimum_cost
          } else NA_real_,
          search_status = result$status,
          termination_reason = result$termination_reason,
          proposal_count = result$proposal_count,
          cache_hit_count = result$cache_hit_count,
          attempted_fits = result$attempted_fit_count,
          successful_fits = result$successful_fit_count,
          failed_fits = result$failed_fit_count,
          verification_fits = verification$verification_fit_count,
          search_seconds = run$timing$search_seconds,
          verification_seconds = run$timing$verification_seconds,
          total_practical_runtime_seconds = run$timing$total_seconds,
          fit_failure_fraction = if (result$attempted_fit_count > 0L) {
            result$failed_fit_count / result$attempted_fit_count
          } else NA_real_,
          stringsAsFactors = FALSE
        )
        at <- at + 1L
      }
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

.ab_sim_empty_summary_row <- function(family_id, seed, realization_id,
                                      execution_status, error = NA_character_) {
  data.frame(
    family_id = as.character(family_id),
    realization_id = as.character(realization_id),
    seed = as.integer(seed),
    generating_effect_category = NA_character_,
    generating_log2_effect_A = NA_real_,
    probability_profile = NA_character_,
    registered_ambiguous_cells = NA_integer_,
    true_label_status = NA_character_,
    true_label_beta = NA_real_,
    original_status = NA_character_,
    original_beta = NA_real_,
    exact_status = NA_character_,
    exact_minimum_cost = NA_integer_,
    exact_best_verified_cost = NA_integer_,
    exact_total_assignments = NA_integer_,
    exact_feasible_assignments = NA_integer_,
    exact_evaluable_assignments = NA_integer_,
    exact_unevaluable_assignments = NA_integer_,
    exact_attempted_fits = NA_integer_,
    exact_failed_fits = NA_integer_,
    exact_search_seconds = NA_real_,
    exact_verification_fits = NA_integer_,
    exact_verification_seconds = NA_real_,
    exact_total_seconds = NA_real_,
    truth_mislabeled_cells = NA_integer_,
    truth_missing_correct_alternatives = NA_integer_,
    execution_status = as.character(execution_status),
    execution_error = as.character(error),
    stringsAsFactors = FALSE
  )
}

.ab_sim_normalize_metric_input <- function(realization_results) {
  if (inherits(realization_results, "ab_simulation_evaluation")) {
    realization_results <- list(realization_results)
  }
  if (!is.list(realization_results)) {
    stop("realization_results must be a list", call. = FALSE)
  }
  if (!all(c("registered_realizations", "latest_attempts") %in%
           names(realization_results))) {
    evaluations <- realization_results
    if (any(!vapply(evaluations, inherits, logical(1L),
                    what = "ab_simulation_evaluation"))) {
      stop("realization_results must contain simulation evaluations", call. = FALSE)
    }
    summaries <- lapply(evaluations, function(evaluation) {
      row <- evaluation$realization_summary
      row$execution_status <- "completed"
      row$execution_error <- NA_character_
      row
    })
    return(list(
      evaluations = evaluations,
      summaries = if (length(summaries)) do.call(rbind, summaries) else data.frame()
    ))
  }

  registered <- realization_results$registered_realizations
  if (!is.data.frame(registered) ||
      !all(c("family_id", "seed") %in% names(registered))) {
    stop("registered_realizations must have family_id and seed columns",
         call. = FALSE)
  }
  if (!"realization_id" %in% names(registered)) {
    registered$realization_id <- paste0(
      registered$family_id, "__seed_", as.integer(registered$seed)
    )
  }
  if (anyDuplicated(registered$realization_id)) {
    stop("registered_realizations contains duplicate realization IDs",
         call. = FALSE)
  }
  attempts <- realization_results$latest_attempts
  if (is.data.frame(attempts)) {
    attempts <- split(attempts, seq_len(nrow(attempts)))
  }
  if (!is.list(attempts)) {
    stop("latest_attempts must be a list or data frame", call. = FALSE)
  }
  attempt_id <- function(attempt) {
    if (inherits(attempt, "ab_simulation_evaluation")) {
      return(as.character(attempt$realization_summary$realization_id[[1L]]))
    }
    if (is.data.frame(attempt) && "realization_id" %in% names(attempt)) {
      return(as.character(attempt$realization_id[[1L]]))
    }
    if (is.list(attempt) && !is.null(attempt$evaluation) &&
        inherits(attempt$evaluation, "ab_simulation_evaluation")) {
      return(as.character(
        attempt$evaluation$realization_summary$realization_id[[1L]]
      ))
    }
    if (is.list(attempt) && !is.null(attempt$realization_id)) {
      return(as.character(attempt$realization_id[[1L]]))
    }
    NA_character_
  }
  ids <- vapply(attempts, attempt_id, character(1L))
  supplied_names <- names(attempts)
  if (!is.null(supplied_names)) {
    missing_id <- is.na(ids) | !nzchar(ids)
    ids[missing_id] <- supplied_names[missing_id]
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      any(!ids %in% registered$realization_id)) {
    stop("latest_attempts must map uniquely to registered realization IDs",
         call. = FALSE)
  }
  names(attempts) <- ids
  evaluations <- list()
  summaries <- vector("list", nrow(registered))
  for (index in seq_len(nrow(registered))) {
    id <- as.character(registered$realization_id[[index]])
    attempt <- attempts[[id]]
    evaluation <- if (inherits(attempt, "ab_simulation_evaluation")) {
      attempt
    } else if (is.list(attempt) &&
               inherits(attempt$evaluation, "ab_simulation_evaluation")) {
      attempt$evaluation
    } else NULL
    if (!is.null(evaluation)) {
      row <- evaluation$realization_summary
      row$execution_status <- "completed"
      row$execution_error <- NA_character_
      summaries[[index]] <- row
      evaluations[[id]] <- evaluation
    } else {
      status <- if (is.null(attempt)) "run_failure" else {
        candidate <- if (is.data.frame(attempt)) {
          if ("execution_status" %in% names(attempt)) attempt$execution_status[[1L]]
          else if ("status" %in% names(attempt)) attempt$status[[1L]]
          else "run_failure"
        } else if (!is.null(attempt$execution_status)) {
          attempt$execution_status[[1L]]
        } else if (!is.null(attempt$status)) attempt$status[[1L]] else "run_failure"
        as.character(candidate)
      }
      if (!status %in% c("generation_failure", "run_failure")) {
        stop("failed latest attempt must be generation_failure or run_failure",
             call. = FALSE)
      }
      error <- if (is.null(attempt)) "no latest attempt recorded" else if (
        is.data.frame(attempt) && "error" %in% names(attempt)
      ) attempt$error[[1L]] else if (is.list(attempt) && !is.null(attempt$error)) {
        attempt$error[[1L]]
      } else NA_character_
      summaries[[index]] <- .ab_sim_empty_summary_row(
        registered$family_id[[index]], registered$seed[[index]], id, status, error
      )
    }
  }
  list(evaluations = unname(evaluations), summaries = do.call(rbind, summaries))
}

.ab_summarize_simulation_metrics <- function(realization_results, metrics_spec) {
  if (!is.list(metrics_spec) ||
      !identical(as.character(metrics_spec$metrics_id),
                 "M3-optimization-biological-utility-v1")) {
    stop("metrics_spec must be M3-optimization-biological-utility-v1",
         call. = FALSE)
  }
  normalized <- .ab_sim_normalize_metric_input(realization_results)
  evaluations <- normalized$evaluations
  summaries <- normalized$summaries
  optimization_rows <- .ab_sim_optimization_rows(evaluations)
  utility_rows <- if (length(evaluations)) {
    do.call(rbind, lapply(evaluations, `[[`, "utility_rows"))
  } else data.frame()
  draw_rows <- if (length(evaluations)) {
    do.call(rbind, lapply(evaluations, `[[`, "probability_draw_rows"))
  } else data.frame()

  registered_status <- as.character(unlist(metrics_spec$status_denominators))
  observed_status <- if (nrow(summaries)) {
    unique(summaries$exact_status[summaries$execution_status == "completed"])
  } else character()
  unknown_status <- setdiff(observed_status[!is.na(observed_status)], registered_status)
  if (length(unknown_status)) {
    stop("unregistered exact status in realization results: ",
         paste(unknown_status, collapse = ", "), call. = FALSE)
  }
  non_result_status <- as.character(unlist(
    metrics_spec$non_result_failure_denominators
  ))
  if (!setequal(non_result_status, c("generation_failure", "run_failure"))) {
    stop("metrics_spec must register generation_failure and run_failure",
         call. = FALSE)
  }
  status_levels <- c(registered_status, non_result_status)
  status_denominators <- data.frame(
    realization_status = status_levels,
    count = vapply(status_levels, function(status) {
      if (status %in% non_result_status) {
        sum(summaries$execution_status == status)
      } else {
        sum(summaries$execution_status == "completed" &
              summaries$exact_status == status, na.rm = TRUE)
      }
    }, integer(1L)),
    total_registered_realizations = nrow(summaries),
    stringsAsFactors = FALSE
  )
  status_denominators$fraction <- if (nrow(summaries)) {
    status_denominators$count / nrow(summaries)
  } else NA_real_

  confidence <- as.numeric(metrics_spec$uncertainty$confidence_level)
  optimization_summary <- data.frame()
  optimization_family_summary <- data.frame()
  if (nrow(optimization_rows)) {
    summarize_optimization_block <- function(block, registered_n,
                                             family_id = NULL) {
      attainable <- block$known_attainable %in% TRUE
      reference_attainable <- block$bidirectional_reference_attainable %in% TRUE
      interval <- .ab_sim_wilson(sum(block$verified_hit & attainable),
                                 sum(attainable), confidence)
      coverage_interval <- .ab_sim_wilson(
        sum(block$verified_hit & reference_attainable),
        sum(reference_attainable), confidence
      )
      data.frame(
        family_id = if (is.null(family_id)) "fixed_registered_mixture" else family_id,
        direction = block$direction[[1L]],
        method = block$method[[1L]],
        comparison_interpretation = block$comparison_interpretation[[1L]],
        registered_realizations = as.integer(registered_n),
        executed_realizations = nrow(block),
        known_attainable_denominator = sum(attainable),
        known_attainable_hits = sum(block$verified_hit & attainable),
        known_attainable_hit_fraction = interval[["estimate"]],
        hit_fraction_wilson_lower = interval[["lower"]],
        hit_fraction_wilson_upper = interval[["upper"]],
        bidirectional_reference_attainable_denominator = sum(reference_attainable),
        bidirectional_vulnerability_coverage_hits = sum(
          block$verified_hit & reference_attainable
        ),
        bidirectional_vulnerability_coverage_fraction = coverage_interval[["estimate"]],
        coverage_wilson_lower = coverage_interval[["lower"]],
        coverage_wilson_upper = coverage_interval[["upper"]],
        verified_hit_count = sum(block$verified_hit),
        median_exact_cost_gap = if (any(is.finite(block$exact_cost_gap))) {
          stats::median(block$exact_cost_gap, na.rm = TRUE)
        } else NA_real_,
        total_search_fits = sum(block$attempted_fits),
        total_verification_fits = sum(block$verification_fits),
        total_failed_fits = sum(block$failed_fits),
        fit_failure_fraction = if (sum(block$attempted_fits) > 0) {
          sum(block$failed_fits) / sum(block$attempted_fits)
        } else NA_real_,
        median_total_runtime_seconds = stats::median(
          block$total_practical_runtime_seconds, na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
    groups <- split(optimization_rows,
                    interaction(optimization_rows$direction,
                                optimization_rows$method, drop = TRUE))
    optimization_summary <- do.call(rbind, lapply(groups, function(block) {
      summarize_optimization_block(block, nrow(summaries))
    }))
    rownames(optimization_summary) <- NULL
    family_groups <- split(
      optimization_rows,
      interaction(optimization_rows$family_id, optimization_rows$direction,
                  optimization_rows$method, drop = TRUE)
    )
    family_registered <- table(summaries$family_id)
    optimization_family_summary <- do.call(rbind, lapply(family_groups,
      function(block) {
        family_id <- as.character(block$family_id[[1L]])
        summarize_optimization_block(
          block, as.integer(family_registered[[family_id]]), family_id
        )
      }))
    rownames(optimization_family_summary) <- NULL
  }

  utility_by_realization <- data.frame()
  biological_utility_summary <- data.frame()
  if (nrow(utility_rows)) {
    utility_groups <- split(utility_rows,
      interaction(utility_rows$family_id, utility_rows$realization_id,
                  utility_rows$method, drop = TRUE))
    utility_by_realization <- do.call(rbind, lapply(utility_groups, function(block) {
      evaluable <- block$fit_status %in% "evaluable" & is.finite(block$beta)
      mean_beta <- if (any(evaluable)) mean(block$beta[evaluable]) else NA_real_
      true_beta <- unique(block$true_label_beta[is.finite(block$true_label_beta)])
      if (length(true_beta) > 1L) {
        stop("utility rows disagree on true-label beta within a realization",
             call. = FALSE)
      }
      true_beta <- if (length(true_beta)) true_beta[[1L]] else NA_real_
      data.frame(
        family_id = block$family_id[[1L]],
        realization_id = block$realization_id[[1L]],
        method = block$method[[1L]],
        method_outcomes = nrow(block),
        evaluable_outcomes = sum(evaluable),
        evaluable_fraction = mean(evaluable),
        mean_evaluable_beta = mean_beta,
        true_label_beta = true_beta,
        mean_beta_minus_true = if (is.finite(mean_beta) && is.finite(true_beta)) {
          mean_beta - true_beta
        } else NA_real_,
        absolute_error_of_mean_beta = if (is.finite(mean_beta) &&
                                          is.finite(true_beta)) {
          abs(mean_beta - true_beta)
        } else NA_real_,
        mean_absolute_beta_error = if (any(evaluable)) {
          mean(block$absolute_beta_error[evaluable], na.rm = TRUE)
        } else NA_real_,
        claim_agreement_fraction = if (any(evaluable)) {
          mean(block$claim_agreement[evaluable], na.rm = TRUE)
        } else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
    rownames(utility_by_realization) <- NULL
    summary_groups <- split(utility_by_realization,
      interaction(utility_by_realization$family_id,
                  utility_by_realization$method, drop = TRUE))
    bootstrap_n <- as.integer(metrics_spec$uncertainty$paired_differences$replicates)
    bootstrap_seed <- as.integer(metrics_spec$uncertainty$paired_differences$seed)
    biological_utility_summary <- do.call(rbind, lapply(seq_along(summary_groups),
      function(group_index) {
        block <- summary_groups[[group_index]]
        total_outcomes <- sum(block$method_outcomes)
        evaluated_outcomes <- sum(block$evaluable_outcomes)
        evaluable_interval <- .ab_sim_bootstrap_mean(
          block$evaluable_fraction, bootstrap_n, confidence,
          bootstrap_seed + 10000L + group_index
        )
        mean_beta_error_interval <- .ab_sim_bootstrap_mean(
          block$absolute_error_of_mean_beta, bootstrap_n, confidence,
          bootstrap_seed + 20000L + group_index
        )
        data.frame(
          family_id = block$family_id[[1L]],
          method = block$method[[1L]],
          independent_realizations = nrow(block),
          method_outcomes = total_outcomes,
          evaluable_outcomes = evaluated_outcomes,
          mean_realization_evaluable_fraction = evaluable_interval[["estimate"]],
          evaluable_fraction_bootstrap_lower = evaluable_interval[["lower"]],
          evaluable_fraction_bootstrap_upper = evaluable_interval[["upper"]],
          interval_unit = "independent_simulation_realization",
          interval_method = "realization_level_bootstrap_mean",
          bootstrap_replicates = bootstrap_n,
          bootstrap_seed = bootstrap_seed + 10000L + group_index,
          mean_absolute_error_of_realization_mean_beta =
            mean_beta_error_interval[["estimate"]],
          mean_beta_error_bootstrap_lower = mean_beta_error_interval[["lower"]],
          mean_beta_error_bootstrap_upper = mean_beta_error_interval[["upper"]],
          median_realization_mean_absolute_beta_error = if (
            any(is.finite(block$mean_absolute_beta_error))
          ) stats::median(block$mean_absolute_beta_error, na.rm = TRUE) else NA_real_,
          mean_realization_claim_agreement = if (
            any(is.finite(block$claim_agreement_fraction))
          ) mean(block$claim_agreement_fraction, na.rm = TRUE) else NA_real_,
          stringsAsFactors = FALSE
        )
      }))
    rownames(biological_utility_summary) <- NULL
  }

  draw_status_denominators <- data.frame()
  if (nrow(draw_rows)) {
    draw_status <- sort(unique(draw_rows$status), method = "radix")
    draw_status_denominators <- data.frame(
      status = draw_status,
      count = vapply(draw_status, function(status) {
        sum(draw_rows$status == status)
      }, integer(1L)),
      total_registered_draws = nrow(draw_rows),
      stringsAsFactors = FALSE
    )
    draw_status_denominators$fraction <-
      draw_status_denominators$count / nrow(draw_rows)
  }

  paired_bootstrap <- data.frame()
  if (nrow(optimization_rows)) {
    all_rows <- optimization_rows[optimization_rows$direction == "all", , drop = FALSE]
    method_pairs <- combn(sort(unique(all_rows$method), method = "radix"), 2L,
                          simplify = FALSE)
    bootstrap_n <- as.integer(metrics_spec$uncertainty$paired_differences$replicates)
    bootstrap_seed <- as.integer(metrics_spec$uncertainty$paired_differences$seed)
    collected <- list()
    at <- 1L
    for (pair_index in seq_along(method_pairs)) {
      pair <- method_pairs[[pair_index]]
      left <- all_rows[all_rows$method == pair[[1L]], , drop = FALSE]
      right <- all_rows[all_rows$method == pair[[2L]], , drop = FALSE]
      matched <- merge(left, right, by = "realization_id", suffixes = c("_left", "_right"))
      difference_values <- list(
        verified_hit = as.numeric(matched$verified_hit_right) -
          as.numeric(matched$verified_hit_left),
        exact_cost_gap = matched$exact_cost_gap_right - matched$exact_cost_gap_left,
        total_practical_runtime_seconds =
          matched$total_practical_runtime_seconds_right -
          matched$total_practical_runtime_seconds_left
      )
      for (metric_index in seq_along(difference_values)) {
        metric <- names(difference_values)[[metric_index]]
        interval <- .ab_sim_bootstrap_mean(
          difference_values[[metric]], bootstrap_n, confidence,
          bootstrap_seed + pair_index * 100L + metric_index
        )
        collected[[at]] <- data.frame(
          left_method = pair[[1L]], right_method = pair[[2L]], metric = metric,
          paired_realizations = as.integer(interval[["n"]]),
          mean_right_minus_left = interval[["estimate"]],
          bootstrap_lower = interval[["lower"]],
          bootstrap_upper = interval[["upper"]],
          bootstrap_replicates = bootstrap_n,
          bootstrap_seed = bootstrap_seed + pair_index * 100L + metric_index,
          stringsAsFactors = FALSE
        )
        at <- at + 1L
      }
    }
    if (length(collected)) paired_bootstrap <- do.call(rbind, collected)
  }

  structure(list(
    realization_summary = summaries,
    status_denominators = status_denominators,
    optimization_rows = optimization_rows,
    optimization_summary = optimization_summary,
    optimization_family_summary = optimization_family_summary,
    utility_rows = utility_rows,
    utility_by_realization = utility_by_realization,
    biological_utility_summary = biological_utility_summary,
    probability_draw_rows = draw_rows,
    probability_draw_status_denominators = draw_status_denominators,
    paired_bootstrap = paired_bootstrap,
    metrics_id = as.character(metrics_spec$metrics_id),
    registered_realization_count = nrow(summaries),
    completed_realization_count = sum(summaries$execution_status == "completed"),
    independent_realization_count = nrow(summaries)
  ), class = c("ab_simulation_metrics", "list"))
}
