.ab_empty_search_log <- function() {
  data.frame(
    proposal_index = integer(),
    proposal_type = character(),
    search_phase = character(),
    source_rank = integer(),
    proposal_rank = integer(),
    fit_rank = integer(),
    proxy_score = numeric(),
    parent_assignment_id = character(),
    assignment_id = character(),
    changes = character(),
    cost = integer(),
    feasible = logical(),
    constraint_failures = character(),
    cache_hit = logical(),
    duplicate_proposal = logical(),
    fit_attempted = logical(),
    fit_resumed = logical(),
    fit_status = character(),
    beta = numeric(),
    margin = numeric(),
    overturn = logical(),
    selected = logical(),
    fit_elapsed_seconds = numeric(),
    warnings = character(),
    failure_reason = character(),
    elapsed_seconds = numeric(),
    evaluation_error = character(),
    stringsAsFactors = FALSE
  )
}

.ab_search_skeleton <- function(status, reason, input_hashes = character()) {
  out <- .ab_result_skeleton(status, reason, input_hashes)
  out$method <- NA_character_
  out$seed <- NA_integer_
  out$direction <- NA_character_
  out$max_fits <- NA_integer_
  out$max_proposals <- NA_real_
  out$max_seconds <- NA_real_
  out$proposal_count <- 0L
  out$duplicate_proposal_count <- 0L
  out$cache_hit_count <- 0L
  out$unique_assignment_count <- 0L
  out$unique_proposal_count <- 0L
  out$infeasible_proposal_count <- 0L
  out$attempted_fit_count <- 0L
  out$successful_fit_count <- 0L
  out$failed_fit_count <- 0L
  out$verification_fit_count <- 0L
  out$candidate_block_size <- NA_integer_
  out$pruning_fit_reserve <- NA_integer_
  out$pruning_chunk_size <- NA_integer_
  out$proxy_target_pseudocount <- NA_real_
  out$proxy_library_pseudocount <- NA_real_
  out$proxy_scale <- NA_real_
  out$ranking_setup_count <- 0L
  out$ranking_setup_seconds <- 0
  out$ranking_peak_bytes <- 0
  out$component_accounting <- data.frame()
  out$witness_trajectory <- data.frame()
  out$first_witness <- NULL
  out$final_witness <- NULL
  out$fits_to_first_witness <- NA_integer_
  out$fits_to_final_witness <- NA_integer_
  out$time_to_first_full_fit_witness_seconds <- NA_real_
  out$time_to_final_full_fit_witness_seconds <- NA_real_
  out$pruning_complete <- NA
  out$search_censored <- FALSE
  out$resumed_fit_count <- 0L
  out$new_fit_count <- 0L
  out$cumulative_fit_count <- 0L
  out$resume_historical_fit_seconds <- 0
  out$resume_historical_search_elapsed_seconds <- 0
  out$progress_errors <- character()
  out$cost_one_early_stop_triggered <- FALSE
  out$cost_one_early_stop_replay_required <- FALSE
  out$cost_one_early_stop <- list(
    eligible_route = FALSE,
    triggered = FALSE,
    independent_replay_required = FALSE,
    independent_replay_status = "not_applicable",
    replay_failure_action = "not_applicable"
  )
  out$candidate_log <- .ab_empty_search_log()
  out$search_path <- data.frame()
  out$elapsed_seconds <- 0
  class(out) <- c("ab_search_result", "ab_result")
  out
}

.ab_direction_ok <- function(data, assignment, direction) {
  if (identical(direction, "all")) return(TRUE)
  changes <- .ab_change_rows(data, assignment)
  if (!nrow(changes)) return(TRUE)
  all(changes$original_label == data$claim$focal_label &
        changes$assigned_label == data$claim$alternative_label)
}

.ab_local_neighbors <- function(data, current, method, direction,
                                max_candidates = Inf) {
  cell_ids <- data$cells$cell_id
  options <- list()
  at <- 1L
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
      options[[at]] <- c(cell_id = cell_id, label = label)
      at <- at + 1L
    }
  }
  candidates <- list()
  truncated <- FALSE
  add_candidate <- function(assignment) {
    if (length(candidates) >= max_candidates) {
      truncated <<- TRUE
      return(FALSE)
    }
    if (.ab_direction_ok(data, assignment, direction)) {
      candidates[[length(candidates) + 1L]] <<- assignment
    }
    TRUE
  }
  if (length(options)) {
    for (move in options) {
      assignment <- current
      assignment[[move[["cell_id"]]]] <- move[["label"]]
      if (!add_candidate(assignment)) break
    }
  }
  if (identical(method, "pair") && length(options) >= 2L && !truncated) {
    for (first_index in seq_len(length(options) - 1L)) {
      for (second_index in seq.int(first_index + 1L, length(options))) {
        first <- options[[first_index]]
        second <- options[[second_index]]
        if (first[["cell_id"]] == second[["cell_id"]]) next
        assignment <- current
        assignment[[first[["cell_id"]]]] <- first[["label"]]
        assignment[[second[["cell_id"]]]] <- second[["label"]]
        if (!add_candidate(assignment)) break
      }
      if (truncated) break
    }
  }
  if (!length(candidates)) {
    return(structure(candidates, truncated = truncated))
  }
  ids <- vapply(candidates, function(assignment) {
    .ab_assignment_id(data, assignment)
  }, character(1L))
  keep <- !duplicated(ids)
  candidates <- candidates[keep]
  ids <- ids[keep]
  structure(candidates[order(ids, method = "radix")], truncated = truncated)
}

.ab_random_assignment <- function(data, direction) {
  original <- data$original_assignment
  eligible <- names(original)[vapply(names(original), function(cell_id) {
    alternatives <- setdiff(data$allowed_by_cell[[cell_id]], original[[cell_id]])
    if (identical(direction, "outgoing_only")) {
      alternatives <- alternatives[
        original[[cell_id]] == data$claim$focal_label &
          alternatives == data$claim$alternative_label
      ]
    }
    length(alternatives) > 0L
  }, logical(1L))]
  max_cost <- min(data$constraints$max_changes_global, length(eligible))
  if (max_cost < 1L) return(NULL)
  cost <- sample.int(max_cost, 1L)
  selected <- sample(eligible, cost, replace = FALSE)
  assignment <- original
  for (cell_id in selected) {
    alternatives <- setdiff(data$allowed_by_cell[[cell_id]], original[[cell_id]])
    if (identical(direction, "outgoing_only")) {
      alternatives <- alternatives[alternatives == data$claim$alternative_label]
    }
    assignment[[cell_id]] <- sample(alternatives, 1L)
  }
  assignment
}

#' Search for a verified annotation counterexample with bounded simple methods
#'
#' @param x A successful `ab_validation` object.
#' @param method Observed-centered proxy-ranked cumulative prefixes by default,
#'   or an explicit legacy/search-development route.
#' @param max_fits Positive full-fit cap including the original assignment.
#' @param seed Local reproducibility seed.
#' @param max_proposals Proposal cap excluding the original fit.
#' @param max_seconds Wall-time cap checked between proposals.
#' @param direction All moves or focal-population outgoing moves only.
#' @param fit_fun Optional custom complete-fit function.
#' @param fit_identity Required versioned label for M2 custom fit functions.
#' @param candidate_block_size Number of ranked candidates fitted before the
#'   ranked route either advances or scans the next block.
#' @param pruning_fit_reserve Full fits kept out of random discovery so they
#'   remain available for deletion pruning. A
#'   `NULL` value uses one quarter of the total cap, bounded by 32.
#' @param pruning_chunk_size Deterministic first-pass deletion chunk size. A
#'   `NULL` value uses 2 when `max_fits <= 32` and 8 otherwise.
#' @param proxy_target_pseudocount Target-count pseudocount in the ranking proxy.
#' @param proxy_library_pseudocount Library-count pseudocount in the ranking proxy.
#' @param proxy_scale Multiplicative CPM scale in the ranking proxy.
#' @param warm_start_assignment Optional complete assignment to refit before
#'   pruning; supported only by `method = "random_prune"`.
#' @param progress_callback Optional function receiving a compact progress record
#'   after every newly executed full fit.
#' @param progress_every Callback interval. The resumable ledger contract
#'   currently requires `1` so no completed full fit is omitted.
#' @param resume_evaluations Optional named list of compact prior fit records,
#'   keyed by canonical assignment ID, for deterministic route reconstruction.
#' @param prefix_sizes Optional strictly increasing positive integer prefix
#'   sizes for `method = "proxy_prefix"`; `NULL` uses the registered default.
#' @return An `ab_search_result` with upper-bound-only evidence semantics.
search_breakpoint <- function(x, method = c(
                                "proxy_prefix", "greedy", "random", "pair", "greedy_random",
                                "random_block", "random_prune", "ranked"
                              ),
                              max_fits, seed, max_proposals = Inf,
                              max_seconds = Inf,
                              direction = c("all", "outgoing_only"),
                              fit_fun = NULL, fit_identity = NULL,
                              candidate_block_size = 8L,
                              pruning_fit_reserve = NULL,
                              pruning_chunk_size = NULL,
                              proxy_target_pseudocount = 0.5,
                              proxy_library_pseudocount = 1,
                              proxy_scale = 1e6,
                              warm_start_assignment = NULL,
                              progress_callback = NULL,
                              progress_every = 1L,
                              resume_evaluations = NULL,
                              prefix_sizes = NULL) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  invalid <- character()
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    invalid <- c(invalid, "x must be a successful ab_validation object")
  }
  if (missing(max_fits) || length(max_fits) != 1L || is.na(max_fits) ||
      !((is.infinite(max_fits) && max_fits > 0) ||
        (.ab_is_whole(max_fits) && max_fits >= 1L))) {
    invalid <- c(invalid, "max_fits must be +Inf or a positive finite integer")
  }
  if (missing(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    invalid <- c(invalid, "seed must be one finite integer")
  }
  valid_cap <- function(value) {
    length(value) == 1L && !is.na(value) &&
      ((is.infinite(value) && value > 0) || .ab_is_whole(value))
  }
  if (!valid_cap(max_proposals)) {
    invalid <- c(invalid, "max_proposals must be +Inf or a nonnegative integer")
  }
  if (length(max_seconds) != 1L || is.na(max_seconds) ||
      (!is.infinite(max_seconds) && (!is.finite(max_seconds) || max_seconds < 0)) ||
      (is.infinite(max_seconds) && max_seconds < 0)) {
    invalid <- c(invalid, "max_seconds must be +Inf or finite and nonnegative")
  }
  valid_positive_integer <- function(value) {
    length(value) == 1L && !is.na(value) && .ab_is_whole(value) && value >= 1L
  }
  if (!valid_positive_integer(candidate_block_size)) {
    invalid <- c(invalid, "candidate_block_size must be a positive finite integer")
  }
  if (!is.null(pruning_fit_reserve) &&
      (length(pruning_fit_reserve) != 1L || is.na(pruning_fit_reserve) ||
       !.ab_is_whole(pruning_fit_reserve) || pruning_fit_reserve < 0L)) {
    invalid <- c(invalid, "pruning_fit_reserve must be NULL or a nonnegative finite integer")
  }
  if (!is.null(pruning_chunk_size) &&
      !valid_positive_integer(pruning_chunk_size)) {
    invalid <- c(invalid, "pruning_chunk_size must be NULL or a positive finite integer")
  }
  valid_positive_number <- function(value) {
    length(value) == 1L && !is.na(value) && is.finite(value) && value > 0
  }
  if (!valid_positive_number(proxy_target_pseudocount) ||
      !valid_positive_number(proxy_library_pseudocount) ||
      !valid_positive_number(proxy_scale)) {
    invalid <- c(
      invalid,
      "ranking proxy pseudocounts and scale must be positive finite numbers"
    )
  }
  if (!is.null(warm_start_assignment) &&
      !identical(method, "random_prune")) {
    invalid <- c(
      invalid,
      "warm_start_assignment is supported only for method = 'random_prune'"
    )
  }
  if (identical(method, "proxy_prefix") && !identical(direction, "all")) {
    invalid <- c(invalid, "proxy_prefix requires direction = 'all'")
  }
  if (!is.null(prefix_sizes) && !identical(method, "proxy_prefix")) {
    invalid <- c(invalid, "prefix_sizes is supported only by method = 'proxy_prefix'")
  }
  if (identical(method, "proxy_prefix") && !is.null(prefix_sizes)) {
    prefix_error <- tryCatch({
      .ab_proxy_prefix_sizes(prefix_sizes)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(prefix_error)) invalid <- c(invalid, prefix_error)
  }
  if (identical(method, "proxy_prefix") &&
      (!identical(as.integer(candidate_block_size), 8L) ||
       !is.null(pruning_fit_reserve) || !is.null(pruning_chunk_size) ||
       !identical(as.numeric(proxy_target_pseudocount), 0.5) ||
       !identical(as.numeric(proxy_library_pseudocount), 1) ||
       !identical(as.numeric(proxy_scale), 1e6) ||
       !is.null(progress_callback) || !is.null(resume_evaluations))) {
    invalid <- c(
      invalid,
      paste0(
        "proxy_prefix uses its fixed registered ranking and does not accept ",
        "block, pruning, proxy, progress-callback, or resume overrides"
      )
    )
  }
  if (!is.null(progress_callback) && !is.function(progress_callback)) {
    invalid <- c(invalid, "progress_callback must be NULL or a function")
  }
  if (!valid_positive_integer(progress_every) || progress_every != 1L) {
    invalid <- c(invalid, "progress_every must be 1 for complete resumable ledgers")
  }
  if (!is.null(resume_evaluations) &&
      (!is.list(resume_evaluations) || is.null(names(resume_evaluations)) ||
       anyNA(names(resume_evaluations)) || any(!nzchar(names(resume_evaluations))) ||
       anyDuplicated(names(resume_evaluations)))) {
    invalid <- c(
      invalid,
      "resume_evaluations must be NULL or a uniquely named list keyed by assignment ID"
    )
  }
  if (length(resume_evaluations) &&
      !method %in% c("greedy_random", "random_block", "random_prune", "ranked")) {
    invalid <- c(invalid, "resume_evaluations is supported only by revised routes")
  }
  hashes <- if (inherits(x, "ab_validation") && isTRUE(x$ok)) {
    .ab_recompute_input_hashes(x)
  } else character()
  if (length(invalid)) {
    result <- .ab_search_skeleton("invalid_input", "invalid_search_arguments", hashes)
    result$errors <- invalid
    return(result)
  }
  candidate_block_size <- as.integer(candidate_block_size)
  effective_pruning_fit_reserve <- if (is.null(pruning_fit_reserve)) {
    if (is.infinite(max_fits)) 0L else min(32L, as.integer(floor(max_fits / 4)))
  } else {
    as.integer(pruning_fit_reserve)
  }
  if (is.finite(max_fits)) {
    effective_pruning_fit_reserve <- min(
      effective_pruning_fit_reserve, max(0L, as.integer(max_fits) - 1L)
    )
  }
  effective_pruning_chunk_size <- if (is.null(pruning_chunk_size)) {
    if (max_fits <= 32L) 2L else 8L
  } else {
    as.integer(pruning_chunk_size)
  }
  data <- x$data
  data$input_hashes <- hashes
  normalized_warm_start <- NULL
  if (!is.null(warm_start_assignment)) {
    normalized_warm_start <- tryCatch(
      .ab_normalize_assignment(data, warm_start_assignment),
      error = function(e) e
    )
    if (inherits(normalized_warm_start, "error")) {
      result <- .ab_search_skeleton(
        "invalid_input", "invalid_warm_start_assignment", hashes
      )
      result$errors <- conditionMessage(normalized_warm_start)
      result$method <- method
      result$seed <- as.integer(seed)
      result$direction <- direction
      return(.ab_attach_contract(result, data))
    }
  }
  resolved <- tryCatch(
    .ab_resolve_fit(x, fit_fun, fit_identity), error = function(e) e
  )
  if (inherits(resolved, "error")) {
    result <- .ab_search_skeleton("invalid_input", "invalid_fit_function", hashes)
    result$errors <- conditionMessage(resolved)
    result$method <- method
    result$seed <- as.integer(seed)
    result$direction <- direction
    return(.ab_attach_contract(result, data))
  }
  if (identical(method, "proxy_prefix")) {
    return(.ab_search_proxy_prefix(
      x, data, hashes, resolved, max_fits, seed, max_proposals, max_seconds,
      prefix_sizes = prefix_sizes
    ))
  }
  resume_fit_by_id <- new.env(hash = TRUE, parent = emptyenv())
  if (length(resume_evaluations)) {
    resume_error <- NULL
    for (assignment_id in names(resume_evaluations)) {
      fit <- tryCatch(
        .ab_resume_fit(
          resume_evaluations[[assignment_id]], assignment_id, resolved$identity
        ),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        resume_error <- paste0(assignment_id, ": ", conditionMessage(fit))
        break
      }
      assign(assignment_id, fit, envir = resume_fit_by_id)
    }
    if (!is.null(resume_error)) {
      result <- .ab_search_skeleton(
        "invalid_input", "invalid_resume_evaluations", hashes
      )
      result$errors <- resume_error
      result$method <- method
      result$seed <- as.integer(seed)
      result$direction <- direction
      return(.ab_attach_contract(result, data, resolved$identity))
    }
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  started <- unname(proc.time()[["elapsed"]])
  random_domain <- NULL
  if (method %in% c("random", "random_prune")) {
    finite_domain <- tryCatch(.ab_candidate_domain(data), error = function(e) e)
    if (inherits(finite_domain, "error")) {
      if (identical(method, "random") && is.infinite(max_proposals)) {
        result <- .ab_search_skeleton(
          "invalid_input", "finite_max_proposals_required_for_large_random_domain",
          hashes
        )
        result$errors <- paste0(
          "random search above the small enumerable domain requires finite max_proposals: ",
          conditionMessage(finite_domain)
        )
        result$method <- method
        result$seed <- as.integer(seed)
        result$direction <- direction
        result$max_fits <- if (is.infinite(max_fits)) Inf else as.integer(max_fits)
        result$max_proposals <- max_proposals
        result$max_seconds <- max_seconds
        return(.ab_attach_contract(result, data, resolved$identity))
      }
    } else {
      eligible_index <- which(
        finite_domain$table$cost > 0L &
          finite_domain$table$cost <= data$constraints$max_changes_global &
          vapply(finite_domain$table$assignment, function(assignment) {
            .ab_direction_ok(data, assignment, direction)
          }, logical(1L))
      )
      random_domain <- finite_domain$table$assignment[eligible_index]
      if (length(random_domain) > 1L) {
        random_domain <- random_domain[sample.int(length(random_domain))]
      }
    }
  }

  cache <- new.env(hash = TRUE, parent = emptyenv())
  seen <- new.env(hash = TRUE, parent = emptyenv())
  state <- new.env(parent = emptyenv())
  lightweight_route <- method %in% c(
    "greedy_random", "random_block", "random_prune", "ranked"
  )
  state$proposal_count <- 0L
  state$duplicate_count <- 0L
  state$cache_hits <- 0L
  state$infeasible <- 0L
  state$attempted <- 0L
  state$successful <- 0L
  state$failed <- 0L
  state$cumulative_successful <- 0L
  state$cumulative_failed <- 0L
  state$candidate_fit_count <- 0L
  state$new_fit_count <- 0L
  state$resumed_fit_count <- 0L
  state$resume_historical_fit_seconds <- 0
  state$resume_historical_search_elapsed_seconds <- 0
  state$current_phase <- "original"
  state$phase_started <- unname(proc.time()[["elapsed"]] - started)
  state$phase_seconds <- stats::setNames(
    numeric(4L), c("original", "discovery", "ranking", "pruning")
  )
  if (method %in% c("random", "random_prune")) {
    state$phase_seconds[["discovery"]] <- state$phase_started
  }
  state$ranking_setup_count <- 0L
  state$ranking_setup_seconds <- 0
  state$ranking_peak_bytes <- 0
  state$progress_errors <- character()
  log_rows <- list()
  witness_events <- list()
  first_witness_record <- NULL
  cost_one_early_stop_record <- NULL

  elapsed <- function() unname(proc.time()[["elapsed"]] - started)
  switch_phase <- function(phase) {
    now <- elapsed()
    old <- state$current_phase
    state$phase_seconds[[old]] <- state$phase_seconds[[old]] +
      now - state$phase_started
    state$current_phase <- phase
    state$phase_started <- now
    invisible(NULL)
  }
  close_phase <- function() {
    now <- elapsed()
    phase <- state$current_phase
    state$phase_seconds[[phase]] <- state$phase_seconds[[phase]] +
      now - state$phase_started
    state$phase_started <- now
    invisible(NULL)
  }
  append_log <- function(assignment, proposal_type, parent_id, cache_hit,
                         duplicate, fit_attempted, fit, failures,
                         selected = FALSE, proposal_index = state$proposal_count,
                         source_rank = NA_integer_, proxy_score = NA_real_,
                         fit_rank = NA_integer_, fit_resumed = FALSE) {
    changes <- .ab_change_rows(data, assignment)
    evaluable <- inherits(fit, "ab_fit") && identical(fit$status, "evaluable")
    row <- data.frame(
      proposal_index = as.integer(proposal_index),
      proposal_type = proposal_type,
      search_phase = state$current_phase,
      source_rank = as.integer(source_rank),
      proposal_rank = as.integer(proposal_index),
      fit_rank = as.integer(fit_rank),
      proxy_score = as.numeric(proxy_score),
      parent_assignment_id = if (is.null(parent_id)) NA_character_ else parent_id,
      assignment_id = .ab_assignment_id(data, assignment),
      changes = .ab_changes_text(changes),
      cost = nrow(changes),
      feasible = !length(failures),
      constraint_failures = paste(failures, collapse = ";"),
      cache_hit = cache_hit,
      duplicate_proposal = duplicate,
      fit_attempted = fit_attempted,
      fit_resumed = fit_resumed,
      fit_status = if (is.null(fit)) NA_character_ else fit$status,
      beta = if (evaluable) fit$beta else NA_real_,
      margin = if (evaluable) fit$margin else NA_real_,
      overturn = if (evaluable) fit$overturn else FALSE,
      selected = selected,
      fit_elapsed_seconds = if (!is.null(fit) &&
                                !is.null(fit$fit_elapsed_seconds)) {
        as.numeric(fit$fit_elapsed_seconds)
      } else NA_real_,
      warnings = if (!is.null(fit) && !is.null(fit$warnings) &&
                     length(fit$warnings)) {
        paste(unique(as.character(fit$warnings)), collapse = " | ")
      } else NA_character_,
      failure_reason = if (!is.null(fit) && !is.null(fit$failure_reason) &&
                           length(fit$failure_reason) == 1L &&
                           !is.na(fit$failure_reason)) {
        as.character(fit$failure_reason)
      } else NA_character_,
      elapsed_seconds = elapsed(),
      evaluation_error = if (is.null(fit) || evaluable) NA_character_ else fit$error,
      stringsAsFactors = FALSE
    )
    row$assignment <- I(list(if (lightweight_route) NULL else assignment))
    log_rows[[length(log_rows) + 1L]] <<- row
    length(log_rows)
  }

  evaluate_assignment <- function(assignment, proposal_type, parent_id = NULL,
                                  is_original = FALSE, source_rank = NA_integer_,
                                  proxy_score = NA_real_) {
    assignment <- .ab_normalize_assignment(data, assignment)
    assignment_id <- .ab_assignment_id(data, assignment)
    key <- .ab_fit_cache_key(data, assignment, resolved$identity)
    duplicate <- exists(assignment_id, envir = seen, inherits = FALSE)
    if (!duplicate) assign(assignment_id, TRUE, envir = seen)
    if (!is_original) {
      state$proposal_count <- state$proposal_count + 1L
      if (duplicate) state$duplicate_count <- state$duplicate_count + 1L
    }
    if (exists(key, envir = cache, inherits = FALSE)) {
      state$cache_hits <- state$cache_hits + as.integer(!is_original)
      record <- get(key, envir = cache, inherits = FALSE)
      record$assignment <- assignment
      if (!is_original && length(record$failures)) {
        state$infeasible <- state$infeasible + 1L
      }
      index <- append_log(
        assignment, proposal_type, parent_id, TRUE, duplicate, FALSE,
        record$fit, record$failures,
        proposal_index = if (is_original) 0L else state$proposal_count,
        source_rank = source_rank, proxy_score = proxy_score
      )
      return(c(record, list(log_index = index, cap = FALSE)))
    }
    failures <- .ab_constraint_failures(data, assignment)
    if (!.ab_direction_ok(data, assignment, direction)) {
      failures <- c(failures, "direction:outgoing_only")
    }
    failures <- unique(failures)
    if (length(failures)) {
      if (!is_original) state$infeasible <- state$infeasible + 1L
      record <- list(fit = NULL, failures = failures, assignment = assignment)
      cached_record <- record
      if (lightweight_route) cached_record$assignment <- NULL
      assign(key, cached_record, envir = cache)
      index <- append_log(
        assignment, proposal_type, parent_id, FALSE, duplicate, FALSE,
        NULL, failures,
        proposal_index = if (is_original) 0L else state$proposal_count,
        source_rank = source_rank, proxy_score = proxy_score
      )
      return(c(record, list(log_index = index, cap = FALSE)))
    }
    if (state$attempted >= max_fits) {
      index <- append_log(
        assignment, proposal_type, parent_id, FALSE, duplicate, FALSE,
        NULL, character(), proposal_index = state$proposal_count,
        source_rank = source_rank, proxy_score = proxy_score
      )
      return(list(
        fit = NULL, failures = character(), assignment = assignment,
        log_index = index, cap = TRUE
      ))
    }
    resumed_fit <- exists(assignment_id, envir = resume_fit_by_id, inherits = FALSE)
    state$attempted <- state$attempted + 1L
    fit_rank <- if (is_original) {
      0L
    } else {
      state$candidate_fit_count <- state$candidate_fit_count + 1L
      state$candidate_fit_count
    }
    if (resumed_fit) {
      fit <- get(assignment_id, envir = resume_fit_by_id, inherits = FALSE)
      state$resumed_fit_count <- state$resumed_fit_count + 1L
      if (!is.null(fit$fit_elapsed_seconds) &&
          length(fit$fit_elapsed_seconds) == 1L &&
          is.finite(fit$fit_elapsed_seconds)) {
        state$resume_historical_fit_seconds <-
          state$resume_historical_fit_seconds + fit$fit_elapsed_seconds
      }
      if (!is.null(fit$search_elapsed_seconds) &&
          length(fit$search_elapsed_seconds) == 1L &&
          is.finite(fit$search_elapsed_seconds)) {
        state$resume_historical_search_elapsed_seconds <- max(
          state$resume_historical_search_elapsed_seconds,
          fit$search_elapsed_seconds
        )
      }
    } else {
      value <- tryCatch(
        .ab_call_fit_fun(resolved$fun, x, assignment), error = function(e) e
      )
      fit <- .ab_standardize_fit(value, data, assignment, assignment_id)
      fit$search_elapsed_seconds <- elapsed()
      fit$historical_fit_index <- as.integer(fit_rank)
      fit$historical_proposal_index <- if (is_original) {
        0L
      } else as.integer(state$proposal_count)
      state$new_fit_count <- state$new_fit_count + 1L
    }
    if (identical(fit$status, "evaluable")) {
      state$cumulative_successful <- state$cumulative_successful + 1L
      if (!resumed_fit) state$successful <- state$successful + 1L
    } else {
      state$cumulative_failed <- state$cumulative_failed + 1L
      if (!resumed_fit) state$failed <- state$failed + 1L
    }
    if (lightweight_route) fit <- .ab_compact_search_fit(fit)
    record <- list(fit = fit, failures = character(), assignment = assignment)
    cached_record <- record
    if (lightweight_route) cached_record$assignment <- NULL
    assign(key, cached_record, envir = cache)
    index <- append_log(
      assignment, proposal_type, parent_id, FALSE, duplicate, !resumed_fit,
      fit, character(),
      proposal_index = if (is_original) 0L else state$proposal_count,
      source_rank = source_rank, proxy_score = proxy_score,
      fit_rank = fit_rank, fit_resumed = resumed_fit
    )
    record <- c(record, list(log_index = index, cap = FALSE))
    if (!is_original && is.null(first_witness_record) &&
        identical(fit$status, "evaluable") && isTRUE(fit$overturn)) {
      first_witness_record <<- record
      witness_events[[length(witness_events) + 1L]] <<- data.frame(
        event_index = as.integer(length(witness_events) + 1L),
        event = "first_full_fit_witness",
        assignment_id = assignment_id,
        cost = as.integer(nrow(.ab_change_rows(data, assignment))),
        beta = as.numeric(fit$beta), margin = as.numeric(fit$margin),
        cumulative_best_cost = as.integer(nrow(.ab_change_rows(data, assignment))),
        attempted_fit_count = as.integer(state$attempted),
        proposal_count = as.integer(state$proposal_count),
        elapsed_seconds = if (!is.null(fit$search_elapsed_seconds) &&
                              is.finite(fit$search_elapsed_seconds)) {
          as.numeric(fit$search_elapsed_seconds)
        } else elapsed(),
        is_final = FALSE,
        stringsAsFactors = FALSE
      )
    }
    if (!resumed_fit && !is.null(progress_callback)) {
      progress_error <- tryCatch({
        progress_callback(list(
          schema = "annotation-breakpoint-search-progress-v1",
          input_hashes = as.list(hashes),
          evaluation_identity = resolved$identity,
          method = method,
          seed = as.integer(seed),
          phase = state$current_phase,
          proposal_count = as.integer(state$proposal_count),
          candidate_fit_count = as.integer(state$candidate_fit_count),
          attempted_fit_count = as.integer(state$attempted),
          new_fit_count = as.integer(state$new_fit_count),
          resumed_fit_count = as.integer(state$resumed_fit_count),
          elapsed_seconds = elapsed(),
          last_assignment_id = assignment_id,
          last_status = fit$status,
          fit = unclass(.ab_compact_search_fit(fit)),
          first_witness_assignment_id = if (is.null(first_witness_record)) {
            NA_character_
          } else .ab_assignment_id(data, first_witness_record$assignment)
        ))
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(progress_error)) {
        state$progress_errors <- c(state$progress_errors, progress_error)
      }
    }
    record
  }
  record_witness_event <- function(record, event) {
    witness_events[[length(witness_events) + 1L]] <<- data.frame(
      event_index = as.integer(length(witness_events) + 1L),
      event = as.character(event),
      assignment_id = .ab_assignment_id(data, record$assignment),
      cost = as.integer(nrow(.ab_change_rows(data, record$assignment))),
      beta = as.numeric(record$fit$beta),
      margin = as.numeric(record$fit$margin),
      cumulative_best_cost = as.integer(min(c(
        vapply(witness_events, function(row) row$cost[[1L]], integer(1L)),
        nrow(.ab_change_rows(data, record$assignment))
      ))),
      attempted_fit_count = as.integer(state$attempted),
      proposal_count = as.integer(state$proposal_count),
      elapsed_seconds = if (!is.null(record$fit$search_elapsed_seconds) &&
                            is.finite(record$fit$search_elapsed_seconds)) {
        as.numeric(record$fit$search_elapsed_seconds)
      } else elapsed(),
      is_final = FALSE,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  move_record_assignment <- function(record, base_assignment) {
    assignment <- base_assignment
    assignment[[record$move$cell_id]] <- record$move$to_label
    assignment
  }

  original_record <- evaluate_assignment(
    data$original_assignment, "original", is_original = TRUE
  )
  original_fit <- original_record$fit
  is_cost_one_early_stop <- function(record, assignment) {
    method %in% c("random_block", "ranked") &&
      inherits(original_fit, "ab_fit") &&
      identical(original_fit$status, "evaluable") &&
      isTRUE(original_fit$claim_holds) &&
      !isTRUE(record$cap) && !length(record$failures) &&
      inherits(record$fit, "ab_fit") &&
      identical(record$fit$status, "evaluable") &&
      isTRUE(record$fit$overturn) &&
      nrow(.ab_change_rows(data, assignment)) == 1L
  }
  result_status <- "no_counterexample_found"
  termination <- ""
  witness_record <- NULL
  pruning_complete <- NA
  pruning_log <- list()
  source_candidate_count <- if (!is.null(random_domain)) {
    as.integer(length(random_domain))
  } else 0L
  path_rows <- list()
  if (is.null(original_fit) || !identical(original_fit$status, "evaluable")) {
    result_status <- "not_assessable"
    termination <- "original_fit_not_assessable"
  } else if (!isTRUE(original_fit$claim_holds)) {
    result_status <- "original_claim_not_met"
    termination <- "original_claim_not_met"
  } else {
    current <- data$original_assignment
    current_fit <- original_fit
    current_id <- .ab_assignment_id(data, current)
    path_rows[[1L]] <- data.frame(
      step = 0L, assignment_id = current_id, cost = 0L,
      beta = current_fit$beta, margin = current_fit$margin,
      stringsAsFactors = FALSE
    )
    log_rows[[original_record$log_index]]$selected <- TRUE
    switch_phase("discovery")

    if (identical(method, "random")) {
      random_index <- 1L
      repeat {
        if (state$proposal_count >= max_proposals) {
          termination <- "max_proposals_reached"
          break
        }
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        if (!is.null(random_domain)) {
          if (random_index > length(random_domain)) {
            termination <- if (length(random_domain)) {
              "random_allowed_domain_exhausted"
            } else "no_eligible_changes"
            break
          }
          proposal <- random_domain[[random_index]]
          random_index <- random_index + 1L
        } else {
          proposal <- .ab_random_assignment(data, direction)
          if (is.null(proposal)) {
            termination <- "no_eligible_changes"
            break
          }
        }
        record <- evaluate_assignment(proposal, "random", current_id)
        if (isTRUE(record$cap)) {
          termination <- "max_fits_reached"
          break
        }
        if (!is.null(record$fit) && identical(record$fit$status, "evaluable") &&
            isTRUE(record$fit$overturn)) {
          witness_record <- record
          log_rows[[record$log_index]]$selected <- TRUE
          result_status <- "verified_counterexample"
          termination <- "verified_counterexample_found"
          break
        }
      }
    } else if (method %in% c("greedy", "pair")) {
      step <- 0L
      repeat {
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        proposal_remaining <- if (is.infinite(max_proposals)) {
          Inf
        } else max(0, max_proposals - state$proposal_count)
        fit_remaining <- max(0, max_fits - state$attempted)
        safe_neighbor_cap <- max(
          1,
          min(10000, fit_remaining + length(ls(cache, all.names = TRUE)))
        )
        neighbor_cap <- min(proposal_remaining, safe_neighbor_cap)
        if (neighbor_cap < 1L) {
          termination <- "max_proposals_reached"
          break
        }
        neighbors <- .ab_local_neighbors(
          data, current, method, direction, max_candidates = neighbor_cap
        )
        neighborhood_truncated <- isTRUE(attr(neighbors, "truncated"))
        if (!length(neighbors)) {
          termination <- "no_eligible_neighbors"
          break
        }
        records <- list()
        cap_reason <- NULL
        for (proposal in neighbors) {
          if (state$proposal_count >= max_proposals) {
            cap_reason <- "max_proposals_reached"
            break
          }
          if (elapsed() >= max_seconds) {
            cap_reason <- "max_seconds_reached"
            break
          }
          record <- evaluate_assignment(
            proposal, paste0(method, "_local"), current_id
          )
          records[[length(records) + 1L]] <- record
          if (isTRUE(record$cap)) {
            cap_reason <- "max_fits_reached"
            break
          }
        }
        evaluable <- which(vapply(records, function(record) {
          !is.null(record$fit) && identical(record$fit$status, "evaluable")
        }, logical(1L)))
        if (!length(evaluable)) {
          termination <- if (is.null(cap_reason)) {
            "no_evaluable_neighbor"
          } else cap_reason
          break
        }
        overturn_records <- evaluable[vapply(records[evaluable], function(record) {
          isTRUE(record$fit$overturn)
        }, logical(1L))]
        if (length(overturn_records)) {
          overturn_cost <- vapply(records[overturn_records], function(record) {
            nrow(.ab_change_rows(data, record$assignment))
          }, integer(1L))
          overturn_id <- vapply(records[overturn_records], function(record) {
            .ab_assignment_id(data, record$assignment)
          }, character(1L))
          selected_order <- order(overturn_cost, overturn_id, method = "radix")
          best <- records[[overturn_records[selected_order[[1L]]]]]
          witness_record <- best
          log_rows[[best$log_index]]$selected <- TRUE
          step <- step + 1L
          path_rows[[length(path_rows) + 1L]] <- data.frame(
            step = step,
            assignment_id = .ab_assignment_id(data, best$assignment),
            cost = nrow(.ab_change_rows(data, best$assignment)),
            beta = best$fit$beta,
            margin = best$fit$margin,
            stringsAsFactors = FALSE
          )
          result_status <- "verified_counterexample"
          termination <- "verified_counterexample_found"
          break
        }
        margins <- vapply(records[evaluable], function(record) record$fit$margin,
                          numeric(1L))
        ids <- vapply(records[evaluable], function(record) {
          .ab_assignment_id(data, record$assignment)
        }, character(1L))
        best_order <- order(margins, ids, method = "radix")
        best <- records[[evaluable[best_order[[1L]]]]]
        if (!(best$fit$margin < current_fit$margin - data$claim$tolerance)) {
          termination <- if (!is.null(cap_reason)) {
            cap_reason
          } else if (neighborhood_truncated) {
            "neighborhood_candidate_cap_reached"
          } else "no_improving_move"
          break
        }
        current <- best$assignment
        current_fit <- best$fit
        current_id <- .ab_assignment_id(data, current)
        step <- step + 1L
        log_rows[[best$log_index]]$selected <- TRUE
        path_rows[[length(path_rows) + 1L]] <- data.frame(
          step = step, assignment_id = current_id,
          cost = nrow(.ab_change_rows(data, current)),
          beta = current_fit$beta, margin = current_fit$margin,
          stringsAsFactors = FALSE
        )
        if (!is.null(cap_reason)) {
          termination <- cap_reason
          break
        }
      }
    } else if (identical(method, "greedy_random")) {
      step <- 0L
      repeat {
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        proposal_remaining <- if (is.infinite(max_proposals)) {
          Inf
        } else max(0, max_proposals - state$proposal_count)
        fit_remaining <- max(0, max_fits - state$attempted)
        safe_neighbor_cap <- max(
          1, fit_remaining + length(ls(cache, all.names = TRUE))
        )
        neighbor_cap <- min(proposal_remaining, safe_neighbor_cap)
        if (neighbor_cap < 1L) {
          termination <- "max_proposals_reached"
          break
        }
        full_source <- .ab_randomized_single_source(data, current, direction)
        source_candidate_count <- max(source_candidate_count, length(full_source))
        neighborhood_truncated <- length(full_source) > neighbor_cap
        neighbors <- if (length(full_source)) {
          full_source[seq_len(min(length(full_source), neighbor_cap))]
        } else list()
        if (!length(neighbors)) {
          termination <- "no_eligible_neighbors"
          break
        }
        records <- list()
        cap_reason <- NULL
        for (candidate in neighbors) {
          if (state$proposal_count >= max_proposals) {
            cap_reason <- "max_proposals_reached"
            break
          }
          if (elapsed() >= max_seconds) {
            cap_reason <- "max_seconds_reached"
            break
          }
          proposal <- current
          proposal[[candidate$cell_id]] <- candidate$to_label
          record <- evaluate_assignment(
            proposal, "greedy_random_local", current_id,
            source_rank = candidate$source_rank
          )
          record$move <- candidate
          record$assignment <- NULL
          records[[length(records) + 1L]] <- record
          if (isTRUE(record$cap)) {
            cap_reason <- "max_fits_reached"
            break
          }
        }
        evaluable <- which(vapply(records, function(record) {
          !is.null(record$fit) && identical(record$fit$status, "evaluable")
        }, logical(1L)))
        if (!length(evaluable)) {
          termination <- if (is.null(cap_reason)) {
            "no_evaluable_neighbor"
          } else cap_reason
          break
        }
        overturn_records <- evaluable[vapply(records[evaluable], function(record) {
          isTRUE(record$fit$overturn)
        }, logical(1L))]
        if (length(overturn_records)) {
          costs <- vapply(records[overturn_records], function(record) {
            nrow(.ab_change_rows(data, move_record_assignment(record, current)))
          }, integer(1L))
          ids <- vapply(records[overturn_records], function(record) {
            .ab_assignment_id(data, move_record_assignment(record, current))
          }, character(1L))
          selected_order <- order(costs, ids, method = "radix")
          best <- records[[overturn_records[selected_order[[1L]]]]]
          best$assignment <- current
          best$assignment[[best$move$cell_id]] <- best$move$to_label
          witness_record <- best
          log_rows[[best$log_index]]$selected <- TRUE
          step <- step + 1L
          path_rows[[length(path_rows) + 1L]] <- data.frame(
            step = step, assignment_id = .ab_assignment_id(data, best$assignment),
            cost = nrow(.ab_change_rows(data, best$assignment)),
            beta = best$fit$beta, margin = best$fit$margin,
            stringsAsFactors = FALSE
          )
          result_status <- "verified_counterexample"
          termination <- "verified_counterexample_found"
          break
        }
        margins <- vapply(records[evaluable], function(record) record$fit$margin,
                          numeric(1L))
        ids <- vapply(records[evaluable], function(record) {
          .ab_assignment_id(data, move_record_assignment(record, current))
        }, character(1L))
        best_order <- order(margins, ids, method = "radix")
        best <- records[[evaluable[best_order[[1L]]]]]
        best$assignment <- current
        best$assignment[[best$move$cell_id]] <- best$move$to_label
        if (!(best$fit$margin < current_fit$margin - data$claim$tolerance)) {
          termination <- if (!is.null(cap_reason)) {
            cap_reason
          } else if (neighborhood_truncated) {
            "neighborhood_candidate_cap_reached"
          } else "no_improving_move"
          break
        }
        current <- best$assignment
        current_fit <- best$fit
        current_id <- .ab_assignment_id(data, current)
        step <- step + 1L
        log_rows[[best$log_index]]$selected <- TRUE
        path_rows[[length(path_rows) + 1L]] <- data.frame(
          step = step, assignment_id = current_id,
          cost = nrow(.ab_change_rows(data, current)),
          beta = current_fit$beta, margin = current_fit$margin,
          stringsAsFactors = FALSE
        )
        if (!is.null(cap_reason)) {
          termination <- cap_reason
          break
        }
      }
    } else if (identical(method, "random_block")) {
      step <- 0L
      repeat {
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        ordered_source <- .ab_randomized_single_source(data, current, direction)
        source_candidate_count <- max(source_candidate_count, length(ordered_source))
        if (!length(ordered_source)) {
          termination <- "no_eligible_neighbors"
          break
        }
        block_start <- 1L
        advanced <- FALSE
        stop_route <- FALSE
        while (block_start <= length(ordered_source)) {
          if (state$proposal_count >= max_proposals) {
            termination <- "max_proposals_reached"
            stop_route <- TRUE
            break
          }
          if (elapsed() >= max_seconds) {
            termination <- "max_seconds_reached"
            stop_route <- TRUE
            break
          }
          block_end <- min(
            length(ordered_source), block_start + candidate_block_size - 1L
          )
          records <- list()
          cap_reason <- NULL
          for (candidate_index in seq.int(block_start, block_end)) {
            if (state$proposal_count >= max_proposals) {
              cap_reason <- "max_proposals_reached"
              break
            }
            if (elapsed() >= max_seconds) {
              cap_reason <- "max_seconds_reached"
              break
            }
            candidate <- ordered_source[[candidate_index]]
            proposal <- current
            proposal[[candidate$cell_id]] <- candidate$to_label
            record <- evaluate_assignment(
              proposal, "randomized_order_block", current_id,
              source_rank = candidate$source_rank
            )
            record$move <- candidate
            if (is_cost_one_early_stop(record, proposal)) {
              cost_one_early_stop_record <- record
              records[[length(records) + 1L]] <- record
              break
            }
            record$assignment <- NULL
            records[[length(records) + 1L]] <- record
            if (isTRUE(record$cap)) {
              cap_reason <- "max_fits_reached"
              break
            }
          }
          if (!is.null(cost_one_early_stop_record)) {
            best <- cost_one_early_stop_record
            witness_record <- best
            log_rows[[best$log_index]]$selected <- TRUE
            step <- step + 1L
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step,
              assignment_id = .ab_assignment_id(data, best$assignment),
              cost = 1L,
              beta = best$fit$beta, margin = best$fit$margin,
              stringsAsFactors = FALSE
            )
            result_status <- "verified_counterexample"
            termination <- "cost_one_early_stop"
            stop_route <- TRUE
            break
          }
          evaluable <- which(vapply(records, function(record) {
            !is.null(record$fit) && identical(record$fit$status, "evaluable")
          }, logical(1L)))
          overturn_records <- evaluable[vapply(records[evaluable], function(record) {
            isTRUE(record$fit$overturn)
          }, logical(1L))]
          if (length(overturn_records)) {
            costs <- vapply(records[overturn_records], function(record) {
              nrow(.ab_change_rows(data, move_record_assignment(record, current)))
            }, integer(1L))
            cells <- vapply(records[overturn_records], function(record) {
              record$move$cell_id
            }, character(1L))
            labels <- vapply(records[overturn_records], function(record) {
              record$move$to_label
            }, character(1L))
            selected_order <- order(costs, cells, labels, method = "radix")
            best <- records[[overturn_records[selected_order[[1L]]]]]
            best$assignment <- move_record_assignment(best, current)
            witness_record <- best
            log_rows[[best$log_index]]$selected <- TRUE
            step <- step + 1L
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step,
              assignment_id = .ab_assignment_id(data, best$assignment),
              cost = nrow(.ab_change_rows(data, best$assignment)),
              beta = best$fit$beta, margin = best$fit$margin,
              stringsAsFactors = FALSE
            )
            result_status <- "verified_counterexample"
            termination <- "verified_counterexample_found"
            stop_route <- TRUE
            break
          }
          improving <- evaluable[vapply(records[evaluable], function(record) {
            record$fit$margin < current_fit$margin - data$claim$tolerance
          }, logical(1L))]
          if (length(improving)) {
            margins <- vapply(records[improving], function(record) {
              record$fit$margin
            }, numeric(1L))
            cells <- vapply(records[improving], function(record) {
              record$move$cell_id
            }, character(1L))
            labels <- vapply(records[improving], function(record) {
              record$move$to_label
            }, character(1L))
            best_order <- order(margins, cells, labels, method = "radix")
            best <- records[[improving[best_order[[1L]]]]]
            best$assignment <- move_record_assignment(best, current)
            current <- best$assignment
            current_fit <- best$fit
            current_id <- .ab_assignment_id(data, current)
            step <- step + 1L
            log_rows[[best$log_index]]$selected <- TRUE
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step, assignment_id = current_id,
              cost = nrow(.ab_change_rows(data, current)),
              beta = current_fit$beta, margin = current_fit$margin,
              stringsAsFactors = FALSE
            )
            advanced <- TRUE
            if (!is.null(cap_reason)) {
              termination <- cap_reason
              stop_route <- TRUE
            }
            break
          }
          if (!is.null(cap_reason)) {
            termination <- cap_reason
            stop_route <- TRUE
            break
          }
          block_start <- block_end + 1L
        }
        if (stop_route) break
        if (advanced) next
        termination <- "no_improving_move_after_randomized_block_scan"
        break
      }
    } else if (identical(method, "ranked")) {
      step <- 0L
      switch_phase("ranking")
      proxy_cache_started <- unname(proc.time()[["elapsed"]])
      proxy_cache <- .ab_proxy_cell_cache(data)
      state$ranking_setup_seconds <- state$ranking_setup_seconds +
        unname(proc.time()[["elapsed"]] - proxy_cache_started)
      state$ranking_peak_bytes <- max(
        state$ranking_peak_bytes, as.numeric(utils::object.size(proxy_cache))
      )
      switch_phase("discovery")
      repeat {
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        switch_phase("ranking")
        ranking_started <- unname(proc.time()[["elapsed"]])
        ranked_source <- .ab_proxy_ranked_source(
          data, current, direction, proxy_cache,
          target_pseudocount = proxy_target_pseudocount,
          library_pseudocount = proxy_library_pseudocount,
          scale = proxy_scale
        )
        ranking_elapsed <- unname(proc.time()[["elapsed"]] - ranking_started)
        state$ranking_setup_count <- state$ranking_setup_count + 1L
        state$ranking_setup_seconds <- state$ranking_setup_seconds + ranking_elapsed
        state$ranking_peak_bytes <- max(
          state$ranking_peak_bytes,
          as.numeric(utils::object.size(proxy_cache)) +
            as.numeric(utils::object.size(ranked_source))
        )
        switch_phase("discovery")
        source_candidate_count <- max(source_candidate_count, length(ranked_source))
        if (!length(ranked_source)) {
          termination <- "no_eligible_neighbors"
          break
        }
        block_start <- 1L
        advanced <- FALSE
        stop_route <- FALSE
        while (block_start <= length(ranked_source)) {
          if (state$proposal_count >= max_proposals) {
            termination <- "max_proposals_reached"
            stop_route <- TRUE
            break
          }
          if (elapsed() >= max_seconds) {
            termination <- "max_seconds_reached"
            stop_route <- TRUE
            break
          }
          block_end <- min(
            length(ranked_source), block_start + candidate_block_size - 1L
          )
          records <- list()
          cap_reason <- NULL
          for (candidate_index in seq.int(block_start, block_end)) {
            if (state$proposal_count >= max_proposals) {
              cap_reason <- "max_proposals_reached"
              break
            }
            if (elapsed() >= max_seconds) {
              cap_reason <- "max_seconds_reached"
              break
            }
            candidate <- ranked_source[[candidate_index]]
            proposal <- current
            proposal[[candidate$cell_id]] <- candidate$to_label
            record <- evaluate_assignment(
              proposal, "ranked_proxy_block", current_id,
              source_rank = candidate$source_rank,
              proxy_score = candidate$proxy_score
            )
            record$move <- candidate
            if (is_cost_one_early_stop(record, proposal)) {
              cost_one_early_stop_record <- record
              records[[length(records) + 1L]] <- record
              break
            }
            record$assignment <- NULL
            records[[length(records) + 1L]] <- record
            if (isTRUE(record$cap)) {
              cap_reason <- "max_fits_reached"
              break
            }
          }
          if (!is.null(cost_one_early_stop_record)) {
            best <- cost_one_early_stop_record
            witness_record <- best
            log_rows[[best$log_index]]$selected <- TRUE
            step <- step + 1L
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step,
              assignment_id = .ab_assignment_id(data, best$assignment),
              cost = 1L,
              beta = best$fit$beta, margin = best$fit$margin,
              stringsAsFactors = FALSE
            )
            result_status <- "verified_counterexample"
            termination <- "cost_one_early_stop"
            stop_route <- TRUE
            break
          }
          evaluable <- which(vapply(records, function(record) {
            !is.null(record$fit) && identical(record$fit$status, "evaluable")
          }, logical(1L)))
          overturn_records <- evaluable[vapply(records[evaluable], function(record) {
            isTRUE(record$fit$overturn)
          }, logical(1L))]
          if (length(overturn_records)) {
            costs <- vapply(records[overturn_records], function(record) {
              nrow(.ab_change_rows(data, move_record_assignment(record, current)))
            }, integer(1L))
            cells <- vapply(records[overturn_records], function(record) {
              record$move$cell_id
            }, character(1L))
            labels <- vapply(records[overturn_records], function(record) {
              record$move$to_label
            }, character(1L))
            selected_order <- order(costs, cells, labels, method = "radix")
            best <- records[[overturn_records[selected_order[[1L]]]]]
            best$assignment <- current
            best$assignment[[best$move$cell_id]] <- best$move$to_label
            witness_record <- best
            log_rows[[best$log_index]]$selected <- TRUE
            step <- step + 1L
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step,
              assignment_id = .ab_assignment_id(data, best$assignment),
              cost = nrow(.ab_change_rows(data, best$assignment)),
              beta = best$fit$beta, margin = best$fit$margin,
              stringsAsFactors = FALSE
            )
            result_status <- "verified_counterexample"
            termination <- "verified_counterexample_found"
            stop_route <- TRUE
            break
          }
          improving <- evaluable[vapply(records[evaluable], function(record) {
            record$fit$margin < current_fit$margin - data$claim$tolerance
          }, logical(1L))]
          if (length(improving)) {
            margins <- vapply(records[improving], function(record) {
              record$fit$margin
            }, numeric(1L))
            cells <- vapply(records[improving], function(record) {
              record$move$cell_id
            }, character(1L))
            labels <- vapply(records[improving], function(record) {
              record$move$to_label
            }, character(1L))
            best_order <- order(margins, cells, labels, method = "radix")
            best <- records[[improving[best_order[[1L]]]]]
            best$assignment <- current
            best$assignment[[best$move$cell_id]] <- best$move$to_label
            current <- best$assignment
            current_fit <- best$fit
            current_id <- .ab_assignment_id(data, current)
            step <- step + 1L
            log_rows[[best$log_index]]$selected <- TRUE
            path_rows[[length(path_rows) + 1L]] <- data.frame(
              step = step, assignment_id = current_id,
              cost = nrow(.ab_change_rows(data, current)),
              beta = current_fit$beta, margin = current_fit$margin,
              stringsAsFactors = FALSE
            )
            advanced <- TRUE
            if (!is.null(cap_reason)) {
              termination <- cap_reason
              stop_route <- TRUE
            }
            break
          }
          if (!is.null(cap_reason)) {
            termination <- cap_reason
            stop_route <- TRUE
            break
          }
          block_start <- block_end + 1L
        }
        if (stop_route) break
        if (advanced) next
        termination <- "no_improving_move_after_ranked_scan"
        break
      }
    } else if (identical(method, "random_prune")) {
      discovery_fit_limit <- max(
        1, max_fits - effective_pruning_fit_reserve
      )
      random_index <- 1L
      if (!is.null(normalized_warm_start)) {
        record <- evaluate_assignment(
          normalized_warm_start, "warm_start_discovery", current_id
        )
        if (isTRUE(record$cap)) {
          termination <- "max_fits_reached_before_warm_start"
        } else if (length(record$failures)) {
          termination <- "warm_start_infeasible"
        } else if (is.null(record$fit) ||
                   !identical(record$fit$status, "evaluable")) {
          termination <- "warm_start_not_assessable"
        } else if (!isTRUE(record$fit$overturn)) {
          termination <- "warm_start_does_not_overturn"
        } else {
          witness_record <- record
          log_rows[[record$log_index]]$selected <- TRUE
          step <- 1L
          path_rows[[length(path_rows) + 1L]] <- data.frame(
            step = step,
            assignment_id = .ab_assignment_id(data, record$assignment),
            cost = nrow(.ab_change_rows(data, record$assignment)),
            beta = record$fit$beta, margin = record$fit$margin,
            stringsAsFactors = FALSE
          )
          result_status <- "verified_counterexample"
        }
      } else repeat {
        if (state$proposal_count >= max_proposals) {
          termination <- "max_proposals_reached"
          break
        }
        if (elapsed() >= max_seconds) {
          termination <- "max_seconds_reached"
          break
        }
        if (state$attempted >= discovery_fit_limit) {
          termination <- "discovery_fit_reserve_reached_no_witness"
          break
        }
        if (!is.null(random_domain)) {
          if (random_index > length(random_domain)) {
            termination <- if (length(random_domain)) {
              "random_allowed_domain_exhausted"
            } else "no_eligible_changes"
            break
          }
          proposal <- random_domain[[random_index]]
          random_index <- random_index + 1L
        } else {
          proposal <- .ab_random_assignment(data, direction)
          if (is.null(proposal)) {
            termination <- "no_eligible_changes"
            break
          }
        }
        record <- evaluate_assignment(proposal, "random_discovery", current_id)
        if (isTRUE(record$cap)) {
          termination <- "max_fits_reached"
          break
        }
        if (!is.null(record$fit) && identical(record$fit$status, "evaluable") &&
            isTRUE(record$fit$overturn)) {
          witness_record <- record
          log_rows[[record$log_index]]$selected <- TRUE
          step <- 1L
          path_rows[[length(path_rows) + 1L]] <- data.frame(
            step = step,
            assignment_id = .ab_assignment_id(data, record$assignment),
            cost = nrow(.ab_change_rows(data, record$assignment)),
            beta = record$fit$beta, margin = record$fit$margin,
            stringsAsFactors = FALSE
          )
          result_status <- "verified_counterexample"
          break
        }
      }
      if (!is.null(witness_record)) {
        switch_phase("pruning")
        pruning_complete <- TRUE
        pruning_stop_reason <- NULL
        can_prune <- function() {
          if (state$proposal_count >= max_proposals) {
            pruning_stop_reason <<- "max_proposals_reached_during_pruning"
            return(FALSE)
          }
          if (elapsed() >= max_seconds) {
            pruning_stop_reason <<- "max_seconds_reached_during_pruning"
            return(FALSE)
          }
          if (state$attempted >= max_fits) {
            pruning_stop_reason <<- "max_fits_reached_during_pruning"
            return(FALSE)
          }
          TRUE
        }
        try_removal <- function(cell_ids, proposal_type) {
          if (!can_prune()) return(NULL)
          proposal <- witness_record$assignment
          proposal[cell_ids] <- data$original_assignment[cell_ids]
          record <- evaluate_assignment(
            proposal, proposal_type,
            .ab_assignment_id(data, witness_record$assignment)
          )
          accepted <- !isTRUE(record$cap) && !length(record$failures) &&
            !is.null(record$fit) && identical(record$fit$status, "evaluable") &&
            isTRUE(record$fit$overturn)
          pruning_log[[length(pruning_log) + 1L]] <<- data.frame(
            attempt = as.integer(length(pruning_log) + 1L),
            proposal_type = proposal_type,
            removed_cells = paste(cell_ids, collapse = ";"),
            assignment_id = .ab_assignment_id(data, proposal),
            cost = as.integer(nrow(.ab_change_rows(data, proposal))),
            feasible = !length(record$failures),
            fit_attempted = isTRUE(log_rows[[record$log_index]]$fit_attempted),
            evaluable = !is.null(record$fit) &&
              identical(record$fit$status, "evaluable"),
            overturn = !is.null(record$fit) && isTRUE(record$fit$overturn),
            accepted = accepted,
            attempted_fit_count = as.integer(state$attempted),
            elapsed_seconds = elapsed(),
            stringsAsFactors = FALSE
          )
          if (accepted) {
            witness_record <<- record
            log_rows[[record$log_index]]$selected <- TRUE
            step <<- step + 1L
            path_rows[[length(path_rows) + 1L]] <<- data.frame(
              step = step,
              assignment_id = .ab_assignment_id(data, record$assignment),
              cost = nrow(.ab_change_rows(data, record$assignment)),
              beta = record$fit$beta, margin = record$fit$margin,
              stringsAsFactors = FALSE
            )
            record_witness_event(record, "pruning_improvement")
          }
          record
        }

        changed <- sort(
          names(witness_record$assignment)[
            witness_record$assignment != data$original_assignment
          ], method = "radix"
        )
        if (length(changed)) {
          groups <- split(
            changed,
            ceiling(seq_along(changed) / effective_pruning_chunk_size)
          )
          for (group in groups) {
            if (!can_prune()) {
              pruning_complete <- FALSE
              break
            }
            try_removal(group, "pruning_chunk_removal")
          }
        }
        if (pruning_complete) {
          changed <- sort(
            names(witness_record$assignment)[
              witness_record$assignment != data$original_assignment
            ], method = "radix"
          )
          for (cell_id in changed) {
            if (!can_prune()) {
              pruning_complete <- FALSE
              break
            }
            try_removal(cell_id, "pruning_single_removal")
          }
        }
        termination <- if (isTRUE(pruning_complete)) {
          "verified_counterexample_pruning_complete"
        } else {
          pruning_stop_reason
        }
      }
    }
  }

  close_phase()
  candidate_log <- if (length(log_rows)) do.call(rbind, log_rows) else .ab_empty_search_log()
  rownames(candidate_log) <- NULL
  result <- .ab_search_skeleton(result_status, termination, hashes)
  result$method <- method
  result$seed <- as.integer(seed)
  result$direction <- direction
  result$max_fits <- if (is.infinite(max_fits)) Inf else as.integer(max_fits)
  result$max_proposals <- max_proposals
  result$max_seconds <- max_seconds
  result$candidate_block_size <- candidate_block_size
  result$pruning_fit_reserve <- effective_pruning_fit_reserve
  result$pruning_chunk_size <- effective_pruning_chunk_size
  result$proxy_target_pseudocount <- as.numeric(proxy_target_pseudocount)
  result$proxy_library_pseudocount <- as.numeric(proxy_library_pseudocount)
  result$proxy_scale <- as.numeric(proxy_scale)
  result$original_fit <- original_fit
  result$proposal_count <- state$proposal_count
  result$duplicate_proposal_count <- state$duplicate_count
  result$cache_hit_count <- state$cache_hits
  result$unique_assignment_count <- length(ls(seen, all.names = TRUE))
  result$unique_proposal_count <- max(0L, result$unique_assignment_count - 1L)
  result$infeasible_proposal_count <- state$infeasible
  result$attempted_fit_count <- state$new_fit_count
  result$cumulative_fit_count <- state$attempted
  result$successful_fit_count <- state$successful
  result$failed_fit_count <- state$failed
  result$ranking_setup_count <- state$ranking_setup_count
  result$ranking_setup_seconds <- state$ranking_setup_seconds
  result$ranking_peak_bytes <- state$ranking_peak_bytes
  result$candidate_log <- candidate_log
  result$candidates <- candidate_log
  result$assignments <- candidate_log
  result$search_path <- if (length(path_rows)) do.call(rbind, path_rows) else data.frame()
  result$elapsed_seconds <- elapsed()
  result$component_accounting <- do.call(rbind, lapply(
    names(state$phase_seconds), function(phase) {
      rows <- candidate_log$search_phase == phase
      data.frame(
        phase = phase,
        proposal_count = sum(rows & candidate_log$proposal_index > 0L),
        attempted_fit_count = sum(rows & candidate_log$fit_attempted),
        resumed_fit_count = sum(rows & candidate_log$fit_resumed),
        cumulative_fit_count = sum(
          rows & (candidate_log$fit_attempted | candidate_log$fit_resumed)
        ),
        new_fit_count = sum(
          rows & candidate_log$fit_attempted & !candidate_log$fit_resumed
        ),
        successful_fit_count = sum(
          rows & candidate_log$fit_attempted &
            candidate_log$fit_status == "evaluable", na.rm = TRUE
        ),
        failed_fit_count = sum(
          rows & candidate_log$fit_attempted &
            candidate_log$fit_status != "evaluable", na.rm = TRUE
        ),
        elapsed_seconds = as.numeric(state$phase_seconds[[phase]]),
        stringsAsFactors = FALSE
      )
    }
  ))
  unique_log <- candidate_log[!duplicated(candidate_log$assignment_id), , drop = FALSE]
  result$total_assignment_count <- nrow(unique_log)
  result$feasible_assignment_count <- sum(unique_log$feasible)
  result$infeasible_assignment_count <- sum(!unique_log$feasible)
  result$evaluable_assignment_count <- sum(
    unique_log$fit_status == "evaluable", na.rm = TRUE
  )
  result$unevaluable_assignment_count <- sum(
    (unique_log$fit_attempted | unique_log$fit_resumed) &
      unique_log$fit_status == "not_assessable",
    na.rm = TRUE
  )
  result$unevaluated_assignment_count <- sum(
    unique_log$feasible & !unique_log$fit_attempted &
      !unique_log$fit_resumed & !unique_log$cache_hit
  )
  result$overturn_assignment_count <- sum(unique_log$overturn)
  result$infeasible_proposal_rate <- if (result$proposal_count) {
    result$infeasible_proposal_count / result$proposal_count
  } else 0
  result$search_complete <- FALSE
  result$fit_accounting <- data.frame(
    proposals = result$proposal_count,
    duplicate_proposals = result$duplicate_proposal_count,
    cache_hits = result$cache_hit_count,
    unique_assignments = result$unique_assignment_count,
    infeasible_proposals = result$infeasible_proposal_count,
    infeasible_proposal_rate = result$infeasible_proposal_rate,
    attempted_fits = result$attempted_fit_count,
    successful_fits = result$successful_fit_count,
    failed_fits = result$failed_fit_count,
    resumed_fits = state$resumed_fit_count,
    new_fits = state$new_fit_count,
    cumulative_route_fits = state$attempted,
    cumulative_successful_fits = state$cumulative_successful,
    cumulative_failed_fits = state$cumulative_failed,
    resume_historical_fit_seconds = state$resume_historical_fit_seconds,
    resume_historical_search_elapsed_seconds =
      state$resume_historical_search_elapsed_seconds,
    original_fits = result$component_accounting$attempted_fit_count[
      result$component_accounting$phase == "original"
    ],
    discovery_fits = result$component_accounting$attempted_fit_count[
      result$component_accounting$phase == "discovery"
    ],
    ranking_fits = result$component_accounting$attempted_fit_count[
      result$component_accounting$phase == "ranking"
    ],
    pruning_fits = result$component_accounting$attempted_fit_count[
      result$component_accounting$phase == "pruning"
    ],
    verification_fits = 0L,
    stringsAsFactors = FALSE
  )
  if (!is.null(witness_record)) {
    assignment <- witness_record$assignment
    changes <- .ab_changed_with_provenance(data, assignment)
    witness_id <- .ab_assignment_id(data, assignment)
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
      beta = witness_record$fit$beta,
      margin = witness_record$fit$margin,
      overturn = TRUE,
      input_hashes = hashes,
      analysis_id = data$claim$estimator_id
    )
    snapshot <- function(record) {
      row <- candidate_log[record$log_index, , drop = FALSE]
      list(
        assignment_id = .ab_assignment_id(data, record$assignment),
        assignment = record$assignment,
        cost = as.integer(nrow(.ab_change_rows(data, record$assignment))),
        beta = as.numeric(record$fit$beta),
        margin = as.numeric(record$fit$margin),
        fit_index = as.integer(row$fit_rank[[1L]]),
        attempted_fit_count = as.integer(row$fit_rank[[1L]] + 1L),
        proposal_index = as.integer(row$proposal_index[[1L]]),
        elapsed_seconds = if (!is.null(record$fit$search_elapsed_seconds) &&
                              is.finite(record$fit$search_elapsed_seconds)) {
          as.numeric(record$fit$search_elapsed_seconds)
        } else as.numeric(row$elapsed_seconds[[1L]])
      )
    }
    result$first_witness <- snapshot(first_witness_record)
    result$final_witness <- snapshot(witness_record)
    result$fits_to_first_witness <- result$first_witness$attempted_fit_count
    result$fits_to_final_witness <- result$final_witness$attempted_fit_count
    result$time_to_first_full_fit_witness_seconds <-
      result$first_witness$elapsed_seconds
    result$time_to_final_full_fit_witness_seconds <-
      result$final_witness$elapsed_seconds
    final_id <- result$final_witness$assignment_id
    if (!length(witness_events) ||
        !identical(witness_events[[length(witness_events)]]$assignment_id[[1L]],
                   final_id)) {
      final <- result$final_witness
      witness_events[[length(witness_events) + 1L]] <- data.frame(
        event_index = as.integer(length(witness_events) + 1L),
        event = "selected_final_witness",
        assignment_id = final$assignment_id, cost = final$cost,
        beta = final$beta, margin = final$margin,
        cumulative_best_cost = as.integer(min(c(
          vapply(witness_events, function(row) row$cost[[1L]], integer(1L)),
          final$cost
        ))),
        attempted_fit_count = final$attempted_fit_count,
        proposal_count = final$proposal_index,
        elapsed_seconds = final$elapsed_seconds, is_final = TRUE,
        stringsAsFactors = FALSE
      )
    } else {
      witness_events[[length(witness_events)]]$is_final <- TRUE
    }
  }
  result$witness_trajectory <- if (length(witness_events)) {
    do.call(rbind, witness_events)
  } else data.frame()
  result$pruning_complete <- pruning_complete
  result$search_censored <- grepl(
    "max_seconds|max_fits|max_proposals|reserve_reached", termination
  )
  result$resumed_fit_count <- state$resumed_fit_count
  result$new_fit_count <- state$new_fit_count
  result$cumulative_fit_count <- state$attempted
  result$resume_historical_fit_seconds <- state$resume_historical_fit_seconds
  result$resume_historical_search_elapsed_seconds <-
    state$resume_historical_search_elapsed_seconds
  result$progress_errors <- unique(state$progress_errors)
  cost_one_early_stop_triggered <- !is.null(cost_one_early_stop_record)
  result$cost_one_early_stop_triggered <- cost_one_early_stop_triggered
  result$cost_one_early_stop_replay_required <- cost_one_early_stop_triggered
  cost_one_early_stop <- list(
    eligible_route = method %in% c("random_block", "ranked"),
    triggered = cost_one_early_stop_triggered,
    assignment_id = if (cost_one_early_stop_triggered) {
      .ab_assignment_id(data, cost_one_early_stop_record$assignment)
    } else NA_character_,
    cost = if (cost_one_early_stop_triggered) 1L else NA_integer_,
    fit_index = if (cost_one_early_stop_triggered) {
      as.integer(candidate_log$fit_rank[
        cost_one_early_stop_record$log_index
      ])
    } else NA_integer_,
    proposal_index = if (cost_one_early_stop_triggered) {
      as.integer(candidate_log$proposal_index[
        cost_one_early_stop_record$log_index
      ])
    } else NA_integer_,
    original_fit_evaluable = identical(original_fit$status, "evaluable"),
    original_claim_holds = identical(original_fit$status, "evaluable") &&
      isTRUE(original_fit$claim_holds),
    candidate_fit_evaluable = cost_one_early_stop_triggered &&
      identical(cost_one_early_stop_record$fit$status, "evaluable"),
    candidate_overturn = cost_one_early_stop_triggered &&
      isTRUE(cost_one_early_stop_record$fit$overturn),
    independent_replay_required = cost_one_early_stop_triggered,
    independent_replay_status = if (cost_one_early_stop_triggered) {
      "pending_external_replay"
    } else "not_applicable",
    replay_failure_action = if (cost_one_early_stop_triggered) {
      "reject_search_result"
    } else "not_applicable"
  )
  result$cost_one_early_stop <- cost_one_early_stop
  pruning_table <- if (length(pruning_log)) do.call(rbind, pruning_log) else data.frame()
  route_id <- switch(method,
    greedy = "canonical_greedy",
    greedy_random = "greedy_random",
    random_block = "randomized_order_block",
    random_prune = if (is.null(normalized_warm_start)) {
      "random_prune"
    } else "random_prune_warm_start",
    ranked = "ranked",
    method
  )
  candidate_order_metadata <- list(
    source_order = switch(method,
      greedy = "legacy_candidate_construction_prefix_then_assignment_id",
      pair = "legacy_candidate_construction_prefix_then_assignment_id",
      greedy_random = "cell_id_assigned_label_then_seeded_full_domain_permutation",
      random_block = "cell_id_assigned_label_then_seeded_full_domain_permutation_in_blocks",
      ranked = "proxy_score_then_cell_id_assigned_label",
      random_prune = if (is.null(normalized_warm_start)) {
        "seeded_random_feasible_assignment_discovery"
      } else "supplied_warm_start_assignment",
      "legacy_seeded_random_feasible_assignment_discovery"
    ),
    full_domain_ordered_before_truncation = method %in% c(
      "greedy_random", "random_block", "ranked"
    ),
    authoritative_selection_tie_break = if (method %in% c(
      "random_block", "ranked"
    )) {
      "cost_then_cell_id_assigned_label"
    } else "cost_then_assignment_id",
    source_rank_field = "source_rank",
    proposal_rank_field = "proposal_rank",
    fit_rank_field = "fit_rank"
  )
  result$efficiency <- list(
    route_id = route_id,
    api_method = method,
    source_candidate_count = as.integer(source_candidate_count),
    ranking_setup_count = as.integer(state$ranking_setup_count),
    ranking_setup_seconds = as.numeric(state$ranking_setup_seconds),
    ranking_peak_bytes = as.numeric(state$ranking_peak_bytes),
    original_fit_count = as.integer(result$fit_accounting$original_fits),
    discovery_fit_count = as.integer(result$fit_accounting$discovery_fits),
    ranking_fit_count = as.integer(result$fit_accounting$ranking_fits),
    pruning_fit_count = as.integer(result$fit_accounting$pruning_fits),
    resumed_fit_count = as.integer(state$resumed_fit_count),
    new_fit_count = as.integer(state$new_fit_count),
    cumulative_fit_count = as.integer(state$attempted),
    resume_historical_fit_seconds = as.numeric(
      state$resume_historical_fit_seconds
    ),
    resume_historical_search_elapsed_seconds = as.numeric(
      state$resume_historical_search_elapsed_seconds
    ),
    current_search_elapsed_seconds = as.numeric(result$elapsed_seconds),
    search_censored = isTRUE(result$search_censored),
    cost_one_early_stop = cost_one_early_stop,
    progress_errors = unique(state$progress_errors),
    resume_evaluations_provided = as.integer(length(resume_evaluations)),
    resume_evaluations_consumed = as.integer(state$resumed_fit_count),
    first_witness = result$first_witness,
    final_witness = result$final_witness,
    pruning_complete = pruning_complete,
    pruning_log = pruning_table,
    candidate_order_metadata = candidate_order_metadata,
    component_accounting = result$component_accounting
  )
  .ab_attach_contract(result, data, resolved$identity)
}
