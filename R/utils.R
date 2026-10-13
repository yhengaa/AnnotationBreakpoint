.ab_sort <- function(x) {
  x <- enc2utf8(as.character(x))
  x[order(x, method = "radix", na.last = TRUE)]
}

.ab_unique_sorted <- function(x) {
  .ab_sort(unique(as.character(x)))
}

.ab_compact <- function(x) {
  x[!is.na(x) & nzchar(x)]
}

.ab_as_character <- function(x) {
  if (is.factor(x)) as.character(x) else as.character(x)
}

.ab_table_names <- function(x) {
  if (is.data.frame(x) || is.list(x)) names(x) else character()
}

.ab_field <- function(x, name, default = NULL) {
  if (!name %in% .ab_table_names(x)) {
    return(default)
  }
  value <- x[[name]]
  if (is.data.frame(x)) {
    if (nrow(x) != 1L) {
      return(default)
    }
    if (is.list(value) && length(value) == 1L) {
      value <- value[[1L]]
    } else {
      value <- value[[1L]]
    }
  }
  value
}

.ab_scalar_character <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) {
    return(NA_character_)
  }
  as.character(x)
}

.ab_scalar_number <- function(x) {
  if (is.null(x) || length(x) != 1L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}

.ab_is_whole <- function(x) {
  length(x) == 1L && is.finite(x) && x >= 0 && x == floor(x)
}

.ab_hash <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 3L)
}

.ab_file_hash <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

.ab_pair_id <- function(left, right) {
  paste0(
    nchar(enc2utf8(left), type = "bytes"), ":", enc2utf8(left),
    nchar(enc2utf8(right), type = "bytes"), ":", enc2utf8(right)
  )
}

.ab_hash_payloads <- function(data) {
  counts_triplet <- Matrix::summary(data$counts)
  counts_payload <- list(
    dim = unname(dim(data$counts)),
    dimnames = unname(dimnames(data$counts)),
    i = counts_triplet$i,
    j = counts_triplet$j,
    x = counts_triplet$x
  )
  cells_payload <- unclass(data$cells)
  specimens_payload <- unclass(data$specimens)
  allowed_payload <- unclass(data$allowed_labels)
  claim_payload <- data$claim
  constraints_payload <- data$constraints

  component_hashes <- c(
    counts = .ab_hash(counts_payload),
    cells = .ab_hash(cells_payload),
    specimens = .ab_hash(specimens_payload),
    allowed_labels = .ab_hash(allowed_payload),
    claim = .ab_hash(claim_payload),
    constraints = .ab_hash(constraints_payload)
  )
  specification_payload <- if (is.null(data$analysis_spec)) {
    list(
      claim = claim_payload,
      estimator = list(
        id = "portable-log-proportion-contrast-v1",
        pseudocount_alpha = 0.5,
        scale = 1e6,
        sample_weighting = "equal"
      )
    )
  } else {
    list(
      claim = claim_payload,
      analysis_spec = data$analysis_spec,
      runtime_policy = .ab_normalize_runtime_policy(data$runtime_policy)
    )
  }
  hashes <- c(
    component_hashes,
    input_hash = .ab_hash(component_hashes[c(
      "counts", "cells", "specimens", "allowed_labels"
    )]),
    specification_hash = .ab_hash(specification_payload),
    constraint_hash = .ab_hash(constraints_payload)
  )
  if (!is.null(data$analysis_spec)) {
    hashes <- c(hashes, analysis_spec_hash = .ab_hash(data$analysis_spec))
  }
  hashes
}

.ab_recompute_input_hashes <- function(x) {
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    stop("x must be a successful ab_validation object", call. = FALSE)
  }
  .ab_hash_payloads(x$data)
}

.ab_assignment_id <- function(data, assignment) {
  assignment <- assignment[.ab_sort(names(assignment))]
  original <- data$original_assignment[names(assignment)]
  changed <- assignment != original
  changed_pairs <- data.frame(
    cell_id = names(assignment)[changed],
    assigned_label = unname(assignment[changed]),
    stringsAsFactors = FALSE
  )
  canonical_json <- jsonlite::toJSON(
    list(
      schema = "annotation-breakpoint-assignment-v1",
      input_hash = unname(data$input_hashes[["input_hash"]]),
      changed_pairs = changed_pairs
    ),
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    digits = NA,
    pretty = FALSE
  )
  digest::digest(canonical_json, algo = "sha256", serialize = FALSE)
}

.ab_change_rows <- function(data, assignment) {
  cells <- data$cells
  assigned <- unname(assignment[cells$cell_id])
  changed <- assigned != cells$original_label
  out <- cells[changed, c(
    "cell_id", "sample_id", "donor_id", "condition", "original_label"
  ), drop = FALSE]
  out$assigned_label <- assigned[changed]
  if (nrow(out)) {
    out <- out[order(out$cell_id, method = "radix"), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.ab_changes_text <- function(changes) {
  if (!nrow(changes)) {
    return("original")
  }
  paste0(changes$cell_id, "=", changes$assigned_label, collapse = ";")
}

.ab_changed_with_provenance <- function(data, assignment) {
  changes <- .ab_change_rows(data, assignment)
  if (!nrow(changes)) {
    changes$evidence_source <- character()
    changes$evidence_version <- character()
    if ("calibrated_confidence" %in% names(data$allowed_labels)) {
      changes$calibrated_confidence <- numeric()
    }
    return(changes)
  }
  selected <- data$allowed_labels[
    match(
      .ab_pair_id(changes$cell_id, changes$assigned_label),
      .ab_pair_id(
        data$allowed_labels$cell_id,
        data$allowed_labels$allowed_label
      )
    ),
    , drop = FALSE
  ]
  changes$evidence_source <- selected$evidence_source
  changes$evidence_version <- selected$evidence_version
  if ("calibrated_confidence" %in% names(data$allowed_labels)) {
    changes$calibrated_confidence <- selected$calibrated_confidence
  }
  changes
}

.ab_normalize_assignment <- function(data, assignment = NULL) {
  cell_ids <- data$cells$cell_id
  if (is.null(assignment)) {
    return(data$original_assignment)
  }
  if (is.factor(assignment)) {
    assignment <- as.character(assignment)
  }
  if (!is.atomic(assignment) || is.null(names(assignment))) {
    stop("assignment must be a named character vector", call. = FALSE)
  }
  if (anyDuplicated(names(assignment))) {
    stop("assignment contains duplicate cell IDs", call. = FALSE)
  }
  assignment_names <- names(assignment)
  assignment <- as.character(assignment)
  names(assignment) <- assignment_names
  if (length(assignment) != length(cell_ids) ||
      !setequal(names(assignment), cell_ids)) {
    stop("assignment must map every input cell ID exactly once", call. = FALSE)
  }
  assignment <- assignment[cell_ids]
  if (anyNA(assignment) || any(!nzchar(assignment))) {
    stop("assignment labels must be nonempty", call. = FALSE)
  }
  allowed_by_cell <- data$allowed_by_cell
  if (!identical(names(allowed_by_cell), cell_ids)) {
    allowed_by_cell <- allowed_by_cell[cell_ids]
  }
  permitted <- vapply(
    seq_along(cell_ids),
    function(index) assignment[[index]] %in% allowed_by_cell[[index]],
    logical(1L)
  )
  if (!all(permitted)) {
    stop(
      paste0(
        "assignment contains impermissible labels for: ",
        paste(cell_ids[!permitted], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  assignment
}

.ab_population_summary <- function(data, assignment) {
  assignment <- .ab_normalize_assignment(data, assignment)
  rows <- vector("list", length(data$sample_ids) * length(data$registered_labels))
  at <- 1L
  for (sample_id in data$sample_ids) {
    sample_index <- data$cells$sample_id == sample_id
    for (label in data$registered_labels) {
      selected_ids <- data$cells$cell_id[sample_index & assignment == label]
      library <- if (length(selected_ids)) {
        sum(data$counts[data$claim$feature_universe, selected_ids, drop = FALSE])
      } else {
        0
      }
      rows[[at]] <- data.frame(
        sample_id = sample_id,
        label = label,
        cell_count = length(selected_ids),
        library_count = as.numeric(library),
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  do.call(rbind, rows)
}

.ab_constraint_failures <- function(data, assignment) {
  assignment <- .ab_normalize_assignment(data, assignment)
  changes <- .ab_change_rows(data, assignment)
  failures <- character()
  constraints <- data$constraints

  if (nrow(changes) > constraints$max_changes_global) {
    failures <- c(failures, "max_changes_global")
  }
  if (nrow(changes)) {
    per_sample <- table(factor(changes$sample_id, levels = data$sample_ids))
    over_sample <- names(per_sample)[
      as.numeric(per_sample) > constraints$max_changes_per_sample[names(per_sample)]
    ]
    if (length(over_sample)) {
      failures <- c(failures, paste0("max_changes_per_sample:", over_sample))
    }
    per_donor <- table(factor(changes$donor_id, levels = data$donor_ids))
    over_donor <- names(per_donor)[
      as.numeric(per_donor) > constraints$max_changes_per_donor[names(per_donor)]
    ]
    if (length(over_donor)) {
      failures <- c(failures, paste0("max_changes_per_donor:", over_donor))
    }
  }

  population <- .ab_population_summary(data, assignment)
  for (idx in seq_len(nrow(population))) {
    row <- population[idx, , drop = FALSE]
    label <- row$label[[1L]]
    if (row$cell_count[[1L]] < constraints$min_cells_per_label[[label]]) {
      failures <- c(failures, paste0("min_cells:", row$sample_id, ":", label))
    }
    if (row$library_count[[1L]] < constraints$min_library_per_label[[label]]) {
      failures <- c(failures, paste0("min_library:", row$sample_id, ":", label))
    }
  }

  focal <- data$claim$focal_label
  focal_valid <- population$label == focal &
    population$cell_count >= constraints$min_cells_per_label[[focal]] &
    population$library_count >= constraints$min_library_per_label[[focal]]
  valid_samples <- population$sample_id[focal_valid]
  valid_specimens <- data$specimens[
    data$specimens$sample_id %in% valid_samples,
    , drop = FALSE
  ]
  for (condition in c(data$claim$control_condition, data$claim$case_condition)) {
    donor_count <- length(unique(valid_specimens$donor_id[
      valid_specimens$condition == condition
    ]))
    if (donor_count < constraints$min_donors_per_condition) {
      failures <- c(failures, paste0("min_donors:", condition, ":", focal))
    }
  }
  .ab_sort(unique(failures))
}

.ab_long_pseudobulk <- function(data, pseudobulk) {
  triplet <- Matrix::summary(pseudobulk$counts)
  if (!nrow(triplet)) {
    return(data.frame(
      feature_id = character(), sample_id = character(), donor_id = character(),
      condition = character(), label = character(), count = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  columns <- pseudobulk$columns[triplet$j, , drop = FALSE]
  out <- data.frame(
    feature_id = rownames(pseudobulk$counts)[triplet$i],
    sample_id = columns$sample_id,
    donor_id = columns$donor_id,
    condition = columns$condition,
    label = columns$label,
    count = as.numeric(triplet$x),
    stringsAsFactors = FALSE
  )
  out <- out[order(
    out$sample_id, out$label, out$feature_id, method = "radix"
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.ab_empty_candidate_table <- function() {
  data.frame(
    rank = integer(),
    assignment_id = character(),
    changes = character(),
    cost = integer(),
    feasible = logical(),
    constraint_failures = character(),
    evaluated = logical(),
    evaluable = logical(),
    beta = numeric(),
    margin = numeric(),
    claim_holds = logical(),
    overturn = logical(),
    fit_elapsed_seconds = numeric(),
    warnings = character(),
    failure_reason = character(),
    evaluation_error = character(),
    stringsAsFactors = FALSE
  )
}

.ab_result_skeleton <- function(status, reason, input_hashes = character()) {
  structure(list(
    status = status,
    original_fit = NULL,
    assignments = .ab_empty_candidate_table(),
    candidates = .ab_empty_candidate_table(),
    effects_by_cost = data.frame(),
    total_assignment_count = 0L,
    feasible_assignment_count = 0L,
    infeasible_assignment_count = 0L,
    evaluable_assignment_count = 0L,
    unevaluable_assignment_count = 0L,
    unevaluated_assignment_count = 0L,
    overturn_assignment_count = 0L,
    best_verified_cost = NA_integer_,
    verified_upper_bound = NA_integer_,
    exact_minimum_cost = NA_integer_,
    witness = NULL,
    witness_id = NULL,
    minimum_cost_overturns = .ab_empty_candidate_table(),
    search_complete = FALSE,
    termination_reason = reason,
    input_hashes = input_hashes,
    errors = character(),
    max_evaluations = NA_real_,
    certified_lower_bound = NA_integer_,
    claim = NULL,
    constraints = NULL,
    analysis_spec = NULL,
    analysis_id = NULL,
    backend_identity = NULL,
    evaluation_identity = NULL,
    result_id = NULL
  ), class = "ab_result")
}

.ab_result_identity <- function(result) {
  candidates <- result$candidates
  if (is.data.frame(candidates)) {
    candidates <- candidates[, setdiff(names(candidates), "assignment"), drop = FALSE]
  }
  canonical <- jsonlite::toJSON(
    list(
      schema = "annotation-breakpoint-result-identity-v1",
      status = result$status,
      input_hashes = as.list(result$input_hashes),
      claim = result$claim,
      constraints = result$constraints,
      analysis_spec = result$analysis_spec,
      analysis_id = result$analysis_id,
      backend_identity = result$backend_identity,
      evaluation_identity = result$evaluation_identity,
      original = if (is.null(result$original_fit)) NULL else list(
        status = result$original_fit$status,
        beta = result$original_fit$beta,
        margin = result$original_fit$margin,
        claim_holds = result$original_fit$claim_holds,
        overturn = result$original_fit$overturn,
        error = result$original_fit$error,
        failure_reason = result$original_fit$failure_reason,
        warnings = result$original_fit$warnings,
        runtime = result$original_fit$runtime
      ),
      candidates = candidates,
      best_verified_cost = result$best_verified_cost,
      verified_upper_bound = result$verified_upper_bound,
      exact_minimum_cost = result$exact_minimum_cost,
      certified_lower_bound = result$certified_lower_bound,
      witness = result$witness,
      witness_id = result$witness_id,
      search_complete = result$search_complete,
      method = result$method,
      seed = result$seed,
      direction = result$direction,
      proxy_prefix_contract = if (is.null(result$proxy_prefix)) NULL else list(
        route = result$proxy_prefix$route,
        registered_prefix_sizes = result$proxy_prefix$registered_prefix_sizes,
        available_prefix_sizes = result$proxy_prefix$available_prefix_sizes,
        ranking_frozen_before_fits = result$proxy_prefix$ranking_frozen_before_fits,
        diagnostics_included = result$proxy_prefix$diagnostics_included
      ),
      fit_accounting = result$fit_accounting,
      search_path = result$search_path,
      termination_reason = result$termination_reason
    ),
    auto_unbox = TRUE,
    dataframe = "rows",
    null = "null",
    na = "null",
    digits = NA,
    pretty = FALSE
  )
  digest::digest(canonical, algo = "sha256", serialize = FALSE)
}

.ab_attach_contract <- function(result, data, evaluation_identity = NULL) {
  result$claim <- data$claim
  result$constraints <- data$constraints
  result$analysis_spec <- data$analysis_spec
  authoritative <- .ab_authoritative_identity(data)
  result$analysis_id <- authoritative$analysis_id
  result$backend_identity <- authoritative
  result$evaluation_identity <- if (is.null(evaluation_identity)) {
    authoritative
  } else {
    evaluation_identity
  }
  result$result_id <- .ab_result_identity(result)
  result
}
