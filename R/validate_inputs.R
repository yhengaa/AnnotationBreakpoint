#' Validate and normalize Annotation Breakpoint inputs
#'
#' @param counts A sparse feature-by-cell `Matrix::dgCMatrix`.
#' @param cells Cell metadata.
#' @param specimens Specimen metadata.
#' @param allowed_labels Cell-specific permitted labels and provenance.
#' @param claim A named one-row object describing the registered claim.
#' @param constraints A named one-row object describing budgets and minima.
#' @param analysis_spec Optional frozen authoritative edgeR analysis specification.
#' @param runtime_policy `"strict"` requires the registered R and edgeR versions;
#'   `"compatible"` permits patch-version differences only after required API checks.
#' @return An object of class `ab_validation`.
#' @export
validate_inputs <- function(counts, cells, specimens, allowed_labels,
                            claim, constraints, analysis_spec = NULL,
                            runtime_policy = c("strict", "compatible")) {
  errors <- character()
  warnings <- character()
  add_error <- function(message) {
    errors <<- c(errors, as.character(message))
  }
  normalized_runtime_policy <- tryCatch(
    .ab_normalize_runtime_policy(runtime_policy), error = function(e) e
  )
  if (inherits(normalized_runtime_policy, "error")) {
    errors <- c(errors, conditionMessage(normalized_runtime_policy))
    normalized_runtime_policy <- "strict"
  }
  normalized_analysis <- .ab_normalize_analysis_spec(analysis_spec)
  if (!normalized_analysis$ok) {
    errors <- c(errors, normalized_analysis$errors)
  }
  if (identical(normalized_runtime_policy, "compatible") &&
      is.null(analysis_spec)) {
    errors <- c(
      errors,
      "runtime_policy = 'compatible' requires an authoritative analysis_spec"
    )
  }

  required_cells <- c(
    "cell_id", "sample_id", "donor_id", "condition", "original_label"
  )
  required_specimens <- c("sample_id", "donor_id", "condition")
  required_allowed <- c(
    "cell_id", "allowed_label", "evidence_source", "evidence_version"
  )

  if (!inherits(counts, "dgCMatrix")) {
    add_error("counts must be a Matrix::dgCMatrix")
  }
  if (!is.data.frame(cells)) {
    add_error("cells must be a data frame")
  } else if (!all(required_cells %in% names(cells))) {
    add_error(paste0(
      "cells is missing required columns: ",
      paste(setdiff(required_cells, names(cells)), collapse = ", ")
    ))
  }
  if (!is.data.frame(specimens)) {
    add_error("specimens must be a data frame")
  } else if (!all(required_specimens %in% names(specimens))) {
    add_error(paste0(
      "specimens is missing required columns: ",
      paste(setdiff(required_specimens, names(specimens)), collapse = ", ")
    ))
  }
  if (!is.data.frame(allowed_labels)) {
    add_error("allowed_labels must be a data frame")
  } else if (!all(required_allowed %in% names(allowed_labels))) {
    add_error(paste0(
      "allowed_labels is missing required columns: ",
      paste(setdiff(required_allowed, names(allowed_labels)), collapse = ", ")
    ))
  }
  if (!(is.list(claim) || is.data.frame(claim))) {
    add_error("claim must be a named one-row list or data frame")
  } else if (is.data.frame(claim) && nrow(claim) != 1L) {
    add_error("claim must contain exactly one row")
  }
  if (!(is.list(constraints) || is.data.frame(constraints))) {
    add_error("constraints must be a named one-row list or data frame")
  } else if (is.data.frame(constraints) && nrow(constraints) != 1L) {
    add_error("constraints must contain exactly one row")
  }

  if (length(errors)) {
    return(structure(list(
      ok = FALSE,
      status = "invalid_input",
      errors = unique(errors),
      warnings = unique(warnings),
      input_hashes = character(),
      data = NULL
    ), class = "ab_validation"))
  }

  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  specimens <- as.data.frame(specimens, stringsAsFactors = FALSE)
  allowed_labels <- as.data.frame(allowed_labels, stringsAsFactors = FALSE)
  for (field in required_cells) cells[[field]] <- .ab_as_character(cells[[field]])
  for (field in required_specimens) {
    specimens[[field]] <- .ab_as_character(specimens[[field]])
  }
  for (field in required_allowed) {
    allowed_labels[[field]] <- .ab_as_character(allowed_labels[[field]])
  }

  check_nonempty <- function(table, fields, table_name) {
    for (field in fields) {
      bad <- is.na(table[[field]]) | !nzchar(table[[field]])
      if (any(bad)) {
        add_error(paste0(table_name, "$", field, " must contain nonempty values"))
      }
    }
  }
  check_nonempty(cells, required_cells, "cells")
  check_nonempty(specimens, required_specimens, "specimens")
  check_nonempty(allowed_labels, required_allowed, "allowed_labels")

  if (anyDuplicated(cells$cell_id)) add_error("cells$cell_id must be unique")
  if (anyDuplicated(specimens$sample_id)) {
    add_error("specimens$sample_id must be unique")
  }
  allowed_pairs <- paste0(
    nchar(allowed_labels$cell_id, type = "bytes"), ":",
    allowed_labels$cell_id,
    nchar(allowed_labels$allowed_label, type = "bytes"), ":",
    allowed_labels$allowed_label
  )
  if (anyDuplicated(allowed_pairs)) {
    add_error("allowed_labels contains duplicate (cell_id, allowed_label) rows")
  }
  if (any(grepl("[;=]", cells$cell_id))) {
    add_error("cell IDs may not contain the reserved delimiters ';' or '='")
  }
  if (any(grepl("::", specimens$sample_id, fixed = TRUE)) ||
      any(grepl("::", allowed_labels$allowed_label, fixed = TRUE))) {
    add_error("sample IDs and labels may not contain the reserved delimiter '::'")
  }
  if (any(grepl("[;=]", allowed_labels$allowed_label))) {
    add_error("labels may not contain the reserved delimiters ';' or '='")
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    add_error("counts must have feature row names and cell-ID column names")
  } else {
    if (anyNA(rownames(counts)) || any(!nzchar(rownames(counts)))) {
      add_error("counts feature IDs must be nonempty")
    }
    if (anyNA(colnames(counts)) || any(!nzchar(colnames(counts)))) {
      add_error("counts cell IDs must be nonempty")
    }
    if (anyDuplicated(rownames(counts))) add_error("counts feature IDs must be unique")
    if (anyDuplicated(colnames(counts))) add_error("counts cell IDs must be unique")
    if (!setequal(colnames(counts), cells$cell_id)) {
      add_error("counts columns must match cells$cell_id exactly")
    }
  }
  if (inherits(counts, "dgCMatrix")) {
    stored <- counts@x
    if (any(!is.finite(stored)) || any(stored < 0) ||
        any(stored != floor(stored))) {
      add_error("counts must be finite nonnegative integers")
    }
  }

  focal_gene <- .ab_scalar_character(
    .ab_field(claim, "focal_gene", .ab_field(claim, "target_gene_id"))
  )
  focal_label <- .ab_scalar_character(.ab_field(claim, "focal_label"))
  case_condition <- .ab_scalar_character(
    .ab_field(claim, "case_condition", .ab_field(claim, "numerator_condition"))
  )
  control_condition <- .ab_scalar_character(
    .ab_field(claim, "control_condition", .ab_field(claim, "denominator_condition"))
  )
  q <- .ab_scalar_number(.ab_field(claim, "q"))
  delta <- .ab_scalar_number(.ab_field(claim, "delta"))
  tolerance <- .ab_scalar_number(.ab_field(claim, "tolerance"))
  declared_estimator <- .ab_scalar_character(.ab_field(claim, "estimator_id"))
  feature_universe <- .ab_field(claim, "feature_universe")
  if (is.list(feature_universe) && length(feature_universe) == 1L) {
    feature_universe <- feature_universe[[1L]]
  }
  if (is.character(feature_universe) && length(feature_universe) == 1L) {
    feature_universe <- strsplit(feature_universe, "|", fixed = TRUE)[[1L]]
  }
  feature_universe <- as.character(feature_universe)

  if (is.na(focal_gene) || !nzchar(focal_gene)) add_error("claim focal_gene is required")
  if (is.na(focal_label) || !nzchar(focal_label)) {
    add_error("claim focal_label is required")
  }
  if (is.na(case_condition) || !nzchar(case_condition) ||
      is.na(control_condition) || !nzchar(control_condition) ||
      identical(case_condition, control_condition)) {
    add_error("claim case_condition and control_condition must be distinct and nonempty")
  }
  if (!length(q) || is.na(q) || !is.finite(q) || !q %in% c(-1, 1)) {
    add_error("claim q must be either -1 or +1")
  }
  if (!length(delta) || is.na(delta) || !is.finite(delta) || delta < 0) {
    add_error("claim delta must be finite and nonnegative")
  }
  if (!length(tolerance) || is.na(tolerance) || !is.finite(tolerance) ||
      tolerance < 0) {
    add_error("claim tolerance must be finite and nonnegative")
  }
  if (!length(feature_universe) || anyNA(feature_universe) ||
      any(!nzchar(feature_universe)) || anyDuplicated(feature_universe)) {
    add_error("claim feature_universe must contain unique nonempty feature IDs")
  } else {
    if (!focal_gene %in% feature_universe) {
      add_error("claim focal_gene must be in feature_universe")
    }
    if (!all(feature_universe %in% rownames(counts))) {
      add_error("every claim feature_universe ID must occur in counts")
    }
  }
  expected_estimator <- if (is.null(normalized_analysis$spec)) {
    "portable-log-proportion-contrast-v1"
  } else {
    normalized_analysis$spec$analysis_id
  }
  if (!is.null(normalized_analysis$spec)) {
    if (is.na(declared_estimator) || !identical(declared_estimator, expected_estimator)) {
      add_error("claim estimator_id must match analysis_spec$analysis_id")
    }
    if (!identical(control_condition,
                   normalized_analysis$spec$design$control_reference) ||
        !identical(case_condition, normalized_analysis$spec$design$case_level)) {
      add_error("claim conditions must match the frozen analysis design levels")
    }
  } else if (!is.na(declared_estimator) &&
             !identical(declared_estimator, expected_estimator)) {
    add_error("claim estimator_id does not match the portable Layer-A analysis")
  }

  constraint_names <- .ab_table_names(constraints)
  constrained_labels <- sub(
    "^min_cells_", "", grep("^min_cells_", constraint_names, value = TRUE)
  )
  constrained_labels <- constrained_labels[
    paste0("min_library_", constrained_labels) %in% constraint_names
  ]
  explicit_alternative <- .ab_scalar_character(
    .ab_field(claim, "alternative_label")
  )
  if (!is.na(explicit_alternative) && nzchar(explicit_alternative)) {
    registered_labels <- c(focal_label, explicit_alternative)
  } else if (length(constrained_labels) == 2L && focal_label %in% constrained_labels) {
    registered_labels <- c(focal_label, setdiff(constrained_labels, focal_label))
  } else {
    allowed_union <- unique(allowed_labels$allowed_label)
    movable_cells <- unique(allowed_labels$cell_id[
      ave(allowed_labels$allowed_label, allowed_labels$cell_id,
          FUN = length) > 1L
    ])
    movable_union <- unique(allowed_labels$allowed_label[
      allowed_labels$cell_id %in% movable_cells
    ])
    candidates <- if (length(movable_union)) movable_union else allowed_union
    if (length(candidates) == 2L && focal_label %in% candidates) {
      registered_labels <- c(focal_label, setdiff(candidates, focal_label))
    } else {
      registered_labels <- character()
      add_error(
        "the two registered labels cannot be derived; provide claim$alternative_label"
      )
    }
  }
  if (length(registered_labels) == 2L &&
      (anyNA(registered_labels) || any(!nzchar(registered_labels)) ||
       registered_labels[[1L]] == registered_labels[[2L]])) {
    add_error("the focal and alternative labels must be distinct and nonempty")
  }

  min_cells <- setNames(numeric(length(registered_labels)), registered_labels)
  min_library <- setNames(numeric(length(registered_labels)), registered_labels)
  if (length(registered_labels) == 2L) {
    for (label in registered_labels) {
      min_cells[[label]] <- .ab_scalar_number(
        .ab_field(constraints, paste0("min_cells_", label))
      )
      min_library[[label]] <- .ab_scalar_number(
        .ab_field(constraints, paste0("min_library_", label))
      )
      if (!.ab_is_whole(min_cells[[label]]) || min_cells[[label]] < 1) {
        add_error(paste0("minimum cells for label ", label, " must be a positive integer"))
      }
      if (!.ab_is_whole(min_library[[label]]) || min_library[[label]] < 1) {
        add_error(paste0(
          "minimum library count for label ", label, " must be a positive integer"
        ))
      }
    }
  }

  global_budget <- .ab_scalar_number(
    .ab_field(constraints, "max_changes_global")
  )
  min_donors <- .ab_scalar_number(
    .ab_field(constraints, "min_donors_per_condition")
  )
  if (!.ab_is_whole(global_budget)) {
    add_error("max_changes_global must be a nonnegative integer")
  }
  if (!.ab_is_whole(min_donors) || min_donors < 2) {
    add_error("min_donors_per_condition must be an integer of at least two")
  }

  sample_ids_raw <- unique(specimens$sample_id)
  donor_ids_raw <- unique(specimens$donor_id)
  parse_budget_map <- function(value, ids, field) {
    if (is.list(value) && length(value) == 1L) value <- value[[1L]]
    value <- suppressWarnings(as.numeric(value))
    original_names <- names(.ab_field(constraints, field))
    if (length(value) == 1L) {
      value <- rep(value, length(ids))
      names(value) <- ids
    } else if (!is.null(original_names) && length(original_names) == length(value)) {
      names(value) <- original_names
      if (anyDuplicated(names(value)) || !setequal(names(value), ids)) {
        add_error(paste0(field, " names must match the registered IDs exactly"))
        return(setNames(rep(NA_real_, length(ids)), ids))
      }
      value <- value[ids]
    } else {
      add_error(paste0(field, " must be a scalar or a named map over registered IDs"))
      return(setNames(rep(NA_real_, length(ids)), ids))
    }
    if (any(!is.finite(value)) || any(value < 0) || any(value != floor(value))) {
      add_error(paste0(field, " values must be nonnegative integers"))
    }
    value
  }
  sample_budgets <- parse_budget_map(
    .ab_field(constraints, "max_changes_per_sample"),
    sample_ids_raw,
    "max_changes_per_sample"
  )
  donor_budgets <- parse_budget_map(
    .ab_field(constraints, "max_changes_per_donor"),
    donor_ids_raw,
    "max_changes_per_donor"
  )

  if (anyDuplicated(specimens$donor_id)) {
    add_error("version 1 requires exactly one sample per donor")
  }
  if (!setequal(unique(cells$sample_id), specimens$sample_id)) {
    add_error("cells and specimens must contain the same sample IDs")
  }
  unknown_allowed <- setdiff(unique(allowed_labels$cell_id), cells$cell_id)
  missing_allowed <- setdiff(cells$cell_id, unique(allowed_labels$cell_id))
  if (length(unknown_allowed)) add_error("allowed_labels contains unknown cell IDs")
  if (length(missing_allowed)) add_error("every cell must have allowed-label rows")

  if (!anyDuplicated(specimens$sample_id) && nrow(cells)) {
    specimen_match <- match(cells$sample_id, specimens$sample_id)
    known <- !is.na(specimen_match)
    matched_donor <- specimens$donor_id[specimen_match]
    matched_condition <- specimens$condition[specimen_match]
    donor_mismatch <- known & (
      is.na(cells$donor_id) | is.na(matched_donor) |
        cells$donor_id != matched_donor
    )
    condition_mismatch <- known & (
      is.na(cells$condition) | is.na(matched_condition) |
        cells$condition != matched_condition
    )
    if (any(donor_mismatch | condition_mismatch, na.rm = TRUE)) {
      add_error("cell sample/donor/condition fields must match specimens")
    }
  }
  expected_conditions <- c(control_condition, case_condition)
  if (!setequal(unique(specimens$condition), expected_conditions)) {
    add_error("specimens must contain exactly the registered two conditions")
  }
  if (!anyDuplicated(specimens$donor_id) && .ab_is_whole(min_donors)) {
    donor_counts <- vapply(
      expected_conditions,
      function(condition) length(unique(specimens$donor_id[
        specimens$condition == condition
      ])),
      integer(1L)
    )
    if (any(donor_counts < min_donors)) {
      add_error("insufficient distinct donors in at least one condition")
    }
  }

  if (!length(errors) && length(registered_labels) == 2L) {
    allowed_by_cell_raw <- split(
      allowed_labels$allowed_label,
      factor(allowed_labels$cell_id, levels = cells$cell_id)
    )
    for (idx in seq_len(nrow(cells))) {
      cell_id <- cells$cell_id[[idx]]
      original <- cells$original_label[[idx]]
      permitted <- unique(allowed_by_cell_raw[[cell_id]])
      if (!original %in% permitted) {
        add_error(paste0("original label is not allowed for cell ", cell_id))
      }
      if (original %in% registered_labels) {
        if (any(!permitted %in% registered_labels)) {
          add_error(paste0(
            "registered-label cell ", cell_id,
            " has an allowed label outside the registered pair"
          ))
        }
      } else if (!identical(.ab_sort(permitted), original)) {
        add_error(paste0(
          "cells outside the registered pair must be immutable: ", cell_id
        ))
      }
    }
  }

  if ("calibrated_confidence" %in% names(allowed_labels)) {
    confidence <- suppressWarnings(as.numeric(allowed_labels$calibrated_confidence))
    supplied <- !is.na(allowed_labels$calibrated_confidence) &
      nzchar(as.character(allowed_labels$calibrated_confidence))
    if (any(supplied & (!is.finite(confidence) | confidence < 0 | confidence > 1))) {
      add_error("calibrated_confidence must be missing or between zero and one")
    }
    allowed_labels$calibrated_confidence <- confidence
  }

  if (length(errors)) {
    return(structure(list(
      ok = FALSE,
      status = "invalid_input",
      errors = unique(errors),
      warnings = unique(warnings),
      input_hashes = character(),
      data = NULL
    ), class = "ab_validation"))
  }

  cell_order <- order(cells$cell_id, method = "radix")
  cells <- cells[cell_order, , drop = FALSE]
  rownames(cells) <- NULL
  specimen_order <- order(specimens$sample_id, method = "radix")
  specimens <- specimens[specimen_order, , drop = FALSE]
  rownames(specimens) <- NULL
  allowed_order <- order(
    allowed_labels$cell_id, allowed_labels$allowed_label, method = "radix"
  )
  allowed_labels <- allowed_labels[allowed_order, , drop = FALSE]
  rownames(allowed_labels) <- NULL
  counts <- counts[
    order(rownames(counts), method = "radix"),
    cells$cell_id,
    drop = FALSE
  ]
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")

  sample_ids <- specimens$sample_id
  donor_ids <- .ab_sort(specimens$donor_id)
  sample_budgets <- sample_budgets[sample_ids]
  donor_budgets <- donor_budgets[donor_ids]
  allowed_by_cell <- split(
    allowed_labels$allowed_label,
    factor(allowed_labels$cell_id, levels = cells$cell_id)
  )
  allowed_by_cell <- lapply(allowed_by_cell, .ab_unique_sorted)
  names(allowed_by_cell) <- cells$cell_id
  original_assignment <- setNames(cells$original_label, cells$cell_id)

  normalized_claim <- list(
    focal_gene = focal_gene,
    focal_label = focal_label,
    alternative_label = registered_labels[[2L]],
    case_condition = case_condition,
    control_condition = control_condition,
    q = as.integer(q),
    delta = as.numeric(delta),
    tolerance = as.numeric(tolerance),
    feature_universe = as.character(feature_universe),
    pseudocount_alpha = 0.5,
    scale = 1e6,
    estimator_id = expected_estimator
  )
  normalized_constraints <- list(
    min_cells_per_label = as.numeric(min_cells),
    min_library_per_label = as.numeric(min_library),
    max_changes_global = as.integer(global_budget),
    max_changes_per_sample = as.numeric(sample_budgets),
    max_changes_per_donor = as.numeric(donor_budgets),
    min_donors_per_condition = as.integer(min_donors)
  )
  names(normalized_constraints$min_cells_per_label) <- names(min_cells)
  names(normalized_constraints$min_library_per_label) <- names(min_library)
  names(normalized_constraints$max_changes_per_sample) <- names(sample_budgets)
  names(normalized_constraints$max_changes_per_donor) <- names(donor_budgets)

  data <- list(
    counts = counts,
    cells = cells,
    specimens = specimens,
    allowed_labels = allowed_labels,
    claim = normalized_claim,
    constraints = normalized_constraints,
    registered_labels = registered_labels,
    sample_ids = sample_ids,
    donor_ids = donor_ids,
    allowed_by_cell = allowed_by_cell,
    original_assignment = original_assignment,
    analysis_spec = normalized_analysis$spec,
    runtime_policy = normalized_runtime_policy
  )
  original_failures <- .ab_constraint_failures(data, original_assignment)
  if (length(original_failures)) {
    return(structure(list(
      ok = FALSE,
      status = "invalid_input",
      errors = paste0(
        "original assignment violates registered feasibility: ",
        paste(original_failures, collapse = ";")
      ),
      warnings = unique(warnings),
      input_hashes = character(),
      data = NULL
    ), class = "ab_validation"))
  }

  data$input_hashes <- .ab_hash_payloads(data)
  structure(list(
    ok = TRUE,
    status = "valid",
    errors = character(),
    warnings = unique(warnings),
    input_hashes = data$input_hashes,
    data = data
  ), class = "ab_validation")
}
