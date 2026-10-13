.ab_fit_failure <- function(data, assignment, pseudobulk, message,
                            sample_expression = data.frame(),
                            pseudobulk_summary = data.frame()) {
  structure(list(
    status = "not_assessable",
    beta = NA_real_,
    margin = NA_real_,
    claim_holds = NA,
    overturn = FALSE,
    sample_expression = sample_expression,
    per_sample = sample_expression,
    pseudobulk_summary = pseudobulk_summary,
    error = as.character(message),
    assignment = assignment,
    assignment_id = if (!is.null(assignment)) {
      .ab_assignment_id(data, assignment)
    } else {
      NA_character_
    },
    pseudobulk = pseudobulk
  ), class = "ab_fit")
}

.ab_invoke_callback <- function(callback, x, assignment, pseudobulk) {
  callback_formals <- tryCatch(names(formals(callback)), error = function(e) NULL)
  arguments <- list(x = x, assignment = assignment, pseudobulk = pseudobulk)
  if (!is.null(callback_formals) && !"..." %in% callback_formals) {
    arguments <- arguments[intersect(names(arguments), callback_formals)]
  }
  do.call(callback, arguments)
}

.ab_callback_beta <- function(value) {
  if (inherits(value, "ab_fit") || is.list(value)) {
    if (!is.null(value$status) && identical(value$status, "not_assessable")) {
      message <- if (!is.null(value$error)) value$error else "callback returned not_assessable"
      stop(as.character(message), call. = FALSE)
    }
    value <- value$beta
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("callback must return one finite numeric beta or a list containing it",
         call. = FALSE)
  }
  as.numeric(value)
}

#' Evaluate the registered Layer-A expression claim
#'
#' @param x A successful `ab_validation` object.
#' @param assignment A named complete assignment, or `NULL` for original labels.
#' @param pseudobulk An optional matching `ab_pseudobulk`.
#' @param callback Optional test callback returning a beta or fit-like list.
#' @return An object of class `ab_fit`.
#' @export
fit_claim <- function(x, assignment = NULL, pseudobulk = NULL, callback = NULL) {
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    empty <- list(
      status = "not_assessable", beta = NA_real_, margin = NA_real_,
      claim_holds = NA, overturn = FALSE, sample_expression = data.frame(),
      per_sample = data.frame(), pseudobulk_summary = data.frame(),
      error = "x must be a successful ab_validation object",
      assignment = NULL, assignment_id = NA_character_, pseudobulk = NULL
    )
    return(structure(empty, class = "ab_fit"))
  }
  data <- x$data
  if (!is.null(data$analysis_spec)) {
    return(.ab_fit_failure(
      data,
      tryCatch(.ab_normalize_assignment(data, assignment), error = function(e) NULL),
      pseudobulk,
      "fit_claim is the portable Layer-A estimator; use fit_edger_claim for this analysis_spec"
    ))
  }
  assignment_error <- NULL
  normalized_assignment <- tryCatch(
    {
      if (is.null(assignment) && inherits(pseudobulk, "ab_pseudobulk")) {
        .ab_normalize_assignment(data, pseudobulk$assignment)
      } else {
        .ab_normalize_assignment(data, assignment)
      }
    },
    error = function(error) {
      assignment_error <<- conditionMessage(error)
      NULL
    }
  )
  if (!is.null(assignment_error)) {
    return(.ab_fit_failure(data, NULL, pseudobulk, assignment_error))
  }

  constraint_failures <- .ab_constraint_failures(data, normalized_assignment)
  if (length(constraint_failures)) {
    return(.ab_fit_failure(
      data, normalized_assignment, pseudobulk,
      paste0("assignment is infeasible: ", paste(constraint_failures, collapse = ";"))
    ))
  }

  if (is.null(pseudobulk)) {
    pseudobulk <- tryCatch(
      build_pseudobulk(x, normalized_assignment, method = "full"),
      error = function(error) error
    )
    if (inherits(pseudobulk, "error")) {
      return(.ab_fit_failure(
        data, normalized_assignment, NULL, conditionMessage(pseudobulk)
      ))
    }
  } else {
    if (!inherits(pseudobulk, "ab_pseudobulk")) {
      return(.ab_fit_failure(
        data, normalized_assignment, pseudobulk,
        "pseudobulk must be an ab_pseudobulk object"
      ))
    }
    current_hashes <- .ab_recompute_input_hashes(x)
    if (!identical(pseudobulk$input_hashes, current_hashes)) {
      return(.ab_fit_failure(
        data, normalized_assignment, pseudobulk,
        "pseudobulk immutable-input hashes do not match x"
      ))
    }
    if (!identical(
      unname(pseudobulk$assignment[data$cells$cell_id]),
      unname(normalized_assignment[data$cells$cell_id])
    )) {
      return(.ab_fit_failure(
        data, normalized_assignment, pseudobulk,
        "pseudobulk assignment does not match assignment"
      ))
    }
  }

  long_summary <- .ab_pseudobulk_columns(data)
  focal_label <- data$claim$focal_label
  focal_rows <- long_summary$label == focal_label
  focal_columns <- long_summary$pseudobulk_id[focal_rows]
  focal_meta <- long_summary[focal_rows, , drop = FALSE]
  target_count <- as.numeric(pseudobulk$counts[
    data$claim$focal_gene, focal_columns, drop = TRUE
  ])
  library_count <- as.numeric(Matrix::colSums(pseudobulk$counts[
    data$claim$feature_universe, focal_columns, drop = FALSE
  ]))
  cell_count <- vapply(
    focal_meta$sample_id,
    function(sample_id) sum(
      data$cells$sample_id == sample_id & normalized_assignment == focal_label
    ),
    integer(1L)
  )
  alpha <- data$claim$pseudocount_alpha
  feature_count <- length(data$claim$feature_universe)
  expression <- log2(
    data$claim$scale * (target_count + alpha) /
      (library_count + alpha * feature_count)
  )
  sample_expression <- data.frame(
    sample_id = focal_meta$sample_id,
    donor_id = focal_meta$donor_id,
    condition = focal_meta$condition,
    focal_label = focal_label,
    cell_count = unname(cell_count),
    target_count = target_count,
    library_count = library_count,
    expression_log2 = expression,
    stringsAsFactors = FALSE
  )
  sample_expression <- sample_expression[
    order(sample_expression$sample_id, method = "radix"), , drop = FALSE
  ]
  rownames(sample_expression) <- NULL
  pseudobulk_summary <- .ab_long_pseudobulk(data, pseudobulk)

  if (any(!is.finite(expression))) {
    return(.ab_fit_failure(
      data, normalized_assignment, pseudobulk,
      "registered estimator produced nonfinite per-sample expression",
      sample_expression, pseudobulk_summary
    ))
  }

  if (is.null(callback)) {
    case_values <- sample_expression$expression_log2[
      sample_expression$condition == data$claim$case_condition
    ]
    control_values <- sample_expression$expression_log2[
      sample_expression$condition == data$claim$control_condition
    ]
    beta <- mean(case_values) - mean(control_values)
  } else {
    callback_value <- tryCatch(
      .ab_invoke_callback(callback, x, normalized_assignment, pseudobulk),
      error = function(error) error
    )
    if (inherits(callback_value, "error")) {
      return(.ab_fit_failure(
        data, normalized_assignment, pseudobulk,
        paste0("callback failed: ", conditionMessage(callback_value)),
        sample_expression, pseudobulk_summary
      ))
    }
    beta_value <- tryCatch(.ab_callback_beta(callback_value), error = function(error) error)
    if (inherits(beta_value, "error")) {
      return(.ab_fit_failure(
        data, normalized_assignment, pseudobulk,
        paste0("callback failed: ", conditionMessage(beta_value)),
        sample_expression, pseudobulk_summary
      ))
    }
    beta <- beta_value
  }

  if (!is.finite(beta)) {
    return(.ab_fit_failure(
      data, normalized_assignment, pseudobulk,
      "registered estimator produced a nonfinite effect",
      sample_expression, pseudobulk_summary
    ))
  }
  margin <- data$claim$q * beta - data$claim$delta
  if (!is.finite(margin)) {
    return(.ab_fit_failure(
      data, normalized_assignment, pseudobulk,
      "registered claim margin is nonfinite",
      sample_expression, pseudobulk_summary
    ))
  }
  claim_holds <- margin > data$claim$tolerance
  overturn <- margin <= data$claim$tolerance

  structure(list(
    status = "evaluable",
    beta = as.numeric(beta),
    margin = as.numeric(margin),
    claim_holds = isTRUE(claim_holds),
    overturn = isTRUE(overturn),
    sample_expression = sample_expression,
    per_sample = sample_expression,
    pseudobulk_summary = pseudobulk_summary,
    error = NA_character_,
    assignment = normalized_assignment,
    assignment_id = .ab_assignment_id(data, normalized_assignment),
    pseudobulk = pseudobulk
  ), class = "ab_fit")
}
