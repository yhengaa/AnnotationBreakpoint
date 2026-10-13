.ab_edger_spec_expected <- function(control_reference = "control",
                                    case_level = "case") {
  list(
    schema_version = "1.0.0",
    analysis_id = "edger-ql-pseudobulk-v1",
    backend = list(package = "edgeR", version = "4.6.3", r_version = "4.5.0"),
    design = list(
      formula = "~ condition",
      control_reference = control_reference,
      case_level = case_level,
      coefficient_role = "case_vs_control",
      donor_fixed_effect = FALSE
    ),
    filter = list(
      `function` = "filterByExpr", min_count = 10, min_total_count = 15,
      large_n = 10, min_prop = 0.7
    ),
    normalization = list(
      `function` = "calcNormFactors", method = "TMM",
      rebuild_library_sizes = TRUE
    ),
    fit = list(
      `function` = "glmQLFit", dispersion = NULL, abundance_trend = TRUE,
      robust = FALSE, legacy = FALSE, top_proportion = NULL,
      prior_count = 0.125
    ),
    test = list(
      `function` = "glmQLFTest", coefficient = "case_vs_control",
      poisson_bound = TRUE
    ),
    endpoint = list(
      field = "logFC", units = "log2", coefficient_check_tolerance = 1e-12
    ),
    requirements = list(
      min_retained_features = 50L, min_donors_per_condition = 2L,
      require_positive_libraries = TRUE
    ),
    failure_policy = list(
      filtered_target = "not_assessable", rank_failure = "not_assessable",
      nonfinite_result = "not_assessable", warning = "record"
    ),
    nuisance_policy_version = "rebuild-all-label-dependent-steps-v1",
    cache_policy_version = "final-fit-canonical-key-v1"
  )
}

.ab_compare_spec <- function(value, expected, path = "analysis_spec") {
  errors <- character()
  if (is.list(expected)) {
    if (!is.list(value) || is.data.frame(value)) {
      return(paste0(path, " must be a named mapping"))
    }
    expected_names <- names(expected)
    observed_names <- names(value)
    if (is.null(observed_names) || !setequal(observed_names, expected_names)) {
      missing <- setdiff(expected_names, observed_names)
      extra <- setdiff(observed_names, expected_names)
      if (length(missing)) {
        errors <- c(errors, paste0(path, " is missing: ", paste(missing, collapse = ",")))
      }
      if (length(extra)) {
        errors <- c(errors, paste0(path, " has unsupported fields: ",
                                   paste(extra, collapse = ",")))
      }
    }
    for (name in intersect(expected_names, observed_names)) {
      errors <- c(errors, .ab_compare_spec(
        value[[name]], expected[[name]], paste0(path, "$", name)
      ))
    }
    return(errors)
  }
  if (is.null(expected)) {
    if (!is.null(value)) errors <- c(errors, paste0(path, " must be null"))
    return(errors)
  }
  if (is.logical(expected)) {
    if (length(value) != 1L || is.na(value) ||
        !identical(as.logical(value), expected)) {
      errors <- c(errors, paste0(path, " does not match the frozen policy"))
    }
  } else if (is.numeric(expected)) {
    number <- suppressWarnings(as.numeric(value))
    if (length(number) != 1L || !is.finite(number) || number != as.numeric(expected)) {
      errors <- c(errors, paste0(path, " does not match the frozen policy"))
    }
  } else if (is.character(expected)) {
    if (length(value) != 1L || is.na(value) ||
        !identical(as.character(value), expected)) {
      errors <- c(errors, paste0(path, " does not match the frozen policy"))
    }
  }
  errors
}

.ab_normalize_analysis_spec <- function(analysis_spec) {
  if (is.null(analysis_spec)) {
    return(list(ok = TRUE, errors = character(), spec = NULL))
  }
  design <- if (is.list(analysis_spec) && !is.data.frame(analysis_spec)) {
    analysis_spec$design
  } else {
    NULL
  }
  control_reference <- if (is.list(design)) design$control_reference else NULL
  case_level <- if (is.list(design)) design$case_level else NULL
  valid_level <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
  }
  level_errors <- character()
  if (!valid_level(control_reference)) {
    level_errors <- c(
      level_errors,
      "analysis_spec$design$control_reference must be one nonempty string"
    )
  }
  if (!valid_level(case_level)) {
    level_errors <- c(
      level_errors,
      "analysis_spec$design$case_level must be one nonempty string"
    )
  }
  if (valid_level(control_reference) && valid_level(case_level) &&
      identical(control_reference, case_level)) {
    level_errors <- c(
      level_errors,
      "analysis_spec design control_reference and case_level must be distinct"
    )
  }
  expected <- if (!length(level_errors)) {
    .ab_edger_spec_expected(control_reference, case_level)
  } else {
    .ab_edger_spec_expected()
  }
  errors <- c(level_errors, .ab_compare_spec(analysis_spec, expected))
  list(ok = !length(errors), errors = unique(errors), spec = if (!length(errors)) {
    expected
  } else {
    NULL
  })
}

.ab_authoritative_identity <- function(data) {
  if (is.null(data$analysis_spec)) {
    return(list(
      type = "authoritative",
      evaluation_id = "portable-log-proportion-contrast-v1",
      analysis_id = "portable-log-proportion-contrast-v1",
      backend = "Layer-A",
      backend_version = NA_character_,
      specification_hash = unname(data$input_hashes[["specification_hash"]])
    ))
  }
  list(
    type = "authoritative",
    evaluation_id = data$analysis_spec$analysis_id,
    analysis_id = data$analysis_spec$analysis_id,
    backend = data$analysis_spec$backend$package,
    backend_version = data$analysis_spec$backend$version,
    nuisance_policy_version = data$analysis_spec$nuisance_policy_version,
    cache_policy_version = data$analysis_spec$cache_policy_version,
    specification_hash = unname(data$input_hashes[["specification_hash"]])
  )
}

.ab_function_hash <- function(fun) {
  payload <- list(
    formals = tryCatch(deparse(formals(fun), width.cutoff = 500L),
                       error = function(e) character()),
    body = tryCatch(deparse(body(fun), width.cutoff = 500L),
                    error = function(e) deparse(fun, width.cutoff = 500L))
  )
  .ab_hash(payload)
}

.ab_resolve_fit <- function(x, fit_fun = NULL, fit_identity = NULL) {
  data <- x$data
  if (is.null(fit_fun)) {
    fun <- if (is.null(data$analysis_spec)) fit_claim else fit_edger_claim
    return(list(fun = fun, identity = .ab_authoritative_identity(data)))
  }
  if (!is.function(fit_fun)) stop("fit_fun must be a function", call. = FALSE)
  if (!is.null(data$analysis_spec) && is.null(fit_identity)) {
    stop("M2 custom fit functions require an explicit versioned fit_identity",
         call. = FALSE)
  }
  function_hash <- .ab_function_hash(fit_fun)
  label <- if (is.null(fit_identity)) {
    "anonymous-custom-fit"
  } else if (is.character(fit_identity) && length(fit_identity) == 1L &&
             !is.na(fit_identity) && nzchar(fit_identity)) {
    fit_identity
  } else {
    stop("fit_identity must be NULL or one nonempty string", call. = FALSE)
  }
  list(
    fun = fit_fun,
    identity = list(
      type = "custom",
      evaluation_id = paste0("custom:", label, ":", function_hash),
      label = label,
      function_hash = function_hash,
      authoritative_analysis_id = .ab_authoritative_identity(data)$analysis_id,
      specification_hash = unname(data$input_hashes[["specification_hash"]])
    )
  )
}

.ab_fit_cache_key <- function(data, assignment, evaluation_identity) {
  .ab_hash(list(
    schema = "annotation-breakpoint-final-fit-cache-v1",
    assignment_id = .ab_assignment_id(data, assignment),
    input_hash = unname(data$input_hashes[["input_hash"]]),
    claim_hash = unname(data$input_hashes[["claim"]]),
    specification_hash = unname(data$input_hashes[["specification_hash"]]),
    constraint_hash = unname(data$input_hashes[["constraint_hash"]]),
    nuisance_policy_version = if (is.null(data$analysis_spec)) {
      "portable-layer-a-v1"
    } else {
      data$analysis_spec$nuisance_policy_version
    },
    evaluation_identity = evaluation_identity
  ))
}
