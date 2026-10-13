.ab_write_tsv <- function(value, path) {
  utils::write.table(
    value,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    eol = "\n",
    fileEncoding = "UTF-8"
  )
}

.ab_write_json <- function(value, path, pretty = TRUE) {
  jsonlite::write_json(
    value,
    path = path,
    auto_unbox = TRUE,
    dataframe = "rows",
    matrix = "rowmajor",
    null = "null",
    na = "null",
    digits = NA,
    pretty = pretty
  )
  write("", file = path, append = TRUE)
}

.ab_yaml_ready <- function(value) {
  if (is.atomic(value) && !is.null(names(value))) {
    return(lapply(as.list(value), .ab_yaml_ready))
  }
  if (is.list(value)) {
    return(lapply(value, .ab_yaml_ready))
  }
  value
}

.ab_named_or_na <- function(value, mode = c("integer", "double", "character")) {
  mode <- match.arg(mode)
  if (is.null(value) || !length(value)) {
    return(switch(mode,
      integer = NA_integer_,
      double = NA_real_,
      character = NA_character_
    ))
  }
  switch(mode,
    integer = as.integer(value[[1L]]),
    double = as.numeric(value[[1L]]),
    character = as.character(value[[1L]])
  )
}

.ab_flat_candidates <- function(candidates) {
  if (!is.data.frame(candidates)) return(data.frame())
  candidates[, setdiff(names(candidates), "assignment"), drop = FALSE]
}

.ab_witness_json <- function(witness) {
  if (is.null(witness)) return(NULL)
  assignment <- witness$assignment
  list(
    witness_id = witness$witness_id,
    assignment_id = witness$assignment_id,
    changes = witness$changes,
    cost = witness$cost,
    beta = witness$beta,
    margin = witness$margin,
    overturn = witness$overturn,
    assignment = data.frame(
      cell_id = names(assignment),
      assigned_label = unname(assignment),
      stringsAsFactors = FALSE
    ),
    changed_cells = witness$changed_cells,
    input_hashes = as.list(witness$input_hashes)
  )
}

.ab_result_json <- function(result) {
  original <- result$original_fit
  list(
    schema_version = "1.0.0",
    status = result$status,
    claim = result$claim,
    constraints = result$constraints,
    analysis_spec = result$analysis_spec,
    analysis_id = result$analysis_id,
    backend_identity = result$backend_identity,
    evaluation_identity = result$evaluation_identity,
    result_id = .ab_result_identity(result),
    original_fit = if (is.null(original)) NULL else list(
      status = original$status,
      beta = original$beta,
      margin = original$margin,
      claim_holds = original$claim_holds,
      overturn = original$overturn,
      error = original$error,
      failure_reason = original$failure_reason,
      warnings = original$warnings,
      runtime = original$runtime,
      sample_expression = original$sample_expression
    ),
    assignment_accounting = list(
      total = result$total_assignment_count,
      feasible = result$feasible_assignment_count,
      infeasible = result$infeasible_assignment_count,
      evaluable = result$evaluable_assignment_count,
      unevaluable = result$unevaluable_assignment_count,
      not_evaluated = result$unevaluated_assignment_count,
      overturn = result$overturn_assignment_count
    ),
    best_verified_cost = result$best_verified_cost,
    exact_minimum_cost = result$exact_minimum_cost,
    certified_lower_bound = result$certified_lower_bound,
    witness_id = result$witness_id,
    witness = .ab_witness_json(result$witness),
    search_complete = result$search_complete,
    max_evaluations = if (is.infinite(result$max_evaluations)) {
      "Inf"
    } else {
      result$max_evaluations
    },
    termination_reason = result$termination_reason,
    method = if (!is.null(result$method)) result$method else "exact_enumeration",
    seed = if (!is.null(result$seed)) result$seed else NULL,
    direction = if (!is.null(result$direction)) result$direction else "all",
    elapsed_seconds = if (!is.null(result$elapsed_seconds)) {
      result$elapsed_seconds
    } else NULL,
    fit_accounting = if (!is.null(result$fit_accounting)) {
      result$fit_accounting
    } else NULL,
    proxy_prefix = if (!is.null(result$proxy_prefix)) {
      result$proxy_prefix
    } else NULL,
    verification_embedded = if (!is.null(result$verification_embedded)) {
      result$verification_embedded
    } else FALSE,
    verification_fit_count = if (!is.null(result$verification_fit_count)) {
      result$verification_fit_count
    } else NULL,
    replay_attempts = if (!is.null(result$replay_attempts)) {
      result$replay_attempts
    } else NULL,
    search_elapsed_seconds = if (!is.null(result$search_elapsed_seconds)) {
      result$search_elapsed_seconds
    } else NULL,
    replay_elapsed_seconds = if (!is.null(result$replay_elapsed_seconds)) {
      result$replay_elapsed_seconds
    } else NULL,
    all_in_elapsed_seconds = if (!is.null(result$all_in_elapsed_seconds)) {
      result$all_in_elapsed_seconds
    } else NULL,
    search_path = if (!is.null(result$search_path)) result$search_path else NULL,
    errors = result$errors,
    input_hashes = as.list(result$input_hashes),
    effects_by_cost = result$effects_by_cost,
    assignments = .ab_flat_candidates(result$candidates)
  )
}

.ab_verification_json <- function(verification) {
  list(
    schema_version = "1.0.0",
    ok = verification$ok,
    errors = verification$errors,
    recomputed_cost = verification$recomputed_cost,
    recomputed_effect = verification$recomputed_effect,
    recomputed_status = verification$recomputed_status,
    witness_id = verification$witness_id,
    changed_cells = verification$changed_cells,
    hash_matches = verification$hash_matches,
    input_hashes = as.list(verification$input_hashes),
    result_identity = verification$result_identity,
    verification_fit_count = verification$verification_fit_count,
    verification_elapsed_seconds = verification$verification_elapsed_seconds,
    evaluation_identity = verification$evaluation_identity,
    before_donor_pseudobulks = verification$before_donor_pseudobulks,
    after_donor_pseudobulks = verification$after_donor_pseudobulks
  )
}

.ab_markdown_value <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[[1L]])) return("NA")
  if (is.numeric(value)) return(format(value[[1L]], digits = 17, scientific = FALSE))
  as.character(value[[1L]])
}

#' Render deterministic result and verification artifacts
#'
#' @param result An `ab_result`.
#' @param verification An `ab_verification` for the same result.
#' @param output_dir Caller-supplied output directory.
#' @param run_metadata Optional named metadata recorded verbatim.
#' @return Paths and SHA-256 hashes, invisibly classed as `ab_report_artifacts`.
#' @export
render_report <- function(result, verification, output_dir,
                          run_metadata = list()) {
  if (!inherits(result, "ab_result")) {
    stop("result must be an ab_result object", call. = FALSE)
  }
  if (!inherits(verification, "ab_verification")) {
    stop("verification must be an ab_verification object", call. = FALSE)
  }
  current_result_identity <- .ab_result_identity(result)
  if (!isTRUE(verification$ok)) {
    stop("verification must be successful before rendering", call. = FALSE)
  }
  if (is.null(verification$result_identity) ||
      !identical(verification$result_identity, current_result_identity)) {
    stop("verification belongs to a different or modified result", call. = FALSE)
  }
  if (!identical(verification$input_hashes, result$input_hashes)) {
    stop("verification immutable hashes do not match result", call. = FALSE)
  }
  if (!identical(verification$recomputed_status, result$status)) {
    stop("verification status does not match result status", call. = FALSE)
  }
  reported_verification_fits <- if (
    isTRUE(result$verification_embedded) &&
      !is.null(result$verification_fit_count)
  ) {
    as.integer(result$verification_fit_count)
  } else {
    as.integer(verification$verification_fit_count)
  }
  result_witness_id <- if (is.null(result$witness)) NULL else result$witness_id
  if (!identical(verification$witness_id, result_witness_id)) {
    stop("verification witness identity does not match result", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path", call. = FALSE)
  }
  if (!is.list(run_metadata)) {
    stop("run_metadata must be a list", call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("failed to create output_dir", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)

  candidates <- .ab_flat_candidates(result$candidates)
  original_beta <- if (is.null(result$original_fit)) {
    NA_real_
  } else {
    result$original_fit$beta
  }
  original_margin <- if (is.null(result$original_fit)) {
    NA_real_
  } else {
    result$original_fit$margin
  }
  original_failure_reason <- if (is.null(result$original_fit) ||
                                 is.null(result$original_fit$failure_reason)) {
    NA_character_
  } else as.character(result$original_fit$failure_reason)
  original_error <- if (is.null(result$original_fit) ||
                        is.null(result$original_fit$error)) {
    NA_character_
  } else as.character(result$original_fit$error)
  original_warnings <- if (is.null(result$original_fit) ||
                           !length(result$original_fit$warnings)) {
    NA_character_
  } else paste(unique(as.character(result$original_fit$warnings)), collapse = " | ")
  runtime <- if (!is.null(result$original_fit) &&
                 !is.null(result$original_fit$runtime)) {
    result$original_fit$runtime
  } else {
    list(
      policy = result$backend_identity$runtime_policy,
      expected_r_version = result$backend_identity$expected_r_version,
      observed_r_version = result$backend_identity$observed_r_version,
      expected_backend_version = result$backend_identity$expected_backend_version,
      observed_backend_version = result$backend_identity$backend_version,
      api_supported = result$backend_identity$runtime_api_supported,
      compatible = result$backend_identity$runtime_compatible
    )
  }
  claim_results <- data.frame(
    status = result$status,
    evidence_label = result$status,
    focal_gene = .ab_named_or_na(result$claim$focal_gene, "character"),
    focal_label = .ab_named_or_na(result$claim$focal_label, "character"),
    case_condition = .ab_named_or_na(result$claim$case_condition, "character"),
    control_condition = .ab_named_or_na(
      result$claim$control_condition, "character"
    ),
    direction_q = .ab_named_or_na(result$claim$q, "integer"),
    threshold_delta_log2 = .ab_named_or_na(result$claim$delta, "double"),
    comparison_tolerance = .ab_named_or_na(
      result$claim$tolerance, "double"
    ),
    estimator_id = .ab_named_or_na(result$claim$estimator_id, "character"),
    analysis_id = .ab_named_or_na(result$analysis_id, "character"),
    backend = .ab_named_or_na(result$backend_identity$backend, "character"),
    backend_version = .ab_named_or_na(
      result$backend_identity$backend_version, "character"
    ),
    evaluation_id = .ab_named_or_na(
      result$evaluation_identity$evaluation_id, "character"
    ),
    method = if (!is.null(result$method)) result$method else "exact_enumeration",
    direction = if (!is.null(result$direction)) result$direction else "all",
    original_beta_log2 = original_beta,
    original_margin_log2 = original_margin,
    original_fit_status = .ab_named_or_na(
      if (is.null(result$original_fit)) NULL else result$original_fit$status,
      "character"
    ),
    original_failure_reason = original_failure_reason,
    original_error = original_error,
    original_warnings = original_warnings,
    runtime_policy = .ab_named_or_na(runtime$policy, "character"),
    expected_r_version = .ab_named_or_na(runtime$expected_r_version, "character"),
    observed_r_version = .ab_named_or_na(runtime$observed_r_version, "character"),
    expected_backend_version = .ab_named_or_na(
      runtime$expected_backend_version, "character"
    ),
    observed_backend_version = .ab_named_or_na(
      runtime$observed_backend_version, "character"
    ),
    runtime_api_supported = if (is.null(runtime$api_supported)) NA else
      isTRUE(runtime$api_supported),
    runtime_compatible = if (is.null(runtime$compatible)) NA else
      isTRUE(runtime$compatible),
    best_verified_cost = .ab_named_or_na(result$best_verified_cost, "integer"),
    exact_minimum_cost = .ab_named_or_na(result$exact_minimum_cost, "integer"),
    certified_lower_bound = .ab_named_or_na(result$certified_lower_bound, "integer"),
    total_assignment_count = result$total_assignment_count,
    feasible_assignment_count = result$feasible_assignment_count,
    infeasible_assignment_count = result$infeasible_assignment_count,
    evaluable_assignment_count = result$evaluable_assignment_count,
    unevaluable_assignment_count = result$unevaluable_assignment_count,
    unevaluated_assignment_count = result$unevaluated_assignment_count,
    overturn_assignment_count = result$overturn_assignment_count,
    witness_id = .ab_named_or_na(result$witness_id, "character"),
    search_complete = result$search_complete,
    termination_reason = result$termination_reason,
    search_fit_count = if (!is.null(result$attempted_fit_count)) {
      result$attempted_fit_count
    } else sum(result$candidates$evaluated),
    verification_fit_count = reported_verification_fits,
    search_elapsed_seconds = if (!is.null(result$elapsed_seconds)) {
      result$elapsed_seconds
    } else NA_real_,
    verification_elapsed_seconds = verification$verification_elapsed_seconds,
    stringsAsFactors = FALSE
  )
  witness_table <- if (is.null(result$witness)) {
    verification$changed_cells[FALSE, , drop = FALSE]
  } else {
    result$witness$changed_cells
  }
  input_hashes <- data.frame(
    hash_name = names(result$input_hashes),
    sha256 = unname(result$input_hashes),
    stringsAsFactors = FALSE
  )
  method_log <- if (!is.null(result$candidate_log) &&
                    is.data.frame(result$candidate_log) &&
                    nrow(result$candidate_log)) {
    .ab_flat_candidates(result$candidate_log)
  } else candidates
  fit_accounting <- if (!is.null(result$fit_accounting)) {
    result$fit_accounting
  } else {
    data.frame(
      proposals = nrow(candidates),
      duplicate_proposals = 0L,
      cache_hits = 0L,
      unique_assignments = nrow(candidates),
      infeasible_proposals = sum(!candidates$feasible),
      attempted_fits = sum(candidates$evaluated),
      successful_fits = sum(candidates$evaluable),
      failed_fits = sum(candidates$evaluated & !candidates$evaluable),
      verification_fits = verification$verification_fit_count,
      stringsAsFactors = FALSE
    )
  }
  fit_accounting$verification_fits <- reported_verification_fits

  paths <- c(
    claim_results = file.path(output_dir, "claim_results.tsv"),
    assignments = file.path(output_dir, "assignments.tsv"),
    effects_by_cost = file.path(output_dir, "effects_by_cost.tsv"),
    witnesses = file.path(output_dir, "witnesses.tsv"),
    before_pseudobulks = file.path(output_dir, "before_pseudobulks.tsv"),
    after_pseudobulks = file.path(output_dir, "after_pseudobulks.tsv"),
    method_log = file.path(output_dir, "method_log.tsv"),
    fit_accounting = file.path(output_dir, "fit_accounting.tsv"),
    search_path = file.path(output_dir, "search_path.tsv"),
    input_hashes = file.path(output_dir, "input_hashes.tsv"),
    analysis_spec_json = file.path(output_dir, "analysis_spec.json"),
    analysis_spec = file.path(output_dir, "analysis_spec.yaml"),
    constraints = file.path(output_dir, "constraints.yaml"),
    result_json = file.path(output_dir, "result.json"),
    verification_json = file.path(output_dir, "verification.json"),
    search_log = file.path(output_dir, "search_log.jsonl"),
    report = file.path(output_dir, "report.md"),
    manifest = file.path(output_dir, "run_manifest.json")
  )
  .ab_write_tsv(claim_results, paths[["claim_results"]])
  .ab_write_tsv(candidates, paths[["assignments"]])
  .ab_write_tsv(result$effects_by_cost, paths[["effects_by_cost"]])
  .ab_write_tsv(witness_table, paths[["witnesses"]])
  .ab_write_tsv(
    verification$before_donor_pseudobulks,
    paths[["before_pseudobulks"]]
  )
  .ab_write_tsv(
    verification$after_donor_pseudobulks,
    paths[["after_pseudobulks"]]
  )
  .ab_write_tsv(method_log, paths[["method_log"]])
  .ab_write_tsv(fit_accounting, paths[["fit_accounting"]])
  .ab_write_tsv(
    if (!is.null(result$search_path)) result$search_path else data.frame(),
    paths[["search_path"]]
  )
  .ab_write_tsv(input_hashes, paths[["input_hashes"]])
  .ab_write_json(result$analysis_spec, paths[["analysis_spec_json"]])
  analysis_yaml <- if (!is.null(result$analysis_spec)) {
    .ab_yaml_ready(result$analysis_spec)
  } else {
    list(
      schema_version = "1.0.0",
      analysis_id = result$claim$estimator_id,
      claim = .ab_yaml_ready(result$claim)
    )
  }
  yaml::write_yaml(
    analysis_yaml,
    paths[["analysis_spec"]]
  )
  yaml::write_yaml(
    list(
      schema_version = "1.0.0",
      constraints = .ab_yaml_ready(result$constraints)
    ),
    paths[["constraints"]]
  )
  .ab_write_json(.ab_result_json(result), paths[["result_json"]])
  .ab_write_json(.ab_verification_json(verification), paths[["verification_json"]])

  search_lines <- character()
  if (nrow(candidates)) {
    search_lines <- vapply(seq_len(nrow(candidates)), function(index) {
      jsonlite::toJSON(
        c(list(event = "candidate"), as.list(candidates[index, , drop = FALSE])),
        auto_unbox = TRUE,
        dataframe = "rows",
        null = "null",
        na = "null",
        digits = NA,
        pretty = FALSE
      )
    }, character(1L))
  }
  search_lines <- c(search_lines, jsonlite::toJSON(
    list(
      event = "termination",
      status = result$status,
      search_complete = result$search_complete,
      termination_reason = result$termination_reason
    ),
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    pretty = FALSE
  ))
  writeLines(search_lines, paths[["search_log"]], useBytes = TRUE)

  markdown <- c(
    "# Annotation Breakpoint report",
    "",
    paste0("Status: `", result$status, "`."),
    "",
    paste0(
      "Registered claim: `", result$claim$focal_gene, "` in population `",
      result$claim$focal_label, "`, `", result$claim$case_condition,
      "` versus `", result$claim$control_condition, "`, q=`",
      result$claim$q, "`, delta=`", result$claim$delta,
      "`, tolerance=`", result$claim$tolerance, "`."
    ),
    paste0("Estimator: `", result$claim$estimator_id, "`."),
    paste0(
      "Witness moves: `",
      if (is.null(result$witness)) "none" else result$witness$changes,
      "`."
    ),
    "",
    paste0("Original beta (log2): `", .ab_markdown_value(original_beta), "`."),
    paste0("Original fit status: `",
           .ab_markdown_value(claim_results$original_fit_status), "`."),
    paste0("Original fit failure reason: `",
           .ab_markdown_value(original_failure_reason), "`."),
    paste0("Original fit error: `", .ab_markdown_value(original_error), "`."),
    paste0("Original fit warnings: `",
           .ab_markdown_value(original_warnings), "`."),
    paste0("Best verified cost: `",
           .ab_markdown_value(result$best_verified_cost), "`."),
    paste0("Exact minimum cost: `",
           .ab_markdown_value(result$exact_minimum_cost), "`."),
    paste0("Termination reason: `", result$termination_reason, "`."),
    "",
    paste0(
      "Assignments: ", result$total_assignment_count,
      " total; ", result$feasible_assignment_count,
      " feasible; ", result$evaluable_assignment_count,
      " evaluable; ", result$unevaluable_assignment_count,
      " unevaluable; ", result$unevaluated_assignment_count,
      " not evaluated."
    ),
    "",
    paste0("Witness verification: `", verification$ok, "`."),
    paste0("Analysis: `", .ab_markdown_value(result$analysis_id), "`; evaluation: `",
           .ab_markdown_value(result$evaluation_identity$evaluation_id), "`."),
    paste0(
      "Runtime policy: `", .ab_markdown_value(runtime$policy),
      "`; R expected/observed: `", .ab_markdown_value(runtime$expected_r_version),
      "` / `", .ab_markdown_value(runtime$observed_r_version),
      "`; edgeR expected/observed: `",
      .ab_markdown_value(runtime$expected_backend_version), "` / `",
      .ab_markdown_value(runtime$observed_backend_version), "`."
    ),
    paste0("Search fits: `", claim_results$search_fit_count,
           "`; independent verification fits: `", verification$verification_fit_count,
           "`."),
    "",
    "A counterexample establishes sensitivity to the declared annotation assumptions. It does not prove misannotation or invalidate the biological effect."
  )
  writeLines(markdown, paths[["report"]], useBytes = TRUE)

  pre_manifest <- paths[names(paths) != "manifest"]
  artifact_hashes <- vapply(pre_manifest, .ab_file_hash, character(1L))
  metadata <- run_metadata
  if (length(metadata) && !is.null(names(metadata))) {
    metadata <- metadata[order(names(metadata), method = "radix")]
  }
  manifest <- list(
    schema_version = "1.0.0",
    status = result$status,
    verification_ok = verification$ok,
    run_metadata = metadata,
    input_hashes = as.list(result$input_hashes),
    artifacts = lapply(names(pre_manifest), function(name) list(
      name = name,
      file = basename(pre_manifest[[name]]),
      sha256 = unname(artifact_hashes[[name]])
    ))
  )
  .ab_write_json(manifest, paths[["manifest"]])
  hashes <- vapply(paths, .ab_file_hash, character(1L))
  structure(list(
    output_dir = output_dir,
    paths = paths,
    hashes = hashes
  ), class = "ab_report_artifacts")
}
