.ab_edger_failure <- function(data, assignment, pseudobulk, reason, message,
                              warnings = character(), diagnostics = list(),
                              started = NULL) {
  elapsed <- if (is.null(started)) NA_real_ else unname(proc.time()[["elapsed"]] - started)
  out <- list(
    status = "not_assessable",
    beta = NA_real_,
    margin = NA_real_,
    claim_holds = NA,
    overturn = FALSE,
    sample_expression = if (!is.null(diagnostics$sample_summary)) {
      diagnostics$sample_summary
    } else data.frame(),
    per_sample = if (!is.null(diagnostics$sample_summary)) {
      diagnostics$sample_summary
    } else data.frame(),
    pseudobulk_summary = if (inherits(pseudobulk, "ab_pseudobulk")) {
      .ab_long_pseudobulk(data, pseudobulk)
    } else data.frame(),
    error = as.character(message),
    failure_reason = as.character(reason),
    warnings = unique(as.character(warnings)),
    assignment = assignment,
    assignment_id = if (!is.null(assignment)) .ab_assignment_id(data, assignment) else NA_character_,
    pseudobulk = pseudobulk,
    analysis_id = data$analysis_spec$analysis_id,
    backend = data$analysis_spec$backend$package,
    backend_version = data$analysis_spec$backend$version,
    evaluation_identity = .ab_authoritative_identity(data),
    fit_elapsed_seconds = elapsed,
    authoritative_rebuild = TRUE
  )
  for (name in names(diagnostics)) out[[name]] <- diagnostics[[name]]
  structure(out, class = "ab_fit")
}

#' Fit the frozen authoritative edgeR quasi-likelihood claim
#'
#' @param x A successful M2 `ab_validation` object.
#' @param assignment A complete assignment, or `NULL` for original labels.
#' @param pseudobulk Optional matching aggregate; the authoritative fit still rebuilds it.
#' @return An `ab_fit` with edgeR diagnostics and registered claim fields.
fit_edger_claim <- function(x, assignment = NULL, pseudobulk = NULL) {
  started <- unname(proc.time()[["elapsed"]])
  if (!inherits(x, "ab_validation") || !isTRUE(x$ok) || is.null(x$data)) {
    return(structure(list(
      status = "not_assessable", beta = NA_real_, margin = NA_real_,
      claim_holds = NA, overturn = FALSE, sample_expression = data.frame(),
      per_sample = data.frame(), pseudobulk_summary = data.frame(),
      error = "x must be a successful ab_validation object",
      failure_reason = "invalid_validation", warnings = character(),
      assignment = NULL, assignment_id = NA_character_, pseudobulk = NULL,
      analysis_id = NA_character_, backend = "edgeR",
      backend_version = NA_character_, evaluation_identity = NULL,
      fit_elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
      authoritative_rebuild = FALSE
    ), class = "ab_fit"))
  }
  data <- x$data
  if (is.null(data$analysis_spec) ||
      !identical(data$analysis_spec$analysis_id, "edger-ql-pseudobulk-v1")) {
    placeholder <- data
    placeholder$analysis_spec <- .ab_edger_spec_expected()
    return(.ab_edger_failure(
      placeholder, NULL, NULL, "m2_analysis_spec_required",
      "fit_edger_claim requires the validated M2 edgeR analysis specification",
      started = started
    ))
  }
  assignment <- tryCatch(
    .ab_normalize_assignment(data, assignment),
    error = function(e) e
  )
  if (inherits(assignment, "error")) {
    return(.ab_edger_failure(
      data, NULL, NULL, "invalid_assignment", conditionMessage(assignment),
      started = started
    ))
  }
  failures <- .ab_constraint_failures(data, assignment)
  if (length(failures)) {
    return(.ab_edger_failure(
      data, assignment, NULL, "infeasible_assignment",
      paste(failures, collapse = ";"), started = started
    ))
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    return(.ab_edger_failure(
      data, assignment, NULL, "backend_unavailable",
      "edgeR is not installed in the active library path", started = started
    ))
  }
  observed_version <- as.character(utils::packageVersion("edgeR"))
  if (!identical(observed_version, data$analysis_spec$backend$version)) {
    return(.ab_edger_failure(
      data, assignment, NULL, "backend_version_mismatch",
      paste0("expected edgeR ", data$analysis_spec$backend$version,
             ", observed ", observed_version), started = started
    ))
  }
  observed_r <- paste(R.version$major, R.version$minor, sep = ".")
  if (!identical(observed_r, data$analysis_spec$backend$r_version)) {
    return(.ab_edger_failure(
      data, assignment, NULL, "r_version_mismatch",
      paste0("expected R ", data$analysis_spec$backend$r_version,
             ", observed ", observed_r), started = started
    ))
  }

  rebuilt <- tryCatch(build_pseudobulk(x, assignment, method = "full"),
                      error = function(e) e)
  if (inherits(rebuilt, "error")) {
    return(.ab_edger_failure(
      data, assignment, NULL, "pseudobulk_rebuild_failure",
      conditionMessage(rebuilt), started = started
    ))
  }
  if (!is.null(pseudobulk)) {
    if (!inherits(pseudobulk, "ab_pseudobulk") ||
        !identical(pseudobulk$input_hashes, rebuilt$input_hashes) ||
        !identical(pseudobulk$assignment, rebuilt$assignment)) {
      return(.ab_edger_failure(
        data, assignment, rebuilt, "supplied_pseudobulk_mismatch",
        "supplied pseudobulk does not match the authoritative rebuild",
        started = started
      ))
    }
  }

  focal <- data$claim$focal_label
  focal_meta <- rebuilt$column_metadata[
    rebuilt$column_metadata$label == focal, , drop = FALSE
  ]
  focal_meta <- focal_meta[match(data$sample_ids, focal_meta$sample_id), , drop = FALSE]
  focal_counts <- as.matrix(rebuilt$counts[
    data$claim$feature_universe, focal_meta$pseudobulk_id, drop = FALSE
  ])
  colnames(focal_counts) <- focal_meta$sample_id
  raw_library_sizes <- colSums(focal_counts)
  sample_summary <- data.frame(
    sample_id = focal_meta$sample_id,
    donor_id = focal_meta$donor_id,
    condition = focal_meta$condition,
    focal_label = focal,
    cell_count = focal_meta$cell_count,
    raw_library_size = as.numeric(raw_library_sizes),
    stringsAsFactors = FALSE
  )
  diagnostics <- list(
    sample_summary = sample_summary,
    raw_library_sizes = setNames(as.numeric(raw_library_sizes), focal_meta$sample_id)
  )
  if (any(!is.finite(raw_library_sizes)) || any(raw_library_sizes <= 0)) {
    return(.ab_edger_failure(
      data, assignment, rebuilt, "nonpositive_library",
      "focal pseudobulk libraries must be finite and positive",
      diagnostics = diagnostics, started = started
    ))
  }

  condition <- factor(
    focal_meta$condition,
    levels = c(data$claim$control_condition, data$claim$case_condition)
  )
  design <- stats::model.matrix(~ condition)
  colnames(design) <- c("(Intercept)", "case_vs_control")
  rownames(design) <- focal_meta$sample_id
  diagnostics$design <- design
  diagnostics$design_formula <- "~ condition"
  diagnostics$control_reference <- data$claim$control_condition
  diagnostics$case_level <- data$claim$case_condition
  diagnostics$coefficient_index <- 2L
  diagnostics$coefficient_name <- "case_vs_control"
  if (qr(design)$rank != ncol(design)) {
    return(.ab_edger_failure(
      data, assignment, rebuilt, "rank_failure",
      "the registered condition design is not full rank",
      diagnostics = diagnostics, started = started
    ))
  }

  captured_warnings <- character()
  fit_value <- tryCatch(
    withCallingHandlers({
      dge <- edgeR::DGEList(counts = focal_counts)
      keep <- edgeR::filterByExpr(
        dge,
        design = design,
        min.count = data$analysis_spec$filter$min_count,
        min.total.count = data$analysis_spec$filter$min_total_count,
        large.n = data$analysis_spec$filter$large_n,
        min.prop = data$analysis_spec$filter$min_prop
      )
      if (!isTRUE(keep[[data$claim$focal_gene]])) {
        stop(structure(list(message = "registered target was removed by filterByExpr"),
                       class = c("ab_target_filtered", "error", "condition")))
      }
      retained <- names(keep)[keep]
      if (length(retained) < data$analysis_spec$requirements$min_retained_features) {
        stop(structure(list(message = paste0(
          "only ", length(retained), " features survived filterByExpr"
        )), class = c("ab_insufficient_features", "error", "condition")))
      }
      dge <- dge[keep, , keep.lib.sizes = FALSE]
      if (any(!is.finite(dge$samples$lib.size)) || any(dge$samples$lib.size <= 0)) {
        stop(structure(list(message = "retained libraries are nonpositive"),
                       class = c("ab_nonpositive_library", "error", "condition")))
      }
      dge <- edgeR::calcNormFactors(
        dge,
        method = data$analysis_spec$normalization$method
      )
      effective <- dge$samples$lib.size * dge$samples$norm.factors
      if (any(!is.finite(dge$samples$norm.factors)) ||
          any(dge$samples$norm.factors <= 0) ||
          any(!is.finite(effective)) || any(effective <= 0)) {
        stop(structure(list(message = "TMM normalization factors are invalid"),
                       class = c("ab_normalization_failure", "error", "condition")))
      }
      fit <- edgeR::glmQLFit(
        dge,
        design = design,
        dispersion = NULL,
        abundance.trend = data$analysis_spec$fit$abundance_trend,
        robust = data$analysis_spec$fit$robust,
        legacy = data$analysis_spec$fit$legacy,
        top.proportion = data$analysis_spec$fit$top_proportion,
        prior.count = data$analysis_spec$fit$prior_count
      )
      test <- edgeR::glmQLFTest(
        fit,
        coef = 2L,
        poisson.bound = data$analysis_spec$test$poisson_bound
      )
      list(dge = dge, keep = keep, retained = retained, fit = fit,
           test = test, effective = effective)
    }, warning = function(warning) {
      captured_warnings <<- c(captured_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }),
    error = function(error) error
  )
  if (inherits(fit_value, "error")) {
    reason <- if (inherits(fit_value, "ab_target_filtered")) {
      "target_filtered"
    } else if (inherits(fit_value, "ab_insufficient_features")) {
      "insufficient_retained_features"
    } else if (inherits(fit_value, "ab_nonpositive_library")) {
      "nonpositive_retained_library"
    } else if (inherits(fit_value, "ab_normalization_failure")) {
      "normalization_failure"
    } else {
      "edger_fit_failure"
    }
    return(.ab_edger_failure(
      data, assignment, rebuilt, reason, conditionMessage(fit_value),
      captured_warnings, diagnostics, started
    ))
  }

  keep <- fit_value$keep
  retained <- fit_value$retained
  fit <- fit_value$fit
  test <- fit_value$test
  table_beta <- unname(test$table[data$claim$focal_gene, "logFC"])
  coefficient_loge <- unname(fit$coefficients[data$claim$focal_gene, 2L])
  coefficient_log2 <- coefficient_loge / log(2)
  difference <- table_beta - coefficient_log2
  diagnostics$filter_keep <- data.frame(
    feature_id = names(keep), kept = unname(as.logical(keep)),
    stringsAsFactors = FALSE
  )
  diagnostics$retained_feature_ids <- retained
  diagnostics$retained_feature_count <- length(retained)
  diagnostics$normalization_factors <- setNames(
    as.numeric(fit_value$dge$samples$norm.factors), focal_meta$sample_id
  )
  diagnostics$retained_library_sizes <- setNames(
    as.numeric(fit_value$dge$samples$lib.size), focal_meta$sample_id
  )
  diagnostics$effective_library_sizes <- setNames(
    as.numeric(fit_value$effective), focal_meta$sample_id
  )
  diagnostics$model_coefficient_loge <- coefficient_loge
  diagnostics$model_coefficient_log2 <- coefficient_log2
  diagnostics$test_table_logFC <- table_beta
  diagnostics$coefficient_difference <- difference
  diagnostics$dispersion <- fit$dispersion
  diagnostics$df_prior <- fit$df.prior
  diagnostics$var_prior <- fit$var.prior
  sample_summary$normalization_factor <- diagnostics$normalization_factors
  sample_summary$retained_library_size <- diagnostics$retained_library_sizes
  sample_summary$effective_library_size <- diagnostics$effective_library_sizes
  diagnostics$sample_summary <- sample_summary

  if (any(!is.finite(c(table_beta, coefficient_loge, coefficient_log2, difference)))) {
    return(.ab_edger_failure(
      data, assignment, rebuilt, "nonfinite_result",
      "edgeR returned a nonfinite target effect or coefficient",
      captured_warnings, diagnostics, started
    ))
  }
  if (abs(difference) > data$analysis_spec$endpoint$coefficient_check_tolerance) {
    return(.ab_edger_failure(
      data, assignment, rebuilt, "coefficient_mismatch",
      "glmQLFTest logFC does not match the registered model coefficient",
      captured_warnings, diagnostics, started
    ))
  }
  margin <- data$claim$q * table_beta - data$claim$delta
  if (!is.finite(margin)) {
    return(.ab_edger_failure(
      data, assignment, rebuilt, "nonfinite_result",
      "registered edgeR margin is nonfinite",
      captured_warnings, diagnostics, started
    ))
  }

  out <- list(
    status = "evaluable",
    beta = as.numeric(table_beta),
    margin = as.numeric(margin),
    claim_holds = isTRUE(margin > data$claim$tolerance),
    overturn = isTRUE(margin <= data$claim$tolerance),
    sample_expression = sample_summary,
    per_sample = sample_summary,
    pseudobulk_summary = .ab_long_pseudobulk(data, rebuilt),
    error = NA_character_,
    failure_reason = NA_character_,
    warnings = unique(captured_warnings),
    assignment = assignment,
    assignment_id = .ab_assignment_id(data, assignment),
    pseudobulk = rebuilt,
    analysis_id = data$analysis_spec$analysis_id,
    backend = data$analysis_spec$backend$package,
    backend_version = observed_version,
    evaluation_identity = .ab_authoritative_identity(data),
    fit_elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
    authoritative_rebuild = TRUE
  )
  for (name in names(diagnostics)) out[[name]] <- diagnostics[[name]]
  structure(out, class = "ab_fit")
}
