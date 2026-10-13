.ab_sim_with_local_seed <- function(seed, fun) {
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be one nonnegative finite integer", call. = FALSE)
  }
  if (!is.function(fun)) stop("fun must be a function", call. = FALSE)
  seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_existed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (seed_existed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  fun()
}

.ab_sim_scalar <- function(x, name, type = c("numeric", "integer", "character"),
                           minimum = NULL) {
  type <- match.arg(type)
  value <- x[[name]]
  if (is.null(value) || length(value) != 1L || is.na(value)) {
    stop("simulation family is missing scalar field: ", name, call. = FALSE)
  }
  value <- switch(type,
    numeric = suppressWarnings(as.numeric(value)),
    integer = suppressWarnings(as.integer(value)),
    character = as.character(value)
  )
  if (type != "character" && !is.finite(value)) {
    stop("simulation family field is not finite: ", name, call. = FALSE)
  }
  if (type == "integer" && as.numeric(x[[name]]) != value) {
    stop("simulation family field must be an integer: ", name, call. = FALSE)
  }
  if (!is.null(minimum) && value < minimum) {
    stop("simulation family field is below its minimum: ", name, call. = FALSE)
  }
  value
}

.ab_sim_registered_mechanisms <- function() {
  c(
    stable_true_effects = "symmetric_correct_uncertainty",
    null_condition_dependent_contamination = "case_B_observed_as_A",
    subtle_effect_with_mislabeling = "asymmetric_bidirectional_mislabeling",
    rare_populations = "rare_A_bidirectional",
    inbound_only_vulnerability = "allowed_B_to_A_only_in_case",
    outbound_only_vulnerability = "allowed_A_to_B_only_in_case",
    balanced_swaps = "within_specimen_opposite_pairs",
    one_influential_donor = "concentrated_case_donor",
    miscalibrated_annotation_probabilities = "asymmetric_bidirectional_mislabeling",
    missing_correct_allowed_alternatives = "locked_wrong_plus_visible_uncertainty",
    continuous_transitional_biology = "six_transition_states"
  )
}

.ab_sim_gamma_shape_scale <- function(cv) {
  if (!is.finite(cv) || cv < 0) stop("gamma CV must be nonnegative", call. = FALSE)
  if (cv == 0) return(c(shape = Inf, scale = 0))
  c(shape = 1 / cv^2, scale = cv^2)
}

.ab_sim_gamma_mean_one <- function(n, cv) {
  parameters <- .ab_sim_gamma_shape_scale(cv)
  if (is.infinite(parameters[["shape"]])) return(rep(1, n))
  stats::rgamma(n, shape = parameters[["shape"]], scale = parameters[["scale"]])
}

.ab_sim_interleaved_sample <- function(cells, label, n, condition = NULL,
                                       sample_id = NULL, exclude = character()) {
  keep <- cells$true_label == label & !cells$cell_id %in% exclude
  if (!is.null(condition)) keep <- keep & cells$condition %in% condition
  if (!is.null(sample_id)) keep <- keep & cells$sample_id %in% sample_id
  candidates <- which(keep)
  if (length(candidates) < n) {
    stop("registered confusion mechanism has too few eligible cells", call. = FALSE)
  }
  by_sample <- split(candidates, cells$sample_id[candidates])
  by_sample <- by_sample[lengths(by_sample) > 0L]
  by_sample <- lapply(by_sample, sample)
  result <- integer()
  while (length(result) < n) {
    progressed <- FALSE
    for (name in names(by_sample)) {
      if (length(by_sample[[name]]) && length(result) < n) {
        result <- c(result, by_sample[[name]][[1L]])
        by_sample[[name]] <- by_sample[[name]][-1L]
        progressed <- TRUE
      }
    }
    if (!progressed) break
  }
  result
}

.ab_sim_confusion_plan <- function(cells, family, ambiguous_cells) {
  family_id <- as.character(family$family_id)
  if (ambiguous_cells != 6L) {
    stop("generator v1 requires exactly six registered ambiguous cells", call. = FALSE)
  }
  observed <- cells$true_label
  registered <- rep(FALSE, nrow(cells))
  locked_wrong <- rep(FALSE, nrow(cells))
  true_state <- rep(NA_real_, nrow(cells))
  selected <- integer()
  flip <- function(index) {
    observed[index] <<- ifelse(cells$true_label[index] == "A", "B", "A")
  }
  mark <- function(index, do_flip = FALSE) {
    registered[index] <<- TRUE
    selected <<- c(selected, index)
    if (do_flip) flip(index)
  }

  if (identical(family_id, "stable_true_effects")) {
    mark(c(.ab_sim_interleaved_sample(cells, "A", 3L),
           .ab_sim_interleaved_sample(cells, "B", 3L)))
  } else if (identical(family_id, "null_condition_dependent_contamination")) {
    mark(.ab_sim_interleaved_sample(cells, "B", 6L, condition = "case"), TRUE)
  } else if (identical(family_id, "subtle_effect_with_mislabeling") ||
             identical(family_id, "miscalibrated_annotation_probabilities")) {
    mark(c(.ab_sim_interleaved_sample(cells, "B", 4L, condition = "case"),
           .ab_sim_interleaved_sample(cells, "A", 2L, condition = "case")), TRUE)
  } else if (identical(family_id, "rare_populations")) {
    mark(c(.ab_sim_interleaved_sample(cells, "B", 3L, condition = "case"),
           .ab_sim_interleaved_sample(cells, "A", 3L, condition = "case")), TRUE)
  } else if (identical(family_id, "inbound_only_vulnerability")) {
    mark(.ab_sim_interleaved_sample(cells, "B", 6L, condition = "case"), FALSE)
  } else if (identical(family_id, "outbound_only_vulnerability")) {
    mark(.ab_sim_interleaved_sample(cells, "A", 6L, condition = "case"), FALSE)
  } else if (identical(family_id, "balanced_swaps")) {
    chosen_samples <- sample(cells$sample_id[!duplicated(cells$sample_id)], 3L)
    for (sample_name in chosen_samples) {
      a <- .ab_sim_interleaved_sample(cells, "A", 1L, sample_id = sample_name)
      b <- .ab_sim_interleaved_sample(cells, "B", 1L, sample_id = sample_name)
      mark(c(a, b), TRUE)
    }
  } else if (identical(family_id, "one_influential_donor")) {
    influential_sample <- unique(cells$sample_id[cells$condition == "case"])[[1L]]
    mark(c(.ab_sim_interleaved_sample(cells, "B", 3L,
                                      sample_id = influential_sample),
           .ab_sim_interleaved_sample(cells, "A", 3L,
                                      sample_id = influential_sample)), TRUE)
  } else if (identical(family_id, "missing_correct_allowed_alternatives")) {
    locked_n <- .ab_sim_scalar(family, "locked_wrong_cells", "integer", 1)
    locked_index <- .ab_sim_interleaved_sample(
      cells, "B", locked_n, condition = "case"
    )
    flip(locked_index)
    locked_wrong[locked_index] <- TRUE
    visible_wrong_b <- .ab_sim_interleaved_sample(
      cells, "B", 2L, condition = "case", exclude = cells$cell_id[locked_index]
    )
    visible_wrong_a <- .ab_sim_interleaved_sample(
      cells, "A", 2L, condition = "case", exclude = cells$cell_id[locked_index]
    )
    visible_correct <- .ab_sim_interleaved_sample(
      cells, "A", 2L, condition = "control",
      exclude = cells$cell_id[c(locked_index, visible_wrong_b, visible_wrong_a)]
    )
    mark(c(visible_wrong_b, visible_wrong_a), TRUE)
    mark(visible_correct, FALSE)
  } else if (identical(family_id, "continuous_transitional_biology")) {
    states <- suppressWarnings(as.numeric(unlist(family$transition_states)))
    if (length(states) != 6L || any(!is.finite(states)) ||
        any(states < 0 | states > 1)) {
      stop("continuous family requires six transition states in [0,1]", call. = FALSE)
    }
    below <- which(states < 0.5)
    above <- which(states >= 0.5)
    if (length(below) + length(above) != 6L) {
      stop("invalid registered transition states", call. = FALSE)
    }
    index <- c(.ab_sim_interleaved_sample(cells, "A", length(below)),
               .ab_sim_interleaved_sample(cells, "B", length(above)))
    selected <- index
    registered[index] <- TRUE
    true_state[index] <- c(states[below], states[above])
    discretized <- ifelse(true_state[index] < 0.5, "A", "B")
    if (!identical(unname(discretized), unname(cells$true_label[index]))) {
      stop("transition-cell binary truth must equal the registered discretization",
           call. = FALSE)
    }
  } else {
    stop("unimplemented registered simulation family: ", family_id, call. = FALSE)
  }
  if (sum(registered) != ambiguous_cells) {
    stop("confusion mechanism did not create exactly six ambiguous cells",
         call. = FALSE)
  }
  list(
    observed_label = observed,
    registered_uncertain = registered,
    locked_wrong = locked_wrong,
    true_state = true_state,
    selected_index = selected,
    influential_sample = if (identical(family_id, "one_influential_donor")) {
      influential_sample
    } else NA_character_
  )
}

.ab_sim_probability_table <- function(cells, profile, generator_id) {
  labels <- lapply(seq_len(nrow(cells)), function(index) {
    if (isTRUE(cells$registered_uncertain[[index]])) c("A", "B") else {
      cells$observed_label[[index]]
    }
  })
  rows <- vector("list", sum(lengths(labels)))
  at <- 1L
  uses_hidden_truth <- FALSE
  for (index in seq_len(nrow(cells))) {
    for (label in labels[[index]]) {
      probability <- if (length(labels[[index]]) == 1L) {
        1
      } else if (identical(profile, "stipulated_moderate")) {
        uses_hidden_truth <- TRUE
        if (identical(label, cells$true_label[[index]])) 0.8 else 0.2
      } else if (identical(profile, "overconfident_wrong")) {
        if (identical(label, cells$observed_label[[index]])) 0.95 else 0.05
      } else if (identical(profile, "state_distance")) {
        uses_hidden_truth <- TRUE
        state <- cells$true_state[[index]]
        if (!is.finite(state)) {
          if (identical(label, cells$true_label[[index]])) 0.9 else 0.1
        } else if (identical(label, "A")) 1 - state else state
      } else {
        stop("unsupported registered probability profile: ", profile,
             call. = FALSE)
      }
      rows[[at]] <- data.frame(
        cell_id = cells$cell_id[[index]],
        label = label,
        probability = probability,
        probability_source = if (identical(profile, "overconfident_wrong")) {
          "intentionally_miscalibrated_stipulation"
        } else if (identical(profile, "state_distance")) {
          "oracle_informed_state_stipulation"
        } else {
          "oracle_informed_label_stipulation"
        },
        probability_version = paste0(generator_id, ":", profile),
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  out <- do.call(rbind, rows)
  totals <- tapply(out$probability, out$cell_id, sum)
  if (any(abs(totals - 1) > 1e-12)) {
    stop("stipulated annotation probabilities are not normalized", call. = FALSE)
  }
  attr(out, "uses_hidden_truth") <- uses_hidden_truth
  out
}

.ab_sim_generating_category <- function(family_id, effect) {
  if (identical(family_id, "null_condition_dependent_contamination")) {
    return("true_null_with_condition_dependent_contamination")
  }
  if (identical(family_id, "stable_true_effects")) {
    return("strong_positive_true_effect")
  }
  if (identical(family_id, "subtle_effect_with_mislabeling")) {
    return("subtle_positive_true_effect")
  }
  if (identical(family_id, "continuous_transitional_biology")) {
    return("continuous_transition_with_binary_analysis_discretization")
  }
  if (effect > 0) "positive_true_effect" else if (effect < 0) {
    "negative_true_effect"
  } else "true_null"
}

.ab_sim_compact_fit <- function(fit) {
  list(
    status = fit$status,
    beta = fit$beta,
    margin = fit$margin,
    claim_holds = fit$claim_holds,
    overturn = fit$overturn,
    failure_reason = fit$failure_reason,
    warnings = fit$warnings,
    error = fit$error,
    fit_elapsed_seconds = fit$fit_elapsed_seconds
  )
}

.ab_sim_make_input <- function(counts, cells, specimens, allowed_labels,
                               claim, constraints, analysis_spec) {
  list(
    counts = counts,
    cells = cells,
    specimens = specimens,
    allowed_labels = allowed_labels,
    claim = claim,
    constraints = constraints,
    analysis_spec = analysis_spec
  )
}

.ab_simulate_case <- function(family, seed, generator_spec, analysis_spec,
                              common = NULL) {
  if (!is.list(family) || is.data.frame(family)) {
    stop("family must be one parsed manifest family mapping", call. = FALSE)
  }
  if (!is.list(generator_spec) ||
      !identical(as.character(generator_spec$generator_id),
                 "gamma-nb-independent-donor-v1")) {
    stop("generator_spec must be gamma-nb-independent-donor-v1", call. = FALSE)
  }
  if (is.null(common)) {
    common <- list(
      target_gene = "G_target", focal_label = "A", alternative_label = "B",
      q = 1, delta = 0.5, tolerance = 1e-12, background_features = 100,
      ambiguous_cells = 6, max_changes_global = 3,
      max_changes_per_sample = 2, max_changes_per_donor = 2,
      min_cells_per_sample_label = 2, min_library_per_sample_label = 50,
      min_donors_per_condition = 2
    )
  }
  if (!is.list(common) || is.data.frame(common)) {
    stop("common must be a parsed manifest common mapping", call. = FALSE)
  }
  registered_value <- function(name, type = c("numeric", "integer", "character"),
                               minimum = NULL) {
    source <- if (!is.null(family[[name]])) family else common
    .ab_sim_scalar(source, name, match.arg(type), minimum)
  }
  normalized_analysis <- .ab_normalize_analysis_spec(analysis_spec)
  if (!isTRUE(normalized_analysis$ok)) {
    stop(paste(normalized_analysis$errors, collapse = "; "), call. = FALSE)
  }
  mechanisms <- .ab_sim_registered_mechanisms()
  family_id <- .ab_sim_scalar(family, "family_id", "character")
  if (!family_id %in% names(mechanisms)) {
    stop("family_id is not registered in generator v1: ", family_id,
         call. = FALSE)
  }
  confusion <- .ab_sim_scalar(family, "confusion", "character")
  if (!identical(unname(mechanisms[[family_id]]), confusion)) {
    stop("family confusion does not match its registered v1 mechanism",
         call. = FALSE)
  }

  .ab_sim_with_local_seed(seed, function() {
    counts_spec <- generator_spec$counts
    donors_per_condition <- .ab_sim_scalar(
      family, "donors_per_condition", "integer", 2
    )
    n_a <- .ab_sim_scalar(family, "cells_A_per_sample", "integer", 1)
    n_b <- .ab_sim_scalar(family, "cells_B_per_sample", "integer", 1)
    effect_a <- .ab_sim_scalar(family, "log2_effect_A", "numeric")
    effect_b <- .ab_sim_scalar(family, "log2_effect_B", "numeric")
    target_a <- .ab_sim_scalar(family, "target_mean_A_control", "numeric", 0)
    target_b <- .ab_sim_scalar(family, "target_mean_B_control", "numeric", 0)
    dispersion <- .ab_sim_scalar(family, "dispersion", "numeric", 0)
    donor_cv <- .ab_sim_scalar(family, "donor_gamma_cv", "numeric", 0)
    cell_cv <- suppressWarnings(as.numeric(counts_spec$cell_size_gamma_cv))
    feature_count <- as.integer(counts_spec$feature_count)
    target_gene <- as.character(counts_spec$target_gene)
    background_count <- feature_count - 1L
    if (feature_count < 50L || background_count < 1L || dispersion <= 0 ||
        target_a <= 0 || target_b <= 0 || !is.finite(cell_cv) || cell_cv < 0) {
      stop("invalid registered count-generation parameters", call. = FALSE)
    }
    if (identical(family_id, "balanced_swaps") && n_a != n_b) {
      stop("balanced_swaps requires equal registered A and B counts", call. = FALSE)
    }

    specimens <- do.call(rbind, lapply(c("control", "case"), function(condition) {
      data.frame(
        sample_id = sprintf("S_%s_%02d", condition,
                            seq_len(donors_per_condition)),
        donor_id = sprintf("D_%s_%02d", condition,
                           seq_len(donors_per_condition)),
        condition = condition,
        stringsAsFactors = FALSE
      )
    }))
    rownames(specimens) <- NULL
    skeleton <- do.call(rbind, lapply(seq_len(nrow(specimens)), function(index) {
      rbind(
        data.frame(sample_id = specimens$sample_id[[index]],
                   donor_id = specimens$donor_id[[index]],
                   condition = specimens$condition[[index]],
                   true_label = "A", stringsAsFactors = FALSE)[rep(1L, n_a), ],
        data.frame(sample_id = specimens$sample_id[[index]],
                   donor_id = specimens$donor_id[[index]],
                   condition = specimens$condition[[index]],
                   true_label = "B", stringsAsFactors = FALSE)[rep(1L, n_b), ]
      )
    }))
    rownames(skeleton) <- NULL
    # The random permutation consumes no hidden cell state and the resulting token
    # contains no label, state, family, condition, or corruption marker.
    skeleton$cell_id <- sample(sprintf("C%07d", seq_len(nrow(skeleton))),
                               nrow(skeleton), replace = FALSE)
    opaque_id_ok <- !any(grepl("A|B|control|case|transition|wrong", skeleton$cell_id,
                              ignore.case = TRUE)) &&
      !any(grepl(family_id, skeleton$cell_id, fixed = TRUE)) &&
      !anyDuplicated(skeleton$cell_id)
    if (!isTRUE(opaque_id_ok)) {
      stop("internal error: optimizer-facing cell IDs encode hidden state",
           call. = FALSE)
    }
    cells_hidden <- skeleton[, c(
      "cell_id", "sample_id", "donor_id", "condition", "true_label"
    )]
    ambiguous_count <- registered_value("ambiguous_cells", "integer", 1)
    plan <- .ab_sim_confusion_plan(cells_hidden, family, ambiguous_count)
    cells_hidden$observed_label <- plan$observed_label
    cells_hidden$registered_uncertain <- plan$registered_uncertain
    cells_hidden$locked_wrong <- plan$locked_wrong
    cells_hidden$true_state <- plan$true_state

    profile <- .ab_sim_scalar(family, "probability_profile", "character")
    if (!profile %in% c("stipulated_moderate", "overconfident_wrong",
                        "state_distance")) {
      stop("unregistered probability profile", call. = FALSE)
    }
    probabilities <- .ab_sim_probability_table(
      cells_hidden, profile, as.character(generator_spec$generator_id)
    )

    features <- c(target_gene, sprintf("%s%03d",
      as.character(counts_spec$background_gene_prefix), seq_len(background_count)
    ))
    base_background <- stats::rgamma(
      background_count,
      shape = as.numeric(counts_spec$background_mean_shape),
      rate = as.numeric(counts_spec$background_mean_rate)
    )
    type_a_background <- base_background * .ab_sim_gamma_mean_one(background_count, 0.15)
    type_b_background <- base_background * .ab_sim_gamma_mean_one(background_count, 0.15)
    donor_effect <- matrix(
      .ab_sim_gamma_mean_one(length(features) * nrow(specimens), donor_cv),
      nrow = length(features), ncol = nrow(specimens),
      dimnames = list(features, specimens$donor_id)
    )
    if (identical(family_id, "one_influential_donor")) {
      multiplier <- .ab_sim_scalar(
        family, "influential_donor_multiplier", "numeric", 0
      )
      influential_donor <- specimens$donor_id[
        match(plan$influential_sample, specimens$sample_id)
      ]
      donor_effect[target_gene, influential_donor] <-
        donor_effect[target_gene, influential_donor] * multiplier
    } else {
      influential_donor <- NA_character_
      multiplier <- NA_real_
    }
    cell_size <- .ab_sim_gamma_mean_one(nrow(cells_hidden), cell_cv)
    names(cell_size) <- cells_hidden$cell_id
    minimum_mean <- as.numeric(counts_spec$minimum_mean)
    mean_matrix <- matrix(NA_real_, nrow = length(features),
                          ncol = nrow(cells_hidden),
                          dimnames = list(features, cells_hidden$cell_id))
    for (index in seq_len(nrow(cells_hidden))) {
      condition_case <- cells_hidden$condition[[index]] == "case"
      state <- cells_hidden$true_state[[index]]
      if (!is.finite(state)) state <- if (cells_hidden$true_label[[index]] == "A") 0 else 1
      mean_a <- c(target_a * 2^(effect_a * condition_case), type_a_background)
      mean_b <- c(target_b * 2^(effect_b * condition_case), type_b_background)
      biological_mean <- (1 - state) * mean_a + state * mean_b
      donor <- cells_hidden$donor_id[[index]]
      mean_matrix[, index] <- pmax(
        minimum_mean,
        biological_mean * donor_effect[, donor] * cell_size[[index]]
      )
    }
    raw_counts <- matrix(
      stats::rnbinom(length(mean_matrix), mu = as.vector(mean_matrix),
                     size = 1 / dispersion),
      nrow = nrow(mean_matrix), ncol = ncol(mean_matrix),
      dimnames = dimnames(mean_matrix)
    )
    if (any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
        any(raw_counts != floor(raw_counts))) {
      stop("negative-binomial generator produced invalid counts", call. = FALSE)
    }
    counts <- methods::as(Matrix::Matrix(raw_counts, sparse = TRUE), "dgCMatrix")

    # Export order is determined only by the opaque token. Because token assignment
    # was independently permuted, neither row order nor count-column order retains
    # the generator's within-sample A-then-B construction blocks.
    export_order <- order(cells_hidden$cell_id, method = "radix")
    cells_hidden <- cells_hidden[export_order, , drop = FALSE]
    rownames(cells_hidden) <- NULL
    counts <- counts[, cells_hidden$cell_id, drop = FALSE]
    probabilities <- probabilities[
      order(probabilities$cell_id, probabilities$label, method = "radix"),
      , drop = FALSE
    ]
    rownames(probabilities) <- NULL
    export_order_ok <- identical(cells_hidden$cell_id, colnames(counts)) &&
      identical(cells_hidden$cell_id,
                sort(cells_hidden$cell_id, method = "radix")) &&
      identical(unique(probabilities$cell_id), cells_hidden$cell_id)
    if (!isTRUE(export_order_ok)) {
      stop("internal error: export-facing cell order is not opaque-ID based",
           call. = FALSE)
    }

    search_cells <- cells_hidden[, c(
      "cell_id", "sample_id", "donor_id", "condition", "observed_label"
    )]
    names(search_cells)[names(search_cells) == "observed_label"] <- "original_label"
    allowed_labels <- do.call(rbind, lapply(seq_len(nrow(cells_hidden)), function(index) {
      allowed <- probabilities$label[probabilities$cell_id == cells_hidden$cell_id[[index]]]
      data.frame(
        cell_id = cells_hidden$cell_id[[index]],
        allowed_label = allowed,
        evidence_source = "simulation_stipulated_uncertainty",
        evidence_version = paste0(generator_spec$generator_id, ":", profile),
        stringsAsFactors = FALSE
      )
    }))
    common_target <- registered_value("target_gene", "character")
    focal_label <- registered_value("focal_label", "character")
    alternative_label <- registered_value("alternative_label", "character")
    background_features <- registered_value("background_features", "integer", 1)
    if (!identical(common_target, target_gene) || background_features != background_count ||
        !identical(c(focal_label, alternative_label), c("A", "B"))) {
      stop("manifest common fields do not match generator v1 labels/features",
           call. = FALSE)
    }
    claim <- data.frame(
      focal_gene = common_target,
      focal_label = focal_label,
      alternative_label = alternative_label,
      case_condition = "case",
      control_condition = "control",
      q = registered_value("q", "integer"),
      delta = registered_value("delta", "numeric", 0),
      tolerance = registered_value("tolerance", "numeric", 0),
      feature_universe = paste(features, collapse = "|"),
      estimator_id = normalized_analysis$spec$analysis_id,
      stringsAsFactors = FALSE
    )
    min_cells <- registered_value("min_cells_per_sample_label", "integer", 1)
    minimum_library <- registered_value(
      "min_library_per_sample_label", "integer", 1
    )
    constraints <- data.frame(
      min_cells_A = min_cells,
      min_cells_B = min_cells,
      min_library_A = minimum_library,
      min_library_B = minimum_library,
      max_changes_global = registered_value("max_changes_global", "integer", 0),
      max_changes_per_sample = registered_value(
        "max_changes_per_sample", "integer", 0
      ),
      max_changes_per_donor = registered_value(
        "max_changes_per_donor", "integer", 0
      ),
      min_donors_per_condition = registered_value(
        "min_donors_per_condition", "integer", 2
      ),
      stringsAsFactors = FALSE
    )
    search_input <- .ab_sim_make_input(
      counts, search_cells, specimens, allowed_labels, claim, constraints,
      normalized_analysis$spec
    )
    expected_search_names <- c(
      "counts", "cells", "specimens", "allowed_labels", "claim", "constraints",
      "analysis_spec"
    )
    if (!identical(names(search_input), expected_search_names)) {
      stop("internal error: search input is not strictly allowlisted", call. = FALSE)
    }
    prohibited <- unlist(generator_spec$truth_separation$prohibited_search_fields)
    search_fields <- unlist(lapply(search_input[c(
      "cells", "specimens", "allowed_labels", "claim", "constraints"
    )], names), use.names = FALSE)
    if (any(prohibited %in% search_fields)) {
      stop("internal error: hidden truth leaked into search tables", call. = FALSE)
    }

    true_cells <- search_cells
    true_cells$original_label <- cells_hidden$true_label
    true_allowed <- data.frame(
      cell_id = true_cells$cell_id,
      allowed_label = true_cells$original_label,
      evidence_source = "simulation_evaluation_truth",
      evidence_version = as.character(generator_spec$generator_id),
      stringsAsFactors = FALSE
    )
    true_input <- .ab_sim_make_input(
      counts, true_cells, specimens, true_allowed, claim, constraints,
      normalized_analysis$spec
    )
    true_validation <- validate_inputs(
      true_input$counts, true_input$cells, true_input$specimens,
      true_input$allowed_labels, true_input$claim, true_input$constraints,
      true_input$analysis_spec
    )
    true_fit <- if (isTRUE(true_validation$ok)) {
      fit_edger_claim(true_validation)
    } else {
      structure(list(
        status = "not_assessable", beta = NA_real_, margin = NA_real_,
        claim_holds = NA, overturn = FALSE,
        failure_reason = "true_label_input_invalid",
        warnings = true_validation$warnings,
        error = paste(true_validation$errors, collapse = "; "),
        fit_elapsed_seconds = NA_real_
      ), class = "ab_fit")
    }

    cell_truth <- cells_hidden
    cell_truth$is_mislabeled <- cells_hidden$observed_label != cells_hidden$true_label
    cell_truth$correct_alternative <- ifelse(
      cell_truth$is_mislabeled, cell_truth$true_label, NA_character_
    )
    cell_truth$correct_alternative_available <- vapply(seq_len(nrow(cell_truth)),
      function(index) {
        any(allowed_labels$cell_id == cell_truth$cell_id[[index]] &
              allowed_labels$allowed_label == cell_truth$true_label[[index]])
      }, logical(1L))
    cell_truth$binary_label_interpretation <- if (
      identical(family_id, "continuous_transitional_biology")
    ) "registered_discretization_not_annotation_error" else "discrete_type"

    donor_effect_table <- data.frame(
      feature_id = rep(rownames(donor_effect), times = ncol(donor_effect)),
      donor_id = rep(colnames(donor_effect), each = nrow(donor_effect)),
      gamma_effect = as.vector(donor_effect),
      stringsAsFactors = FALSE
    )
    cell_size_table <- data.frame(
      cell_id = names(cell_size), gamma_size_factor = unname(cell_size),
      stringsAsFactors = FALSE
    )
    generating_effects <- data.frame(
      label = c("A", "B"),
      target_control_mean = c(target_a, target_b),
      target_log2_condition_effect = c(effect_a, effect_b),
      stringsAsFactors = FALSE
    )
    category <- .ab_sim_generating_category(family_id, effect_a)
    uses_truth <- isTRUE(attr(probabilities, "uses_hidden_truth"))
    generation_record <- list(
      schema_version = "1.0.0",
      generator_id = as.character(generator_spec$generator_id),
      family_id = family_id,
      seed = as.integer(seed),
      confusion_mechanism = confusion,
      probability_profile = profile,
      probability_source = unique(probabilities$probability_source),
      probability_construction_uses_hidden_generating_values = uses_truth,
      probability_calibration_status = if (identical(profile, "overconfident_wrong")) {
        "intentionally_miscalibrated"
      } else if (uses_truth) {
        "oracle_informed_stipulated_not_empirically_calibrated"
      } else {
        "stipulated_not_empirically_calibrated"
      },
      probabilities_frozen_before_search = TRUE,
      optimizer_cell_id_rule = "random permutation of opaque numeric tokens; hidden state is not an input",
      optimizer_cell_ids_state_independent = opaque_id_ok,
      export_order_uses_opaque_cell_id_only = export_order_ok,
      independent_donors = TRUE,
      one_specimen_per_donor = TRUE,
      donor_feature_gamma_effects_independent_draws = TRUE,
      donor_gamma_cv = donor_cv,
      cell_size_gamma_cv = cell_cv,
      negative_binomial_dispersion = dispersion,
      feature_count = length(features),
      cell_count = nrow(cells_hidden),
      registered_ambiguous_cell_count = sum(cells_hidden$registered_uncertain),
      mislabeled_cell_count = sum(cell_truth$is_mislabeled),
      locked_wrong_cell_count = sum(cells_hidden$locked_wrong),
      missing_correct_alternative_count = sum(
        cell_truth$is_mislabeled & !cell_truth$correct_alternative_available
      ),
      balanced_swap_invariant = if (identical(family_id, "balanced_swaps")) {
        all(table(search_cells$sample_id, search_cells$original_label)[, c("A", "B")] == n_a)
      } else NA,
      influential_donor = influential_donor,
      influential_donor_multiplier = multiplier,
      generating_effect_category = category,
      limitation = as.character(family$limitation),
      count_diagnostics = list(
        all_integer_nonnegative = all(raw_counts >= 0 & raw_counts == floor(raw_counts)),
        zero_fraction = mean(raw_counts == 0),
        minimum_library_size = min(colSums(raw_counts)),
        maximum_library_size = max(colSums(raw_counts))
      ),
      latent_factor_diagnostics = list(
        donor_effect_min = min(donor_effect),
        donor_effect_max = max(donor_effect),
        donor_effect_mean = mean(donor_effect),
        cell_size_min = min(cell_size),
        cell_size_max = max(cell_size),
        cell_size_mean = mean(cell_size)
      )
    )
    if (identical(family_id, "balanced_swaps") &&
        (!isTRUE(generation_record$balanced_swap_invariant) ||
         constraints$min_cells_A != n_a || constraints$min_cells_B != n_b)) {
      stop("balanced swap invariant or per-label minima were not enforced",
           call. = FALSE)
    }

    truth <- list(
      schema_version = "1.0.0",
      family_id = family_id,
      seed = as.integer(seed),
      generating_effect_category = category,
      generating_log2_effect_A = effect_a,
      probability_construction = list(
        profile = profile,
        uses_hidden_generating_values = uses_truth,
        calibration_status = generation_record$probability_calibration_status,
        frozen_before_search = TRUE
      ),
      generating_effects = generating_effects,
      cell_truth = cell_truth,
      donor_feature_effects = donor_effect_table,
      cell_size_factors = cell_size_table,
      true_label_input = true_input,
      true_label_validation = list(
        ok = true_validation$ok,
        status = true_validation$status,
        errors = true_validation$errors,
        warnings = true_validation$warnings
      ),
      true_label_fit = true_fit,
      true_label_fit_summary = .ab_sim_compact_fit(true_fit)
    )
    structure(list(
      search_input = search_input,
      probabilities = probabilities,
      truth = truth,
      generation_record = generation_record
    ), class = c("ab_simulation_case", "list"))
  })
}

.ab_sim_write_table <- function(value, path) {
  utils::write.table(
    value, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", eol = "\n", fileEncoding = "UTF-8"
  )
}

.ab_sim_write_input <- function(input, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  counts_table <- data.frame(
    gene_id = rownames(input$counts), as.matrix(input$counts),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  .ab_sim_write_table(counts_table, file.path(directory, "counts.tsv"))
  for (name in c("cells", "specimens", "allowed_labels", "claim", "constraints")) {
    .ab_sim_write_table(input[[name]], file.path(directory, paste0(name, ".tsv")))
  }
  yaml::write_yaml(.ab_yaml_ready(input$analysis_spec),
                   file.path(directory, "analysis_spec.yaml"))
  invisible(directory)
}

.ab_write_simulation_case <- function(simulation, output_dir) {
  if (!inherits(simulation, "ab_simulation_case")) {
    stop("simulation must be an ab_simulation_case", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path", call. = FALSE)
  }
  required_search <- c(
    "counts", "cells", "specimens", "allowed_labels", "claim", "constraints",
    "analysis_spec"
  )
  if (!identical(names(simulation$search_input), required_search)) {
    stop("simulation search_input is not strictly allowlisted", call. = FALSE)
  }
  root <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  search_dir <- file.path(root, "search_input")
  method_dir <- file.path(root, "method_input")
  truth_dir <- file.path(root, "evaluation_truth")
  true_input_dir <- file.path(truth_dir, "true_label_input")
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)
  .ab_sim_write_input(simulation$search_input, search_dir)
  .ab_sim_write_table(simulation$probabilities,
                      file.path(method_dir, "annotation_probabilities.tsv"))
  .ab_write_json(list(
    schema_version = "1.0.0",
    family_id = simulation$generation_record$family_id,
    seed = simulation$generation_record$seed,
    generator_id = simulation$generation_record$generator_id,
    probability_profile = simulation$generation_record$probability_profile,
    probabilities_frozen_before_search = TRUE
  ), file.path(method_dir, "case_metadata.json"))
  .ab_sim_write_table(simulation$truth$cell_truth,
                      file.path(truth_dir, "cell_truth.tsv"))
  .ab_sim_write_table(simulation$truth$generating_effects,
                      file.path(truth_dir, "generating_effects.tsv"))
  .ab_sim_write_table(simulation$truth$donor_feature_effects,
                      file.path(truth_dir, "donor_feature_effects.tsv"))
  .ab_sim_write_table(simulation$truth$cell_size_factors,
                      file.path(truth_dir, "cell_size_factors.tsv"))
  .ab_write_json(simulation$truth$true_label_fit_summary,
                 file.path(truth_dir, "true_label_fit.json"))
  .ab_write_json(list(
    schema_version = simulation$truth$schema_version,
    family_id = simulation$truth$family_id,
    seed = simulation$truth$seed,
    generating_effect_category = simulation$truth$generating_effect_category
  ), file.path(truth_dir, "generating_truth.json"))
  .ab_write_json(simulation$generation_record,
                 file.path(truth_dir, "generation_record.json"))
  .ab_sim_write_input(simulation$truth$true_label_input, true_input_dir)

  paths <- list(
    search_input = search_dir,
    probabilities = file.path(method_dir, "annotation_probabilities.tsv"),
    method_metadata = file.path(method_dir, "case_metadata.json"),
    truth = truth_dir,
    true_label_input = true_input_dir
  )
  files <- sort(list.files(root, recursive = TRUE, full.names = TRUE),
                method = "radix")
  hashes <- setNames(vapply(files, .ab_file_hash, character(1L)),
                     substring(files, nchar(root) + 2L))
  structure(list(root = root, paths = paths, hashes = hashes),
            class = c("ab_simulation_artifacts", "list"))
}
