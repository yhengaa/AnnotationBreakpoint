.ab_scaling_with_local_seed <- function(seed, fun) {
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be one nonnegative finite integer", call. = FALSE)
  }
  if (!is.function(fun)) stop("fun must be a function", call. = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  fun()
}

.ab_scaling_scalar <- function(x, name, type = c("numeric", "integer", "character"),
                               minimum = NULL) {
  type <- match.arg(type)
  value <- x[[name]]
  if (is.null(value) || length(value) != 1L || is.na(value)) {
    stop("missing scaling field: ", name, call. = FALSE)
  }
  parsed <- switch(type,
    numeric = suppressWarnings(as.numeric(value)),
    integer = suppressWarnings(as.integer(value)),
    character = as.character(value)
  )
  if (type != "character" && (!is.finite(parsed) ||
      (type == "integer" && as.numeric(value) != parsed))) {
    stop("invalid scaling field: ", name, call. = FALSE)
  }
  if (!is.null(minimum) && parsed < minimum) {
    stop("scaling field is below its minimum: ", name, call. = FALSE)
  }
  parsed
}

.ab_scaling_gamma_one <- function(n, cv) {
  if (!is.finite(cv) || cv < 0) stop("gamma CV must be nonnegative", call. = FALSE)
  if (cv == 0) return(rep(1, n))
  stats::rgamma(n, shape = 1 / cv^2, scale = cv^2)
}

.ab_scaling_config <- function(project_root = ".") {
  paths <- file.path(normalizePath(project_root, mustWork = TRUE), "configs", "scaling",
                     c("generator_v1.yaml", "development_manifest_v1.yaml",
                       "search_grid_v1.yaml"))
  analysis_path <- file.path(normalizePath(project_root, mustWork = TRUE),
                             "configs", "m2_analysis_spec.yaml")
  if (any(!file.exists(c(paths, analysis_path)))) {
    stop("one or more registered scaling configuration files are missing", call. = FALSE)
  }
  out <- list(
    generator = yaml::read_yaml(paths[[1L]]),
    manifest = yaml::read_yaml(paths[[2L]]),
    search = yaml::read_yaml(paths[[3L]]),
    analysis = yaml::read_yaml(analysis_path),
    paths = c(generator = paths[[1L]], manifest = paths[[2L]],
              search = paths[[3L]], analysis = analysis_path),
    hashes = setNames(vapply(c(paths, analysis_path), function(path) {
      digest::digest(file = path, algo = "sha256", serialize = FALSE)
    }, character(1L)), c("generator", "manifest", "search", "analysis"))
  )
  parent_seeds <- as.integer(unlist(out$manifest$parent_seeds))
  offsets <- as.integer(unlist(out$generator$seed_stream_offsets))
  role_names <- names(out$generator$seed_stream_offsets)
  stream_table <- expand.grid(parent_seed = parent_seeds, role = role_names,
                              stringsAsFactors = FALSE)
  stream_table$stream_seed <- stream_table$parent_seed +
    offsets[match(stream_table$role, role_names)]
  if (anyDuplicated(stream_table$stream_seed)) {
    collision <- stream_table$stream_seed[duplicated(stream_table$stream_seed) |
      duplicated(stream_table$stream_seed, fromLast = TRUE)]
    stop("registered RNG stream collision: ", paste(sort(unique(collision)), collapse = ","),
         call. = FALSE)
  }
  out$rng_stream_table <- stream_table
  out
}

.ab_scaling_evaluation_config <- function(project_root = ".") {
  root <- normalizePath(project_root, mustWork = TRUE)
  paths <- c(
    generator = file.path(root, "configs", "scaling", "heldout_generator_v1.yaml"),
    manifest = file.path(root, "configs", "scaling",
                         "heldout_evaluation_manifest_draft_v1.yaml"),
    analysis = file.path(root, "configs", "m2_analysis_spec.yaml")
  )
  if (any(!file.exists(paths))) {
    stop("one or more held-out scaling registration files are missing", call. = FALSE)
  }
  out <- list(
    generator = yaml::read_yaml(paths[["generator"]]),
    manifest = yaml::read_yaml(paths[["manifest"]]),
    analysis = yaml::read_yaml(paths[["analysis"]]),
    paths = paths,
    hashes = setNames(vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE
    ), character(1L)), names(paths))
  )
  seeds <- vapply(out$manifest$parents, function(x) as.integer(x$seed), integer(1L))
  offsets <- as.integer(unlist(out$generator$seed_stream_offsets))
  roles <- names(out$generator$seed_stream_offsets)
  streams <- expand.grid(parent_seed = seeds, role = roles,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  streams$stream_seed <- streams$parent_seed + offsets[match(streams$role, roles)]
  if (length(seeds) != 24L || anyDuplicated(seeds) ||
      anyDuplicated(streams$stream_seed)) {
    stop("held-out parent seeds or parent-by-role RNG streams are not unique",
         call. = FALSE)
  }
  out$rng_stream_table <- streams
  out
}

.ab_scaling_evaluation_parent_spec <- function(config, parent) {
  parents <- config$manifest$parents
  if (is.numeric(parent) && length(parent) == 1L && is.finite(parent)) {
    index <- as.integer(parent)
  } else {
    ids <- vapply(parents, function(x) as.character(x$parent_id), character(1L))
    index <- match(as.character(parent), ids)
  }
  if (is.na(index) || index < 1L || index > length(parents)) {
    stop("unknown registered held-out parent", call. = FALSE)
  }
  registered <- parents[[index]]
  recipe_id <- as.character(registered$recipe_id)
  recipe <- config$generator$recipes[[recipe_id]]
  if (is.null(recipe)) stop("held-out parent references an unknown recipe", call. = FALSE)
  out <- utils::modifyList(recipe, registered)
  out$scenario <- recipe_id
  out$stratum <- as.character(recipe$parent_category)
  out
}

.ab_scaling_evaluation_registry <- function(config) {
  parents <- lapply(seq_len(length(config$manifest$parents)), function(index) {
    .ab_scaling_evaluation_parent_spec(config, index)
  })
  parent_table <- do.call(rbind, lapply(seq_along(parents), function(index) {
    x <- parents[[index]]
    data.frame(
      parent_index = index, parent_id = as.character(x$parent_id),
      seed = as.integer(x$seed), parent_category = as.character(x$parent_category),
      recipe_id = as.character(x$recipe_id), stringsAsFactors = FALSE
    )
  }))
  exact <- expand.grid(
    parent_index = seq_len(nrow(parent_table)), k = c(6L, 12L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  exact$method <- "exact"
  exact$max_fits <- NA_integer_
  exact$tier <- "exact_reference"
  exact$task_id <- seq.int(0L, nrow(exact) - 1L)
  practical <- expand.grid(
    parent_index = seq_len(nrow(parent_table)),
    k = c(6L, 12L, 24L, 48L, 96L, 192L),
    method = c("greedy", "random", "pair"), max_fits = c(32L, 128L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  practical$task_id <- seq.int(nrow(exact), nrow(exact) + nrow(practical) - 1L)
  practical$tier <- "primary_equal_budget"
  high_cap <- expand.grid(
    parent_index = seq_len(nrow(parent_table)), k = 192L,
    method = c("greedy", "random", "pair"), max_fits = 1024L,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  high_cap$tier <- "coverage_resource_diagnostic"
  high_cap$task_id <- seq.int(
    nrow(exact) + nrow(practical),
    nrow(exact) + nrow(practical) + nrow(high_cap) - 1L
  )
  tasks <- rbind(exact[, names(practical)], practical, high_cap[, names(practical)])
  tasks <- merge(tasks, parent_table, by = "parent_index", sort = FALSE)
  tasks <- tasks[order(tasks$task_id), , drop = FALSE]
  rownames(tasks) <- NULL
  tasks$run_id <- ifelse(
    tasks$method == "exact",
    paste(tasks$parent_id, paste0("K", tasks$k), "exact", sep = "__"),
    paste(tasks$parent_id, paste0("K", tasks$k), tasks$method,
          paste0("fits", tasks$max_fits), sep = "__")
  )
  recipe_counts <- table(parent_table$recipe_id)
  valid <- nrow(parent_table) == 24L && nrow(exact) == 48L &&
    nrow(practical) == 864L && nrow(high_cap) == 72L && nrow(tasks) == 984L &&
    !anyDuplicated(tasks$task_id) && !anyDuplicated(tasks$run_id) &&
    length(recipe_counts) == 6L && all(recipe_counts == 4L) &&
    identical(tasks$task_id, 0:983)
  if (!valid) stop("held-out scaling task registry invariant failed", call. = FALSE)
  list(parents = parents, parent_table = parent_table, tasks = tasks,
       rng_stream_table = config$rng_stream_table)
}

.ab_scaling_parent_spec <- function(config, parent) {
  parents <- config$manifest$parents
  if (is.numeric(parent) && length(parent) == 1L && is.finite(parent)) {
    index <- as.integer(parent)
  } else {
    ids <- vapply(parents, function(x) as.character(x$parent_id), character(1L))
    index <- match(as.character(parent), ids)
  }
  if (is.na(index) || index < 1L || index > length(parents)) {
    stop("unknown registered development parent", call. = FALSE)
  }
  parents[[index]]
}

.ab_scaling_rank_cells <- function(cells, seed, k_max) {
  if (k_max != 192L) stop("generator v1 requires maximum K=192", call. = FALSE)
  first_plan <- data.frame(
    rank = 1:6,
    donor_id = c("D_control_01", "D_case_01", "D_control_02", "D_case_02",
                 "D_control_03", "D_case_03"),
    original_label = c("B", "A", "B", "A", "B", "A"),
    stringsAsFactors = FALSE
  )
  random_key <- .ab_scaling_with_local_seed(seed, function() stats::runif(nrow(cells)))
  working <- cells[, c("cell_id", "sample_id", "donor_id", "condition",
                       "original_label")]
  working$random_key <- random_key
  first <- do.call(rbind, lapply(seq_len(nrow(first_plan)), function(index) {
    eligible <- which(working$donor_id == first_plan$donor_id[[index]] &
                        working$original_label == first_plan$original_label[[index]])
    eligible <- eligible[order(working$random_key[eligible], method = "radix")]
    if (!length(eligible)) stop("cannot construct registered K6 prefix", call. = FALSE)
    working[eligible[[1L]], , drop = FALSE]
  }))
  if (anyDuplicated(first$cell_id)) stop("K6 prefix is not unique", call. = FALSE)

  # Six cells per donor and label gives 192 cells and exact donor/direction balance.
  strata <- split(seq_len(nrow(working)),
                  interaction(working$donor_id, working$original_label, drop = TRUE))
  selected <- lapply(strata, function(index) {
    fixed <- index[working$cell_id[index] %in% first$cell_id]
    other <- setdiff(index, fixed)
    other <- other[order(working$random_key[other], method = "radix")]
    c(fixed, utils::head(other, 6L - length(fixed)))
  })
  if (any(lengths(selected) != 6L)) {
    stop("cannot construct donor/direction-balanced K192 domain", call. = FALSE)
  }
  selected_ids <- unname(unlist(lapply(selected, function(index) working$cell_id[index])))
  remaining_ids <- setdiff(selected_ids, first$cell_id)
  remaining <- working[match(remaining_ids, working$cell_id), , drop = FALSE]
  # A seeded key orders the remainder without using expression counts or planted status.
  remaining <- remaining[order(remaining$random_key, remaining$cell_id,
                               method = "radix"), , drop = FALSE]
  ranked <- rbind(first, remaining)
  ranked$ambiguity_rank <- seq_len(nrow(ranked))
  ranked$direction <- ifelse(ranked$original_label == "A", "A_to_B", "B_to_A")
  ranked <- ranked[, c("ambiguity_rank", "cell_id", "sample_id", "donor_id",
                       "condition", "original_label", "direction")]
  rownames(ranked) <- NULL
  if (nrow(ranked) != k_max || anyDuplicated(ranked$cell_id) ||
      !identical(ranked$direction[1:6],
                 c("B_to_A", "A_to_B", "B_to_A", "A_to_B", "B_to_A", "A_to_B"))) {
    stop("registered ambiguity ranking invariant failed", call. = FALSE)
  }
  balance <- table(ranked$donor_id, ranked$original_label)
  if (!all(balance[, c("A", "B"), drop = FALSE] == 6L)) {
    stop("K192 donor/direction balance invariant failed", call. = FALSE)
  }
  ranked
}

.ab_scaling_generate_parent <- function(parent_spec, generator_spec,
                                        analysis_spec,
                                        registration_hashes = character()) {
  seed <- .ab_scaling_scalar(parent_spec, "seed", "integer", 0)
  parent_id <- .ab_scaling_scalar(parent_spec, "parent_id", "character")
  dimensions <- generator_spec$dimensions
  count_spec <- generator_spec$counts
  constraints_spec <- generator_spec$constraints
  donors_per_condition <- .ab_scaling_scalar(dimensions, "donors_per_condition",
                                              "integer", 1)
  cells_field <- if ("cells_per_population_per_donor" %in% names(dimensions)) {
    "cells_per_population_per_donor"
  } else "cells_per_population_per_donor_before_registered_mislabeling"
  cells_per_population <- .ab_scaling_scalar(dimensions, cells_field, "integer", 1)
  feature_count <- .ab_scaling_scalar(dimensions, "feature_count", "integer", 2)
  if (donors_per_condition != 8L || cells_per_population != 128L ||
      feature_count != 101L) {
    stop("generator dimensions differ from the frozen 16-donor/4096-cell/101-feature design",
         call. = FALSE)
  }
  offsets <- lapply(generator_spec$seed_stream_offsets, as.integer)
  specimens <- do.call(rbind, lapply(c("control", "case"), function(condition) {
    data.frame(
      sample_id = sprintf("S_%s_%02d", condition, seq_len(donors_per_condition)),
      donor_id = sprintf("D_%s_%02d", condition, seq_len(donors_per_condition)),
      condition = condition,
      stringsAsFactors = FALSE
    )
  }))
  rownames(specimens) <- NULL
  cells <- do.call(rbind, lapply(seq_len(nrow(specimens)), function(index) {
    data.frame(
      sample_id = specimens$sample_id[[index]],
      donor_id = specimens$donor_id[[index]],
      condition = specimens$condition[[index]],
      generating_label = rep(c("A", "B"), each = cells_per_population),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cells) <- NULL
  # IDs are a random permutation of opaque tokens drawn in a stream that never sees
  # labels, challenge membership, expression counts, or endpoint values.
  cells$cell_id <- .ab_scaling_with_local_seed(seed + offsets$opaque_cell_ids,
    function() sample(sprintf("C%07d", seq_len(nrow(cells))), nrow(cells), FALSE)
  )
  if (anyDuplicated(cells$cell_id) || any(grepl(
      "A|B|case|control|cost|challenge", cells$cell_id, ignore.case = TRUE))) {
    stop("opaque optimizer cell-ID invariant failed", call. = FALSE)
  }
  cells$original_label <- cells$generating_label
  original_label_recipe <- if (is.null(parent_spec$original_label_recipe)) {
    "generating_label"
  } else as.character(parent_spec$original_label_recipe)
  contamination <- data.frame(
    cell_id = character(), sample_id = character(), donor_id = character(),
    condition = character(), generating_label = character(),
    original_label = character(), stringsAsFactors = FALSE
  )
  if (identical(original_label_recipe, "case_B_observed_as_A")) {
    if (is.null(offsets$unconditional_recipe)) {
      stop("condition-dependent contamination requires its registered RNG stream",
           call. = FALSE)
    }
    contaminated_per_donor <- .ab_scaling_scalar(
      parent_spec, "contaminated_B_cells_per_case_donor", "integer", 1
    )
    contamination_key <- .ab_scaling_with_local_seed(
      seed + offsets$unconditional_recipe, function() stats::runif(nrow(cells))
    )
    chosen <- unlist(lapply(specimens$donor_id[specimens$condition == "case"],
      function(donor_id) {
        eligible <- which(cells$donor_id == donor_id & cells$generating_label == "B")
        if (length(eligible) < contaminated_per_donor) {
          stop("too few case-B cells for the registered contamination recipe",
               call. = FALSE)
        }
        utils::head(eligible[order(contamination_key[eligible], method = "radix")],
                    contaminated_per_donor)
      }), use.names = FALSE)
    cells$original_label[chosen] <- "A"
    contamination <- cells[chosen, c("cell_id", "sample_id", "donor_id", "condition",
                                      "generating_label", "original_label"), drop = FALSE]
    rownames(contamination) <- NULL
  } else if (!identical(original_label_recipe, "generating_label")) {
    stop("unsupported registered original-label recipe", call. = FALSE)
  }
  ranking <- .ab_scaling_rank_cells(
    cells, seed + offsets$ambiguity_ranking,
    max(as.integer(unlist(generator_spec$ambiguity$nested_k)))
  )
  challenge_cost <- .ab_scaling_scalar(parent_spec, "challenge_cost", "integer", 0)
  positions <- if (challenge_cost == 0L) integer() else {
    as.integer(unlist(generator_spec$ambiguity$challenge_core_positions[[
      paste0("cost_", challenge_cost)
    ]]))
  }
  challenge <- ranking[ranking$ambiguity_rank %in% positions, , drop = FALSE]
  if (nrow(challenge) != challenge_cost || any(challenge$ambiguity_rank == 1L) ||
      any(challenge$original_label != "A") || any(challenge$condition != "case") ||
      anyDuplicated(challenge$donor_id)) {
    stop("registered nonleading challenge-core invariant failed", call. = FALSE)
  }

  target_gene <- as.character(dimensions$target_gene)
  features <- c(target_gene, sprintf("%s%03d",
    as.character(dimensions$background_gene_prefix), seq_len(feature_count - 1L)
  ))
  base_background <- .ab_scaling_with_local_seed(
    seed + offsets$background_means,
    function() stats::rgamma(feature_count - 1L,
      shape = as.numeric(count_spec$background_mean_shape),
      rate = as.numeric(count_spec$background_mean_rate))
  )
  type_effect <- .ab_scaling_with_local_seed(seed + offsets$type_effects,
    function() list(
      A = .ab_scaling_gamma_one(feature_count - 1L, 0.15),
      B = .ab_scaling_gamma_one(feature_count - 1L, 0.15)
    )
  )
  donor_cv <- as.numeric(count_spec$donor_gamma_cv)
  donor_effect <- .ab_scaling_with_local_seed(seed + offsets$donor_effects,
    function() matrix(.ab_scaling_gamma_one(feature_count * nrow(specimens), donor_cv),
      nrow = feature_count, ncol = nrow(specimens),
      dimnames = list(features, specimens$donor_id))
  )
  cell_size <- .ab_scaling_with_local_seed(seed + offsets$cell_sizes,
    function() .ab_scaling_gamma_one(nrow(cells), as.numeric(count_spec$cell_size_gamma_cv))
  )
  names(cell_size) <- cells$cell_id
  effect_a <- .ab_scaling_scalar(parent_spec, "baseline_log2_effect_A", "numeric")
  multiplier <- .ab_scaling_scalar(parent_spec, "challenge_target_mean_multiplier",
                                   "numeric", 0)
  mean_matrix <- matrix(NA_real_, nrow = feature_count, ncol = nrow(cells),
                        dimnames = list(features, cells$cell_id))
  for (index in seq_len(nrow(cells))) {
    label <- cells$generating_label[[index]]
    case <- cells$condition[[index]] == "case"
    target_mean <- if (label == "A") {
      as.numeric(count_spec$target_mean_A_control) * 2^(effect_a * case)
    } else {
      as.numeric(count_spec$target_mean_B_control) *
        2^(as.numeric(count_spec$log2_effect_B) * case)
    }
    if (cells$cell_id[[index]] %in% challenge$cell_id) {
      target_mean <- target_mean * multiplier
    }
    biological_mean <- c(target_mean, base_background * type_effect[[label]])
    mean_matrix[, index] <- pmax(as.numeric(count_spec$minimum_mean),
      biological_mean * donor_effect[, cells$donor_id[[index]]] * cell_size[[index]])
  }
  dispersion <- as.numeric(count_spec$negative_binomial_dispersion)
  raw_counts <- .ab_scaling_with_local_seed(seed + offsets$negative_binomial_counts,
    function() matrix(stats::rnbinom(length(mean_matrix), mu = as.vector(mean_matrix),
      size = 1 / dispersion), nrow = nrow(mean_matrix), ncol = ncol(mean_matrix),
      dimnames = dimnames(mean_matrix))
  )
  if (any(!is.finite(raw_counts)) || any(raw_counts < 0) ||
      any(raw_counts != floor(raw_counts))) {
    stop("negative-binomial generator produced invalid counts", call. = FALSE)
  }
  counts <- methods::as(Matrix::Matrix(raw_counts, sparse = TRUE), "dgCMatrix")
  export_order <- order(cells$cell_id, method = "radix")
  cells <- cells[export_order, , drop = FALSE]
  hidden_recipe_truth <- cells[, c("cell_id", "generating_label", "original_label"),
                               drop = FALSE]
  hidden_recipe_truth$is_registered_contamination <-
    hidden_recipe_truth$cell_id %in% contamination$cell_id
  cells <- cells[, c("cell_id", "sample_id", "donor_id", "condition",
                     "original_label"), drop = FALSE]
  rownames(cells) <- NULL
  rownames(hidden_recipe_truth) <- NULL
  counts <- counts[, cells$cell_id, drop = FALSE]
  if (!identical(cells$cell_id, sort(cells$cell_id, method = "radix")) ||
      !identical(colnames(counts), cells$cell_id)) {
    stop("export order is not determined solely by opaque cell ID", call. = FALSE)
  }
  challenge$optimizer_cell_order_rank <- match(challenge$cell_id, cells$cell_id)
  k6_ids <- ranking$cell_id[ranking$ambiguity_rank <= 6L]
  challenge$optimizer_k6_candidate_rank <- match(
    challenge$cell_id, sort(k6_ids, method = "radix")
  )
  claim_spec <- generator_spec$claim
  claim <- data.frame(
    focal_gene = target_gene,
    focal_label = as.character(claim_spec$focal_label),
    alternative_label = as.character(claim_spec$alternative_label),
    case_condition = "case", control_condition = "control",
    q = as.numeric(claim_spec$q), delta = as.numeric(claim_spec$delta),
    tolerance = as.numeric(claim_spec$tolerance),
    feature_universe = paste(features, collapse = "|"),
    estimator_id = as.character(analysis_spec$analysis_id),
    stringsAsFactors = FALSE
  )
  constraints <- data.frame(
    min_cells_A = as.integer(constraints_spec$min_cells_per_sample_label),
    min_cells_B = as.integer(constraints_spec$min_cells_per_sample_label),
    min_library_A = as.integer(constraints_spec$min_library_per_sample_label),
    min_library_B = as.integer(constraints_spec$min_library_per_sample_label),
    max_changes_global = as.integer(constraints_spec$max_changes_global),
    max_changes_per_sample = as.integer(constraints_spec$max_changes_per_sample),
    max_changes_per_donor = as.integer(constraints_spec$max_changes_per_donor),
    min_donors_per_condition = as.integer(constraints_spec$min_donors_per_condition),
    stringsAsFactors = FALSE
  )
  evaluation_only <- isTRUE(generator_spec$evaluation_only)
  record <- list(
    schema_version = "1.0.0",
    generator_id = as.character(generator_spec$generator_id),
    evaluation_only = evaluation_only,
    registration_state = if (evaluation_only) {
      "freeze_ready_pending_owner_acceptance"
    } else "frozen_before_any_development_fit",
    development_only = !evaluation_only,
    parent_id = parent_id, seed = seed,
    stratum = as.character(parent_spec$stratum),
    scenario = as.character(parent_spec$scenario),
    baseline_log2_effect_A = effect_a,
    challenge_cost = challenge_cost,
    challenge_target_mean_multiplier = multiplier,
    parent_category = if (is.null(parent_spec$parent_category)) {
      as.character(parent_spec$stratum)
    } else as.character(parent_spec$parent_category),
    original_label_recipe = original_label_recipe,
    registered_contamination_cell_count = nrow(contamination),
    challenge_core = challenge,
    ambiguity_ranking_rule = as.character(generator_spec$ambiguity$construction),
    challenge_core_nonleading_in_K6 = !nrow(challenge) ||
      all(challenge$ambiguity_rank %in% 2:6),
    ambiguity_uses_expression_or_endpoint = FALSE,
    optimizer_order_rule = as.character(generator_spec$ambiguity$optimizer_order),
    optimizer_order_uses_hidden_state = FALSE,
    nested_k = as.integer(unlist(generator_spec$ambiguity$nested_k)),
    counts_and_claim_fixed_across_k = TRUE,
    independent_donors = TRUE,
    one_specimen_per_donor = TRUE,
    rng_streams = data.frame(
      role = names(offsets),
      stream_seed = seed + as.integer(unlist(offsets)),
      stringsAsFactors = FALSE
    ),
    registration_hashes = registration_hashes,
    count_diagnostics = list(
      all_integer_nonnegative = all(raw_counts >= 0 & raw_counts == floor(raw_counts)),
      zero_fraction = mean(raw_counts == 0),
      minimum_cell_library = min(colSums(raw_counts)),
      maximum_cell_library = max(colSums(raw_counts))
    )
  )
  structure(list(
    parent_id = parent_id, seed = seed, counts = counts, cells = cells,
    specimens = specimens, claim = claim, constraints = constraints,
    analysis_spec = analysis_spec, ambiguity_ranking = ranking,
    challenge_core = challenge, hidden_recipe_truth = hidden_recipe_truth,
    generation_record = record
  ), class = c("ab_scaling_parent", "list"))
}

.ab_scaling_make_domain <- function(parent, k) {
  if (!inherits(parent, "ab_scaling_parent")) {
    stop("parent must be an ab_scaling_parent", call. = FALSE)
  }
  k <- as.integer(k)
  registered_k <- as.integer(parent$generation_record$nested_k)
  if (!length(registered_k)) registered_k <- c(6L, 48L, 192L)
  if (length(k) != 1L || is.na(k) || !k %in% registered_k) {
    stop("k must be one of the parent generator's registered nested values",
         call. = FALSE)
  }
  ambiguous <- parent$ambiguity_ranking$cell_id[
    parent$ambiguity_ranking$ambiguity_rank <= k
  ]
  rows <- lapply(seq_len(nrow(parent$cells)), function(index) {
    cell_id <- parent$cells$cell_id[[index]]
    original <- parent$cells$original_label[[index]]
    labels <- if (cell_id %in% ambiguous) c("A", "B") else original
    generator_id <- as.character(parent$generation_record$generator_id)
    evidence_source <- if (isTRUE(parent$generation_record$evaluation_only)) {
      "scaling_evaluation_nested_ambiguity"
    } else "scaling_development_nested_ambiguity"
    data.frame(
      cell_id = cell_id, allowed_label = labels,
      evidence_source = evidence_source,
      evidence_version = paste0(generator_id, ":K", k),
      stringsAsFactors = FALSE
    )
  })
  allowed <- do.call(rbind, rows)
  rownames(allowed) <- NULL
  input <- list(
    counts = parent$counts, cells = parent$cells, specimens = parent$specimens,
    allowed_labels = allowed, claim = parent$claim, constraints = parent$constraints,
    analysis_spec = parent$analysis_spec
  )
  structure(list(parent_id = parent$parent_id, seed = parent$seed, k = k,
                 search_input = input), class = c("ab_scaling_domain", "list"))
}

.ab_scaling_validate_domains <- function(parent, domains, expected_k) {
  expected_k <- as.integer(expected_k)
  if (!inherits(parent, "ab_scaling_parent") || length(domains) != length(expected_k)) {
    stop("parent and all registered domains are required", call. = FALSE)
  }
  domains <- domains[order(vapply(domains, `[[`, integer(1L), "k"))]
  if (!identical(vapply(domains, `[[`, integer(1L), "k"), sort(expected_k))) {
    stop("domains do not match the registered nested K values", call. = FALSE)
  }
  invariant_names <- c("counts", "cells", "specimens", "claim", "constraints",
                       "analysis_spec")
  hashes <- vapply(invariant_names, function(name) {
    values <- lapply(domains, function(domain) domain$search_input[[name]])
    length(unique(vapply(values, digest::digest, character(1L), algo = "sha256",
                         serialize = TRUE, serializeVersion = 3L))) == 1L
  }, logical(1L))
  ambiguous <- lapply(domains, function(domain) {
    tab <- table(domain$search_input$allowed_labels$cell_id)
    sort(names(tab)[tab > 1L], method = "radix")
  })
  nested <- all(vapply(seq_len(2L), function(index) {
    all(ambiguous[[index]] %in% ambiguous[[index + 1L]])
  }, logical(1L)))
  counts <- vapply(ambiguous, length, integer(1L))
  challenge_in_k6 <- all(parent$challenge_core$cell_id %in% ambiguous[[1L]])
  ok <- all(hashes) && nested && identical(counts, sort(expected_k)) &&
    challenge_in_k6
  structure(list(ok = ok, invariant_hash_equal = hashes, nested = nested,
                 ambiguous_counts = counts, challenge_in_k6 = challenge_in_k6),
            class = c("ab_scaling_nested_check", "list"))
}

.ab_scaling_validate_nested <- function(parent, domains) {
  .ab_scaling_validate_domains(parent, domains, c(6L, 48L, 192L))
}

.ab_scaling_write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA", eol = "\n", fileEncoding = "UTF-8")
}

.ab_scaling_write_parent <- function(parent, output_dir) {
  if (!inherits(parent, "ab_scaling_parent")) {
    stop("parent must be an ab_scaling_parent", call. = FALSE)
  }
  root <- file.path(normalizePath(output_dir, mustWork = FALSE), parent$parent_id)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(parent, file.path(root, "parent.rds"), version = 3L)
  .ab_scaling_write_table(parent$ambiguity_ranking,
                          file.path(root, "ambiguity_ranking.tsv"))
  .ab_scaling_write_table(parent$challenge_core,
                          file.path(root, "challenge_core.tsv"))
  jsonlite::write_json(parent$generation_record,
                       file.path(root, "generation_record.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 17)
  registered_k <- as.integer(parent$generation_record$nested_k)
  if (!length(registered_k)) registered_k <- c(6L, 48L, 192L)
  domains <- lapply(registered_k, function(k) .ab_scaling_make_domain(parent, k))
  for (domain in domains) saveRDS(domain, file.path(root, paste0("domain_K", domain$k, ".rds")),
                                  version = 3L)
  check <- .ab_scaling_validate_domains(parent, domains, registered_k)
  jsonlite::write_json(unclass(check), file.path(root, "nested_domain_check.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  if (!isTRUE(check$ok)) stop("nested-domain check failed", call. = FALSE)
  files <- sort(list.files(root, recursive = TRUE, full.names = TRUE), method = "radix")
  list(root = root, files = files,
       hashes = setNames(vapply(files, function(path) digest::digest(
         file = path, algo = "sha256", serialize = FALSE), character(1L)), basename(files)))
}

.ab_scaling_assignment_distance <- function(first, second) {
  if (is.null(first) || is.null(second)) return(NA_integer_)
  ids <- union(names(first), names(second))
  if (!length(ids) || any(!ids %in% names(first)) || any(!ids %in% names(second))) {
    return(NA_integer_)
  }
  as.integer(sum(first[ids] != second[ids]))
}

.ab_scaling_neighbor_source <- function(data, current, direction = "all") {
  rows <- list()
  assignments <- list()
  at <- 1L
  for (cell_id in data$cells$cell_id) {
    alternatives <- setdiff(data$allowed_by_cell[[cell_id]], current[[cell_id]])
    if (identical(direction, "outgoing_only")) {
      alternatives <- alternatives[
        current[[cell_id]] == data$original_assignment[[cell_id]] &&
          data$original_assignment[[cell_id]] == data$claim$focal_label &&
          alternatives == data$claim$alternative_label
      ]
    }
    for (label in alternatives) {
      assignment <- current
      assignment[[cell_id]] <- label
      if (!.ab_direction_ok(data, assignment, direction)) next
      failures <- .ab_constraint_failures(data, assignment)
      cell <- data$cells[match(cell_id, data$cells$cell_id), , drop = FALSE]
      rows[[at]] <- data.frame(
        source_position = at, cell_id = cell_id,
        donor_id = cell$donor_id[[1L]], sample_id = cell$sample_id[[1L]],
        original_label = data$original_assignment[[cell_id]],
        parent_label = current[[cell_id]], proposed_label = label,
        direction = paste0(current[[cell_id]], "_to_", label),
        assignment_id = .ab_assignment_id(data, assignment),
        feasible = !length(failures),
        constraint_failures = paste(failures, collapse = ";"),
        stringsAsFactors = FALSE
      )
      assignments[[at]] <- assignment
      at <- at + 1L
    }
  }
  table <- if (length(rows)) do.call(rbind, rows) else data.frame(
    source_position = integer(), cell_id = character(), donor_id = character(),
    sample_id = character(), original_label = character(), parent_label = character(),
    proposed_label = character(), direction = character(), assignment_id = character(),
    feasible = logical(), constraint_failures = character(), stringsAsFactors = FALSE
  )
  table$assignment <- I(assignments)
  table
}

.ab_scaling_log_assignments <- function(data, candidate_log) {
  assignments <- list()
  original_id <- .ab_assignment_id(data, data$original_assignment)
  assignments[[original_id]] <- data$original_assignment
  if (nrow(candidate_log) && "assignment" %in% names(candidate_log)) {
    for (index in seq_len(nrow(candidate_log))) {
      value <- candidate_log$assignment[[index]]
      if (!is.null(value)) assignments[[candidate_log$assignment_id[[index]]]] <- value
    }
  }
  assignments
}

.ab_scaling_instrument_search <- function(validation, result, parent, k,
                                          outer_timing = NULL,
                                          vm_hwm_before_kb = NA_real_,
                                          vm_hwm_after_kb = NA_real_) {
  if (!inherits(validation, "ab_validation") || !isTRUE(validation$ok)) {
    stop("validation must be successful", call. = FALSE)
  }
  if (!inherits(parent, "ab_scaling_parent")) {
    stop("parent must be an ab_scaling_parent", call. = FALSE)
  }
  data <- validation$data
  log <- result$candidate_log
  if (is.null(log) || !is.data.frame(log)) log <- data.frame()
  map <- .ab_scaling_log_assignments(data, log)
  original_id <- .ab_assignment_id(data, data$original_assignment)
  original_source <- .ab_scaling_neighbor_source(
    data, data$original_assignment,
    if (is.null(result$direction)) "all" else result$direction
  )
  ambiguous <- names(data$allowed_by_cell)[
    lengths(data$allowed_by_cell) > 1L
  ]
  challenge <- parent$challenge_core
  challenge_canonical_rank <- match(
    challenge$cell_id,
    original_source$cell_id
  )

  if (nrow(log)) {
    parent_assignment <- lapply(log$parent_assignment_id, function(id) {
      if (is.na(id) || !nzchar(id)) NULL else map[[id]]
    })
    log$move_distance_from_parent <- vapply(seq_len(nrow(log)), function(index) {
      if (identical(log$proposal_type[[index]], "original")) return(0L)
      .ab_scaling_assignment_distance(parent_assignment[[index]],
                                      log$assignment[[index]])
    }, integer(1L))
    log$changed_cell_ids_from_parent <- vapply(seq_len(nrow(log)), function(index) {
      first <- parent_assignment[[index]]
      second <- log$assignment[[index]]
      if (is.null(first) || is.null(second)) return(NA_character_)
      ids <- names(first)[first != second[names(first)]]
      paste(ids, collapse = ";")
    }, character(1L))
  } else {
    log$move_distance_from_parent <- integer()
    log$changed_cell_ids_from_parent <- character()
  }

  path <- result$search_path
  state_rows <- list()
  if (is.data.frame(path) && nrow(path)) {
    for (index in seq_len(nrow(path))) {
      state_id <- path$assignment_id[[index]]
      assignment <- map[[state_id]]
      if (is.null(assignment) && identical(state_id, original_id)) {
        assignment <- data$original_assignment
      }
      source <- .ab_scaling_neighbor_source(data, assignment, result$direction)
      at_state <- log$parent_assignment_id %in% state_id
      logged_ids <- unique(log$assignment_id[at_state])
      fitted_ids <- unique(log$assignment_id[at_state & log$fit_attempted])
      successful_ids <- unique(log$assignment_id[
        at_state & log$fit_attempted & log$fit_status %in% "evaluable"
      ])
      prefix <- 0L
      while (prefix < nrow(source) &&
             source$assignment_id[[prefix + 1L]] %in% logged_ids) prefix <- prefix + 1L
      state_rows[[index]] <- data.frame(
        step = as.integer(path$step[[index]]), assignment_id = state_id,
        state_cost = as.integer(path$cost[[index]]),
        feasible_one_cell_neighbors = sum(source$feasible),
        proposed_one_cell_neighbors = sum(source$assignment_id %in% logged_ids),
        fitted_one_cell_neighbors = sum(source$assignment_id %in% fitted_ids),
        successfully_fitted_one_cell_neighbors = sum(
          source$assignment_id %in% successful_ids
        ),
        proposed_feasible_fraction = if (sum(source$feasible)) {
          sum(source$feasible & source$assignment_id %in% logged_ids) /
            sum(source$feasible)
        } else NA_real_,
        fitted_feasible_fraction = if (sum(source$feasible)) {
          sum(source$feasible & source$assignment_id %in% fitted_ids) /
            sum(source$feasible)
        } else NA_real_,
        source_prefix_generated = prefix,
        first_unvisited_source_position = if (prefix < nrow(source)) prefix + 1L else NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }
  state_coverage <- if (length(state_rows)) do.call(rbind, state_rows) else data.frame()

  cell_rows <- lapply(seq_len(nrow(original_source)), function(index) {
    source_id <- original_source$assignment_id[[index]]
    logged <- which(log$parent_assignment_id %in% original_id &
                      log$assignment_id %in% source_id)
    data.frame(
      source_position = original_source$source_position[[index]],
      cell_id = original_source$cell_id[[index]],
      donor_id = original_source$donor_id[[index]],
      sample_id = original_source$sample_id[[index]],
      direction = original_source$direction[[index]],
      ambiguity_rank = parent$ambiguity_ranking$ambiguity_rank[
        match(original_source$cell_id[[index]], parent$ambiguity_ranking$cell_id)
      ],
      challenge_core = original_source$cell_id[[index]] %in% challenge$cell_id,
      feasible = original_source$feasible[[index]],
      proposed = length(logged) > 0L,
      fit_attempted = length(logged) > 0L && any(log$fit_attempted[logged]),
      fit_evaluable = length(logged) > 0L && any(
        log$fit_status[logged] %in% "evaluable"
      ),
      stringsAsFactors = FALSE
    )
  })
  cell_coverage <- if (length(cell_rows)) do.call(rbind, cell_rows) else data.frame()
  by_donor_direction <- if (nrow(cell_coverage)) {
    keys <- interaction(cell_coverage$donor_id, cell_coverage$direction, drop = TRUE)
    do.call(rbind, lapply(split(seq_len(nrow(cell_coverage)), keys), function(index) {
      data.frame(
        donor_id = cell_coverage$donor_id[index[[1L]]],
        direction = cell_coverage$direction[index[[1L]]],
        eligible_cells = length(unique(cell_coverage$cell_id[index])),
        proposed_cells = length(unique(cell_coverage$cell_id[index][cell_coverage$proposed[index]])),
        successfully_fitted_cells = length(unique(cell_coverage$cell_id[index][
          cell_coverage$fit_evaluable[index]
        ])), stringsAsFactors = FALSE
      )
    }))
  } else data.frame()
  if (nrow(by_donor_direction)) rownames(by_donor_direction) <- NULL

  selected <- if (nrow(log)) log[log$selected, , drop = FALSE] else log
  fitted <- if (nrow(log)) log[log$fit_attempted, , drop = FALSE] else log
  practical_id <- paste(parent$parent_id, paste0("K", k), result$method,
                        paste0("fits", result$max_fits), sep = "__")
  timing <- if (is.null(outer_timing)) rep(NA_real_, 5L) else as.numeric(outer_timing)
  names(timing) <- if (length(timing) == 5L) {
    c("user.self", "sys.self", "elapsed", "user.child", "sys.child")
  } else names(timing)
  summary <- data.frame(
    run_id = practical_id, parent_id = parent$parent_id, parent_seed = parent$seed,
    k = as.integer(k), method = as.character(result$method),
    search_seed = as.integer(result$seed), max_fits = as.integer(result$max_fits),
    max_proposals = as.numeric(result$max_proposals), max_seconds = as.numeric(result$max_seconds),
    status = as.character(result$status),
    termination = as.character(result$termination_reason),
    eligible_cells = length(ambiguous), eligible_moves = nrow(original_source),
    original_feasible_one_cell_neighbors = sum(original_source$feasible),
    distinct_cells_proposed = if (nrow(log)) length(unique(unlist(strsplit(
      log$changed_cell_ids_from_parent[!is.na(log$changed_cell_ids_from_parent) &
        nzchar(log$changed_cell_ids_from_parent)], ";", fixed = TRUE)))) else 0L,
    distinct_cells_successfully_fitted = if (nrow(log)) length(unique(unlist(strsplit(
      log$changed_cell_ids_from_parent[log$fit_status %in% "evaluable" &
        !is.na(log$changed_cell_ids_from_parent) & nzchar(log$changed_cell_ids_from_parent)],
      ";", fixed = TRUE)))) else 0L,
    one_move_proposals = sum(log$move_distance_from_parent == 1L, na.rm = TRUE),
    two_move_proposals = sum(log$move_distance_from_parent == 2L, na.rm = TRUE),
    one_move_fit_attempts = sum(log$move_distance_from_parent == 1L & log$fit_attempted,
                                na.rm = TRUE),
    two_move_fit_attempts = sum(log$move_distance_from_parent == 2L & log$fit_attempted,
                                na.rm = TRUE),
    maximum_hamming_depth_fitted = if (nrow(fitted)) max(fitted$cost, na.rm = TRUE) else NA_integer_,
    maximum_hamming_depth_reached = if (is.data.frame(path) && nrow(path)) {
      max(path$cost, na.rm = TRUE)
    } else NA_integer_,
    maximum_path_steps = if (is.data.frame(path) && nrow(path)) max(path$step) else 0L,
    selected_nonoriginal_steps = sum(selected$proposal_type != "original"),
    proposal_count = as.integer(result$proposal_count),
    attempted_fit_count = as.integer(result$attempted_fit_count),
    successful_fit_count = as.integer(result$successful_fit_count),
    failed_fit_count = as.integer(result$failed_fit_count),
    internal_search_elapsed_seconds = as.numeric(result$elapsed_seconds),
    search_call_wall_seconds = if (length(timing)) unname(timing[["elapsed"]]) else NA_real_,
    search_call_user_cpu_seconds = if (length(timing)) unname(timing[["user.self"]]) else NA_real_,
    search_call_system_cpu_seconds = if (length(timing)) unname(timing[["sys.self"]]) else NA_real_,
    search_setup_and_wrapper_wall_seconds = if (length(timing)) {
      max(0, unname(timing[["elapsed"]]) - as.numeric(result$elapsed_seconds))
    } else NA_real_,
    process_cumulative_vmhwm_before_run_kb = vm_hwm_before_kb,
    process_cumulative_vmhwm_after_run_kb = vm_hwm_after_kb,
    memory_measurement_limitation = "VmHWM is cumulative within the worker process; use fresh-process MaxRSS for per-run peak",
    challenge_core_ambiguity_ranks = paste(challenge$ambiguity_rank, collapse = ";"),
    challenge_core_canonical_candidate_ranks = paste(challenge_canonical_rank,
                                                      collapse = ";"),
    interpretation = "executed implementation under registered cell order; not a generic method comparison",
    stringsAsFactors = FALSE
  )
  list(summary = summary, candidate_log = log, state_coverage = state_coverage,
       cell_coverage = cell_coverage, donor_direction_coverage = by_donor_direction,
       original_neighbor_source = original_source)
}

.ab_scaling_read_vm_hwm_kb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  value <- readLines("/proc/self/status", warn = FALSE)
  value <- value[startsWith(value, "VmHWM:")]
  if (!length(value)) return(NA_real_)
  suppressWarnings(as.numeric(strsplit(trimws(value[[1L]]), "[[:space:]]+")[[1L]][[2L]]))
}
