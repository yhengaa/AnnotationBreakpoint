.ab_public_hash <- function(namespace, seed, ...) {
  fields <- c(namespace, as.character(as.integer(seed)), as.character(c(...)))
  digest::digest(paste(fields, collapse = "|"), algo = "sha256",
                 serialize = FALSE)
}

.ab_public_hash_rank <- function(values, namespace, seed, ...) {
  extra <- list(...)
  hashes <- vapply(values, function(value) {
    do.call(.ab_public_hash, c(list(namespace, seed), extra, list(value)))
  }, character(1L))
  order(hashes, values, method = "radix")
}

.ab_public_condition_table <- function(spec) {
  donors <- as.character(unlist(spec$selection$evaluation_donors))
  ranked <- donors[.ab_public_hash_rank(
    donors, spec$artificial_conditions$namespace,
    spec$artificial_conditions$seed
  )]
  expected_case <- sort(ranked[seq_len(3L)], method = "radix")
  declared_case <- sort(as.character(unlist(
    spec$artificial_conditions$case_donors
  )), method = "radix")
  if (!identical(expected_case, declared_case)) {
    stop("declared artificial case donors do not match the registered hash rule",
         call. = FALSE)
  }
  data.frame(
    donor_id = donors,
    sample_id = paste0("Hao5P_day0_", donors),
    condition = ifelse(donors %in% declared_case, "case", "control"),
    artificial_group = TRUE,
    stringsAsFactors = FALSE
  )
}

.ab_public_select_gene_panel <- function(snapshot, spec) {
  calibration_donors <- as.character(unlist(
    spec$selection$calibration_donors
  ))
  metadata <- snapshot$selected_cells
  index <- which(metadata$donor %in% calibration_donors)
  counts <- snapshot$rna_counts[, index, drop = FALSE]
  blocks <- paste(metadata$donor[index], metadata$celltype.l2[index], sep = "::")
  block_levels <- as.vector(outer(
    calibration_donors, as.character(unlist(spec$selection$populations)),
    paste, sep = "::"
  ))
  incidence <- Matrix::sparseMatrix(
    i = seq_along(index), j = match(blocks, block_levels), x = 1,
    dims = c(length(index), length(block_levels)),
    dimnames = list(colnames(counts), block_levels)
  )
  block_counts <- counts %*% incidence
  totals <- as.numeric(Matrix::rowSums(counts))
  detected <- as.numeric(Matrix::rowSums(counts > 0))
  minimum_block <- apply(as.matrix(block_counts), 1L, min)
  feature_ids <- rownames(counts)
  prefixes <- as.character(unlist(spec$gene_panel$exclusions$prefixes))
  exact_exclusions <- as.character(unlist(spec$gene_panel$exclusions$exact))
  prefix_excluded <- Reduce(
    `|`, lapply(prefixes, function(prefix) startsWith(feature_ids, prefix)),
    init = rep(FALSE, length(feature_ids))
  )
  eligible <- !prefix_excluded & !feature_ids %in% exact_exclusions &
    totals >= as.numeric(
      spec$gene_panel$eligibility$minimum_total_calibration_count
    ) & minimum_block >= as.numeric(
      spec$gene_panel$eligibility$minimum_count_each_donor_population_block
    ) & detected >= as.numeric(
      spec$gene_panel$eligibility$minimum_detected_calibration_cells
    )
  eligible_table <- data.frame(
    feature_id = feature_ids[eligible],
    calibration_total = totals[eligible],
    calibration_detected_cells = detected[eligible],
    minimum_donor_population_count = minimum_block[eligible],
    stringsAsFactors = FALSE
  )
  eligible_table <- eligible_table[order(
    eligible_table$calibration_total, eligible_table$feature_id,
    method = "radix"
  ), , drop = FALSE]
  rownames(eligible_table) <- NULL
  if (nrow(eligible_table) < 5L) {
    stop("calibration-only rule selected fewer than five eligible genes",
         call. = FALSE)
  }
  quantiles <- c(
    primary_no_injection = as.numeric(
      spec$gene_panel$selection$primary_no_injection_quantile
    ),
    primary_subtle = as.numeric(
      spec$gene_panel$selection$primary_subtle_quantile
    ),
    primary_strong = as.numeric(
      spec$gene_panel$selection$primary_strong_quantile
    ),
    control_1 = as.numeric(unlist(
      spec$gene_panel$selection$non_injected_control_quantiles
    ))[[1L]],
    control_2 = as.numeric(unlist(
      spec$gene_panel$selection$non_injected_control_quantiles
    ))[[2L]]
  )
  nominal <- round(1 + quantiles * (nrow(eligible_table) - 1L))
  selected <- integer()
  for (position in nominal) {
    candidates <- setdiff(seq_len(nrow(eligible_table)), selected)
    distance <- abs(candidates - position)
    selected <- c(selected, candidates[order(distance, candidates)[[1L]]])
  }
  panel <- eligible_table[selected, , drop = FALSE]
  panel$panel_role <- names(quantiles)
  panel$selection_quantile <- unname(quantiles)
  panel$eligible_rank <- selected
  panel$injected <- panel$panel_role %in% c(
    "primary_subtle", "primary_strong"
  )
  panel$nominal_log2_fold_change <- c(
    as.numeric(spec$injection$cases$no_injection$nominal_log2_fold_change),
    as.numeric(spec$injection$cases$subtle_injection$nominal_log2_fold_change),
    as.numeric(spec$injection$cases$strong_injection$nominal_log2_fold_change),
    0, 0
  )
  panel$selection_information <- "P1/P2 capped calibration RNA counts only"
  panel <- panel[, c(
    "panel_role", "feature_id", "injected", "nominal_log2_fold_change",
    "selection_quantile", "eligible_rank", "calibration_total",
    "calibration_detected_cells", "minimum_donor_population_count",
    "selection_information"
  )]
  list(panel = panel, eligible_gene_count = nrow(eligible_table),
       eligible = eligible_table)
}

.ab_public_ambiguity <- function(selected_cells, spec) {
  donors <- as.character(unlist(spec$selection$evaluation_donors))
  population_rank <- donors[.ab_public_hash_rank(
    donors, spec$ambiguity$namespace, spec$ambiguity$seed
  )]
  chosen_population <- setNames(
    c(rep("CD14 Mono", 3L), rep("CD16 Mono", 3L)), population_rank
  )
  rows <- lapply(donors, function(donor) {
    population <- chosen_population[[donor]]
    candidates <- selected_cells[
      selected_cells$donor == donor &
        selected_cells$celltype.l2 == population, , drop = FALSE
    ]
    rank <- .ab_public_hash_rank(
      candidates$cell_id, "hao5p-ambiguity-cell-v1",
      spec$ambiguity$seed, donor, population
    )
    candidates[rank[[1L]], , drop = FALSE]
  })
  chosen <- do.call(rbind, rows)
  corruption_rank <- donors[.ab_public_hash_rank(
    donors, spec$ambiguity$corruption$namespace,
    spec$ambiguity$corruption$seed
  )]
  corrupted <- sort(corruption_rank[seq_len(
    as.integer(spec$ambiguity$corruption$count)
  )], method = "radix")
  declared <- sort(as.character(unlist(
    spec$ambiguity$corruption$donors
  )), method = "radix")
  if (!identical(corrupted, declared)) {
    stop("declared corruption donors do not match the registered hash rule",
         call. = FALSE)
  }
  reference <- chosen$celltype.l2
  observed <- reference
  is_corrupted <- chosen$donor %in% corrupted
  observed[is_corrupted] <- ifelse(
    reference[is_corrupted] == "CD14 Mono", "CD16 Mono", "CD14 Mono"
  )
  data.frame(
    cell_id = chosen$cell_id,
    donor_id = chosen$donor,
    sample_id = paste0("Hao5P_day0_", chosen$donor),
    source_lane = chosen$lane,
    reference_label = reference,
    observed_label = observed,
    corrupted = is_corrupted,
    alternative_label = ifelse(
      observed == "CD14 Mono", "CD16 Mono", "CD14 Mono"
    ),
    stringsAsFactors = FALSE
  )
}

.ab_public_inject_counts <- function(source_counts, reference_cells, panel,
                                     conditions, spec) {
  counts <- source_counts
  injection_rows <- list()
  at <- 1L
  affected_donors <- conditions$donor_id[conditions$condition == "case"]
  targets <- panel[panel$injected, , drop = FALSE]
  for (target_index in seq_len(nrow(targets))) {
    gene <- targets$feature_id[[target_index]]
    theta <- targets$nominal_log2_fold_change[[target_index]]
    for (donor in affected_donors) {
      ids <- reference_cells$cell_id[
        reference_cells$donor == donor &
          reference_cells$celltype.l2 == spec$injection$affected_reference_population
      ]
      old <- as.numeric(source_counts[gene, ids, drop = TRUE])
      old_total <- sum(old)
      target_total <- round(old_total * 2^theta)
      addition_total <- as.integer(target_total - old_total)
      weights <- old + 1
      quota <- addition_total * weights / sum(weights)
      addition <- floor(quota)
      remainder <- addition_total - sum(addition)
      tie_hash <- vapply(ids, function(id) .ab_public_hash(
        "hao5p-injection-v1", spec$injection$seed, gene, id
      ), character(1L))
      if (remainder > 0L) {
        order_index <- order(-(quota - addition), tie_hash, ids,
                             method = "radix")
        addition[order_index[seq_len(remainder)]] <-
          addition[order_index[seq_len(remainder)]] + 1
      }
      new <- old + addition
      counts[gene, ids] <- new
      injection_rows[[at]] <- data.frame(
        feature_id = gene, panel_role = targets$panel_role[[target_index]],
        donor_id = donor, cell_id = ids, source_count = old,
        added_count = addition, frozen_count = new,
        nominal_log2_fold_change = theta,
        source_pseudobulk_total = old_total,
        frozen_pseudobulk_total = sum(new),
        target_pseudobulk_total = target_total,
        stringsAsFactors = FALSE
      )
      at <- at + 1L
    }
  }
  counts <- methods::as(Matrix::drop0(counts), "dgCMatrix")
  injection <- if (length(injection_rows)) {
    do.call(rbind, injection_rows)
  } else data.frame()
  list(counts = counts, injection = injection)
}

.ab_prepare_public_benchmark <- function(snapshot, spec, analysis_spec) {
  if (!inherits(snapshot, "ab_hao5p_snapshot")) {
    stop("snapshot must be an ab_hao5p_snapshot", call. = FALSE)
  }
  conditions <- .ab_public_condition_table(spec)
  gene_selection <- .ab_public_select_gene_panel(snapshot, spec)
  evaluation_donors <- conditions$donor_id
  evaluation_index <- which(snapshot$selected_cells$donor %in% evaluation_donors)
  evaluation_index <- evaluation_index[order(
    snapshot$selected_cells$cell_id[evaluation_index], method = "radix"
  )]
  reference_cells <- snapshot$selected_cells[evaluation_index, , drop = FALSE]
  source_counts <- snapshot$rna_counts[, evaluation_index, drop = FALSE]
  if (!identical(colnames(source_counts), reference_cells$cell_id)) {
    stop("evaluation source counts and cell metadata are not aligned",
         call. = FALSE)
  }
  ambiguity <- .ab_public_ambiguity(reference_cells, spec)
  observed <- setNames(reference_cells$celltype.l2, reference_cells$cell_id)
  observed[ambiguity$cell_id] <- ambiguity$observed_label
  cells <- data.frame(
    cell_id = reference_cells$cell_id,
    sample_id = paste0("Hao5P_day0_", reference_cells$donor),
    donor_id = reference_cells$donor,
    condition = conditions$condition[
      match(reference_cells$donor, conditions$donor_id)
    ],
    original_label = unname(observed[reference_cells$cell_id]),
    stringsAsFactors = FALSE
  )
  specimens <- conditions[, c("sample_id", "donor_id", "condition")]
  ambiguous_ids <- ambiguity$cell_id
  allowed_rows <- lapply(seq_len(nrow(cells)), function(index) {
    cell_id <- cells$cell_id[[index]]
    label <- cells$original_label[[index]]
    if (cell_id %in% ambiguous_ids) {
      data.frame(
        cell_id = cell_id,
        allowed_label = c(label, ifelse(
          label == "CD14 Mono", "CD16 Mono", "CD14 Mono"
        )),
        evidence_source = c(
          "semi_synthetic_observed_label", "curated_related_population_pair"
        ),
        evidence_version = "Hao5P-WNN-semisynthetic-v1",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        cell_id = cell_id, allowed_label = label,
        evidence_source = "published_WNN_label",
        evidence_version = "Hao5P-GSE164378-celltype.l2",
        stringsAsFactors = FALSE
      )
    }
  })
  allowed_labels <- do.call(rbind, allowed_rows)
  probabilities <- do.call(rbind, lapply(ambiguous_ids, function(cell_id) {
    labels <- allowed_labels$allowed_label[allowed_labels$cell_id == cell_id]
    data.frame(cell_id = cell_id, label = labels, probability = 0.5,
               probability_source = "stipulated_uniform_v1",
               stringsAsFactors = FALSE)
  }))
  injected <- .ab_public_inject_counts(
    source_counts, reference_cells, gene_selection$panel, conditions, spec
  )
  search_input <- list(
    counts = injected$counts,
    cells = cells,
    specimens = specimens,
    allowed_labels = allowed_labels,
    analysis_spec = analysis_spec
  )
  evaluation_truth <- list(
    source_counts = source_counts,
    reference_cells = reference_cells,
    ambiguity_truth = ambiguity,
    injection_truth = injected$injection,
    published_reference_assignment = setNames(
      reference_cells$celltype.l2, reference_cells$cell_id
    )
  )
  structure(list(
    schema_version = "1.0.0",
    benchmark_id = spec$benchmark_id,
    source_snapshot_id = snapshot$snapshot_id,
    conditions = conditions,
    gene_panel = gene_selection$panel,
    eligible_gene_count = gene_selection$eligible_gene_count,
    search_input = search_input,
    method_input = list(probabilities = probabilities),
    evaluation_truth = evaluation_truth,
    spec = spec
  ), class = c("ab_public_benchmark", "list"))
}

.ab_public_claim <- function(benchmark, feature_id) {
  panel <- benchmark$gene_panel
  if (length(feature_id) != 1L || is.na(feature_id) ||
      !feature_id %in% panel$feature_id) {
    stop("feature_id must identify one registered public-panel gene",
         call. = FALSE)
  }
  spec <- benchmark$spec
  list(
    focal_gene = as.character(feature_id),
    focal_label = as.character(spec$claim$focal_label),
    alternative_label = as.character(spec$claim$alternative_label),
    case_condition = as.character(spec$claim$case_condition),
    control_condition = as.character(spec$claim$control_condition),
    q = as.integer(spec$claim$q),
    delta = as.numeric(spec$claim$delta),
    tolerance = as.numeric(spec$claim$tolerance),
    estimator_id = as.character(spec$claim$estimator_id),
    feature_universe = rownames(benchmark$search_input$counts)
  )
}

.ab_public_constraints <- function(benchmark, variant_id) {
  variants <- benchmark$spec$constraint_variants
  ids <- vapply(variants, function(value) as.character(value$id), character(1L))
  index <- match(variant_id, ids)
  if (is.na(index)) stop("unknown public constraint variant", call. = FALSE)
  variant <- variants[[index]]
  labels <- c(
    as.character(benchmark$spec$claim$focal_label),
    as.character(benchmark$spec$claim$alternative_label)
  )
  out <- list(
    max_changes_global = as.integer(variant$max_changes_global),
    max_changes_per_sample = as.integer(variant$max_changes_per_sample),
    max_changes_per_donor = as.integer(variant$max_changes_per_donor),
    min_donors_per_condition = 2L
  )
  out[[paste0("min_cells_", labels[[1L]])]] <-
    as.integer(variant$min_cells_per_sample_label)
  out[[paste0("min_cells_", labels[[2L]])]] <-
    as.integer(variant$min_cells_per_sample_label)
  out[[paste0("min_library_", labels[[1L]])]] <-
    as.integer(variant$min_library_count_per_sample_label)
  out[[paste0("min_library_", labels[[2L]])]] <-
    as.integer(variant$min_library_count_per_sample_label)
  out
}

.ab_public_validation <- function(benchmark, feature_id, variant_id) {
  input <- benchmark$search_input
  validate_inputs(
    input$counts, input$cells, input$specimens, input$allowed_labels,
    .ab_public_claim(benchmark, feature_id),
    .ab_public_constraints(benchmark, variant_id), input$analysis_spec
  )
}

.ab_public_map_assignment <- function(benchmark) {
  input <- benchmark$search_input
  probabilities <- benchmark$method_input$probabilities
  assignment <- setNames(input$cells$original_label, input$cells$cell_id)
  rows <- split(
    probabilities,
    factor(probabilities$cell_id, levels = unique(probabilities$cell_id))
  )
  for (cell_id in names(rows)) {
    block <- rows[[cell_id]]
    maximum <- max(block$probability)
    candidates <- block$label[block$probability == maximum]
    observed <- assignment[[cell_id]]
    assignment[[cell_id]] <- if (observed %in% candidates) {
      observed
    } else sort(candidates, method = "radix")[[1L]]
  }
  assignment
}
