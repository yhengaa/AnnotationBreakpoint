# Outcome-blind helpers for the bounded COMBAT public-data eligibility screen.
# These functions use donor, sample, cell-label, and assay-support metadata only.

.pu_require_columns <- function(x, required, table_name) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      table_name, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pu_read_combat_metadata <- function(clinical_path, composition_path) {
  clinical <- utils::read.delim(
    clinical_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
  )
  composition <- utils::read.csv(
    composition_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  .pu_require_columns(
    clinical,
    c(
      "COMBAT_participant_timepoint_ID", "COMBAT_ID", "scRNASeq_sample_ID",
      "Source", "Age", "Sex", "Institute", "TimeSinceOnset"
    ),
    "clinical"
  )
  .pu_require_columns(
    composition,
    c("scRNASeq_sample_ID", "MinorSubset", "MinorSubset_Count"),
    "composition"
  )
  list(clinical = clinical, composition = composition)
}

pu_combat_one_sample_candidates <- function(
    clinical,
    composition,
    labels = c("CD4.NAIVE", "CD4.TCM"),
    minimum_cells = 20L,
    conditions = c("COVID_SEV", "COVID_CRIT")) {
  stopifnot(length(labels) == 2L, length(conditions) == 2L)
  clinical <- clinical[
    !is.na(clinical$scRNASeq_sample_ID) & clinical$Source %in% conditions,
    ,
    drop = FALSE
  ]
  if (!nrow(clinical)) stop("no qualifying CITE-seq clinical rows", call. = FALSE)

  severity_rank <- stats::setNames(seq_along(conditions), conditions)
  by_donor <- split(clinical, clinical$COMBAT_ID)
  selected <- lapply(by_donor, function(rows) {
    rank <- unname(severity_rank[rows$Source])
    rows <- rows[rank == max(rank), , drop = FALSE]
    rows <- rows[order(
      rows$COMBAT_participant_timepoint_ID,
      rows$scRNASeq_sample_ID,
      method = "radix"
    ), , drop = FALSE]
    rows[1L, , drop = FALSE]
  })
  selected <- do.call(rbind, selected)
  rownames(selected) <- NULL

  keep <- composition$MinorSubset %in% labels
  counts <- stats::reshape(
    composition[keep, c(
      "scRNASeq_sample_ID", "MinorSubset", "MinorSubset_Count"
    )],
    idvar = "scRNASeq_sample_ID",
    timevar = "MinorSubset",
    direction = "wide"
  )
  names(counts) <- sub("^MinorSubset_Count\\.", "", names(counts))
  selected <- merge(
    selected,
    counts,
    by = "scRNASeq_sample_ID",
    all.x = TRUE,
    sort = FALSE
  )
  selected <- selected[order(selected$COMBAT_ID, method = "radix"), , drop = FALSE]
  selected$eligible_population_counts <-
    !is.na(selected[[labels[[1L]]]]) &
    !is.na(selected[[labels[[2L]]]]) &
    selected[[labels[[1L]]]] >= minimum_cells &
    selected[[labels[[2L]]]] >= minimum_cells
  selected
}

.pu_choose_time_matched_pairs <- function(case_rows, control_rows) {
  order_rows <- function(x) {
    time <- suppressWarnings(as.numeric(x$TimeSinceOnset))
    time[is.na(time)] <- Inf
    x[order(time, x$COMBAT_ID, x$scRNASeq_sample_ID, method = "radix"), , drop = FALSE]
  }
  case_rows <- order_rows(case_rows)
  control_rows <- order_rows(control_rows)
  n <- min(nrow(case_rows), nrow(control_rows))
  if (!n) return(NULL)

  permutations <- function(values) {
    if (length(values) <= 1L) return(list(values))
    unlist(lapply(seq_along(values), function(i) {
      lapply(permutations(values[-i]), function(rest) c(values[[i]], rest))
    }), recursive = FALSE)
  }
  case_subsets <- if (nrow(case_rows) == n) {
    list(seq_len(n))
  } else {
    utils::combn(seq_len(nrow(case_rows)), n, simplify = FALSE)
  }
  control_subsets <- if (nrow(control_rows) == n) {
    list(seq_len(n))
  } else {
    utils::combn(seq_len(nrow(control_rows)), n, simplify = FALSE)
  }
  alternatives <- list()
  at <- 1L
  for (case_idx in case_subsets) {
    for (control_idx in control_subsets) {
      for (control_order in permutations(control_idx)) {
        crows <- case_rows[case_idx, , drop = FALSE]
        rrows <- control_rows[control_order, , drop = FALSE]
        batch_mismatch <- if (all(c("Pool_ID") %in% names(crows)) &&
                              all(c("Pool_ID") %in% names(rrows))) {
          sum(crows$Pool_ID != rrows$Pool_ID)
        } else 0L
        time_difference <- sum(abs(
          suppressWarnings(as.numeric(crows$TimeSinceOnset)) -
            suppressWarnings(as.numeric(rrows$TimeSinceOnset))
        ))
        if (!is.finite(time_difference)) time_difference <- Inf
        signature <- paste(
          paste(crows$COMBAT_ID, rrows$COMBAT_ID, sep = ":"),
          collapse = "|"
        )
        alternatives[[at]] <- list(
          case = crows,
          control = rrows,
          batch_mismatch = batch_mismatch,
          time_difference = time_difference,
          signature = signature
        )
        at <- at + 1L
      }
    }
  }
  score <- data.frame(
    batch = vapply(alternatives, `[[`, integer(1L), "batch_mismatch"),
    time = vapply(alternatives, `[[`, numeric(1L), "time_difference"),
    signature = vapply(alternatives, `[[`, character(1L), "signature"),
    stringsAsFactors = FALSE
  )
  best <- alternatives[[order(
    score$batch, score$time, score$signature, method = "radix"
  )[[1L]]]]
  case_rows <- best$case
  control_rows <- best$control
  data.frame(
    match_pair_id = sprintf("P%02d", seq_len(n)),
    case_donor_id = case_rows$COMBAT_ID,
    case_sample_id = case_rows$scRNASeq_sample_ID,
    control_donor_id = control_rows$COMBAT_ID,
    control_sample_id = control_rows$scRNASeq_sample_ID,
    age_band = as.character(case_rows$Age),
    sex_code = as.character(case_rows$Sex),
    case_time_since_onset = suppressWarnings(as.numeric(case_rows$TimeSinceOnset)),
    control_time_since_onset = suppressWarnings(as.numeric(control_rows$TimeSinceOnset)),
    case_pool_id = if ("Pool_ID" %in% names(case_rows)) case_rows$Pool_ID else NA_character_,
    control_pool_id = if ("Pool_ID" %in% names(control_rows)) control_rows$Pool_ID else NA_character_,
    pool_matched = if ("Pool_ID" %in% names(case_rows)) {
      case_rows$Pool_ID == control_rows$Pool_ID
    } else NA,
    stringsAsFactors = FALSE
  )
}

pu_match_combat_panel <- function(
    candidates,
    case_condition = "COVID_CRIT",
    control_condition = "COVID_SEV",
    target_per_condition = 12L) {
  candidates <- candidates[candidates$eligible_population_counts, , drop = FALSE]
  candidates$stratum <- paste(candidates$Age, candidates$Sex, sep = "|")
  strata <- sort(intersect(
    unique(candidates$stratum[candidates$Source == case_condition]),
    unique(candidates$stratum[candidates$Source == control_condition])
  ), method = "radix")
  pairs <- lapply(strata, function(stratum) {
    case_rows <- candidates[
      candidates$Source == case_condition & candidates$stratum == stratum,
      , drop = FALSE
    ]
    control_rows <- candidates[
      candidates$Source == control_condition & candidates$stratum == stratum,
      , drop = FALSE
    ]
    out <- .pu_choose_time_matched_pairs(case_rows, control_rows)
    if (!is.null(out)) out$stratum <- stratum
    out
  })
  pairs <- do.call(rbind, pairs)
  if (is.null(pairs) || nrow(pairs) < target_per_condition) {
    stop("fewer than target exact age-band/sex matched donor pairs", call. = FALSE)
  }
  pairs <- pairs[order(
    abs(pairs$case_time_since_onset - pairs$control_time_since_onset),
    pairs$stratum,
    pairs$case_donor_id,
    pairs$control_donor_id,
    method = "radix"
  ), , drop = FALSE]
  pairs <- pairs[seq_len(target_per_condition), , drop = FALSE]
  pairs <- pairs[order(pairs$stratum, pairs$case_donor_id, method = "radix"), , drop = FALSE]
  pairs$match_pair_id <- sprintf("P%02d", seq_len(nrow(pairs)))
  rownames(pairs) <- NULL
  pairs
}

pu_assign_development_partition <- function(pairs, development_pairs = 4L,
                                            seed = "combat-public-v1") {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest is required for deterministic donor partitioning", call. = FALSE)
  }
  if (development_pairs <= 0L || development_pairs >= nrow(pairs)) {
    stop("development_pairs must leave nonempty development and evaluation sets", call. = FALSE)
  }
  pair_hash <- vapply(seq_len(nrow(pairs)), function(i) {
    digest::digest(
      paste(
        seed,
        pairs$case_donor_id[[i]], pairs$case_sample_id[[i]],
        pairs$control_donor_id[[i]], pairs$control_sample_id[[i]],
        sep = "|"
      ),
      algo = "sha256",
      serialize = FALSE
    )
  }, character(1L))
  pairs$partition <- "evaluation"
  alternatives <- utils::combn(
    seq_len(nrow(pairs)), development_pairs, simplify = FALSE
  )
  distribution_score <- function(index, field) {
    levels <- sort(unique(pairs[[field]]), method = "radix")
    total <- table(factor(pairs[[field]], levels = levels))
    chosen <- table(factor(pairs[[field]][index], levels = levels))
    sum(abs(as.numeric(chosen) - development_pairs * as.numeric(total) / nrow(pairs)))
  }
  batch_fields <- intersect(
    c("case_pool_id", "control_pool_id"), names(pairs)
  )
  batch_score <- vapply(alternatives, function(index) {
    if (!length(batch_fields) || anyNA(pairs[, batch_fields, drop = FALSE])) return(0)
    sum(vapply(batch_fields, function(field) {
      distribution_score(index, field)
    }, numeric(1L)))
  }, numeric(1L))
  stratum_score <- vapply(
    alternatives,
    distribution_score,
    numeric(1L),
    field = "stratum"
  )
  signatures <- vapply(alternatives, function(index) {
    digest::digest(
      paste(seed, paste(sort(pair_hash[index]), collapse = "|"), sep = "|"),
      algo = "sha256",
      serialize = FALSE
    )
  }, character(1L))
  best <- alternatives[[order(
    batch_score, stratum_score, signatures, method = "radix"
  )[[1L]]]]
  pairs$partition[best] <- "development"
  stopifnot(sum(pairs$partition == "development") == development_pairs)
  pairs$partition_seed <- seed
  pairs$partition_hash <- pair_hash
  pairs
}

pu_percentage_budgets <- function(ambiguous_by_donor,
                                  percentages = c(0.01, 0.05)) {
  donor_names <- names(ambiguous_by_donor)
  ambiguous_by_donor <- as.integer(ambiguous_by_donor)
  names(ambiguous_by_donor) <- donor_names
  if (is.null(names(ambiguous_by_donor)) || any(ambiguous_by_donor < 0L)) {
    stop("ambiguous_by_donor must be a named nonnegative integer vector", call. = FALSE)
  }
  total <- sum(ambiguous_by_donor)
  do.call(rbind, lapply(percentages, function(p) {
    if (!is.finite(p) || p <= 0 || p > 1) stop("invalid budget percentage", call. = FALSE)
    data.frame(
      budget_name = paste0(format(100 * p, trim = TRUE), "pct"),
      percentage = p,
      ambiguous_cells = total,
      global_max_cost = if (total) max(1L, ceiling(p * total)) else 0L,
      donor_id = names(ambiguous_by_donor),
      donor_ambiguous_cells = ambiguous_by_donor,
      donor_max_cost = ifelse(
        ambiguous_by_donor > 0L,
        pmax(1L, ceiling(p * ambiguous_by_donor)),
        0L
      ),
      stringsAsFactors = FALSE
    )
  }))
}

.pu_require_rhdf5 <- function() {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop(
      "COMBAT H5AD adapters require the optional rhdf5 package",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.pu_h5ad_read_column <- function(path, frame, column) {
  .pu_require_rhdf5()
  base <- paste0(frame, "/", column)
  listing <- rhdf5::h5ls(path, recursive = TRUE, datasetinfo = FALSE)
  direct <- listing$group == frame & listing$name == column
  if (!any(direct)) stop("H5AD column is absent: ", base, call. = FALSE)
  children <- rhdf5::h5ls(path, recursive = TRUE, datasetinfo = FALSE)
  available <- paste0(children$group, "/", children$name)
  legacy_category_path <- paste0(frame, "/__categories/", column)
  if (legacy_category_path %in% available) {
    categories <- rhdf5::h5read(path, legacy_category_path)
    codes <- as.integer(rhdf5::h5read(path, base))
    out <- rep(NA_character_, length(codes))
    keep <- codes >= 0L
    out[keep] <- as.character(categories[codes[keep] + 1L])
    return(out)
  }
  if (listing$otype[which(direct)[[1L]]] == "H5I_DATASET") {
    return(rhdf5::h5read(path, base))
  }
  category_path <- paste0(base, "/categories")
  code_path <- paste0(base, "/codes")
  if (!all(c(category_path, code_path) %in% available)) {
    stop("unsupported H5AD column encoding: ", base, call. = FALSE)
  }
  categories <- rhdf5::h5read(path, category_path)
  codes <- as.integer(rhdf5::h5read(path, code_path))
  out <- rep(NA_character_, length(codes))
  keep <- codes >= 0L
  out[keep] <- as.character(categories[codes[keep] + 1L])
  out
}

pu_h5ad_frame_columns <- function(path, frame = "/obs") {
  .pu_require_rhdf5()
  listing <- rhdf5::h5ls(path, recursive = TRUE, datasetinfo = FALSE)
  sort(listing$name[listing$group == frame], method = "radix")
}

.pu_stream_csr_rows <- function(path, matrix_group, selected_rows,
                                n_columns, block_rows = 10000L,
                                keep_columns = NULL) {
  .pu_require_rhdf5()
  indptr <- as.numeric(rhdf5::h5read(path, paste0(matrix_group, "/indptr")))
  n_rows <- length(indptr) - 1L
  selected_rows <- as.integer(selected_rows)
  if (!length(selected_rows) || any(selected_rows < 1L) ||
      any(selected_rows > n_rows)) {
    stop("selected_rows are outside the H5AD matrix", call. = FALSE)
  }
  if (anyDuplicated(selected_rows)) {
    stop("selected_rows must be unique", call. = FALSE)
  }
  selected_lookup <- integer(n_rows)
  selected_lookup[selected_rows] <- seq_along(selected_rows)
  p <- integer(length(selected_rows) + 1L)
  j_parts <- vector("list", length(selected_rows))
  x_parts <- vector("list", length(selected_rows))

  starts <- seq.int(1L, n_rows, by = as.integer(block_rows))
  for (block_start in starts) {
    block_end <- min(n_rows, block_start + as.integer(block_rows) - 1L)
    chosen <- selected_rows[selected_rows >= block_start & selected_rows <= block_end]
    if (!length(chosen)) next
    data_start <- indptr[[block_start]] + 1
    data_end <- indptr[[block_end + 1L]]
    if (data_end < data_start) next
    index <- as.integer(rhdf5::h5read(
      path,
      paste0(matrix_group, "/indices"),
      start = data_start,
      count = data_end - data_start + 1
    ))
    value <- as.numeric(rhdf5::h5read(
      path,
      paste0(matrix_group, "/data"),
      start = data_start,
      count = data_end - data_start + 1
    ))
    offset <- indptr[[block_start]]
    for (row in chosen) {
      at <- selected_lookup[[row]]
      lo <- indptr[[row]] - offset + 1
      hi <- indptr[[row + 1L]] - offset
      if (hi < lo) next
      row_index <- index[seq.int(lo, hi)] + 1L
      row_value <- value[seq.int(lo, hi)]
      if (!is.null(keep_columns)) {
        keep <- row_index %in% keep_columns
        row_index <- row_index[keep]
        row_value <- row_value[keep]
      }
      j_parts[[at]] <- row_index - 1L
      x_parts[[at]] <- row_value
    }
    rm(index, value)
  }
  lengths <- lengths(j_parts)
  p[-1L] <- cumsum(lengths)
  j <- as.integer(unlist(j_parts, use.names = FALSE))
  x <- as.numeric(unlist(x_parts, use.names = FALSE))
  if (!is.null(keep_columns)) {
    mapped <- match(j + 1L, keep_columns)
    if (anyNA(mapped)) stop("internal selected-column mapping failed", call. = FALSE)
    j <- as.integer(mapped - 1L)
    n_columns <- length(keep_columns)
  }
  methods::new(
    "dgRMatrix",
    p = as.integer(p),
    j = j,
    x = x,
    Dim = as.integer(c(length(selected_rows), n_columns))
  )
}

pu_h5ad_extract_raw_counts <- function(path, selected_rows, feature_ids,
                                       rna_columns, block_rows = 10000L) {
  if (anyDuplicated(feature_ids) || anyDuplicated(rna_columns)) {
    stop("feature IDs and selected RNA columns must be unique", call. = FALSE)
  }
  csr <- .pu_stream_csr_rows(
    path,
    "/layers/raw",
    selected_rows,
    n_columns = length(feature_ids),
    block_rows = block_rows,
    keep_columns = rna_columns
  )
  counts <- Matrix::t(methods::as(csr, "CsparseMatrix"))
  if (!inherits(counts, "dgCMatrix")) {
    counts <- methods::as(counts, "dgCMatrix")
  }
  rownames(counts) <- feature_ids[rna_columns]
  counts
}

pu_h5ad_extract_columns <- function(path, matrix_group, selected_rows,
                                    all_feature_count, selected_columns,
                                    block_rows = 10000L) {
  csr <- .pu_stream_csr_rows(
    path,
    matrix_group,
    selected_rows,
    n_columns = all_feature_count,
    block_rows = block_rows,
    keep_columns = selected_columns
  )
  as.matrix(csr)
}

pu_equal_waterfill_quota <- function(block_sizes, total_limit,
                                     minimum_each = 20L) {
  block_names <- names(block_sizes)
  block_sizes <- as.integer(block_sizes)
  names(block_sizes) <- block_names
  if (is.null(names(block_sizes)) || any(block_sizes < minimum_each)) {
    stop("all named sampling blocks must meet the minimum", call. = FALSE)
  }
  if (sum(block_sizes) <= total_limit) return(block_sizes)
  if (length(block_sizes) * minimum_each > total_limit) {
    stop("total limit cannot preserve block minima", call. = FALSE)
  }
  quota <- rep(as.integer(minimum_each), length(block_sizes))
  remaining <- as.integer(total_limit - sum(quota))
  while (remaining > 0L) {
    open <- which(quota < block_sizes)
    if (!length(open)) break
    add <- min(remaining, length(open))
    ordered <- open[order(quota[open], names(block_sizes)[open], method = "radix")]
    quota[ordered[seq_len(add)]] <- quota[ordered[seq_len(add)]] + 1L
    remaining <- remaining - add
  }
  names(quota) <- names(block_sizes)
  quota
}
