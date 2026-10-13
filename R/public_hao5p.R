.ab_hao5p_files <- function(raw_dir) {
  relative <- c(
    metadata = "GSE164378_sc.meta.data_5P.csv.gz",
    rna_matrix = "GSM5008740_RNA_5P-matrix.mtx.gz",
    rna_barcodes = "GSM5008740_RNA_5P-barcodes.tsv.gz",
    rna_features = "GSM5008740_RNA_5P-features.tsv.gz",
    adt_matrix = "GSM5008741_ADT_5P-matrix.mtx.gz",
    adt_barcodes = "GSM5008741_ADT_5P-barcodes.tsv.gz",
    adt_features = "GSM5008741_ADT_5P-features.tsv.gz"
  )
  paths <- setNames(file.path(raw_dir, unname(relative)), names(relative))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Hao5P raw directory is missing: ",
         paste(basename(missing), collapse = ", "), call. = FALSE)
  }
  paths
}

.ab_hao5p_matrix_header <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  banner <- readLines(connection, n = 1L, warn = FALSE)
  expected_banner <- "%%MatrixMarket matrix coordinate integer general"
  if (!identical(banner, expected_banner)) {
    stop("unsupported Hao5P Matrix Market banner in ", basename(path),
         call. = FALSE)
  }
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("missing Matrix Market dimensions", call. = FALSE)
    if (!startsWith(line, "%")) break
  }
  tokens <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  values <- suppressWarnings(as.numeric(tokens))
  if (length(values) != 3L || any(!is.finite(values)) ||
      any(values < 0) || any(values != floor(values))) {
    stop("invalid Matrix Market dimensions", call. = FALSE)
  }
  list(
    banner = banner,
    n_features = as.integer(values[[1L]]),
    n_cells = as.integer(values[[2L]]),
    declared_entries = values[[3L]]
  )
}

.ab_hao5p_read_vector <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  values <- readLines(connection, warn = FALSE)
  values <- sub("\r$", "", enc2utf8(values))
  if (!length(values) || anyNA(values) || any(!nzchar(values))) {
    stop("identifier file contains missing or empty rows: ", basename(path),
         call. = FALSE)
  }
  values
}

.ab_hao5p_read_features <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  features <- utils::read.delim(
    connection, header = FALSE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, check.names = FALSE,
    col.names = c("feature_id", "feature_name", "declared_feature_type")
  )
  if (ncol(features) != 3L || !nrow(features) ||
      anyNA(features$feature_id) || any(!nzchar(features$feature_id))) {
    stop("feature file must contain three nonempty identifier columns: ",
         basename(path), call. = FALSE)
  }
  features[] <- lapply(features, enc2utf8)
  features
}

.ab_hao5p_triplet_qc <- function(matrix, chunk_size = 1000000L) {
  if (!inherits(matrix, "TsparseMatrix")) {
    stop("Matrix Market reader did not return a triplet sparse matrix",
         call. = FALSE)
  }
  entries <- length(matrix@x)
  if (!entries) {
    return(list(
      coordinate_order = "column_then_row_strict",
      coordinate_order_valid = TRUE,
      duplicate_coordinate_count = 0,
      finite_entries = TRUE,
      nonnegative_entries = TRUE,
      integer_entries = TRUE,
      explicit_zero_entries = 0
    ))
  }
  ordered <- TRUE
  duplicate_count <- 0
  finite <- TRUE
  nonnegative <- TRUE
  integer <- TRUE
  zeros <- 0
  previous_i <- NA_integer_
  previous_j <- NA_integer_
  starts <- seq.int(1, entries, by = as.integer(chunk_size))
  for (start in starts) {
    end <- min(entries, start + as.integer(chunk_size) - 1L)
    index <- start:end
    current_i <- matrix@i[index]
    current_j <- matrix@j[index]
    current_x <- matrix@x[index]
    prior_i <- c(previous_i, current_i[-length(current_i)])
    prior_j <- c(previous_j, current_j[-length(current_j)])
    compare <- !is.na(prior_j)
    if (any(current_j[compare] < prior_j[compare] |
            (current_j[compare] == prior_j[compare] &
             current_i[compare] < prior_i[compare]))) {
      ordered <- FALSE
    }
    duplicate_count <- duplicate_count + sum(
      compare & current_j == prior_j & current_i == prior_i
    )
    finite <- finite && all(is.finite(current_x))
    nonnegative <- nonnegative && all(current_x >= 0)
    integer <- integer && all(current_x == floor(current_x))
    zeros <- zeros + sum(current_x == 0)
    previous_i <- current_i[[length(current_i)]]
    previous_j <- current_j[[length(current_j)]]
  }
  list(
    coordinate_order = "column_then_row_strict",
    coordinate_order_valid = ordered,
    duplicate_coordinate_count = as.numeric(duplicate_count),
    finite_entries = finite,
    nonnegative_entries = nonnegative,
    integer_entries = integer,
    explicit_zero_entries = as.numeric(zeros)
  )
}

.ab_read_hao5p_sparse_assay <- function(raw_dir, assay = c("RNA", "ADT"),
                                        expected_sha256 = NULL,
                                        verify_hash = FALSE) {
  assay <- match.arg(assay)
  paths <- .ab_hao5p_files(raw_dir)
  prefix <- tolower(assay)
  matrix_path <- unname(paths[[paste0(prefix, "_matrix")]])
  barcode_path <- unname(paths[[paste0(prefix, "_barcodes")]])
  feature_path <- unname(paths[[paste0(prefix, "_features")]])
  header <- .ab_hao5p_matrix_header(matrix_path)
  barcodes <- .ab_hao5p_read_vector(barcode_path)
  features <- .ab_hao5p_read_features(feature_path)
  if (anyDuplicated(barcodes)) {
    stop(assay, " barcode identifiers are duplicated", call. = FALSE)
  }
  if (anyDuplicated(features$feature_id)) {
    stop(assay, " feature identifiers are duplicated", call. = FALSE)
  }
  if (header$n_cells != length(barcodes) ||
      header$n_features != nrow(features)) {
    stop(assay, " Matrix Market dimensions do not match identifier files",
         call. = FALSE)
  }
  observed_hash <- NA_character_
  if (isTRUE(verify_hash)) {
    observed_hash <- .ab_file_hash(matrix_path)
    if (!is.null(expected_sha256) &&
        !identical(observed_hash, as.character(expected_sha256))) {
      stop(assay, " matrix SHA-256 differs from the accepted source hash",
           call. = FALSE)
    }
  }

  connection <- gzfile(matrix_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  started <- unname(proc.time()[["elapsed"]])
  triplet <- Matrix::readMM(connection)
  read_seconds <- unname(proc.time()[["elapsed"]] - started)
  close(connection)
  on.exit(NULL, add = FALSE)
  if (!identical(as.integer(dim(triplet)),
                 c(header$n_features, header$n_cells))) {
    stop(assay, " parsed matrix dimensions differ from its header", call. = FALSE)
  }
  if (length(triplet@x) != header$declared_entries) {
    stop(assay, " parsed entry count differs from the declared count",
         call. = FALSE)
  }
  entry_qc <- .ab_hao5p_triplet_qc(triplet)
  if (!isTRUE(entry_qc$coordinate_order_valid)) {
    stop(assay,
         " coordinates are not in canonical column/row order; low-memory duplicate validation cannot certify this file",
         call. = FALSE)
  }
  if (entry_qc$duplicate_coordinate_count != 0 ||
      !isTRUE(entry_qc$finite_entries) ||
      !isTRUE(entry_qc$nonnegative_entries) ||
      !isTRUE(entry_qc$integer_entries)) {
    stop(assay, " matrix has duplicate or invalid count entries", call. = FALSE)
  }
  started_coercion <- unname(proc.time()[["elapsed"]])
  matrix <- methods::as(triplet, "CsparseMatrix")
  matrix <- Matrix::drop0(matrix)
  coercion_seconds <- unname(proc.time()[["elapsed"]] - started_coercion)
  rm(triplet)
  rownames(matrix) <- features$feature_id
  colnames(matrix) <- barcodes
  stored_entries <- length(matrix@x)
  expected_stored <- header$declared_entries - entry_qc$explicit_zero_entries
  if (stored_entries != expected_stored) {
    stop(assay, " compressed sparse entry count changed unexpectedly",
         call. = FALSE)
  }
  cell_totals <- as.numeric(Matrix::colSums(matrix))
  names(cell_totals) <- barcodes
  total_slots <- as.numeric(header$n_features) * as.numeric(header$n_cells)
  declared_feature_types <- sort(unique(features$declared_feature_type),
                                 method = "radix")
  adt_anomaly <- identical(assay, "ADT") &&
    identical(declared_feature_types, "Gene Expression")
  qc <- list(
    schema_version = "1.0.0",
    assay = assay,
    assay_identity_basis = if (identical(assay, "ADT")) {
      "GSM5008741 ADT accession and filename; not the generic feature-type field"
    } else {
      "GSM5008740 RNA accession and filename"
    },
    matrix_path = normalizePath(matrix_path, mustWork = TRUE),
    matrix_sha256 = observed_hash,
    matrix_sha256_verified = isTRUE(verify_hash),
    dimensions = c(features = header$n_features, cells = header$n_cells),
    declared_entries = header$declared_entries,
    stored_nonzero_entries = as.numeric(stored_entries),
    declared_sparsity = 1 - header$declared_entries / total_slots,
    stored_sparsity = 1 - stored_entries / total_slots,
    duplicate_barcode_count = 0L,
    duplicate_feature_id_count = 0L,
    duplicate_feature_name_count = sum(duplicated(features$feature_name)),
    identifiers_match_dimensions = TRUE,
    entry_validation = entry_qc,
    cell_total_minimum = min(cell_totals),
    cell_total_maximum = max(cell_totals),
    zero_total_cell_count = sum(cell_totals == 0),
    declared_feature_type_values = declared_feature_types,
    feature_type_anomaly = adt_anomaly,
    feature_type_anomaly_record = if (adt_anomaly) {
      paste(
        "The ADT feature file labels its third column Gene Expression.",
        "Assay identity is retained from the verified ADT accession/source,",
        "and ADT counts remain separate from RNA normalization."
      )
    } else NA_character_,
    read_seconds = read_seconds,
    sparse_coercion_seconds = coercion_seconds
  )
  structure(list(
    assay = assay,
    matrix = matrix,
    barcodes = barcodes,
    features = features,
    cell_totals = cell_totals,
    qc = qc
  ), class = c("ab_hao5p_assay", "list"))
}

.ab_read_hao5p_metadata <- function(path, expected_barcodes = NULL) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  metadata <- utils::read.csv(
    connection, header = TRUE, row.names = 1L, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata$cell_id <- enc2utf8(rownames(metadata))
  rownames(metadata) <- NULL
  required <- c(
    "nCount_ADT", "nFeature_ADT", "nCount_RNA", "nFeature_RNA",
    "lane", "donor", "time", "celltype.l1", "celltype.l2", "celltype.l3",
    "Batch", "cell_id"
  )
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop("Hao5P metadata is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(metadata$cell_id) || anyNA(metadata$cell_id) ||
      any(!nzchar(metadata$cell_id))) {
    stop("Hao5P metadata cell identifiers are invalid or duplicated",
         call. = FALSE)
  }
  if (!is.null(expected_barcodes) &&
      !identical(metadata$cell_id, expected_barcodes)) {
    stop("metadata and matrix barcode order differ", call. = FALSE)
  }
  for (name in c("nCount_ADT", "nFeature_ADT", "nCount_RNA", "nFeature_RNA")) {
    metadata[[name]] <- suppressWarnings(as.numeric(metadata[[name]]))
    if (any(!is.finite(metadata[[name]]) | metadata[[name]] < 0)) {
      stop("metadata QC field is invalid: ", name, call. = FALSE)
    }
  }
  for (name in c("lane", "donor", "time", "celltype.l1", "celltype.l2",
                 "celltype.l3", "Batch")) {
    metadata[[name]] <- enc2utf8(as.character(metadata[[name]]))
    if (anyNA(metadata[[name]]) || any(!nzchar(metadata[[name]]))) {
      stop("metadata identity field is empty: ", name, call. = FALSE)
    }
  }
  metadata
}

.ab_hao5p_selection_hash <- function(cell_id, seed,
                                     namespace = "hao5p-cell-selection-v1") {
  vapply(cell_id, function(id) digest::digest(
    paste(namespace, as.integer(seed), id, sep = "|"),
    algo = "sha256", serialize = FALSE
  ), character(1L))
}

.ab_select_hao5p_day0_monocytes <- function(metadata, seed,
                                            cap_per_population_donor = 50L,
                                            namespace = "hao5p-cell-selection-v1",
                                            time = "0",
                                            populations = c("CD14 Mono", "CD16 Mono"),
                                            donors = paste0("P", 1:8),
                                            calibration_donors = c("P1", "P2"),
                                            evaluation_donors = paste0("P", 3:8)) {
  if (!is.data.frame(metadata) ||
      !all(c("cell_id", "donor", "time", "lane", "celltype.l2") %in%
           names(metadata))) {
    stop("metadata lacks Hao5P selection fields", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed != floor(seed)) {
    stop("selection seed must be one finite integer", call. = FALSE)
  }
  if (length(cap_per_population_donor) != 1L ||
      is.na(cap_per_population_donor) || cap_per_population_donor < 1L ||
      cap_per_population_donor != floor(cap_per_population_donor)) {
    stop("selection cap must be a positive integer", call. = FALSE)
  }
  if (length(namespace) != 1L || is.na(namespace) || !nzchar(namespace)) {
    stop("selection hash namespace must be one nonempty string", call. = FALSE)
  }
  if (length(populations) != 2L || anyDuplicated(populations) ||
      any(!nzchar(populations)) || length(donors) < 2L || anyDuplicated(donors) ||
      !setequal(c(calibration_donors, evaluation_donors), donors) ||
      length(intersect(calibration_donors, evaluation_donors))) {
    stop("invalid registered selection populations or donor split", call. = FALSE)
  }
  baseline <- metadata$time == as.character(time)
  eligible <- baseline & metadata$donor %in% donors &
    metadata$celltype.l2 %in% populations
  coverage <- table(
    factor(metadata$donor[eligible], levels = donors),
    factor(metadata$celltype.l2[eligible], levels = populations)
  )
  if (any(coverage < cap_per_population_donor)) {
    stop("at least one day-0 donor/population group is below the registered cap",
         call. = FALSE)
  }
  provenance <- metadata
  provenance$matrix_column_index <- seq_len(nrow(metadata))
  provenance$eligible_day0_pair <- eligible
  provenance$selection_hash <- NA_character_
  provenance$selection_rank <- NA_integer_
  provenance$selected <- FALSE
  provenance$sample_id <- ifelse(
    baseline, paste0("Hao5P_day0_", metadata$donor), NA_character_
  )
  provenance$donor_role <- ifelse(
    metadata$donor %in% calibration_donors, "calibration",
    ifelse(metadata$donor %in% evaluation_donors, "evaluation", NA_character_)
  )
  eligible_index <- which(eligible)
  provenance$selection_hash[eligible_index] <- .ab_hao5p_selection_hash(
    provenance$cell_id[eligible_index], seed, namespace
  )
  selected <- integer()
  for (donor in donors) {
    for (population in populations) {
      index <- which(eligible & metadata$donor == donor &
                       metadata$celltype.l2 == population)
      order_index <- order(provenance$selection_hash[index],
                           provenance$cell_id[index], method = "radix")
      ranked <- index[order_index]
      provenance$selection_rank[ranked] <- seq_along(ranked)
      chosen <- ranked[seq_len(cap_per_population_donor)]
      provenance$selected[chosen] <- TRUE
      selected <- c(selected, chosen)
    }
  }
  selected <- selected[order(
    match(metadata$donor[selected], donors),
    match(metadata$celltype.l2[selected], populations),
    provenance$selection_rank[selected], method = "radix"
  )]
  specimens <- do.call(rbind, lapply(donors, function(donor) {
    donor_baseline <- metadata[baseline & metadata$donor == donor, , drop = FALSE]
    data.frame(
      sample_id = paste0("Hao5P_day0_", donor),
      donor_id = donor,
      source_time = "0",
      donor_role = if (donor %in% calibration_donors) "calibration" else "evaluation",
      pooled_technical_lanes = paste(
        sort(unique(donor_baseline$lane), method = "radix"), collapse = ";"
      ),
      technical_lane_count = length(unique(donor_baseline$lane)),
      biological_specimen_count = 1L,
      stringsAsFactors = FALSE
    )
  }))
  selected_cells <- provenance[selected, , drop = FALSE]
  if (!all(table(selected_cells$donor, selected_cells$celltype.l2) ==
           cap_per_population_donor) ||
      nrow(specimens) != 8L || anyDuplicated(specimens$donor_id)) {
    stop("deterministic selection violated its donor/population invariant",
         call. = FALSE)
  }
  structure(list(
    selected_cell_ids = selected_cells$cell_id,
    selected_cells = selected_cells,
    specimens = specimens,
    full_cell_provenance = provenance,
    coverage_before_cap = unclass(coverage),
    selection_seed = as.integer(seed),
    selection_hash_namespace = namespace,
    cap_per_population_donor = as.integer(cap_per_population_donor),
    technical_lane_policy = paste(
      "Pool lanes only within the same published donor and day-0 time point;",
      "retain lane in cell provenance; one specimen per real donor."
    )
  ), class = c("ab_hao5p_selection", "list"))
}

.ab_hao5p_source_hashes <- function(raw_dir) {
  manifest_path <- file.path(raw_dir, "acquisition_manifest.json")
  if (!file.exists(manifest_path)) {
    stop("accepted acquisition_manifest.json is missing", call. = FALSE)
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  entries <- manifest$files
  setNames(
    vapply(entries, function(entry) as.character(entry$sha256), character(1L)),
    vapply(entries, function(entry) as.character(entry$file), character(1L))
  )
}

.ab_hao5p_inspect <- function(raw_dir, verify_hashes = TRUE) {
  raw_dir <- normalizePath(raw_dir, mustWork = TRUE)
  paths <- .ab_hao5p_files(raw_dir)
  accepted_hashes <- .ab_hao5p_source_hashes(raw_dir)
  started <- unname(proc.time()[["elapsed"]])
  observed_hashes <- if (isTRUE(verify_hashes)) {
    setNames(vapply(paths, .ab_file_hash, character(1L)), basename(paths))
  } else setNames(rep(NA_character_, length(paths)), basename(paths))
  if (isTRUE(verify_hashes)) {
    expected <- accepted_hashes[names(observed_hashes)]
    if (anyNA(expected) || !identical(unname(observed_hashes), unname(expected))) {
      stop("raw Hao5P files differ from accepted acquisition hashes",
           call. = FALSE)
    }
  }
  rna <- .ab_read_hao5p_sparse_assay(raw_dir, "RNA")
  adt <- .ab_read_hao5p_sparse_assay(raw_dir, "ADT")
  rna$qc$matrix_sha256 <- observed_hashes[[basename(paths[["rna_matrix"]])]]
  rna$qc$matrix_sha256_verified <- isTRUE(verify_hashes)
  adt$qc$matrix_sha256 <- observed_hashes[[basename(paths[["adt_matrix"]])]]
  adt$qc$matrix_sha256_verified <- isTRUE(verify_hashes)
  if (!identical(rna$barcodes, adt$barcodes)) {
    stop("RNA and ADT barcode orders differ", call. = FALSE)
  }
  metadata <- .ab_read_hao5p_metadata(paths[["metadata"]], rna$barcodes)
  rna_difference <- rna$cell_totals - metadata$nCount_RNA
  adt_difference <- adt$cell_totals - metadata$nCount_ADT
  total_audit <- list(
    rna_reported_total_match_count = sum(rna_difference == 0),
    rna_reported_total_mismatch_count = sum(rna_difference != 0),
    rna_maximum_absolute_total_difference = max(abs(rna_difference)),
    adt_reported_total_match_count = sum(adt_difference == 0),
    adt_reported_total_mismatch_count = sum(adt_difference != 0),
    adt_maximum_absolute_total_difference = max(abs(adt_difference)),
    caveat = paste(
      "Metadata totals are audited against every supplied feature row.",
      "Any mismatch is retained because metadata may reflect an upstream",
      "feature/QC convention rather than the stored feature universe."
    )
  )
  elapsed <- unname(proc.time()[["elapsed"]] - started)
  structure(list(
    raw_dir = raw_dir,
    paths = paths,
    accepted_sha256 = accepted_hashes,
    observed_sha256 = observed_hashes,
    hashes_verified = isTRUE(verify_hashes),
    rna = rna,
    adt = adt,
    metadata = metadata,
    qc = list(
      schema_version = "1.0.0",
      rna = rna$qc,
      adt = adt$qc,
      ordered_barcode_match_rna_adt_metadata = TRUE,
      total_audit = total_audit,
      assays_kept_separate = TRUE,
      endpoint_analysis_performed = FALSE,
      total_inspection_seconds = elapsed
    )
  ), class = c("ab_hao5p_inspection", "list"))
}

.ab_hao5p_select <- function(inspected, design) {
  if (!inherits(inspected, "ab_hao5p_inspection")) {
    stop("inspected must be an ab_hao5p_inspection", call. = FALSE)
  }
  if (!is.list(design) || is.null(design$selection) ||
      !identical(as.character(design$benchmark_id),
                 "hao5p-day0-cd14-cd16-semisynthetic-v1")) {
    stop("design must be the registered Hao5P benchmark mapping",
         call. = FALSE)
  }
  selection_spec <- design$selection
  sampling <- selection_spec$sampling
  if (!identical(as.character(sampling$algorithm),
                 "ascending_sha256_of_namespace_seed_barcode")) {
    stop("unsupported deterministic selection algorithm", call. = FALSE)
  }
  selection <- .ab_select_hao5p_day0_monocytes(
    inspected$metadata,
    seed = as.integer(sampling$seed),
    cap_per_population_donor = as.integer(
      selection_spec$cells_per_population_per_donor
    ),
    namespace = as.character(sampling$namespace),
    time = as.character(selection_spec$time),
    populations = as.character(unlist(selection_spec$populations)),
    donors = as.character(unlist(selection_spec$donors)),
    calibration_donors = as.character(unlist(selection_spec$calibration_donors)),
    evaluation_donors = as.character(unlist(selection_spec$evaluation_donors))
  )
  column_index <- match(selection$selected_cell_ids, inspected$rna$barcodes)
  if (anyNA(column_index)) stop("selected barcodes are absent from matrices",
                                call. = FALSE)
  selected_metadata <- selection$selected_cells
  selected_metadata$rna_total_from_matrix <- inspected$rna$cell_totals[column_index]
  selected_metadata$adt_total_from_matrix <- inspected$adt$cell_totals[column_index]
  structure(list(
    schema_version = "1.0.0",
    snapshot_id = paste0(design$benchmark_id, "-selection-v1"),
    benchmark_id = as.character(design$benchmark_id),
    source = list(
      accession = as.character(design$source$study),
      rna_accession = as.character(design$source$assay_policy$rna_accession),
      adt_accession = as.character(design$source$assay_policy$adt_accession),
      accepted_sha256 = inspected$accepted_sha256,
      observed_sha256 = inspected$observed_sha256,
      hashes_verified = inspected$hashes_verified
    ),
    rna_counts = inspected$rna$matrix[, column_index, drop = FALSE],
    adt_counts = inspected$adt$matrix[, column_index, drop = FALSE],
    rna_features = inspected$rna$features,
    adt_features = inspected$adt$features,
    selected_cells = selected_metadata,
    specimens = selection$specimens,
    full_cell_provenance = selection$full_cell_provenance,
    selection = selection[c(
      "coverage_before_cap", "selection_seed", "selection_hash_namespace",
      "cap_per_population_donor", "technical_lane_policy"
    )],
    qc = inspected$qc
  ), class = c("ab_hao5p_snapshot", "list"))
}

.ab_prepare_hao5p_snapshot <- function(raw_dir, selection_seed,
                                       cap_per_population_donor = 50L,
                                       verify_matrix_hashes = FALSE) {
  raw_dir <- normalizePath(raw_dir, mustWork = TRUE)
  paths <- .ab_hao5p_files(raw_dir)
  hashes <- .ab_hao5p_source_hashes(raw_dir)
  rna_expected <- hashes[[basename(paths[["rna_matrix"]])]]
  adt_expected <- hashes[[basename(paths[["adt_matrix"]])]]
  rna <- .ab_read_hao5p_sparse_assay(
    raw_dir, "RNA", rna_expected, verify_hash = verify_matrix_hashes
  )
  adt <- .ab_read_hao5p_sparse_assay(
    raw_dir, "ADT", adt_expected, verify_hash = verify_matrix_hashes
  )
  if (!identical(rna$barcodes, adt$barcodes)) {
    stop("RNA and ADT barcode orders differ", call. = FALSE)
  }
  metadata <- .ab_read_hao5p_metadata(paths[["metadata"]], rna$barcodes)
  rna_difference <- rna$cell_totals - metadata$nCount_RNA
  adt_difference <- adt$cell_totals - metadata$nCount_ADT
  total_audit <- list(
    rna_reported_total_match_count = sum(rna_difference == 0),
    rna_reported_total_mismatch_count = sum(rna_difference != 0),
    rna_maximum_absolute_total_difference = max(abs(rna_difference)),
    adt_reported_total_match_count = sum(adt_difference == 0),
    adt_reported_total_mismatch_count = sum(adt_difference != 0),
    adt_maximum_absolute_total_difference = max(abs(adt_difference)),
    caveat = paste(
      "Metadata totals are audited against every supplied feature row.",
      "Any mismatch must be retained because the metadata may reflect an upstream",
      "feature/QC convention rather than the stored matrix feature universe."
    )
  )
  selection <- .ab_select_hao5p_day0_monocytes(
    metadata, selection_seed, cap_per_population_donor
  )
  column_index <- match(selection$selected_cell_ids, rna$barcodes)
  if (anyNA(column_index)) stop("selected barcodes are absent from matrices",
                                call. = FALSE)
  selected_metadata <- selection$selected_cells
  selected_metadata$rna_total_from_matrix <- rna$cell_totals[column_index]
  selected_metadata$adt_total_from_matrix <- adt$cell_totals[column_index]
  structure(list(
    schema_version = "1.0.0",
    snapshot_id = paste0(
      "hao5p-day0-cd14-cd16-cap", as.integer(cap_per_population_donor),
      "-seed", as.integer(selection_seed), "-v1"
    ),
    source = list(
      accession = "GSE164378",
      rna_accession = "GSM5008740",
      adt_accession = "GSM5008741",
      accepted_sha256 = hashes
    ),
    rna_counts = rna$matrix[, column_index, drop = FALSE],
    adt_counts = adt$matrix[, column_index, drop = FALSE],
    rna_features = rna$features,
    adt_features = adt$features,
    selected_cells = selected_metadata,
    specimens = selection$specimens,
    full_cell_provenance = selection$full_cell_provenance,
    selection = selection[c(
      "coverage_before_cap", "selection_seed", "cap_per_population_donor",
      "technical_lane_policy"
    )],
    qc = list(
      rna = rna$qc,
      adt = adt$qc,
      ordered_barcode_match_rna_adt_metadata = TRUE,
      total_audit = total_audit,
      assays_kept_separate = TRUE,
      endpoint_analysis_performed = FALSE
    )
  ), class = c("ab_hao5p_snapshot", "list"))
}

.ab_hao5p_write_tsv <- function(value, path) {
  utils::write.table(
    value, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", eol = "\n", fileEncoding = "UTF-8"
  )
  invisible(path)
}

.ab_hao5p_write_tsv_gz <- function(value, path) {
  connection <- gzfile(path, open = "wt", compression = 6L)
  on.exit(close(connection), add = TRUE)
  utils::write.table(
    value, connection, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", eol = "\n", fileEncoding = "UTF-8"
  )
  invisible(path)
}

.ab_write_hao5p_snapshot <- function(snapshot, output_dir) {
  if (!inherits(snapshot, "ab_hao5p_snapshot")) {
    stop("snapshot must be an ab_hao5p_snapshot", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("output_dir must be one nonempty path", call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    source_snapshot = file.path(output_dir, "source_snapshot.rds"),
    rna_counts = file.path(output_dir, "rna_counts_sparse.rds"),
    adt_counts = file.path(output_dir, "adt_counts_sparse.rds"),
    rna_features = file.path(output_dir, "rna_features.tsv.gz"),
    adt_features = file.path(output_dir, "adt_features.tsv.gz"),
    selected_cells = file.path(output_dir, "selected_cells.tsv.gz"),
    specimens = file.path(output_dir, "specimens.tsv"),
    full_cell_provenance = file.path(output_dir, "full_cell_provenance.tsv.gz"),
    qc = file.path(output_dir, "qc.json")
  )
  saveRDS(snapshot, paths[["source_snapshot"]], compress = "xz", version = 3L)
  saveRDS(snapshot$rna_counts, paths[["rna_counts"]], compress = "xz",
          version = 3L)
  saveRDS(snapshot$adt_counts, paths[["adt_counts"]], compress = "xz",
          version = 3L)
  .ab_hao5p_write_tsv_gz(snapshot$rna_features, paths[["rna_features"]])
  .ab_hao5p_write_tsv_gz(snapshot$adt_features, paths[["adt_features"]])
  .ab_hao5p_write_tsv_gz(snapshot$selected_cells, paths[["selected_cells"]])
  .ab_hao5p_write_tsv(snapshot$specimens, paths[["specimens"]])
  .ab_hao5p_write_tsv_gz(snapshot$full_cell_provenance,
                         paths[["full_cell_provenance"]])
  .ab_write_json(snapshot$qc, paths[["qc"]])
  hashes <- vapply(paths, .ab_file_hash, character(1L))
  manifest_path <- file.path(output_dir, "snapshot_manifest.json")
  manifest <- list(
    schema_version = "1.0.0",
    snapshot_id = snapshot$snapshot_id,
    source = list(
      accession = snapshot$source$accession,
      rna_accession = snapshot$source$rna_accession,
      adt_accession = snapshot$source$adt_accession,
      accepted_sha256 = as.list(snapshot$source$accepted_sha256),
      observed_sha256 = as.list(snapshot$source$observed_sha256),
      hashes_verified = snapshot$source$hashes_verified
    ),
    sparse_storage = TRUE,
    assays_kept_separate = TRUE,
    selected_cell_count = ncol(snapshot$rna_counts),
    real_donor_count = nrow(snapshot$specimens),
    calibration_donors = c("P1", "P2"),
    evaluation_donors = paste0("P", 3:8),
    endpoint_analysis_performed = FALSE,
    artifacts = lapply(names(paths), function(name) list(
      path = basename(paths[[name]]),
      bytes = unname(file.info(paths[[name]])$size),
      sha256 = unname(hashes[[name]])
    ))
  )
  names(manifest$artifacts) <- names(paths)
  .ab_write_json(manifest, manifest_path)
  structure(list(
    output_dir = output_dir,
    manifest = manifest_path,
    paths = paths,
    hashes = c(hashes, snapshot_manifest = .ab_file_hash(manifest_path))
  ), class = c("ab_hao5p_snapshot_artifacts", "list"))
}
