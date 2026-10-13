copy_combat_io_fixture <- function() {
  source_dir <- system.file(
    "extdata", "examples", "exact_toy",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  output_dir <- tempfile("ab-combat-io-")
  dir.create(output_dir)
  files <- list.files(source_dir, full.names = TRUE)
  copied <- file.copy(files, output_dir)
  if (!all(copied)) stop("failed to copy the packaged toy input")
  output_dir
}

read_combat_io_fixture <- function(path) {
  getFromNamespace("read_input_directory", "AnnotationBreakpoint")(path)
}

validate_combat_io_fixture <- function(raw) {
  do.call(
    validate_inputs,
    raw[c(
      "counts", "cells", "specimens", "allowed_labels", "claim",
      "constraints", "analysis_spec"
    )]
  )
}

read_combat_io_table <- function(path) {
  utils::read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"
  )
}

write_combat_io_table <- function(value, path) {
  utils::write.table(
    value, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
}

set_combat_constraint_token <- function(path, field, token) {
  constraints_path <- file.path(path, "constraints.tsv")
  constraints <- read_combat_io_table(constraints_path)
  constraints[[field]] <- token
  write_combat_io_table(constraints, constraints_path)
}

test_that("toy counts.tsv behavior is unchanged", {
  path <- copy_combat_io_fixture()
  raw <- read_combat_io_fixture(path)
  count_table <- read_combat_io_table(file.path(path, "counts.tsv"))
  expected <- as.matrix(count_table[, -1L, drop = FALSE])
  storage.mode(expected) <- "double"
  rownames(expected) <- count_table$gene_id

  expect_s4_class(raw$counts, "dgCMatrix")
  expect_equal(as.matrix(raw$counts), expected, tolerance = 0)
  expect_identical(
    basename(attr(raw, "source_paths")[["counts"]]),
    "counts.tsv"
  )
  expect_identical(
    unname(attr(raw, "source_hashes")[["counts"]]),
    unname(getFromNamespace(".ab_file_hash", "AnnotationBreakpoint")(
      file.path(path, "counts.tsv")
    ))
  )
  expect_identical(
    names(raw$cells),
    names(read_combat_io_table(file.path(path, "cells.tsv")))
  )
  expect_true(all(vapply(raw$constraints, is.character, logical(1L))))
  expect_true(validate_combat_io_fixture(raw)$ok)
})

test_that("normal input requires exactly one count representation", {
  both <- copy_combat_io_fixture()
  baseline <- read_combat_io_fixture(both)
  saveRDS(baseline$counts, file.path(both, "counts.rds"))
  expect_error(
    read_combat_io_fixture(both),
    "exactly one of counts.tsv or counts.rds, not both",
    fixed = TRUE
  )

  missing <- copy_combat_io_fixture()
  unlink(file.path(missing, "counts.tsv"))
  expect_error(
    read_combat_io_fixture(missing),
    "counts.tsv or counts.rds",
    fixed = TRUE
  )
})

test_that("counts.rds accepts sparse matrices without densifying", {
  path <- copy_combat_io_fixture()
  baseline <- read_combat_io_fixture(path)
  saveRDS(baseline$counts, file.path(path, "counts.rds"))
  unlink(file.path(path, "counts.tsv"))

  observed <- read_combat_io_fixture(path)
  expect_s4_class(observed$counts, "dgCMatrix")
  expect_identical(observed$counts, baseline$counts)
  expect_identical(
    basename(attr(observed, "source_paths")[["counts"]]),
    "counts.rds"
  )
  expect_identical(
    unname(attr(observed, "source_hashes")[["counts"]]),
    unname(getFromNamespace(".ab_file_hash", "AnnotationBreakpoint")(
      file.path(path, "counts.rds")
    ))
  )
  expect_true(validate_combat_io_fixture(observed)$ok)
})

test_that("counts.rds rejects dense and unnamed payloads", {
  dense_path <- copy_combat_io_fixture()
  dense <- read_combat_io_fixture(dense_path)$counts
  unlink(file.path(dense_path, "counts.tsv"))
  saveRDS(as.matrix(dense), file.path(dense_path, "counts.rds"))
  expect_error(
    read_combat_io_fixture(dense_path),
    "sparse Matrix object",
    fixed = TRUE
  )

  unnamed_path <- copy_combat_io_fixture()
  unnamed <- read_combat_io_fixture(unnamed_path)$counts
  rownames(unnamed) <- NULL
  unlink(file.path(unnamed_path, "counts.tsv"))
  saveRDS(unnamed, file.path(unnamed_path, "counts.rds"))
  expect_error(
    read_combat_io_fixture(unnamed_path),
    "nonempty row and column names",
    fixed = TRUE
  )
})

test_that("JSON cap objects decode to named numeric maps and validate", {
  path <- copy_combat_io_fixture()
  specimens <- read_combat_io_table(file.path(path, "specimens.tsv"))
  sample_caps <- stats::setNames(c(1, 2, 1, 2), specimens$sample_id)
  donor_caps <- stats::setNames(c(2, 1, 2, 1), specimens$donor_id)
  set_combat_constraint_token(
    path, "max_changes_per_sample",
    jsonlite::toJSON(as.list(sample_caps), auto_unbox = TRUE)
  )
  set_combat_constraint_token(
    path, "max_changes_per_donor",
    jsonlite::toJSON(as.list(donor_caps), auto_unbox = TRUE)
  )

  raw <- read_combat_io_fixture(path)
  observed_sample <- raw$constraints$max_changes_per_sample[[1L]]
  observed_donor <- raw$constraints$max_changes_per_donor[[1L]]
  expect_type(observed_sample, "double")
  expect_type(observed_donor, "double")
  expect_identical(observed_sample, sample_caps)
  expect_identical(observed_donor, donor_caps)

  validation <- validate_combat_io_fixture(raw)
  expect_true(validation$ok)
  expect_identical(
    validation$data$constraints$max_changes_per_sample,
    sample_caps[validation$data$sample_ids]
  )
  expect_identical(
    validation$data$constraints$max_changes_per_donor,
    donor_caps[validation$data$donor_ids]
  )
})

test_that("malformed, array, duplicate, and wrong-name cap encodings fail", {
  malformed <- copy_combat_io_fixture()
  set_combat_constraint_token(
    malformed, "max_changes_per_sample", '{"S_C1":1,'
  )
  expect_error(
    read_combat_io_fixture(malformed),
    "max_changes_per_sample must be valid JSON",
    fixed = TRUE
  )

  duplicate <- copy_combat_io_fixture()
  set_combat_constraint_token(
    duplicate, "max_changes_per_sample", '{"S_C1":1,"S_C1":2}'
  )
  expect_error(
    read_combat_io_fixture(duplicate),
    "unique nonempty ID keys",
    fixed = TRUE
  )

  array <- copy_combat_io_fixture()
  set_combat_constraint_token(
    array, "max_changes_per_sample", "[1,2,1,2]"
  )
  array_validation <- validate_combat_io_fixture(
    read_combat_io_fixture(array)
  )
  expect_false(array_validation$ok)
  expect_match(
    paste(array_validation$errors, collapse = " | "),
    "max_changes_per_sample values must be nonnegative integers",
    fixed = TRUE
  )

  wrong_name <- copy_combat_io_fixture()
  set_combat_constraint_token(
    wrong_name, "max_changes_per_sample",
    '{"S_C1":1,"S_C2":2,"S_K1":1,"WRONG":2}'
  )
  wrong_name_validation <- validate_combat_io_fixture(
    read_combat_io_fixture(wrong_name)
  )
  expect_false(wrong_name_validation$ok)
  expect_match(
    paste(wrong_name_validation$errors, collapse = " | "),
    "max_changes_per_sample names must match the registered IDs exactly",
    fixed = TRUE
  )
})

test_that("optional adt_ambiguity values decode to logical TRUE and FALSE", {
  path <- copy_combat_io_fixture()
  cells_path <- file.path(path, "cells.tsv")
  cells <- read_combat_io_table(cells_path)
  expected <- rep(c(TRUE, FALSE), length.out = nrow(cells))
  cells$adt_ambiguity <- ifelse(expected, "TRUE", "FALSE")
  write_combat_io_table(cells, cells_path)

  raw <- read_combat_io_fixture(path)
  expect_type(raw$cells$adt_ambiguity, "logical")
  expect_identical(raw$cells$adt_ambiguity, expected)
  validation <- validate_combat_io_fixture(raw)
  expect_true(validation$ok)
  expect_type(validation$data$cells$adt_ambiguity, "logical")
  expect_identical(validation$data$cells$adt_ambiguity, expected)

  invalid_path <- copy_combat_io_fixture()
  invalid_cells_path <- file.path(invalid_path, "cells.tsv")
  invalid_cells <- read_combat_io_table(invalid_cells_path)
  invalid_cells$adt_ambiguity <- rep("FALSE", nrow(invalid_cells))
  invalid_cells$adt_ambiguity[[1L]] <- "unknown"
  write_combat_io_table(invalid_cells, invalid_cells_path)
  expect_error(
    read_combat_io_fixture(invalid_path),
    "adt_ambiguity must contain only TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("prefix sizes reject values above the R integer range", {
  prefix_sizes <- getFromNamespace(
    ".ab_proxy_prefix_sizes", "AnnotationBreakpoint"
  )
  expect_error(
    prefix_sizes(c(1, as.double(.Machine$integer.max) + 1)),
    "strictly increasing list of unique positive integers",
    fixed = TRUE
  )
})
