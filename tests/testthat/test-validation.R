test_that("static validation rejects malformed identities, counts, and allowed labels", {
  base <- load_raw_fixture("fragile")

  bad <- base
  bad$cells <- rbind(bad$cells, bad$cells[1L, , drop = FALSE])
  expect_invalid(bad, "duplicate|cell")

  bad <- base
  rownames(bad$counts)[2L] <- rownames(bad$counts)[1L]
  expect_invalid(bad, "feature|duplicate")

  bad <- base
  bad$counts[1L, 1L] <- -1
  expect_invalid(bad, "negative|nonnegative|count")

  bad <- base
  bad$counts[1L, 1L] <- 0.5
  expect_invalid(bad, "integer|fractional|count")

  bad <- base
  bad$counts[1L, 1L] <- Inf
  expect_invalid(bad, "finite|count")

  bad <- base
  bad$allowed_labels <- rbind(bad$allowed_labels, bad$allowed_labels[1L, , drop = FALSE])
  expect_invalid(bad, "duplicate|allowed")

  bad <- base
  bad$allowed_labels$evidence_source[1L] <- ""
  expect_invalid(bad, "evidence|provenance|source")

  bad <- base
  bad$allowed_labels <- bad$allowed_labels[
    !(bad$allowed_labels$cell_id == "F_C1_A" & bad$allowed_labels$allowed_label == "A"),
    , drop = FALSE
  ]
  expect_invalid(bad, "original|allowed")

  bad <- base
  bad$cells$donor_id[bad$cells$cell_id == "F_C1_A"] <- "wrong_donor"
  expect_invalid(bad, "donor|specimen|metadata")

  bad <- base
  bad$cells$donor_id[bad$cells$cell_id == "F_C1_A"] <- NA_character_
  expect_invalid(bad, "donor|missing|nonempty|metadata")
})

test_that("design, replication, budgets, and both-label minima are validated", {
  base <- load_raw_fixture("fragile")

  bad <- base
  bad$specimens$donor_id[bad$specimens$sample_id == "S_C2"] <- "D_C1"
  bad$cells$donor_id[bad$cells$sample_id == "S_C2"] <- "D_C1"
  expect_invalid(bad, "donor|repeated|replication")

  inadequate <- drop_sample(base, "S_C2")
  inadequate$constraints$min_donors_per_condition <- 1L
  expect_invalid(inadequate, "donor|replication|condition")

  for (field in c("max_changes_global", "max_changes_per_sample", "max_changes_per_donor")) {
    bad <- base
    bad$constraints[[field]] <- -1L
    expect_invalid(bad, "budget|changes|nonnegative")

    bad <- base
    bad$constraints[[field]] <- 0.5
    expect_invalid(bad, "budget|changes|integer")
  }

  for (field in c("min_cells_A", "min_cells_B")) {
    bad <- base
    bad$constraints[[field]] <- 2L
    expect_invalid(bad, "minimum|min_cells|feasible")
  }
  for (field in c("min_library_A", "min_library_B")) {
    bad <- base
    bad$constraints[[field]] <- 101L
    expect_invalid(bad, "minimum|min_library|feasible")
  }
})

test_that("assignment input rejects missing, duplicate, conflicting, and impermissible labels", {
  raw <- load_raw_fixture("fragile")
  x <- validate_raw_fixture(raw)
  expect_true(x$ok)
  assignment <- original_assignment(raw)

  expect_error(build_pseudobulk(x, assignment[-1L]), "assignment|cell|missing")
  expect_error(build_pseudobulk(x, c(assignment, F_C1_in = "A")), "duplicate|assignment|cell")

  bad <- assignment
  bad[["F_C1_in"]] <- "C"
  expect_error(build_pseudobulk(x, bad), "allowed|label|assignment")

  conflict <- c(assignment, F_C1_in = "A")
  names(conflict)[length(conflict) - 1L] <- "F_C1_in"
  expect_error(build_pseudobulk(x, conflict), "duplicate|conflict|assignment")
})
