test_that("all twelve frozen assignment effects match independent raw-count calculations", {
  expected_status <- c(
    fragile = "exact_minimum",
    stable = "exhaustively_stable",
    coordinated = "exact_minimum"
  )
  expected_best <- c(fragile = 1, stable = NA, coordinated = 2)
  expected_exact <- c(fragile = 1, stable = NA, coordinated = 2)

  for (case_id in names(expected_status)) {
    raw <- load_raw_fixture(case_id)
    x <- validate_raw_fixture(raw)
    expect_true(x$ok, info = case_id)
    result <- enumerate_breakpoint(x)

    hand_effects <- vapply(
      raw$expected$changes,
      function(changes) independent_beta(raw, assignment_from_changes(raw, changes)),
      numeric(1)
    )
    expect_equal(unname(hand_effects), unname(raw$expected$beta), tolerance = 1e-12, info = case_id)
    expect_identical(result$assignments$changes, raw$expected$changes, info = case_id)
    expect_equal(unname(result$assignments$beta), unname(hand_effects), tolerance = 1e-12, info = case_id)
    expect_identical(result$assignments$cost, raw$expected$cost, info = case_id)
    expect_identical(result$assignments$feasible, rep(TRUE, 4L), info = case_id)
    expect_identical(result$assignments$overturn, raw$expected$overturn, info = case_id)

    expect_identical(result$total_assignment_count, 4L, info = case_id)
    expect_identical(result$feasible_assignment_count, 4L, info = case_id)
    expect_identical(result$infeasible_assignment_count, 0L, info = case_id)
    expect_identical(result$evaluable_assignment_count, 4L, info = case_id)
    expect_identical(result$unevaluable_assignment_count, 0L, info = case_id)
    expect_identical(result$status, unname(expected_status[[case_id]]), info = case_id)
    expect_equal(result$best_verified_cost, unname(expected_best[[case_id]]), info = case_id)
    expect_equal(result$exact_minimum_cost, unname(expected_exact[[case_id]]), info = case_id)
  }
})

test_that("fragile and coordinated witnesses replay while stable is exhaustive", {
  fragile_raw <- load_raw_fixture("fragile")
  fragile_x <- validate_raw_fixture(fragile_raw)
  fragile <- enumerate_breakpoint(fragile_x)
  expect_identical(fragile$witness$changes, "F_K1_in=A")
  expect_equal(fragile$witness$beta, 0.7854578635213407, tolerance = 1e-12)
  fragile_check <- verify_witness(fragile_x, fragile)
  expect_true(fragile_check$ok)

  stable_raw <- load_raw_fixture("stable")
  stable <- enumerate_breakpoint(validate_raw_fixture(stable_raw))
  expect_identical(stable$status, "exhaustively_stable")
  expect_null(stable$witness)
  expect_equal(min(stable$assignments$beta), 2.5914993511034403, tolerance = 1e-12)

  coordinated_raw <- load_raw_fixture("coordinated")
  coordinated_x <- validate_raw_fixture(coordinated_raw)
  coordinated <- enumerate_breakpoint(coordinated_x)
  expect_identical(coordinated$witness$changes, "CO_K1_in=A;CO_K1_out=B")
  expect_identical(coordinated$witness$cost, 2L)
  expect_equal(coordinated$witness$beta, 0.759702897877574, tolerance = 1e-12)
  expect_true(verify_witness(coordinated_x, coordinated)$ok)

  single_rows <- coordinated$assignments[coordinated$assignments$cost == 1L, , drop = FALSE]
  expect_true(all(single_rows$feasible))
  expect_false(any(single_rows$overturn))
  original <- original_assignment(coordinated_raw)
  witness_assignment <- coordinated$witness$assignment
  moved <- names(original)[original != witness_assignment[names(original)]]
  expect_setequal(moved, c("CO_K1_in", "CO_K1_out"))
  expect_setequal(
    paste0(original[moved], "->", witness_assignment[moved]),
    c("B->A", "A->B")
  )
})

test_that("pseudobulks conserve counts and incremental forward/reverse moves equal full sums", {
  raw <- load_raw_fixture("fragile")
  x <- validate_raw_fixture(raw)
  original <- original_assignment(raw)
  moved <- original
  moved[["F_K1_in"]] <- "A"

  baseline <- build_pseudobulk(x, assignment = original, method = "full")
  full <- build_pseudobulk(x, assignment = moved, method = "full")
  incremental <- build_pseudobulk(
    x,
    method = "incremental",
    base = baseline,
    moves = c(F_K1_in = "A")
  )
  expect_equal(as.matrix(incremental$counts), as.matrix(full$counts), tolerance = 0)
  expect_identical(incremental$assignment, full$assignment)

  reversed <- build_pseudobulk(
    x,
    method = "incremental",
    base = incremental,
    moves = data.frame(
      cell_id = "F_K1_in",
      from_label = "A",
      to_label = "B",
      stringsAsFactors = FALSE
    )
  )
  expect_equal(as.matrix(reversed$counts), as.matrix(baseline$counts), tolerance = 0)
  expect_identical(reversed$assignment, baseline$assignment)

  for (sample_id in raw$specimens$sample_id) {
    raw_cells <- raw$cells$cell_id[raw$cells$sample_id == sample_id]
    raw_total <- Matrix::rowSums(raw$counts[, raw_cells, drop = FALSE])
    pb_columns <- which(full$column_metadata$sample_id == sample_id)
    observed_total <- Matrix::rowSums(full$counts[, pb_columns, drop = FALSE])
    expect_equal(
      unname(observed_total[names(raw_total)]),
      unname(raw_total),
      tolerance = 0
    )
  }
  expect_true(full$all_conserved)
  expect_gt(nrow(full$conservation), 0L)
  expect_true(all(full$conservation$conserved))
})

test_that("row and column reordering preserves effects and canonical witness selection", {
  raw <- load_raw_fixture("fragile")
  baseline <- enumerate_breakpoint(validate_raw_fixture(raw))

  reordered <- raw
  reordered$counts <- reordered$counts[rev(seq_len(nrow(reordered$counts))),
                                        rev(seq_len(ncol(reordered$counts))), drop = FALSE]
  reordered$cells <- reordered$cells[rev(seq_len(nrow(reordered$cells))), , drop = FALSE]
  reordered$specimens <- reordered$specimens[rev(seq_len(nrow(reordered$specimens))), , drop = FALSE]
  reordered$allowed_labels <- reordered$allowed_labels[
    rev(seq_len(nrow(reordered$allowed_labels))), , drop = FALSE
  ]
  observed <- enumerate_breakpoint(validate_raw_fixture(reordered))

  baseline_effects <- stats::setNames(baseline$assignments$beta, baseline$assignments$changes)
  observed_effects <- stats::setNames(observed$assignments$beta, observed$assignments$changes)
  expect_equal(observed_effects[names(baseline_effects)], baseline_effects, tolerance = 1e-12)
  expect_identical(observed$witness$changes, baseline$witness$changes)
  expect_identical(observed$witness_id, baseline$witness_id)
})
