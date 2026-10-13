test_that("empty ambiguity and zero budgets reproduce the baseline", {
  raw <- load_raw_fixture("fragile")
  original_effect <- independent_beta(raw)

  no_ambiguity <- raw
  original <- original_assignment(raw)
  no_ambiguity$allowed_labels <- no_ambiguity$allowed_labels[
    no_ambiguity$allowed_labels$allowed_label == original[no_ambiguity$allowed_labels$cell_id],
    , drop = FALSE
  ]
  no_ambiguity_result <- enumerate_breakpoint(validate_raw_fixture(no_ambiguity))
  expect_identical(no_ambiguity_result$total_assignment_count, 1L)
  expect_identical(no_ambiguity_result$feasible_assignment_count, 1L)
  expect_equal(no_ambiguity_result$original_fit$beta, original_effect, tolerance = 1e-12)
  expect_identical(no_ambiguity_result$status, "exhaustively_stable")

  zero_budget <- raw
  zero_budget$constraints$max_changes_global <- 0L
  zero_budget$constraints$max_changes_per_sample <- 0L
  zero_budget$constraints$max_changes_per_donor <- 0L
  zero_result <- enumerate_breakpoint(validate_raw_fixture(zero_budget))
  expect_identical(zero_result$feasible_assignment_count, 1L)
  expect_equal(zero_result$original_fit$beta, original_effect, tolerance = 1e-12)
  expect_identical(zero_result$status, "exhaustively_stable")
})

test_that("global, sample, and donor budgets independently exclude the cost-two exchange", {
  for (field in c("max_changes_global", "max_changes_per_sample", "max_changes_per_donor")) {
    raw <- load_raw_fixture("coordinated")
    raw$constraints[[field]] <- 1L
    result <- enumerate_breakpoint(validate_raw_fixture(raw))
    pair <- candidate_row(result, "CO_K1_in=A;CO_K1_out=B")
    expect_false(pair$feasible, info = field)
    expect_false(pair$overturn, info = field)
    expect_identical(result$status, "exhaustively_stable", info = field)
  }
})

test_that("candidate minima apply to A and B cells and libraries", {
  inbound <- load_raw_fixture("fragile")
  for (sample_id in c("S_C2", "S_K2")) {
    inbound <- add_fixed_cell(
      inbound,
      paste0("PAD_B_", sample_id),
      sample_id,
      "B",
      c(G_target = 0, G_ref1 = 50, G_ref2 = 50)
    )
  }
  inbound$constraints$min_cells_B <- 2L
  inbound$constraints$min_library_B <- 150L
  inbound_result <- enumerate_breakpoint(validate_raw_fixture(inbound))
  expect_false(candidate_row(inbound_result, "F_K1_in=A")$feasible)

  outbound <- load_raw_fixture("coordinated")
  for (sample_id in c("S_C1", "S_C2", "S_K2")) {
    outbound <- add_fixed_cell(
      outbound,
      paste0("PAD_A_", sample_id),
      sample_id,
      "A",
      c(G_target = 10, G_ref1 = 20, G_ref2 = 20)
    )
  }
  outbound$constraints$min_cells_A <- 2L
  outbound$constraints$min_library_A <- 125L
  outbound_result <- enumerate_breakpoint(validate_raw_fixture(outbound))
  expect_false(candidate_row(outbound_result, "CO_K1_out=B")$feasible)
})

test_that("invalid input, unassessable original, and unmet original claim stay distinct", {
  raw <- load_raw_fixture("fragile")
  invalid <- raw
  invalid$counts[1L, 1L] <- -1
  invalid_result <- enumerate_breakpoint(validate_raw_fixture(invalid))
  expect_identical(invalid_result$status, "invalid_input")

  x <- validate_raw_fixture(raw)
  not_assessable <- enumerate_breakpoint(x, fit_fun = function(...) stop("forced original failure"))
  expect_identical(not_assessable$status, "not_assessable")
  expect_true(is.na(not_assessable$best_verified_cost))
  expect_true(is.na(not_assessable$exact_minimum_cost))

  unmet <- raw
  unmet$claim$delta <- 10
  unmet_result <- enumerate_breakpoint(validate_raw_fixture(unmet))
  expect_identical(unmet_result$status, "original_claim_not_met")
  expect_true(is.na(unmet_result$best_verified_cost))
  expect_true(is.na(unmet_result$exact_minimum_cost))
})

test_that("failed callbacks block exactness and exhaustive stability", {
  coordinated_raw <- load_raw_fixture("coordinated")
  coordinated_x <- validate_raw_fixture(coordinated_raw)
  coordinated_fit <- function(x, assignment = NULL, ...) {
    key <- independent_changes(coordinated_raw, assignment)
    if (identical(key, "CO_K1_in=A")) return(failed_fit())
    fit_claim(x, assignment = assignment)
  }
  incomplete <- enumerate_breakpoint(coordinated_x, fit_fun = coordinated_fit)
  expect_identical(incomplete$status, "exhaustive_search_incomplete")
  expect_identical(incomplete$best_verified_cost, 2L)
  expect_true(is.na(incomplete$exact_minimum_cost))
  expect_identical(incomplete$witness$changes, "CO_K1_in=A;CO_K1_out=B")
  expect_false(candidate_row(incomplete, "CO_K1_in=A")$overturn)

  stable_raw <- load_raw_fixture("stable")
  stable_x <- validate_raw_fixture(stable_raw)
  stable_fit <- function(x, assignment = NULL, ...) {
    key <- independent_changes(stable_raw, assignment)
    if (identical(key, "ST_K1_out=B")) return(failed_fit())
    fit_claim(x, assignment = assignment)
  }
  incomplete_stable <- enumerate_breakpoint(stable_x, fit_fun = stable_fit)
  expect_identical(incomplete_stable$status, "exhaustive_search_incomplete")
  expect_true(is.na(incomplete_stable$exact_minimum_cost))

  nonfinite_fit <- function(x, assignment = NULL, ...) {
    key <- independent_changes(stable_raw, assignment)
    if (identical(key, "ST_K2_in=A")) {
      return(fit_claim(x, assignment = assignment, callback = function(...) Inf))
    }
    fit_claim(x, assignment = assignment)
  }
  nonfinite_stable <- enumerate_breakpoint(stable_x, fit_fun = nonfinite_fit)
  expect_identical(nonfinite_stable$status, "exhaustive_search_incomplete")
  expect_false(candidate_row(nonfinite_stable, "ST_K2_in=A")$overturn)
  expect_true(is.na(nonfinite_stable$exact_minimum_cost))
})

test_that("an unrelated same-cost failure does not block an exact witness", {
  raw <- load_raw_fixture("fragile")
  x <- validate_raw_fixture(raw)
  same_cost_failure <- function(x, assignment = NULL, ...) {
    key <- independent_changes(raw, assignment)
    if (identical(key, "F_C1_in=A")) return(failed_fit())
    fit_claim(x, assignment = assignment)
  }
  result <- enumerate_breakpoint(
    x, fit_fun = same_cost_failure, fit_identity = "same-cost-failure-v1"
  )
  expect_identical(result$status, "exact_minimum")
  expect_identical(result$best_verified_cost, 1L)
  expect_identical(result$exact_minimum_cost, 1L)
  expect_identical(result$witness$changes, "F_K1_in=A")
  expect_false(verify_witness(x, result)$ok)
  expect_true(verify_witness(
    x, result, fit_fun = same_cost_failure,
    fit_identity = "same-cost-failure-v1"
  )$ok)
})

test_that("truncation is unresolved and budget envelopes are monotone", {
  stable <- validate_raw_fixture(load_raw_fixture("stable"))
  truncated <- enumerate_breakpoint(stable, max_evaluations = 2L)
  expect_identical(truncated$status, "no_counterexample_found")
  expect_match(truncated$termination_reason, "limit|truncat|budget", ignore.case = TRUE)
  expect_true(is.na(truncated$exact_minimum_cost))

  minima <- numeric(3L)
  for (budget in 0:2) {
    raw <- load_raw_fixture("stable")
    raw$constraints$max_changes_global <- budget
    raw$constraints$max_changes_per_sample <- budget
    raw$constraints$max_changes_per_donor <- budget
    result <- enumerate_breakpoint(validate_raw_fixture(raw))
    evaluated <- result$assignments$feasible & result$assignments$evaluated
    minima[[budget + 1L]] <- min(result$assignments$beta[evaluated])
  }
  expect_true(all(diff(minima) <= 1e-12))
})

test_that("max_evaluations accepts positive infinity and rejects invalid limits", {
  x <- validate_raw_fixture(load_raw_fixture("fragile"))
  unbounded <- enumerate_breakpoint(x, max_evaluations = Inf)
  expect_identical(unbounded$status, "exact_minimum")

  for (limit in list(-Inf, 0, 1.5)) {
    result <- enumerate_breakpoint(x, max_evaluations = limit)
    expect_identical(result$status, "invalid_input", info = as.character(limit))
    expect_match(result$termination_reason, "max_evaluations|invalid", ignore.case = TRUE)
  }
})

test_that("enumeration capacity is a non-scientific not-assessable outcome", {
  raw <- load_raw_fixture("fragile")
  raw <- add_fixed_cell(
    raw, "CAP_1", "S_C1", "A",
    c(G_target = 0, G_ref1 = 0, G_ref2 = 0)
  )
  raw <- add_fixed_cell(
    raw, "CAP_2", "S_C2", "B",
    c(G_target = 0, G_ref1 = 0, G_ref2 = 0)
  )
  raw <- add_fixed_cell(
    raw, "CAP_3", "S_K2", "A",
    c(G_target = 0, G_ref1 = 0, G_ref2 = 0)
  )
  allowed_counts <- table(raw$allowed_labels$cell_id)
  fixed_ids <- names(allowed_counts)[allowed_counts == 1L]
  originals <- original_assignment(raw)
  raw$allowed_labels <- rbind(raw$allowed_labels, data.frame(
    cell_id = fixed_ids,
    allowed_label = ifelse(originals[fixed_ids] == "A", "B", "A"),
    evidence_source = "synthetic_fixture_design",
    evidence_version = "v1",
    stringsAsFactors = FALSE
  ))

  x <- validate_raw_fixture(raw)
  expect_true(x$ok)
  result <- enumerate_breakpoint(x)
  expect_identical(result$status, "not_assessable")
  expect_identical(result$termination_reason, "enumeration_capacity_exceeded")
  expect_true(is.na(result$best_verified_cost))
  expect_true(is.na(result$exact_minimum_cost))
})

test_that("boundary tolerance and alternate q/delta are supported or explicitly rejected", {
  raw <- load_raw_fixture("fragile")
  baseline <- independent_beta(raw)

  boundary <- raw
  boundary$claim$delta <- baseline - 0.5e-12
  boundary_validation <- validate_raw_fixture(boundary)
  if (boundary_validation$ok) {
    expect_identical(enumerate_breakpoint(boundary_validation)$status, "original_claim_not_met")
  } else {
    expect_match(paste(boundary_validation$errors, collapse = " "), "delta|unsupported", ignore.case = TRUE)
  }

  above_boundary <- raw
  above_boundary$claim$delta <- baseline - 2e-12
  above_validation <- validate_raw_fixture(above_boundary)
  if (above_validation$ok) {
    expect_true(enumerate_breakpoint(above_validation)$original_fit$claim_holds)
  } else {
    expect_match(paste(above_validation$errors, collapse = " "), "delta|unsupported", ignore.case = TRUE)
  }

  delta_zero <- raw
  delta_zero$claim$delta <- 0
  delta_zero_validation <- validate_raw_fixture(delta_zero)
  if (delta_zero_validation$ok) {
    expect_true(enumerate_breakpoint(delta_zero_validation)$original_fit$claim_holds)
  } else {
    expect_match(paste(delta_zero_validation$errors, collapse = " "), "delta|unsupported", ignore.case = TRUE)
  }

  negative_direction <- load_raw_fixture("stable")
  negative_direction$claim$q <- -1L
  negative_direction$cells$condition <- ifelse(
    negative_direction$cells$condition == "case", "control", "case"
  )
  negative_direction$specimens$condition <- ifelse(
    negative_direction$specimens$condition == "case", "control", "case"
  )
  negative_direction$claim$delta <- 0
  q_validation <- validate_raw_fixture(negative_direction)
  if (q_validation$ok) {
    expect_true(enumerate_breakpoint(q_validation)$original_fit$claim_holds)
  } else {
    expect_match(paste(q_validation$errors, collapse = " "), "q|direction|unsupported", ignore.case = TRUE)
  }
})
