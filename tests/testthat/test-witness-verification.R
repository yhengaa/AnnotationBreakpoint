test_that("witness verification rejects cost, effect, label, hash, and immutable-data tampering", {
  raw <- load_raw_fixture("fragile")
  x <- validate_raw_fixture(raw)
  result <- enumerate_breakpoint(x)
  expect_true(verify_witness(x, result)$ok)

  wrong_cost <- result$witness
  wrong_cost$cost <- wrong_cost$cost + 1L
  expect_false(verify_witness(x, result, witness = wrong_cost)$ok)

  fractional_cost <- result$witness
  fractional_cost$cost <- 1.5
  expect_false(verify_witness(x, result, witness = fractional_cost)$ok)

  wrong_effect <- result$witness
  wrong_effect$beta <- wrong_effect$beta + 0.25
  expect_false(verify_witness(x, result, witness = wrong_effect)$ok)

  wrong_label <- result$witness
  wrong_label$assignment[["F_K1_in"]] <- "C"
  expect_false(verify_witness(x, result, witness = wrong_label)$ok)

  wrong_hashes <- result$input_hashes
  wrong_hashes[[1L]] <- paste0("tampered-", wrong_hashes[[1L]])
  expect_false(verify_witness(
    x,
    result,
    expected_input_hashes = wrong_hashes
  )$ok)

  wrong_original <- result
  wrong_original$original_fit$beta <- wrong_original$original_fit$beta + 0.25
  expect_false(verify_witness(x, wrong_original)$ok)

  unmet_original <- result
  unmet_original$status <- "verified_counterexample"
  unmet_original$original_fit$claim_holds <- FALSE
  expect_false(verify_witness(x, unmet_original)$ok)

  changed_identity <- raw
  changed_identity$specimens$donor_id[
    changed_identity$specimens$sample_id == "S_K1"
  ] <- "D_K1_CHANGED"
  changed_identity$cells$donor_id[changed_identity$cells$sample_id == "S_K1"] <- "D_K1_CHANGED"
  changed_x <- validate_raw_fixture(changed_identity)
  expect_true(changed_x$ok)
  expect_false(verify_witness(
    changed_x,
    result,
    expected_input_hashes = result$input_hashes
  )$ok)
})

test_that("reserved display delimiters are rejected as structured invalid input", {
  raw <- load_raw_fixture("coordinated")
  raw <- rename_cell(raw, "CO_K1_in", "x")
  raw <- rename_cell(raw, "CO_K1_out", "y")
  raw <- add_fixed_cell(
    raw,
    "x=A;y",
    "S_K2",
    "A",
    c(G_target = 0, G_ref1 = 0, G_ref2 = 0)
  )
  raw$allowed_labels <- rbind(raw$allowed_labels, data.frame(
    cell_id = "x=A;y",
    allowed_label = "B",
    evidence_source = "synthetic_fixture_design",
    evidence_version = "v1",
    stringsAsFactors = FALSE
  ))
  expect_invalid(raw, "delimiter|reserved|cell")
})
