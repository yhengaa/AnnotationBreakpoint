test_that("renderer rejects failed or cross-result verification", {
  fragile_x <- validate_raw_fixture(load_raw_fixture("fragile"))
  fragile <- enumerate_breakpoint(fragile_x)
  fragile_verification <- verify_witness(fragile_x, fragile)
  expect_true(fragile_verification$ok)

  coordinated_x <- validate_raw_fixture(load_raw_fixture("coordinated"))
  coordinated <- enumerate_breakpoint(coordinated_x)
  expect_error(
    render_report(coordinated, fragile_verification, tempfile("ab-cross-case-")),
    "verification|hash|match"
  )

  tampered <- fragile$witness
  tampered$cost <- 1.5
  failed_verification <- verify_witness(fragile_x, fragile, witness = tampered)
  expect_false(failed_verification$ok)
  expect_error(
    render_report(fragile, failed_verification, tempfile("ab-failed-verification-")),
    "verification|failed"
  )
})

test_that("witness-free unresolved and unmet statuses render no fictitious witness", {
  stable_raw <- load_raw_fixture("stable")
  stable_x <- validate_raw_fixture(stable_raw)
  no_hit <- enumerate_breakpoint(stable_x, max_evaluations = 2L)

  unmet_raw <- load_raw_fixture("fragile")
  unmet_raw$claim$delta <- 10
  unmet_x <- validate_raw_fixture(unmet_raw)
  unmet <- enumerate_breakpoint(unmet_x)

  incomplete_fit <- function(x, assignment = NULL, ...) {
    key <- independent_changes(stable_raw, assignment)
    if (identical(key, "ST_K1_out=B")) return(failed_fit())
    fit_claim(x, assignment = assignment)
  }
  incomplete <- enumerate_breakpoint(
    stable_x, fit_fun = incomplete_fit,
    fit_identity = "render-incomplete-fit-v1"
  )

  cases <- list(
    list(x = stable_x, result = no_hit, fit_fun = NULL, fit_identity = NULL),
    list(x = unmet_x, result = unmet, fit_fun = NULL, fit_identity = NULL),
    list(x = stable_x, result = incomplete, fit_fun = incomplete_fit,
         fit_identity = "render-incomplete-fit-v1")
  )
  expected_statuses <- c(
    "no_counterexample_found", "original_claim_not_met", "exhaustive_search_incomplete"
  )
  for (index in seq_along(cases)) {
    item <- cases[[index]]
    expect_identical(item$result$status, expected_statuses[[index]])
    expect_null(item$result$witness)
    verification <- verify_witness(
      item$x, item$result, fit_fun = item$fit_fun,
      fit_identity = item$fit_identity
    )
    expect_true(verification$ok)
    expect_identical(verification$recomputed_status, expected_statuses[[index]])
    output_dir <- tempfile(paste0("ab-witnessless-", index, "-"))
    artifacts <- render_report(item$result, verification, output_dir)
    witness_table <- utils::read.delim(
      artifacts$paths[["witnesses"]], sep = "\t", stringsAsFactors = FALSE
    )
    claim_table <- utils::read.delim(
      artifacts$paths[["claim_results"]], sep = "\t", stringsAsFactors = FALSE
    )
    expect_identical(nrow(witness_table), 0L)
    expect_identical(claim_table$status[[1L]], expected_statuses[[index]])
    expect_true(is.na(claim_table$witness_id[[1L]]))
  }
})

