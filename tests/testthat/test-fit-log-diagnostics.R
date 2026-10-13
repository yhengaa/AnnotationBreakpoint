test_that("enumeration and search logs preserve fit warnings and failure reasons", {
  raw <- load_raw_fixture("fragile")
  x <- validate_raw_fixture(raw)
  controlled_fit <- function(x, assignment = NULL, ...) {
    changes <- independent_changes(raw, assignment)
    if (identical(changes, "F_C1_in=A")) {
      return(list(
        status = "not_assessable",
        error = "controlled candidate fit failed",
        failure_reason = "controlled_candidate_failure",
        warnings = "controlled failure warning"
      ))
    }
    list(
      status = "evaluable",
      beta = 2,
      warnings = paste0("controlled success warning: ", changes)
    )
  }

  exhaustive <- enumerate_breakpoint(
    x, fit_fun = controlled_fit, fit_identity = "controlled-log-fit-v1"
  )
  original <- exhaustive$candidates[
    exhaustive$candidates$changes == "original", , drop = FALSE
  ]
  successful <- exhaustive$candidates[
    exhaustive$candidates$changes == "F_K1_in=A", , drop = FALSE
  ]
  failed <- exhaustive$candidates[
    exhaustive$candidates$changes == "F_C1_in=A", , drop = FALSE
  ]
  expect_match(original$warnings, "controlled success warning: original", fixed = TRUE)
  expect_match(successful$warnings, "controlled success warning: F_K1_in=A",
               fixed = TRUE)
  expect_identical(failed$failure_reason, "controlled_candidate_failure")
  expect_match(failed$warnings, "controlled failure warning", fixed = TRUE)
  expect_match(failed$evaluation_error, "controlled candidate fit failed", fixed = TRUE)

  practical <- search_breakpoint(
    x, method = "greedy", max_fits = 4L, seed = 97L,
    max_proposals = 10L, fit_fun = controlled_fit,
    fit_identity = "controlled-log-fit-v1"
  )
  original_log <- practical$candidate_log[
    practical$candidate_log$changes == "original", , drop = FALSE
  ]
  successful_log <- practical$candidate_log[
    practical$candidate_log$changes == "F_K1_in=A", , drop = FALSE
  ]
  failed_log <- practical$candidate_log[
    practical$candidate_log$changes == "F_C1_in=A", , drop = FALSE
  ]
  expect_match(original_log$warnings, "controlled success warning: original",
               fixed = TRUE)
  expect_match(successful_log$warnings, "controlled success warning: F_K1_in=A",
               fixed = TRUE)
  expect_identical(failed_log$failure_reason, "controlled_candidate_failure")
  expect_match(failed_log$warnings, "controlled failure warning", fixed = TRUE)
  expect_match(failed_log$evaluation_error, "controlled candidate fit failed",
               fixed = TRUE)
})
