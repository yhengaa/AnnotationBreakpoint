closeout_validation <- function() {
  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  example_dir <- system.file(
    "extdata", "examples", "exact_toy", package = "AnnotationBreakpoint",
    mustWork = TRUE
  )
  raw <- reader(example_dir)
  do.call(
    validate_inputs,
    raw[c("counts", "cells", "specimens", "allowed_labels", "claim", "constraints")]
  )
}

test_that("proxy_prefix reproduces the frozen fixed-prefix nomination order", {
  validation <- closeout_validation()
  expect_true(validation$ok)
  plan <- getFromNamespace(".ab_proxy_prefix_plan", "AnnotationBreakpoint")(
    validation$data
  )
  expect_identical(plan$route,
                   "observed_centered_proxy_ranked_cumulative_prefixes")
  expect_identical(plan$ranking$cell_id, c("F_K1_in", "F_C1_in"))
  expect_identical(plan$ranking$to_label, c("A", "A"))
  expect_identical(plan$registered_prefix_sizes,
                   c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L))
  expect_identical(plan$available_prefix_sizes, c(1L, 2L))
  expect_identical(
    names(which(plan$assignments$k1 != validation$data$original_assignment)),
    "F_K1_in"
  )
  expect_identical(
    names(which(plan$assignments$k2 != validation$data$original_assignment)),
    c("F_C1_in", "F_K1_in")
  )
})

test_that("default proxy_prefix returns only an independently replayed upper bound", {
  validation <- closeout_validation()
  fit_by_cost <- function(x, assignment) {
    2 - sum(assignment != x$data$original_assignment)
  }
  result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L, fit_fun = fit_by_cost
  )
  expect_identical(result$method, "proxy_prefix")
  expect_identical(result$status, "verified_counterexample")
  expect_identical(result$best_verified_cost, 1L)
  expect_true(result$independent_verification$ok)
  expect_true(result$verification_embedded)
  expect_identical(result$verification_fit_count, 2L)
  expect_identical(result$fit_accounting$verification_fits, 2L)
  expect_lte(result$search_elapsed_seconds, result$all_in_elapsed_seconds)
  expect_false(result$search_complete)
  expect_true(is.na(result$exact_minimum_cost))
  expect_true(is.na(result$certified_lower_bound))

  legacy <- search_breakpoint(
    validation, method = "greedy", max_fits = 4L, seed = 1L,
    max_proposals = 4L, fit_fun = fit_by_cost
  )
  expect_identical(legacy$method, "greedy")
})

test_that("rendered replay counts reuse embedded totals and preserve legacy counts", {
  validation <- closeout_validation()
  fit_by_cost <- function(x, assignment) {
    2 - sum(assignment != x$data$original_assignment)
  }
  proxy <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L, fit_fun = fit_by_cost
  )
  proxy_dir <- tempfile("proxy-prefix-report-")
  render_report(proxy, proxy$independent_verification, proxy_dir)
  proxy_accounting <- utils::read.delim(
    file.path(proxy_dir, "fit_accounting.tsv"), check.names = FALSE
  )
  expect_identical(proxy_accounting$verification_fits, 2L)
  proxy_log <- utils::read.delim(
    file.path(proxy_dir, "method_log.tsv"), check.names = FALSE
  )
  expect_identical(proxy_log$prefix_size, c(1L, 2L))
  proxy_json <- jsonlite::read_json(file.path(proxy_dir, "result.json"))
  expect_identical(
    proxy_json$proxy_prefix$route,
    "observed_centered_proxy_ranked_cumulative_prefixes"
  )
  expect_true(proxy_json$verification_embedded)
  expect_equal(proxy_json$verification_fit_count, 2)
  expect_length(proxy_json$replay_attempts, 1L)
  expect_lte(proxy_json$search_elapsed_seconds,
             proxy_json$all_in_elapsed_seconds)

  legacy <- search_breakpoint(
    validation, method = "greedy", max_fits = 4L, seed = 1L,
    max_proposals = 4L, fit_fun = fit_by_cost
  )
  legacy_verification <- verify_witness(
    validation, legacy, fit_fun = fit_by_cost,
    fit_identity = "anonymous-custom-fit"
  )
  expect_true(legacy_verification$ok)
  legacy_dir <- tempfile("legacy-report-")
  render_report(legacy, legacy_verification, legacy_dir)
  legacy_accounting <- utils::read.delim(
    file.path(legacy_dir, "fit_accounting.tsv"), check.names = FALSE
  )
  expect_identical(
    legacy_accounting$verification_fits,
    legacy_verification$verification_fit_count
  )
})

test_that("proxy_prefix no-hit remains unresolved and replays the original", {
  validation <- closeout_validation()
  no_hit <- function(x, assignment) 2
  result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L, fit_fun = no_hit
  )
  expect_identical(result$status, "no_counterexample_found")
  expect_null(result$witness)
  expect_false(result$search_complete)
  expect_identical(result$termination_reason,
                   "registered_prefixes_exhausted_no_verified_counterexample")
  expect_true(result$independent_verification$ok)
  expect_identical(result$verification_fit_count, 1L)
  expect_identical(result$fit_accounting$verification_fits, 1L)
})

test_that("an empty ambiguity domain stays empty and unresolved", {
  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  example_dir <- system.file(
    "extdata", "examples", "exact_toy", package = "AnnotationBreakpoint",
    mustWork = TRUE
  )
  raw <- reader(example_dir)
  original <- stats::setNames(raw$cells$original_label, raw$cells$cell_id)
  raw$allowed_labels <- raw$allowed_labels[
    raw$allowed_labels$allowed_label == original[raw$allowed_labels$cell_id],
    , drop = FALSE
  ]
  validation <- do.call(
    validate_inputs,
    raw[c("counts", "cells", "specimens", "allowed_labels", "claim", "constraints")]
  )
  expect_true(validation$ok)
  result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L,
    fit_fun = function(x, assignment) 2
  )
  expect_identical(result$status, "no_counterexample_found")
  expect_identical(nrow(result$candidates), 0L)
  expect_identical(result$proposal_count, 0L)
  expect_false(result$search_complete)
})

test_that("an inadmissible exact prefix is rejected without replacement", {
  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  example_dir <- system.file(
    "extdata", "examples", "exact_toy", package = "AnnotationBreakpoint",
    mustWork = TRUE
  )
  raw <- reader(example_dir)
  copied <- c(F_C2_B_extra = "F_C2_B", F_K2_B_extra = "F_K2_B")
  additions <- raw$counts[, unname(copied), drop = FALSE]
  colnames(additions) <- names(copied)
  raw$counts <- cbind(raw$counts, additions)
  source_rows <- match(unname(copied), raw$cells$cell_id)
  cells <- raw$cells[source_rows, , drop = FALSE]
  cells$cell_id <- names(copied)
  raw$cells <- rbind(raw$cells, cells)
  raw$allowed_labels <- rbind(
    raw$allowed_labels,
    data.frame(
      cell_id = names(copied), allowed_label = "B",
      evidence_source = "software_closeout_fixture",
      evidence_version = "v1", stringsAsFactors = FALSE
    )
  )
  raw$constraints$min_cells_B <- 2L
  validation <- do.call(
    validate_inputs,
    raw[c("counts", "cells", "specimens", "allowed_labels", "claim", "constraints")]
  )
  expect_true(validation$ok)
  result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L,
    fit_fun = function(x, assignment) 2
  )
  expect_identical(result$status, "no_counterexample_found")
  expect_identical(result$candidates$prefix_size, c(1L, 2L))
  expect_false(any(result$candidates$feasible))
  expect_false(any(result$candidates$evaluated))
  expect_match(result$candidates$changes[[1L]], "F_K1_in=A", fixed = TRUE)
  expect_match(result$candidates$changes[[2L]], "F_C1_in=A", fixed = TRUE)
  expect_match(result$candidates$changes[[2L]], "F_K1_in=A", fixed = TRUE)
  expect_identical(result$attempted_fit_count, 1L)
})

test_that("proxy_prefix refuses an unverified fitted counterexample", {
  validation <- closeout_validation()
  calls <- 0L
  unstable <- function(x, assignment) {
    calls <<- calls + 1L
    if (calls <= 3L) {
      2 - sum(assignment != x$data$original_assignment)
    } else {
      99
    }
  }
  expect_error(
    search_breakpoint(
      validation, max_fits = 9L, seed = 1L, fit_fun = unstable
    ),
    "failed independent replay"
  )
})

test_that("proxy_prefix rejects route options that it does not implement", {
  validation <- closeout_validation()
  result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L,
    candidate_block_size = 4L,
    fit_fun = function(x, assignment) 2
  )
  expect_identical(result$status, "invalid_input")
  expect_match(paste(result$errors, collapse = " | "),
               "fixed registered ranking")
  resume_result <- search_breakpoint(
    validation, max_fits = 9L, seed = 1L,
    resume_evaluations = list(),
    fit_fun = function(x, assignment) 2
  )
  expect_identical(resume_result$status, "invalid_input")
})
