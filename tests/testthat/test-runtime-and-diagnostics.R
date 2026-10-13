edger_example_inputs <- function(runtime_policy = "strict") {
  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  example_dir <- system.file(
    "extdata", "examples", "edger_witness",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  spec <- yaml::read_yaml(file.path(example_dir, "analysis_spec.yaml"))
  raw <- reader(example_dir, analysis_spec = spec)
  validate_inputs(
    raw$counts, raw$cells, raw$specimens, raw$allowed_labels,
    raw$claim, raw$constraints, analysis_spec = spec,
    runtime_policy = runtime_policy
  )
}

test_that("runtime policy is validated, persisted, and hashed", {
  strict <- edger_example_inputs("strict")
  compatible <- edger_example_inputs("compatible")

  expect_true(strict$ok)
  expect_true(compatible$ok)
  expect_identical(strict$data$runtime_policy, "strict")
  expect_identical(compatible$data$runtime_policy, "compatible")
  expect_identical(
    strict$input_hashes[["analysis_spec_hash"]],
    compatible$input_hashes[["analysis_spec_hash"]]
  )
  expect_false(identical(
    strict$input_hashes[["specification_hash"]],
    compatible$input_hashes[["specification_hash"]]
  ))

  raw <- load_raw_fixture("fragile")
  invalid <- validate_inputs(
    raw$counts, raw$cells, raw$specimens, raw$allowed_labels,
    raw$claim, raw$constraints, runtime_policy = "unsupported"
  )
  expect_false(invalid$ok)
  expect_match(paste(invalid$errors, collapse = " | "), "runtime_policy|arg")
})

test_that("runtime helpers require the registered edgeR API and patch compatibility", {
  version_allowed <- getFromNamespace(
    ".ab_runtime_version_allowed", "AnnotationBreakpoint"
  )
  expect_true(version_allowed("4.5.0", "4.5.0", "strict"))
  expect_false(version_allowed("4.5.1", "4.5.0", "strict"))
  expect_true(version_allowed("4.5.1", "4.5.0", "compatible"))
  expect_false(version_allowed("4.6.0", "4.5.0", "compatible"))
  expect_true(version_allowed("4.6.4", "4.6.3", "compatible"))
  expect_false(version_allowed("4.4.0", "4.6.3", "compatible"))

  api <- getFromNamespace(".ab_edger_api_support", "AnnotationBreakpoint")()
  expect_true(api$ok, info = paste(c(api$missing, api$details), collapse = " | "))

  observed_r <- paste(R.version$major, R.version$minor, sep = ".")
  observed_edger <- as.character(utils::packageVersion("edgeR"))
  current_data <- list(
    runtime_policy = "strict",
    analysis_spec = list(backend = list(
      r_version = observed_r,
      version = observed_edger
    ))
  )
  contract <- getFromNamespace(
    ".ab_runtime_contract", "AnnotationBreakpoint"
  )(current_data)
  expect_true(contract$compatible)
  expect_true(contract$api_supported)
  expect_identical(contract$expected_r_version, observed_r)
  expect_identical(contract$observed_r_version, observed_r)
  expect_identical(contract$expected_backend_version, observed_edger)
  expect_identical(contract$observed_backend_version, observed_edger)
})

test_that("CLI exposes compatible mode only for analysis-spec input", {
  executable <- system.file(
    "exec", "annotation-breakpoint",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  help_stdout <- tempfile("ab-cli-help-")
  help_stderr <- tempfile("ab-cli-help-error-")
  help_status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(executable), "--help"),
    stdout = help_stdout, stderr = help_stderr
  )
  expect_identical(help_status, 0L)
  expect_match(
    paste(readLines(help_stdout, warn = FALSE), collapse = "\n"),
    "--allow-version-mismatch", fixed = TRUE
  )
  expect_match(
    paste(readLines(help_stdout, warn = FALSE), collapse = "\n"),
    "--prefix-sizes", fixed = TRUE
  )

  input_dir <- system.file(
    "extdata", "examples", "exact_toy",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  refusal_stdout <- tempfile("ab-cli-refusal-")
  refusal_stderr <- tempfile("ab-cli-refusal-error-")
  refusal_status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "--vanilla", shQuote(executable),
      "--input-dir", shQuote(input_dir),
      "--layer-a", "--allow-version-mismatch",
      "--output-dir", shQuote(tempfile("ab-cli-refusal-output-"))
    ),
    stdout = refusal_stdout, stderr = refusal_stderr
  ))
  expect_true(refusal_status != 0L)
  expect_match(
    paste(readLines(refusal_stderr, warn = FALSE), collapse = "\n"),
    "--allow-version-mismatch requires normal input with --analysis-spec",
    fixed = TRUE
  )

  prefix_stdout <- tempfile("ab-cli-prefix-refusal-")
  prefix_stderr <- tempfile("ab-cli-prefix-refusal-error-")
  prefix_status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      "--vanilla", shQuote(executable),
      "--input-dir", shQuote(input_dir), "--layer-a",
      "--method", "exact", "--prefix-sizes", "1,2,3",
      "--output-dir", shQuote(tempfile("ab-cli-prefix-refusal-output-"))
    ),
    stdout = prefix_stdout, stderr = prefix_stderr
  ))
  expect_true(prefix_status != 0L)
  expect_match(
    paste(readLines(prefix_stderr, warn = FALSE), collapse = "\n"),
    "--prefix-sizes is available only with --method proxy_prefix",
    fixed = TRUE
  )
})

test_that("strict edgeR default proxy prefix replays the cost-one witness", {
  x <- edger_example_inputs("strict")
  expect_true(x$ok)
  input_hashes <- x$input_hashes
  runtime <- getFromNamespace(
    ".ab_runtime_contract", "AnnotationBreakpoint"
  )(x$data)
  expect_true(runtime$compatible)
  expect_true(runtime$api_supported)
  expect_identical(runtime$policy, "strict")
  expect_identical(runtime$expected_r_version, "4.5.0")
  expect_identical(runtime$observed_r_version, "4.5.0")
  expect_identical(runtime$expected_backend_version, "4.6.3")
  expect_identical(runtime$observed_backend_version, "4.6.3")

  result <- search_breakpoint(x, max_fits = 9L, seed = 82901L)
  expect_identical(result$method, "proxy_prefix")
  expect_identical(result$status, "verified_counterexample")
  expect_identical(result$best_verified_cost, 1L)
  expect_identical(result$witness$cost, 1L)
  expect_identical(result$original_fit$status, "evaluable")
  expect_true(result$original_fit$authoritative_rebuild)
  expect_identical(result$original_fit$backend, "edgeR")
  expect_identical(result$original_fit$backend_version, "4.6.3")
  expect_identical(result$evaluation_identity$type, "authoritative")

  expect_true(result$verification_embedded)
  expect_s3_class(result$independent_verification, "ab_verification")
  expect_true(result$independent_verification$ok)
  expect_identical(result$verification_fit_count, 2L)
  expect_identical(
    result$independent_verification$result_identity,
    getFromNamespace(".ab_result_identity", "AnnotationBreakpoint")(result)
  )

  expect_identical(result$input_hashes, input_hashes)
  expect_identical(result$independent_verification$input_hashes, input_hashes)
  expect_identical(
    getFromNamespace(".ab_recompute_input_hashes", "AnnotationBreakpoint")(x),
    input_hashes
  )
})

test_that("rendered artifacts preserve original failure diagnostics and runtime identity", {
  x <- edger_example_inputs("compatible")
  expect_true(x$ok)

  controlled_failure <- function(x, assignment = NULL, ...) {
    runtime <- getFromNamespace(
      ".ab_runtime_contract", "AnnotationBreakpoint"
    )(x$data)
    structure(list(
      status = "not_assessable",
      beta = NA_real_,
      margin = NA_real_,
      claim_holds = NA,
      overturn = FALSE,
      sample_expression = data.frame(),
      per_sample = data.frame(),
      pseudobulk_summary = data.frame(),
      error = "controlled original fit error",
      failure_reason = "controlled_original_failure",
      warnings = c("controlled warning one", "controlled warning two"),
      runtime = runtime,
      assignment = assignment,
      assignment_id = NA_character_,
      pseudobulk = NULL
    ), class = "ab_fit")
  }

  result <- enumerate_breakpoint(
    x,
    fit_fun = controlled_failure,
    fit_identity = "controlled-runtime-diagnostics-v1"
  )
  expect_identical(result$status, "not_assessable")
  verification <- verify_witness(
    x, result,
    fit_fun = controlled_failure,
    fit_identity = "controlled-runtime-diagnostics-v1"
  )
  expect_true(verification$ok)

  output_dir <- tempfile("ab-runtime-diagnostics-")
  artifacts <- render_report(result, verification, output_dir)
  result_json <- jsonlite::read_json(
    artifacts$paths[["result_json"]], simplifyVector = TRUE
  )
  claim_results <- utils::read.delim(
    artifacts$paths[["claim_results"]], sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  report <- paste(readLines(artifacts$paths[["report"]], warn = FALSE),
                  collapse = "\n")

  expect_identical(
    result_json$original_fit$failure_reason,
    "controlled_original_failure"
  )
  expect_identical(
    result_json$original_fit$error,
    "controlled original fit error"
  )
  expect_identical(
    result_json$original_fit$warnings,
    c("controlled warning one", "controlled warning two")
  )
  expect_identical(result_json$original_fit$runtime$policy, "compatible")
  expect_identical(
    result_json$original_fit$runtime$expected_r_version,
    x$data$analysis_spec$backend$r_version
  )
  expect_identical(
    result_json$original_fit$runtime$observed_r_version,
    paste(R.version$major, R.version$minor, sep = ".")
  )
  expect_identical(
    result_json$original_fit$runtime$expected_backend_version,
    x$data$analysis_spec$backend$version
  )
  expect_identical(
    result_json$original_fit$runtime$observed_backend_version,
    as.character(utils::packageVersion("edgeR"))
  )

  expect_identical(
    claim_results$original_failure_reason,
    "controlled_original_failure"
  )
  expect_identical(claim_results$original_error, "controlled original fit error")
  expect_match(claim_results$original_warnings, "controlled warning one")
  expect_identical(claim_results$runtime_policy, "compatible")
  expect_identical(
    claim_results$expected_r_version,
    x$data$analysis_spec$backend$r_version
  )
  expect_identical(
    claim_results$observed_r_version,
    paste(R.version$major, R.version$minor, sep = ".")
  )
  expect_identical(
    claim_results$expected_backend_version,
    x$data$analysis_spec$backend$version
  )
  expect_identical(
    claim_results$observed_backend_version,
    as.character(utils::packageVersion("edgeR"))
  )

  expect_match(report, "controlled_original_failure", fixed = TRUE)
  expect_match(report, "controlled original fit error", fixed = TRUE)
  expect_match(report, "controlled warning one", fixed = TRUE)
  expect_match(report, "Runtime policy: `compatible`", fixed = TRUE)
  expect_match(report, x$data$analysis_spec$backend$r_version, fixed = TRUE)
  expect_match(report, x$data$analysis_spec$backend$version, fixed = TRUE)
})
