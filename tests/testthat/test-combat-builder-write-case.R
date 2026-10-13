combat_builder_functions <- function() {
  builder <- system.file(
    "scripts", "build_combat_inputs.R",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  wanted <- c(".sha256", ".write_tsv", ".map_json", ".write_case")
  environment <- new.env(parent = globalenv())
  expressions <- parse(builder)
  found <- character()
  for (expression in expressions) {
    if (!is.call(expression) || !identical(expression[[1L]], as.name("<-")) ||
        !is.symbol(expression[[2L]])) next
    name <- as.character(expression[[2L]])
    if (!name %in% wanted) next
    eval(expression, envir = environment)
    found <- c(found, name)
  }
  expect_setequal(found, wanted)
  environment
}

combat_builder_toy <- function() {
  input_dir <- system.file(
    "extdata", "examples", "edger_witness",
    package = "AnnotationBreakpoint", mustWork = TRUE
  )
  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  input <- reader(
    input_dir,
    analysis_spec = yaml::read_yaml(file.path(input_dir, "analysis_spec.yaml"))
  )

  labels <- c(A = "CD4.TCM", B = "CD4.NAIVE")
  conditions <- c(control = "COVID_SEV", case = "COVID_CRIT")
  input$cells$original_label <- unname(labels[input$cells$original_label])
  input$cells$condition <- unname(conditions[input$cells$condition])
  input$cells$adt_ambiguity <- input$cells$cell_id %in% c("W_C1_in", "W_K1_in")
  input$specimens$condition <- unname(conditions[input$specimens$condition])
  input$allowed_labels$allowed_label <- unname(
    labels[input$allowed_labels$allowed_label]
  )
  input$claim$focal_label <- "CD4.TCM"
  input$claim$alternative_label <- "CD4.NAIVE"
  input$claim$case_condition <- "COVID_CRIT"
  input$claim$control_condition <- "COVID_SEV"
  names(input$constraints)[names(input$constraints) == "min_cells_A"] <-
    "min_cells_CD4.TCM"
  names(input$constraints)[names(input$constraints) == "min_cells_B"] <-
    "min_cells_CD4.NAIVE"
  names(input$constraints)[names(input$constraints) == "min_library_A"] <-
    "min_library_CD4.TCM"
  names(input$constraints)[names(input$constraints) == "min_library_B"] <-
    "min_library_CD4.NAIVE"
  sample_ids <- input$specimens$sample_id
  donor_ids <- input$specimens$donor_id
  sample_caps <- stats::setNames(rep(2, length(sample_ids)), sample_ids)
  donor_caps <- stats::setNames(rep(2, length(donor_ids)), donor_ids)
  input$constraints$max_changes_per_sample <- I(list(sample_caps))
  input$constraints$max_changes_per_donor <- I(list(donor_caps))
  input$analysis_spec <- getFromNamespace(
    ".ab_edger_spec_expected", "AnnotationBreakpoint"
  )("COVID_SEV", "COVID_CRIT")

  validation <- do.call(validate_inputs, c(
    input[c(
      "counts", "cells", "specimens", "allowed_labels", "claim",
      "constraints", "analysis_spec"
    )],
    list(runtime_policy = "strict")
  ))
  expect_true(
    validation$ok,
    info = paste(validation$errors, collapse = " | ")
  )
  if (!isTRUE(validation$ok)) {
    stop("toy COMBAT fixture did not validate: ",
         paste(validation$errors, collapse = "; "), call. = FALSE)
  }

  old_components <- c(
    "counts", "cells", "specimens", "allowed_labels", "claim",
    "constraints", "input_hash", "analysis_spec_hash"
  )
  historical_specification_hash <- getFromNamespace(
    ".ab_hash", "AnnotationBreakpoint"
  )(list(
    claim = validation$data$claim,
    analysis_spec = validation$data$analysis_spec
  ))
  expected <- data.frame(
    case_id = "TOY_COMBAT",
    component = c(old_components, "old_specification_hash"),
    sha256 = c(
      unname(validation$input_hashes[old_components]),
      historical_specification_hash
    ),
    stringsAsFactors = FALSE
  )
  config <- list(
    claim = list(
      focal_label = "CD4.TCM",
      alternative_label = "CD4.NAIVE",
      case_condition = "COVID_CRIT",
      control_condition = "COVID_SEV",
      q = as.integer(validation$data$claim$q),
      delta = validation$data$claim$delta,
      tolerance = validation$data$claim$tolerance,
      estimator_id = validation$data$claim$estimator_id
    ),
    constraints = list(
      min_cells_CD4.TCM = 1,
      min_cells_CD4.NAIVE = 1,
      min_library_CD4.TCM = 1,
      min_library_CD4.NAIVE = 1,
      max_changes_global = 2,
      max_changes_per_sample = as.list(sample_caps),
      max_changes_per_donor = as.list(donor_caps),
      min_donors_per_condition = 4
    )
  )
  list(
    input = input, validation = validation, expected = expected,
    config = config, old_components = old_components
  )
}

test_that("shipped COMBAT writer round-trips sparse inputs and hash receipts", {
  observed_r <- paste(R.version$major, R.version$minor, sep = ".")
  observed_edgeR <- as.character(utils::packageVersion("edgeR"))
  skip_if_not(identical(observed_r, "4.5.0"))
  skip_if_not(identical(observed_edgeR, "4.6.3"))

  functions <- combat_builder_functions()
  toy <- combat_builder_toy()
  output_root <- tempfile("ab-combat-writer-")
  dir.create(output_root)
  gene <- toy$validation$data$claim$focal_gene
  receipt <- functions$.write_case(
    case_id = "TOY_COMBAT",
    gene = gene,
    counts = toy$input$counts,
    cells = toy$input$cells,
    specimens = toy$input$specimens,
    allowed_labels = toy$input$allowed_labels,
    config = toy$config,
    analysis_spec = toy$input$analysis_spec,
    output_root = output_root,
    expected = toy$expected
  )

  input_dir <- file.path(output_root, "inputs", gene)
  expect_setequal(
    list.files(input_dir),
    c(
      "counts.rds", "cells.tsv", "specimens.tsv", "allowed_labels.tsv",
      "claim.tsv", "constraints.tsv"
    )
  )
  expect_s4_class(readRDS(file.path(input_dir, "counts.rds")), "dgCMatrix")

  reader <- getFromNamespace("read_input_directory", "AnnotationBreakpoint")
  rebuilt <- reader(
    input_dir,
    analysis_spec = yaml::read_yaml(
      file.path(output_root, "analysis_specs", paste0(gene, ".yaml"))
    )
  )
  rebuilt_validation <- do.call(validate_inputs, c(
    rebuilt[c(
      "counts", "cells", "specimens", "allowed_labels", "claim",
      "constraints", "analysis_spec"
    )],
    list(runtime_policy = "strict")
  ))
  expect_true(rebuilt_validation$ok)
  expect_type(rebuilt$cells$adt_ambiguity, "logical")
  expect_identical(
    names(rebuilt$constraints$max_changes_per_sample[[1L]]),
    toy$input$specimens$sample_id
  )
  expect_identical(
    names(rebuilt$constraints$max_changes_per_donor[[1L]]),
    toy$input$specimens$donor_id
  )

  historical <- receipt$component != "strict_specification_hash"
  expect_true(all(receipt$matches[historical]))
  expect_identical(
    receipt$observed_sha256[match(
      toy$old_components, receipt$component
    )],
    unname(rebuilt_validation$input_hashes[toy$old_components])
  )
  expect_identical(
    receipt$observed_sha256[receipt$component == "old_specification_hash"],
    toy$expected$sha256[toy$expected$component == "old_specification_hash"]
  )
  strict_row <- receipt$component == "strict_specification_hash"
  expect_true(is.na(receipt$expected_sha256[strict_row]))
  expect_true(is.na(receipt$matches[strict_row]))
  expect_match(receipt$observed_sha256[strict_row], "^[0-9a-f]{64}$")
  expect_identical(
    receipt$translation[strict_row],
    "0.3.0_strict_runtime_policy_included"
  )
})
