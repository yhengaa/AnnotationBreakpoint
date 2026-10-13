reference_normalize_assignment <- function(data, assignment = NULL) {
  cell_ids <- data$cells$cell_id
  if (is.null(assignment)) return(data$original_assignment)
  if (is.factor(assignment)) assignment <- as.character(assignment)
  if (!is.atomic(assignment) || is.null(names(assignment))) {
    stop("assignment must be a named character vector", call. = FALSE)
  }
  if (anyDuplicated(names(assignment))) {
    stop("assignment contains duplicate cell IDs", call. = FALSE)
  }
  assignment_names <- names(assignment)
  assignment <- as.character(assignment)
  names(assignment) <- assignment_names
  if (length(assignment) != length(cell_ids) ||
      !setequal(names(assignment), cell_ids)) {
    stop("assignment must map every input cell ID exactly once", call. = FALSE)
  }
  assignment <- assignment[cell_ids]
  if (anyNA(assignment) || any(!nzchar(assignment))) {
    stop("assignment labels must be nonempty", call. = FALSE)
  }
  permitted <- vapply(
    cell_ids,
    function(cell_id) assignment[[cell_id]] %in% data$allowed_by_cell[[cell_id]],
    logical(1L)
  )
  if (!all(permitted)) {
    stop(
      paste0(
        "assignment contains impermissible labels for: ",
        paste(cell_ids[!permitted], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  assignment
}

test_that("positional assignment normalization preserves registered semantics", {
  data <- list(
    cells = data.frame(cell_id = c("cell_3", "cell_1", "cell_2")),
    original_assignment = c(cell_3 = "A", cell_1 = "B", cell_2 = "A"),
    allowed_by_cell = list(
      cell_3 = c("A", "B"), cell_1 = "B", cell_2 = c("A", "B")
    )
  )
  reordered <- c(cell_2 = "B", cell_3 = "A", cell_1 = "B")
  expected <- reference_normalize_assignment(data, reordered)
  observed <- AnnotationBreakpoint:::.ab_normalize_assignment(data, reordered)
  expect_identical(observed, expected)
  expect_identical(names(observed), data$cells$cell_id)

  invalid <- reordered
  invalid[["cell_1"]] <- "A"
  expected_error <- tryCatch(
    reference_normalize_assignment(data, invalid),
    error = conditionMessage
  )
  observed_error <- tryCatch(
    AnnotationBreakpoint:::.ab_normalize_assignment(data, invalid),
    error = conditionMessage
  )
  expect_identical(observed_error, expected_error)

  missing <- reordered[names(reordered) != "cell_2"]
  expected_missing <- tryCatch(
    reference_normalize_assignment(data, missing),
    error = conditionMessage
  )
  observed_missing <- tryCatch(
    AnnotationBreakpoint:::.ab_normalize_assignment(data, missing),
    error = conditionMessage
  )
  expect_identical(observed_missing, expected_missing)
})
