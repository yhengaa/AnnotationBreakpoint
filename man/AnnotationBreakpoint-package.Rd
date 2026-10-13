\name{AnnotationBreakpoint-package}
\alias{AnnotationBreakpoint}
\alias{AnnotationBreakpoint-package}
\title{Annotation reassignment sensitivity analysis}
\description{
AnnotationBreakpoint validates registered single-cell count and metadata
inputs, constructs donor-level pseudobulks, searches permitted label
assignments, independently replays results, and writes machine-readable
reports.
}
\details{
The registered effect claim holds only when
\eqn{q\beta > \delta + tolerance}. Practical search can establish a
replay-verified upper bound on breakpoint cost. Exact enumeration can establish
an exact minimum or stability within a declared finite domain.

The default \code{proxy_prefix} route fixes one observed-label ranking and
by default fully refits available cumulative prefixes of sizes 1, 2, 4, 8,
16, 32, 64, and 128. A validated optional prefix schedule selects different
exact prefixes from the same frozen ranking. See
\code{vignette("annotation-breakpoint")} for a complete workflow.
}
\seealso{
\code{\link{validate_inputs}}, \code{\link{search_breakpoint}},
\code{\link{enumerate_breakpoint}}, \code{\link{verify_witness}}
}
