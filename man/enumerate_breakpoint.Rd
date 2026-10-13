\name{enumerate_breakpoint}
\alias{enumerate_breakpoint}
\title{Enumerate a finite assignment domain exactly}
\description{Evaluate every feasible assignment needed for an exact breakpoint result.}
\usage{
enumerate_breakpoint(x, fit_fun = NULL, max_evaluations = Inf,
  fit_identity = NULL)
}
\arguments{
  \item{x}{A successful \code{ab_validation} object.}
  \item{fit_fun}{An optional complete-fit function. The authoritative fit for
  \code{x} is used by default.}
  \item{max_evaluations}{Maximum feasible fits, including the original, or
  \code{Inf}.}
  \item{fit_identity}{A nonempty versioned identity required when a custom fit
  function is used with the authoritative analysis.}
}
\value{
An \code{ab_result} containing the finite assignment table, fit accounting,
original fit, status, optional witness, immutable hashes, evaluation identity,
and independent verification.
}
\details{
The exact route accepts at most 12 ambiguous cells and 4,096 allowed
assignments. \code{exact_minimum} requires a verified crossing and successful
evaluation of every cheaper feasible assignment. \code{exhaustively_stable}
requires successful evaluation of every feasible assignment with no crossing.
Failures or a finite evaluation cap prevent those certifications.
}
\seealso{\code{\link{search_breakpoint}}, \code{\link{verify_witness}}}
