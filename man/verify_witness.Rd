\name{verify_witness}
\alias{verify_witness}
\title{Replay and verify an Annotation Breakpoint result}
\description{Independently replay a reported assignment and verify its result semantics.}
\usage{
verify_witness(x, result, witness = result$witness,
  expected_input_hashes = result$input_hashes, fit_fun = NULL,
  fit_identity = NULL)
}
\arguments{
  \item{x}{The successful \code{ab_validation} used for the result.}
  \item{result}{An exact \code{ab_result} or practical
  \code{ab_search_result}.}
  \item{witness}{Witness to replay; defaults to the result witness.}
  \item{expected_input_hashes}{Expected immutable input and specification
  hashes.}
  \item{fit_fun}{Matching explicit replay function for a custom-fit result.}
  \item{fit_identity}{Matching versioned identity for that custom function.}
}
\value{
An \code{ab_verification} containing \code{ok}, replay errors, recomputed
status, cost and effect, hash comparisons, changed-cell provenance,
before/after donor pseudobulks, fit count, elapsed time, and evaluation
identity.
}
\details{
Verification recomputes immutable hashes and checks result/evaluation identity,
allowed-label membership, all constraints, unit cost, and fresh original and
witness fits. Exact and stability claims receive the additional finite-domain
replay required by their certification semantics. A successful verification
is tied to the specific result identity and cannot be reused after mutation.
}
\seealso{\code{\link{render_report}}, \code{\link{enumerate_breakpoint}},
\code{\link{search_breakpoint}}}
