\name{validate_inputs}
\alias{validate_inputs}
\title{Validate and normalize Annotation Breakpoint inputs}
\description{Validate the six-table schema, registered claim, and analysis contract.}
\usage{
validate_inputs(counts, cells, specimens, allowed_labels, claim,
  constraints, analysis_spec = NULL,
  runtime_policy = c("strict", "compatible"))
}
\arguments{
  \item{counts}{A feature-by-cell \code{Matrix::dgCMatrix} of finite,
  nonnegative integer raw counts with unique feature and cell names.}
  \item{cells}{Cell metadata with \code{cell_id}, \code{sample_id},
  \code{donor_id}, \code{condition}, and \code{original_label}.}
  \item{specimens}{One row per sample with matching sample, donor, and
  condition fields.}
  \item{allowed_labels}{One row per cell/allowed-label pair, including
  nonempty \code{evidence_source} and \code{evidence_version}.}
  \item{claim}{A named one-row object containing the focal gene and labels,
  two-condition contrast, direction, threshold, tolerance, feature universe,
  and estimator identity.}
  \item{constraints}{A named one-row object containing global, sample, and
  donor change limits; cell and library minima; and donor replication.}
  \item{analysis_spec}{The frozen authoritative edgeR specification.
  \code{NULL} selects the transparent Layer-A estimator.}
  \item{runtime_policy}{\code{"strict"} requires the expected R and edgeR
  versions exactly. \code{"compatible"} permits patch-level differences only
  after required edgeR API capability checks.}
}
\value{
An object of class \code{ab_validation}. Successful objects have
\code{ok = TRUE}, normalized data, runtime policy and identity, and immutable
input/specification hashes. Invalid inputs return \code{ok = FALSE} with
structured validation errors.
}
\details{
Validation enforces cell and specimen identity, label provenance, count
integrity, the registered two-label problem, condition replication, budgets,
minimum cells and libraries, and the analysis policy. Compatible runtime mode
does not change any scientific analysis option. Major or minor version
differences and unsupported edgeR APIs remain invalid for authoritative fits.
}
\seealso{\code{\link{fit_claim}}, \code{\link{fit_edger_claim}}}
