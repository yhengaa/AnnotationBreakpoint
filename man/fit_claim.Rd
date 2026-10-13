\name{fit_claim}
\alias{fit_claim}
\title{Evaluate the registered Layer-A claim}
\description{Evaluate the transparent Layer-A estimator for one complete assignment.}
\usage{
fit_claim(x, assignment = NULL, pseudobulk = NULL, callback = NULL)
}
\arguments{
  \item{x}{A successful Layer-A \code{ab_validation} object.}
  \item{assignment}{A named complete assignment; \code{NULL} uses original
  labels.}
  \item{pseudobulk}{An optional matching \code{ab_pseudobulk}.}
  \item{callback}{An optional test callback returning one beta or a fit-like
  list.}
}
\value{
An \code{ab_fit} containing status, beta, registered margin, claim and
overturn indicators, sample-level expression, pseudobulk summary, warnings,
failure fields, assignment identity, and elapsed time.
}
\details{
Layer A computes the portable log-proportion contrast used by the installed
small fixture. It is separate from the authoritative edgeR backend. The claim
holds only when \eqn{q\beta > \delta + tolerance}.
}
\seealso{\code{\link{fit_edger_claim}}, \code{\link{validate_inputs}}}
