\name{fit_edger_claim}
\alias{fit_edger_claim}
\title{Fit the registered edgeR quasi-likelihood claim}
\description{Run the registered authoritative edgeR analysis for one assignment.}
\usage{
fit_edger_claim(x, assignment = NULL, pseudobulk = NULL)
}
\arguments{
  \item{x}{A successful authoritative \code{ab_validation} object.}
  \item{assignment}{A named complete assignment; \code{NULL} uses original
  labels.}
  \item{pseudobulk}{An optional matching aggregate. The authoritative fit
  rebuilds all label-dependent analysis steps.}
}
\value{
An \code{ab_fit}. Evaluable fits contain the edgeR log2 effect and diagnostics.
Unavailable fits retain a structured \code{failure_reason}, backend
\code{error}, warnings, and expected/observed R and edgeR versions.
}
\details{
The validated specification fixes filtering, TMM normalization, design,
quasi-likelihood fit, test, endpoint, and failure policy. Strict mode requires
R 4.5.0 and edgeR 4.6.3 exactly. Compatible mode permits patch-level runtime
differences only after checking the required edgeR API, including the
non-legacy quasi-likelihood route. Runtime policy and observed identity are
part of result, cache, and replay identity.
}
\seealso{\code{\link{validate_inputs}}, \code{\link{fit_claim}}}
