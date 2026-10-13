\name{build_pseudobulk}
\alias{build_pseudobulk}
\title{Construct sample-by-label pseudobulks}
\description{Aggregate immutable cell counts by specimen and candidate label.}
\usage{
build_pseudobulk(x, assignment = NULL,
  method = c("full", "incremental"), base = NULL, moves = NULL)
}
\arguments{
  \item{x}{A successful \code{ab_validation} object.}
  \item{assignment}{A named complete cell-to-label assignment. \code{NULL}
  uses the original labels.}
  \item{method}{\code{"full"} reconstructs all aggregates;
  \code{"incremental"} updates \code{base} and verifies the result against a
  full reconstruction.}
  \item{base}{An \code{ab_pseudobulk} from the same immutable inputs, required
  for incremental updates.}
  \item{moves}{A named cell-to-new-label vector or move table used in
  incremental mode. Supply either \code{moves} or \code{assignment}.}
}
\value{
An \code{ab_pseudobulk} containing sparse counts, column metadata, assignment,
move records, per-sample/per-feature conservation checks, immutable hashes, and
the canonical assignment ID.
}
\details{
Counts remain attached to their original sample and donor. A label move only
changes the sample-by-label aggregate that receives that cell's count vector.
The function stops if conservation fails or an incremental update differs from
independent full reconstruction.
}
\seealso{\code{\link{fit_claim}}, \code{\link{fit_edger_claim}}}
