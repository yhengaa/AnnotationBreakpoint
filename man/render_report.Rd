\name{render_report}
\alias{render_report}
\title{Render deterministic result and verification artifacts}
\description{Write machine-readable and human-readable artifacts for a completed run.}
\usage{
render_report(result, verification, output_dir, run_metadata = list())
}
\arguments{
  \item{result}{An exact or practical \code{ab_result}.}
  \item{verification}{A successful \code{ab_verification} for the unchanged
  result.}
  \item{output_dir}{Destination directory. An existing directory is replaced
  through the function's staged rendering process.}
  \item{run_metadata}{Optional named caller metadata recorded in the report.}
}
\value{
An invisible \code{ab_report_artifacts} object containing the normalized
output directory, artifact paths, and SHA-256 hashes.
}
\details{
Artifacts include the claim result, candidate or prefix log, fit accounting,
changed-cell witness, donor pseudobulks, hashes, constraints, analysis
specification, JSON records, and a Markdown report. Original-fit status,
structured failure reason, backend error, warnings, and expected/observed
runtime versions are propagated to the human- and machine-readable summaries.
Rendering stops if verification is unsuccessful or belongs to a different or
modified result.
}
\seealso{\code{\link{verify_witness}}}
