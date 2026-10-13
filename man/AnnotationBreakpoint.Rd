\name{AnnotationBreakpoint}
\alias{AnnotationBreakpoint}
\alias{validate_inputs}
\alias{build_pseudobulk}
\alias{fit_claim}
\alias{fit_edger_claim}
\alias{enumerate_breakpoint}
\alias{search_breakpoint}
\alias{verify_witness}
\alias{render_report}
\title{Constrained annotation-reassignment sensitivity analysis}
\description{
Validates immutable single-cell count and metadata inputs, builds donor-level
pseudobulks, evaluates a registered expression-effect claim, searches the
allowed label assignments, independently replays a result, and writes a
machine-readable report.
}
\usage{
validate_inputs(counts, cells, specimens, allowed_labels, claim,
  constraints, analysis_spec = NULL)

build_pseudobulk(x, assignment = NULL, method = c("full", "incremental"),
  base = NULL, moves = NULL)

fit_claim(x, assignment = NULL, pseudobulk = NULL, callback = NULL)

fit_edger_claim(x, assignment = NULL, pseudobulk = NULL)

enumerate_breakpoint(x, fit_fun = NULL, max_evaluations = Inf,
  fit_identity = NULL)

search_breakpoint(x, method = c("proxy_prefix", "greedy", "random", "pair",
  "greedy_random", "random_block", "random_prune", "ranked"), max_fits,
  seed, max_proposals = Inf, max_seconds = Inf,
  direction = c("all", "outgoing_only"), fit_fun = NULL,
  fit_identity = NULL, candidate_block_size = 8L,
  pruning_fit_reserve = NULL, pruning_chunk_size = NULL,
  proxy_target_pseudocount = 0.5, proxy_library_pseudocount = 1,
  proxy_scale = 1e6, warm_start_assignment = NULL,
  progress_callback = NULL, progress_every = 1L,
  resume_evaluations = NULL)

verify_witness(x, result, witness = result$witness,
  expected_input_hashes = result$input_hashes, fit_fun = NULL,
  fit_identity = NULL)

render_report(result, verification, output_dir, run_metadata = list())
}
\arguments{
  \item{counts}{Sparse feature-by-cell integer count matrix.}
  \item{cells}{Cell metadata with immutable specimen, donor, condition, and original-label fields.}
  \item{specimens}{One row per independent specimen.}
  \item{allowed_labels}{Cell-specific allowed labels with evidence provenance.}
  \item{claim}{Registered focal gene, population, contrast, direction, threshold, and tolerance.}
  \item{constraints}{Registered move budget and minimum cell/library constraints.}
  \item{analysis_spec}{Optional authoritative backend specification; \code{NULL} selects Layer A.}
  \item{x}{A successful validation object.}
  \item{assignment}{A named complete label assignment; \code{NULL} uses original labels.}
  \item{method}{Aggregation mode or search method, as applicable. Practical
  search defaults to observed-centered proxy-ranked cumulative prefixes.}
  \item{base}{Existing pseudobulk for an incremental update.}
  \item{moves}{Explicit cell label moves for an incremental update.}
  \item{pseudobulk}{Optional matching aggregate; authoritative fitting rebuilds label-dependent steps.}
  \item{callback}{Optional Layer-A test callback.}
  \item{fit_fun}{Optional complete-fit callback.}
  \item{fit_identity}{Required versioned identity for an authoritative custom callback.}
  \item{max_evaluations}{Maximum exact feasible fits.}
  \item{max_fits}{Maximum practical-search fits, including the original, or \code{Inf}.}
  \item{seed}{Deterministic search seed.}
  \item{max_proposals}{Maximum proposals after the original.}
  \item{max_seconds}{Wall-time cap checked between proposals.}
  \item{direction}{All allowed moves or focal-population outgoing moves.}
  \item{candidate_block_size}{Number of fully fitted candidates per continuation block.}
  \item{pruning_fit_reserve}{For finite fit caps, the number withheld from random discovery for deletion pruning.}
  \item{pruning_chunk_size}{Initial deterministic deletion chunk size.}
  \item{proxy_target_pseudocount}{Target-count pseudocount for proposal ranking.}
  \item{proxy_library_pseudocount}{Library-count pseudocount for proposal ranking.}
  \item{proxy_scale}{Multiplicative library scale for proposal ranking.}
  \item{warm_start_assignment}{Optional complete assignment for a pruning-only diagnostic.}
  \item{progress_callback}{Optional callback receiving compact fit progress for deterministic resumption.}
  \item{progress_every}{Number of newly executed full fits between progress callbacks.}
  \item{resume_evaluations}{Optional named compact prior fit records keyed by assignment ID.}
  \item{result}{An exact or practical result.}
  \item{witness}{Witness to replay; defaults to the result witness.}
  \item{expected_input_hashes}{Expected immutable input hashes.}
  \item{verification}{Verification for the same result.}
  \item{output_dir}{New or replaceable report directory.}
  \item{run_metadata}{Named caller metadata recorded in the report.}
}
\value{
Typed validation, pseudobulk, fit, result, verification, or report-artifact
objects. See the installed tutorial for status and evidence semantics.
}
\details{
The registered claim holds only when \eqn{q\beta > \delta + tolerance}; the
boundary and values below it overturn the claim. Practical search returns at
most a verified upper bound. Only exhaustive evidence can establish stability,
and an exact minimum requires excluding every cheaper feasible assignment.
The default \code{proxy_prefix} route fixes its observed-label proxy ranking
before fitting, evaluates only available registered prefix sizes 1, 2, 4, 8,
16, 32, 64, and 128, and rejects an inadmissible exact prefix without replacing
its cells. Every returned witness includes an independent replay. A no-hit is
unresolved and is not a stability certificate.
}
\examples{
example_dir <- system.file("extdata", "examples", "exact_toy",
                           package = "AnnotationBreakpoint")
if (nzchar(example_dir)) {
  # The installed tutorial shows how to read these ordinary TSV inputs.
  list.files(example_dir)
}
}
