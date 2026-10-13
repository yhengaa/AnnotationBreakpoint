\name{search_breakpoint}
\alias{search_breakpoint}
\title{Search for a replay-verified annotation counterexample}
\description{Search a constrained reassignment domain and replay any reported witness.}
\usage{
search_breakpoint(x, method = c("proxy_prefix", "greedy", "random", "pair",
  "greedy_random", "random_block", "random_prune", "ranked"), max_fits,
  seed, max_proposals = Inf, max_seconds = Inf,
  direction = c("all", "outgoing_only"), fit_fun = NULL,
  fit_identity = NULL, candidate_block_size = 8L,
  pruning_fit_reserve = NULL, pruning_chunk_size = NULL,
  proxy_target_pseudocount = 0.5, proxy_library_pseudocount = 1,
  proxy_scale = 1e6, warm_start_assignment = NULL,
  progress_callback = NULL, progress_every = 1L,
  resume_evaluations = NULL, prefix_sizes = NULL)
}
\arguments{
  \item{x}{A successful \code{ab_validation} object.}
  \item{method}{Search route. \code{"proxy_prefix"} is the supported default.
  Other named values are retained as explicit R-API compatibility routes.}
  \item{max_fits}{Positive full-fit cap including the original, or
  \code{Inf}.}
  \item{seed}{One integer used without changing the caller's random state.}
  \item{max_proposals}{Maximum proposals after the original assignment.}
  \item{max_seconds}{Wall-time cap checked between proposals.}
  \item{direction}{All permitted moves or focal-population outgoing moves.
  The default proxy-prefix route requires \code{"all"}.}
  \item{fit_fun}{An optional custom complete-fit function.}
  \item{fit_identity}{A versioned identity required for an authoritative
  custom fit function.}
  \item{candidate_block_size}{Full fits per continuation block for the retained
  ranked compatibility route.}
  \item{pruning_fit_reserve}{Fits reserved for deletion pruning by the retained
  random-prune compatibility route.}
  \item{pruning_chunk_size}{Initial deletion chunk size for that compatibility
  route.}
  \item{proxy_target_pseudocount}{Target-count pseudocount for fixed proxy
  ranking.}
  \item{proxy_library_pseudocount}{Library-count pseudocount for fixed proxy
  ranking.}
  \item{proxy_scale}{Library scale for fixed proxy ranking.}
  \item{warm_start_assignment}{Optional complete assignment for a pruning-only
  compatibility diagnostic.}
  \item{progress_callback}{Optional function receiving a compact record after
  each newly executed full fit.}
  \item{progress_every}{Callback interval. Resumable ledgers require 1.}
  \item{resume_evaluations}{Optional named compact prior fit records keyed by
  canonical assignment ID.}
  \item{prefix_sizes}{Optional unique, strictly increasing vector of positive
  integer prefix sizes for \code{method = "proxy_prefix"}. \code{NULL} uses
  \code{c(1, 2, 4, 8, 16, 32, 64, 128)}.}
}
\value{
An \code{ab_search_result} with original-fit diagnostics, proposal and fit
accounting, candidate log, optional verified witness, runtime and evaluation
identity, immutable hashes, and independent verification.
}
\details{
The default route ranks permitted single-cell changes once at the observed
assignment and fully refits available cumulative prefixes of sizes 1, 2, 4,
8, 16, 32, 64, and 128. It does not rerank or repair an inadmissible prefix.
The smallest evaluated crossing is replayed independently. A returned witness
is an upper bound; a no-hit is unresolved. A caller-supplied
\code{prefix_sizes} vector changes only the exact prefixes evaluated from the
same single ranking. The schedule is part of result identity and the prefix
log, and constraint rejection remains exact-prefix rejection without repair.
The result records the requested schedule in
\code{proxy_prefix$registered_prefix_sizes} and its domain-available subset in
\code{proxy_prefix$available_prefix_sizes}.
}
\seealso{\code{\link{enumerate_breakpoint}}, \code{\link{verify_witness}}}
