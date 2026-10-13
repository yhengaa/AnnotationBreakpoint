# Annotation Breakpoint

**Project closed.** Version 0.2.0 is a software utility with Method B (`proxy_prefix`) as its default. The methods claim is dropped; no further experiments or development are planned.

Annotation Breakpoint measures how sensitive a predeclared donor-level
expression claim is to a finite set of plausible cell-type label changes. It
returns a replayable changed assignment when one overturns the claim, and it
states whether the search established an exact minimum, an upper bound, finite-
domain stability, or only an unresolved no-hit result.

A counterexample shows that an expression claim is sensitive to cell-type
labels under the supplied assumptions. It does not show that a moved cell is
mislabeled, identify the biologically correct annotation, or disprove the
underlying biological effect.

## Software disposition and evidence

Annotation Breakpoint is released as a sensitivity utility. The registered
comparison ended the methods claim: across 13 fixed cases the frozen
proxy-prefix procedure found 13 valid counterexamples and canonical greedy
found 11, while greedy had no exclusive discovery or cost advantage of two or
more cells. This does not demonstrate a unique solver contribution, biological
utility, or held-out performance. The publication snapshot includes the accepted
[result note](reports/DECISIVE-COMPARATOR-001/RESULT_NOTE.md),
[complete table](reports/DECISIVE-COMPARATOR-001/results_13x3.tsv), and
[accepted table review](reports/DECISIVE-COMPARATOR-001/accepted_table_review.json),
which preserve the registered decision, independent replay, unresolved
outcomes, and failed-attempt accounting.

The publication-safe [final evaluation evidence summary](reports/FINAL-EVALUATION-001/EVIDENCE_SUMMARY.md)
summarizes all 984 frozen scaling tasks and 96 public-data searches, with the
earlier [negative Hao result](reports/M4-PUBLIC-001/NEGATIVE_HAO.md) preserved.

Small-domain exact-reference results support the implementation. The scaling
decision remained inconclusive because it had fewer than the required 12
independent reference cases, and observed K=192 coverage depended strongly on
the fit allowance. On COMBAT, one permitted cell reassignment crossed the
near-threshold IL18R1 effect claim; added discovery and broad practical or
biological utility were not demonstrated. The earlier Hao result was negative.

The accepted historical 0.1.0 package evidence and prior correctness checks
remain unchanged. The new 0.2.0 `proxy_prefix` integration is checked and
packaged separately during software closeout. The evaluation reports
distinguish unresolved searches, finite-domain proofs, and separate search
allowances; they do not turn no-hit results into robustness certificates.

## Default practical engine

The default API and CLI method is `proxy_prefix`, labeled **observed-centered
proxy-ranked cumulative prefixes**. It deterministically ranks permitted
single-cell label changes with the frozen observed-centered proxy, then
evaluates only cumulative assignments at the registered prefix sizes
`1, 2, 4, 8, 16, 32, 64, 128`. A smaller domain evaluates only registered
prefixes that are available. If an exact prefix assignment violates a global,
sample, donor, cell-count, library, or replication constraint, that prefix is
rejected; the engine does not skip cells or fill the prefix with lower-ranked
changes.

Every available admissible registered prefix receives a complete authoritative
refit, even after an earlier fitted prefix crosses the threshold. The engine
then considers prefixes in registered order and returns the earliest crossing
prefix that passes an independent replay from immutable input. The normal
engine does not run the comparator's singleton-correlation diagnostics, which
were benchmark-only and could not create a discovery. A witness is a verified
upper bound on breakpoint cost. A completed prefix search with no replayed
witness is `no_counterexample_found` and remains unresolved.

`proxy_prefix` is separate from the adaptive historical `ranked` engine. The R
API retains explicit experimental legacy selections, including `ranked`, for
reproducibility. The CLI exposes only `proxy_prefix`, `exact`, `random`,
`greedy`, and `pair`. It never silently substitutes `ranked` for the default.
`proxy_prefix` is bidirectional and rejects `direction="outgoing_only"`; that
restricted direction applies only to compatible explicitly selected legacy
practical methods.

## Accepted saved examples

The final public examples use the COMBAT CD4.TCM expression contrast
COVID_CRIT versus COVID_SEV, with 16 independent donors (8 per group),
CD4.TCM/CD4.NAIVE alternatives, and `K=8,898` ambiguous cells. They retain the
original 5% domain: global budget 445, its frozen per-donor caps, and effect
threshold 0.5. They are saved development examples, not new evaluation claims.

| Gene | Original beta | Replayed witness beta | Cost | Meaning |
| --- | ---: | ---: | ---: | --- |
| FAM118A | 0.8106752561679119 | 0.4731965904097605 | 64 | verified upper bound |
| PDXK | 0.7274641768393035 | 0.3373883736018561 | 32 | verified upper bound |

Inspect the accepted rows without fitting a model:

```sh
awk -F '\t' 'NR == 1 || (($1 == "REAL_FAM118A" || $1 == "REAL_PDXK") && $3 == "B")' \
  reports/DECISIVE-COMPARATOR-001/results_13x3.tsv
```

The exact registered prefix effects are in
[`method_b_prefixes.tsv`](reports/DECISIVE-COMPARATOR-001/method_b_prefixes.tsv).
The accepted saved-table reconciliation, including independent replay status,
is in [`accepted_table_review.json`](reports/DECISIVE-COMPARATOR-001/accepted_table_review.json).
Detailed immutable-input replay receipts remain in the accepted local evidence
archive and are excluded here because they contain local storage paths.
Neither cost is a proven minimum. Crossing 0.5 concerns the fitted expression
effect, not statistical significance, annotation truth, or a biological
discovery.

## Registered target

For fitted log2 effect `beta`, registered direction `q` in `{-1, +1}`, effect
threshold `delta >= 0`, and nonnegative numerical `tolerance`, the original
claim is exactly

```text
q * beta > delta + tolerance
```

An evaluable changed assignment overturns it exactly when

```text
q * beta <= delta + tolerance
```

Equivalently, the software stores `margin = q * beta - delta`, requires
`margin > tolerance` for the claim, and uses `margin <= tolerance` for the
overturn boundary. This is an effect-threshold sensitivity target. It is not a
p-value threshold, a test of annotation correctness, or a cell-level
mislabeling probability.

The primary cost is the number of cells whose assigned label differs from its
original label. A two-cell swap costs two. Every moved cell remains in its
original sample, donor, and condition, and its counts are neither modified nor
dropped.

## Requirements and installation

The package is an R 4.5.0 workflow. The authoritative backend requires exactly
edgeR 4.6.3; a version mismatch is reported as `not_assessable`. The package
also imports `Matrix`, `digest`, `jsonlite`, and `yaml`.

Install from the source tree:

```sh
R CMD INSTALL .
```

Or install the versioned source tarball produced by the release workflow:

```sh
R CMD INSTALL AnnotationBreakpoint_0.2.0.tar.gz
```

The source tarball does not contain the cached Hao5P public counts.

## Quick start from the source tree

The shipped `exact_toy` example uses the transparent Layer A estimator and has
no answer table:

```sh
Rscript --vanilla exec/annotation-breakpoint \
  --input-dir inst/extdata/examples/exact_toy \
  --layer-a \
  --output-dir ab-output/exact-toy \
  --method exact
```

The shipped `edger_witness` example rebuilds the authoritative edgeR analysis
for every evaluated assignment. It uses the default `proxy_prefix` engine:

```sh
Rscript --vanilla exec/annotation-breakpoint \
  --input-dir inst/extdata/examples/edger_witness \
  --analysis-spec inst/extdata/examples/edger_witness/analysis_spec.yaml \
  --output-dir ab-output/edger-witness \
  --method proxy_prefix \
  --max-fits 32 \
  --max-proposals 256 \
  --max-seconds 60 \
  --seed 82901 \
  --direction all
```

The normal final line has this machine-readable shape:

```text
status=<status> best_verified_cost=<integer-or-NA> exact_minimum_cost=<integer-or-NA> verification_ok=TRUE method=<method> output_dir=<path>
```

`exact_toy` is expected to return `status=exact_minimum` with exact cost 1.
`edger_witness` with the command above is expected to return
`status=verified_counterexample` with best verified cost 1 and
`exact_minimum_cost=NA`: the practical result is a verified upper bound, not an
exact-minimum claim.

## Quick start from an installed package

The executable and example files are installed inside the package. Resolve
them with `system.file()` rather than assuming a working directory:

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"
AB_EXACT="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "exact_toy", package="AnnotationBreakpoint"))')"
AB_EDGER="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "edger_witness", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EXACT" \
  --layer-a \
  --output-dir ab-output/installed-exact \
  --method exact

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EDGER" \
  --analysis-spec "$AB_EDGER/analysis_spec.yaml" \
  --output-dir ab-output/installed-edger \
  --method proxy_prefix \
  --max-fits 32 \
  --max-proposals 256 \
  --max-seconds 60 \
  --seed 82901 \
  --direction all
```

When `method` is omitted, both the R API and CLI use `proxy_prefix`. The CLI
also accepts `exact`, `random`, `greedy`, and `pair`; adaptive `ranked` is an
R-API-only experimental legacy selection. Exact enumeration accepts
`--max-evaluations`; practical methods accept `--max-fits`,
`--max-proposals`, `--max-seconds`, and `--seed`. Practical search can use
`--direction outgoing_only`; exact enumeration is bidirectional only.

## Input files

Normal file input is a directory containing six tab-separated files. Select
exactly one analysis target: pass `--layer-a`, or pass the authoritative YAML
with `--analysis-spec`. Supplying neither or both is refused.

| File | Required content |
| --- | --- |
| `counts.tsv` | `gene_id` followed by cell-ID columns; finite nonnegative integer raw counts |
| `cells.tsv` | `cell_id`, `sample_id`, `donor_id`, `condition`, `original_label` |
| `specimens.tsv` | `sample_id`, `donor_id`, `condition` |
| `allowed_labels.tsv` | `cell_id`, `allowed_label`, `evidence_source`, `evidence_version`; optional `calibrated_confidence` |
| `claim.tsv` | one row declaring gene, focal and alternative labels, two condition names, `q`, `delta`, `tolerance`, feature universe, and estimator |
| `constraints.tsv` | one row declaring global/sample/donor budgets, per-sample label minima, library minima, and donor replication minimum |

Each cell needs at least one allowed-label row, and its original label must be
allowed. A cell with one allowed label is immutable. `evidence_source` and
`evidence_version` must identify how each allowed alternative was constructed;
they are copied into a verified witness. Optional `calibrated_confidence` is
validated to lie in `[0, 1]`, but version 1 search uses unit Hamming cost and
does not interpret that field as a probability.

Version 1 requires one sample per biological donor and exactly two distinct,
nonempty declared conditions, with at least two donors per condition. It uses
an intercept plus case-versus-control condition coefficient and does not add a
donor fixed effect when condition is constant within donor. Paired or repeated
specimens are unsupported; do not flatten them into independent donors.

The constraints enforce:

- maximum changed cells globally, per sample, and per donor;
- minimum cells in each registered label within each sample;
- minimum library count in each registered label within each sample; and
- enough valid focal-label donors in both conditions.

The minimum-library rule is computed over the registered feature universe
before backend `filterByExpr` filtering. Backend filtering then rebuilds its
own retained library sizes. Cell identities and feature counts remain
conserved within their original specimens; only their allocation between the
two registered labels changes.

## Analysis policies

With `--layer-a`, Layer A computes each donor/sample's focal-label
expression as

```text
log2(1e6 * (target_count + 0.5) /
     (feature_universe_library_count + 0.5 * number_of_features))
```

and uses the equally weighted case mean minus control mean as `beta`. This is a
transparent exact-oracle instrument for tiny examples.

The authoritative `edger-ql-pseudobulk-v1` policy uses R 4.5.0 and edgeR 4.6.3.
For every assignment it freshly rebuilds focal-population donor pseudobulks,
the `~ condition` design, `filterByExpr`, post-filter library sizes, TMM
normalization, dispersion and `glmQLFit`, then obtains the case-versus-control
log2 effect from `glmQLFTest`. It retains at least 50 filtered features and the
target gene. No offset, normalization factor, dispersion, or fitted model is
reused across assignments. The exact options are stored in the supplied and
rendered analysis specification.

## Search and result interpretation

Exact enumeration is limited to at most 12 ambiguous cells and 4,096 allowed
assignments. `exact_minimum` means a verified witness exists and every cheaper
feasible assignment was ruled out for the same counts, claim, constraints, and
analysis. It does not turn an evaluated-assignment effect range into a general
full-domain bound.

Practical engines return upper-bound evidence only. `proxy_prefix` is the
software default under the final registered DROP disposition. It fully fits
every available admissible registered prefix, then returns the earliest
replay-passing crossing in registered order. Exact enumeration remains
available for small-domain proof. The CLI retains explicit greedy, random, and
pair selections; the R API also retains adaptive ranked for historical
reproducibility. A practical no-hit never establishes stability.

The important statuses are:

| Status | Interpretation |
| --- | --- |
| `exact_minimum` | Verified witness plus exact minimum within the declared finite problem |
| `verified_counterexample` | Verified witness and upper bound; no minimum or stability certificate |
| `exhaustively_stable` | Every feasible assignment in the declared finite domain was evaluated successfully and none overturned |
| `no_counterexample_found` | No hit within the reported search/evaluation budget; unresolved |
| `exhaustive_search_incomplete` | A processed exact domain contained failed feasible fits, so minimum/stability proof is unavailable |
| `original_claim_not_met` | The original assignment did not satisfy the registered strict claim; no zero-cost discovery is reported |
| `not_assessable` | The requested analysis could not yield an evaluable result |
| `invalid_input` | Static schema, identity, provenance, design, or original-feasibility validation failed |

A heuristic no-hit is never evidence of stability. A per-gene witness belongs
only to that gene's registered problem. Witnesses for different genes can use
incompatible assignments and must not be combined into one purported dataset;
a joint multi-gene claim requires a separately registered common-assignment
optimization.

## Outputs and replay

Every successful CLI invocation writes `claim_results.tsv`, assignment and
method logs, `fit_accounting.tsv`, `witnesses.tsv`, before/after donor
pseudobulks, input hashes, rendered constraints and analysis specifications,
`result.json`, `verification.json`, `search_log.jsonl`, `report.md`, and
`run_manifest.json`.

The CLI automatically performs an independent witness replay from immutable
inputs before writing outputs. `verification_ok=TRUE` means the assignment,
allowed-label membership, unit cost, hashes, and fresh fitted effect agreed.
Replay the same command into another output directory and compare status,
costs, witness ID, changed rows, input hashes, and recomputed effect. Runtime
fields and files containing them can differ across runs.

See [the tutorial](inst/doc/annotation-breakpoint-tutorial.md) for R API use,
failure examples, a replay check, and detailed status interpretation. The
project release record also carries the full methods/failure-mode report.

## Publication evidence

This GitHub snapshot contains selected accepted result notes and compact tables
under [`reports/`](reports/PUBLICATION/README.md). The publication manifest
records every included file and its SHA-256 digest. Multi-gigabyte raw inputs,
scratch caches, failed partial artifacts, private execution records, and the
full local-history closeout archive are intentionally excluded. The accepted
local project retains those materials separately; this repository makes no
claim that the excluded scratch storage is permanent or backed up.

## Evidence limits

The accepted M3 benchmark used 132 held-out synthetic realizations across 11
predeclared families, six ambiguous cells, and cost at most three. The public
M4 demonstration used 600 day-0 Hao5P monocytes from six donors, artificial
balanced condition labels, five genes, and six stipulated ambiguous cells. It
found no public witness: NME2 and KLF2 were stable only in their six-cell finite
domains, while LILRA5, LILRB1, and FYB1 did not meet the original claim. These
results do not establish broad biological robustness or population-wide search
performance.

The readiness claim is limited to a local sensitivity implementation. The
accepted comparator did not establish a unique solver contribution, and the
scaling and public-utility decisions remained inconclusive. The negative Hao
example remains negative. No broad practical superiority, scaling guarantee,
real biological utility, or further experiment is claimed.

S1/drop-cell paper-specific code, scANVI, and scDEcrypter were inspected as
methodological context but were not executed. The release includes matched
outgoing-only and removal/probability baselines; it does not claim to reproduce
those methods.
