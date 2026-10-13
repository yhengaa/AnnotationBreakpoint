# Annotation Breakpoint tutorial

This tutorial runs the two small examples shipped with the installed package.
`exact_toy` exercises transparent Layer A exact enumeration. `edger_witness`
uses the authoritative edgeR 4.6.3 backend and the default `proxy_prefix`
search. Neither
example needs files from the developer's benchmark directories.

## 1. Resolve installed example paths

In R:

```r
library(AnnotationBreakpoint)

exact_dir <- system.file(
  "extdata", "examples", "exact_toy",
  package = "AnnotationBreakpoint", mustWork = TRUE
)
edger_dir <- system.file(
  "extdata", "examples", "edger_witness",
  package = "AnnotationBreakpoint", mustWork = TRUE
)
cli <- system.file(
  "exec", "annotation-breakpoint",
  package = "AnnotationBreakpoint", mustWork = TRUE
)
```

For a source-tree run, replace those values with
`inst/extdata/examples/exact_toy`, `inst/extdata/examples/edger_witness`, and
`exec/annotation-breakpoint`.

## 2. Run exact Layer A from the CLI

Installed package:

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"
AB_EXACT="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "exact_toy", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EXACT" \
  --layer-a \
  --output-dir ab-output/tutorial-exact \
  --method exact
```

Source tree:

```sh
Rscript --vanilla exec/annotation-breakpoint \
  --input-dir inst/extdata/examples/exact_toy \
  --layer-a \
  --output-dir ab-output/tutorial-exact-source \
  --method exact
```

The example should report `status=exact_minimum`,
`best_verified_cost=1`, `exact_minimum_cost=1`, and
`verification_ok=TRUE`. The original strict claim holds; one permitted label
change crosses onto the closed overturn side. Every cheaper feasible assignment
has therefore been ruled out for this exact problem.

Layer A uses only the registered feature universe. For each focal-population
sample it computes

```text
log2(1e6 * (target + 0.5) / (library + 0.5 * feature_count))
```

and defines `beta` as the equally weighted case mean minus control mean. It is
a transparent small-domain oracle, not the authoritative RNA-seq backend.

## 3. Run the authoritative backend

The authoritative example requires exactly R 4.5.0 and edgeR 4.6.3. The YAML
file freezes every numerical analysis option.

Installed package:

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"
AB_EDGER="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "edger_witness", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EDGER" \
  --analysis-spec "$AB_EDGER/analysis_spec.yaml" \
  --output-dir ab-output/tutorial-edger \
  --method proxy_prefix \
  --max-fits 32 \
  --max-proposals 256 \
  --max-seconds 60 \
  --seed 82901 \
  --direction all
```

Source tree:

```sh
Rscript --vanilla exec/annotation-breakpoint \
  --input-dir inst/extdata/examples/edger_witness \
  --analysis-spec inst/extdata/examples/edger_witness/analysis_spec.yaml \
  --output-dir ab-output/tutorial-edger-source \
  --method proxy_prefix \
  --max-fits 32 \
  --max-proposals 256 \
  --max-seconds 60 \
  --seed 82901 \
  --direction all
```

This should report `status=verified_counterexample`,
`best_verified_cost=1`, `exact_minimum_cost=NA`, and
`verification_ok=TRUE`. `NA` is required for the exact minimum: the practical
engine found a
verified upper bound but did not enumerate the finite domain.

For every attempted assignment the backend:

1. rebuilds focal-label donor pseudobulks from immutable raw counts;
2. rebuilds the intercept-plus-condition design with the declared control and
   case strings;
3. calls `filterByExpr` with `min.count=10`, `min.total.count=15`,
   `large.n=10`, and `min.prop=0.7`;
4. requires the focal gene and at least 50 features to survive;
5. resets post-filter library sizes and recalculates TMM factors;
6. refits `glmQLFit` with abundance trend enabled, `robust=FALSE`,
   `legacy=FALSE`, and `prior.count=0.125`; and
7. obtains the case-versus-control log2 effect from `glmQLFTest` with
   `poisson.bound=TRUE`.

No fitted nuisance quantity is shared between assignments. Warnings are
recorded. A backend error, nonfinite output, filtered target, insufficient
features, bad normalization, or design failure is not an overturn.

## 4. Run the same examples through the R API

The file readers are internal convenience functions used by the CLI. The
analysis API begins with exported `validate_inputs()`.

```r
library(AnnotationBreakpoint)

read_input_directory <- getFromNamespace(
  "read_input_directory", "AnnotationBreakpoint"
)

exact_inputs <- read_input_directory(exact_dir)
exact_validation <- do.call(
  validate_inputs,
  exact_inputs[c(
    "counts", "cells", "specimens", "allowed_labels", "claim",
    "constraints", "analysis_spec"
  )]
)
stopifnot(exact_validation$ok)

exact_result <- enumerate_breakpoint(exact_validation)
stopifnot(
  identical(exact_result$status, "exact_minimum"),
  identical(exact_result$exact_minimum_cost, 1L)
)

exact_replay <- verify_witness(exact_validation, exact_result)
stopifnot(exact_replay$ok)

exact_artifacts <- render_report(
  exact_result, exact_replay, "ab-output/r-api-exact",
  run_metadata = list(example = "exact_toy", interface = "R")
)
```

For edgeR:

```r
library(yaml)

edger_spec <- yaml::read_yaml(file.path(edger_dir, "analysis_spec.yaml"))
edger_inputs <- read_input_directory(edger_dir, analysis_spec = edger_spec)
edger_validation <- do.call(
  validate_inputs,
  edger_inputs[c(
    "counts", "cells", "specimens", "allowed_labels", "claim",
    "constraints", "analysis_spec"
  )]
)
stopifnot(edger_validation$ok)

edger_result <- search_breakpoint(
  edger_validation,
  method = "proxy_prefix",
  max_fits = 32,
  seed = 82901L,
  max_proposals = 256,
  max_seconds = 60,
  direction = "all"
)
stopifnot(
  identical(edger_result$status, "verified_counterexample"),
  identical(edger_result$best_verified_cost, 1L),
  is.na(edger_result$exact_minimum_cost)
)

edger_replay <- verify_witness(edger_validation, edger_result)
stopifnot(
  edger_replay$ok,
  identical(edger_replay$recomputed_cost, 1L)
)

edger_artifacts <- render_report(
  edger_result, edger_replay, "ab-output/r-api-edger",
  run_metadata = list(example = "edger_witness", interface = "R")
)
```

`verify_witness()` bypasses the practical-search cache. It checks immutable
input hashes, allowed-label membership, all budgets and minima, assignment
identity, cost, and a fresh original/witness fit. It returns the changed-cell
provenance and before/after donor pseudobulks.

## 5. Understand the input contract

`counts.tsv` starts with `gene_id`; every remaining header is a cell ID.
Counts must be finite nonnegative integers, with unique feature and cell IDs.
In the R API, counts are a feature-by-cell `Matrix::dgCMatrix` with matching
dimnames.

`cells.tsv` has:

```text
cell_id  sample_id  donor_id  condition  original_label
```

`specimens.tsv` has:

```text
sample_id  donor_id  condition
```

Every cell's immutable sample/donor/condition tuple must match the specimen
table. Version 1 accepts one sample per donor and two declared conditions. Each
condition needs at least two independent donors; paired/repeated specimens are
not supported.

`allowed_labels.tsv` has one row per cell/allowed-label pair:

```text
cell_id  allowed_label  evidence_source  evidence_version  [calibrated_confidence]
```

The original label must appear. A permitted alternative needs a concrete,
versioned provenance record. It should be built without looking at breakpoint
outputs. A confidence value is optional and checked for `[0,1]`, but it does
not alter version 1 unit cost. Search results remain conditional on the allowed
set even when the biologically correct label is missing from it.

The claim row names the focal gene, focal and alternative labels, case and
control strings, `q`, `delta`, `tolerance`, the pipe-delimited feature universe,
and the estimator ID. The original claim is strict:

```text
q * beta > delta + tolerance
```

Equality belongs to the overturn side. If the original assignment does not
satisfy the strict inequality, the status is `original_claim_not_met` rather
than a zero-cost witness.

The constraints row supplies `min_cells_<label>`, `min_library_<label>`,
`max_changes_global`, `max_changes_per_sample`, `max_changes_per_donor`, and
`min_donors_per_condition`. Scalar sample/donor budgets apply to every
registered ID through the file interface; the R API can also use named maps.

The minimum-library constraint is checked per sample/label over the registered
feature universe before backend filtering. It prevents the optimizer from
manufacturing an empty or near-empty aggregate. `filterByExpr` then defines the
backend's retained universe and post-filter library sizes. These are different
checks.

## 6. Interpret search modes

- `exact` considers the declared finite allowed-label domain, up to 12
  ambiguous cells and 4,096 assignments. `exact_minimum` rules out every
  cheaper feasible assignment. `exhaustively_stable` requires every feasible
  assignment to fit successfully without an overturn.
- `proxy_prefix`, the software default, is labeled **observed-centered
  proxy-ranked cumulative prefixes**. It deterministically ranks permitted
  changes and fully refits exact cumulative prefixes of
  `1, 2, 4, 8, 16, 32, 64, 128`. A smaller domain uses only available
  registered prefixes. An inadmissible exact prefix is rejected without
  skipping or filling. Every available admissible prefix is fully fitted,
  including those after an earlier crossing; the earliest crossing prefix in
  registered order that passes replay is an upper bound. A no-hit remains
  unresolved. Normal use omits benchmark-only singleton diagnostics.
- `ranked` is the separate adaptive historical engine. It is available only by
  explicit R API selection for reproducibility and is not a CLI method.
- `greedy` evaluates feasible single-label neighbors and takes the cheapest
  overturn or a strictly margin-improving step. It remains available only by
  explicit selection for reproducibility.
- `pair` adds two-cell local moves and swaps. It retained no hit/cost advantage
  over greedy in M3 and used about twice as many search fits.
- `random` samples allowed integer assignments with a local seed. A hit is a
  verified upper bound; a no-hit remains unresolved.
- `proxy_prefix` rejects `outgoing_only`. For compatible explicitly selected
  legacy methods, `outgoing_only` permits only focal-to-alternative changes. It
  is restricted coverage of a different feasible set, not a bidirectional
  exact comparison. The CLI method choices are `proxy_prefix`, `exact`,
  `random`, `greedy`, and `pair`.

Every method uses Hamming distance from the original assignment, so a reversal
can lower current cost. Practical caps include the original fit in
`max_fits`. Duplicate proposals and fresh full fits are recorded separately.

## 7. Inspect the accepted public examples without refitting

The accepted saved comparator artifacts contain two public `proxy_prefix`
witnesses for the COMBAT CD4.TCM expression contrast COVID_CRIT versus
COVID_SEV, with 16 donors (8 per group), CD4.TCM/CD4.NAIVE alternatives, and
`K=8,898`. They use the original global budget 445, frozen 5% per-donor caps,
and effect threshold 0.5:

| Gene | Original beta | Replayed beta | Cost |
| --- | ---: | ---: | ---: |
| FAM118A | 0.8106752561679119 | 0.4731965904097605 | 64 |
| PDXK | 0.7274641768393035 | 0.3373883736018561 | 32 |

From a source checkout:

```sh
awk -F '\t' 'NR == 1 || (($1 == "REAL_FAM118A" || $1 == "REAL_PDXK") && $3 == "B")' \
  reports/DECISIVE-COMPARATOR-001/results_13x3.tsv
```

The prefix-level effects are in `method_b_prefixes.tsv`, and the independent
replays are in the owner receipts named in the repository README. These costs
are verified upper bounds, not exact minima. Crossing 0.5 concerns the fitted
expression effect; it is not a significance result, annotation truth, or a
biological discovery. The accepted project disposition found no demonstrated
unique solver contribution, retained the negative Hao result, and left the
scaling and broad public-utility questions inconclusive.

## 8. Replay and compare outputs

The CLI always calls `verify_witness()` before rendering. For a second replay,
rerun the same command with a new output directory. Compare semantic output in
R:

```r
left <- read.delim("ab-output/tutorial-edger/claim_results.tsv")
right <- read.delim("ab-output/tutorial-edger-replay/claim_results.tsv")

fields <- c(
  "status", "original_beta_log2", "best_verified_cost",
  "exact_minimum_cost", "witness_id"
)
stopifnot(identical(left[fields], right[fields]))

left_hashes <- read.delim("ab-output/tutorial-edger/input_hashes.tsv")
right_hashes <- read.delim("ab-output/tutorial-edger-replay/input_hashes.tsv")
stopifnot(identical(left_hashes, right_hashes))

left_moves <- read.delim("ab-output/tutorial-edger/witnesses.tsv")
right_moves <- read.delim("ab-output/tutorial-edger-replay/witnesses.tsv")
stopifnot(identical(left_moves, right_moves))
```

Search and verification timing fields may differ. Consequently, whole-file
hashes of timing-bearing JSON or TSV artifacts need not match even when the
registered result replays exactly.

`witnesses.tsv` preserves both `original_label` and `assigned_label`; it never
overwrites the supplied annotation. `verification.json` records recomputed
cost/effect/status and hash matches. `before_pseudobulks.tsv` and
`after_pseudobulks.tsv` show the donor-level count transfer.

## 9. Invalid input and failures

Validation is inspectable without running a search:

```r
bad <- exact_inputs
bad$allowed_labels <- bad$allowed_labels[
  !(bad$allowed_labels$cell_id == bad$cells$cell_id[[1L]] &
      bad$allowed_labels$allowed_label == bad$cells$original_label[[1L]]),
]
bad_validation <- do.call(
  validate_inputs,
  bad[c(
    "counts", "cells", "specimens", "allowed_labels", "claim",
    "constraints", "analysis_spec"
  )]
)
stopifnot(
  !bad_validation$ok,
  identical(bad_validation$status, "invalid_input")
)
bad_validation$errors
```

The CLI exits with code 2 for validated `invalid_input` and code 3 for failed
witness verification. Parse, missing-file, dependency, or artifact-write errors
are ordinary nonzero R process failures. Scientific statuses such as
`original_claim_not_met`, `not_assessable`, `no_counterexample_found`, and
`exhaustive_search_incomplete` are written as machine-readable results and do
not become counterexamples.

An unevaluable feasible assignment is retained in the assignment/method log
with its warning, reason, or error. It is never counted as an overturn. A
failed cheaper assignment blocks an exact-minimum conclusion; any failed
feasible assignment blocks unqualified exhaustive stability.

## 10. Limits carried into interpretation

Results are conditional on the supplied raw counts, the upstream QC/cell set,
the registered feature universe, the allowed-label construction, constraints,
and analysis policy. The tool cannot recover cells removed before input or
guarantee that the correct alternative is allowed.

Per-gene searches may return different assignments. Do not merge those
witnesses into one alternative annotation. A joint gene-set claim needs a new
common-assignment registration.

The held-out benchmark had 132 synthetic realizations, six ambiguous cells,
and cost at most three. The public demonstration used six Hao5P donors and an
artificial condition split; no public witness was found. S1, scANVI, and
scDEcrypter were not executed, and this workflow does not reproduce them.
