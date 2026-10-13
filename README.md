# Annotation Breakpoint

Annotation Breakpoint is an R utility for measuring how sensitive a predeclared donor-level expression-effect claim is to allowed cell-type label changes. It returns a replay-verified changed assignment when one crosses the specified effect threshold and distinguishes exact minima, finite-domain stability, verified upper bounds, and unresolved no-hit searches.

A witness shows sensitivity under the supplied labels, allowed alternatives, constraints, and analysis specification. It does not identify the biologically correct annotation, prove that a moved cell was mislabeled, or establish a change in statistical significance.

## Install

Download the 0.2.0 source package from the [GitHub release](https://github.com/yhengaa/AnnotationBreakpoint/releases/tag/v0.2.0):

```sh
curl -L -O https://github.com/yhengaa/AnnotationBreakpoint/releases/download/v0.2.0/AnnotationBreakpoint_0.2.0.tar.gz
echo "17e57844e6c8abad4147032bf87e1d3d0d0e03776f5226e5992be56eba0d04d3  AnnotationBreakpoint_0.2.0.tar.gz" | sha256sum -c -
R CMD INSTALL AnnotationBreakpoint_0.2.0.tar.gz
```

The package requires R 4.5.0 or later. Its authoritative pseudobulk backend requires exactly edgeR 4.6.3; a version mismatch is returned as `not_assessable`. Other dependencies are listed in [`DESCRIPTION`](DESCRIPTION).

## Quick start

The installed `exact_toy` fixture is a small transparent example:

```sh
AB_CLI="$(Rscript --vanilla -e 'cat(system.file("exec", "annotation-breakpoint", package="AnnotationBreakpoint"))')"
AB_EXACT="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "exact_toy", package="AnnotationBreakpoint"))')"

Rscript --vanilla "$AB_CLI" \
  --input-dir "$AB_EXACT" \
  --layer-a \
  --output-dir ab-output/exact-toy \
  --method exact
```

The example should report `status=exact_minimum`, `best_verified_cost=1`, `exact_minimum_cost=1`, and `verification_ok=TRUE`. The [tutorial](inst/doc/annotation-breakpoint-tutorial.md) also covers the authoritative edgeR fixture, the R API, failure handling, and replay.

## Default practical search

The default `proxy_prefix` method uses an observed-centered, proxy-ranked prefix search:

1. Rank each permitted single-cell label change once using a finite donor-logCPM change score at the observed annotation.
2. Evaluate cumulative prefixes of sizes `1, 2, 4, 8, 16, 32, 64, 128` without reranking.
3. Reject an inadmissible exact prefix rather than repairing it or filling it with lower-ranked cells.
4. Subject to the caller's fit and time limits, rebuild the complete label-dependent analysis and perform a full authoritative fit for every admissible available prefix.
5. Independently replay the smallest evaluated prefix that crosses the threshold and report it only when replay passes.

The proxy is an AnnotationBreakpoint finite-change adaptation. It is not AMIP, a formal first-order edgeR influence approximation, or a reproduction of a published annotation method. A replayed result from this practical search is an upper bound on breakpoint cost; it is not a certified minimum. A no-hit remains unresolved.

Exact enumeration is available for domains with at most 12 ambiguous cells and 4,096 allowed assignments. Explicit legacy `random`, `greedy`, and `pair` methods remain available in the CLI; the separate adaptive `ranked` method remains available through the R API for reproducibility.

## Claim and input contract

For fitted log2 effect `beta`, direction `q` in `{-1, +1}`, threshold `delta`, and numerical `tolerance`, the original claim is

```text
q * beta > delta + tolerance
```

An evaluable changed assignment overturns the claim when

```text
q * beta <= delta + tolerance
```

This is an effect-threshold target, not a p-value threshold or a test of annotation correctness. Cost is the number of cells whose assigned label differs from its original label; counts and specimen membership never change.

A normal input directory contains six tab-separated files:

| File | Content |
| --- | --- |
| `counts.tsv` | Integer raw counts, with `gene_id` followed by cell columns |
| `cells.tsv` | Cell, sample, donor, condition, and original-label metadata |
| `specimens.tsv` | Sample-to-donor and condition mapping |
| `allowed_labels.tsv` | Allowed labels and versioned evidence for each cell |
| `claim.tsv` | Target gene, populations, conditions, direction, threshold, feature universe, and estimator |
| `constraints.tsv` | Global/sample/donor change limits, population and library minima, and donor replication minimum |

Choose one analysis target explicitly: `--layer-a` for the transparent fixture estimator, or `--analysis-spec` for the authoritative edgeR policy. The full schema is in [Tutorial §5](inst/doc/annotation-breakpoint-tutorial.md#5-understand-the-input-contract).

## Outputs and replay

The CLI writes the claim result, fitted prefix or search log, fit accounting, changed-cell witness, before/after donor pseudobulks, input hashes, rendered constraints and analysis specification, JSON records, and a Markdown report. It performs an independent replay from immutable inputs before reporting a witness.

`verification_ok=TRUE` means the assignment identity, allowed-label membership, constraints, unit cost, input hashes, and fresh fitted effect agreed. See [Tutorial §8](inst/doc/annotation-breakpoint-tutorial.md#8-replay-and-compare-outputs) for a replay comparison.

Important result statuses are:

| Status | Meaning |
| --- | --- |
| `exact_minimum` | A verified witness exists and every cheaper feasible assignment was ruled out |
| `verified_counterexample` | A replayed witness gives an upper bound; no minimum is certified |
| `exhaustively_stable` | Every feasible assignment in the declared finite domain fitted successfully and none crossed the threshold |
| `no_counterexample_found` | No witness was found within the reported practical search; unresolved |
| `original_claim_not_met` | The original labels did not satisfy the strict claim |
| `not_assessable` | The requested analysis could not produce an evaluable result |
| `invalid_input` | Static schema, identity, provenance, design, or feasibility checks failed |

## Saved public examples

The repository includes two replay-verified `proxy_prefix` examples from the COMBAT CD4.TCM expression contrast `COVID_CRIT` versus `COVID_SEV`. The inputs contain 16 independent donors (8 per group), CD4.TCM/CD4.NAIVE alternatives, and 8,898 ambiguous cells. Both examples keep the 0.5 effect threshold, the 445-cell global limit, and the original per-donor constraints.

| Gene | Original beta | Replayed witness beta | Cost | Interpretation |
| --- | ---: | ---: | ---: | --- |
| FAM118A | 0.8106752561679119 | 0.4731965904097605 | 64 | verified upper bound |
| PDXK | 0.7274641768393035 | 0.3373883736018561 | 32 | verified upper bound |

The complete prefix effects are in [`method_b_prefixes.tsv`](reports/DECISIVE-COMPARATOR-001/method_b_prefixes.tsv), and the 39-row comparison is in [`results_13x3.tsv`](reports/DECISIVE-COMPARATOR-001/results_13x3.tsv). Neither reported cost is a proven minimum. These development examples do not establish biological utility or annotation truth.

## Evaluation and limits

A fixed 13-case development comparison used ten eligibility-screened synthetic cases and three COMBAT development claims. `proxy_prefix` found replayed witnesses in 13/13 cases and canonical greedy in 11/13. Canonical greedy had zero exclusive discoveries and no cost advantage of at least two cells. This result did not support a unique search-method contribution. It also does not estimate held-out performance: the synthetic cases were eligibility-screened, and the public claims reused development data. The detailed [comparator result note](reports/DECISIVE-COMPARATOR-001/RESULT_NOTE.md) preserves the pre-specified decision and accounting.

Earlier evaluations remain deliberately visible:

- In the scaling evaluation, the primary K=192 comparison had only ten exact-attainable parents, below its pre-specified minimum denominator of 12, so its conclusion was **inconclusive**. Search coverage changed substantially with the fit allowance.
- In the COMBAT public-utility evaluation, only three claims originally met the effect threshold, below its minimum denominator of ten. There were zero strict added-information claims, so the utility conclusion was **inconclusive**.
- The [Hao5P demonstration](reports/M4-PUBLIC-001/NEGATIVE_HAO.md) found no allowed counterexample among the six originally holding NME2/KLF2 finite domains. That negative result applies only to the registered six-cell domains and is not a population-wide robustness claim.

See the compact [evaluation evidence summary](reports/FINAL-EVALUATION-001/EVIDENCE_SUMMARY.md) and [evidence index](reports/PUBLICATION/README.md). These results support the implementation and its replay/accounting behavior. They do not establish broad practical superiority, a scaling guarantee, or biological utility.

## Scope and license

The current package supports one sample per biological donor, two declared conditions, and an intercept-plus-condition design. Paired or repeated specimens are not supported. Each witness is specific to one gene and one declared problem; witnesses for different genes may be incompatible and must not be combined into a joint annotation.

The repository excludes raw public counts, caches, multi-gigabyte model objects, and execution logs. Compact numerical evidence is included under [`reports/`](reports/PUBLICATION/README.md).

AnnotationBreakpoint is distributed under the [MIT License](LICENSE.md). Dependency and attribution details are in [`inst/NOTICE`](inst/NOTICE).
