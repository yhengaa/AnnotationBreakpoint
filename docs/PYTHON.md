# Python interface

`annotation-breakpoint-python` 0.1.0 is a standard-library Python interface to
the installed `AnnotationBreakpoint` 0.3.0 R engine. The R package still reads
the inputs, fits edgeR models, searches assignments, and independently replays
witnesses. Python launches the installed R CLI without a shell and validates
the fresh result artifacts before returning them.

## Requirements and installation

- Python 3.10 or later
- R 4.5.x and edgeR 4.6.3
- the AnnotationBreakpoint R package 0.3.0

Install the R package as described in the main [README](../README.md), then
install the Python subdirectory from GitHub:

```sh
python -m pip install \
  "annotation-breakpoint-python @ git+https://github.com/yhengaa/AnnotationBreakpoint.git#subdirectory=python"
```

The Python package has no third-party runtime dependencies.

## Command line

Run the bundled portable toy example immediately after installing both
packages:

```sh
AB_EXAMPLE="$(Rscript --vanilla -e 'cat(system.file("extdata", "examples", "exact_toy", package="AnnotationBreakpoint"))')"
annotation-breakpoint-py \
  --input-dir "$AB_EXAMPLE" \
  --layer-a \
  --output-dir ab-output/python-quickstart
```

The default fixed-prefix search uses the same defaults as the R CLI:

```sh
annotation-breakpoint-py \
  --input-dir inputs/FAM118A \
  --analysis-spec analysis_specs/FAM118A.yaml \
  --output-dir results/FAM118A
```

Use an explicit prefix schedule or exact enumeration as follows:

```sh
annotation-breakpoint-py \
  --input-dir inputs/FAM118A \
  --analysis-spec analysis_specs/FAM118A.yaml \
  --output-dir results/FAM118A-custom \
  --prefix-sizes 1,2,4,8,16,32,64,128,256 \
  --max-fits 1024 --max-proposals Inf --max-seconds Inf

annotation-breakpoint-py \
  --input-dir inputs/tiny \
  --layer-a \
  --output-dir results/tiny-exact \
  --method exact --max-evaluations Inf
```

For an authoritative analysis on a compatible R 4.5.x/edgeR 4.6.x patch
runtime, add `--allow-version-mismatch`. This is forwarded to the R engine and
recorded in its runtime identity. It does not substitute another estimator for
an unsupported edgeR API.

Discovery normally uses `Rscript` from `PATH` and the CLI installed inside the
R package. Use `--rscript`, `--r-cli`,
`ANNOTATION_BREAKPOINT_RSCRIPT`, or `ANNOTATION_BREAKPOINT_R_CLI` for explicit
paths. Paths containing spaces are passed as individual arguments.

`--timeout` is a Python-side wall-time limit. The R search limits
`--max-fits`, `--max-proposals`, and `--max-seconds` remain part of the saved
engine record. A timeout or nonzero R exit returns a nonzero Python exit and
includes captured backend output on standard error.

The command prints a JSON summary. The complete, unmodified R documents remain
in `result.json` and `verification.json` within the output directory.

## Python API

```python
from annotation_breakpoint import AnnotationBreakpointBackend

backend = AnnotationBreakpointBackend()
run = backend.run(
    input_dir="inputs/FAM118A",
    analysis_spec="analysis_specs/FAM118A.yaml",
    output_dir="results/FAM118A",
    method="proxy_prefix",
    prefix_sizes=(1, 2, 4, 8, 16, 32, 64, 128),
    max_fits=1024,
    timeout=3600,
)

print(run.status)
print(run.witness_found)
print(run.result)
print(run.verification)
print(run.artifacts)
```

The caller must choose a new output path for every invocation. Existing paths
are rejected before R starts, so older reports cannot be mistaken for the
current run. On success the bridge checks every artifact hash in
`run_manifest.json`. Status and input hashes must agree across the result,
verification, and manifest; result and witness identities must agree between
the result and its independent replay.

The `verification_ok` property reports successful independent replay of the
saved status. The `witness_found` property is true only for a replayed witness.
Consequently, valid `no_counterexample_found`, `original_claim_not_met`,
`not_assessable`, `exhaustive_search_incomplete`, and `exhaustively_stable`
results stay distinct. An incomplete exact result can still retain a replayed
upper-bound witness when failed feasible fits block a minimum certificate.

Exceptions preserve useful failure information:

- `BackendDiscoveryError`: Rscript or the installed R CLI cannot be found.
- `UnsupportedBackendError`: the R package is not version 0.3.0.
- `BackendProcessError`: R exits nonzero; stdout, stderr, command, and partial
  output path remain available.
- `BackendTimeoutError`: the Python wall-time limit expires, with any captured
  streams and partial output path retained.
- `StaleOutputError`: the requested output path already exists.
- `ArtifactError`: fresh output is missing, malformed, hash-invalid, or has
  inconsistent identities.

The convenience function `run_annotation_breakpoint()` accepts the same run
arguments plus optional `rscript`, `r_cli`, `env`, and `discovery_timeout`
backend settings.
