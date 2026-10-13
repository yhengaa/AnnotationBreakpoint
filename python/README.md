# Annotation Breakpoint for Python

This package is a Python interface to the installed
`AnnotationBreakpoint` 0.3.0 R package. It launches the validated R CLI with an
argument array, preserves its JSON and TSV outputs, and checks the emitted
artifact hashes and result/replay identities before returning.

It requires Python 3.10 or later, R 4.5.x, edgeR 4.6.3, and the
AnnotationBreakpoint R package 0.3.0. It does not reimplement the estimator in
Python.

Install from the repository subdirectory:

```sh
python -m pip install \
  "annotation-breakpoint-python @ git+https://github.com/yhengaa/AnnotationBreakpoint.git#subdirectory=python"
```

Run the authoritative edgeR workflow:

```sh
annotation-breakpoint-py \
  --input-dir inputs/FAM118A \
  --analysis-spec analysis_specs/FAM118A.yaml \
  --output-dir results/FAM118A
```

Or call the API:

```python
from annotation_breakpoint import run_annotation_breakpoint

run = run_annotation_breakpoint(
    input_dir="inputs/FAM118A",
    analysis_spec="analysis_specs/FAM118A.yaml",
    output_dir="results/FAM118A",
)
print(run.status, run.witness_found)
print(run.result)        # unchanged result.json document
print(run.verification)  # unchanged verification.json document
```

`witness_found` is separate from process and replay success. A valid
`no_counterexample_found`, `original_claim_not_met`, or `not_assessable` run
therefore remains distinct from a verified witness.

See the repository [Python guide](../docs/PYTHON.md) for exact search,
compatible runtime policy, explicit executable paths, limits, and error
handling.
