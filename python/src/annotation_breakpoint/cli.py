"""Command-line entry point for the Python bridge."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Sequence

from .bridge import AnnotationBreakpointBackend, BridgeError


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="annotation-breakpoint-py",
        description="Run the installed Annotation Breakpoint 0.3.0 R engine.",
    )
    parser.add_argument("--input-dir", required=True)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--analysis-spec")
    target.add_argument("--layer-a", action="store_true")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method", choices=("proxy_prefix", "exact"), default="proxy_prefix")
    parser.add_argument("--allow-version-mismatch", action="store_true")
    parser.add_argument("--prefix-sizes", type=_parse_prefix_sizes)
    parser.add_argument("--max-evaluations", type=_parse_number, default=float("inf"))
    parser.add_argument("--max-fits", type=_parse_number, default=1024)
    parser.add_argument("--max-proposals", type=_parse_number, default=float("inf"))
    parser.add_argument("--max-seconds", type=_parse_number, default=float("inf"))
    parser.add_argument("--seed", type=int, default=82901)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--rscript")
    parser.add_argument("--r-cli")
    return parser


def _parse_number(value: str) -> int | float:
    if value.lower() in {"inf", "infinity"}:
        return float("inf")
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a nonnegative integer or Inf") from error
    if number < 0:
        raise argparse.ArgumentTypeError("expected a nonnegative integer or Inf")
    return number


def _parse_prefix_sizes(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated positive integers") from error
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return values


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    backend = AnnotationBreakpointBackend(rscript=args.rscript, r_cli=args.r_cli)
    try:
        run = backend.run(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            analysis_spec=args.analysis_spec,
            layer_a=args.layer_a,
            method=args.method,
            allow_version_mismatch=args.allow_version_mismatch,
            prefix_sizes=args.prefix_sizes,
            max_evaluations=args.max_evaluations,
            max_fits=args.max_fits,
            max_proposals=args.max_proposals,
            max_seconds=args.max_seconds,
            seed=args.seed,
            timeout=args.timeout,
        )
    except (BridgeError, ValueError) as error:
        payload = error.as_dict() if isinstance(error, BridgeError) else {
            "error_type": type(error).__name__,
            "message": str(error),
        }
        print(json.dumps(payload, sort_keys=True), file=sys.stderr)
        return 2
    print(json.dumps(run.summary(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
