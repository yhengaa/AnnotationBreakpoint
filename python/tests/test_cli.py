from __future__ import annotations

import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

from annotation_breakpoint.bridge import BackendInfo, RunResult
from annotation_breakpoint.cli import main


class CliTests(unittest.TestCase):
    def test_cli_emits_structured_summary_and_passes_limits(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = RunResult(
                result={
                    "status": "no_counterexample_found",
                    "result_id": "result-id",
                    "witness_id": None,
                    "witness": None,
                    "best_verified_cost": None,
                    "exact_minimum_cost": None,
                    "input_hashes": {"input_hash": "hash"},
                    "backend_identity": {"backend": "edgeR"},
                    "original_fit": {"runtime": {"policy": "strict"}},
                },
                verification={"ok": True},
                manifest={},
                artifacts={"result_json": root / "result.json"},
                backend=BackendInfo(root / "Rscript", root / "cli", "0.3.0"),
                command=("Rscript",),
                stdout="",
                stderr="",
                output_dir=root / "output",
            )
            stdout = io.StringIO()
            with mock.patch(
                "annotation_breakpoint.cli.AnnotationBreakpointBackend.run",
                return_value=result,
            ) as invoked, redirect_stdout(stdout):
                code = main(
                    [
                        "--input-dir",
                        str(root),
                        "--analysis-spec",
                        str(root / "spec.yaml"),
                        "--output-dir",
                        str(root / "output"),
                        "--prefix-sizes",
                        "1,3,8",
                        "--max-fits",
                        "21",
                    ]
                )
            self.assertEqual(code, 0)
            payload = json.loads(stdout.getvalue())
            self.assertEqual(payload["status"], "no_counterexample_found")
            self.assertFalse(payload["witness_found"])
            arguments = invoked.call_args.kwargs
            self.assertEqual(arguments["prefix_sizes"], (1, 3, 8))
            self.assertEqual(arguments["max_fits"], 21)

    def test_cli_serializes_bridge_errors_to_stderr(self) -> None:
        stderr = io.StringIO()
        with tempfile.TemporaryDirectory() as directory, mock.patch(
            "annotation_breakpoint.cli.AnnotationBreakpointBackend.run",
            side_effect=ValueError("controlled input error"),
        ), redirect_stderr(stderr):
            code = main(
                [
                    "--input-dir",
                    directory,
                    "--layer-a",
                    "--output-dir",
                    str(Path(directory) / "output"),
                ]
            )
        self.assertEqual(code, 2)
        payload = json.loads(stderr.getvalue())
        self.assertEqual(payload["error_type"], "ValueError")
        self.assertIn("controlled input error", payload["message"])


if __name__ == "__main__":
    unittest.main()
