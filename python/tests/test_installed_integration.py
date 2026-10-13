from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from annotation_breakpoint import AnnotationBreakpointBackend


@unittest.skipUnless(
    os.environ.get("ANNOTATION_BREAKPOINT_RUN_INTEGRATION") == "1",
    "set ANNOTATION_BREAKPOINT_RUN_INTEGRATION=1 for the tiny edgeR parity test",
)
class InstalledEdgeRParityTests(unittest.TestCase):
    def test_python_matches_direct_r_cli(self) -> None:
        rscript = shutil.which("Rscript")
        self.assertIsNotNone(rscript)
        discovery = subprocess.run(
            [
                rscript,
                "--vanilla",
                "-e",
                'cat(system.file("exec","annotation-breakpoint",package="AnnotationBreakpoint"),"\\n"); cat(system.file("extdata","examples","edger_witness",package="AnnotationBreakpoint"),"\\n")',
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        r_cli, fixture = (line.strip() for line in discovery.stdout.splitlines()[:2])
        self.assertTrue(r_cli)
        self.assertTrue(fixture)
        analysis_spec = str(Path(fixture) / "analysis_spec.yaml")
        with tempfile.TemporaryDirectory(prefix="ab-python-installed-parity-") as directory:
            root = Path(directory)
            direct_dir = root / "direct R output"
            python_dir = root / "Python output"
            common = [
                "--input-dir",
                fixture,
                "--analysis-spec",
                analysis_spec,
                "--method",
                "proxy_prefix",
                "--max-fits",
                "32",
                "--seed",
                "82901",
            ]
            direct = subprocess.run(
                [rscript, "--vanilla", r_cli, *common, "--output-dir", str(direct_dir)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(direct.returncode, 0, direct.stderr or direct.stdout)
            python_run = AnnotationBreakpointBackend(rscript=rscript, r_cli=r_cli).run(
                input_dir=fixture,
                analysis_spec=analysis_spec,
                output_dir=python_dir,
                method="proxy_prefix",
                max_fits=32,
                seed=82901,
            )
            direct_result = json.loads((direct_dir / "result.json").read_text())
            direct_verification = json.loads((direct_dir / "verification.json").read_text())
            self.assertEqual(python_run.status, direct_result["status"])
            self.assertEqual(
                python_run.result["best_verified_cost"], direct_result["best_verified_cost"]
            )
            self.assertEqual(python_run.result["witness"], direct_result["witness"])
            self.assertEqual(python_run.result["input_hashes"], direct_result["input_hashes"])
            self.assertEqual(
                python_run.result["original_fit"]["beta"], direct_result["original_fit"]["beta"]
            )
            self.assertEqual(
                python_run.verification["recomputed_effect"],
                direct_verification["recomputed_effect"],
            )
            self.assertEqual(
                python_run.verification["result_identity"], python_run.result["result_id"]
            )
            self.assertEqual(
                direct_verification["result_identity"], direct_result["result_id"]
            )
            self.assertEqual(
                python_run.verification["witness_id"], direct_verification["witness_id"]
            )
            self.assertEqual(
                python_run.verification["evaluation_identity"],
                direct_verification["evaluation_identity"],
            )
            self.assertTrue(python_run.verification_ok)
            self.assertTrue(python_run.witness_found)


if __name__ == "__main__":
    unittest.main()
