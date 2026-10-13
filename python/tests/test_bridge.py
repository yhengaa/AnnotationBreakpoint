from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from annotation_breakpoint import (
    AnnotationBreakpointBackend,
    ArtifactError,
    BackendDiscoveryError,
    BackendProcessError,
    BackendTimeoutError,
    StaleOutputError,
    UnsupportedBackendError,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _documents(
    status: str = "verified_counterexample", *, force_witness: bool | None = None
) -> tuple[dict, dict]:
    witness = (
        status in {"verified_counterexample", "exact_minimum"}
        if force_witness is None
        else force_witness
    )
    result = {
        "status": status,
        "result_id": "result-identity",
        "witness_id": "witness-identity" if witness else None,
        "witness": {"cost": 1, "assignment": [{"cell_id": "C 1"}]} if witness else None,
        "best_verified_cost": 1 if witness else None,
        "exact_minimum_cost": 1 if status == "exact_minimum" else None,
        "input_hashes": {"counts": "counts-hash", "input_hash": "input-hash"},
        "backend_identity": {"backend": "edgeR", "backend_version": "4.6.3"},
        "original_fit": {"runtime": {"policy": "strict", "observed_r_version": "4.5.0"}},
    }
    verification = {
        "ok": True,
        "recomputed_status": status,
        "result_identity": "result-identity",
        "witness_id": "witness-identity" if witness else None,
        "input_hashes": dict(result["input_hashes"]),
        "hash_matches": [{"hash_name": "counts", "matches": True}],
    }
    return result, verification


def _write_output(
    target: Path,
    *,
    status: str = "verified_counterexample",
    mutate=None,
    malformed_result: bool = False,
    force_witness: bool | None = None,
) -> None:
    target.mkdir(parents=True, exist_ok=True)
    result, verification = _documents(status, force_witness=force_witness)
    if mutate:
        mutate(result, verification)
    payloads = {
        "result_json": ("result.json", "{" if malformed_result else json.dumps(result)),
        "verification_json": ("verification.json", json.dumps(verification)),
        "claim_results": ("claim_results.tsv", "status\n" + status + "\n"),
        "report": ("report.md", "# Report\n"),
    }
    artifacts = []
    for name, (filename, content) in payloads.items():
        path = target / filename
        path.write_text(content, encoding="utf-8")
        artifacts.append({"name": name, "file": filename, "sha256": _sha256(path)})
    manifest = {
        "status": status,
        "verification_ok": True,
        "run_metadata": {"package_version": "0.3.0"},
        "input_hashes": dict(result["input_hashes"]),
        "artifacts": artifacts,
    }
    (target / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


class BackendTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="ab python bridge ")
        self.root = Path(self.temp.name)
        self.cli = self.root / "R CLI with spaces"
        self.cli.write_text("fixture", encoding="utf-8")
        self.input_dir = self.root / "input with spaces"
        self.input_dir.mkdir()
        self.spec = self.input_dir / "analysis spec.yaml"
        self.spec.write_text("schema_version: 1.0.0\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _backend(self) -> AnnotationBreakpointBackend:
        return AnnotationBreakpointBackend(rscript=sys.executable, r_cli=self.cli)

    def _successful_side_effect(
        self, *, status="verified_counterexample", mutate=None, force_witness=None
    ):
        def invoke(command, **kwargs):
            self.assertIs(kwargs["shell"], False)
            if "-e" in command:
                return subprocess.CompletedProcess(
                    command, 0, f"0.3.0\n{self.cli}\n", ""
                )
            target = Path(command[command.index("--output-dir") + 1])
            _write_output(
                target,
                status=status,
                mutate=mutate,
                force_witness=force_witness,
            )
            return subprocess.CompletedProcess(command, 0, "backend stdout", "backend stderr")

        return invoke

    def test_shell_free_argv_paths_and_structured_witness(self) -> None:
        output = self.root / "output with spaces"
        calls = []

        def invoke(command, **kwargs):
            calls.append((command, kwargs))
            return self._successful_side_effect()(command, **kwargs)

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            run = self._backend().run(
                input_dir=self.input_dir,
                analysis_spec=self.spec,
                output_dir=output,
                prefix_sizes=(1, 3, 8),
                max_fits=12,
            )
        self.assertEqual(run.status, "verified_counterexample")
        self.assertTrue(run.witness_found)
        self.assertTrue(run.verification_ok)
        self.assertEqual(run.result["input_hashes"], run.verification["input_hashes"])
        command = calls[1][0]
        self.assertIn(str(self.input_dir.resolve()), command)
        self.assertIn(str(self.spec.resolve()), command)
        self.assertIn("1,3,8", command)
        self.assertFalse(calls[1][1]["shell"])

    def test_valid_nonwitness_statuses_are_distinct(self) -> None:
        for status in (
            "no_counterexample_found",
            "original_claim_not_met",
            "not_assessable",
            "exhaustively_stable",
            "exhaustive_search_incomplete",
        ):
            with self.subTest(status=status):
                output = self.root / ("output-" + status)
                with mock.patch(
                    "annotation_breakpoint.bridge.subprocess.run",
                    side_effect=self._successful_side_effect(status=status),
                ):
                    run = self._backend().run(
                        input_dir=self.input_dir,
                        analysis_spec=self.spec,
                        output_dir=output,
                    )
                self.assertEqual(run.status, status)
                self.assertTrue(run.verification_ok)
                self.assertFalse(run.witness_found)

    def test_incomplete_exact_result_can_retain_verified_upper_bound(self) -> None:
        with mock.patch(
            "annotation_breakpoint.bridge.subprocess.run",
            side_effect=self._successful_side_effect(
                status="exhaustive_search_incomplete", force_witness=True
            ),
        ):
            run = self._backend().run(
                input_dir=self.input_dir,
                analysis_spec=self.spec,
                output_dir=self.root / "incomplete-with-witness",
                method="exact",
            )
        self.assertEqual(run.status, "exhaustive_search_incomplete")
        self.assertTrue(run.witness_found)

    def test_existing_output_is_rejected_before_backend_launch(self) -> None:
        output = self.root / "existing"
        output.mkdir()
        (output / "result.json").write_text("stale", encoding="utf-8")
        with mock.patch("annotation_breakpoint.bridge.subprocess.run") as invoked:
            with self.assertRaises(StaleOutputError):
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=output,
                )
        invoked.assert_not_called()

    def test_concurrent_output_claim_is_rejected_before_analysis(self) -> None:
        output = self.root / "claimed-during-discovery"
        calls = 0

        def invoke(command, **kwargs):
            nonlocal calls
            calls += 1
            output.mkdir()
            return subprocess.CompletedProcess(command, 0, f"0.3.0\n{self.cli}\n", "")

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            with self.assertRaisesRegex(StaleOutputError, "another caller"):
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=output,
                )
        self.assertEqual(calls, 1)

    def test_nonzero_exit_preserves_backend_streams(self) -> None:
        calls = 0

        def invoke(command, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 1:
                return subprocess.CompletedProcess(command, 0, f"0.3.0\n{self.cli}\n", "")
            return subprocess.CompletedProcess(command, 7, "partial stdout", "backend failed")

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            with self.assertRaises(BackendProcessError) as caught:
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=self.root / "failed",
                )
        self.assertEqual(caught.exception.returncode, 7)
        self.assertEqual(caught.exception.stdout, "partial stdout")
        self.assertEqual(caught.exception.stderr, "backend failed")

    def test_timeout_preserves_partial_streams(self) -> None:
        calls = 0

        def invoke(command, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 1:
                return subprocess.CompletedProcess(command, 0, f"0.3.0\n{self.cli}\n", "")
            raise subprocess.TimeoutExpired(command, 0.01, output=b"partial", stderr=b"timed out")

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            with self.assertRaises(BackendTimeoutError) as caught:
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=self.root / "timed-out",
                    timeout=0.01,
                )
        self.assertEqual(caught.exception.stdout, "partial")
        self.assertEqual(caught.exception.stderr, "timed out")

    def test_malformed_json_never_returns_success(self) -> None:
        def invoke(command, **kwargs):
            if "-e" in command:
                return subprocess.CompletedProcess(command, 0, f"0.3.0\n{self.cli}\n", "")
            target = Path(command[command.index("--output-dir") + 1])
            _write_output(target, malformed_result=True)
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            with self.assertRaises(ArtifactError):
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=self.root / "malformed",
                )

    def test_identity_mismatch_never_returns_success(self) -> None:
        def mutate(result, verification):
            verification["result_identity"] = "another-result"

        with mock.patch(
            "annotation_breakpoint.bridge.subprocess.run",
            side_effect=self._successful_side_effect(mutate=mutate),
        ):
            with self.assertRaisesRegex(ArtifactError, "different result identity"):
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=self.root / "identity-mismatch",
                )

    def test_artifact_hash_mismatch_never_returns_success(self) -> None:
        def invoke(command, **kwargs):
            if "-e" in command:
                return subprocess.CompletedProcess(command, 0, f"0.3.0\n{self.cli}\n", "")
            target = Path(command[command.index("--output-dir") + 1])
            _write_output(target)
            (target / "report.md").write_text("changed", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, "", "")

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            with self.assertRaisesRegex(ArtifactError, "hash mismatch"):
                self._backend().run(
                    input_dir=self.input_dir,
                    analysis_spec=self.spec,
                    output_dir=self.root / "hash-mismatch",
                )

    def test_wrong_package_version_is_rejected(self) -> None:
        with mock.patch(
            "annotation_breakpoint.bridge.subprocess.run",
            return_value=subprocess.CompletedProcess([], 0, f"0.2.0\n{self.cli}\n", ""),
        ):
            with self.assertRaises(UnsupportedBackendError):
                self._backend().discover()

    def test_missing_rscript_has_useful_error(self) -> None:
        with mock.patch("annotation_breakpoint.bridge.shutil.which", return_value=None):
            with self.assertRaisesRegex(BackendDiscoveryError, "Rscript was not found"):
                AnnotationBreakpointBackend().discover()

    def test_exact_route_and_argument_validation(self) -> None:
        calls = []

        def invoke(command, **kwargs):
            calls.append(command)
            return self._successful_side_effect(status="exact_minimum")(command, **kwargs)

        with mock.patch("annotation_breakpoint.bridge.subprocess.run", side_effect=invoke):
            run = self._backend().run(
                input_dir=self.input_dir,
                analysis_spec=self.spec,
                output_dir=self.root / "exact",
                method="exact",
                max_evaluations=42,
            )
        self.assertEqual(run.status, "exact_minimum")
        self.assertIn("42", calls[1])
        self.assertNotIn("--max-fits", calls[1])
        with self.assertRaisesRegex(ValueError, "prefix_sizes"):
            self._backend().run(
                input_dir=self.input_dir,
                analysis_spec=self.spec,
                output_dir=self.root / "never",
                method="exact",
                prefix_sizes=(1, 2),
            )


if __name__ == "__main__":
    unittest.main()
