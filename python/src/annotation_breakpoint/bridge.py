"""Subprocess bridge to the validated Annotation Breakpoint R engine."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

EXPECTED_R_PACKAGE_VERSION = "0.3.0"
REQUIRED_WITNESS_STATUSES = frozenset({"exact_minimum", "verified_counterexample"})
OPTIONAL_WITNESS_STATUSES = frozenset({"exhaustive_search_incomplete"})
SUPPORTED_STATUSES = frozenset(
    {
        "exact_minimum",
        "verified_counterexample",
        "exhaustively_stable",
        "exhaustive_search_incomplete",
        "no_counterexample_found",
        "original_claim_not_met",
        "not_assessable",
    }
)
_DISCOVERY_EXPRESSION = (
    'cat(as.character(packageVersion("AnnotationBreakpoint")), "\\n", '
    'system.file("exec", "annotation-breakpoint", '
    'package="AnnotationBreakpoint"), "\\n", sep="")'
)


class BridgeError(RuntimeError):
    """Base class for Python bridge failures."""

    def as_dict(self) -> dict[str, Any]:
        return {"error_type": type(self).__name__, "message": str(self)}


class BackendDiscoveryError(BridgeError):
    """The required R executable or installed package could not be found."""


class UnsupportedBackendError(BridgeError):
    """The discovered R package is not the supported release."""


class StaleOutputError(BridgeError):
    """The requested output path already exists."""


class ArtifactError(BridgeError):
    """Fresh backend output is missing, malformed, or internally inconsistent."""


class BackendProcessError(BridgeError):
    """The R backend exited unsuccessfully."""

    def __init__(
        self,
        message: str,
        *,
        returncode: int,
        command: Sequence[str],
        stdout: str,
        stderr: str,
        output_dir: Path,
    ) -> None:
        super().__init__(message)
        self.returncode = returncode
        self.command = tuple(command)
        self.stdout = stdout
        self.stderr = stderr
        self.output_dir = output_dir

    def as_dict(self) -> dict[str, Any]:
        value = super().as_dict()
        value.update(
            {
                "returncode": self.returncode,
                "command": list(self.command),
                "stdout": self.stdout,
                "stderr": self.stderr,
                "output_dir": str(self.output_dir),
            }
        )
        return value


class BackendTimeoutError(BridgeError):
    """The R backend exceeded the Python-side wall-time limit."""

    def __init__(
        self,
        message: str,
        *,
        timeout: float,
        command: Sequence[str],
        stdout: str,
        stderr: str,
        output_dir: Path,
    ) -> None:
        super().__init__(message)
        self.timeout = timeout
        self.command = tuple(command)
        self.stdout = stdout
        self.stderr = stderr
        self.output_dir = output_dir

    def as_dict(self) -> dict[str, Any]:
        value = super().as_dict()
        value.update(
            {
                "timeout": self.timeout,
                "command": list(self.command),
                "stdout": self.stdout,
                "stderr": self.stderr,
                "output_dir": str(self.output_dir),
            }
        )
        return value


@dataclass(frozen=True)
class BackendInfo:
    """Discovered backend executables and installed R package version."""

    rscript: Path
    r_cli: Path
    package_version: str


@dataclass(frozen=True)
class RunResult:
    """Fresh structured output from one R backend invocation."""

    result: Mapping[str, Any]
    verification: Mapping[str, Any]
    manifest: Mapping[str, Any]
    artifacts: Mapping[str, Path]
    backend: BackendInfo
    command: tuple[str, ...]
    stdout: str
    stderr: str
    output_dir: Path

    @property
    def status(self) -> str:
        return str(self.result["status"])

    @property
    def verification_ok(self) -> bool:
        return self.verification.get("ok") is True

    @property
    def witness_found(self) -> bool:
        return (
            self.verification_ok
            and self.result.get("witness") is not None
            and self.result.get("best_verified_cost") is not None
        )

    def summary(self) -> dict[str, Any]:
        """Return a compact JSON-ready summary without altering raw documents."""

        return {
            "status": self.status,
            "witness_found": self.witness_found,
            "verification_ok": self.verification_ok,
            "best_verified_cost": self.result.get("best_verified_cost"),
            "exact_minimum_cost": self.result.get("exact_minimum_cost"),
            "result_id": self.result.get("result_id"),
            "witness_id": self.result.get("witness_id"),
            "input_hashes": self.result.get("input_hashes"),
            "backend_identity": self.result.get("backend_identity"),
            "runtime": (self.result.get("original_fit") or {}).get("runtime"),
            "artifacts": {key: str(value) for key, value in self.artifacts.items()},
            "output_dir": str(self.output_dir),
            "r_package_version": self.backend.package_version,
        }


def _text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _number(value: int | float, name: str) -> str:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric")
    if math.isnan(float(value)) or float(value) < 0:
        raise ValueError(f"{name} must be nonnegative")
    if math.isinf(float(value)):
        return "Inf"
    return str(value)


def _prefix_schedule(values: Sequence[int] | None) -> tuple[int, ...] | None:
    if values is None:
        return None
    converted: list[int] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError("prefix_sizes must contain positive integers")
        converted.append(value)
    if not converted or any(a >= b for a, b in zip(converted, converted[1:])):
        raise ValueError("prefix_sizes must be nonempty and strictly increasing")
    return tuple(converted)


class AnnotationBreakpointBackend:
    """Discover and invoke the installed AnnotationBreakpoint 0.3.0 R CLI."""

    def __init__(
        self,
        *,
        rscript: str | os.PathLike[str] | None = None,
        r_cli: str | os.PathLike[str] | None = None,
        env: Mapping[str, str] | None = None,
        discovery_timeout: float = 30.0,
    ) -> None:
        self._rscript = rscript
        self._r_cli = r_cli
        self._env = dict(env or {})
        self._discovery_timeout = discovery_timeout
        self._backend_info: BackendInfo | None = None

    def _environment(self) -> dict[str, str]:
        environment = os.environ.copy()
        environment.update(self._env)
        return environment

    def discover(self) -> BackendInfo:
        """Locate Rscript and the installed R CLI, enforcing package 0.3.0."""

        if self._backend_info is not None:
            return self._backend_info
        supplied = self._rscript or os.environ.get("ANNOTATION_BREAKPOINT_RSCRIPT")
        located = str(supplied) if supplied is not None else shutil.which("Rscript")
        if not located:
            raise BackendDiscoveryError(
                "Rscript was not found; install R or provide rscript/"
                "ANNOTATION_BREAKPOINT_RSCRIPT"
            )
        rscript = Path(located).expanduser()
        if not rscript.is_file():
            raise BackendDiscoveryError(f"Rscript does not exist: {rscript}")
        command = [str(rscript), "--vanilla", "-e", _DISCOVERY_EXPRESSION]
        try:
            completed = subprocess.run(
                command,
                shell=False,
                capture_output=True,
                text=True,
                timeout=self._discovery_timeout,
                env=self._environment(),
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise BackendDiscoveryError(
                f"could not query the AnnotationBreakpoint R package: {error}"
            ) from error
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise BackendDiscoveryError(
                "could not query the AnnotationBreakpoint R package"
                + (f": {detail}" if detail else "")
            )
        lines = completed.stdout.splitlines()
        if len(lines) < 2 or not lines[0].strip():
            raise BackendDiscoveryError("R package discovery returned malformed output")
        package_version = lines[0].strip()
        if package_version != EXPECTED_R_PACKAGE_VERSION:
            raise UnsupportedBackendError(
                "annotation-breakpoint-python 0.1.0 requires the "
                f"AnnotationBreakpoint R package {EXPECTED_R_PACKAGE_VERSION}; "
                f"found {package_version}"
            )
        supplied_cli = self._r_cli or os.environ.get("ANNOTATION_BREAKPOINT_R_CLI")
        cli_value = str(supplied_cli) if supplied_cli is not None else lines[1].strip()
        if not cli_value:
            raise BackendDiscoveryError(
                "the installed AnnotationBreakpoint package has no R CLI"
            )
        r_cli = Path(cli_value).expanduser()
        if not r_cli.is_file():
            raise BackendDiscoveryError(f"Annotation Breakpoint R CLI does not exist: {r_cli}")
        self._backend_info = BackendInfo(
            rscript=rscript.resolve(),
            r_cli=r_cli.resolve(),
            package_version=package_version,
        )
        return self._backend_info

    def run(
        self,
        *,
        input_dir: str | os.PathLike[str],
        output_dir: str | os.PathLike[str],
        analysis_spec: str | os.PathLike[str] | None = None,
        layer_a: bool = False,
        method: str = "proxy_prefix",
        allow_version_mismatch: bool = False,
        prefix_sizes: Sequence[int] | None = None,
        max_evaluations: int | float = math.inf,
        max_fits: int | float = 1024,
        max_proposals: int | float = math.inf,
        max_seconds: int | float = math.inf,
        seed: int = 82901,
        timeout: float | None = None,
    ) -> RunResult:
        """Run one fresh analysis and validate all emitted identities.

        A successful process and verification do not imply that a witness was
        found. Inspect ``status`` or ``witness_found`` on the returned object.
        """

        if method not in {"proxy_prefix", "exact"}:
            raise ValueError("method must be 'proxy_prefix' or 'exact'")
        schedule = _prefix_schedule(prefix_sizes)
        if method == "exact" and schedule is not None:
            raise ValueError("prefix_sizes is available only with proxy_prefix")
        if bool(analysis_spec) == bool(layer_a):
            raise ValueError("provide exactly one of analysis_spec or layer_a=True")
        if allow_version_mismatch and (layer_a or analysis_spec is None):
            raise ValueError(
                "allow_version_mismatch requires an authoritative analysis_spec"
            )
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise ValueError("seed must be an integer")
        if timeout is not None and timeout <= 0:
            raise ValueError("timeout must be positive")

        source = Path(input_dir).expanduser().resolve()
        if not source.is_dir():
            raise ValueError(f"input_dir is not a directory: {source}")
        spec = None if analysis_spec is None else Path(analysis_spec).expanduser().resolve()
        if spec is not None and not spec.is_file():
            raise ValueError(f"analysis_spec is not a file: {spec}")
        target = Path(output_dir).expanduser().resolve()
        if target.exists():
            raise StaleOutputError(
                f"output_dir already exists; choose a fresh path: {target}"
            )
        target.parent.mkdir(parents=True, exist_ok=True)

        backend = self.discover()
        command = [
            str(backend.rscript),
            "--vanilla",
            str(backend.r_cli),
            "--input-dir",
            str(source),
        ]
        if spec is not None:
            command.extend(["--analysis-spec", str(spec)])
        else:
            command.append("--layer-a")
        if allow_version_mismatch:
            command.append("--allow-version-mismatch")
        command.extend(["--output-dir", str(target), "--method", method])
        if method == "exact":
            command.extend(["--max-evaluations", _number(max_evaluations, "max_evaluations")])
        else:
            command.extend(
                [
                    "--max-fits",
                    _number(max_fits, "max_fits"),
                    "--max-proposals",
                    _number(max_proposals, "max_proposals"),
                    "--max-seconds",
                    _number(max_seconds, "max_seconds"),
                    "--seed",
                    str(seed),
                ]
            )
            if schedule is not None:
                command.extend(["--prefix-sizes", ",".join(map(str, schedule))])
        try:
            target.mkdir(exist_ok=False)
        except FileExistsError as error:
            raise StaleOutputError(
                f"output_dir was claimed by another caller: {target}"
            ) from error
        try:
            completed = subprocess.run(
                command,
                shell=False,
                capture_output=True,
                text=True,
                timeout=timeout,
                env=self._environment(),
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            raise BackendTimeoutError(
                f"R backend exceeded the {timeout}-second Python timeout",
                timeout=float(timeout),
                command=command,
                stdout=_text(error.stdout),
                stderr=_text(error.stderr),
                output_dir=target,
            ) from error
        except OSError as error:
            raise BackendProcessError(
                f"could not start the R backend: {error}",
                returncode=-1,
                command=command,
                stdout="",
                stderr=str(error),
                output_dir=target,
            ) from error
        if completed.returncode != 0:
            raise BackendProcessError(
                f"R backend exited with status {completed.returncode}",
                returncode=completed.returncode,
                command=command,
                stdout=completed.stdout,
                stderr=completed.stderr,
                output_dir=target,
            )
        return self._load_fresh_output(
            target=target,
            backend=backend,
            command=command,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )

    def _load_fresh_output(
        self,
        *,
        target: Path,
        backend: BackendInfo,
        command: Sequence[str],
        stdout: str,
        stderr: str,
    ) -> RunResult:
        manifest_path = target / "run_manifest.json"
        if not manifest_path.is_file():
            raise ArtifactError(f"fresh run did not produce {manifest_path}")
        manifest = _read_json(manifest_path, "run manifest")
        records = manifest.get("artifacts")
        if not isinstance(records, list):
            raise ArtifactError("run manifest has no artifact list")
        artifacts: dict[str, Path] = {"manifest": manifest_path}
        seen: set[str] = set()
        for record in records:
            if not isinstance(record, dict):
                raise ArtifactError("run manifest contains a malformed artifact record")
            name, filename, expected = (
                record.get("name"),
                record.get("file"),
                record.get("sha256"),
            )
            if not all(isinstance(item, str) and item for item in (name, filename, expected)):
                raise ArtifactError("run manifest contains an incomplete artifact record")
            if name in seen or Path(filename).name != filename:
                raise ArtifactError("run manifest contains duplicate or unsafe artifact paths")
            seen.add(name)
            path = target / filename
            if not path.is_file():
                raise ArtifactError(f"fresh run artifact is missing: {filename}")
            observed = _sha256(path)
            if observed != expected:
                raise ArtifactError(
                    f"artifact hash mismatch for {filename}: expected {expected}, found {observed}"
                )
            artifacts[name] = path
        for required in ("result_json", "verification_json", "claim_results", "report"):
            if required not in artifacts:
                raise ArtifactError(f"run manifest is missing required artifact {required}")
        result = _read_json(artifacts["result_json"], "result")
        verification = _read_json(artifacts["verification_json"], "verification")
        self._validate_documents(result, verification, manifest, backend)
        return RunResult(
            result=result,
            verification=verification,
            manifest=manifest,
            artifacts=artifacts,
            backend=backend,
            command=tuple(command),
            stdout=stdout,
            stderr=stderr,
            output_dir=target,
        )

    @staticmethod
    def _validate_documents(
        result: Mapping[str, Any],
        verification: Mapping[str, Any],
        manifest: Mapping[str, Any],
        backend: BackendInfo,
    ) -> None:
        status = result.get("status")
        if status not in SUPPORTED_STATUSES:
            raise ArtifactError(f"unsupported or missing result status: {status!r}")
        if manifest.get("status") != status or verification.get("recomputed_status") != status:
            raise ArtifactError("result, verification, and manifest statuses do not agree")
        if verification.get("ok") is not True or manifest.get("verification_ok") is not True:
            raise ArtifactError("backend did not emit a successful independent verification")
        result_id = result.get("result_id")
        if not isinstance(result_id, str) or not result_id:
            raise ArtifactError("result identity is missing")
        if verification.get("result_identity") != result_id:
            raise ArtifactError("verification belongs to a different result identity")
        hashes = result.get("input_hashes")
        if not isinstance(hashes, dict) or not hashes:
            raise ArtifactError("result input hashes are missing")
        if verification.get("input_hashes") != hashes or manifest.get("input_hashes") != hashes:
            raise ArtifactError("result, verification, and manifest input hashes do not agree")
        result_witness = result.get("witness_id")
        verification_witness = verification.get("witness_id")
        if result_witness != verification_witness:
            raise ArtifactError("result and verification witness identities do not agree")
        for match in verification.get("hash_matches") or []:
            if not isinstance(match, dict) or match.get("matches") is not True:
                raise ArtifactError("verification contains a failed input-hash comparison")
        run_metadata = manifest.get("run_metadata")
        if not isinstance(run_metadata, dict):
            raise ArtifactError("run metadata are missing")
        if str(run_metadata.get("package_version")) != backend.package_version:
            raise ArtifactError("run metadata report a different R package version")
        has_witness = result.get("witness") is not None
        if status in REQUIRED_WITNESS_STATUSES and not has_witness:
            raise ArtifactError(f"witness status {status!r} has no witness record")
        if status not in REQUIRED_WITNESS_STATUSES | OPTIONAL_WITNESS_STATUSES and has_witness:
            raise ArtifactError(f"non-witness status {status!r} unexpectedly has a witness")


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ArtifactError(f"could not read {label} JSON at {path}: {error}") from error
    if not isinstance(value, dict):
        raise ArtifactError(f"{label} JSON must be an object")
    return value


def run_annotation_breakpoint(**kwargs: Any) -> RunResult:
    """Convenience function using default backend discovery."""

    backend_keys = {"rscript", "r_cli", "env", "discovery_timeout"}
    backend_options = {key: kwargs.pop(key) for key in list(kwargs) if key in backend_keys}
    return AnnotationBreakpointBackend(**backend_options).run(**kwargs)
