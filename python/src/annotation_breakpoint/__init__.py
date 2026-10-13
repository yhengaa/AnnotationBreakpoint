"""Python interface to the installed Annotation Breakpoint R engine."""

from .bridge import (
    AnnotationBreakpointBackend,
    ArtifactError,
    BackendDiscoveryError,
    BackendInfo,
    BackendProcessError,
    BackendTimeoutError,
    BridgeError,
    RunResult,
    StaleOutputError,
    UnsupportedBackendError,
    run_annotation_breakpoint,
)

__all__ = [
    "AnnotationBreakpointBackend",
    "ArtifactError",
    "BackendDiscoveryError",
    "BackendInfo",
    "BackendProcessError",
    "BackendTimeoutError",
    "BridgeError",
    "RunResult",
    "StaleOutputError",
    "UnsupportedBackendError",
    "run_annotation_breakpoint",
]

__version__ = "0.1.0"
