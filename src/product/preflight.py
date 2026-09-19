"""Headless preflight adapter for orchestration and future batch callers."""

from __future__ import annotations

from dataclasses import dataclass

from .capabilities import Capability, CapabilityReport, detect_capabilities
from .datasets import DataRegistry
from .status import ProductStatus


@dataclass(frozen=True, slots=True)
class AnalysisRequest:
    target: str
    tissue: str
    engine: str = "classic"


@dataclass(frozen=True, slots=True)
class PreflightResult:
    request: AnalysisRequest
    status: ProductStatus
    missing: tuple[Capability, ...]
    capabilities: CapabilityReport


ENGINE_REQUIREMENTS = {
    "classic": (Capability.HAS_PPI, Capability.HAS_TISSUE_CONTEXT, Capability.HAS_CANONICAL_MAPPING),
    "directed": (
        Capability.HAS_PPI, Capability.HAS_TISSUE_CONTEXT,
        Capability.HAS_CANONICAL_MAPPING, Capability.HAS_DIRECTION,
    ),
    "evidence": (Capability.HAS_PPI, Capability.HAS_EVIDENCE_CHANNELS),
}


def preflight(request: AnalysisRequest, registry: DataRegistry) -> PreflightResult:
    engine = request.engine.casefold()
    if engine not in ENGINE_REQUIREMENTS:
        raise ValueError(f"unknown engine: {request.engine}")
    report = detect_capabilities(registry)
    missing = tuple(
        capability for capability in ENGINE_REQUIREMENTS[engine]
        if not report.available(capability)
    )
    return PreflightResult(
        request=request,
        status=ProductStatus.DATA_UNAVAILABLE if missing else ProductStatus.READY,
        missing=missing,
        capabilities=report,
    )
