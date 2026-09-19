"""Run provenance contract shared by headless, UI and export adapters."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Mapping

from .capabilities import CapabilityReport
from .versions import ENGINE_VERSIONS, EXPORT_SCHEMA_VERSION, PREPROCESSING_VERSIONS, SOPHIARK_VERSION


def build_run_manifest(
    *, target: str, tissue: str, dataset_versions: Mapping[str, str],
    dataset_fingerprints: Mapping[str, str], capabilities: CapabilityReport,
    engine_config_fingerprints: Mapping[str, str] | None = None,
    timestamp: datetime | None = None,
) -> dict[str, Any]:
    now = timestamp or datetime.now(timezone.utc)
    return {
        "sophiark_version": SOPHIARK_VERSION,
        "target": target,
        "tissue": tissue,
        "engine_versions": dict(ENGINE_VERSIONS),
        "engine_config_fingerprints": dict(engine_config_fingerprints or {}),
        "datasets": {
            name: {"version": dataset_versions.get(name), "fingerprint": fingerprint}
            for name, fingerprint in dataset_fingerprints.items()
        },
        "preprocessing_versions": dict(PREPROCESSING_VERSIONS),
        "timestamp_utc": now.astimezone(timezone.utc).isoformat(),
        "availability": {key.value: value.value for key, value in capabilities.statuses.items()},
        "missing_optional_layers": sorted(
            name for name, health in capabilities.datasets.items()
            if not health.required and health.status is not ProductStatus.READY
        ),
        "export_schema_version": EXPORT_SCHEMA_VERSION,
    }


from .status import ProductStatus  # kept last to keep the manifest imports visually grouped
