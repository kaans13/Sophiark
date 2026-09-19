"""Offline orchestration boundary for the read-only Research Context stack.

The normal constructor is dependency-injected and performs no I/O.  The
explicit :meth:`ResearchContextService.from_project_root` constructor performs
read-only local index discovery/loading and is therefore named separately.
Neither path enables an external provider or touches the scientific core.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime
from enum import Enum
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time
from typing import Any

import pandas as pd

from .config import ResearchConfig
from .entity_resolution import EntityResolver
from .local_context import LocalResearchContext, build_local_research_context
from .models import (
    EvidenceAssertion,
    ObservationType,
    ProvenanceRecord,
    SimulationResultSnapshot,
    immutable_mapping,
)
from .observations import (
    DetectedObservation,
    ObservationDetectionResult,
    detect_observations,
)
from .providers import LocalEvidenceAdapters
from .relationships import RelationshipBuildResult, build_target_response_relationships
from .semantic_adapter import SemanticSimulationResult, SemanticResultView, adapt_snapshot
from .snapshot import RESULT_SCHEMA_VERSION, build_snapshot


_EXTERNAL_CONFIG_FIELDS = (
    "external_context_enabled",
    "provider_uniprot_enabled",
    "provider_interpro_enabled",
    "provider_quickgo_enabled",
    "provider_reactome_enabled",
    "provider_europepmc_enabled",
)


class ServiceInitializationMode(str, Enum):
    INJECTED_NO_IO = "INJECTED_NO_IO"
    PROJECT_ROOT_READ_ONLY_IO = "PROJECT_ROOT_READ_ONLY_IO"


@dataclass(frozen=True, slots=True)
class ResearchTimings:
    snapshot_build_ms: float = 0.0
    semantic_adapter_ms: float = 0.0
    local_context_ms: float = 0.0
    relationship_ms: float = 0.0
    observation_ms: float = 0.0
    signaling_ms: float = 0.0
    total_ms: float = 0.0

    def __post_init__(self) -> None:
        values = tuple(float(value) for value in asdict(self).values())
        if any(not math.isfinite(value) or value < 0 for value in values):
            raise ValueError("research timing values must be finite and non-negative")
        component_total = sum(values[:-1])
        if self.total_ms + 1e-9 < component_total:
            raise ValueError("total_ms cannot be shorter than its component timings")

    def as_mapping(self) -> Mapping[str, float]:
        return immutable_mapping(asdict(self))


@dataclass(frozen=True, slots=True)
class OfflineResearchBundle:
    """One immutable, snapshot-bound offline Research Context result."""

    bundle_id: str
    snapshot: SimulationResultSnapshot
    semantic_result: SemanticSimulationResult
    local_context: LocalResearchContext
    relationships: RelationshipBuildResult
    observation_result: ObservationDetectionResult
    timings: ResearchTimings
    config: ResearchConfig
    notices: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    signaling_result: Any = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "notices", tuple(map(str, self.notices)))
        object.__setattr__(self, "errors", tuple(map(str, self.errors)))
        snapshot_id = self.snapshot.snapshot_id
        simulation_id = self.snapshot.simulation_id
        components = (
            ("semantic_result", self.semantic_result.snapshot_id, self.semantic_result.simulation_id),
            ("local_context", self.local_context.snapshot_id, self.local_context.simulation_id),
            ("relationships", self.relationships.snapshot_id, self.relationships.simulation_id),
            ("observation_result", self.observation_result.snapshot_id, self.observation_result.simulation_id),
        )
        for name, component_snapshot_id, component_simulation_id in components:
            if component_snapshot_id != snapshot_id or component_simulation_id != simulation_id:
                raise ValueError(f"{name} is stale or belongs to a different snapshot")
        expected_prefix = f"offline:{snapshot_id[:16]}:"
        if not self.bundle_id.startswith(expected_prefix):
            raise ValueError("bundle_id is not bound to the supplied snapshot")
        _validate_offline_config(self.config)

    @property
    def snapshot_id(self) -> str:
        return self.snapshot.snapshot_id

    @property
    def simulation_id(self) -> str:
        return self.snapshot.simulation_id

    @property
    def healthy(self) -> bool:
        return not self.errors

    def observations_of_type(
        self, observation_type: ObservationType | str,
    ) -> tuple[DetectedObservation, ...]:
        return self.observation_result.of_type(observation_type)

    @property
    def family_observations(self) -> tuple[DetectedObservation, ...]:
        return self.observations_of_type(ObservationType.FAMILY_COOCCURRENCE)

    @property
    def functional_observations(self) -> tuple[DetectedObservation, ...]:
        return self.observations_of_type(ObservationType.FUNCTIONAL_COOCCURRENCE)

    @property
    def relationship_observations(self) -> tuple[DetectedObservation, ...]:
        return self.observations_of_type(ObservationType.TARGET_RESPONSE_RELATIONSHIP)

    def overview(self) -> Mapping[str, Any]:
        """Return a lightweight immutable projection; no new interpretation."""

        return immutable_mapping({
            "bundle_id": self.bundle_id,
            "snapshot_id": self.snapshot_id,
            "simulation_id": self.simulation_id,
            "species": self.snapshot.species,
            "taxon_id": self.snapshot.taxon_id,
            "tissue": self.snapshot.tissue,
            "selection_mode": self.snapshot.selection_mode,
            "target_count": len(self.semantic_result.target_ids),
            "resolved_entity_count": len(self.local_context.entities),
            "unresolved_entity_count": len(self.local_context.unresolved_entities),
            "relationship_count": len(self.relationships.assertions),
            "observation_count": len(self.observation_result.observations),
            "family_observation_count": len(self.family_observations),
            "functional_observation_count": len(self.functional_observations),
            "signaling_status": getattr(getattr(self, "signaling_result", None), "status", None),
            "healthy": self.healthy,
            "notice_count": len(self.notices),
            "error_count": len(self.errors),
            "timings_ms": self.timings.as_mapping(),
        })


def _validate_offline_config(config: ResearchConfig) -> None:
    if not isinstance(config, ResearchConfig):
        raise TypeError("config must be a ResearchConfig")
    if not config.research_context_enabled:
        raise ValueError("research_context_enabled must be true for explicit offline research")
    enabled_external = tuple(
        field_name for field_name in _EXTERNAL_CONFIG_FIELDS
        if bool(getattr(config, field_name))
    )
    if enabled_external:
        raise ValueError(
            "offline ResearchContextService requires all external flags to remain false: "
            + ", ".join(enabled_external)
        )


def _deduplicated_text(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values if str(value).strip()))


def _elapsed_ms(started: float) -> float:
    return max(0.0, (time.perf_counter() - started) * 1000.0)


def _empty_semantic_result(
    snapshot: SimulationResultSnapshot, message: str,
) -> SemanticResultView:
    return SemanticResultView(
        snapshot_id=snapshot.snapshot_id,
        simulation_id=snapshot.simulation_id,
        scope=snapshot.scope,
        target_ids=snapshot.targets,
        target_records=(),
        positive_redistribution=(),
        network_losses=(),
        fdr_supported=(),
        go_terms=(),
        kegg_pathways=(),
        unclassified_enrichment=(),
        notices=(message,),
    )


def _empty_local_context(
    snapshot: SimulationResultSnapshot,
    semantic_result: SemanticSimulationResult,
    message: str,
) -> LocalResearchContext:
    return LocalResearchContext(
        snapshot_id=snapshot.snapshot_id,
        simulation_id=snapshot.simulation_id,
        semantic_result=semantic_result,
        entities=(),
        unresolved_entities=(),
        roles={},
        annotations={},
        functional_context=semantic_result.functional_rows,
        source_statuses=(),
        provenance=snapshot.provenance,
        notices=(message,),
        errors=(message,),
    )


def _empty_relationships(
    snapshot: SimulationResultSnapshot, message: str,
) -> RelationshipBuildResult:
    return RelationshipBuildResult(
        simulation_id=snapshot.simulation_id,
        snapshot_id=snapshot.snapshot_id,
        assertions=(),
        notices=(message,),
    )


def _empty_observations(
    snapshot: SimulationResultSnapshot, message: str,
) -> ObservationDetectionResult:
    return ObservationDetectionResult(
        simulation_id=snapshot.simulation_id,
        snapshot_id=snapshot.snapshot_id,
        observations=(),
        notices=(message,),
    )


def _bundle_id(
    snapshot: SimulationResultSnapshot,
    relationships: RelationshipBuildResult,
    observations: ObservationDetectionResult,
) -> str:
    payload = {
        "snapshot_id": snapshot.snapshot_id,
        "relationships": [item.assertion_id for item in relationships.assertions],
        "observations": [item.observation_id for item in observations.observations],
    }
    encoded = json.dumps(
        payload, ensure_ascii=False, allow_nan=False,
        sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:16]
    return f"offline:{snapshot.snapshot_id[:16]}:{digest}"


@dataclass(frozen=True, slots=True)
class ResearchContextService:
    """Reusable offline orchestrator with explicitly injected local dependencies."""

    resolver: EntityResolver
    providers: LocalEvidenceAdapters
    config: ResearchConfig
    initialization_mode: ServiceInitializationMode = ServiceInitializationMode.INJECTED_NO_IO
    project_root: str | None = None
    signaling_graphs: Mapping[str, Any] = field(default_factory=dict)
    signaling_load_errors: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.resolver, EntityResolver):
            raise TypeError("resolver must be an EntityResolver")
        if not isinstance(self.providers, LocalEvidenceAdapters):
            raise TypeError("providers must be LocalEvidenceAdapters")
        _validate_offline_config(self.config)
        object.__setattr__(self, "initialization_mode", ServiceInitializationMode(self.initialization_mode))
        object.__setattr__(self, "signaling_graphs", immutable_mapping(self.signaling_graphs))
        object.__setattr__(self, "signaling_load_errors", immutable_mapping(self.signaling_load_errors))

    @classmethod
    def from_project_root(
        cls,
        project_root: str | Path,
        *,
        config: ResearchConfig,
        include_human: bool = True,
        include_mouse: bool = True,
    ) -> "ResearchContextService":
        """Explicit read-only local-I/O constructor for app-level caching."""

        _validate_offline_config(config)
        root = Path(project_root).resolve()
        resolver = EntityResolver.from_project_data(
            root,
            include_human=include_human,
            include_mouse=include_mouse,
            external_fallback=None,
        )
        providers = LocalEvidenceAdapters.from_project_root(root)
        signaling_graphs: dict[str, Any] = {}
        signaling_errors: dict[str, str] = {}
        if config.directed_signaling_enabled:
            from src.signaling import (
                DirectedSignalingGraph, load_omnipath_signaling_dataset, load_signaling_cache,
            )
            species_flags = ((9606, include_human), (10090, include_mouse))
            for taxon_id, enabled in species_flags:
                if not enabled:
                    continue
                path = root / "data" / "regulatory" / "omnipath" / f"omnipath_{taxon_id}.tsv"
                cache_path = root / "data" / "processed" / f"omnipath_signaling_{taxon_id}.sqlite"
                try:
                    try:
                        dataset = load_signaling_cache(cache_path, source_snapshot=path)
                    except (FileNotFoundError, ValueError, OSError, sqlite3.DatabaseError):
                        dataset = load_omnipath_signaling_dataset(path, resolver=resolver, taxon_id=taxon_id)
                    signaling_graphs[str(taxon_id)] = DirectedSignalingGraph(
                        dataset
                    )
                except Exception as exc:
                    signaling_errors[str(taxon_id)] = f"{type(exc).__name__}: {exc}"
        return cls(
            resolver=resolver,
            providers=providers,
            config=config,
            initialization_mode=ServiceInitializationMode.PROJECT_ROOT_READ_ONLY_IO,
            project_root=str(root),
            signaling_graphs=signaling_graphs,
            signaling_load_errors=signaling_errors,
        )

    def build_offline_bundle(
        self,
        snapshot: SimulationResultSnapshot,
        *,
        explicit_relationships: Iterable[EvidenceAssertion] = (),
        snapshot_build_ms: float = 0.0,
    ) -> OfflineResearchBundle:
        """Run the offline pipeline; every stage remains snapshot-bound."""

        if not isinstance(snapshot, SimulationResultSnapshot):
            raise TypeError("snapshot must be a SimulationResultSnapshot")
        _validate_offline_config(self.config)
        if not math.isfinite(float(snapshot_build_ms)) or float(snapshot_build_ms) < 0:
            raise ValueError("snapshot_build_ms must be finite and non-negative")
        explicit = tuple(explicit_relationships)
        errors: list[str] = []

        started = time.perf_counter()
        try:
            semantic_result = adapt_snapshot(snapshot)
        except Exception as exc:
            message = f"semantic_adapter: {type(exc).__name__}: {exc}"
            errors.append(message)
            semantic_result = _empty_semantic_result(snapshot, message)
        semantic_ms = _elapsed_ms(started)

        started = time.perf_counter()
        try:
            local_context = build_local_research_context(
                snapshot,
                resolver=self.resolver,
                providers=self.providers,
                semantic_result=semantic_result,
            )
        except Exception as exc:
            message = f"local_context: {type(exc).__name__}: {exc}"
            errors.append(message)
            local_context = _empty_local_context(snapshot, semantic_result, message)
        local_ms = _elapsed_ms(started)
        errors.extend(local_context.errors)

        started = time.perf_counter()
        try:
            relationships = build_target_response_relationships(
                local_context,
                explicit_evidence=explicit,
            )
        except Exception as exc:
            message = f"relationships: {type(exc).__name__}: {exc}"
            errors.append(message)
            relationships = _empty_relationships(snapshot, message)
        relationship_ms = _elapsed_ms(started)

        started = time.perf_counter()
        try:
            observation_result = detect_observations(
                local_context,
                relationships,
                config=self.config,
                created_at=snapshot.created_at,
            )
        except Exception as exc:
            message = f"observations: {type(exc).__name__}: {exc}"
            errors.append(message)
            observation_result = _empty_observations(snapshot, message)
        observation_ms = _elapsed_ms(started)

        started = time.perf_counter()
        signaling_result = None
        if self.config.directed_signaling_enabled:
            from src.signaling.analysis import safe_analyze_signaling_context
            graph = self.signaling_graphs.get(str(snapshot.taxon_id))
            report_frame = snapshot.report.to_frame()
            tissue_support = None
            if str(snapshot.tissue).strip().casefold() not in {"", "none"}:
                tissue_support = tuple(report_frame.get("gene", ())) + tuple(snapshot.targets)
            signaling_result = safe_analyze_signaling_context(
                graph, report_frame, targets=snapshot.targets,
                max_depth=self.config.directed_signaling_max_depth,
                tissue_supported_entities=tissue_support,
            )
            load_error = self.signaling_load_errors.get(str(snapshot.taxon_id))
            if load_error and signaling_result.error is None:
                signaling_result = type(signaling_result)(
                    status=signaling_result.status, rows=signaling_result.rows,
                    convergence=signaling_result.convergence, metadata=signaling_result.metadata,
                    notices=(*signaling_result.notices, f"OmniPath load error: {load_error}"),
                    error=load_error,
                )
        signaling_ms = _elapsed_ms(started)

        snapshot_ms = float(snapshot_build_ms)
        timings = ResearchTimings(
            snapshot_build_ms=snapshot_ms,
            semantic_adapter_ms=semantic_ms,
            local_context_ms=local_ms,
            relationship_ms=relationship_ms,
            observation_ms=observation_ms,
            signaling_ms=signaling_ms,
            total_ms=(snapshot_ms + semantic_ms + local_ms + relationship_ms + observation_ms + signaling_ms),
        )
        notices = _deduplicated_text((
            *semantic_result.notices,
            *local_context.notices,
            *relationships.notices,
            *observation_result.notices,
            *(tuple(getattr(signaling_result, "notices", ())) if signaling_result is not None else ()),
        ))
        return OfflineResearchBundle(
            bundle_id=_bundle_id(snapshot, relationships, observation_result),
            snapshot=snapshot,
            semantic_result=semantic_result,
            local_context=local_context,
            relationships=relationships,
            observation_result=observation_result,
            timings=timings,
            config=self.config,
            notices=notices,
            errors=_deduplicated_text(errors),
            signaling_result=signaling_result,
        )

    # Concise alias used by non-UI callers.
    build_bundle = build_offline_bundle

    def build_from_results(
        self,
        report: pd.DataFrame,
        signed_presentation_result: pd.DataFrame | None,
        enrichment: pd.DataFrame | None,
        *,
        species: str,
        taxon_id: int,
        tissue: str,
        targets: Sequence[str],
        attenuation: float | None,
        selection_mode: str,
        test_limit: int | None,
        tested_count: int | None,
        returned_count: int,
        top_n: int | None = None,
        threshold_parameters: Mapping[str, Any] | None = None,
        result_schema_version: int = RESULT_SCHEMA_VERSION,
        provenance: Iterable[ProvenanceRecord] | None = None,
        simulation_id: str | None = None,
        created_at: datetime | str | None = None,
        explicit_relationships: Iterable[EvidenceAssertion] = (),
    ) -> OfflineResearchBundle:
        started = time.perf_counter()
        snapshot = build_research_snapshot(
            report,
            signed_presentation_result,
            enrichment,
            species=species,
            taxon_id=taxon_id,
            tissue=tissue,
            targets=targets,
            attenuation=attenuation,
            selection_mode=selection_mode,
            test_limit=test_limit,
            tested_count=tested_count,
            returned_count=returned_count,
            top_n=top_n,
            threshold_parameters=threshold_parameters,
            result_schema_version=result_schema_version,
            provenance=provenance,
            simulation_id=simulation_id,
            created_at=created_at,
        )
        snapshot_ms = _elapsed_ms(started)
        return self.build_offline_bundle(
            snapshot,
            explicit_relationships=explicit_relationships,
            snapshot_build_ms=snapshot_ms,
        )


def build_research_snapshot(
    report: pd.DataFrame,
    signed_presentation_result: pd.DataFrame | None,
    enrichment: pd.DataFrame | None,
    *,
    species: str,
    taxon_id: int,
    tissue: str,
    targets: Sequence[str],
    attenuation: float | None,
    selection_mode: str,
    test_limit: int | None,
    tested_count: int | None,
    returned_count: int,
    top_n: int | None = None,
    threshold_parameters: Mapping[str, Any] | None = None,
    result_schema_version: int = RESULT_SCHEMA_VERSION,
    provenance: Iterable[ProvenanceRecord] | None = None,
    simulation_id: str | None = None,
    created_at: datetime | str | None = None,
) -> SimulationResultSnapshot:
    """Freeze exact, already-bounded results without filtering or re-ranking."""

    return build_snapshot(
        report,
        signed_presentation_result,
        enrichment,
        species=species,
        taxon_id=taxon_id,
        tissue=tissue,
        targets=targets,
        attenuation=attenuation,
        selection_mode=selection_mode,
        test_limit=test_limit,
        tested_count=tested_count,
        returned_count=returned_count,
        top_n=top_n,
        threshold_parameters=threshold_parameters,
        result_schema_version=result_schema_version,
        provenance=provenance,
        simulation_id=simulation_id,
        created_at=created_at,
    )


def build_offline_research_bundle(
    service: ResearchContextService,
    snapshot: SimulationResultSnapshot,
    *,
    explicit_relationships: Iterable[EvidenceAssertion] = (),
) -> OfflineResearchBundle:
    """Pure facade-friendly delegation for an already-created service."""

    if not isinstance(service, ResearchContextService):
        raise TypeError("service must be a ResearchContextService")
    return service.build_offline_bundle(
        snapshot,
        explicit_relationships=explicit_relationships,
    )


__all__ = [
    "OfflineResearchBundle",
    "ResearchContextService",
    "ResearchTimings",
    "ServiceInitializationMode",
    "build_offline_research_bundle",
    "build_research_snapshot",
]
