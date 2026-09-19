"""Read-only Research Context package.

Public objects are imported lazily so ``import src.research`` performs no
dataset load, database initialization, Streamlit work, or external request.
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any


_LAZY_EXPORTS = {
    "AliasType": ("src.research.models", "AliasType"),
    "Entity": ("src.research.models", "Entity"),
    "EntityAlias": ("src.research.models", "EntityAlias"),
    "EntityResolution": ("src.research.models", "EntityResolution"),
    "EntityType": ("src.research.models", "EntityType"),
    "ResolutionStatus": ("src.research.models", "ResolutionStatus"),
    "CacheTTLs": ("src.research.config", "CacheTTLs"),
    "DEFAULT_RESEARCH_CONFIG": ("src.research.config", "DEFAULT_RESEARCH_CONFIG"),
    "PresentationLimits": ("src.research.config", "PresentationLimits"),
    "ProviderTimeout": ("src.research.config", "ProviderTimeout"),
    "ProviderTimeouts": ("src.research.config", "ProviderTimeouts"),
    "ResearchConfig": ("src.research.config", "ResearchConfig"),
    "RetrievalBudget": ("src.research.config", "RetrievalBudget"),
    "EvidenceProvider": ("src.research.provider_runtime", "EvidenceProvider"),
    "ProviderCapability": ("src.research.provider_runtime", "ProviderCapability"),
    "ProviderExecutionPolicy": ("src.research.provider_runtime", "ProviderExecutionPolicy"),
    "ProviderExecutor": ("src.research.provider_runtime", "ProviderExecutor"),
    "ProviderRequest": ("src.research.provider_runtime", "ProviderRequest"),
    "ProviderResult": ("src.research.provider_runtime", "ProviderResult"),
    "ProviderStatus": ("src.research.provider_runtime", "ProviderStatus"),
    "SessionCircuitBreaker": ("src.research.provider_runtime", "SessionCircuitBreaker"),
    "CacheEvidenceClass": ("src.research.provider_cache", "CacheEvidenceClass"),
    "CacheFirstProviderClient": ("src.research.provider_cache", "CacheFirstProviderClient"),
    "ProviderCache": ("src.research.provider_cache", "ProviderCache"),
    "ProviderCacheEntry": ("src.research.provider_cache", "ProviderCacheEntry"),
    "ProviderCacheKey": ("src.research.provider_cache", "ProviderCacheKey"),
    "ProviderCacheStats": ("src.research.provider_cache", "ProviderCacheStats"),
    "ProviderTTLPolicy": ("src.research.provider_cache", "ProviderTTLPolicy"),
    "ExternalBiologyLoad": ("src.research.external_biology", "ExternalBiologyLoad"),
    "ExternalBiologyLoader": ("src.research.external_biology", "ExternalBiologyLoader"),
    "ExternalBiologyPlan": ("src.research.external_biology", "ExternalBiologyPlan"),
    "plan_external_biology": ("src.research.external_biology", "plan_external_biology"),
    "LiteratureLoad": ("src.research.literature", "LiteratureLoad"),
    "LiteratureLoader": ("src.research.literature", "LiteratureLoader"),
    "LiteraturePlan": ("src.research.literature", "LiteraturePlan"),
    "LiteratureQuery": ("src.research.literature", "LiteratureQuery"),
    "LiteratureQueryBuilder": ("src.research.literature", "LiteratureQueryBuilder"),
    "LiteratureRepository": ("src.research.literature", "LiteratureRepository"),
    "LiteratureSearchHit": ("src.research.literature", "LiteratureSearchHit"),
    "LiteratureSearchResult": ("src.research.literature", "LiteratureSearchResult"),
    "Publication": ("src.research.literature", "Publication"),
    "PublicationMatch": ("src.research.literature", "PublicationMatch"),
    "PublicationStoreResult": ("src.research.literature", "PublicationStoreResult"),
    "RelevantPassage": ("src.research.literature", "RelevantPassage"),
    "extract_relevant_passages": ("src.research.literature", "extract_relevant_passages"),
    "plan_literature_requests": ("src.research.literature", "plan_literature_requests"),
    "publications_from_result": ("src.research.literature", "publications_from_result"),
    "ExternalContextSnapshot": ("src.research.external_context", "ExternalContextSnapshot"),
    "ExternalContextSnapshotStore": ("src.research.external_context", "ExternalContextSnapshotStore"),
    "ExternalSnapshotStatus": ("src.research.external_context", "ExternalSnapshotStatus"),
    "SnapshotRestoreResult": ("src.research.external_context", "SnapshotRestoreResult"),
    "SnapshotRestoreStatus": ("src.research.external_context", "SnapshotRestoreStatus"),
    "SnapshotSaveResult": ("src.research.external_context", "SnapshotSaveResult"),
    "create_external_context_snapshot": ("src.research.external_context", "create_external_context_snapshot"),
    "refresh_external_context_snapshot": ("src.research.external_context", "refresh_external_context_snapshot"),
    "ExternalContextUIState": ("src.research.presenter", "ExternalContextUIState"),
    "ExternalDisplayState": ("src.research.presenter", "ExternalDisplayState"),
    "PresenterDeck": ("src.research.presenter", "PresenterDeck"),
    "PresenterSection": ("src.research.presenter", "PresenterSection"),
    "PresenterSectionKind": ("src.research.presenter", "PresenterSectionKind"),
    "build_external_context_ui_state": ("src.research.presenter", "build_external_context_ui_state"),
    "build_presenter_deck": ("src.research.presenter", "build_presenter_deck"),
    "render_external_context_status": ("src.research.presenter", "render_external_context_status"),
    "render_presenter_mode": ("src.research.presenter", "render_presenter_mode"),
    "BriefFormat": ("src.research.brief", "BriefFormat"),
    "ResearchBrief": ("src.research.brief", "ResearchBrief"),
    "build_research_brief": ("src.research.brief", "build_research_brief"),
    "render_research_brief_download": ("src.research.brief", "render_research_brief_download"),
    "ResearchExportFormat": ("src.research.ui", "ResearchExportFormat"),
    "build_research_export_tables": ("src.research.ui", "build_research_export_tables"),
    "research_explorer_download_payload": ("src.research.ui", "research_explorer_download_payload"),
    "render_research_explorer_download": ("src.research.ui", "render_research_explorer_download"),
    "DiseaseBankService": ("src.research.disease_bank.service", "DiseaseBankService"),
    "DiseaseLoad": ("src.research.disease_bank.service", "DiseaseLoad"),
    "DiseaseBenchmark": ("src.research.disease_bank.benchmark", "DiseaseBenchmark"),
    "benchmark_disease_reference": ("src.research.disease_bank.benchmark", "benchmark_disease_reference"),
    "build_interpreter_context": ("src.research.interpreter", "build_interpreter_context"),
    "DeterministicInterpreter": ("src.research.interpreter", "DeterministicInterpreter"),
    "InterpreterContext": ("src.research.interpreter", "InterpreterContext"),
    "InterpreterRender": ("src.research.interpreter", "InterpreterRender"),
    "ResearchInterpreter": ("src.research.interpreter", "ResearchInterpreter"),
    "CURRENT_RESULT_ORDER_NOTE": ("src.research.query_planner", "CURRENT_RESULT_ORDER_NOTE"),
    "LookupIntent": ("src.research.query_planner", "LookupIntent"),
    "LookupOperation": ("src.research.query_planner", "LookupOperation"),
    "QUERY_PLANNER_VERSION": ("src.research.query_planner", "QUERY_PLANNER_VERSION"),
    "QueryPlan": ("src.research.query_planner", "QueryPlan"),
    "QueryPlanner": ("src.research.query_planner", "QueryPlanner"),
    "RetrievalLevel": ("src.research.query_planner", "RetrievalLevel"),
    "plan_research_queries": ("src.research.query_planner", "plan_research_queries"),
    "EntityResolver": ("src.research.entity_resolution", "EntityResolver"),
    "ExternalCandidateSet": ("src.research.entity_resolution", "ExternalCandidateSet"),
    "ProjectEntityIndex": ("src.research.entity_resolution", "ProjectEntityIndex"),
    "empty_project_index": ("src.research.entity_resolution", "empty_project_index"),
    "external_cache_key": ("src.research.entity_resolution", "external_cache_key"),
    "load_human_project_index": ("src.research.entity_resolution", "load_human_project_index"),
    "load_mouse_project_index": ("src.research.entity_resolution", "load_mouse_project_index"),
}

__all__ = sorted(_LAZY_EXPORTS)


def __getattr__(name: str) -> Any:
    target = _LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attribute_name = target
    value = getattr(import_module(module_name), attribute_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))


if TYPE_CHECKING:
    from .config import (
        CacheTTLs,
        PresentationLimits,
        ProviderTimeout,
        ProviderTimeouts,
        ResearchConfig,
        RetrievalBudget,
    )
    from .entity_resolution import (
        EntityResolver,
        ExternalCandidateSet,
        ProjectEntityIndex,
    )
    from .provider_cache import (
        CacheEvidenceClass,
        CacheFirstProviderClient,
        ProviderCache,
        ProviderCacheEntry,
        ProviderCacheKey,
        ProviderCacheStats,
        ProviderTTLPolicy,
    )
    from .provider_runtime import (
        EvidenceProvider,
        ProviderCapability,
        ProviderExecutionPolicy,
        ProviderExecutor,
        ProviderRequest,
        ProviderResult,
        ProviderStatus,
        SessionCircuitBreaker,
    )
    from .external_biology import (
        ExternalBiologyLoad,
        ExternalBiologyLoader,
        ExternalBiologyPlan,
        plan_external_biology,
    )
    from .literature import (
        LiteratureLoad,
        LiteratureLoader,
        LiteraturePlan,
        LiteratureQuery,
        LiteratureQueryBuilder,
        LiteratureRepository,
        LiteratureSearchHit,
        LiteratureSearchResult,
        Publication,
        PublicationMatch,
        PublicationStoreResult,
        RelevantPassage,
        extract_relevant_passages,
        plan_literature_requests,
        publications_from_result,
    )
    from .external_context import (
        ExternalContextSnapshot,
        ExternalContextSnapshotStore,
        ExternalSnapshotStatus,
        SnapshotRestoreResult,
        SnapshotRestoreStatus,
        SnapshotSaveResult,
        create_external_context_snapshot,
        refresh_external_context_snapshot,
    )
    from .presenter import (
        ExternalContextUIState,
        ExternalDisplayState,
        PresenterDeck,
        PresenterSection,
        PresenterSectionKind,
        build_external_context_ui_state,
        build_presenter_deck,
        render_external_context_status,
        render_presenter_mode,
    )
    from .brief import (
        BriefFormat,
        ResearchBrief,
        build_research_brief,
        render_research_brief_download,
    )
    from .ui import (
        ResearchExportFormat,
        build_research_export_tables,
        research_explorer_download_payload,
        render_research_explorer_download,
    )
    from .interpreter import (
        DeterministicInterpreter,
        InterpreterContext,
        InterpreterRender,
        ResearchInterpreter,
        build_interpreter_context,
    )
    from .query_planner import (
        CURRENT_RESULT_ORDER_NOTE,
        LookupIntent,
        LookupOperation,
        QUERY_PLANNER_VERSION,
        QueryPlan,
        QueryPlanner,
        RetrievalLevel,
        plan_research_queries,
    )
    from .models import (
        AliasType,
        Entity,
        EntityAlias,
        EntityResolution,
        EntityType,
        ResolutionStatus,
    )
