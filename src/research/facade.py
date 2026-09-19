"""Single application import surface for the optional Research Context UI.

Importing this module performs no project-data load and imports no UI toolkit.
Phase 5B render symbols are resolved lazily only when a caller requests them.
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .config import DEFAULT_RESEARCH_CONFIG, PresentationLimits, ResearchConfig
from .local_context import LocalResearchContext
from .models import (
    EvidenceAssertion, ObservationType, ProvenanceKind, ProvenanceRecord,
    SimulationResultSnapshot,
)
from .service import (
    OfflineResearchBundle,
    ResearchContextService,
    ResearchTimings,
    ServiceInitializationMode,
    build_offline_research_bundle,
    build_research_snapshot,
)


_LAZY_UI_EXPORTS = {
    "ExplorerSection": ("src.research.ui", "ExplorerSection"),
    "ResearchExportFormat": ("src.research.ui", "ResearchExportFormat"),
    "ResearchExplorerView": ("src.research.ui", "ResearchExplorerView"),
    "build_research_explorer_view": ("src.research.ui", "build_research_explorer_view"),
    "build_research_export_tables": ("src.research.ui", "build_research_export_tables"),
    "research_explorer_download_payload": ("src.research.ui", "research_explorer_download_payload"),
    "render_research_explorer": ("src.research.ui", "render_research_explorer"),
    "render_research_explorer_download": ("src.research.ui", "render_research_explorer_download"),
    "DiseaseBankService": ("src.research.disease_bank.service", "DiseaseBankService"),
    "DiseaseLoad": ("src.research.disease_bank.service", "DiseaseLoad"),
    "DiseaseBenchmark": ("src.research.disease_bank.benchmark", "DiseaseBenchmark"),
    "benchmark_disease_reference": ("src.research.disease_bank.benchmark", "benchmark_disease_reference"),
    "render_disease_context": ("src.research.disease_bank.ui", "render_disease_context"),
    "build_presenter_deck": ("src.research.presenter", "build_presenter_deck"),
    "build_external_context_ui_state": ("src.research.presenter", "build_external_context_ui_state"),
    "render_presenter_mode": ("src.research.presenter", "render_presenter_mode"),
    "render_external_context_status": ("src.research.presenter", "render_external_context_status"),
    "build_research_brief": ("src.research.brief", "build_research_brief"),
    "render_research_brief_download": ("src.research.brief", "render_research_brief_download"),
    "build_interpreter_context": ("src.research.interpreter", "build_interpreter_context"),
    "DeterministicInterpreter": ("src.research.interpreter", "DeterministicInterpreter"),
    "ResearchInterpreter": ("src.research.interpreter", "ResearchInterpreter"),
}


def __getattr__(name: str) -> Any:
    target = _LAZY_UI_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attribute_name = target
    value = getattr(import_module(module_name), attribute_name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))


__all__ = [
    "DEFAULT_RESEARCH_CONFIG",
    "DiseaseBankService",
    "DiseaseBenchmark",
    "DiseaseLoad",
    "DeterministicInterpreter",
    "EvidenceAssertion",
    "ExplorerSection",
    "LocalResearchContext",
    "ObservationType",
    "ProvenanceKind",
    "ProvenanceRecord",
    "OfflineResearchBundle",
    "PresentationLimits",
    "ResearchConfig",
    "ResearchContextService",
    "ResearchInterpreter",
    "ResearchExportFormat",
    "ResearchExplorerView",
    "ResearchTimings",
    "ServiceInitializationMode",
    "SimulationResultSnapshot",
    "build_offline_research_bundle",
    "build_external_context_ui_state",
    "build_interpreter_context",
    "build_presenter_deck",
    "build_research_brief",
    "build_research_export_tables",
    "build_research_explorer_view",
    "build_research_snapshot",
    "benchmark_disease_reference",
    "render_research_explorer",
    "render_research_explorer_download",
    "render_external_context_status",
    "render_disease_context",
    "render_presenter_mode",
    "render_research_brief_download",
    "research_explorer_download_payload",
]


if TYPE_CHECKING:
    from .ui import (
        ExplorerSection,
        ResearchExportFormat,
        ResearchExplorerView,
        build_research_export_tables,
        research_explorer_download_payload,
        render_research_explorer_download,
    )
    from .presenter import (
        build_external_context_ui_state,
        build_presenter_deck,
        render_external_context_status,
        render_presenter_mode,
    )
    from .brief import build_research_brief, render_research_brief_download
    from .interpreter import (
        DeterministicInterpreter,
        ResearchInterpreter,
        build_interpreter_context,
    )
