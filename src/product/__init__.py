"""Lightweight public facade around Sophiark's frozen scientific engines.

Exports resolve lazily so a focused product import does not load Unified,
Streamlit, datasets, or persistence resources as a side effect.
"""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "Capability": ("capabilities", "Capability"),
    "CapabilityReport": ("capabilities", "CapabilityReport"),
    "detect_capabilities": ("capabilities", "detect_capabilities"),
    "DataRegistry": ("datasets", "DataRegistry"),
    "DatasetHealth": ("datasets", "DatasetHealth"),
    "DatasetSpec": ("datasets", "DatasetSpec"),
    "default_registry": ("datasets", "default_registry"),
    "inspect_build_report": ("datasets", "inspect_build_report"),
    "MappingReport": ("entities", "MappingReport"),
    "normalize_ensp": ("entities", "normalize_ensp"),
    "normalize_string_protein_id": ("entities", "normalize_string_protein_id"),
    "ProjectPaths": ("paths", "ProjectPaths"),
    "AnalysisBundleReader": ("analysis_bundle", "AnalysisBundleReader"),
    "AnalysisBundleWriter": ("analysis_bundle", "AnalysisBundleWriter"),
    "BundleError": ("analysis_bundle", "BundleError"),
    "BundleStatus": ("analysis_bundle", "BundleStatus"),
    "SavedUnifiedResult": ("analysis_bundle", "SavedUnifiedResult"),
    "HistoryService": ("history", "HistoryService"),
    "ImportStatus": ("history", "ImportStatus"),
    "rebuild_history_index": ("history", "rebuild_history_index"),
    "ProductStatus": ("status", "ProductStatus"),
    "SOPHIARK_VERSION": ("versions", "SOPHIARK_VERSION"),
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    try:
        module_name, attribute = _EXPORTS[name]
    except KeyError as error:
        raise AttributeError(name) from error
    value = getattr(import_module(f"{__name__}.{module_name}"), attribute)
    globals()[name] = value
    return value
