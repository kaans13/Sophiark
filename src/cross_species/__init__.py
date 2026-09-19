from .comparison import ComparisonConfig, compare_ortholog_responses
from .entity_resolution import MouseEntityResolution, resolve_mouse_entity
from .orthology import OrthologRecord, OrthologyClass, OrthologyIndex
from .proxy import FunctionalProxyCandidate, discover_functional_proxies
from .suitability import SuitabilityResult, SuitabilityStatus, check_mouse_suitability

__all__ = [
    "ComparisonConfig",
    "FunctionalProxyCandidate",
    "MouseEntityResolution",
    "OrthologRecord",
    "OrthologyClass",
    "OrthologyIndex",
    "SuitabilityResult",
    "SuitabilityStatus",
    "check_mouse_suitability",
    "compare_ortholog_responses",
    "discover_functional_proxies",
    "resolve_mouse_entity",
]
