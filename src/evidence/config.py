from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import hashlib
import json
from pathlib import Path


ENGINE_VERSION = "evidence-beta-1.0.0"
STRING_VERSION = "12.0"
STRING_TAXON = 9606


class EvidenceMode(str, Enum):
    PARITY = "parity"
    E1_WEIGHT_ONLY = "e1_weight_only"
    E2_TOPOLOGY_WEIGHT = "e2_topology_weight"


class PhysicalPolicy(str, Enum):
    CONTEXT_ONLY = "context_only"
    BOUNDED = "bounded"


class CorumPolicy(str, Enum):
    CLASSIC_PAIRWISE = "classic_pairwise"  # parity control only
    COMPLEX_RESPONSE_ONLY = "complex_response_only"


@dataclass(frozen=True, slots=True)
class EvidenceConfig:
    mode: EvidenceMode = EvidenceMode.E1_WEIGHT_ONLY
    s_nontext_policy: str = "validated_string_combiner_excluding_textmining_direct_and_transferred"
    string_prior: float = 0.041
    text_gate: str = "G1"
    hill_n: float = 2.0
    hill_m: float = 2.0
    hill_k_text: float = 0.6
    hill_k_evidence: float = 0.6
    lambda_text: float = 0.05
    physical_policy: PhysicalPolicy = PhysicalPolicy.CONTEXT_ONLY
    lambda_physical: float = 0.0
    corum_policy: CorumPolicy = CorumPolicy.COMPLEX_RESPONSE_ONLY
    string_threshold: int = 500
    engine_version: str = ENGINE_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "mode", EvidenceMode(self.mode))
        object.__setattr__(self, "physical_policy", PhysicalPolicy(self.physical_policy))
        object.__setattr__(self, "corum_policy", CorumPolicy(self.corum_policy))
        if self.mode is EvidenceMode.PARITY and self.corum_policy is not CorumPolicy.CLASSIC_PAIRWISE:
            object.__setattr__(self, "corum_policy", CorumPolicy.CLASSIC_PAIRWISE)
        if self.physical_policy is PhysicalPolicy.CONTEXT_ONLY and self.lambda_physical != 0:
            raise ValueError("context_only physical policy requires lambda_physical=0")
        for name in ("lambda_text", "lambda_physical"):
            if not 0 <= float(getattr(self, name)) <= 0.10:
                raise ValueError(f"{name} must be in [0, 0.10]")
        if self.text_gate not in {"G1", "G2", "G3"}:
            raise ValueError("text_gate must be G1, G2, or G3")

    def identity(self) -> str:
        payload = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()

    def to_dict(self) -> dict[str, object]:
        values = asdict(self)
        values["mode"] = self.mode.value
        values["physical_policy"] = self.physical_policy.value
        values["corum_policy"] = self.corum_policy.value
        return values


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FULL_LINKS_PATH = PROJECT_ROOT / "data" / "raw" / "9606.protein.links.full.v12.0.txt.gz"
EVIDENCE_DB_PATH = PROJECT_ROOT / "data" / "processed" / "string_evidence_v12_9606.sqlite"
CLASSIC_DB_PATH = PROJECT_ROOT / "data" / "raw" / "hinterland_core.db"
PHYSICAL_DB_PATH = PROJECT_ROOT / "data" / "raw" / "hinterland_core_physical.db"
CORUM_PATH = PROJECT_ROOT / "data" / "raw" / "coreComplexes.txt"
