"""Perturbasyon senaryolarının merkezi ve açık tanımı."""

from __future__ import annotations

from dataclasses import dataclass


SEVERITY_FRACTIONS: dict[str, float] = {
    "Mild": 0.50,
    "Moderate": 0.20,
    "Strong": 0.05,
    "Near-complete": 0.001,
}

NULL_MODEL_PRESETS: dict[str, int] = {
    "Fast": 19,
    "Standard": 99,
    "Rigorous": 499,
}


def severity_fraction(name: str) -> float:
    try:
        return SEVERITY_FRACTIONS[name]
    except KeyError as exc:
        raise ValueError(f"Bilinmeyen perturbasyon şiddeti: {name!r}") from exc


def severity_name(fraction: float, *, tolerance: float = 1e-12) -> str:
    for name, configured in SEVERITY_FRACTIONS.items():
        if abs(float(fraction) - configured) <= tolerance:
            return name
    return "Custom"


@dataclass(frozen=True)
class PerturbationConfig:
    attenuation_fraction: float = SEVERITY_FRACTIONS["Near-complete"]
    severity: str = "Near-complete"
    redistribution_mode: str = "null_fdr"
    exploratory_top_n: int = 250
    percentile: float = 99.0
    null_iterations: int = NULL_MODEL_PRESETS["Standard"]
    fdr_alpha: float = 0.05
    random_seed: int = 42
    efficiency_sample_sources: int = 256
    local_efficiency_sample_nodes: int = 128

    def __post_init__(self) -> None:
        if not 0 < self.attenuation_fraction <= 1:
            raise ValueError("attenuation_fraction (0, 1] aralığında olmalıdır")
        if self.redistribution_mode not in {"null_fdr", "top_n", "percentile"}:
            raise ValueError("redistribution_mode null_fdr, top_n veya percentile olmalıdır")
        if self.exploratory_top_n < 1:
            raise ValueError("exploratory_top_n pozitif olmalıdır")
        if not 0 < self.percentile <= 100:
            raise ValueError("percentile (0, 100] aralığında olmalıdır")
        if self.null_iterations < 0:
            raise ValueError("null_iterations negatif olamaz")
