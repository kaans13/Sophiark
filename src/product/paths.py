"""Portable project paths with explicit, optional environment overrides."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


def _resolved(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


@dataclass(frozen=True, slots=True)
class ProjectPaths:
    root: Path
    data: Path
    raw: Path
    processed: Path
    derived: Path
    cache: Path
    outputs: Path
    history: Path
    config: Path

    @classmethod
    def discover(cls, root: str | Path | None = None) -> "ProjectPaths":
        project_root = _resolved(
            root
            or os.environ.get("SOPHIARK_ROOT")
            or Path(__file__).resolve().parents[2]
        )
        data = _resolved(os.environ.get("SOPHIARK_DATA_DIR", project_root / "data"))
        outputs = _resolved(os.environ.get("SOPHIARK_OUTPUT_DIR", project_root / "outputs"))
        cache = _resolved(os.environ.get("SOPHIARK_CACHE_DIR", outputs / "cache"))
        history = _resolved(os.environ.get("SOPHIARK_HISTORY_DIR", outputs / "history"))
        return cls(
            root=project_root,
            data=data,
            raw=data / "raw",
            processed=data / "processed",
            derived=data / "derived",
            cache=cache,
            outputs=outputs,
            history=history,
            config=project_root / "config",
        )

    def ensure_runtime_dirs(self) -> None:
        """Create only Sophiark-owned mutable directories; raw data stays immutable."""
        for path in (self.processed, self.derived, self.cache, self.outputs, self.history):
            path.mkdir(parents=True, exist_ok=True)
