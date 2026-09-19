"""File-level drift guard for the scientific reference implementation."""

from __future__ import annotations

import json
from pathlib import Path

from .fingerprints import sha256_file
from .paths import ProjectPaths


class EngineFreezeViolation(RuntimeError):
    pass


def verify_engine_freeze(
    *, root: str | Path | None = None, manifest_path: str | Path | None = None
) -> dict[str, str]:
    paths = ProjectPaths.discover(root)
    manifest_file = Path(manifest_path) if manifest_path else paths.config / "engine-freeze.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    entries = payload.get("files")
    if not isinstance(entries, list) or not entries:
        raise EngineFreezeViolation("scientific engine freeze manifest has no file entries")
    declared = [item.get("path") for item in entries if isinstance(item, dict)]
    if len(declared) != len(entries) or any(not isinstance(path, str) or not path for path in declared):
        raise EngineFreezeViolation("scientific engine freeze manifest contains an invalid file entry")
    duplicates = sorted({path for path in declared if declared.count(path) > 1})
    if duplicates:
        raise EngineFreezeViolation(
            "scientific engine freeze manifest contains duplicate entries: " + ", ".join(duplicates)
        )
    violations: dict[str, str] = {}
    for item in entries:
        relative = item["path"]
        relative_path = Path(relative)
        if (
            relative_path.is_absolute()
            or ".." in relative_path.parts
            or not relative_path.parts
            or relative_path.parts[0] != "src"
            or any(part.casefold() in {"cache", "outputs", "output", "generated", "__pycache__"}
                   for part in relative_path.parts)
        ):
            violations[relative] = "unsafe/non-source manifest path"
            continue
        expected_hash = item.get("sha256")
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            violations[relative] = "invalid sha256 declaration"
            continue
        source = paths.root / relative_path
        if not source.is_file():
            violations[relative] = "missing"
            continue
        actual = sha256_file(source)
        if actual != expected_hash:
            violations[relative] = actual
    if violations:
        detail = ", ".join(f"{path}={value}" for path, value in violations.items())
        raise EngineFreezeViolation(f"scientific engine freeze violated: {detail}")
    return {item["path"]: item["sha256"] for item in entries}
