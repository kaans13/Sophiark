"""Print a machine-readable Sophiark dataset/capability health report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.product.capabilities import detect_capabilities
from src.product.datasets import default_registry
from src.product.engine_freeze import verify_engine_freeze
from src.product.paths import ProjectPaths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", help="Project root; defaults to SOPHIARK_ROOT or auto-discovery")
    parser.add_argument("--fingerprint", action="store_true", help="Hash full dataset files (may be slow)")
    args = parser.parse_args()
    paths = ProjectPaths.discover(args.root)
    registry = default_registry(paths)
    report = detect_capabilities(registry).as_dict()
    if args.fingerprint:
        report["datasets"] = {
            key: value.as_dict() for key, value in registry.health_report(fingerprint=True).items()
        }
    report["overall_status"] = registry.overall_status().value
    report["engine_freeze_files"] = len(verify_engine_freeze(root=paths.root))
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
