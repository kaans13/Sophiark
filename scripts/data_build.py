"""Build biological datasets into deterministic shadow/staging locations only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.product.builders import (
    SourceBuildError, build_corum_shadow, build_hpa_shadow,
    build_omnipath_shadow, build_trrust_shadow,
)
from src.product.paths import ProjectPaths


DEFAULT_SOURCES = {
    "hpa": "data/raw/rna_tissue_hpa.tsv",
    "omnipath": "data/regulatory/omnipath/omnipath_9606.tsv",
    "corum": "data/raw/coreComplexes.txt",
    "trrust": "data/raw/trrust_human.tsv",
}


def _build(dataset: str, *, source: Path, paths: ProjectPaths, output_root: Path):
    if dataset == "hpa":
        return build_hpa_shadow(
            source, mapping_path=paths.processed / "ensp_to_ensg_map.csv",
            output_root=output_root,
        )
    if dataset == "omnipath":
        return build_omnipath_shadow(
            source, project_root=paths.root, output_root=output_root, taxon_id=9606,
        )
    if dataset == "corum":
        return build_corum_shadow(
            source, mapping_path=paths.processed / "ensp_with_symbols.csv",
            output_root=output_root,
        )
    if dataset == "trrust":
        return build_trrust_shadow(
            source, mapping_path=paths.processed / "ensp_with_symbols.csv",
            output_root=output_root,
        )
    raise ValueError(f"unsupported dataset: {dataset}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=(*DEFAULT_SOURCES, "all"), required=True)
    parser.add_argument("--source", type=Path, help="Explicit raw source; valid only for a single dataset")
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", type=Path, help="Defaults to data/derived/staging")
    args = parser.parse_args()
    if args.dataset == "all" and args.source is not None:
        parser.error("--source cannot be combined with --dataset all")
    paths = ProjectPaths.discover(args.project_root)
    output_root = (args.output_root or paths.derived / "staging").resolve()
    selected = tuple(DEFAULT_SOURCES) if args.dataset == "all" else (args.dataset,)
    reports: dict[str, object] = {}
    exit_code = 0
    for dataset in selected:
        source = args.source.resolve() if args.source is not None else paths.root / DEFAULT_SOURCES[dataset]
        try:
            reports[dataset] = _build(
                dataset, source=source, paths=paths, output_root=output_root,
            ).as_dict()
        except (SourceBuildError, FileNotFoundError, ValueError, OSError) as exc:
            reports[dataset] = {
                "dataset_id": dataset, "status": "BUILD_FAILED",
                "errors": [f"{type(exc).__name__}: {exc}"], "source": str(source),
            }
            exit_code = 1
    aggregate = {
        "status": "READY" if exit_code == 0 else "PARTIAL",
        "promotion_performed": False,
        "output_root": str(output_root),
        "datasets": reports,
    }
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "aggregate-build-report.json").write_text(
        json.dumps(aggregate, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8",
    )
    print(json.dumps(aggregate, ensure_ascii=False, indent=2, sort_keys=True))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
