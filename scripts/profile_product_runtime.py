"""Measure product/orchestration boundaries without instrumenting frozen engines."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from time import perf_counter
from unittest.mock import patch
import uuid

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.product.analysis_bundle import AnalysisBundleReader, AnalysisBundleWriter
from src.unified.service import clear_unified_cache, run_unified_analysis


TARGETS = {
    "CFTR": "ENSP00000003084",
    "TGFBR1": "ENSP00000447297",
    "ANO2": "ENSP00000348453",
}


def _rss_mb() -> float:
    try:
        import psutil
        return psutil.Process().memory_info().rss / (1024 * 1024)
    except ImportError:
        # Keep the benchmark usable in a minimal declared test environment.
        import ctypes
        from ctypes import wintypes

        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        handle = ctypes.windll.kernel32.OpenProcess(0x1000 | 0x0400, False, os.getpid())
        try:
            if not ctypes.windll.psapi.GetProcessMemoryInfo(
                handle, ctypes.byref(counters), counters.cb,
            ):
                return float("nan")
            return counters.WorkingSetSize / (1024 * 1024)
        finally:
            ctypes.windll.kernel32.CloseHandle(handle)


class BoundaryTimer:
    def __init__(self) -> None:
        self.values: dict[str, float] = {}

    def wrap(self, target: object, attribute: str, label: str, stack: ExitStack) -> None:
        original = getattr(target, attribute)

        def measured(*args, **kwargs):
            started = perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                self.values[label] = self.values.get(label, 0.0) + perf_counter() - started

        stack.enter_context(patch.object(target, attribute, measured))


def _run(symbol: str, output: Path, *, use_cache: bool) -> dict[str, object]:
    import src.biology_logic as classic_logic
    import src.main as classic_main
    import src.services.directed_analysis as directed
    import src.unified.service as unified

    timer = BoundaryTimer()
    memory_before = _rss_mb()
    with ExitStack() as stack:
        timer.wrap(classic_main, "arayuz_icin_motoru_hazirla", "tissue_graph_preparation", stack)
        timer.wrap(classic_logic, "run_infection_simulation", "classic_calculation", stack)
        timer.wrap(directed, "load_signaling_cache", "omnipath_load", stack)
        timer.wrap(directed, "prepare_effective_directed_source_graph", "directed_source_preparation", stack)
        timer.wrap(directed, "_cached_hybrid_graph", "directed_overlay", stack)
        timer.wrap(directed, "safe_run_directed_calculation", "directed_calculation", stack)
        timer.wrap(unified, "_synthesize", "unified_synthesis", stack)
        timer.wrap(unified, "_target_edge_provenance", "evidence_provenance", stack)
        timer.wrap(unified, "_complex_context", "corum_context", stack)
        timer.wrap(unified, "_provenance", "run_provenance", stack)
        started = perf_counter()
        report = run_unified_analysis(
            target=TARGETS[symbol], tissue="Lung", use_cache=use_cache,
            bc_sample_sources=4,
        )
        wall = perf_counter() - started
    memory_after_analysis = _rss_mb()
    bundle = output / f"{symbol}_{uuid.uuid4().hex[:8]}.sophiark"
    started = perf_counter()
    AnalysisBundleWriter().write(
        report, bundle, created_at=datetime.now(timezone.utc),
    )
    bundle_write = perf_counter() - started
    started = perf_counter()
    loaded = AnalysisBundleReader().load(bundle)
    bundle_load = perf_counter() - started
    return {
        "symbol": symbol,
        "target": TARGETS[symbol],
        "status": report.status.value,
        "wall_seconds": wall,
        "boundaries_seconds": timer.values,
        "unattributed_seconds": max(0.0, wall - sum(timer.values.values())),
        "rss_before_mb": memory_before,
        "rss_after_analysis_mb": memory_after_analysis,
        "rss_delta_mb": memory_after_analysis - memory_before,
        "bundle_write_seconds": bundle_write,
        "bundle_load_seconds": bundle_load,
        "bundle_bytes": bundle.stat().st_size,
        "bundle_load_reran_engines": False,
        "roundtrip_status_equal": loaded.report.status is report.status,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    clear_unified_cache()
    runs = []
    for index, symbol in enumerate(TARGETS):
        row = _run(symbol, args.output, use_cache=False)
        row["scenario"] = "cold" if index == 0 else "same_tissue_next_target"
        runs.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
    clear_unified_cache()
    miss = _run("CFTR", args.output, use_cache=True)
    miss["scenario"] = "cache_miss"
    runs.append(miss)
    hit = _run("CFTR", args.output, use_cache=True)
    hit["scenario"] = "cache_hit"
    runs.append(hit)
    payload = {
        "measured_at": datetime.now(timezone.utc).isoformat(),
        "runs": runs,
    }
    (args.output / "profile.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
