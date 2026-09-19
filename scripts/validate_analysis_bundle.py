"""Run representative live Unified analyses and measure `.sophiark` roundtrips."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from time import perf_counter
import uuid

from pandas.testing import assert_frame_equal

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.product.analysis_bundle import AnalysisBundleReader, AnalysisBundleWriter, _portable  # noqa: E402
from src.unified import run_unified_analysis  # noqa: E402


TARGETS = {
    "CFTR": "ENSP00000003084",
    "TGFBR1": "ENSP00000447297",
    "ANO2": "ENSP00000348453",
}


def main() -> int:
    output = ROOT / "outputs" / "analysis_bundle_validation"
    output.mkdir(parents=True, exist_ok=True)
    writer = AnalysisBundleWriter()
    reader = AnalysisBundleReader()
    runs = []
    for symbol, target in TARGETS.items():
        started = perf_counter()
        live = run_unified_analysis(target=target, tissue="Lung", bc_sample_sources=4, use_cache=False)
        live_seconds = perf_counter() - started
        path = output / f"{symbol}_Lung.sophiark"
        started = perf_counter()
        written = writer.write(live, path, analysis_id=str(uuid.uuid5(uuid.NAMESPACE_URL, f"sophiark:{symbol}:Lung")))
        write_seconds = perf_counter() - started
        started = perf_counter()
        loaded = reader.load(path)
        load_seconds = perf_counter() - started
        assert loaded.report.target == live.target and loaded.report.tissue == live.tissue
        assert loaded.report.status is live.status
        assert loaded.report.classic_status == live.classic_status
        assert loaded.report.directed_status == live.directed_status
        assert loaded.report.agreement == live.agreement
        assert loaded.report.context_availability == _portable(live.context_availability)
        assert loaded.report.provenance == _portable(live.provenance)
        assert loaded.report.errors == _portable(live.errors)
        for name in ("candidates", "classic_raw", "directed_raw", "complex_context"):
            assert_frame_equal(getattr(loaded.report, name), getattr(live, name), check_exact=True)
        runs.append({
            "symbol": symbol, "target": target, "status": live.status.value,
            "candidate_rows": len(live.candidates), "classic_rows": len(live.classic_raw),
            "directed_rows": len(live.directed_raw), "bundle_bytes": written.bytes_written,
            "live_seconds": live_seconds, "write_seconds": write_seconds,
            "load_seconds": load_seconds, "roundtrip_exact": True,
            "content_fingerprint": written.content_fingerprint,
        })
        print(json.dumps(runs[-1]))
    summary = {"status": "COMPLETE", "runs": runs, "load_reran_engines": False}
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
