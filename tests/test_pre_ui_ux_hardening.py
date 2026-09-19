from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import io
import re
import tomllib
from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pandas as pd

from src.product.runtime_resources import (
    DetachedPreparedContextCache, clear_runtime_resource_caches, symbol_mapping,
)
from src.ui.metric_presentation import (
    display_status, format_metric, metric_label, ordered_columns, prepare_display_frame,
)
from src.runtime_logging import ensure_utf8_stream
from src.product.analysis_bundle import AnalysisBundleReader, AnalysisBundleWriter
from src.unified.service import UnifiedResearchReport, UnifiedStatus


class _Graph:
    def __init__(self, names: tuple[str, ...] = ("A", "B")) -> None:
        self.vs = {"name": list(names)}

    def copy(self):
        return _Graph(tuple(self.vs["name"]))


def test_detached_context_cache_reuses_builder_without_sharing_working_state():
    cache = DetachedPreparedContextCache(maxsize=2)
    calls = 0

    def build():
        nonlocal calls
        calls += 1
        return _Graph(), pd.DataFrame({"value": [1.0]})

    first_graph, first_scores = cache.get_or_build(("Lung", "fingerprint"), build)
    first_graph.vs["name"].append("MUTATED")
    first_scores.loc[0, "value"] = 99.0
    second_graph, second_scores = cache.get_or_build(("Lung", "fingerprint"), build)
    assert calls == 1
    assert second_graph.vs["name"] == ["A", "B"]
    assert second_scores.loc[0, "value"] == 1.0


def test_detached_context_cache_resets_engine_rng_before_each_new_build(monkeypatch):
    import src.product.runtime_resources as resources

    resets = []
    monkeypatch.setattr(resources, "_reset_igraph_rng", lambda seed=42: resets.append(seed))
    cache = DetachedPreparedContextCache(maxsize=2)

    def build():
        return _Graph(), pd.DataFrame({"value": [1.0]})

    cache.get_or_build(("Lung", "v1"), build)
    cache.get_or_build(("Lung", "v1"), build)
    cache.get_or_build(("Liver", "v1"), build)

    assert resets == [42, 42]


def test_detached_context_cache_coalesces_concurrent_same_key_builds():
    cache = DetachedPreparedContextCache(maxsize=2)
    entered = Event()
    release = Event()
    calls = 0

    def build():
        nonlocal calls
        calls += 1
        entered.set()
        assert release.wait(2)
        return _Graph(), pd.DataFrame({"value": [1.0]})

    with ThreadPoolExecutor(max_workers=2) as pool:
        first = pool.submit(cache.get_or_build, ("Adipose",), build)
        assert entered.wait(1)
        second = pool.submit(cache.get_or_build, ("Adipose",), build)
        release.set()
        first_graph, _ = first.result(timeout=2)
        second_graph, _ = second.result(timeout=2)

    assert calls == 1
    assert first_graph is not second_graph


def test_detached_context_cache_serializes_distinct_igraph_builds():
    cache = DetachedPreparedContextCache(maxsize=2)
    first_entered = Event()
    release_first = Event()
    second_entered = Event()

    def first_build():
        first_entered.set()
        assert release_first.wait(2)
        return _Graph(), pd.DataFrame({"value": [1.0]})

    def second_build():
        second_entered.set()
        return _Graph(), pd.DataFrame({"value": [2.0]})

    with ThreadPoolExecutor(max_workers=2) as pool:
        first = pool.submit(cache.get_or_build, ("Lung",), first_build)
        assert first_entered.wait(1)
        second = pool.submit(cache.get_or_build, ("Liver",), second_build)
        assert not second_entered.wait(.1)
        release_first.set()
        first.result(timeout=2)
        _, scores = second.result(timeout=2)

    assert second_entered.is_set()
    assert scores.loc[0, "value"] == 2.0


def test_symbol_mapping_cache_invalidates_on_file_identity_change(tmp_path: Path):
    clear_runtime_resource_caches()
    path = tmp_path / "mapping.csv"
    path.write_text("gene,Symbol\nENSP1,A\n", encoding="utf-8")
    assert symbol_mapping(path) == {"ENSP1": "A"}
    path.write_text("gene,Symbol\nENSP1,B\nENSP2,C\n", encoding="utf-8")
    assert symbol_mapping(path) == {"ENSP1": "B", "ENSP2": "C"}


def test_presentation_keeps_machine_schema_and_distinguishes_semantics():
    source = pd.DataFrame({
        "Evidence_Experimental_Support_Present": pd.Series([True, False, pd.NA], dtype="boolean"),
        "classic_response": [0.0, -1.25, 2.5],
        "entity_id": ["E1", "E2", "E3"],
        "symbol": ["A", "B", "C"],
        "finding_class": ["ROBUST_CROSS_MODEL", "DIRECTION_SENSITIVE", pd.NA],
    })
    before = source.copy(deep=True)
    display = prepare_display_frame(source)
    assert list(display.columns[:2]) == ["symbol", "entity_id"]
    assert list(display["Evidence_Experimental_Support_Present"]) == [
        "Evet", "Hayır", "Kullanılamıyor",
    ]
    assert display.loc[2, "finding_class"] == "—"
    pd.testing.assert_frame_equal(source, before)
    assert list(source.columns) != ordered_columns(source.columns)


def test_percent_and_status_formatting_are_explicit_not_magnitude_inferred():
    assert metric_label("symbol") == "Gen"
    assert metric_label("entity_id") == "ENSP"
    assert metric_label("Symbol") == "Gen"
    assert metric_label("gene") == "ENSP"
    assert format_metric("classic_response", 0.42) == "+0.4200%"
    assert format_metric("sign_agreement_rate", 0.42) == "42.0%"
    assert format_metric("Evidence_Experimental_Support_Present", False) == "Hayır"
    assert format_metric("Evidence_Experimental_Support_Present", pd.NA) == "Kullanılamıyor"
    assert display_status("ENGINE_PARTIAL") == "Motor kısmi"
    assert display_status("UNSUPPORTED_BUNDLE_VERSION") == "Desteklenmeyen bundle sürümü"
    assert display_status("ID_CONFLICT") == "Kimlik çakışması"


def test_cp1254_stream_is_safely_reconfigured_for_unicode_logs():
    raw = io.BytesIO()
    stream = io.TextIOWrapper(raw, encoding="cp1254", errors="strict")
    ensure_utf8_stream(stream)
    stream.write("doku ağırlığı → hazır ✓")
    stream.flush()
    assert raw.getvalue().decode("utf-8") == "doku ağırlığı → hazır ✓"


def test_analysis_service_import_is_streamlit_free():
    import src.services.analysis_service as service

    assert "streamlit" not in service.__dict__
    assert isinstance(service._PREPARED_NETWORKS, DetachedPreparedContextCache)


def test_unified_same_tissue_targets_reuse_detached_prepared_context(monkeypatch):
    import src.biology_logic as classic_logic
    import src.main as classic_main
    import src.services.directed_analysis as directed_service
    import src.unified.service as unified

    unified.clear_unified_cache()
    preparations = 0

    def prepare(**_):
        nonlocal preparations
        preparations += 1
        return _Graph(("ENSP_A", "ENSP_B")), pd.DataFrame({"gene": ["ENSP_A", "ENSP_B"]})

    classic = pd.DataFrame({
        "gene": ["ENSP_A", "ENSP_B"], "Delta_PageRank_Pct": [1.0, -1.0],
    })
    directed = pd.DataFrame({
        "entity": ["ENSP_A", "ENSP_B"], "Directed_Redistribution_Pct": [1.0, -1.0],
    })
    monkeypatch.setattr(classic_main, "arayuz_icin_motoru_hazirla", prepare)
    monkeypatch.setattr(classic_logic, "run_infection_simulation", lambda *_, **__: None)
    monkeypatch.setattr(classic_logic, "latest_signed_redistribution", lambda: classic.copy())
    monkeypatch.setattr(
        directed_service, "run_optional_directed_engine",
        lambda **_: SimpleNamespace(calculation=SimpleNamespace(
            report=directed.copy(), metadata={}, status=SimpleNamespace(value="AVAILABLE"), error=None,
        )),
    )
    unified._default_engine_runs(
        target="ENSP_A", tissue="Lung", block_strength=.001, damping=.85, bc_sample_sources=4,
    )
    unified._default_engine_runs(
        target="ENSP_B", tissue="Lung", block_strength=.001, damping=.85, bc_sample_sources=4,
    )
    assert preparations == 1


def test_existing_v1_bundle_remains_machine_canonical_and_displayable(tmp_path: Path):
    candidates = pd.DataFrame({
        "entity_id": ["ENSP_A"], "symbol": ["A"],
        "classic_response": [0.42], "finding_class": ["ROBUST_CROSS_MODEL"],
    })
    report = UnifiedResearchReport(
        target="ENSP_A", tissue="Lung", status=UnifiedStatus.COMPLETE,
        classic_status="COMPLETE", directed_status="AVAILABLE", candidates=candidates,
    )
    path = tmp_path / "existing-v1.sophiark"
    AnalysisBundleWriter().write(report, path)
    loaded = AnalysisBundleReader().load(path)
    assert loaded.manifest["bundle_schema_version"] == 1
    assert list(loaded.report.candidates.columns) == list(candidates.columns)
    display = prepare_display_frame(loaded.report.candidates)
    assert list(display.columns[:2]) == ["symbol", "entity_id"]
    assert display.loc[0, "finding_class"] == "Modeller arası sağlam"
    pd.testing.assert_frame_equal(loaded.report.candidates, candidates)


def test_active_runtime_contains_no_user_specific_absolute_path():
    root = Path(__file__).resolve().parents[1]
    candidates = [root / "app.py", *sorted((root / "src").rglob("*.py"))]
    pattern = re.compile(r"(?i)([A-Z]:[\\/]Users[\\/]|OneDrive[\\/])")
    hits = [str(path.relative_to(root)) for path in candidates if pattern.search(path.read_text(encoding="utf-8"))]
    assert hits == []


def test_pyproject_declares_core_and_test_tooling():
    root = Path(__file__).resolve().parents[1]
    payload = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    dependencies = payload["project"]["dependencies"]
    test_dependencies = payload["project"]["optional-dependencies"]["test"]
    assert any(item.startswith("numpy") for item in dependencies)
    assert any(item.startswith("pandas") for item in dependencies)
    assert any(item.startswith("pytest") for item in test_dependencies)
