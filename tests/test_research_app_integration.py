"""Static and fake-state acceptance tests for the minimal app integration."""

from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest

from src import state as app_state
from src.core import cache_manager


ROOT = Path(__file__).resolve().parents[1]
APP_PATH = ROOT / "app.py"


def _app_source() -> str:
    return APP_PATH.read_text(encoding="utf-8")


def _app_tree() -> ast.Module:
    return ast.parse(_app_source())


def _function(tree: ast.Module, name: str) -> ast.FunctionDef:
    return next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    )


def _calls(node: ast.AST, name: str) -> tuple[ast.Call, ...]:
    return tuple(
        item for item in ast.walk(node)
        if isinstance(item, ast.Call)
        and (
            isinstance(item.func, ast.Name) and item.func.id == name
            or isinstance(item.func, ast.Attribute) and item.func.attr == name
        )
    )


class ResearchAppIntegrationTests(unittest.TestCase):
    def test_app_imports_only_the_research_facade(self) -> None:
        imports = tuple(
            node for node in _app_tree().body
            if isinstance(node, ast.ImportFrom)
            and node.module is not None
            and node.module.startswith("src.research")
        )
        self.assertTrue(imports)
        self.assertEqual({node.module for node in imports}, {"src.research.facade"})
        imported_names = {alias.name for node in imports for alias in node.names}
        self.assertTrue({
            "ResearchConfig", "ResearchContextService", "build_research_snapshot",
            "render_research_explorer",
        } <= imported_names)

    def test_research_surface_is_between_deep_interpretation_and_raw_audit(self) -> None:
        source = _app_source()
        deep = source.index("_deep_interpretation_surface = st.container()")
        research = source.index("_research_surface = st.container()")
        audit = source.index("_data_audit_surface = st.container()")
        self.assertLess(deep, research)
        self.assertLess(research, audit)
        self.assertLess(source.index("with _research_surface:"), source.index("with _data_audit_surface:"))

    def test_snapshot_uses_engine_specific_signed_source_and_explicit_scope(self) -> None:
        calls = _calls(_app_tree(), "build_research_snapshot")
        self.assertEqual(len(calls), 1)
        call = calls[0]
        self.assertGreaterEqual(len(call.args), 3)
        self.assertIsInstance(call.args[0], ast.Name)
        self.assertEqual(call.args[0].id, "df")
        self.assertIsInstance(call.args[1], ast.Name)
        self.assertEqual(call.args[1].id, "_research_signed_input")
        self.assertNotIn("signed_redistribution_result", ast.unparse(call))
        source = _app_source()
        self.assertIn(
            '_active_flow.signed_response\n            if _active_engine == "directed"\n            else _negative_loss_df',
            source,
        )
        keyword_nodes = {keyword.arg: keyword.value for keyword in call.keywords}
        keywords = set(keyword_nodes)
        self.assertTrue({
            "species", "taxon_id", "tissue", "targets", "attenuation",
            "selection_mode", "test_limit", "tested_count", "returned_count", "top_n",
            "threshold_parameters",
        } <= keywords)
        self.assertIsInstance(keyword_nodes["test_limit"], ast.Constant)
        self.assertIsNone(keyword_nodes["test_limit"].value)
        self.assertEqual(ast.unparse(keyword_nodes["tested_count"]), "_research_tested_count")
        self.assertEqual(
            ast.unparse(keyword_nodes["returned_count"]),
            "int(_research_unique_gene_count(stresli_df) or 0)",
        )
        self.assertEqual(ast.unparse(keyword_nodes["top_n"]), "_research_top_n")
        self.assertEqual(
            ast.unparse(keyword_nodes["threshold_parameters"]),
            "_research_threshold_parameters",
        )

    def test_simulation_time_top_n_is_frozen_and_not_read_from_results_widget(self) -> None:
        source = _app_source()
        execute_call = _calls(_app_tree(), "execute_simulation")[0]
        execute_keywords = {keyword.arg: keyword.value for keyword in execute_call.keywords}
        scientific_options = execute_keywords["scientific_options"]
        scientific_options_source = ast.unparse(scientific_options)
        self.assertIn("'exploratory_top_n': _simulation_exploratory_top_n", scientific_options_source)
        self.assertIn("'tissue': hedef_doku if doku_aktif else 'None'", scientific_options_source)
        self.assertIn("'tissue_normalization_mode'", scientific_options_source)
        self.assertIn("app_state.store_research_simulation_scope(", source)
        self.assertIn("top_n=(\n                        _simulation_exploratory_top_n", source)
        research_block = source[source.index("with _research_surface:"):source.index("df_goster =")]
        self.assertNotIn("redistribution_candidate_limit", research_block)
        self.assertIn("_research_top_n_from_report(df)", research_block)
        self.assertIn('"scientific_top_n_operator": ">"', source)
        self.assertIn('"scientific_top_n_positive_response_only": True', source)
        self.assertIn('_research_unique_gene_count(_simulation_signed)', source)
        self.assertIn('_research_run_scope.get("tested_count")', source)
        self.assertIn('"loss_presentation_operator": "< -abs(threshold)"', source)
        self.assertIn('"loss_presentation_only": True', research_block)

    def test_legacy_summary_reads_only_canonical_efficiency_fields(self) -> None:
        source = _app_source()
        readout = source[source.index("def _safe_float"):source.index("if stresli_genler:")]
        self.assertIn('"Global_Efficiency_Change_Pct"', readout)
        self.assertIn('"Mean_Local_Efficiency_Change_Pct"', readout)
        self.assertIn('"Target_Local_Efficiency_Change_Pct"', readout)
        self.assertNotIn('"Efficiency_Kayip_Pct"', readout)
        self.assertNotIn('"Local_Efficiency_Kayip_Pct"', readout)

    def test_service_is_resource_cached_and_bundle_auto_loads_once_per_snapshot(self) -> None:
        tree = _app_tree()
        cached = _function(tree, "_cached_research_service")
        self.assertTrue(any(
            isinstance(decorator, ast.Call)
            and isinstance(decorator.func, ast.Attribute)
            and decorator.func.attr == "cache_resource"
            for decorator in cached.decorator_list
        ))

        renderer = _function(tree, "_render_research_surface")
        renderer_source = ast.unparse(renderer)
        self.assertIn("if bundle is None", renderer_source)
        self.assertIn("_cached_research_service", renderer_source)
        self.assertIn("build_offline_bundle", renderer_source)
        self.assertIn("store_research_bundle", renderer_source)
        self.assertNotIn("explore_research_context", renderer_source)
        self.assertTrue(_calls(renderer, "render_research_explorer"))
        self.assertNotIn("requests", renderer_source.casefold())

    def test_snapshot_change_evicts_stale_bundle_but_same_snapshot_reuses_it(self) -> None:
        original_st = app_state.st
        app_state.st = SimpleNamespace(session_state={})
        try:
            app_state.initialize()
            self.assertFalse(app_state.bind_research_snapshot(
                snapshot_id="snapshot-one", simulation_id="simulation-one",
            ))
            bundle = SimpleNamespace(snapshot_id="snapshot-one", simulation_id="simulation-one")
            app_state.store_research_bundle(bundle)
            app_state.st.session_state["research_external_snapshot_id"] = "external-one"

            self.assertFalse(app_state.bind_research_snapshot(
                snapshot_id="snapshot-one", simulation_id="simulation-one",
            ))
            self.assertIs(app_state.research_bundle_for("snapshot-one"), bundle)

            self.assertTrue(app_state.bind_research_snapshot(
                snapshot_id="snapshot-two", simulation_id="simulation-two",
            ))
            self.assertIsNone(app_state.research_bundle_for("snapshot-one"))
            self.assertIsNone(app_state.st.session_state["research_offline_bundle"])
            self.assertIsNone(app_state.st.session_state["research_external_snapshot_id"])
            self.assertFalse(app_state.st.session_state["research_explorer_open"])
            with self.assertRaises(ValueError):
                app_state.store_research_bundle(bundle)
        finally:
            app_state.st = original_st

    def test_simulation_scope_is_defensively_stored_and_cleared(self) -> None:
        original_st = app_state.st
        app_state.st = SimpleNamespace(session_state={})
        thresholds = {"threshold": 0.05, "scientific_top_n_operator": ">"}
        try:
            app_state.initialize()
            app_state.store_research_simulation_scope(
                species="Homo sapiens",
                taxon_id=9606,
                tissue="Pancreas",
                targets=("ENSP_TARGET",),
                attenuation=0.001,
                selection_mode="top_n",
                test_limit=None,
                tested_count=199,
                returned_count=20,
                top_n=20,
                threshold_parameters=thresholds,
            )
            thresholds["threshold"] = 99.0
            scope = app_state.research_simulation_scope()
            self.assertEqual(scope["top_n"], 20)
            self.assertIsNone(scope["test_limit"])
            self.assertEqual(scope["tested_count"], 199)
            self.assertEqual(scope["returned_count"], 20)
            self.assertEqual(scope["threshold_parameters"]["threshold"], 0.05)
            scope["threshold_parameters"]["threshold"] = 88.0
            self.assertEqual(
                app_state.research_simulation_scope()["threshold_parameters"]["threshold"],
                0.05,
            )
            app_state.clear_research_session()
            self.assertIsNone(app_state.research_simulation_scope())
        finally:
            app_state.st = original_st

    def test_reset_clears_research_session_but_preserves_persistent_sqlite(self) -> None:
        research_keys = set(app_state.RESEARCH_SESSION_KEYS)
        self.assertTrue(research_keys <= set(cache_manager._STATE_KEYS_TO_PURGE))
        original_st = cache_manager.st
        original_force_gc = cache_manager.force_gc
        counters = {"data": 0, "resource": 0}

        class _Cache:
            def __init__(self, key: str) -> None:
                self.key = key

            def clear(self) -> None:
                counters[self.key] += 1

        fake_session = {key: object() for key in research_keys}
        cache_manager.st = SimpleNamespace(
            session_state=fake_session,
            cache_data=_Cache("data"),
            cache_resource=_Cache("resource"),
        )
        cache_manager.force_gc = lambda: 0
        try:
            with TemporaryDirectory() as directory:
                sqlite_path = Path(directory) / "research_evidence.sqlite3"
                sqlite_path.write_bytes(b"persistent fixture")
                purged, collected = cache_manager.scorched_earth_reset()
                self.assertEqual(purged, len(research_keys))
                self.assertEqual(collected, 0)
                self.assertEqual(counters, {"data": 1, "resource": 1})
                self.assertTrue(sqlite_path.is_file())
                self.assertEqual(sqlite_path.read_bytes(), b"persistent fixture")
        finally:
            cache_manager.st = original_st
            cache_manager.force_gc = original_force_gc

    def test_existing_raw_export_surface_remains_unchanged(self) -> None:
        source = _app_source()
        self.assertEqual(source.count("render_raw_data(df_goster)"), 1)
        self.assertIn("df_goster = uygula_symbol_map_lokal(df, _gene_to_symbol)", source)
        self.assertNotIn("research", source[source.index("with _data_audit_surface:"):source.index("if _interpretation_result")].casefold())


if __name__ == "__main__":
    unittest.main()
