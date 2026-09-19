"""Reusing presentation projections must preserve every exported field."""
from tests.test_research_ui import _bundle
from src.research import ui


def test_export_projection_reuse_is_exactly_equivalent(monkeypatch):
    bundle = _bundle()
    optimized = ui.build_research_export_tables(bundle)
    original_functions = ui._function_items
    original_relationships = ui._relationship_items

    def repeated_functions(*args, **kwargs):
        kwargs.pop("_prepared_rows", None)
        kwargs.pop("_prepared_summaries", None)
        return original_functions(*args, **kwargs)

    def repeated_relationships(*args, **kwargs):
        kwargs.pop("_prepared_summaries", None)
        return original_relationships(*args, **kwargs)

    monkeypatch.setattr(ui, "_function_items", repeated_functions)
    monkeypatch.setattr(ui, "_relationship_items", repeated_relationships)
    assert ui.build_research_export_tables(bundle) == optimized


def test_export_builds_function_rows_once(monkeypatch):
    calls = []
    original = ui._function_rows

    def counted(context):
        calls.append(context)
        return original(context)

    monkeypatch.setattr(ui, "_function_rows", counted)
    ui.build_research_export_tables(_bundle())
    assert len(calls) == 1
